#!/usr/bin/env python3
"""Classify seed inputs, candidate loss and the executed operations on final paths."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from local_truth import envelope, category
from pilot import archive_records, event_groups, read_records, process_event, pnorm, match, good_track
from seed_geometry import taubin
from artifact_io import publish_gzip


def finite_q(value):
    return value is not None and math.isfinite(value) and value != 0


def ratio(q, truth_q):
    return abs(truth_q/q) if finite_q(q) else None


def seed_diagnostics(seed, circle, truth, records):
    p = pnorm(seed)
    pt = np.linalg.norm(truth['momentum_F_MeV'])*.001
    barcode = truth['barcode']
    pdg = truth['pdg']
    result = dict(seed=seed['seed'], p_ratio=p/pt if p is not None else None,
                  seed_majority_other=match(seed) != barcode,
                  stations=circle.get('stations'), covariance_diagonal=np.diag(seed['cov']).tolist())
    bounds = envelope(records, barcode, pdg, seed['position_A_mm'][0])
    result['local_class'] = category(p, bounds)
    segments = circle.get('segments', [])
    result['segment_majority_other'] = any(not s['matches'] or s['matches'][0][0] != barcode for s in segments)
    result['segment_mixed_deposits'] = any(len(s['matches']) > 1 for s in segments)
    result['seed_ambiguous_70pct'] = len(seed['matches']) > 1 and seed['matches'][1][1] >= .7*len(seed['hits'])
    slope_errors = []
    center_errors = []
    truth_fake = []
    for s in segments:
        station = s['station']
        direction = truth['station_direction_F'][station]
        position = truth['station_position_F_mm'][station]
        if direction is None or position is None or not direction[2]:
            continue
        center = np.asarray(s['position_F_mm'])
        local_truth = np.asarray(position) + np.asarray(direction)*((center[2]-position[2])/direction[2])
        center_errors.append(float(center[1]-local_truth[1]))
        momentum = s['momentum_F']
        if momentum[2]:
            slope_errors.append(float(momentum[1]/momentum[2]-direction[1]/direction[2]))
        truth_fake.extend([center+offset*np.asarray(direction) for offset in [0., -30., 30.]])
    result['segment_y_slope_errors'] = slope_errors
    result['segment_y_center_errors_mm'] = center_errors
    points = circle.get('fake_positions_F_mm')
    if points is not None and len(points) >= 3 and circle.get('stations', 0) > 1:
        fit = taubin(np.asarray(points)[:, [2, 1]])
        result.update(circle_valid=fit['valid'], circle_condition=fit.get('covariance_condition'),
                      normalized_determinant=fit.get('normalized_determinant'))
        production = circle.get('radius_mm')
        if fit['valid'] and production and math.isfinite(production):
            result['radius_recompute_relative_error'] = fit['radius_mm']/production-1
        for label, values in [('centers_only', [s['position_F_mm'] for s in segments]),
                              ('truth_directions', truth_fake)]:
            if len(values) < 3:
                continue
            alternative = taubin(np.asarray(values)[:, [2, 1]])
            if alternative['valid']:
                result[label+'_p_ratio'] = .000165*alternative['radius_mm']/pt
                result[label+'_charge'] = 1 if alternative['center'][1] < 0 else -1
    elif circle.get('stations') == 1:
        result['single_station_default'] = True
    if p is None:
        result['primary_class'] = 'invalid_seed_momentum'
    elif result['segment_majority_other']:
        result['primary_class'] = 'segment_majority_other'
    elif result['local_class'] == 'excluded_local_tail':
        result['primary_class'] = 'excluded_local_tail'
    elif result.get('single_station_default'):
        result['primary_class'] = 'single_station_default'
    elif result.get('circle_valid') is False:
        result['primary_class'] = 'invalid_circle'
    elif result.get('truth_directions_p_ratio', 0.) >= .5 and p/pt < .5:
        result['primary_class'] = 'sensitive_to_segment_direction'
    elif result.get('centers_only_p_ratio', 0.) >= .5 and p/pt < .5:
        result['primary_class'] = 'sensitive_to_fake_point_weight'
    else:
        result['primary_class'] = result['local_class']
    return result


def operations(trace, initial_q, prefix):
    steps = trace['steps']
    result = []
    previous = initial_q
    for s in steps:
        if 'predicted' in s:
            q = s['predicted'][4]
            result.append(dict(stage=prefix, operation='prediction', before=previous, after=q,
                               state=s['state'], z_F_mm=s['z_F_mm'], hit=s.get('hit')))
            previous = q
        if s['measurement_flag'] and 'predicted' in s and 'filtered' in s:
            q = s['filtered'][4]
            result.append(dict(stage=prefix, operation='update', before=s['predicted'][4], after=q,
                               state=s['state'], z_F_mm=s['z_F_mm'], hit=s.get('hit'),
                               innovation=s.get('update_residual_mm'), gain_qop=s.get('update_gain_qop')))
        if 'filtered' in s:
            previous = s['filtered'][4]
    for s in reversed(steps):
        if 'filtered' in s and 'smoothed' in s:
            result.append(dict(stage=prefix, operation='smoothing', before=s['filtered'][4], after=s['smoothed'][4],
                               state=s['state'], z_F_mm=s['z_F_mm'], hit=s.get('hit')))
    return result


def analyze_particle(row, records):
    final = row['final']
    seed_id, candidate_id = final['seed'], final['candidate']
    path = {r['stage']: r for r in records if r['seed'] == seed_id and
            (r['stage'] in ['seed', 'seed_circle'] or r['candidate'] == candidate_id)}
    truth = next(r for r in records if r['stage']=='truth' and r['barcode']==row['barcode'])
    seed = path['seed']
    truth_q = row['truth_charge']/row['truth_p_GeV']
    seed_q = seed['bound'][4]
    result = dict(dsid=row['dsid'], file_index=row['file_index'], source_entry=row['source_entry'],
                  barcode=row['barcode'], truth_p_GeV=row['truth_p_GeV'], truth_charge=row['truth_charge'],
                  final_p_ratio=final['p_ratio'], final_wrong=final['wrong'], seed=seed_id, candidate=candidate_id,
                  seed_wrong=finite_q(seed_q) and seed_q*truth_q<0,
                  final_hit_ids=path['final']['hits'], selected_hit_ids=path.get('selected',{}).get('hits',[]),
                  seed_diagnostics=seed_diagnostics(seed, path.get('seed_circle', {}), truth, records), traces={})
    bounds = envelope(records, row['barcode'], truth['pdg'], final['position_A_mm'][0])
    result['local_class'] = category(final['p_GeV'], bounds)
    result['local_bounds'] = bounds
    all_ops = []
    for stage, initial in [('selected', seed_q), ('refit', path.get('refit_input', seed)['bound'][4])]:
        trace = path.get(stage+'_state_trace')
        if trace is None:
            continue
        steps = trace['steps']
        all_ops.extend(operations(trace, initial, stage))
        totals = dict(states=len(steps), material_X0=sum(s.get('material_full_X0',0.)*s.get('material_path_correction',1.) for s in steps),
                      material_thickness_mm=sum(s.get('material_full_thickness_mm',0.)*s.get('material_path_correction',1.) for s in steps),
                      max_slab_X0=max((s.get('material_full_X0',0.) for s in steps),default=0.),
                      max_gain_qop=max((abs(v) for s in steps for v in s.get('update_gain_qop',[])),default=0.),
                      max_update_error=max((abs(s.get('update_identity_error',0.)) for s in steps),default=0.),
                      max_smoothing_error=max((abs(s.get('smoothing_identity_error',0.)) for s in steps),default=0.),
                      invalid_matrix_checks=sum(bool(s.get('invalid_matrix_check')) for s in steps))
        result['traces'][stage] = totals
        totals['max_update_error_over_truth_qop']=totals['max_update_error']/abs(truth_q)
        totals['max_smoothing_error_over_truth_qop']=totals['max_smoothing_error']/abs(truth_q)
        first = steps[0] if steps else {}
        target = path.get('selected' if stage=='selected' else 'refit')
        if target and pnorm(target):
            qt = target.get('bound',[0.,0.,0.,0.,target.get('charge',1.)/pnorm(target)])[4]
            smoothed = first.get('smoothed',first.get('filtered'))
            if smoothed:
                all_ops.append(dict(stage=stage,operation='target',before=smoothed[4],after=qt,
                                    state=-1,z_F_mm=target['position_A_mm'][0],hit=None))
        if stage=='selected' and 'converted' in path and 'refit_input' in path:
            for a,b in [('selected','converted'),('converted','refit_input')]:
                pa,pb=pnorm(path[a]),pnorm(path[b])
                if pa and pb:
                    qa=path[a].get('charge',math.copysign(1.,path[a].get('bound',[0]*5)[4]))/pa
                    qb=path[b].get('charge',math.copysign(1.,path[b].get('bound',[0]*5)[4]))/pb
                    all_ops.append(dict(stage='interfaces',operation=a+'_to_'+b,before=qa,after=qb,state=-1,
                                        z_F_mm=path[b]['position_A_mm'][0],hit=None))
    crossings = [o for o in all_ops if finite_q(o['before']) and finite_q(o['after']) and
                 ratio(o['before'],truth_q)>=.5 and ratio(o['after'],truth_q)<.5]
    flips = [o for o in all_ops if finite_q(o['before']) and finite_q(o['after']) and o['before']*o['after']<0]
    result['first_new_tail_operation'] = crossings[0] if crossings else None
    result['first_sign_flip_operation'] = flips[0] if flips else None
    valid_ops = [o for o in all_ops if finite_q(o['before']) and finite_q(o['after'])]
    result['largest_qop_change'] = max(valid_ops,key=lambda o:abs(o['after']-o['before'])) if valid_ops else None
    if result['largest_qop_change']:
        result['largest_qop_change']['delta_over_truth_qop'] = abs(result['largest_qop_change']['after']-result['largest_qop_change']['before'])/abs(truth_q)
    seed_r = ratio(seed_q,truth_q)
    lineage = row['final_lineage']
    if final['p_ratio'] < .5:
        first_stage = next((s for s in ['seed','found','converted','refit_input','refit','final']
                            if s in lineage and lineage[s]['p_ratio']<.5), 'invalid_seed_or_target')
        crossing_stage = {'found':'selected','converted':'interfaces','refit_input':'interfaces',
                          'refit':'refit','final':'refit'}.get(first_stage)
        relevant = [o for o in crossings if o['stage']==crossing_stage]
        if seed_r is not None and seed_r < .5:
            result['tail_class'] = 'seed:'+result['seed_diagnostics']['primary_class']
        elif relevant:
            op = relevant[0]
            result['tail_class'] = op['stage']+':'+op['operation']
        else:
            result['tail_class'] = 'unresolved_recorded_path'
        result['first_tail_stage'] = first_stage
        result['tail_entry_operation'] = relevant[0] if relevant else None
    shared = path.get('selected_shared_state_trace')
    if shared:
        result['shared_smoothing_nonclosure'] = any(abs(s.get('smoothing_identity_error',0.))>1e-6*abs(truth_q) for s in shared['steps'])
    result['conversion_relative_p_change'] = abs(pnorm(path['selected'])/pnorm(path['converted'])-1) if 'selected' in path and 'converted' in path and pnorm(path['converted']) else None
    result['refit_input_relative_p_change'] = abs(pnorm(path['converted'])/pnorm(path['refit_input'])-1) if 'converted' in path and 'refit_input' in path and pnorm(path['refit_input']) else None
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--nominal',type=Path,help='Original compact audit, for exact physics-observer pairing')
    parser.add_argument('--limit-tasks',type=int)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    tasks=json.loads((args.input/'tasks.json').read_text())
    if args.limit_tasks:tasks=tasks[:args.limit_tasks]
    nominal_index=None
    if args.nominal:
        archive=(args.nominal/'index.json').exists()
        entries=json.loads((args.nominal/('index.json' if archive else 'tasks.json')).read_text())
        entries=entries['tasks'] if archive else entries
        nominal_index={(t['dsid'],t['file_index'],t['skip']):t for t in entries}
    counts=defaultdict(Counter)
    values=defaultdict(lambda:defaultdict(list))
    seed_counts=defaultdict(Counter)
    tail_classes=defaultdict(Counter)
    candidate_counts=defaultdict(Counter)
    examples=defaultdict(list)
    same=Counter()
    selection=['seed','found','selected','converted','refit_input','refit','final']
    fields=['seed','candidate','momentum_A_GeV','position_A_mm','charge','bound','bound_trk','variance_qop',
            'hits','matches','chi2','ndof','measurements','layers','rmax_mm','refit']
    layer_values=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    local_particles=Path.cwd()/'operation_particles.jsonl.gz'
    with gzip.open(local_particles,'wt') as output:
        for task in tasks:
            folder=args.input/task['shard']
            assert (folder/'exit_code.txt').read_text().strip()=='0',folder
            events=event_groups(read_records(folder/'stages.jsonl.gz'))
            nominal_events=None
            if nominal_index is not None:
                nt=nominal_index[(task['dsid'],task['file_index'],task['skip'])]
                nr=archive_records(args.nominal,nt) if archive else read_records(args.nominal/nt['shard']/'stages.jsonl.gz')
                nominal_events=event_groups(nr)
            for ordinal,records in enumerate(events):
                dsid=task['dsid'];counter=counts[dsid];counter['events']+=1
                if nominal_events is not None:
                    original=next(nominal_events)
                    for stage in selection:
                        left=[{k:r[k] for k in fields if k in r} for r in original if r['stage']==stage]
                        right=[{k:r[k] for k in fields if k in r} for r in records if r['stage']==stage]
                        same[stage+'_events']+=1
                        if left!=right:
                            same[stage+'_different']+=1
                            if len(examples['observer_difference'])<5:
                                examples['observer_difference'].append(dict(dsid=dsid,file_index=task['file_index'],source_entry=task['skip']+ordinal,stage=stage))
                header=records[0];rows,flow=process_event(tuple(header[k] for k in ['run','lumi','event']),records,dsid)
                for key,value in flow.items():
                    if key!='events':counter['cutflow_'+key]+=value
                for row in rows:
                    row.update(file_index=task['file_index'],source_entry=task['skip']+ordinal)
                    if not (row['fiducial'] and row['truth_H18']):continue
                    counter['fiducial_H18']+=1
                    found=[r for r in records if r['stage']=='found' and match(r)==row['barcode'] and pnorm(r)]
                    selected=[r for r in records if r['stage']=='selected' and match(r)==row['barcode'] and pnorm(r)]
                    accurate=lambda r: .8<=pnorm(r)/row['truth_p_GeV']<=1.2 and r.get('charge',1)*row['truth_charge']>0 and r.get('measurements',0)>=12
                    cf=bool(any(accurate(r) for r in found));cs=bool(any(accurate(r) for r in selected))
                    cc=candidate_counts[dsid];cc['truth']+=1;cc['has_accurate_found']+=cf;cc['accurate_found_no_accurate_selected']+=cf and not cs
                    cc['accurate_found_no_matched_selected']+=cf and not selected
                    if not row['final'] or not row['final']['good95']:
                        cc['accurate_found_no_good95_final']+=cf
                        continue
                    result=analyze_particle(row,records)
                    counter['good95']+=1
                    sd=result['seed_diagnostics']
                    sc=seed_counts[dsid];sc['final_ancestry_seeds']+=1
                    sc['invalid_seed_momentum']+=sd['p_ratio'] is None
                    sc['single_station_default']+=bool(sd.get('single_station_default'))
                    sc['seed_majority_other']+=sd['seed_majority_other']
                    sc['segment_majority_other']+=sd['segment_majority_other']
                    sc['seed_ambiguous_70pct']+=sd['seed_ambiguous_70pct']
                    if sd['p_ratio'] is not None and sd['p_ratio']<.5:
                        sc['production_seed_tail']+=1;sc['tail:'+sd['primary_class']]+=1
                    for key in ['radius_recompute_relative_error','circle_condition','normalized_determinant','centers_only_p_ratio','truth_directions_p_ratio']:
                        if sd.get(key) is not None and math.isfinite(sd[key]):values[dsid]['seed_'+key].append(sd[key])
                    values[dsid]['segment_abs_y_slope_error'].extend(abs(v) for v in sd['segment_y_slope_errors'])
                    values[dsid]['segment_abs_y_center_error_mm'].extend(abs(v) for v in sd['segment_y_center_errors_mm'])
                    if result['final_p_ratio']<.5:
                        counter['production_tail']+=1
                        counter['production_tail_local:'+result['local_class']]+=1
                        tail_classes[dsid][result['tail_class']]+=1
                        counter['first_tail_stage:'+result['first_tail_stage']]+=1
                        label=result['tail_class']
                        if len(examples[dsid+':'+label])<3:examples[dsid+':'+label].append(result)
                    if result['final_wrong']:
                        counter['wrong']+=1
                        flip=result['first_sign_flip_operation']
                        counter['wrong_first_flip:'+('seed' if result['seed_wrong'] else (flip['stage']+':'+flip['operation'] if flip else 'unrecorded'))]+=1
                    counter['shared_smoothing_nonclosure']+=result.get('shared_smoothing_nonclosure',False)
                    for key in ['conversion_relative_p_change','refit_input_relative_p_change']:
                        if result[key] is not None:values[dsid][key].append(result[key])
                    for stage,metrics in result['traces'].items():
                        for key,value in metrics.items():values[dsid][stage+'_'+key].append(value)
                        for name in ['update','smoothing']:
                            key='max_'+name+'_error_over_truth_qop'
                            counter[stage+'_'+name+'_nonclosure_over_1e-6']+=metrics[key]>1e-6
                            label=dsid+':'+stage+':max_'+name+'_error'
                            previous=examples[label]
                            if not previous or metrics[key]>previous[0]['value']:
                                examples[label]=[dict(dsid=dsid,file_index=row['file_index'],source_entry=row['source_entry'],
                                                      barcode=row['barcode'],value=metrics[key])]
                    for trace in (r for r in records if r['stage'] in ['selected_state_trace','refit_state_trace'] and r['seed']==result['seed'] and r['candidate']==result['candidate']):
                        for step in trace['steps']:
                            if not step['measurement_flag']:continue
                            z=round(step['z_F_mm'],2)
                            for component in ['predicted','filtered','smoothed']:
                                if component in step:
                                    p_ratio=ratio(step[component][4],row['truth_charge']/row['truth_p_GeV'])
                                    if p_ratio is not None:layer_values[dsid][trace['trace_stage']][(z,component)].append(p_ratio)
                    output.write(json.dumps(result,separators=(',',':'),allow_nan=False)+'\n')
            if nominal_events is not None:
                assert next(nominal_events,None) is None,task['shard']
            print('Completed',task['dsid'],task['file_index'],task['skip'],flush=True)
    publish_gzip(local_particles,args.output/'particles.jsonl.gz')
    def describe(x):
        return dict(n=len(x),quantiles=np.quantile(x,[0,.16,.5,.84,1]).tolist()) if x else dict(n=0)
    result=dict(samples={d:dict(counts=dict(c),seed_counts=dict(seed_counts[d]),tail_classes=dict(tail_classes[d]),
                                candidate_counts=dict(candidate_counts[d]),metrics={k:describe(v) for k,v in values[d].items()})
                         for d,c in counts.items()},observer_pair=dict(same))
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (args.output/'examples.json').write_text(json.dumps(dict(examples),indent=2)+'\n')
    for dsid,stages in layer_values.items():
        fig,axes=plt.subplots(1,2,figsize=(12,4),sharey=True)
        for ax,(stage,data) in zip(axes,sorted(stages.items())):
            for component in ['predicted','filtered','smoothed']:
                zs=sorted(z for z,c in data if c==component)
                medians=[np.median(data[(z,component)]) for z in zs]
                ax.plot(zs,medians,marker='.',label=component)
            ax.axhline(1.,color='grey',lw=1);ax.axhline(.5,color='grey',lw=1,ls='--')
            ax.set(xlabel='Physical z (FASER frame) [mm]',title=stage,ylim=(0,2));ax.legend()
        axes[0].set_ylabel('Median p / production truth p')
        fig.suptitle('DSID '+dsid+': final Good95 ancestry');fig.tight_layout()
        fig.savefig(args.output/(dsid+'_layers.png'),dpi=160);plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
