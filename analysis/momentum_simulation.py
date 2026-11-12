#!/usr/bin/env python3
"""Analyze deterministic same-plane truth injection on a paired G4 mother sample."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

if __package__:
    from .pilot import event_groups, read_records, process_event, wilson
    from .artifact_io import publish_gzip
else:
    from pilot import event_groups, read_records, process_event, wilson
    from artifact_io import publish_gzip


def sensor_bounds(path):
    groups = defaultdict(list)
    with gzip.open(path,'rt') as stream:
        reader=csv.DictReader(stream)
        assert {'event','source_entry','g4_event'}.issubset(reader.fieldnames),'Sensor truth requires Athena event identity (patch 0005)'
        for row in reader:
            if int(row['parent']) != 0 or int(row['track']) != 1:
                continue
            key=tuple(int(row[k]) for k in ['event','station','layer','row','module','sensor'])
            groups[key].append(row)
    result=defaultdict(list)
    for key,rows in groups.items():
        momenta=[math.sqrt(sum(float(row[prefix+component+'_MeV'])**2 for component in ['px','py','pz']))*.001
                 for row in rows for prefix in ['', 'post_']]
        result[key[0]].append(dict(z_entry_mm=min(float(row['z_mm']) for row in rows),
                                  lower_GeV=min(momenta),upper_GeV=max(momenta),sensor=key[1:]))
    return result


def local_bound(sensors,z):
    if not sensors:return None
    sensor=min(sensors,key=lambda r:abs(r['z_entry_mm']+.1425-z))
    return sensor if abs(sensor['z_entry_mm']+.1425-z)<.6 else None


def metrics(rows, total):
    good=[r for r in rows if r['final'] and r['final']['good95'] and r['fiducial'] and r['truth_H18']]
    ratio=np.array([r['final']['p_ratio'] for r in good])
    tails={str(t):int(np.sum(ratio<t)) for t in [.8,.5,.2]}
    return dict(mother_events=total,truth=len(rows),fiducial_H18=sum(r['fiducial'] and r['truth_H18'] for r in rows),
                matched_final=sum(bool(r['final']) for r in rows),good95=len(good),good95_efficiency=len(good)/total if total else None,
                tails=tails,tail_half_wilson95=wilson(tails['0.5'],len(good)),
                wrong=sum(r['final']['wrong'] for r in good),
                p_ratio_quantiles=np.quantile(ratio,[.16,.5,.84]).tolist() if len(good) else [],
                actually_injected_final=sum(r.get('actually_injected',False) for r in good),
                sensor_local_tail_classes=dict(Counter(r.get('sensor_local_class','unavailable') for r in good)))


def t1_selection(rows,injections):
    seeds=defaultdict(set)
    for item in injections.values():
        if item['eligible'] and item.get('good_seed_T1',False):
            seeds[(item['event_ordinal'],item['barcode'])].add(item['seed'])
    mother=[];path_rows=[];counts=Counter()
    for row in rows:
        eligible=seeds.get((row['source_entry'],row['barcode']))
        if not eligible:continue
        mother.append(row)
        final=row['final'];same=bool(final and final['seed'] in eligible)
        counts['final_from_good_seed']+=same
        counts['final_from_other_seed']+=bool(final and not same)
        counts['no_matched_final']+=not final
        path_rows.append(row if same else dict(row,final=None))
    result=metrics(path_rows,len(mother))
    result.update(definition='Mother truth has at least one good T0 seed before finding; final must descend from one of those seeds',
                  lineage_counts=dict(counts),any_seed_final=metrics(mother,len(mother)))
    return result


def paired(a,b):
    a={(r['source_entry'],r['barcode']):r for r in a};b={(r['source_entry'],r['barcode']):r for r in b}
    assert a.keys()==b.keys()
    matrix=[];counts=Counter()
    for key,x in a.items():
        y=b[key]
        assert x['truth_p_GeV']==y['truth_p_GeV']
        if not (x['fiducial'] and x['truth_H18']):continue
        gx=bool(x['final'] and x['final']['good95']);gy=bool(y['final'] and y['final']['good95'])
        tx=bool(gx and x['final']['p_ratio']<.5);ty=bool(gy and y['final']['p_ratio']<.5)
        counts[f'good_{int(gx)}{int(gy)}']+=1
        if gx and gy:counts[f'tail_{int(tx)}{int(ty)}']+=1
        matrix.append([1,gx,gy,tx,ty])
    if not matrix:return dict(counts=dict(counts))
    matrix=np.asarray(matrix,dtype=float)
    def rates(x):
        with np.errstate(divide='ignore',invalid='ignore'):
            return np.stack([(x[...,2]-x[...,1])/x[...,0],x[...,4]/x[...,2]-x[...,3]/x[...,1]],axis=-1)
    rng=np.random.default_rng(80501)
    bootstrap=matrix[rng.integers(0,len(matrix),size=(1000,len(matrix)))].sum(axis=1)
    intervals=np.nanquantile(rates(bootstrap),[.025,.975],axis=0)
    return dict(counts=dict(counts),differences={label:dict(value=float(v),event_bootstrap95=intervals[:,i].tolist())
           for i,(label,v) in enumerate(zip(['good95_efficiency','tail_half_fraction'],rates(matrix.sum(axis=0))))})


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--tasks',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    summary={};all_rows=[];checks=Counter();layer_values=defaultdict(lambda:defaultdict(list))
    for task in json.loads(args.tasks.read_text()):
        folder=Path(task['output'])
        assert (folder/'exit_code.txt').read_text().strip()=='0',folder
        with gzip.open(folder/'injections.jsonl.gz','rt') as stream:injections=[json.loads(line) for line in stream]
        injections={(r['event_ordinal'],r['seed']):r for r in injections}
        sensors=sensor_bounds(Path(task.get('sensor_truth',folder/'sensor_truth.csv.gz')))
        assert len(sensors)==task['events'],('Incomplete sensor event identity',task['shard'],len(sensors))
        mode_rows={};mode_cutflow={};seed_cov={}
        for mode in task['modes']:
            rows=[];counts=Counter()
            for ordinal,records in enumerate(event_groups(read_records(folder/mode/'stages.jsonl.gz'))):
                header=records[0]
                event_rows,flow=process_event(tuple(header[k] for k in ['run','lumi','event']),records,task['dsid'])
                counts.update(flow)
                for seed in [r for r in records if r['stage']=='seed']:
                    key=(ordinal,seed['seed']);inj=injections.get(key)
                    if mode=='T0':seed_cov[key]=seed['cov']
                    elif mode in ['T2','T3','T4','T3_current','T4_current'] and inj and inj['eligible']:
                        assert seed['cov']==seed_cov[key],('Covariance changed under mean injection',task['shard'],mode,key)
                        expected=list(inj['truth_mean'])
                        if mode=='T2':
                            expected=list(inj['original_mean']);expected[4]=inj['truth_mean'][4]
                        assert np.allclose(seed['bound'],expected,rtol=0,atol=1e-13),(mode,key,seed['bound'],expected)
                        checks['verified_actual_seed_mean_and_covariance']+=1
                for row in event_rows:
                    row.update(source_entry=ordinal,file_index=task['file_index'],mode=mode,grid=task['shard'])
                    final=row['final'];inj=injections.get((ordinal,final['seed'])) if final else None
                    row['actually_injected']=bool(inj and inj['eligible'] and mode.startswith(('T2','T3','T4')))
                    row['good_seed_T1']=bool(inj and inj.get('good_seed_T1',False))
                    if final:
                        bound=local_bound(sensors[header['event']],final['position_A_mm'][0])
                        row['sensor_bound']=bound
                        if bound:
                            row['sensor_local_class']='guaranteed_local_tail' if final['p_GeV']<.5*bound['lower_GeV'] else 'excluded_local_tail' if final['p_GeV']>=.5*bound['upper_GeV'] else 'interval_unresolved'
                        for trace in [r for r in records if r['stage'] in ['selected_state_trace','refit_state_trace'] and r['seed']==final['seed'] and r['candidate']==final['candidate']]:
                            for step in trace['steps']:
                                if not step['measurement_flag']:continue
                                bound=local_bound(sensors[header['event']],step['z_F_mm'])
                                if not bound:checks['state_sensor_bound_missing']+=1;continue
                                for component in ['predicted','filtered','smoothed']:
                                    if component not in step:continue
                                    q=step[component][4]
                                    if q is None or not math.isfinite(q) or q==0:checks['invalid_state_qop']+=1;continue
                                    # The pre/post interval is retained; the plot uses its upper endpoint.
                                    layer_values[task['shard']+':'+mode+':'+trace['trace_stage']][(round(step['z_F_mm'],2),component)].append(1/abs(q)/bound['upper_GeV'])
                    rows.append(row)
                checks['processed_mode_events']+=1
            assert counts['events']==task['events']
            mode_rows[mode]=rows;mode_cutflow[mode]=dict(counts);all_rows.extend(rows)
        results={mode:metrics(rows,task['events']) for mode,rows in mode_rows.items()}
        results['T1']=t1_selection(mode_rows['T0'],injections)
        results['T1_current']=t1_selection(mode_rows['T0_current'],injections)
        comparisons={}
        for a,b in [('T0','T2'),('T2','T3'),('T3','T4'),('T0_current','T0'),('T3_current','T3'),
                    ('T4_current','T4'),('T4_current','T4_no_scatt_current'),('T4_current','T4_no_eloss_current'),
                    ('T4','T4_electron'),('T4','T4_seedC01'),('T4','T4_interfaces')]:
            if a in mode_rows and b in mode_rows:comparisons[a+'->'+b]=paired(mode_rows[a],mode_rows[b])
        summary[task['shard']]=dict(pdg=task['pdg'],p_GeV=task['p_GeV'],events=task['events'],sensor_truth_events=len(sensors),
            injection_seeds=len(injections),eligible_seeds=sum(i['eligible'] for i in injections.values()),
            injection_plane_distance_max_mm=max((abs(i['plane_distance_mm']) for i in injections.values() if i['eligible']),default=None),
            results=results,cutflow=mode_cutflow,comparisons=comparisons)
        print('Completed',task['shard'],flush=True)
    (args.output/'summary.json').write_text(json.dumps(dict(grids=summary,checks=dict(checks)),indent=2)+'\n')
    local_particles=Path.cwd()/'simulation_particles.jsonl.gz'
    with gzip.open(local_particles,'wt') as stream:
        for row in all_rows:stream.write(json.dumps(row,separators=(',',':'),allow_nan=False)+'\n')
    publish_gzip(local_particles,args.output/'particles.jsonl.gz')
    for grid,item in summary.items():
        labels=['T0_current','T1_current','T0','T1','T2','T3','T4','T4_current','T4_no_scatt_current','T4_no_eloss_current']
        labels=[m for m in labels if m in item['results']]
        fig,axes=plt.subplots(1,2,figsize=(13,4))
        for ax,quantity in zip(axes,['tail','efficiency']):
            ys=[];lo=[];hi=[]
            for mode in labels:
                r=item['results'][mode]
                n=r['good95'] if quantity=='tail' else r['mother_events']
                k=r['tails']['0.5'] if quantity=='tail' else r['good95']
                v=k/n if n else float('nan');interval=wilson(k,n) if n else [0.,0.]
                ys.append(v);lo.append(v-interval[0] if n else 0.);hi.append(interval[1]-v if n else 0.)
            ax.errorbar(range(len(labels)),ys,yerr=[lo,hi],fmt='o')
            ax.set_xticks(range(len(labels)));ax.set_xticklabels(labels,rotation=55,ha='right')
            ax.set_ylabel('Tail fraction (p ratio < 0.5)' if quantity=='tail' else 'Good95 / mother sample')
            ax.grid(axis='y',alpha=.3)
        fig.suptitle(f"PDG {item['pdg']}, p = {item['p_GeV']} GeV; {item['events']} G4 events")
        fig.tight_layout();fig.savefig(args.output/(grid+'_injections.png'),dpi=160);plt.close(fig)
    # Full numerical per-layer summaries remain exportable without repeating matrices.
    layers={name:[dict(z_F_mm=z,component=c,n=len(v),p_over_sensor_upper_quantiles=np.quantile(v,[.16,.5,.84]).tolist())
                  for (z,c),v in sorted(data.items())] for name,data in layer_values.items()}
    (args.output/'sensor_layers.json').write_text(json.dumps(layers,indent=2)+'\n')
    print(json.dumps(dict(checks=dict(checks),grids=len(summary))))


if __name__=='__main__':main()
