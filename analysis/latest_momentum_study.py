#!/usr/bin/env python3
"""Attribute native Calypso q/p changes, retaining the original momentum plot."""
import argparse
from collections import Counter, defaultdict
import json
import os
from pathlib import Path

import numpy as np


def describe(values):
    values = np.asarray(values)
    finite = values[np.isfinite(values)]
    if not len(finite):
        return dict(raw=len(values), finite=0)
    return dict(raw=len(values), finite=len(finite),
        mean=float(np.mean(finite)), std=float(np.std(finite)),
        quantiles=dict(zip(['p05','p16','p50','p84','p95'],np.quantile(finite,[.05,.16,.5,.84,.95]).tolist())),
        below_minus01=float(np.mean(finite<-.1)), above_plus01=float(np.mean(finite>.1)),
        below_minus05=float(np.mean(finite<-.5)))


def p_residual(q, p):
    with np.errstate(divide='ignore',invalid='ignore'):
        return 1/(np.abs(q)*p)-1


def native_formula_residual(q, truth_q):
    # ROOT TTreeFormula's division-by-zero protection is part of the reference.
    def divide(a,b):
        a,b=np.broadcast_arrays(a,b)
        with np.errstate(divide='ignore',invalid='ignore'):
            return np.divide(a,b,out=np.zeros_like(a,dtype=float),where=b!=0)
    p=divide(1.,np.abs(q));truth=divide(1.,np.abs(truth_q))
    return divide(p-truth,truth)


def describe_p(values):
    result=describe(values)
    edges=np.linspace(-1.2,1.2,101)
    counts,_=np.histogram(np.asarray(values)[np.isfinite(values)],edges)
    if np.any(counts):
        index=int(np.argmax(counts))
        result['modal_bin']=[float(edges[index]),float(edges[index+1])]
    return result


def subset_files(data, files):
    """Restrict all tables together so figures cannot mix input coverage."""
    files=set(files)
    state_mask=np.isin(data['values'][:,0],list(files))
    result={**data,'values':data['values'][state_mask],
        'hits':data['hits'][state_mask],'geometry':data['geometry'][state_mask],
        'tracks':data['tracks'][np.isin(data['tracks'][:,0],list(files))],
        'seeds':data['seeds'][np.isin(data['seeds'][:,0],list(files))],
        'coverage':[x for x in data['coverage'] if x['file_index'] in files]}
    assert {x['file_index'] for x in result['coverage']}==files
    return result


def event_interval(success, paired, source_events, bounded=True):
    keys=paired[:,:2].astype(int)
    _,inverse=np.unique(keys,axis=0,return_inverse=True)
    count=np.bincount(inverse)
    passed=np.bincount(inverse,weights=np.asarray(success,dtype=float))
    ratio=float(np.sum(passed)/np.sum(count))
    variance=source_events/(source_events-1)*np.sum((passed-ratio*count)**2)/np.sum(count)**2
    half=1.96*np.sqrt(variance)
    interval=[ratio-half,ratio+half]
    if bounded:interval=[max(0.,interval[0]),min(1.,interval[1])]
    return dict(fraction=ratio,ci95=interval,
                unit='source-event cluster delta interval')


def load_group(root, name, files=None):
    paths = sorted((root/'runs'/name).glob('file*'))
    if files is not None:
        paths=[p for p in paths if int(p.name[4:]) in files]
    if name=='baseline' and (root/'baseline_replacements.json').exists():
        replacements=json.loads((root/'baseline_replacements.json').read_text())
        paths=[Path(replacements.get(p.name,str(p))) for p in paths]
    values, tracks, seeds, identities, hits = [], [], [], [], []
    coverage = []
    fields = track_fields = None
    for path in paths:
        task = json.loads((path/'source_task.json').read_text())
        assert (path/'exit_code.txt').read_text().strip()=='0', path
        summary=json.loads((path/'summary.json').read_text())
        assert summary['source_events']==task['events']
        run=json.loads((path/'run.json').read_text())
        assert run['success'] and run['isolate_smoothing']
        assert run['args']['material']==str(root/'material-latest-rowmajor.json')
        assert run['calypso_dir'].endswith('/geoRotation_calypso_latest/run')
        assert any('/software_v4/lib/libFaserActsKalmanFilter.so' in p for p in run['loaded_acts_libraries'])
        data = np.load(path/'states.npz')
        writer = np.load(path/'writer.npz')
        assert fields is None or list(fields)==list(data['fields'])
        fields, track_fields = data['fields'], data['track_fields']
        # First column makes input-file identity independent of production event IDs.
        values.append(np.column_stack([np.full(len(data['values']), task['file_index']),data['values']]))
        tracks.append(np.column_stack([np.full(len(data['tracks']),task['file_index']),data['tracks']]))
        seeds.append(np.column_stack([np.full(len(writer['seeds']),task['file_index']),writer['seeds']]))
        identities.append(data['geometry_ids']); hits.append(data['hit_ids'])
        coverage.append(dict(file=task['input'],file_index=task['file_index'],skip=task['skip'],events=task['events']))
    assert len({c['file'] for c in coverage})==len(coverage)
    return dict(values=np.concatenate(values),fields=['file',*fields], tracks=np.concatenate(tracks),
        track_fields=['file',*track_fields],seeds=np.concatenate(seeds),
        geometry=np.concatenate(identities),hits=np.concatenate(hits),paths=paths,coverage=coverage)


def transitions(data, require_good=True, same_truth=False):
    v = data['values']; f = {n:i for i,n in enumerate(data['fields'])}
    keep = (v[:,f['github']]>0)&(v[:,f['purity']]>.7)
    if require_good:
        keep &= v[:,f['good_seed']]>0
    if same_truth:
        keep &= v[:,f['same_seed_truth']]>0
    indices = np.flatnonzero(keep)
    # Input is already in forward path order for each track.
    boundaries = np.r_[0,1+np.flatnonzero(np.any(np.diff(v[indices,:3],axis=0)!=0,axis=1)),len(indices)]
    output, paths = [], []
    for start,end in zip(boundaries[:-1],boundaries[1:]):
        ii=indices[start:end]; a=v[ii]
        measurement=a[:,f['measurement']]>0
        s1=np.flatnonzero(measurement&(a[:,f['station']]==1))
        s2=np.flatnonzero(measurement&(a[:,f['station']]==2))
        if not len(s1) or not len(s2):
            continue
        first2,last1=int(s2[0]),int(s1[-1])
        row1,row2=a[last1],a[first2]
        p=row2[f['p_production_GeV']]
        tq=row2[f['truth_charge']]/p
        seed=row2[f['seed_qop']]
        s1f=row1[f['q_filtered']]; s2p=row2[f['q_predicted']]
        accumulated=a[:first2+1]
        trans=float(np.sum(accumulated[:,f['dq_transport']]))
        updates=float(np.sum(accumulated[:-1,f['dq_update']]))
        closure=(s2p-seed)-(trans+updates)
        assert abs(closure)<1e-10
        output.append([*row2[:3],row2[f['seed']],p,tq,seed,s1f,s2p,
            row2[f['q_smoothed']],trans,updates,closure,
            row2[f['sigma_over_truth']],row2[f['same_seed_truth']],row2[f['track_barcode']],
            row2[f['q_filtered']],a[s2[1],f['q_predicted']] if len(s2)>1 else np.nan,
            row2[f['measurement_residual_mm']],row2[f['gain_qop_per_mm']],
            row2[f['layer']],a[s2[1],f['layer']] if len(s2)>1 else np.nan,
            row2[f['innovation_mm2']],row2[f['corr_qop_phi']],
            row2[f['var_predicted']],row2[f['var_filtered']]])
        paths.append((tuple(v[ii[measurement],f['station']].astype(int)),
                      tuple(int(x) for x in data['hits'][ii[measurement]])))
    columns=['file','entry','candidate','seed','p','truth_qop','seed_qop','s1_filtered',
        's2_predicted','s2_smoothed','sum_transport','sum_updates','closure',
        'sigma_over_truth','same_truth','track_barcode','s2_first_filtered','s2_second_predicted',
        's2_first_residual_mm','s2_first_gain','s2_first_layer','s2_second_layer',
        's2_first_innovation_mm2','s2_first_corr_qop_phi','s2_first_var_predicted','s2_first_var_filtered']
    return np.asarray(output).reshape(-1,len(columns)),columns,paths


def summarize(data):
    v=data['values'];f={n:i for i,n in enumerate(data['fields'])}
    t=data['tracks'];tf={n:i for i,n in enumerate(data['track_fields'])}
    selected=(t[:,tf['github']]>0)&(t[:,tf['purity']]>.7)&(t[:,tf['good_seed']]>0)
    result=dict(coverage=data['coverage'],source_events=sum(x['events'] for x in data['coverage']),
        all_tracks=len(t),good_seed_tracks=int(selected.sum()),
        seed_final_majority_different=int(np.sum(selected&(t[:,tf['seed_barcode']]!=t[:,tf['track_barcode']]))))
    result['layers']={}
    for same in [False,True]:
        base=(v[:,f['github']]>0)&(v[:,f['purity']]>.7)&(v[:,f['good_seed']]>0)&(v[:,f['measurement']]>0)
        if same:base&=v[:,f['same_seed_truth']]>0
        bylayer={}
        for station in [1,2,3]:
            for layer in [0,1,2]:
                rows=v[base&(v[:,f['station']]==station)&(v[:,f['layer']]==layer)]
                bylayer[f'S{station}L{layer}']={}
                for component in ['predicted','filtered','smoothed']:
                    # Native writer truth availability, values in GeV^-1 here.
                    qtruth=rows[:,f['native_truth_qop_GeV_inv']]
                    residual=native_formula_residual(rows[:,f['q_'+component]],qtruth)
                    bylayer[f'S{station}L{layer}'][component]=describe_p(residual)
        result['layers']['same_truth' if same else 'all']=bylayer
    paired,cols,_=transitions(data)
    c={n:i for i,n in enumerate(cols)}
    p=paired[:,c['p']]
    res={n:p_residual(paired[:,c[n]],p) for n in ['seed_qop','s1_filtered','s2_predicted',
        's2_first_filtered','s2_second_predicted','s2_smoothed']}
    result['transition']=dict(count=len(paired),components={n:describe_p(x) for n,x in res.items()},
        relative_p_change_s1f_to_s2p=describe(np.abs(paired[:,c['s1_filtered']]/paired[:,c['s2_predicted']])-1),
        signed_normalized_transport=describe(paired[:,c['sum_transport']]/paired[:,c['truth_qop']]),
        signed_normalized_updates=describe(paired[:,c['sum_updates']]/paired[:,c['truth_qop']]),
        curvature_error=describe(paired[:,c['s2_predicted']]/paired[:,c['truth_qop']]-1),
        sigma_over_truth=describe(paired[:,c['sigma_over_truth']]),
        s2_left=int(np.sum(res['s2_predicted']<-.1)),
        s2_left_already_s1=int(np.sum((res['s2_predicted']<-.1)&(res['s1_filtered']<-.1))))
    result['transition']['s2_two_sensors']=dict(
        count=len(paired),both_layer0=int(np.sum((paired[:,c['s2_first_layer']]==0)&(paired[:,c['s2_second_layer']]==0))),
        first_update_over_truth=describe((paired[:,c['s2_first_filtered']]-paired[:,c['s2_predicted']])/paired[:,c['truth_qop']]),
        first_residual_mm=describe(paired[:,c['s2_first_residual_mm']]),
        first_gain=describe(paired[:,c['s2_first_gain']]),
        first_innovation_sigma_mm=describe(np.sqrt(paired[:,c['s2_first_innovation_mm2']])),
        first_residual_pull=describe(paired[:,c['s2_first_residual_mm']]/np.sqrt(paired[:,c['s2_first_innovation_mm2']])),
        first_corr_qop_phi=describe(paired[:,c['s2_first_corr_qop_phi']]),
        first_filtered_to_second_predicted=describe(np.abs(paired[:,c['s2_first_filtered']]/paired[:,c['s2_second_predicted']])-1))
    result['transition']['left_fraction_interval']=event_interval(res['s2_predicted']<-.1,paired,result['source_events'])
    result['transition']['charge_wrong_fraction']={name:float(np.mean(paired[:,c[name]]*paired[:,c['truth_qop']]<0))
        for name in res}
    seed_phi={(int(s[0]),int(s[1]),int(s[2])):abs(s[6]-s[12]) for s in data['seeds']}
    phi_good=np.asarray([seed_phi[(int(row[0]),int(row[1]),int(row[c['seed']]))]<1e-4 for row in paired])
    result['transition_subsets']={}
    for name,mask in [('same_truth',paired[:,c['same_truth']]>0),
                      ('same_truth_and_both_S2_layer0',(paired[:,c['same_truth']]>0)&
                       (paired[:,c['s2_first_layer']]==0)&(paired[:,c['s2_second_layer']]==0)),
                      ('same_truth_and_phi_1e4',(paired[:,c['same_truth']]>0)&phi_good)]:
        result['transition_subsets'][name]=dict(count=int(mask.sum()),
            components={n:describe_p(x[mask]) for n,x in res.items()},
            curvature={n:describe(paired[mask,c[n]]/paired[mask,c['truth_qop']]-1)
                       for n in res},
            charge_wrong={n:float(np.mean(paired[mask,c[n]]*paired[mask,c['truth_qop']]<0))
                          for n in res} if mask.any() else {})
    ungated,ungated_cols,_=transitions(data,require_good=False,same_truth=True)
    uc={n:i for i,n in enumerate(ungated_cols)}
    result['all_seeds_same_truth_transition']=dict(count=len(ungated),
        components={name:describe_p(p_residual(ungated[:,uc[name]],ungated[:,uc['p']]))
                    for name in res})
    result['seed_gate_ablation']=seed_gate_ablation(data,ungated,ungated_cols)
    seeds=data['seeds']
    result['seeds']=dict(count=len(seeds),
        momentum_residual=describe(np.abs(seeds[:,14]/seeds[:,8])-1),
        phi_difference_mrad=describe((seeds[:,6]-seeds[:,12])*1e3),
        theta_difference_mrad=describe((seeds[:,7]-seeds[:,13])*1e3))
    seed_by_key={tuple(s[:3].astype(int)):s for s in seeds}
    selected_seeds=np.asarray([seed_by_key[(int(row[0]),int(row[tf['entry']]),int(row[tf['seed']]))]
                              for row in t[selected]])
    result['seeds']['selected_track_ancestors']=dict(track_weighted_count=len(selected_seeds),
        momentum_residual=describe_p(np.abs(selected_seeds[:,14]/selected_seeds[:,8])-1),
        phi_difference_mrad=describe((selected_seeds[:,6]-selected_seeds[:,12])*1e3),
        theta_difference_mrad=describe((selected_seeds[:,7]-selected_seeds[:,13])*1e3))
    native_reference=v[:,f['native_truth_qop_GeV_inv']]*v[:,f['p_production_GeV']]/v[:,f['truth_charge']]-1
    result['native_truth_vs_production_curvature']=describe(native_reference)
    error=v[:,f['update_identity_error']];error=error[np.isfinite(error)]
    result['max_update_identity_error']=float(np.max(np.abs(error)))
    cov_error=v[:,f['variance_update_identity_error']]
    cov_error=cov_error[np.isfinite(cov_error)]
    result['max_variance_update_identity_error']=float(np.max(np.abs(cov_error)))
    result['covariance_health']={name:dict(
        below_minus_1e8=int(np.sum(v[:,f[name]] < -1e-8)),
        nonfinite=int(np.sum(~np.isfinite(v[:,f[name]]))),minimum=float(np.nanmin(v[:,f[name]])))
        for name in ['predicted_corr_min_eigenvalue','filtered_corr_min_eigenvalue']}
    result['updates_by_layer']={}
    mask=(v[:,f['github']]>0)&(v[:,f['purity']]>.7)&(v[:,f['good_seed']]>0)&(v[:,f['measurement']]>0)
    result['selected_update_checks']=dict(
        max_mean_identity_error_GeV_inv=float(np.nanmax(np.abs(v[mask,f['update_identity_error']]))),
        max_covariance_relative_identity_error=float(np.nanmax(np.abs(v[mask,f['variance_update_identity_error']])/
            np.maximum(np.abs(v[mask,f['var_predicted']]),np.abs(v[mask,f['var_filtered']])))))
    for station in [1,2,3]:
        for layer in [0,1,2]:
            a=v[mask&(v[:,f['station']]==station)&(v[:,f['layer']]==layer)]
            result['updates_by_layer'][f'S{station}L{layer}']={name:describe(a[:,f[name]])
                for name in ['measurement_residual_mm','gain_qop_per_mm','corr_qop_phi','var_predicted','var_filtered']}
            result['updates_by_layer'][f'S{station}L{layer}']['dq_over_truth']=describe(
                a[:,f['dq_update']]*a[:,f['p_production_GeV']]/a[:,f['truth_charge']])
    result['first_tail_operations']=first_tail_operations(data)
    return result


def seed_gate_values(data,rows,columns):
    c={name:i for i,name in enumerate(columns)}
    lookup={tuple(s[:3].astype(int)):s for s in data['seeds']}
    s=np.asarray([lookup[(int(row[0]),int(row[1]),int(row[c['seed']]))] for row in rows])
    with np.errstate(divide='ignore',invalid='ignore'):
        residual=np.abs(s[:,14]/s[:,8])-1
    momentum=np.isfinite(residual)&(np.abs(residual)<.1)
    charge=s[:,8]*s[:,14]>0
    theta=np.abs(s[:,7]-s[:,13])<1e-4
    cuts=dict(all=np.ones(len(rows),bool),p_only=momentum,charge_only=charge,theta_only=theta,
        p_and_charge=momentum&charge,charge_and_theta=charge&theta,
        original=momentum&charge&theta)
    return residual,cuts


def seed_gate_ablation(data,rows,columns):
    """Selections of one frozen native reconstruction, with no rerun of tracking."""
    c={name:i for i,name in enumerate(columns)}
    residual,cuts=seed_gate_values(data,rows,columns)
    output={}
    for name,mask in cuts.items():
        output[name]=dict(tracks=int(mask.sum()),
            seed_residual=describe_p(residual[mask]),
            s2_first_predicted=describe_p(p_residual(rows[mask,c['s2_predicted']],rows[mask,c['p']])),
            s2_smoothed=describe_p(p_residual(rows[mask,c['s2_smoothed']],rows[mask,c['p']])))
    return output


def model_gate_validation(data, closure):
    """Validate a seed-model-only scale on files excluded from its determination."""
    rows,columns,_=transitions(data,require_good=False,same_truth=True)
    c={name:i for i,name in enumerate(columns)}
    residual,cuts=seed_gate_values(data,rows,columns)
    factor=1+closure['residual']['quantiles']['p50']
    evaluation=rows[:,c['file']]!=0
    model_gate=(np.abs((1+residual)/factor-1)<.1)&cuts['charge_and_theta']
    masks=dict(no_seed_gate=evaluation,original=evaluation&cuts['original'],
        native_model_calibrated=evaluation&model_gate)
    output=dict(calibration_file=0,calibration_source_events=closure['source_events'],
        evaluation_source_events=sum(x['events'] for x in data['coverage'] if x['file_index']!=0),
        factor=factor,definition='Analysis selection only; frozen reconstruction unchanged; file 0 excluded from evaluation')
    output['selections']={name:dict(tracks=int(mask.sum()),
        seed_residual=describe_p(residual[mask]),
        s2_smoothed=describe_p(p_residual(rows[mask,c['s2_smoothed']],rows[mask,c['p']]))
        ) for name,mask in masks.items()}
    return output,rows,c,masks


def first_tail_operations(data):
    """First crossing labels describe a trajectory's history, not causation."""
    v=data['values'];f={n:i for i,n in enumerate(data['fields'])}
    keep=(v[:,f['github']]>0)&(v[:,f['purity']]>.7)&(v[:,f['good_seed']]>0)&(v[:,f['same_seed_truth']]>0)
    ii=np.flatnonzero(keep)
    boundaries=np.r_[0,1+np.flatnonzero(np.any(np.diff(v[ii,:3],axis=0)!=0,axis=1)),len(ii)]
    counters=defaultdict(Counter)
    for begin,end in zip(boundaries[:-1],boundaries[1:]):
        a=v[ii[begin:end]]
        s2=np.flatnonzero((a[:,f['station']]==2)&(a[:,f['measurement']]>0))
        if not len(s2):continue
        first2=s2[0];p=a[first2,f['p_production_GeV']]
        for threshold,name in [(-.1,'minus01'),(-.5,'minus05')]:
            if p_residual(a[first2,f['q_predicted']],p)>=threshold:continue
            counter=counters[name];counter['s2_left']+=1
            history=[('seed',p_residual(a[0,f['seed_qop']],p))]
            for j,row in enumerate(a[:first2+1]):
                label=f"S{int(row[f['station']])}L{int(row[f['layer']])}"
                history.append((label+'_transport',p_residual(row[f['q_predicted']],p)))
                if j<first2:
                    history.append((label+'_update',p_residual(row[f['q_filtered']],p)))
            counter[next(label for label,value in history if value<threshold)]+=1
    return {name:dict(counts) for name,counts in counters.items()}


def compare_paths(baseline, control):
    a, fields, paths_a = transitions(baseline)
    b, _, paths_b = transitions(control)
    c = {n:i for i,n in enumerate(fields)}
    def keyed(rows, paths):
        out = {}
        for row,path in zip(rows,paths):
            key=tuple(int(row[c[k]]) for k in ['file','entry','seed','track_barcode'])+(path,)
            assert key not in out, 'Ambiguous path matching'
            out[key]=row
        return out
    aa,bb=keyed(a,paths_a),keyed(b,paths_b)
    keys=sorted(aa.keys()&bb.keys())
    result=dict(baseline_eligible=len(aa),control_eligible=len(bb),same_seed_hits_truth=len(keys))
    if keys:
        x,y=np.asarray([aa[k] for k in keys]),np.asarray([bb[k] for k in keys])
        components=['s2_predicted','s2_first_filtered','s2_second_predicted','s2_smoothed']
        for component in components:
            rx=p_residual(x[:,c[component]],x[:,c['p']]);ry=p_residual(y[:,c[component]],y[:,c['p']])
            result[component]=dict(baseline=describe(rx),control=describe(ry),
                difference=describe(ry-rx),
                baseline_left=int(np.sum(rx<-.1)),control_left=int(np.sum(ry<-.1)),
                baseline_wrong_charge=int(np.sum(x[:,c[component]]*x[:,c['truth_qop']]<0)),
                control_wrong_charge=int(np.sum(y[:,c[component]]*y[:,c['truth_qop']]<0)),
                qop_difference_over_truth=describe((y[:,c[component]]-x[:,c[component]])/x[:,c['truth_qop']]),
                left_fraction_difference=event_interval((ry<-.1).astype(float)-(rx<-.1),x,
                    sum(item['events'] for item in baseline['coverage']),bounded=False))
        same=x[:,c['same_truth']]>0
        result['same_seed_track_truth']=dict(count=int(same.sum()),
            components={component:dict(
                baseline=describe_p(p_residual(x[same,c[component]],x[same,c['p']])),
                control=describe_p(p_residual(y[same,c[component]],y[same,c['p']])))
                for component in components})
    return result


def recorded_truth_bounds(root, data):
    chains=defaultdict(list)
    for path in sorted((root/'truth_vertices').glob('file*')):
        if int(path.name[4:]) not in [item['file_index'] for item in data['coverage']]:
            continue
        assert (path/'exit_code.txt').read_text().strip()=='0'
        file=int(path.name[4:])
        for entry,barcode,pdg,z,p,pz in np.load(path/'vertices.npz')['values']:
            chains[(file,int(entry),int(barcode)%1000000,int(pdg))].append((z,p,pz))
    assert chains
    v=data['values'];f={n:i for i,n in enumerate(data['fields'])}
    keep=(v[:,f['github']]>0)&(v[:,f['purity']]>.7)&(v[:,f['good_seed']]>0)&(v[:,f['measurement']]>0)
    counters=defaultdict(Counter)
    residual_intervals=defaultdict(list)
    for row in v[keep]:
        pdg=int(-11*row[f['truth_charge']])
        key=(int(row[0]),int(row[f['entry']]),int(row[f['track_barcode']])%1000000,pdg)
        chain=sorted(chains.get(key,[]))
        valid=bool(chain) and all(point[2]>0 for point in chain) and all(
            b[1]<=a[1]*(1+1e-6) for a,b in zip(chain,chain[1:]))
        before=[point for point in chain if point[0]<=row[f['z_F_mm']]]
        after=[point for point in chain if point[0]>row[f['z_F_mm']]]
        for component in ['predicted','filtered','smoothed']:
            label=f"S{int(row[f['station']])}L{int(row[f['layer']])}_{component}"
            counter=counters[label];counter['all']+=1
            p=1/abs(row[f['q_'+component]])
            if not valid or not before:
                counter['invalid_chain']+=1;continue
            upper=before[-1][1];lower=after[0][1] if after else 0.
            counter['valid_chain']+=1
            if lower>0 and row[f['same_seed_truth']]>0:
                # True residual is bracketed, without assuming an interpolation.
                residual_intervals[label].append([p/upper-1,p/lower-1,
                    upper/row[f['p_production_GeV']],lower/upper,
                    row[f['purity']],row[f['truth_charge']],row[f['p_production_GeV']]])
            if upper<.9*row[f['p_production_GeV']]:counter['recorded_upstream_loss_above_10pct']+=1
            for threshold,name in [(.9,'minus01'),(.5,'minus05')]:
                if p>=threshold*row[f['p_production_GeV']]:continue
                counter['production_left_'+name]+=1
                category='guaranteed_local' if p<threshold*lower else 'excluded_local' if p>=threshold*upper else 'unresolved'
                counter[name+'_'+category]+=1
    narrow={}
    for label,rows in residual_intervals.items():
        rows=np.asarray(rows)
        narrow[label]={}
        for ratio in [.99,.95]:
            selected=rows[rows[:,3]>=ratio]
            narrow[label][str(ratio)]=dict(states=len(selected),
                minimum_local_residual=describe_p(selected[:,0]),
                maximum_local_residual=describe_p(selected[:,1]),
                upper_bound_over_production=describe(selected[:,2]),
                purity95_minimum=describe_p(selected[selected[:,4]>=.95,0]),
                purity95_maximum=describe_p(selected[selected[:,4]>=.95,1]),
                charge={str(charge):dict(states=int(np.sum(selected[:,5]==charge)),
                    minimum=describe_p(selected[selected[:,5]==charge,0]),
                    maximum=describe_p(selected[selected[:,5]==charge,1])) for charge in [-1,1]},
                momentum={f'{low}_{high}':dict(states=int(np.sum((selected[:,6]>=low)&(selected[:,6]<high))),
                    minimum=describe_p(selected[(selected[:,6]>=low)&(selected[:,6]<high),0]),
                    maximum=describe_p(selected[(selected[:,6]>=low)&(selected[:,6]<high),1]))
                    for low,high in [(0,500),(500,1000),(1000,2000),(2000,7000)]})
    # A tight local interval alone does not exclude radiation later in the fit.
    # Require a bound spanning the entire measured trajectory separately.
    selected=v[keep&(v[:,f['same_seed_truth']]>0)]
    boundaries=np.r_[0,1+np.flatnonzero(np.any(np.diff(selected[:,:3],axis=0)!=0,axis=1)),len(selected)]
    span=defaultdict(list)
    for begin,end in zip(boundaries[:-1],boundaries[1:]):
        a=selected[begin:end]
        first2=np.flatnonzero(a[:,f['station']]==2)
        last3=np.flatnonzero(a[:,f['station']]==3)
        if not len(first2) or not len(last3):continue
        row=a[first2[0]];last=a[last3[-1]]
        key=(int(row[0]),int(row[f['entry']]),int(row[f['track_barcode']])%1000000,
             int(-11*row[f['truth_charge']]))
        chain=sorted(chains.get(key,[]))
        if not chain or not all(point[2]>0 for point in chain) or not all(
                b[1]<=a[1]*(1+1e-6) for a,b in zip(chain,chain[1:])):continue
        before=[point for point in chain if point[0]<=a[0,f['z_F_mm']]]
        after=[point for point in chain if point[0]>last[f['z_F_mm']]]
        if not before:continue
        upper=before[-1][1];lower=after[0][1] if after else 0.
        p=row[f['p_production_GeV']]
        values=[p_residual(row[f['q_predicted']],p),p_residual(row[f['q_smoothed']],p),
                p_residual(last[f['q_smoothed']],p),upper/p,lower/upper,row[f['purity']],
                p_residual(row[f['q_smoothed']],upper),
                p_residual(row[f['q_smoothed']],lower) if lower>0 else np.nan]
        span['all_valid'].append(values)
        for ratio in [.99,.95,.9]:
            if lower/upper>=ratio:span[f'entire_path_lower_over_upper_ge_{ratio}'].append(values)
    span_summary={name:dict(tracks=len(rows),
        s2_predicted=describe_p(np.asarray(rows)[:,0]),s2_smoothed=describe_p(np.asarray(rows)[:,1]),
        s3_last_smoothed=describe_p(np.asarray(rows)[:,2]),
        p_at_start_upper_over_production=describe(np.asarray(rows)[:,3]),
        path_lower_over_upper=describe(np.asarray(rows)[:,4]),
        s2_smoothed_minimum_local_residual=describe_p(np.asarray(rows)[:,6]),
        s2_smoothed_maximum_local_residual=describe_p(np.asarray(rows)[:,7])) for name,rows in span.items()}
    return dict(definition='Forward monotone same-PDG HepMC continuation bounds; zero lower bound without a later vertex; no interpolation or exact sensor momentum',
        per_layer={k:dict(v) for k,v in counters.items()},
        same_seed_truth_narrow_intervals=narrow,same_seed_truth_tracking_span=span_summary)


def plot(root, output, groups, summaries):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    output.mkdir(exist_ok=False)
    bins={}
    histogram_metadata={}
    def hist(ax,values,edges,key,**kwargs):
        values=np.asarray(values)
        finite=values[np.isfinite(values)]
        counts,_=np.histogram(finite,edges)
        bins[key+'_counts']=counts;bins[key+'_edges']=np.asarray(edges)
        bins[key+'_sumw2']=counts.copy()
        histogram_metadata[key]=dict(raw=len(values),finite=len(finite),
            underflow=int(np.sum(finite<edges[0])),overflow=int(np.sum(finite>edges[-1])))
        ax.stairs(counts,edges,**kwargs)
    def save(fig,name,pdf):
        fig.savefig(output/(name+'.png'),dpi=160)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig);plt.close(fig)
    data=groups['baseline']
    paired,cols,_=transitions(data);c={n:i for i,n in enumerate(cols)}
    scope=f"110051 | {len(data['coverage'])} files, {summaries['baseline']['source_events']} source events | latest Calypso"
    edges=np.linspace(-1.2,1.2,101)
    with PdfPages(output/'report.pdf') as pdf:
        fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
        for name,color in [('seed_qop','black'),('s1_filtered','red'),('s2_predicted','green'),('s2_smoothed','blue')]:
            residual=p_residual(paired[:,c[name]],paired[:,c['p']])
            hist(axes[0],residual,edges,name,label=name.replace('_',' '),color=color)
        axes[0].set(xlabel='(p - p production) / p production',ylabel='Tracks',title='Good-seed tracks with both S1 and S2')
        axes[0].legend(fontsize=8)
        residual1=p_residual(paired[:,c['s1_filtered']],paired[:,c['p']])
        residual2=p_residual(paired[:,c['s2_predicted']],paired[:,c['p']])
        h,xe,ye,_=axes[1].hist2d(residual1,residual2,bins=80,range=[[-1,1],[-1,1]],cmin=1)
        bins.update(transition2d_counts=np.nan_to_num(h),transition2d_xedges=xe,transition2d_yedges=ye)
        axes[1].plot([-1,1],[-1,1],color='red',lw=.8)
        axes[1].set(xlabel='Last S1 filtered p residual',ylabel='First S2 predicted p residual',title=f'Paired tracks: {len(paired)}')
        fig.suptitle(scope+'\nGitHub cuts + final purity > 0.7 + original good-seed gate; production truth')
        save(fig,'station2_origin',pdf)
        fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
        for same,ax in zip([False,True],axes):
            mask=paired[:,c['same_truth']]>0 if same else np.ones(len(paired),bool)
            for name,color,label in [('s2_predicted','green','First S2 sensor predicted'),
                    ('s2_first_filtered','red','First S2 sensor filtered'),
                    ('s2_second_predicted','orange','Next S2 sensor predicted'),
                    ('s2_smoothed','blue','First S2 sensor smoothed')]:
                residual=p_residual(paired[mask,c[name]],paired[mask,c['p']])
                hist(ax,residual,edges,('same_truth_' if same else 'all_')+name+'_sensors',
                    color=color,label=label)
            ax.set(xlabel='(p - p production) / p production',ylabel='Tracks',
                title=('Same seed / track truth' if same else 'Original good-seed cohort')+f' (N={mask.sum()})')
            ax.legend(fontsize=7)
        fig.suptitle(scope+'\nSeparate measurement sensors in forward path order; GitHub + purity > 0.7 + good seed')
        save(fig,'station2_sensors',pdf)
        fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
        for name,color in [('sum_transport','orange'),('sum_updates','purple')]:
            values=paired[:,c[name]]/paired[:,c['truth_qop']]
            hist(axes[0],values,np.linspace(-10,10,101),name,label=name.replace('_',' '),color=color)
        axes[0].set(xlabel='Accumulated change / signed truth q/p',ylabel='Tracks',title='Seed to first S2 predicted')
        axes[0].legend(fontsize=8)
        for component,color,label in [('s2_predicted','green','First S2 predicted'),
                ('s2_first_filtered','red','First S2 filtered'),('s2_smoothed','blue','First S2 smoothed')]:
            hist(axes[1],paired[:,c[component]]/paired[:,c['truth_qop']]-1,np.linspace(-20,20,101),
                component+'_curvature',color=color,label=label)
        axes[1].set(xlabel='Signed q/p / q/p truth - 1',ylabel='Tracks',title='Curvature before the absolute inverse')
        axes[1].legend(fontsize=7)
        hist(axes[2],paired[:,c['sigma_over_truth']],np.linspace(0,100,101),'s2_sigma')
        axes[2].set(xlabel='Predicted sigma(q/p) / |q/p truth|',ylabel='Tracks',title='First S2 curvature uncertainty')
        fig.suptitle(scope+'\nGood-seed gate; signed truth q/p from majority production particle')
        save(fig,'curvature_operations',pdf)
        seeds=data['seeds']
        # file, entry, seed, barcode, reco[6], truth[6]
        finite=np.all(np.isfinite(seeds[:,[6,7,8,12,13,14]]),axis=1)
        s=seeds[finite]
        good=(s[:,8]*s[:,14]>0)&(np.abs(np.abs(s[:,14]/s[:,8])-1)<.1)&(np.abs(s[:,7]-s[:,13])<1e-4)
        fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
        for mask,label in [(np.ones(len(s),bool),'All seeds with valid truth'),(good,'Original good-seed gate')]:
            prefix='good_seed' if label.startswith('Original') else 'all_seeds'
            hist(axes[0],np.abs(s[mask,14]/s[mask,8])-1,edges,prefix+'_p',label=label)
            hist(axes[1],(s[mask,6]-s[mask,12])*1e3,np.linspace(-10,10,101),prefix+'_phi',label=label)
            hist(axes[2],(s[mask,7]-s[mask,13])*1e3,np.linspace(-1,1,101),prefix+'_theta',label=label)
        tt=data['tracks'];tf={name:i for i,name in enumerate(data['track_fields'])}
        mask=(tt[:,tf['github']]>0)&(tt[:,tf['purity']]>.7)&(tt[:,tf['good_seed']]>0)
        lookup={tuple(x[:3].astype(int)):x for x in seeds}
        ancestors=np.asarray([lookup[(int(x[0]),int(x[tf['entry']]),int(x[tf['seed']]))] for x in tt[mask]])
        label="Selected tracks' ancestor seeds"
        hist(axes[0],np.abs(ancestors[:,14]/ancestors[:,8])-1,edges,'selected_seed_p',label=label)
        hist(axes[1],(ancestors[:,6]-ancestors[:,12])*1e3,np.linspace(-10,10,101),'selected_seed_phi',label=label)
        hist(axes[2],(ancestors[:,7]-ancestors[:,13])*1e3,np.linspace(-1,1,101),'selected_seed_theta',label=label)
        for ax,title,xlabel in zip(axes,['Momentum','Bending direction phi_A','Non-bending direction theta_A'],
                ['(p seed - p truth) / p truth','Seed - propagated truth [mrad]','Seed - propagated truth [mrad]']):
            ax.set(title=title,xlabel=xlabel,ylabel='Seeds');ax.legend(fontsize=7)
        fig.suptitle(scope+'\nNative RootSeedWriter: production particle propagated to seed surface')
        save(fig,'seed_residuals',pdf)
        ungated,ucols,_=transitions(data,require_good=False,same_truth=True)
        uc={name:i for i,name in enumerate(ucols)}
        seed_r,cuts=seed_gate_values(data,ungated,ucols)
        smooth_r=p_residual(ungated[:,uc['s2_smoothed']],ungated[:,uc['p']])
        fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
        for name,label in [('all','No seed gate'),('p_and_charge','Seed p + charge'),
                           ('charge_and_theta','Seed charge + theta'),('original','Original seed gate')]:
            mask=cuts[name];values=smooth_r[mask]
            counts,_=np.histogram(values[np.isfinite(values)],edges)
            key='seed_selection_'+name
            bins[key+'_counts']=counts;bins[key+'_edges']=edges;bins[key+'_sumw2']=counts.copy()
            histogram_metadata[key]=dict(raw=len(values),finite=int(np.isfinite(values).sum()),
                underflow=int(np.sum(values<edges[0])),overflow=int(np.sum(values>=edges[-1])))
            axes[0].stairs(counts/max(1,len(values)),edges,label=f'{label} (N={len(values)})')
        axes[0].set(xlabel='First S2 smoothed p residual to production truth',ylabel='Fraction of raw tracks / bin',
                    title='Selections of the same native reconstruction')
        axes[0].legend(fontsize=7)
        h,xe,ye,_=axes[1].hist2d(seed_r,smooth_r,bins=80,range=[[-.5,.5],[-.5,.5]],cmin=1)
        bins.update(seed_smooth2d_counts=np.nan_to_num(h),seed_smooth2d_xedges=xe,seed_smooth2d_yedges=ye)
        for cut in [-.1,.1]:axes[1].axvline(cut,color='red',lw=.8)
        axes[1].set(xlabel='Seed p residual to propagated seed truth',ylabel='First S2 smoothed p residual',
                    title='Same seed / track truth; before seed gate')
        fig.suptitle(scope+'\nGitHub + purity > 0.7; same seed / track truth; production p; isolation ON')
        save(fig,'seed_gate_selection',pdf)
        def controls_figure(names,files,filename):
            fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
            for name in names:
                d=subset_files(groups[name],files)
                a,columns,_=transitions(d);cc={n:i for i,n in enumerate(columns)}
                for ax,component in zip(axes,['s2_predicted','s2_smoothed']):
                    hist(ax,p_residual(a[:,cc[component]],a[:,cc['p']]),edges,filename+'_'+name+'_'+component,
                         label=f'{name} (N={len(a)})')
                    ax.set(title=component.replace('_',' '),xlabel='(p - p production) / p production',ylabel='Tracks')
            for ax in axes:ax.legend(fontsize=7)
            events=sum(x['events'] for x in subset_files(groups['baseline'],files)['coverage'])
            fig.suptitle(f'110051 | {len(files)} files, {events} source events per configuration | latest Calypso\nOriginal seed gate; each configuration selects its own tracks; production truth')
            save(fig,filename,pdf)
        if len(groups)>1:
            names=[name for name,d in groups.items() if any(x['file_index']==0 for x in d['coverage'])]
            if len(names)>1:controls_figure(names,{0},'controls_file0')
            complete=[name for name,d in groups.items() if len(d['coverage'])==len(data['coverage'])]
            if len(complete)>1 and len(data['coverage'])>1:
                controls_figure(complete,{x['file_index'] for x in data['coverage']},'controls_full')
        closure=summaries['baseline'].get('native_seed_field_closure')
        if closure:
            rows=json.loads((root/'seed_field_closure_v3/closure.json').read_text())['rows']
            particles=[p for event in rows for p in event['particles'] if p['success']]
            fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
            hist(axes[0],[p['native_circle_residual'] for p in particles],
                np.linspace(-.2,.4,101),'seed_model_residual',color='black')
            hist(axes[1],[p['effective_circle_field_T'] for p in particles],
                np.linspace(.3,.7,101),'seed_model_effective_field',color='purple')
            axes[0].set(xlabel='(p native circle - p input) / p input',ylabel='Particles',
                title='Native CircleFit with +/-30 mm points and B = 0.55 T')
            axes[1].set(xlabel='p input / (0.3 R) [T]',ylabel='Particles',
                title='Effective field for the circle estimator')
            fig.suptitle(f"110051 file 0 | {closure['source_events']} source events | {closure['successful_particles']} particles\nNative Calypso full-field propagation; no measurement noise, energy loss or scattering; p > 10 GeV")
            save(fig,'seed_model_closure',pdf)
            if any(x['file_index']!=0 for x in data['coverage']):
                validation,rows,vc,masks=model_gate_validation(data,closure)
                fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
                for name,color,label in [('no_seed_gate','black','No seed gate'),
                        ('original','red','Original seed p gate'),
                        ('native_model_calibrated','blue','Native-model calibrated seed p gate')]:
                    residual=p_residual(rows[masks[name],vc['s2_smoothed']],rows[masks[name],vc['p']])
                    finite=residual[np.isfinite(residual)]
                    counts,_=np.histogram(finite,edges)
                    key='model_validation_'+name
                    bins[key+'_counts']=counts;bins[key+'_edges']=edges;bins[key+'_sumw2']=counts.copy()
                    ax.stairs(counts/len(finite),edges,color=color,label=label+f' (N={len(residual)})')
                    histogram_metadata[key]=dict(raw=len(residual),finite=len(finite),
                        underflow=int(np.sum(finite<edges[0])),overflow=int(np.sum(finite>edges[-1])),normalization='Raw finite tracks')
                ax.set(xlim=(-.5,.5),xlabel='First S2 smoothed (p - p production) / p production',
                       ylabel='Fraction of selected tracks per bin')
                ax.legend(fontsize=8)
                fig.suptitle(f"110051 | {validation['evaluation_source_events']} source events excluding file 0 | same seed/track truth\nAnalysis cuts only; scale {validation['factor']:.5f} from 100-event native seed-model closure")
                save(fig,'seed_model_gate_validation',pdf)
    np.savez_compressed(output/'bin_counts.npz',edges=edges,**bins)
    (output/'histogram_metadata.json').write_text(json.dumps(histogram_metadata,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--groups',nargs='+',default=['baseline'])
    parser.add_argument('--files',type=int,nargs='+',help='Explicitly labelled subset for mechanism pilots')
    parser.add_argument('--control-files',type=int,nargs='+',help='Explicit control coverage when baseline uses more input files')
    parser.add_argument('--summary-only',action='store_true')
    args=parser.parse_args()
    groups={name:load_group(args.root,name,args.files if name=='baseline' or args.control_files is None
                           else args.control_files) for name in args.groups}
    frozen=json.loads((args.root.parent.parent/'stage04_seed_purity/20261113_seed_good_v2/inputs.json').read_text())
    expected={x['path']:x['entries'] for x in frozen['samples'][0]['files']}
    coverage=groups['baseline']['coverage']
    assert len(expected)==20
    assert len(coverage)==(len(expected) if args.files is None else len(args.files))
    assert all(x['skip']==0 and x['events']==expected[x['file']] for x in coverage)
    summaries={name:summarize(data) for name,data in groups.items()}
    summaries['baseline']['recorded_truth_bounds']=recorded_truth_bounds(args.root,groups['baseline'])
    closure_path=args.root/'seed_field_closure_v3/closure.json'
    if closure_path.exists():
        closure=json.loads(closure_path.read_text())
        particles=[p for event in closure['rows'] for p in event['particles']]
        successful=[p for p in particles if p['success']]
        summaries['baseline']['native_seed_field_closure']=dict(
            source=str(closure_path),source_events=closure['events'],particles=len(particles),
            successful_particles=len(successful),definition=closure['definition'],
            residual=describe_p([p['native_circle_residual'] for p in successful]),
            effective_circle_field_T=describe([p['effective_circle_field_T'] for p in successful]))
        if any(x['file_index']!=0 for x in groups['baseline']['coverage']):
            summaries['baseline']['model_gate_validation']=model_gate_validation(
                groups['baseline'],summaries['baseline']['native_seed_field_closure'])[0]
    for name,data in groups.items():
        if name!='baseline':
            # Restrict baseline to the control's frozen source coverage.
            base=groups['baseline']
            files=[x['file_index'] for x in data['coverage']]
            subset=subset_files(base,files)
            summaries[name]['matched_paths']=compare_paths(subset,data)
    if args.summary_only:
        args.output.mkdir(exist_ok=False)
        (args.output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
        return
    plot(args.root,args.output,groups,summaries)
    (args.output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    import ROOT
    ROOT.gROOT.SetBatch(True)
    import subprocess
    reference_figures=[]
    raw_entries={};reference_bins={}
    labels=dict(baseline='latest baseline',truth_qop='truth seed q/p',truth_mean='truth seed mean',
        qop_prior_100='q/p prior x100',qop_prior_1='q/p prior x1',truth_qop_prior1='truth q/p + prior x1')
    for group,data in groups.items():
        counts=args.output/(group+'_counts_3station_forward.root')
        prefix=args.output/('110051_latest_good_seed' if group=='baseline' else '110051_'+group+'_nine_layers')
        # Worker files also contain canvases. Merge only the 27 reference
        # histograms; ROOT canvas ownership makes generic TFileMerger unsafe.
        merged_histograms={};expected_entries={};expected_counts={}
        names=[f'h_{component}_{station}_{layer}' for component in range(3)
               for station in [1,2,3] for layer in [0,1,2]]
        for path in data['paths']:
            source_file=ROOT.TFile.Open(str(path/'native_counts.root'))
            for name in names:
                source_hist=source_file.Get(name)
                assert source_hist and source_hist.GetNbinsX()==100
                values=np.asarray([source_hist.GetBinContent(i) for i in range(102)])
                if name not in merged_histograms:
                    merged_histograms[name]=source_hist.Clone(name)
                    merged_histograms[name].SetDirectory(0)
                    expected_entries[name]=source_hist.GetEntries()
                    expected_counts[name]=values
                else:
                    assert merged_histograms[name].Add(source_hist)
                    expected_entries[name]+=source_hist.GetEntries()
                    expected_counts[name]+=values
            source_file.Close()
        merged_file=ROOT.TFile(str(counts),'RECREATE')
        for name,histogram in merged_histograms.items():
            assert histogram.GetEntries()==expected_entries[name]
            np.testing.assert_equal([histogram.GetBinContent(i) for i in range(102)],expected_counts[name])
            histogram.Write(name)
        merged_file.Close()
        env={**os.environ,'REFERENCE_SCOPE_LABEL':f"110051 | {len(data['coverage'])} files / {summaries[group]['source_events']} source events | ACTS A | HepMC p | GitHub + good seed + purity>0.7"}
        invocation=str(args.root/'reference/ckf_residual_perLayer.C')+'('+','.join(
            json.dumps(x) for x in [str(counts),str(prefix),labels.get(group,group),'t_eQOP'])+',true)'
        subprocess.run(['root','-l','-b','-q',invocation],env=env,check=True)
        reference_figures.append(prefix.with_suffix('.pdf'))
        actual=ROOT.TFile.Open(str(prefix.with_suffix('.root')))
        source=ROOT.TFile.Open(str(counts))
        raw_entries[group]={}
        for component in range(3):
            for station in [1,2,3]:
                for layer in [0,1,2]:
                    name=f'h_{component}_{station}_{layer}'
                    a,b=actual.Get(name),source.Get(name)
                    assert a.GetEntries()==b.GetEntries()
                    assert all(a.GetBinContent(i)==b.GetBinContent(i) for i in range(102))
                    raw_entries[group][name]=a.GetEntries()
                    reference_bins[group+'_'+name+'_counts_with_flow']=np.asarray([a.GetBinContent(i) for i in range(102)])
                    reference_bins[group+'_'+name+'_sumw2_with_flow']=np.asarray([a.GetBinError(i)**2 for i in range(102)])
        actual.Close();source.Close()
    np.savez_compressed(args.output/'reference_bin_counts.npz',edges=np.linspace(-1.2,1.2,101),**reference_bins)
    import sys
    sys.path.insert(0,str(args.root.parent.parent/'stage04_seed_purity/20261113_seed_good_v2/software_tools'))
    from pypdf import PdfReader, PdfWriter
    from PIL import Image
    (args.output/'report.pdf').rename(args.output/'diagnostics.pdf')
    merged=PdfWriter()
    for path in reference_figures:merged.append(str(path))
    merged.append(str(args.output/'diagnostics.pdf'))
    with (args.output/'report.pdf').open('wb') as stream:
        merged.write(stream)
    index=[]
    for png in sorted(args.output.glob('*.png')):
        with Image.open(png) as image:
            image.load();size=list(image.size)
        pdf=png.with_suffix('.pdf')
        assert len(PdfReader(pdf).pages)==1
        index.append(dict(png=png.name,pdf=pdf.name,pixels=size))
    (args.output/'figure_validation.json').write_text(json.dumps(dict(
        figure_index=index,pdf_pages=len(PdfReader(args.output/'report.pdf').pages),
        reference_histograms=27*len(groups),reference_bins_with_flow=2754*len(groups),
        reference_raw_entries=raw_entries,reference_histograms_exact=True),indent=2)+'\n')


if __name__=='__main__':
    main()
