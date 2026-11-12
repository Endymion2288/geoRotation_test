#!/usr/bin/env python3
"""Quantify native Calypso angle definitions, structure and held-out response."""
import argparse
import json
from pathlib import Path

import numpy as np

from slope_tables import angle, event_first_mask, histogram, join_truth, merge_tables, physical, slope, write_json
from slope_statistics import conditional_response, discover_windows, distribution_structure, peak_statistics, event_ratio_interval


EDGES = np.linspace(-.01,.01,101)
PLANES = ['st1','VetoNu','sensor']


def derive(tables, variant, dsid):
    tracks, truth = tables['tracks'],tables['truth']
    index, matched = join_truth(tables)
    long = tracks['station_bits'] & 14 == 14
    p = np.linalg.norm(tracks['upstream_F'][:,3:],axis=1)
    prod_p = np.linalg.norm(tracks['production_momentum_F_GeV'],axis=1)
    pdg = -13 if dsid in ['100044','100116'] else 13 if dsid == '100117' else 11
    signal_truth = (truth['status']==1) & ((np.abs(truth['pdg'])==11) if abs(pdg)==11 else (truth['pdg']==pdg))
    pos = truth['production_position_F_mm']
    station = truth['truth_station_position_F_mm'][:,1:,:2]
    fiducial = (np.hypot(pos[:,0],pos[:,1])<100) & np.all(np.linalg.norm(station,axis=2)<100,axis=1)
    fiducial &= np.sum(truth['truth_hit_counts'][:,1:],axis=1)>=18
    quality = (tracks['n_layers']>=7) & (tracks['ndof']+5>=12) & (tracks['rmax_mm']<95)
    with np.errstate(divide='ignore',invalid='ignore'):
        quality &= tracks['chi2']/tracks['ndof']<25
    pure = tracks['purity']>=.7
    primary = matched & signal_truth[index]
    high_quality = long & quality & primary & pure & fiducial[index]
    pwindow = [50.,100.] if abs(pdg)==11 else [50.,250.]
    selections = dict(all=long,good95=long & quality,matched_pure=high_quality,
        reco_p_window=long & (p>pwindow[0]) & (p<pwindow[1]),
        truth_p_window=long & primary & (prod_p>pwindow[0]) & (prod_p<pwindow[1]),
        legacy_first_pz=event_first_mask(tracks['sample_entry'],tracks['upstream_raw'][:,5],long))
    truth_direction = truth['truth_model_momentum_F_GeV']
    model = np.stack([slope(truth_direction[:,i]) for i in range(2)],axis=1)
    model = np.where(truth['truth_model_valid']==1,model,np.nan)
    sensor = slope(truth['truth_sensor_direction_F'])
    truth_slopes = np.column_stack([model,sensor])
    reconstructed = np.stack([slope(tracks['projection'][:,i,3:6]) for i in range(3)],axis=1)
    reconstructed = np.where(tracks['projection_valid']==1,reconstructed,np.nan)
    residuals = reconstructed-truth_slopes[index]
    residuals[~matched,:] = np.nan
    with np.errstate(divide='ignore',invalid='ignore'):
        raw_st1 = np.arctan(tracks['upstream_raw'][:,4]/tracks['upstream_raw'][:,5])
        p_ratio = p/prod_p
    legacy_veto = np.where(tracks['projection_valid'][:,1]==1,tracks['projection'][:,1,7],999999.)
    # Literal enabled GitHub macro: before theta_x with 100<raw pz<250;
    # after uses projected y/sqrt(x^2+y^2), with no p selection.
    if variant=='beforeGeo':
        with np.errstate(divide='ignore',invalid='ignore'):
            literal = np.arctan(tracks['projection'][:,1,3]/tracks['projection'][:,1,5])
        literal = np.where(tracks['projection_valid'][:,1]==1,literal,999999.)
        literal_mask = long & (tracks['upstream_raw'][:,5]>100) & (tracks['upstream_raw'][:,5]<250)
    else:
        # The rotated writer stores X=z_A=-x_F and Y=y_A=y_F.
        # Its X sign cancels inside this particular squared denominator.
        xy = tracks['projection'][:,1,:2]
        with np.errstate(divide='ignore',invalid='ignore'):
            literal = xy[:,1]/np.linalg.norm(xy,axis=1)
        literal = np.where(tracks['projection_valid'][:,1]==1,literal,1/np.sqrt(2))
        literal_mask = long
    return dict(index=index,matched=matched,primary=primary,p=p,prod_p=prod_p,p_ratio=p_ratio,
        long=long,quality=quality,pure=pure,high_quality=high_quality,signal_truth=signal_truth,
        fiducial_truth=fiducial & signal_truth,selections=selections,pwindow=pwindow,
        truth_slopes=truth_slopes,reconstructed=reconstructed,residuals=residuals,
        raw_st1=raw_st1,legacy_veto=legacy_veto,literal=literal,literal_mask=literal_mask,
        wrong_charge=primary & (tracks['charge']*np.where(truth['pdg'][index]>0,-1,1)<0),
        momentum_tail=primary & (p_ratio<.5),rotated=variant!='beforeGeo')


def analyse_variant(tables, derived, coverage, windows, split, dsid):
    tracks,truth = tables['tracks'],tables['truth']
    d = derived
    long = d['long']
    end = int(tables['events']['sample_entry'].max(initial=-1)+1)
    begin = int(tables['events']['sample_entry'].min(initial=end))
    counts = {}
    for plane in range(3):
        for name,selected in d['selections'].items():
            counts[f'{PLANES[plane]}|{name}|physical'] = histogram(d['reconstructed'][selected,plane],EDGES)
            if name in ['all','matched_pure','truth_p_window']:
                counts[f'{PLANES[plane]}|{name}|residual'] = histogram(d['residuals'][selected,plane],np.linspace(-.005,.005,201))
        values = d['reconstructed'][long,plane]
        for width in [1e-4,2e-4,4e-4]:
            for offset in [0.,.5]:
                edges = np.arange(-.01,.01000001,width)+offset*width
                counts[f'{PLANES[plane]}|bin_{width:g}_offset_{offset:g}'] = histogram(values,edges)
    for name,values in [('raw_st1',d['raw_st1']),('legacy_veto',d['legacy_veto'])]:
        for selection in ['all','legacy_first_pz','reco_p_window']:
            counts[f'{name}|{selection}'] = histogram(values[d['selections'][selection]],EDGES)
    literal_edges = np.linspace(-.01,.01,101 if not d['rotated'] else 51)
    counts['literal_github_macro'] = histogram(d['literal'][d['literal_mask']],literal_edges)
    for plane in range(3):
        values = d['truth_slopes'][d['fiducial_truth'],plane]
        counts[f'{PLANES[plane]}|truth_fiducial'] = histogram(values,EDGES)
    native_f = physical(tracks['upstream_raw'][:,3:],d['rotated'])
    raw = tracks['upstream_raw'][:,3:]
    with np.errstate(divide='ignore',invalid='ignore'):
        alternate = raw[:,1]/raw[:,0] if d['rotated'] else raw[:,1]/raw[:,2]
    common = np.isfinite(alternate) & np.isfinite(slope(native_f))
    frame_error = np.max(np.abs(alternate[common]-slope(native_f)[common]),initial=0)
    metrics = dict(source_events=coverage['source_events'],native_tracks=len(long),long_tracks=int(long.sum()),
        primary_truth=int(d['signal_truth'].sum()),fiducial_H18_truth=int(d['fiducial_truth'].sum()),
        high_quality_matched_tracks=int(d['high_quality'].sum()),frame_formula_max_absolute_error=float(frame_error),
        raw_F_momentum_max_error_GeV=float(np.max(np.abs(native_f-tracks['upstream_F'][:,3:]),initial=0)),
        upstream_surface_differs_from_st1=int(np.count_nonzero(long & (np.abs(tracks['upstream_F'][:,2])>1e-6))),
        legacy_st1_outside_plot=int(np.count_nonzero(long & (np.abs(d['raw_st1'])>=.01))),
        legacy_first_pz_selected_tracks=int(d['selections']['legacy_first_pz'].sum()),
        recoil_window_GeV=d['pwindow'],production_truth_is_local_G4=False,
        covariance_mode='actual transformed track covariance' if coverage.get('actual_covariance') else 'identity (native Ntuple convention)')
    correct_charge = d['primary'] & ~d['wrong_charge']
    for name,mask in [('wrong_charge',d['wrong_charge']),('momentum_tail',d['momentum_tail'])]:
        metrics[name] = event_ratio_interval(tracks['sample_entry'],mask & long,d['primary'] & long,(begin,end))
    efficient = np.bincount(d['index'][d['high_quality']],minlength=len(truth['barcode']))>0
    metrics['truth_efficiency'] = event_ratio_interval(truth['sample_entry'],efficient & d['fiducial_truth'],d['fiducial_truth'],(begin,end))
    metrics['propagation'] = {}
    for i,plane in enumerate(PLANES):
        metrics['propagation'][plane] = dict(
            success=int(np.count_nonzero(long & (tracks['projection_valid'][:,i]==1))),
            failed=int(np.count_nonzero(long & (tracks['projection_valid'][:,i]==0))),
            missing_reference=int(np.count_nonzero(long & (tracks['projection_valid'][:,i]<0))),
            invalid_slope=int(np.count_nonzero(long & ~np.isfinite(d['reconstructed'][:,i]))))
    structures,peaks,response = {},{},{ }
    validation_range = (max(begin,split),end)
    for i,plane in enumerate(PLANES):
        for selection in ['all','reco_p_window','truth_p_window','matched_pure']:
            mask = d['selections'][selection]
            key = plane+'|'+selection
            structures[key] = distribution_structure(d['reconstructed'][mask,i])
            if validation_range[1]>validation_range[0]:
                peaks[key] = peak_statistics(d['reconstructed'][:,i],tracks['sample_entry'],mask,windows.get(key,windows[plane]),validation_range)
        if i<2 and begin<split<end:
            for selection in ['all','reco_p_window','truth_p_window']:
                key = plane if selection=='all' else plane+'|'+selection
                truth_population = d['signal_truth'].copy()
                if selection=='truth_p_window':
                    tp = np.linalg.norm(truth['production_momentum_F_GeV'],axis=1)
                    truth_population &= (tp>d['pwindow'][0]) & (tp<d['pwindow'][1])
                response[key] = conditional_response(d['truth_slopes'][:,i],
                    np.linalg.norm(truth['production_momentum_F_GeV'],axis=1),
                    truth['truth_sensor_position_F_mm'][:,0],truth['truth_sensor_position_F_mm'][:,1],
                    truth['sample_entry'],d['reconstructed'][:,i],d['index'],
                    d['primary'] & d['selections'][selection] & d['pure'],tracks['sample_entry'],split,EDGES,
                    windows=windows.get(plane+'|'+selection,windows[plane]),truth_selected=truth_population)
    groups = {
        'production_p_GeV':(d['prod_p'],[0,50,100,250,500,1000,np.inf]),
        'module_eta':(tracks['module_eta'],[-2,-1,0,1,2]),
        'module_phi':(tracks['module_phi'],[-.5,.5,1.5,2.5,3.5,4.5]),
        'number_of_hits':(tracks['n_hits'],[0,12,15,18,21,1000]),
        'cluster_width':(tracks['cluster_width_mean'],[0,1.01,1.5,2,3,1000]),
        'purity':(tracks['purity'],[0,.5,.7,.9,1.01]),
        'production_slope':(slope(tracks['production_momentum_F_GeV']),[-np.inf,-.004,-.002,0,.002,.004,np.inf]),
        'production_x_mm':(truth['production_position_F_mm'][d['index'],0],[-np.inf,-50,0,50,np.inf]),
        'true_sensor_y_mm':(truth['truth_sensor_position_F_mm'][d['index'],1],[-np.inf,-50,0,50,np.inf]),
    }
    for group,(values,edges) in groups.items():
        for j,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
            selected = long & d['primary'] & (values>=lo) & (values<hi)
            for i,plane in enumerate(PLANES):
                counts[f'{plane}|group_{group}_{j}'] = histogram(d['reconstructed'][selected,i],EDGES)
                counts[f'{plane}|reco_p_window|group_{group}_{j}'] = histogram(d['reconstructed'][selected & d['selections']['reco_p_window'],i],EDGES)
    for label,mask in [('momentum_tail',d['momentum_tail']),('momentum_normal',d['primary'] & ~d['momentum_tail']),
                       ('wrong_charge',d['wrong_charge']),('correct_charge',correct_charge)]:
        for i,plane in enumerate(PLANES):
            counts[f'{plane}|{label}'] = histogram(d['reconstructed'][long & mask,i],EDGES)
    maps = {}
    for i,plane in enumerate(PLANES):
        for name,values,edges in [('production_p_GeV',d['prod_p'],[0,25,50,75,100,150,250,400,700,1000,2000,5000]),
                                 ('true_sensor_y_mm',truth['truth_sensor_position_F_mm'][d['index'],1],np.linspace(-100,100,51)),
                                 ('module_eta',tracks['module_eta'],[-2,-1,0,1,2])]:
            mask = long & d['primary'] & np.isfinite(d['reconstructed'][:,i]) & np.isfinite(values)
            h,x,y = np.histogram2d(d['reconstructed'][mask,i],values[mask],bins=[EDGES,edges])
            maps[plane+'|'+name] = dict(slope_edges=x.tolist(),condition_edges=y.tolist(),counts=h.astype(int).tolist(),
                entries=int(mask.sum()),out_of_map=int(mask.sum()-h.sum()))
        mask = long & d['primary'] & np.isfinite(d['residuals'][:,i])
        variance = tracks['projection'][:,i,6]
        valid = mask & np.isfinite(variance) & (variance>0)
        counts[f'{plane}|matched|pull'] = histogram(d['residuals'][valid,i]/np.sqrt(variance[valid]),np.linspace(-10,10,201))
        metrics['propagation'][plane]['pull_invalid_variance'] = int(np.count_nonzero(mask & ~valid))
    return dict(metrics=metrics,histograms=counts,structure=structures,peaks=peaks,response=response,maps=maps,
                groups={k:[None if not np.isfinite(x) else x for x in edges] for k,(_,edges) in groups.items()})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--dsid',default='100044')
    parser.add_argument('--variants',nargs='+',choices=['beforeGeo','afterGeo'],default=['beforeGeo','afterGeo'])
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    manifest = json.loads((args.root/'manifest.json').read_text())
    split = manifest['held_out_split_sample_entry']
    data,derived,coverage = {},{},{}
    for variant in args.variants:
        data[variant],coverage[variant] = merge_tables(args.root,variant)
        derived[variant] = derive(data[variant],variant,args.dsid)
        first_task = json.loads((args.root/(variant+'_tasks.json')).read_text())[0]
        source = json.loads((Path(first_task['output'])/'source_task.json').read_text())
        coverage[variant]['actual_covariance'] = '--actual-covariance' in source['options']
    frozen = args.output/'peak_windows.json'
    if frozen.exists():
        windows = json.loads(frozen.read_text())['windows']
    else:
        if 'afterGeo' not in data:
            raise ValueError('Discover and freeze windows on afterGeo before analysing beforeGeo alone')
        d,tracks = derived['afterGeo'],data['afterGeo']['tracks']
        windows = {}
        for i,plane in enumerate(PLANES):
            for selection in ['all','reco_p_window','truth_p_window','matched_pure']:
                selected = d['selections'][selection] & (tracks['sample_entry']<split)
                key = plane if selection=='all' else plane+'|'+selection
                windows[key] = discover_windows(d['reconstructed'][selected,i])
        write_json(frozen,dict(source=str(args.root),definition='AfterGeo discovery entries only; all native long tracks',
            held_out_boundary=split,windows=windows))
    for variant in data:
        result = analyse_variant(data[variant],derived[variant],coverage[variant],windows,split,args.dsid)
        write_json(args.output/(variant+'_summary.json'),result)
        d,tracks = derived[variant],data[variant]['tracks']
        peak_rows = []
        for i,plane in enumerate(PLANES):
            for w,window in enumerate(windows[plane]):
                distance = np.abs(d['reconstructed'][:,i]-window['center'])/window['half_width']
                mask = d['long'] & (tracks['sample_entry']>=split) & (distance<3)
                rows = np.flatnonzero(mask)
                peak_rows.extend(np.column_stack([tracks['sample_entry'][rows],tracks['barcode'][rows],tracks['ordinal'][rows],
                                                  np.full(len(rows),i),np.full(len(rows),w),(distance[rows]>=1).astype(int)]))
        np.savez_compressed(args.output/(variant+'_peak_events.npz'),
            rows=np.asarray(peak_rows,dtype=np.int64).reshape(-1,6),
            columns=np.array(['sample_entry','barcode','ordinal','plane','window','sideband']))
    for value in coverage.values():
        value['held_out_split_sample_entry'] = split
    if (args.output/'coverage.json').exists():
        coverage = json.loads((args.output/'coverage.json').read_text()) | coverage
    write_json(args.output/'coverage.json',coverage)
    print(json.dumps({v:coverage[v]['source_events'] for v in coverage}),flush=True)


if __name__ == '__main__':
    main()
