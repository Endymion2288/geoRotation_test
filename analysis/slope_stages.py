#!/usr/bin/env python3
"""Compare audited native stages after propagation by the actual Calypso tool."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_tables import histogram, join_truth, slope
from slope_statistics import event_ratio_interval
from slope_controls import archive, save


LABELS = {1:'Seed',2:'Selected CKF target',3:'CKF predicted',4:'CKF filtered',5:'CKF smoothed',
          6:'ACTS to Trk',7:'Refit input',8:'Refit predicted',9:'Refit filtered',10:'Refit smoothed',11:'Final'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--states',type=Path,required=True)
    parser.add_argument('--transport',type=Path,required=True)
    parser.add_argument('--windows',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    windows = json.loads(args.windows.read_text())['windows']
    observations,rows,features = [],[],[]
    closure = []
    offset = 0
    for shard in sorted((args.states/'shards').glob('s*')):
        target=args.transport/'shards'/shard.name
        assert (shard/'exit_code.txt').read_text().strip()=='0' and (target/'exit_code.txt').read_text().strip()=='0'
        state,transport = archive(shard/'states.npz'),archive(target/'transport.npz')
        r,t = archive(shard/'tracks.npz'),archive(shard/'truth.npz')
        f = archive(shard/'candidate_features.npz')
        assert np.array_equal(state['entry'],transport['source_entry']) and np.array_equal(state['stage'],transport['stage'])
        r['sample_entry'],t['sample_entry'] = r['source_entry'],t['source_entry']
        index,matched=join_truth({'tracks':r,'truth':t})
        finals=np.flatnonzero(transport['stage']==11)
        assert len(finals)==len(r['source_entry']) and np.array_equal(transport['source_entry'][finals],r['source_entry'])
        common=np.isfinite(transport['projection_F'][finals]) & np.isfinite(r['projection'][:,:,:6])
        closure.append(float(np.max(np.abs(transport['projection_F'][finals][common]-r['projection'][:,:,:6][common]),initial=0)))
        assert closure[-1]<1e-10
        keys = list(zip(transport['source_entry'][finals],transport['seed'][finals],transport['candidate'][finals]))
        all_events=archive(shard/'events.npz')['source_entry']
        dense={int(entry):offset+i for i,entry in enumerate(all_events)}
        offset+=len(all_events)
        pure=(r['station_bits']&14==14) & matched & (t['pdg'][index]==-13) & (r['purity']>=.7)
        final_p=np.linalg.norm(r['upstream_F'][:,3:],axis=1)
        by_candidate={key:j for j,key in enumerate(keys)}
        feature_index={(e,s,c):i for i,(e,s,c) in enumerate(zip(f['entry'],f['seed'],f['candidate']))}
        representatives={}
        for k in range(len(transport['stage'])):
            key=(transport['source_entry'][k],transport['seed'][k],transport['candidate'][k],transport['stage'][k])
            old=representatives.get(key)
            if old is None or transport['input_position_F_mm'][k,2]<transport['input_position_F_mm'][old,2]:representatives[key]=k
        for key,j in by_candidate.items():
            if not pure[j]:continue
            ti=index[j]
            truth_slopes=np.r_[slope(t['truth_model_momentum_F_GeV'][ti]),slope(t['truth_sensor_direction_F'][ti])]
            projection=np.where(transport['projection_valid'][finals[j]]==1,slope(transport['projection_F'][finals[j],:,3:]),np.nan)
            fi=feature_index[key]
            observation=dict(entry=int(key[0]),seed=int(key[1]),candidate=int(key[2]),barcode=int(r['barcode'][j]),dense_entry=dense[int(key[0])],
                final_p_GeV=float(final_p[j]),truth_p_GeV=float(np.linalg.norm(t['production_momentum_F_GeV'][ti])),
                charge=float(r['charge'][j]),hits=int(r['n_hits'][j]),module_eta=int(r['module_eta'][j]),module_phi=int(r['module_phi'][j]),
                true_x_mm=float(t['truth_sensor_position_F_mm'][ti,0]),true_y_mm=float(t['truth_sensor_position_F_mm'][ti,1]),
                final_slopes=projection,truth_slopes=truth_slopes,
                features={k:f[k][fi].item() for k in f if k not in ['entry','seed','candidate','barcode']})
            observations.append(observation)
            for code in LABELS:
                representative=representatives.get((*key,code)) if code!=1 else representatives.get((key[0],key[1],-1,1))
                if representative is None:continue
                k=representative
                values=np.where(transport['projection_valid'][k]==1,slope(transport['projection_F'][k,:,3:]),np.nan)
                rows.append(dict(dense_entry=observation['dense_entry'],entry=int(key[0]),seed=int(key[1]),candidate=int(key[2]),stage=code,
                    reco_p_window=50<final_p[j]<250,values=values,truth=truth_slopes))
        selected=np.isin(state['stage'],[3,4,5])
        cohort_keys={key for key,j in by_candidate.items() if pure[j]}
        cohort=np.array([(e,s,c) in cohort_keys for e,s,c in zip(state['entry'],state['seed'],state['candidate'])])
        for label,codes in [('CKF',[3,4,5]),('refit',[8,9,10])]:
            mask=cohort & (state['stage']==codes[-1])
            features.append(dict(shard=shard.name,stage=label,states=int(mask.sum()),
                update_angular_identity_max=float(np.nanmax(np.abs(state['update_error'][mask,2:4]))),
                smoothing_angular_identity_max=float(np.nanmax(np.abs(state['smoothing_error'][mask,2:4]))),
                smoothing_angular_identity_above_1e_minus_8=int(np.count_nonzero(np.any(np.abs(state['smoothing_error'][mask,2:4])>1e-8,axis=1)))))
    stages=np.array([r['stage'] for r in rows]);entry=np.array([r['dense_entry'] for r in rows])
    values=np.asarray([r['values'] for r in rows]);truth=np.asarray([r['truth'] for r in rows])
    pwindow=np.array([r['reco_p_window'] for r in rows])
    result=dict(source_events=offset,matched_pure_native_long_tracks=len(observations),
        final_projection_max_absolute_closure=max(closure),stage_definitions=LABELS,
        measurement_state_definition='Most upstream measurement state on each candidate, P/F/S and refit P/F/S; all projected to common physical planes',
        truth_reference='SiHit chord at sensor; Calypso-transported model at st1 and VetoNu',
        matrix_identities=features,planes={},matched_controls=[],candidate_groups={})
    report=PdfPages(args.output/'stage_report.pdf')
    for i,plane in enumerate(['st1','VetoNu','sensor']):
        for selection in ['all','reco_p_window']:
            selected=np.ones(len(rows),dtype=bool) if selection=='all' else pwindow
            hist,peak = {},{}
            mode_windows=windows.get(plane+'|'+selection,windows[plane])
            for code,label in LABELS.items():
                keep=selected & (stages==code)
                hist[label]=dict(slope=histogram(values[keep,i],np.linspace(-.01,.01,101)),
                                 residual=histogram(values[keep,i]-truth[keep,i],np.linspace(-.005,.005,201)))
                peak[label]=[]
                for w in mode_windows:
                    valid=keep & np.isfinite(values[:,i])
                    in_peak=valid & (np.abs(values[:,i]-w['center'])<w['half_width'])
                    peak[label].append(event_ratio_interval(entry,in_peak,valid,(0,offset)))
            result['planes'][plane+'|'+selection]=dict(histograms=hist,peak_fractions=peak,windows=mode_windows)
            fig,axes=plt.subplots(2,2,figsize=(12,8))
            for code in [1,2,5,6,7,10,11]:
                label=LABELS[code];h=hist[label]['slope'];r=hist[label]['residual']
                axes[0,0].stairs(h['counts'],h['edges'],label=label)
                axes[0,1].stairs(r['counts'],r['edges'],label=label)
            for k,w in enumerate(mode_windows):
                estimate=[peak[label][k]['estimate'] for label in LABELS.values()]
                axes[1,0].plot(list(LABELS),estimate,marker='o',label=f"Window {w['center']:.6f}")
            for code in [3,4,5,8,9,10]:
                label=LABELS[code];h=hist[label]['residual'];axes[1,1].stairs(h['counts'],h['edges'],label=label)
            axes[0,0].set(xlabel='Physical sy',ylabel='Entries / bin')
            for ax in [axes[0,1],axes[1,1]]:ax.set(xlabel='Reco sy - same-plane truth sy',ylabel='Entries / bin')
            axes[1,0].set(xlabel='Native stage code',ylabel='Frozen-window fraction')
            for ax in axes.flat:ax.legend(fontsize=7);ax.grid(alpha=.2)
            fig.suptitle(f'Native Calypso stages | {plane} | {selection} | {offset:,} source events\nFixed final-track p selection; missing stages and invalid projections are retained in counts')
            save(fig,args.output,f'{plane}_{selection}_stages',report)
    # The matched controls retain charge, module and hit multiplicity, then
    # require close truth momentum and transverse truth position.
    for i,plane in enumerate(['st1','VetoNu']):
        for wi,w in enumerate(windows.get(plane+'|reco_p_window',windows[plane])):
            available=[o for o in observations if o['entry']>=250000 and 50<o['final_p_GeV']<250 and np.isfinite(o['final_slopes'][i])]
            peaks=[o for o in available if abs(o['final_slopes'][i]-w['center'])<w['half_width']]
            side=[o for o in available if w['half_width']<=abs(o['final_slopes'][i]-w['center'])<3*w['half_width']]
            pairs=[]
            for p in peaks:
                eligible=[c for c in side if c['entry']!=p['entry'] and c['charge']==p['charge'] and c['module_eta']==p['module_eta'] and c['module_phi']==p['module_phi'] and abs(c['hits']-p['hits'])<=2 and abs(np.log(c['truth_p_GeV']/p['truth_p_GeV']))<.2 and np.hypot(c['true_x_mm']-p['true_x_mm'],c['true_y_mm']-p['true_y_mm'])<30]
                if eligible:
                    c=min(eligible,key=lambda c:(np.log(c['truth_p_GeV']/p['truth_p_GeV'])/.1)**2+((c['true_x_mm']-p['true_x_mm'])/10)**2+((c['true_y_mm']-p['true_y_mm'])/10)**2)
                    pairs.append(dict(peak=[p[k] for k in ['entry','seed','candidate']],control=[c[k] for k in ['entry','seed','candidate']]))
            result['matched_controls'].append(dict(plane=plane,window=wi,peak_tracks=len(peaks),sideband_tracks=len(side),matched_peaks=len(pairs),pairs=pairs,
                matching='Same charge and first module; hit multiplicity within two; truth p within exp(0.2); true sensor separation below 30 mm; controls may be reused'))
    for group,edges in [('archived_full_slab_X0',[0,.2,1,10,200]),('seed_stations',[.5,1.5,2.5,3.5]),
                        ('refit',[-.5,.5,1.5]),('shared_hit_fraction',[-.01,.001,.1,1.01])]:
        fig,axes=plt.subplots(1,2,figsize=(11,4.5))
        result['candidate_groups'][group]=[]
        for lo,hi in zip(edges[:-1],edges[1:]):
            selected=[o for o in observations if 50<o['final_p_GeV']<250 and lo<=o['features'][group]<hi]
            counts={}
            for i,plane in enumerate(['st1','VetoNu']):
                h=histogram(np.array([o['final_slopes'][i] for o in selected]),np.linspace(-.01,.01,101))
                counts[plane]=h
                axes[i].stairs(h['counts'],h['edges'],label=f'{lo:g} to {hi:g}; N={len(selected)}')
                axes[i].set(title=plane,xlabel='Physical final sy',ylabel='Entries / bin')
            result['candidate_groups'][group].append(dict(lower=lo,upper=hi,tracks=len(selected),histograms=counts))
        for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.2)
        fig.suptitle(f'Native candidate contributions by {group} | final p = 50-250 GeV\nArchived slab X0 sums cover recorded surfaces only; source events = {offset:,}')
        save(fig,args.output,'candidate_'+group,report)
    report.close()
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(args.output/'summary.json')


if __name__=='__main__':main()
