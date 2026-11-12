#!/usr/bin/env python3
"""Locate angular structures using native Calypso ROOT writers and propagation."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_controls import archive, save
from slope_tables import join_truth, slope, histogram, write_json
from slope_statistics import event_ratio_interval


LABELS={3:'CKF predicted',4:'CKF filtered',5:'CKF smoothed',8:'Refit predicted',
        9:'Refit filtered',10:'Refit smoothed',11:'Final'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--states',type=Path,required=True)
    parser.add_argument('--transport',type=Path,required=True)
    parser.add_argument('--variant',required=True)
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    rows=[];seeds=[];source_events=0;closures=[];changes=[];observer_checks=[]
    for p in sorted((args.states/'shards'/args.variant).glob('s*')):
        q=args.transport/'shards'/args.variant/p.name
        assert (p/'exit_code.txt').read_text().strip()=='0'
        assert (q/'exit_code.txt').read_text().strip()=='0'
        checked={}
        for table in ['events','truth','tracks','segments']:
            actual=archive(p/(table+'.npz'))
            original=archive(args.baseline/'shards'/args.variant/p.name/(table+'.npz'))
            assert set(actual)==set(original)
            for key in actual:assert np.array_equal(actual[key],original[key],equal_nan=True),(p.name,table,key)
            checked[table]=dict(branches=len(actual),all_saved_values_equal=True)
        observer_checks.append(dict(shard=p.name,tables=checked))
        r,t=archive(p/'tracks.npz'),archive(p/'truth.npz')
        state,projection=archive(p/'native_states.npz'),archive(q/'transport.npz')
        assert np.array_equal(state['entry'],projection['source_entry'])
        assert np.array_equal(state['stage'],projection['stage'])
        assert np.array_equal(state['track_index'],projection['seed'])
        r['sample_entry'],t['sample_entry']=r['source_entry'],t['source_entry']
        index,matched=join_truth({'tracks':r,'truth':t})
        pure=(r['station_bits']&14==14)&matched&(t['pdg'][index]==-13)&(r['purity']>=.7)
        lookup={(int(e),int(i)):j for j,(e,i) in enumerate(zip(r['source_entry'],r['ordinal']))}
        by_stage={(int(e),int(i),int(s)):k for k,(e,i,s) in enumerate(zip(state['entry'],state['track_index'],state['stage']))}
        all_entries=archive(p/'events.npz')['source_entry']
        dense={int(e):source_events+j for j,e in enumerate(all_entries)}
        source_events+=len(all_entries)
        for (entry,ordinal),j in lookup.items():
            final=by_stage[(entry,ordinal,11)]
            common=(projection['projection_valid'][final]==1)&(r['projection_valid'][j]==1)
            delta=np.abs(projection['projection_F'][final,common]-r['projection'][j,common,:6])
            closures.append(float(np.max(delta,initial=0)))
            assert closures[-1]<1e-10
            if not pure[j]:continue
            codes={code:by_stage.get((entry,ordinal,code)) for code in LABELS}
            complete=all(k is not None for k in codes.values())
            ti=index[j]
            truth=np.r_[slope(t['truth_model_momentum_F_GeV'][ti]),slope(t['truth_sensor_direction_F'][ti])]
            p_reco=np.linalg.norm(r['upstream_F'][j,3:])
            p_truth=np.linalg.norm(t['production_momentum_F_GeV'][ti])
            values={}
            for code,k in codes.items():
                if k is None:continue
                value=np.where(projection['projection_valid'][k]==1,slope(projection['projection_F'][k,:,3:]),np.nan)
                values[code]=value
                rows.append(dict(entry=entry,dense_entry=dense[entry],ordinal=ordinal,barcode=int(r['barcode'][j]),
                    stage=code,values=value,truth=truth,complete=complete,reco_p_window=50<p_reco<250,
                    truth_p_window=50<p_truth<250))
            if complete and entry>=250000 and 50<p_reco<250 and abs(values[11][0])<.0003:
                for a,b in [(3,4),(4,5),(5,10),(10,11)]:
                    changes.append(dict(entry=entry,ordinal=ordinal,barcode=int(r['barcode'][j]),
                        stages=[a,b],reco_p_GeV=float(p_reco),truth_p_GeV=float(p_truth),
                        before=values[a].tolist(),after=values[b].tolist(),
                        final_peak=bool(abs(values[11][0])<.0001)))
        truth_index={(int(e),int(b)):i for i,(e,b) in enumerate(zip(t['source_entry'],t['barcode']))}
        for k in np.flatnonzero(state['stage']==1):
            ti=truth_index.get((int(state['entry'][k]),int(state['barcode'][k])))
            if ti is None or t['pdg'][ti]!=-13:continue
            value=np.where(projection['projection_valid'][k]==1,slope(projection['projection_F'][k,:,3:]),np.nan)
            seeds.append(dict(entry=int(state['entry'][k]),values=value,
                truth=np.r_[slope(t['truth_model_momentum_F_GeV'][ti]),slope(t['truth_sensor_direction_F'][ti])],
                truth_p_GeV=float(np.linalg.norm(t['production_momentum_F_GeV'][ti]))))
    result=dict(variant=args.variant,source_events=source_events,final_projection_max_absolute_closure=max(closures),
        selected_tracks=sum(r['stage']==11 for r in rows),stage_rows=len(rows),stage_definitions=LABELS,
        pairing='Same final track and complete CKF/refit stages; selected by final momentum, not intermediate momentum',
        seed_definition='Separate matched-primary population, not associated to final ancestor seeds',
        native_writer_observer_closure=observer_checks,planes={},peak_transitions=changes)
    stage=np.array([r['stage'] for r in rows]);entry=np.array([r['dense_entry'] for r in rows])
    validation=np.array([r['entry']>=250000 for r in rows]);complete=np.array([r['complete'] for r in rows])
    value=np.array([r['values'] for r in rows]);truth=np.array([r['truth'] for r in rows])
    result['availability']={}
    for selection in ['all','reco_p_window','truth_p_window']:
        select=np.ones(len(rows),bool) if selection=='all' else np.array([r[selection] for r in rows])
        final=(stage==11)&select&validation&np.isfinite(value[:,0])
        near=final&(np.abs(value[:,0])<.0001)
        result['availability'][selection]=dict(final_tracks=int(final.sum()),
            complete_tracks=int((final&complete).sum()),final_near_zero_tracks=int(near.sum()),
            complete_near_zero_tracks=int((near&complete).sum()),
            per_stage_rows={LABELS[c]:int(((stage==c)&select&validation).sum()) for c in LABELS})
    edges=np.linspace(-.003,.003,241)
    with PdfPages(args.output/'native_stage_report.pdf') as report:
        for i,plane in enumerate(['st1','VetoNu','sensor']):
            for selection in ['all','reco_p_window','truth_p_window']:
                select=np.ones(len(rows),bool) if selection=='all' else np.array([r[selection] for r in rows])
                fig,axes=plt.subplots(2,2,figsize=(12,8))
                item={}
                for code,label in LABELS.items():
                    mask=(stage==code)&complete&select&validation
                    h=histogram(value[mask,i],edges);res=histogram(value[mask,i]-truth[mask,i],edges)
                    valid=mask&np.isfinite(value[:,i])
                    near=valid&(np.abs(value[:,i])<.0001)
                    fraction=event_ratio_interval(entry,near,valid,(0,source_events))
                    item[label]=dict(slope=h,residual=res,near_zero_fraction=fraction)
                    axes[0,0].stairs(h['counts'],edges,label=label)
                    axes[0,1].stairs(res['counts'],edges,label=label)
                    if code in [3,4,5,8,9,10]:axes[1,1].stairs(h['counts'],edges,label=label)
                axes[0,0].set(xlabel='Physical sy',ylabel='Tracks / 25 microradians')
                axes[0,1].set(xlabel='Reco sy - same-plane truth sy',ylabel='Tracks / bin')
                estimates=np.array([item[l]['near_zero_fraction']['estimate'] for l in LABELS.values()])
                ci=np.array([item[l]['near_zero_fraction']['interval_95'] for l in LABELS.values()]).T
                axes[1,0].errorbar(list(LABELS),estimates,yerr=np.maximum(0,np.array([estimates-ci[0],ci[1]-estimates])),fmt='o-')
                axes[1,0].set(xlabel='Native stage code',ylabel='Fraction within |sy| < 100 microradians')
                axes[1,1].set(xlabel='Physical sy',ylabel='Tracks / bin',xlim=(-.0005,.0005))
                for ax in axes.flat:
                    ax.grid(alpha=.2)
                    if ax is not axes[1,0]:ax.legend(fontsize=7)
                fig.suptitle(f'{args.variant} | {plane} | {selection} | {source_events:,} source events\nValidation entry >= 250,000; same primary mu+ tracks, purity >= 0.7 and available CKF/refit stages')
                save(fig,args.output,plane+'_'+selection+'_native_stages',report)
                result['planes'][plane+'|'+selection]=item
        fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharey=True)
        result['seed_population']={}
        for ax,(i,plane) in zip(axes,enumerate(['st1','VetoNu'])):
            selected=[s for s in seeds if s['entry']>=250000 and 50<s['truth_p_GeV']<250]
            for k,label in [('values','Seed population'),('truth','Same particles: truth model')]:
                h=histogram(np.array([s[k][i] for s in selected]),edges)
                ax.stairs(h['counts'],edges,label=label)
                result['seed_population'][plane+'|'+k]=h
            ax.set(title=plane,xlabel='Physical sy',ylabel='Seeds / bin');ax.legend();ax.grid(alpha=.2)
        fig.suptitle(f'{args.variant} | Separate seed population | Production truth p = 50-250 GeV\nValidation events; no claim of the same final ancestor seed')
        save(fig,args.output,'seed_population',report)
    # JSON has no non-finite placeholders; failed projections remain in histogram
    # invalid counts and in the compact source arrays.
    for r in result['peak_transitions']:
        for key in ['before','after']:r[key]=[x if np.isfinite(x) else None for x in r[key]]
    write_json(args.output/'summary.json',result)
    print(args.output/'summary.json')


if __name__=='__main__':main()
