#!/usr/bin/env python3
"""Compare native Calypso refit scales while holding the CKF candidates fixed."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_centroids import collect
from slope_controls import archive, save
from slope_statistics import event_ratio_interval
from slope_tables import histogram, join_truth, slope, write_json


def paired_fraction_change(entries, before, after, selected, count):
    denominator=np.bincount(entries,weights=selected,minlength=count)
    numerator=np.bincount(entries,weights=selected*(after.astype(int)-before.astype(int)),minlength=count)
    rng=np.random.default_rng(50572);values=[]
    for _ in range(1000):
        weights=rng.poisson(1.,count)
        if weights@denominator:values.append((weights@numerator)/(weights@denominator))
    return dict(estimate=float(numerator.sum()/denominator.sum()),
                interval_95=np.quantile(values,[.025,.975]).tolist(),
                definition='After minus native near-zero fraction, same candidate and same selected denominator',
                unit='source event',bootstrap_repeats=1000)


def ckf_closure(paths, baseline):
    result=[]
    def source_key(path):
        task=json.loads((path/'source_task.json').read_text())
        return task['input'],task['skip'],task['events']
    original={source_key(path):path for path in (baseline/paths[0].parent.name).glob('s*')}
    for p in paths:
        actual=archive(p/'native_states.npz')
        reference=original[source_key(p)]
        native=archive(reference/'native_states.npz')
        ia,ib=np.isin(actual['stage'],[1,3,4,5]),np.isin(native['stage'],[1,3,4,5])
        for key in ['entry','track_index','stage','values','charge','station','layer']:
            assert np.array_equal(actual[key][ia],native[key][ib],equal_nan=True),(p,key)
        result.append(dict(shard=p.name,original_shard=str(reference),seed_CKF_values_and_identity_equal=True,rows=int(ia.sum())))
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--variant',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    baseline_paths=sorted((args.baseline/'shards'/args.variant).glob('s*'))
    original,count=collect(baseline_paths)
    truth,base=original['truth'],original['tracks']
    index,matched=join_truth(original)
    pure=matched&(truth['pdg'][index]==-13)&(base['station_bits']&14==14)&(base['purity']>=.7)
    validation=base['source_entry']>=250000
    p=np.linalg.norm(base['upstream_F'][:,3:],axis=1)
    pt=np.linalg.norm(truth['production_momentum_F_GeV'][index],axis=1)
    vb=np.where(base['projection_valid'][:,0]==1,slope(base['projection'][:,0,3:6]),np.nan)
    reference=slope(truth['truth_model_momentum_F_GeV'][index,0])
    result=dict(variant=args.variant,source_events=count,actual_native_refit_scale=10,
        definition='Same original CKF candidate ordinal and majority truth; p windows and purity fixed to native scale10. Only native KalmanFitterTool1.SeedCovarianceScale changes.',
        original_baseline=str(args.baseline),scales={})
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for scale in [1,10,100]:
        paths=sorted((args.root/f'scale{scale}/shards'/args.variant).glob('s*'))
        for path in paths:assert (path/'exit_code.txt').read_text().strip()=='0'
        tables,n=collect(paths);assert n==count
        for name in ['events','truth','segments']:
            assert set(tables[name])==set(original[name])
            for key in tables[name]:
                assert np.array_equal(tables[name][key],original[name][key],equal_nan=True),(scale,name,key)
        r=tables['tracks']
        for key in ['source_entry','ordinal']:assert np.array_equal(r[key],base[key])
        closure=ckf_closure(paths,args.baseline/'shards')
        if scale==10:
            for key in r:
                equal=np.array_equal(r[key],base[key]) if r[key].dtype.kind in 'US' else np.array_equal(r[key],base[key],equal_nan=True)
                assert equal,key
        va=np.where(r['projection_valid'][:,0]==1,slope(r['projection'][:,0,3:6]),np.nan)
        same=r['hit_signatures']==base['hit_signatures']
        c=dict(seed_CKF_closure=closure,source_truth_and_segments_equal=True,
            final_rows=len(va),same_complete_hit_sets=int(same.sum()),
            majority_identity_changed=int(np.count_nonzero(r['barcode']!=base['barcode'])),cohorts={},angles={})
        for label,window in [('all',np.ones(len(p),bool)),('native_reco_p', (50<p)&(p<250)),('truth_p',(50<pt)&(pt<250))]:
            mask=pure&validation&window
            for same_hits in [False,True]:
                select=mask&same if same_hits else mask
                finite=select&np.isfinite(va)&np.isfinite(vb)
                near=np.abs(va)<.0001
                c['cohorts'][label+('|same_hits' if same_hits else '')]=dict(
                    selected_candidates=int(select.sum()),common_finite=int(finite.sum()),
                    near_zero=event_ratio_interval(r['sample_entry'],finite&near,finite,(0,count)),
                    paired_near_zero_change=paired_fraction_change(r['sample_entry'],np.abs(vb)<.0001,near,finite,count),
                    wrong_charge=event_ratio_interval(r['sample_entry'],finite&(r['charge']!=1),finite,(0,count)),
                    momentum_tail=event_ratio_interval(r['sample_entry'],finite&(np.linalg.norm(r['upstream_F'][:,3:],axis=1)<.5*pt),finite,(0,count)))
            if label=='native_reco_p':
                for ax,values,edges in [(axes[0,0],va,np.linspace(-.005,.005,251)),
                    (axes[0,1],va-reference,np.linspace(-.001,.001,201)),
                    (axes[1,0],va,np.linspace(-.0005,.0005,101))]:
                    h=histogram(values[mask],edges)
                    ax.stairs(h['counts'],edges,label=f'Refit scale {scale}'+(' (native)' if scale==10 else ''))
                c['fixed_native_p_histogram']=histogram(va[mask],np.linspace(-.005,.005,251))
        # Flag raw angle chart crossings before additional transport. A flag is
        # diagnostic evidence, not a stand-alone estimate of algorithmic bias.
        lookup={(int(e),int(o)):i for i,(e,o) in enumerate(zip(r['source_entry'],r['ordinal']))}
        flags={name:np.zeros(len(r['ordinal']),bool) for name in ['phi_outside','theta_outside','large_phi_difference']}
        for path in paths:
            a=archive(path/'native_angles.npz');m=a['writer']==1
            row=np.array([lookup[(int(e),int(o))] for e,o in zip(a['entry'],a['track_index'])])
            for name,bad in [('phi_outside',np.abs(a['ePHI_smt'])>np.pi),
                ('theta_outside',(a['eTHETA_smt']<0)|(a['eTHETA_smt']>np.pi)),
                ('large_phi_difference',np.abs(a['ePHI_smt']-a['ePHI_prt'])>np.pi)]:
                flags[name][row[m&bad]]=True
        selected=pure&validation&(50<p)&(p<250)
        for name,flag in flags.items():
            c['angles'][name]=dict(selected_tracks=int(selected.sum()),flagged=int((selected&flag).sum()),
                near_zero=int((selected&(np.abs(va)<.0001)).sum()),
                near_zero_flagged=int((selected&flag&(np.abs(va)<.0001)).sum()))
        result['scales'][str(scale)]=c
    axes[0,0].set(xlabel='Physical sy at st1',ylabel='Fixed native candidates / 40 microradians')
    axes[0,1].set(xlabel='Reco sy - same-plane truth-model sy',ylabel='Fixed native candidates / bin')
    axes[1,0].set(xlabel='Physical sy at st1 (zoom)',ylabel='Fixed native candidates / 10 microradians')
    scales=[1,10,100]
    values=[result['scales'][str(s)]['cohorts']['native_reco_p']['near_zero'] for s in scales]
    y=np.array([r['estimate'] for r in values]);ci=np.array([r['interval_95'] for r in values]).T
    axes[1,1].errorbar(scales,y,yerr=np.maximum(0,[y-ci[0],ci[1]-y]),fmt='o-')
    axes[1,1].set(xscale='log',xlabel='Native refit covariance scale (default 10)',ylabel='Fraction within |sy| < 100 microradians')
    for ax in axes.flat:
        ax.grid(alpha=.2)
        if ax is not axes[1,1]:ax.legend(fontsize=8)
    fig.suptitle(f'{args.variant} | Actual Calypso refit control | {count:,} source events\nValidation entry >= 250,000; original native pure long candidates and p = 50-250 GeV')
    with PdfPages(args.output/'refit_covariance_report.pdf') as report:save(fig,args.output,'refit_covariance_control',report)
    write_json(args.output/'summary.json',result)
    print(args.output/'summary.json')


if __name__=='__main__':main()
