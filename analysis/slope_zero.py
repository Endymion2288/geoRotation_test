#!/usr/bin/env python3
"""A beam-centred, predefined window audit on the complete native S05 sample."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_analysis import derive
from slope_tables import merge_tables, histogram, write_json
from slope_statistics import conditional_response
from slope_controls import save


def event_interval(entries, numerator, denominator, total):
    """Ratio delta interval with source events as the independent observations."""
    x = np.bincount(entries, weights=np.asarray(numerator, float), minlength=total)
    y = np.bincount(entries, weights=np.asarray(denominator, float), minlength=total)
    estimate = x.sum()/y.sum() if y.sum() else None
    if estimate is None:
        return dict(numerator=0, denominator=0, estimate=None, interval_95=None)
    score = x-estimate*y
    se = np.sqrt(total/(total-1)*np.sum(score*score))/y.sum()
    return dict(numerator=int(x.sum()), denominator=int(y.sum()), estimate=float(estimate),
                interval_95=[max(0.,float(estimate-1.96*se)), min(1.,float(estimate+1.96*se))],
                method='Source-event ratio delta interval')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--variant',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    tables,coverage=merge_tables(args.root,args.variant)
    d=derive(tables,args.variant,'100044')
    r,t=tables['tracks'],tables['truth'];total=coverage['source_events']
    split=json.loads((args.root/'manifest.json').read_text())['held_out_split_sample_entry']
    validation=r['sample_entry']>=split
    result=dict(variant=args.variant,source_events=total,validation_entries=[split,total],
        window_origin='Beam axis sy=0; fixed half-widths 25, 50, 100, 200 microradians, without fitting the validation sample',
        reference='SiHit chord at sensor; production-magnitude/SiHit-chord Calypso transport model at fixed planes',
        planes={},response={})
    edges=np.linspace(-.002,.002,161)
    with PdfPages(args.output/'near_zero_report.pdf') as report:
        for i,plane in enumerate(['st1','VetoNu','sensor']):
            fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True)
            for ax,selection in zip(axes.flat,['all','reco_p_window','truth_p_window','matched_pure']):
                cohort=d['selections'][selection]&d['primary']&d['pure']&validation
                value=d['reconstructed'][:,i];truth=d['truth_slopes'][d['index'],i]
                finite=cohort&np.isfinite(value)&np.isfinite(truth)
                rec=histogram(value[finite],edges);tru=histogram(truth[finite],edges)
                ax.stairs(rec['counts'],edges,label='Native final reco')
                ax.stairs(tru['counts'],edges,label='Same selected particles: truth',linestyle='--')
                ax.set(title=selection,xlabel='Physical sy',ylabel='Tracks / 25 microradians',xlim=(-.002,.002));ax.legend(fontsize=8);ax.grid(alpha=.2)
                item=dict(reco=rec,same_cohort_truth=tru,windows=[])
                for half in [25e-6,50e-6,100e-6,200e-6]:
                    centre=finite&(np.abs(value)<half)
                    side=finite&(np.abs(value)>=half)&(np.abs(value)<3*half)
                    local=np.bincount(r['sample_entry'],weights=centre.astype(float)-.5*side,minlength=total)
                    excess=float(local.sum())
                    variance=total/(total-1)*(float(local@local)-excess*excess/total)
                    se=float(np.sqrt(max(0.,variance)))
                    groups={}
                    for name,group in [('momentum_tail',d['momentum_tail']),('wrong_charge',d['wrong_charge']),
                                       ('truth_p_above_250',d['prod_p']>=250),('truth_p_below_50',d['prod_p']<=50),
                                       ('truth_within_same_zero_window',np.abs(truth)<half)]:
                        groups[name]=event_interval(r['sample_entry'],centre&group,centre,total)
                    item['windows'].append(dict(half_width=half,
                        reco_fraction=event_interval(r['sample_entry'],centre,finite,total),
                        truth_fraction=event_interval(r['sample_entry'],finite&(np.abs(truth)<half),finite,total),
                        sideband_tracks=int(side.sum()),linear_sideband_excess=excess,
                        linear_sideband_excess_interval_95=[excess-1.96*se,excess+1.96*se],groups=groups))
                result['planes'][plane+'|'+selection]=item
            fig.suptitle(f'{args.variant} | {plane} | DSID 100044 | {total:,} source events\nValidation >= {split:,}; fixed beam-centred windows; same primary truth and purity >= 0.7')
            save(fig,args.output,plane+'_near_zero',report)
        # The old discovery windows were defined on afterGeo and omitted the
        # beam-centred beforeGeo structure. This additional physical window is
        # independent of both histograms; old windows and results stay frozen.
        response=conditional_response(d['truth_slopes'][:,0],np.linalg.norm(t['production_momentum_F_GeV'],axis=1),
            t['truth_sensor_position_F_mm'][:,0],t['truth_sensor_position_F_mm'][:,1],t['sample_entry'],
            d['reconstructed'][:,0],d['index'],d['primary']&d['selections']['reco_p_window']&d['pure'],
            r['sample_entry'],split,np.linspace(-.01,.01,101),
            windows=[dict(center=0.,half_width=.0001)],truth_selected=d['signal_truth'],phase_period=.08/63.)
        result['response']['st1|reco_p_window']=response
    write_json(args.output/'summary.json',result)
    print(args.output/'summary.json')


if __name__=='__main__':main()
