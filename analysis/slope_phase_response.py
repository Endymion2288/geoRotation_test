#!/usr/bin/env python3
"""Held-out Calypso response conditioned on the measured strip/lever-arm phase."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_tables import merge_tables,write_json
from slope_analysis import derive,EDGES,PLANES
from slope_statistics import conditional_response
from slope_controls import save


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--variant',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    tables,coverage=merge_tables(args.root,args.variant)
    d=derive(tables,args.variant,'100044')
    truth,tracks=tables['truth'],tables['tracks']
    split=json.loads((args.root/'manifest.json').read_text())['held_out_split_sample_entry']
    windows=json.loads((args.reference/'peak_windows.json').read_text())['windows']
    coarse=json.loads((args.reference/(args.variant+'_summary.json')).read_text())['response']
    result=dict(source_events=coverage['source_events'],variant=args.variant,
                period_definition='Actual geometry strip pitch 0.08 mm / three-layer lever arm 63 mm',
                phase_period=.08/63.,responses={})
    with PdfPages(args.output/'phase_response.pdf') as report:
        for i,plane in enumerate(PLANES[:2]):
            key=plane+'|reco_p_window'
            response=conditional_response(d['truth_slopes'][:,i],np.linalg.norm(truth['production_momentum_F_GeV'],axis=1),
                truth['truth_sensor_position_F_mm'][:,0],truth['truth_sensor_position_F_mm'][:,1],truth['sample_entry'],
                d['reconstructed'][:,i],d['index'],d['primary']&d['selections']['reco_p_window']&d['pure'],
                tracks['sample_entry'],split,EDGES,windows=windows[key],truth_selected=d['signal_truth'],phase_period=.08/63.)
            result['responses'][key]=response
            fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
            edges=np.array(response['edges']);observed=np.array(response['observed_counts']);prediction=np.array(response['predicted_counts'])
            axes[0].stairs(observed,edges,label='Held-out native reconstruction',color='#1565c0')
            axes[0].stairs(coarse[key]['predicted_counts'],edges,label='Coarse response',color='#999999',linestyle=':')
            axes[0].stairs(prediction,edges,label='Geometry-phase response',color='#222222',linestyle='--')
            low,high=np.array(response['predicted_counts_interval_95'])
            axes[0].fill_between(edges,np.r_[low,low[-1]],np.r_[high,high[-1]],step='post',alpha=.2,color='#444444')
            axes[1].stairs(observed-prediction,edges,label='Observed - phase prediction')
            low,high=np.array(response['observed_minus_predicted_interval_95'])
            axes[1].fill_between(edges,np.r_[low,low[-1]],np.r_[high,high[-1]],step='post',alpha=.2,color='#444444',label='Event bootstrap 95%')
            axes[1].axhline(0,color='k',lw=.6)
            for w in windows[key]:
                axes[0].axvspan(w['center']-w['half_width'],w['center']+w['half_width'],alpha=.07,color='green')
            for ax in axes:ax.set(xlim=(-.006,.004),ylabel='Entries / bin');ax.legend(fontsize=8);ax.grid(alpha=.2)
            axes[1].set_xlabel('Physical sy')
            fig.suptitle(f"{args.variant} | {plane} | p = 50-250 GeV | {coverage['source_events']:,} source events\nDiscovery <{split:,}; validation >= {split:,}; geometry phase fixed at 1.26984 mrad")
            save(fig,args.output,plane+'_phase_response',report)
    write_json(args.output/'summary.json',result)
    print(args.output/'summary.json')


if __name__=='__main__':main()
