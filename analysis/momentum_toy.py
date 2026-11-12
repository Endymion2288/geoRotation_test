#!/usr/bin/env python3
"""Summarize actual ACTS ideal-fit closure, separating noiseless replicas."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binomtest


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    data=np.genfromtxt(args.input,names=True,delimiter=',')
    summary=[]
    for p in np.unique(data['p_truth_GeV']):
        for q in [-1,1]:
            for mode in range(3):
                rows=data[(data['p_truth_GeV']==p)&(data['charge']==q)&(data['mode']==mode)]
                noise=rows[rows['replica']>0];good=noise[noise['success']==1]
                ratios=good['p_reco_GeV']/p
                qtruth=q/p;pull=(good['qop_reco']-qtruth)/good['qop_sigma']
                n=len(good);k=int(np.sum(ratios<.5))
                interval=binomtest(k,n).proportion_ci() if n else None
                det=rows[(rows['replica']==0)&(rows['success']==1)]
                summary.append(dict(p_GeV=float(p),charge=q,mode=mode,trials=len(noise),success=n,
                    noiseless_relative_p_error=float(det['p_reco_GeV'][0]/p-1) if len(det) else None,
                    qop_bias=float(np.mean(good['qop_reco']-qtruth)),pull_mean=float(np.mean(pull)),
                    pull_std=float(np.std(pull,ddof=1)),mean_pull_standard_error=float(np.std(pull,ddof=1)/np.sqrt(n)),
                    p_ratio_quantiles=np.quantile(ratios,[.16,.5,.84]).tolist(),
                    tails={str(t):int(np.sum(ratios<t)) for t in [.8,.5,.2]},
                    wrong=int(np.sum(good['qop_reco']*q<0)),tail_half_ci95=[float(interval.low),float(interval.high)] if n else None))
    (args.output/'summary.json').write_text(json.dumps(dict(fits=len(data),groups=summary),indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for mode in range(3):
        rows=[r for r in summary if r['mode']==mode and r['charge']==1]
        x=[r['p_GeV'] for r in rows]
        axes[0].errorbar(x,[r['pull_mean'] for r in rows],yerr=[1.96*r['mean_pull_standard_error'] for r in rows],marker='o',label=f'Mode {mode}')
        axes[1].plot(x,[r['pull_std'] for r in rows],marker='o',label=f'Mode {mode}')
    for ax in axes:
        ax.set_xscale('log');ax.set_xlabel('Truth momentum [GeV]');ax.legend()
    axes[0].axhline(0,color='grey',lw=1);axes[0].set_ylabel('Mean q/p pull (95% mean interval)')
    axes[1].axhline(1,color='grey',lw=1);axes[1].set_ylabel('q/p pull standard deviation')
    fig.tight_layout();fig.savefig(args.output/'toy_pulls.png',dpi=160);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
