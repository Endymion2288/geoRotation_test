#!/usr/bin/env python3
"""Plot paired material validation and compare material/gain on common track paths."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read(path):
    return json.loads(path.read_text())


def core_summary(path):
    events=defaultdict(lambda:defaultdict(list))
    with gzip.open(path,'rt') as stream:
        for line in stream:
            row=json.loads(line)
            if not (row['fiducial'] and row['truth_H18']):continue
            final=row['final']
            value=final['p_ratio'] if final and final['good95'] else float('nan')
            events[row['dsid']][(row['file_index'],row['source_entry'])].append(value)
    result={}
    rng=np.random.default_rng(40126)
    for dsid,groups in events.items():
        width=max(map(len,groups.values()));matrix=np.full((len(groups),width),np.nan)
        for i,values in enumerate(groups.values()):matrix[i,:len(values)]=values
        good=matrix[np.isfinite(matrix)]
        matrix[(matrix<.8)|(matrix>1.2)]=np.nan
        core=matrix[np.isfinite(matrix)]
        bootstrap=[]
        for _ in range(16):
            replicas=matrix[rng.integers(0,len(matrix),size=(32,len(matrix)))].reshape(32,-1)
            quantiles=np.nanquantile(replicas,[.16,.5,.84],axis=1).T
            bootstrap.extend(np.stack([quantiles[:,1],(quantiles[:,2]-quantiles[:,0])*.5],axis=-1))
        q=np.quantile(core,[.16,.5,.84]);ci=np.quantile(bootstrap,[.025,.975],axis=0)
        result[dsid]=dict(good95=len(good),core_n=len(core),definition='0.8 <= p ratio <= 1.2 within Good95; event bootstrap 512 replicas',
            core_median=float(q[1]),core_half_68width=float((q[2]-q[0])*.5),
            median_bootstrap95=ci[:,0].tolist(),half_68width_bootstrap95=ci[:,1].tolist(),
            tails={str(t):int(np.sum(good<t)) for t in [.8,.5,.2]})
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--operations',action='store_true')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    r=args.root
    cohorts={
        'Discovery':(read(r.parent/'stage03_baseline/merged_current/summary.json'),read(r/'merged_consistent/summary.json'),read(r/'paired_material.json')),
        'Independent':(read(r/'independent_merged/current/summary.json'),read(r/'independent_merged/strict/summary.json'),read(r/'paired_independent.json')),
    }
    dsids=['100044','100116','100117','110033','110051'];x=np.arange(len(dsids))
    fig,axes=plt.subplots(2,3,figsize=(14,8))
    for row,(cohort,(base,alt,pair)) in enumerate(cohorts.items()):
        for ax,key,label in zip(axes[row],['tail_fraction_difference','misid_fraction_difference','good95_efficiency_difference'],
                                ['Tail fraction change [percentage points]','Charge misID change [percentage points]','Good95 efficiency change [percentage points]']):
            values=[pair[d]['rate_differences'][key]['value']*100 for d in dsids]
            intervals=np.array([pair[d]['rate_differences'][key]['event_bootstrap95'] for d in dsids]).T*100
            ax.errorbar(x,values,yerr=[np.asarray(values)-intervals[0],intervals[1]-np.asarray(values)],fmt='o',capsize=3)
            ax.axhline(0,color='grey',lw=1);ax.set_xticks(x);ax.set_xticklabels(dsids,rotation=40)
            ax.set(ylabel=label,title=cohort);ax.grid(axis='y',alpha=.3)
    fig.suptitle('Strictly rotated current22 material minus unrotated current22; event bootstrap 95% intervals')
    fig.tight_layout();fig.savefig(args.output/'paired_material_validation.png',dpi=180);plt.close(fig)
    compact={}
    for name,(a,b,pair) in cohorts.items():
        compact[name]={d:dict(current=a['samples'][d],strict=b['samples'][d],paired=pair[d]) for d in dsids}
    (args.output/'validation.json').write_text(json.dumps(compact,indent=2)+'\n')
    if not args.operations:return
    cores={variant:core_summary(r/'independent_merged'/variant/'particles.jsonl.gz') for variant in ['current','strict']}
    (args.output/'independent_core.json').write_text(json.dumps(cores,indent=2)+'\n')
    def load(variant):
        with gzip.open(r/'operations_analysis'/variant/'particles.jsonl.gz','rt') as stream:
            rows=[json.loads(line) for line in stream]
        return {(a['dsid'],a['file_index'],a['source_entry'],a['barcode']):a for a in rows}
    a,b=load('current'),load('strict')
    results={};fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    for d in dsids:
        pairs=[(a[k],b[k]) for k in a.keys()&b.keys() if k[0]==d]
        counts=Counter();values={k:[] for k in ['X0_current','X0_strict','gain_current','gain_strict','final_p_ratio_current','final_p_ratio_strict']}
        for left,right in pairs:
            same=left['selected_hit_ids']==right['selected_hit_ids'] and left['seed']==right['seed']
            counts['common_good95']+=1;counts['same_selected_hits_and_seed']+=same
            if not same:continue
            counts['tail_recovered']+=left['final_p_ratio']<.5 and right['final_p_ratio']>=.5
            counts['tail_introduced']+=left['final_p_ratio']>=.5 and right['final_p_ratio']<.5
            for name,item in [('current',left),('strict',right)]:
                values['X0_'+name].append(item['traces']['selected']['material_X0'])
                values['gain_'+name].append(item['traces']['selected']['max_gain_qop'])
                values['final_p_ratio_'+name].append(item['final_p_ratio'])
        results[d]=dict(counts=dict(counts),same_path_quantiles={k:np.quantile(v,[0,.16,.5,.84,1]).tolist() if v else [] for k,v in values.items()})
        if values['X0_current']:
            axes[0].scatter(values['X0_current'],values['X0_strict'],s=3,alpha=.15,label=d,rasterized=True)
            axes[1].scatter(values['final_p_ratio_current'],values['final_p_ratio_strict'],s=3,alpha=.15,label=d,rasterized=True)
    axes[0].set(xlabel='Selected path X/X0, unrotated map',ylabel='Selected path X/X0, strictly rotated map',xscale='log',yscale='log')
    axes[1].set(xlabel='Final p ratio, unrotated map',ylabel='Final p ratio, strictly rotated map',xlim=(0,2),ylim=(0,2))
    axes[1].plot([0,2],[0,2],color='grey',lw=1)
    for ax in axes:ax.legend(markerscale=3)
    fig.suptitle('Common Good95 tracks with identical selected hit sequence and seed')
    fig.tight_layout();fig.savefig(args.output/'same_path_material.png',dpi=180);plt.close(fig)
    (args.output/'paired_operations.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
