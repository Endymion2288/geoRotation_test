#!/usr/bin/env python3
"""Compare the same truth particles, resampling complete source events."""
import argparse
from collections import Counter,defaultdict
import gzip
import json
from pathlib import Path

import numpy as np


def load(path):
    with gzip.open(path,'rt') as stream:
        rows=[json.loads(line) for line in stream]
    return {(r['dsid'],r['file_index'],r['source_entry'],r['barcode']):r for r in rows}


def good(row):
    return bool(row['fiducial'] and row['truth_H18'] and row['final'] and row['final']['good95'])


def compare(a,b):
    assert a.keys()==b.keys(),'Paired truth coverage differs'
    result={}
    for dsid in sorted({k[0] for k in a}):
        pairs=[(a[k],b[k]) for k in a if k[0]==dsid]
        counts=Counter();events=defaultdict(lambda:np.zeros(9,dtype=int));common=[]
        for x,y in pairs:
            assert x['truth_p_GeV']==y['truth_p_GeV'] and x['fiducial']==y['fiducial'] and x['truth_H18']==y['truth_H18']
            if not x['fiducial'] or not x['truth_H18']:continue
            gx,gy=good(x),good(y);fx,fy=x['final'],y['final']
            tx=bool(gx and fx['p_ratio']<.5);ty=bool(gy and fy['p_ratio']<.5)
            wx=bool(gx and fx['wrong']);wy=bool(gy and fy['wrong'])
            events[(x['file_index'],x['source_entry'])]+=[1,gx,gy,tx,ty,wx,wy,gx and not wx,gy and not wy]
            counts[f'good_{int(gx)}{int(gy)}']+=1
            if gx and gy:
                common.append((x,y));counts[f'tail_{int(tx)}{int(ty)}']+=1;counts[f'wrong_{int(wx)}{int(wy)}']+=1
                counts['winner_seed_changed']+=fx['seed']!=fy['seed']
                counts['winner_candidate_changed']+=fx['candidate']!=fy['candidate']
        matrix=np.stack(list(events.values()));total=matrix.sum(axis=0)
        def rates(t):
            with np.errstate(divide='ignore',invalid='ignore'):
                return np.stack([t[...,2]/t[...,0]-t[...,1]/t[...,0],
                    t[...,4]/t[...,2]-t[...,3]/t[...,1],t[...,6]/t[...,2]-t[...,5]/t[...,1],
                    (t[...,8]-t[...,7])/t[...,0],(t[...,4]-t[...,3])/t[...,0]],axis=-1)
        rng=np.random.default_rng(40041);samples=[]
        for _ in range(10):
            boot=matrix[rng.integers(0,len(matrix),size=(100,len(matrix)))].sum(axis=1);samples.append(rates(boot))
        ci=np.nanquantile(np.concatenate(samples),[.025,.975],axis=0)
        labels=['good95_efficiency_difference','tail_fraction_difference','misid_fraction_difference','correct_efficiency_difference','tail_yield_per_truth_difference']
        stats={label:dict(value=float(v),event_bootstrap95=ci[:,i].tolist()) for i,(label,v) in enumerate(zip(labels,rates(total)))}
        distributions={}
        for label,index in [('baseline',0),('alternative',1)]:
            selected=[p[index] for p in pairs if good(p[index])]
            values=np.array([r['final']['p_ratio'] for r in selected])
            slopes=np.array([r['final']['slope_y']-r['truth_slope_y'] for r in selected])
            distributions[label]=dict(n=len(selected),p_ratio_quantiles=np.quantile(values,[.16,.5,.84]).tolist(),
                slope_residual_quantiles=np.quantile(slopes,[.16,.5,.84]).tolist(),
                tails={str(t):int(np.sum(values<t)) for t in [.8,.5,.2]})
        result[dsid]=dict(counts=dict(counts),paired_events=len(matrix),fiducial_H18_truth=int(total[0]),
            rate_differences=stats,distributions=distributions,common_good95=len(common),
            common_p_ratio_shift_quantiles=np.quantile([y['final']['p_ratio']-x['final']['p_ratio'] for x,y in common],[.16,.5,.84]).tolist() if common else [])
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--baseline',type=Path,required=True);parser.add_argument('--alternative',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=compare(load(args.baseline),load(args.alternative))
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
