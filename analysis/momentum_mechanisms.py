#!/usr/bin/env python3
"""Cross-tabulate operation labels and local-truth bounds on the selected cohort."""
import argparse
from collections import Counter,defaultdict
import gzip
import json
from pathlib import Path

import numpy as np


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    results={}
    for variant in ['current','strict']:
        counts=defaultdict(Counter);tables=defaultdict(Counter);groups=defaultdict(lambda:defaultdict(list));examples={}
        with gzip.open(args.root/'operations_analysis'/variant/'particles.jsonl.gz','rt') as stream:
            for line in stream:
                r=json.loads(line);d=r['dsid'];tail=r['final_p_ratio']<.5
                counts[d]['good95']+=1;counts[d]['tail']+=tail
                event=(r['file_index'],r['source_entry'])
                groups[d][event].append([1,tail,tail and r['local_class']=='guaranteed_local_tail',
                                        tail and r['local_class']=='unresolved_local_truth_interval'])
                if not tail:continue
                label=r['tail_class'];tables[d][label+'|'+r['local_class']]+=1
                op=r.get('tail_entry_operation')
                if op and variant+':'+d+':'+label not in examples:
                    qtruth=r['truth_charge']/r['truth_p_GeV']
                    examples[variant+':'+d+':'+label]=dict(dsid=d,file_index=r['file_index'],source_entry=r['source_entry'],barcode=r['barcode'],
                        operation=op,ratio_before=abs(qtruth/op['before']),ratio_after=abs(qtruth/op['after']),final_ratio=r['final_p_ratio'],
                        selected_material_X0=r['traces']['selected']['material_X0'],local_class=r['local_class'])
        statistics={};rng=np.random.default_rng(40412)
        for d,events in groups.items():
            matrix=np.stack([np.sum(rows,axis=0) for rows in events.values()]);total=matrix.sum(axis=0)
            boot=[]
            for _ in range(10):
                sample=matrix[rng.integers(0,len(matrix),size=(100,len(matrix)))].sum(axis=1)
                boot.append(sample[:,1:]/sample[:,0,None])
            ci=np.quantile(np.concatenate(boot),[.025,.975],axis=0)
            statistics[d]=dict(counts=dict(counts[d]),tail_operation_by_local_bound=dict(tables[d]),
                fractions={label:dict(value=float(total[i+1]/total[0]),event_bootstrap95=ci[:,i].tolist())
                           for i,label in enumerate(['production_tail','guaranteed_local_tail','unresolved_local_interval_tail'])})
        results[variant]=dict(samples=statistics,examples=examples)
    (args.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
