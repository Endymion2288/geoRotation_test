#!/usr/bin/env python3
"""Compare frozen representative particles under one-factor reconstruction controls."""
import argparse
from collections import Counter,defaultdict
from copy import deepcopy
import gzip
import json
import math
from pathlib import Path

import numpy as np

from pilot import process_event,pnorm
from trace_momentum import analyze


def summarize_steps(steps,truth):
    result={};counts=Counter()
    for i,step in enumerate(steps):
        for component in ['predicted','filtered','smoothed']:
            if component in step:
                v=step[component]
                if v['p_ratio']<.5:counts[component+'_tail_steps']+=1
                if v['wrong']:counts[component+'_wrong_steps']+=1
        if 'update' in step:
            update=step['update'];delta=update['delta_qop']
            if 'largest_update' not in result or abs(delta)>abs(result['largest_update']['delta_qop']):
                result['largest_update']=dict(z_F_mm=step['x_A_mm'],state=step['state'],**update)
        for component in ['predicted','filtered']:
            if component in step and step[component]['p_ratio']<.5 and 'first_forward_tail' not in result:
                result['first_forward_tail']=dict(z_F_mm=step['x_A_mm'],state=step['state'],operation=component,p_ratio=step[component]['p_ratio'])
            if component in step and step[component]['wrong'] and 'first_forward_wrong' not in result:
                result['first_forward_wrong']=dict(z_F_mm=step['x_A_mm'],state=step['state'],operation=component,qop=step[component]['qop'])
    result['counts']=dict(counts)
    for kind in ['update','smoothing']:
        values=[s[kind]['normalized_qop_identity_error'] for s in steps if kind in s]
        result[kind+'_identity']=dict(links=len(values),nonclosing=sum(v>1e-6 for v in values),maximum=max(values,default=0.))
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    original=json.loads((args.root/'initial_tasks.json').read_text())[:30]
    variants=sorted(p.name for p in (args.root/'controls_cases').iterdir() if p.is_dir())
    results=[];summary=defaultdict(Counter)
    for task in original:
        name=Path(task['output']).name
        baseline_records=[json.loads(l) for l in gzip.open(args.root/'case_baseline'/name/'stages.jsonl.gz','rt')]
        baseline_finals={(r['seed'],r['candidate']):set(r['hits']) for r in baseline_records if r['stage']=='final'}
        paths={'baseline':args.root/'controls_smoke_v2/observer'/name}
        paths.update({v:args.root/'controls_cases'/v/name for v in variants})
        for variant,folder in paths.items():
            assert (folder/'exit_code.txt').read_text().strip()=='0',folder
            records=[json.loads(line) for line in gzip.open(folder/'stages.jsonl.gz','rt')]
            header=next(r for r in records if r['stage']=='event');key=tuple(header[k] for k in ['run','lumi','event'])
            particles,cutflow=process_event(key,records,task['dsid']);by_barcode={r['barcode']:r for r in particles}
            for case in task['cases']:
                source=case['particle'];particle=by_barcode[source['barcode']]
                result=dict(variant=variant,dsid=task['dsid'],file_index=task['file_index'],source_entry=task['skip'],
                    barcode=source['barcode'],category=case['category'],truth_p_GeV=source['truth_p_GeV'],
                    final=particle['final'],lineage=particle['final_lineage'])
                summary[variant]['particles']+=1
                final=particle['final']
                if final:
                    summary[variant]['matched']+=1
                    summary[variant]['good95']+=final['good95']
                    summary[variant]['tail_half']+=final['p_ratio']<.5
                    summary[variant]['wrong']+=final['wrong']
                    summary[variant]['winner_seed_changed']+=final['seed']!=source['final']['seed']
                    final_raw=next(r for r in records if r['stage']=='final' and r['seed']==final['seed'] and r['candidate']==final['candidate'])
                    baseline_hits=baseline_finals[(source['final']['seed'],source['final']['candidate'])]
                    result['final_hit_set_changed']=set(final_raw['hits'])!=baseline_hits
                    summary[variant]['final_hit_set_changed']+=result['final_hit_set_changed']
                    c=deepcopy(case);c['particle']['final']=final
                    traced=analyze(records,c)
                    result['operations']={stage:summarize_steps(steps,traced['truth_qop_production']) for stage,steps in traced['stages'].items()}
                results.append(result)
    (args.output/'cases.json').write_text(json.dumps(results,indent=2)+'\n')
    (args.output/'summary.json').write_text(json.dumps({k:dict(v) for k,v in summary.items()},indent=2)+'\n')
    print(json.dumps({k:dict(v) for k,v in summary.items()},indent=2))


if __name__=='__main__':main()
