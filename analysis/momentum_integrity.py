#!/usr/bin/env python3
"""Compare observer replays field by field, including selection and fallback."""
import argparse
from collections import Counter
import json
from pathlib import Path

from pilot import archive_records, event_groups, read_records


def good95(row):
    return bool(row.get('ndof',0)>0 and row.get('chi2') is not None
                and row['chi2']/row['ndof']<25 and row.get('layers',0)>=7
                and row['ndof']+5>=12 and row.get('rmax_mm',float('inf'))<95)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--nominal',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    archive=(args.nominal/'index.json').exists()
    source=json.loads((args.nominal/('index.json' if archive else 'tasks.json')).read_text())
    tasks=source['tasks'] if archive else source
    index={(t['dsid'],t['file_index'],t['skip']):t for t in tasks}
    stages={'seed','found','selected','converted','refit_input','refit','final'}
    counts=Counter();differences=Counter();examples={}
    fields=['seed','candidate','momentum_A_GeV','position_A_mm','charge','bound','bound_trk',
            'variance_qop','hits','matches','chi2','ndof','measurements','layers','rmax_mm','refit']
    for task in json.loads((args.input/'tasks.json').read_text()):
        nominal=index[(task['dsid'],task['file_index'],task['skip'])]
        old=archive_records(args.nominal,nominal) if archive else read_records(args.nominal/nominal['shard']/'stages.jsonl.gz')
        new=read_records(args.input/task['shard']/'stages.jsonl.gz')
        left=event_groups(old);right=event_groups(new)
        for ordinal,b in enumerate(right):
            a=next(left);counts['events']+=1
            for stage in stages:
                aa=[r for r in a if r['stage']==stage];bb=[r for r in b if r['stage']==stage]
                assert len(aa)==len(bb),(task,ordinal,stage)
                counts[stage+'_tracks']+=len(aa)
                for x,y in zip(aa,bb):
                    for field in fields:
                        if x.get(field)!=y.get(field):
                            key=stage+':'+field;differences[key]+=1
                            if key not in examples:examples[key]=dict(dsid=task['dsid'],file_index=task['file_index'],source_entry=task['skip']+ordinal,
                                seed=x['seed'],candidate=x['candidate'],nominal=x.get(field),observer=y.get(field))
                    if stage=='final':
                        counts['final_good95_changed']+=good95(x)!=good95(y)
                        counts['final_momentum_or_charge_changed']+=any(x.get(k)!=y.get(k) for k in ['momentum_A_GeV','charge'])
        assert next(left,None) is None
        print('Completed',task['dsid'],task['file_index'],task['skip'],flush=True)
    (args.output/'summary.json').write_text(json.dumps(dict(counts=dict(counts),differences=dict(differences),examples=examples),indent=2)+'\n')


if __name__=='__main__':main()
