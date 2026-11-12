#!/usr/bin/env python3
"""Refine T1 using pre-finding seed eligibility, retaining failed/pruned truth."""
import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path

from momentum_simulation import t1_selection


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--tasks',type=Path,required=True)
    parser.add_argument('--particles',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();groups=defaultdict(list)
    with gzip.open(args.particles,'rt') as stream:
        for line in stream:
            row=json.loads(line)
            if row['mode'] in ['T0','T0_current']:groups[(row['grid'],row['mode'])].append(row)
    result={}
    for task in json.loads(args.tasks.read_text()):
        with gzip.open(Path(task['output'])/'injections.jsonl.gz','rt') as stream:
            injections={(r['event_ordinal'],r['seed']):r for r in map(json.loads,stream)}
        result[task['shard']]={label:t1_selection(groups[(task['shard'],mode)],injections)
                              for label,mode in [('T1','T0'),('T1_current','T0_current')]}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
