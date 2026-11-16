#!/usr/bin/env python3
"""Run the latest Calypso and reduce ROOT/audit output on worker scratch."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args()
    task=json.loads(args.task.read_text())
    destination=Path(task['output'])
    destination.mkdir(parents=True,exist_ok=True)
    code=Path(__file__).resolve().parents[1]
    command=['python3',str(code/'configs/reference_tracking.py'),'--input',task['input'],
        '--events',str(task['events']),'--skip',str(task['skip']),'--material',task['material'],
        '--variant','afterGeo','--native-final-purity','--seed-unbiased','--isolate-smoothing',
        *task.get('options',[])]
    if task.get('audit',True):
        command.append('--momentum-audit')
    with Path('tracking.log').open('w') as stream:
        status=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT).returncode
    with Path('tracking.log').open('rb') as source,gzip.open(destination/'tracking.log.gz','wb') as target:
        shutil.copyfileobj(source,target)
    if Path('run.json').exists():
        shutil.copy2('run.json',destination/'run.json')
    (destination/'source_task.json').write_text(json.dumps(task,indent=2)+'\n')
    if status:
        return status
    if task.get('retain_audit'):
        with Path('momentum_audit.jsonl').open('rb') as source,gzip.open(destination/'audit.jsonl.gz','wb') as target:
            shutil.copyfileobj(source,target)
    sys.path.insert(0,str(code/'analysis'))
    from latest_momentum import reduce_run
    reduce_run(Path.cwd(),destination,task)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
