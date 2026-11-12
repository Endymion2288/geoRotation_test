#!/usr/bin/env python3
"""Run native Calypso transport and export compact arrays from worker scratch."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import uproot


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    task = json.loads(args.task.read_text())
    output = Path(task['output'])
    command = [sys.executable, task['runner'], '--input', task['input'],
               '--material', task['material'], '--variant', task['variant'],
               '--events', str(task.get('events', 1)), '--skip', str(task.get('skip', 0))]
    if task.get('direction_scan'):
        command.append('--direction-scan')
    if task.get('states'):
        command.extend(['--states', task['states']])
    with Path('transport.log').open('w') as log:
        status = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    with Path('transport.log').open('rb') as source, gzip.open(output/'transport.log.gz', 'wb') as target:
        shutil.copyfileobj(source, target)
    if status:
        raise RuntimeError(f'Native transport failed: {status}; see {output}')
    with uproot.open('transport_probe.root') as file:
        data = file['transport'].arrays(library='np')
    compact = Path('transport.npz')
    np.savez_compressed(compact, **data)
    with np.load(compact,allow_pickle=False) as check:
        assert np.array_equal(check['source_entry'],data['source_entry'])
    shutil.copy2(compact,output/'transport.npz')
    (output/'source_task.json').write_text(json.dumps(task, indent=2)+'\n')
    print(f"Exported {len(data['source_entry'])} actual Calypso transport inputs", flush=True)


if __name__ == '__main__':
    main()
