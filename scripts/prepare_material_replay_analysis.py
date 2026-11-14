#!/usr/bin/env python3
"""Freeze independent small and complete material charge-audit reductions."""
import argparse
import json
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--revision', default='replay_analysis_v1')
    args = parser.parse_args()
    root = args.root.resolve()
    replay = root/'material_encoding_control_v1/replay_v1'
    assert json.loads((replay/'smoke_validation.json').read_text())['all_exact']
    destination = replay.parent/args.revision
    destination.mkdir(exist_ok=False)
    shutil.copytree(Path(__file__).resolve().parents[1]/'analysis', destination/'code/analysis',
                    ignore=shutil.ignore_patterns('__pycache__'))
    for scope in ['smoke', 'full']:
        command = ['python3', str(destination/'code/analysis/material_encoding_replay.py'),
                   '--root', str(root), '--output', str(destination/scope)]
        if scope == 'smoke':
            command += ['--jobs', str(replay/'smoke_tasks.json')]
        jobs = [dict(kind='analysis', command=command, output=str(destination/'jobs'/scope))]
        (destination/(scope+'_tasks.json')).write_text(json.dumps(jobs, indent=2)+'\n')
    print(destination)


if __name__ == '__main__':
    main()
