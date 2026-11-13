#!/usr/bin/env python3
"""Require exact native reconstruction before retaining only refit counters."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from validate_charge_audit import compare_tables


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    task = json.loads(args.task.read_text())
    subprocess.run([sys.executable, str(Path(__file__).with_name('run_slope_task.py')),
                    '--task', str(args.task)], check=True)
    output = Path(task['output'])
    fields = compare_tables(output, Path(task['original_output']))
    (output/'counter_validation.json').write_text(json.dumps(dict(
        all_58_fields_exact=len(fields) == 58, fields=len(fields), original=task['original_output'],
        only_final_refit_counters_retained=True), indent=2)+'\n')
    assert len(fields) == 58
    for name in ['events', 'truth', 'tracks', 'segments']:
        (output/(name+'.npz')).unlink()


if __name__ == '__main__':
    main()
