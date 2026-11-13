#!/usr/bin/env python3
"""Run actual Calypso audit and require exact agreement with the source shard."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from validate_charge_audit import compare_tables, inspect_audit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    task = json.loads(args.task.read_text())
    subprocess.run([sys.executable, str(Path(__file__).with_name('run_slope_task.py')),
                    '--task', str(args.task)], check=True)
    output = Path(task['output'])
    result = dict(original=task['original_output'], source_entry=task['skip'],
                  sample_entry=task['sample_entry'], targets=task['replay_targets'])
    try:
        checked = compare_tables(output, Path(task['original_output']), task['skip'])
        result.update(all_reconstruction_arrays_identical=True, fields_exact=len(checked),
                      **inspect_audit(output))
    except AssertionError as error:
        result.update(all_reconstruction_arrays_identical=False, error=str(error))
    (output/'replay_validation.json').write_text(json.dumps(result, indent=2)+'\n')
    assert result['all_reconstruction_arrays_identical'], result


if __name__ == '__main__':
    main()
