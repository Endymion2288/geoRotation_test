#!/usr/bin/env python3
"""Record batch Python dependencies and check unchanged native Calypso output."""
import argparse
import importlib
import json
from pathlib import Path
import subprocess
import sys

from validate_charge_alignment import equal_tables


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    task = json.loads(args.task.read_text())
    output = Path(task['output'])
    modules = {}
    for name in ['numpy', 'awkward', 'uproot', 'scipy']:
        module = importlib.import_module(name)
        modules[name] = dict(version=module.__version__, path=module.__file__)
    record = dict(executable=sys.executable, version=sys.version, modules=modules,
                  user_site_dependency_made_explicit=True)
    (output/'python_environment.json').write_text(json.dumps(record, indent=2)+'\n')
    subprocess.run([sys.executable, str(Path(__file__).with_name('run_slope_task.py')),
                    '--task', str(args.task)], check=True)
    fields = equal_tables(output, Path(task['original_output']))
    assert len(fields) == 58
    (output/'validation.json').write_text(json.dumps(dict(all_58_fields_exact=True,
        fields=fields, python=record), indent=2)+'\n')


if __name__ == '__main__':
    main()
