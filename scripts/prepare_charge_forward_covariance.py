#!/usr/bin/env python3
"""Freeze read-only forward-covariance analysis of the actual Calypso audit."""
import argparse
import json
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--revision', type=int, default=1)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.root.resolve()/'common_prior_ckf_chart_v1'
    source = root/'smoothing_operations_v1'
    output = root/'forward_covariance_v1'
    output.mkdir(exist_ok=args.revision > 1)
    revision = 'v'+str(args.revision)
    code = output/('code_'+revision)
    code.mkdir(exist_ok=False)
    paths = ['scripts/validate_charge_components.py'] + ['analysis/'+name for name in
        ['charge_forward_covariance.py', 'charge_component_refit.py', 'charge_refit_input.py',
         'charge_smoothing_operations.py', 'slope_tables.py']]
    for relative in paths:
        target = code/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, target)
    tasks = [dict(kind='analysis', command=['python3', str(code/'analysis/charge_forward_covariance.py'),
        '--input', str(source), '--output', str(output/('analysis_'+revision))], output=str(output/'jobs'/('analysis_'+revision)))]
    (output/('analysis_tasks_'+revision+'.json')).write_text(json.dumps(tasks, indent=2)+'\n')
    (output/('manifest_'+revision+'.json')).write_text(json.dumps(dict(source_manifest=str(source/'manifest.json'),
        source_validation=str(source/'analysis_v2/validation.json'), code=str(code),
        frozen_files=paths, isolated_smoothing=True,
        scope='Read-only native selected-state analysis, including passive/material surfaces; '
              'no reconstruction changes and no standalone transport.'), indent=2)+'\n')
    print(output)


if __name__ == '__main__':
    main()
