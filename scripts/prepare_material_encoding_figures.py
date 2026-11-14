#!/usr/bin/env python3
"""Freeze complete-sample material-control plots and report assembly."""
import argparse
import json
from pathlib import Path
import shutil


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--revision', default='figures_v1')
    args = parser.parse_args()
    root = args.root.resolve()
    control = root/'material_encoding_control_v1'
    assert json.loads((control/'full_validation.json').read_text())['all_exact']
    destination = control/args.revision
    destination.mkdir(exist_ok=False)
    code = destination/'code'
    shutil.copytree(Path(__file__).resolve().parents[1]/'analysis', code/'analysis',
                    ignore=shutil.ignore_patterns('__pycache__'))
    output = destination/'results'
    tools = root.parents[1]/'stage05_geometry_comparison/software_tools'
    for phase in ['plot', 'validate']:
        dump(destination/(phase+'_tasks.json'), [dict(kind='analysis', command=['python3',
            str(code/'analysis/material_encoding_figures.py'), '--phase', phase,
            '--input', str(args.analysis.resolve()), '--figures', str(output),
            '--pdf-tools', str(tools)], output=str(destination/'jobs'/phase))])
    dump(destination/'report_tasks.json', [dict(kind='analysis', command=['python3',
        str(code/'analysis/charge_report.py'), '--figures', str(root/'report_v5'), str(output),
        '--output', str(root/'report_v6'), '--pdf-tools', str(tools)],
        output=str(destination/'jobs/report'))])
    dump(destination/'manifest.json', dict(analysis=str(args.analysis.resolve()),
        material_manifest=str(control/'manifest.json'), figures=str(output), pages=12,
        scope='All three 500k-event material encoding controls; native Calypso reconstruction with isolated smoothing ON'))
    print(destination)


if __name__ == '__main__':
    main()
