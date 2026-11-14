#!/usr/bin/env python3
"""Freeze mechanism figures for the complete corrected-material replay."""
import argparse
import json
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--revision', default='replay_figures_v1')
    args = parser.parse_args()
    root = args.root.resolve()
    analysis = args.analysis.resolve()
    assert json.loads((analysis/'full/summary.json').read_text())['events'] == 119
    destination = root/'material_encoding_control_v1'/args.revision
    destination.mkdir(exist_ok=False)
    shutil.copytree(Path(__file__).resolve().parents[1]/'analysis', destination/'code/analysis',
                    ignore=shutil.ignore_patterns('__pycache__'))
    output = destination/'results'
    pdf_tools = root.parents[1]/'stage05_geometry_comparison/software_tools'
    for phase in ['plot', 'validate']:
        jobs = [dict(kind='analysis', command=['python3',
            str(destination/'code/analysis/material_encoding_replay_figures.py'), '--phase', phase,
            '--input', str(analysis/'full'), '--figures', str(output), '--pdf-tools', str(pdf_tools)],
            output=str(destination/'jobs'/phase))]
        (destination/(phase+'_tasks.json')).write_text(json.dumps(jobs, indent=2)+'\n')
    jobs = [dict(kind='analysis', command=['python3', str(destination/'code/analysis/charge_report.py'),
        '--figures', str(root/'report_v6'), str(output), '--output', str(root/'report_v7'),
        '--pdf-tools', str(pdf_tools)], output=str(destination/'jobs/report'))]
    (destination/'report_tasks.json').write_text(json.dumps(jobs, indent=2)+'\n')
    print(destination)


if __name__ == '__main__':
    main()
