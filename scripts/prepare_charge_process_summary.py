#!/usr/bin/env python3
"""Freeze reduction and four report pages for the resolved Calypso matrix."""
import argparse
import json
from pathlib import Path
import shutil

from prepare_cluster_measurements import dump


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    control = root/'fixed_hit_process_v1'
    resolved = control/'resolved_tasks_v1.json'
    assert len(json.loads(resolved.read_text())) == 792
    repo = Path(__file__).resolve().parents[1]
    code = control/'analysis_code_v2'
    for folder in ['scripts', 'analysis']:
        shutil.copytree(repo/folder, code/folder, ignore=shutil.ignore_patterns('__pycache__'))
    analysis = control/'analysis_v2/full'
    figures = control/'figures_v2'
    figures.mkdir()
    pdf_tools = root.parent.parent/'stage05_geometry_comparison/software_tools'
    tasks = []

    def task(name, command, destination):
        path = control/(name+'_tasks_v2.json')
        dump(path, [dict(kind='analysis', command=command, output=str(destination))])
        tasks.append(str(path))

    task('full_analysis', ['python3', str(code/'analysis/charge_process_noise.py'),
        '--root', str(root), '--jobs', str(resolved), '--output', str(analysis)],
        control/'jobs/analysis_v2_full')
    for phase in ['plot', 'validate']:
        task('figures_'+phase, ['python3', str(code/'analysis/charge_process_noise_figures.py'),
            '--phase', phase, '--input', str(analysis), '--output', str(figures/'results'),
            '--pdf-tools', str(pdf_tools)], figures/'jobs'/phase)
    task('report', ['python3', str(code/'analysis/charge_report.py'), '--figures', str(root/'report_v9'),
        str(figures/'results'), '--output', str(root/'report_v10'), '--pdf-tools', str(pdf_tools)],
        figures/'jobs/report')
    dump(control/'sequence_tasks_v2.json', [dict(kind='analysis', command=['python3',
        str(code/'scripts/run_analysis_sequence.py'), '--task-files', *tasks,
        '--worker', str(code/'scripts/momentum_task.py')], output=str(control/'jobs/sequence_v2'))])


if __name__ == '__main__':
    main()
