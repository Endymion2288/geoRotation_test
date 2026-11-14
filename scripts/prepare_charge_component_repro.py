#!/usr/bin/env python3
"""Make one-event, one-Gaussian-component native Calypso reproducers.

The input component files are indexed by local event.  A selected source
event is remapped to local event zero so that the native CKF component hook
can be run with ``SkipEvents`` while preserving the original RDO event.
"""
import argparse
import json
from pathlib import Path
import shutil


def select_rows(source, entry, component):
    rows = []
    for line in source.read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        if int(fields[0]) == entry and int(fields[2]) == component:
            fields[0] = '0'
            rows.append(' '.join(fields))
    if len(rows) != 1:
        raise RuntimeError(f'{source}: expected one row for entry={entry}, component={component}; got {len(rows)}')
    return rows[0] + '\n'


def prepare(args):
    root = args.root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    code = output / 'code'
    (code / 'scripts').mkdir(parents=True)
    (code / 'configs').mkdir()
    repo = Path(__file__).resolve().parents[1]
    for name in ('run_slope_task.py', 'prepare_cluster_measurements.py'):
        shutil.copy2(repo / 'scripts' / name, code / 'scripts' / name)
    shutil.copy2(repo / 'configs' / 'slope_tracking.py', code / 'configs' / 'slope_tracking.py')
    tasks = []
    for variant in ('beforeGeo', 'afterGeo'):
        source_task = json.loads((root / 'common_prior_components_v2' / 'smoke_v2' /
                                  '100117' / f'{variant}_source.json').read_text())
        source_components = Path(source_task['component_file'])
        for component in range(args.first, args.last + 1):
            group = output / variant / f'component_{component:03d}'
            group.mkdir(parents=True)
            component_file = group / 'component.txt'
            component_file.write_text(select_rows(source_components, args.entry, component))
            task = dict(source_task,
                        skip=args.entry, events=1,
                        runner=str(code / 'configs' / 'slope_tracking.py'),
                        output=str(group / 'result'), component_id=component,
                        component_file=str(component_file))
            options = ['--no-seed-writer', '--actual-covariance', '--isolate-smoothing',
                       '--prior-component-audit', '--seed-component-input-file', str(component_file)]
            task['options'] = options
            task_path = group / 'source.json'
            task_path.write_text(json.dumps(task, indent=2) + '\n')
            tasks.append(dict(kind='analysis', command=['python3', str(code / 'scripts' / 'run_slope_task.py'),
                                                        '--task', str(task_path)], output=task['output']))
    (output / 'tasks.json').write_text(json.dumps(tasks, indent=2) + '\n')
    for variant in ('beforeGeo', 'afterGeo'):
        selected = [task for task in tasks if f'/{variant}/' in task['output']]
        (output / f'{variant}_tasks.json').write_text(json.dumps(selected, indent=2) + '\n')
    (output / 'manifest.json').write_text(json.dumps(dict(
        dsid='100117', source_entry=args.entry, components=[args.first, args.last],
        scope='one finite Gaussian component per native Calypso CKF process; local event remapped to zero',
        tasks=len(tasks)), indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True,
                        help='stage06 isolated result root containing common_prior_components_v2')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--entry', type=int, default=3)
    parser.add_argument('--first', type=int, default=0)
    parser.add_argument('--last', type=int, default=31)
    args = parser.parse_args()
    if args.first < 0 or args.last < args.first:
        parser.error('invalid component range')
    prepare(args)


if __name__ == '__main__':
    main()
