#!/usr/bin/env python3
"""Freeze single-component native logger replays of incomplete observation paths."""
import argparse
import json
from pathlib import Path
import shutil

from prepare_cluster_measurements import dump


CASES = {
    'complete': (0, 0),
    'missing_early_hits': (0, 16),
    'step_limit': (1, 4),
    'terminal_only': (17, 0),
}


def prepare(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'navigation_log_v1'
    output.mkdir(exist_ok=False)
    (output/'tasks').mkdir()
    code = output/'code'
    for name in ('configs/slope_tracking.py', 'scripts/run_slope_task.py'):
        target = code/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/name, target)
    points, tasks = [], {v:[] for v in ('beforeGeo', 'afterGeo')}
    for variant in tasks:
        reference = root/'runs/100044'/variant/'fixed_observation'
        assert (reference/'exit_code.txt').read_text().strip() == '0'
        original = json.loads((reference/'source_task.json').read_text())
        for name, (entry, component) in CASES.items():
            source = Path(original['component_file'])
            lines = [line for line in source.read_text().splitlines()
                     if (int(line.split()[0]), int(line.split()[2])) == (entry, component)]
            assert len(lines) == 1
            inputs = output/'inputs'/variant/(name+'.txt')
            inputs.parent.mkdir(parents=True, exist_ok=True)
            inputs.write_text(lines[0]+'\n')
            options = original['options'].copy()
            options[options.index('--seed-component-input-file')+1] = str(inputs)
            options += ['--ckf-acts-logging', 'VERBOSE']
            folder = output/'runs'/name/variant
            task = dict(original, runner=str(code/'configs/slope_tracking.py'),
                component_file=str(inputs), options=options, output=str(folder),
                reference=str(reference), navigation_case=name,
                navigation_source_entry=entry, navigation_component=component)
            path = output/'tasks'/(name+'_'+variant+'.json')
            dump(path, task)
            tasks[variant].append(dict(kind='analysis', command=['python3',
                str(code/'scripts/run_slope_task.py'), '--task', str(path)], output=str(folder)))
            points.append(dict(case=name, variant=variant, entry=entry,
                component=component, output=str(folder), reference=str(reference)))
    for variant, jobs in tasks.items():
        dump(output/(variant+'_tasks.json'), jobs)
    dump(output/'manifest.json', dict(points=points, dsid='100044', source_events_per_replay=20,
        distinct_source_events=20, reference=str(root), isolated_smoothing=True,
        scope='Four diagnostic components in each actual native Calypso version; exact archived component '
              'inputs, frozen observation IDs and libraries. Only ACTS logger verbosity and the number '
              'of independent component calls change. Component audit equality must be verified. '
              'These selected examples do not estimate population fractions.'))


def analysis(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'navigation_log_v1'
    code = output/'analysis_code_v1'
    shutil.copytree(root/'pipeline_code_v1', code,
        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    shutil.copy2(repo/'analysis/charge_fitter_navigation.py', code/'analysis/charge_fitter_navigation.py')
    shutil.copy2(repo/'scripts/validate_charge_components.py', code/'scripts/validate_charge_components.py')
    dump(output/'analysis_tasks_v1.json', [dict(kind='analysis', command=['python3',
        str(code/'analysis/charge_fitter_navigation.py'), '--input', str(output),
        '--output', str(output/'analysis_v1')], output=str(output/'jobs/analysis_v1'))])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--phase', choices=['runs', 'analysis'], default='runs')
    args = parser.parse_args()
    (prepare if args.phase == 'runs' else analysis)(args.input.resolve())
