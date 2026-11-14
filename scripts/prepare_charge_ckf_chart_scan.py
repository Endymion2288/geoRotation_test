#!/usr/bin/env python3
"""Freeze the existing finite-prior width scan with native regular CKF controls."""
import argparse
import json
from pathlib import Path
import shutil


SAMPLES = ('100044', '100116', '100117')
SERIES = ('matched_n32_f256_v1', 'matched_n32_f1024_v1',
          'matched_n64_f1024_v1', 'matched_n64_f4096_v1')


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def prepare(root):
    repo = Path(__file__).resolve().parents[1]
    chart = root/'common_prior_ckf_chart_v1'
    assert json.loads((chart/'gate_v1/validation.json').read_text())['all_exact']
    output = chart/'width_scan_v1'
    output.mkdir(exist_ok=False)
    code = output/'run_code'
    for relative in ('configs/slope_tracking.py', 'scripts/run_slope_task.py'):
        target = code/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, target)
    jobs, points = {'beforeGeo': [], 'afterGeo': []}, []
    for series in SERIES:
        for dsid in SAMPLES:
            point = dict(series=series, dsid=dsid)
            for variant, suffix in [('beforeGeo', 'regular_both'), ('afterGeo', 'native')]:
                name = variant+'_'+suffix
                if series == SERIES[0]:
                    folder = (chart/series/dsid/'both' if variant == 'beforeGeo' else
                              root/'common_prior_refit_input_v1/scale1_v1'/series/dsid/name)
                    assert (folder/'exit_code.txt').read_text().strip() == '0'
                else:
                    source = root/'common_prior_component_refit_v1/tasks'/(series+'_'+dsid+'_'+name+'.json')
                    task = json.loads(source.read_text())
                    assert '--refit-covariance-scale' not in task['options']
                    options = task['options']+['--refit-covariance-scale', '1.0']
                    if variant == 'beforeGeo':
                        options += ['--regular-ckf-update-chart', '--regular-ckf-smoothing-chart']
                    folder = output/'runs'/series/dsid/name
                    task.update(runner=str(code/'configs/slope_tracking.py'), options=options,
                                output=str(folder), reference_component_task=str(source),
                                ckf_chart='both' if variant == 'beforeGeo' else 'native',
                                refit_input_audit=True)
                    task_path = output/'tasks'/(series+'_'+dsid+'_'+name+'.json')
                    write(task_path, task)
                    jobs[variant].append(dict(kind='analysis', command=['python3',
                        str(code/'scripts/run_slope_task.py'), '--task', str(task_path)], output=str(folder)))
                point[variant] = str(folder)
            points.append(point)
    for variant, tasks in jobs.items():
        write(output/(variant+'_tasks.json'), tasks)
    write(output/'manifest.json', dict(points=points, series=list(SERIES), samples=list(SAMPLES),
        reused_baseline_series=SERIES[0], new_reconstructions=sum(map(len, jobs.values())),
        before_library=str(chart/'overlay_beforeGeo/lib'),
        after_library=str(root/'common_prior_component_refit_v1/overlay_afterGeo/lib'),
        isolated_smoothing=True, refit_diagonal_scale=1.,
        scope='Actual native CKF+conversion+refit; before regular update/smoothing in both fits; '
              'after native; fixed finite-prior input families, no mixture posterior aggregation.'))
    print(json.dumps(dict(output=str(output), jobs={k:len(v) for k,v in jobs.items()}), indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root)


if __name__ == '__main__':
    main()
