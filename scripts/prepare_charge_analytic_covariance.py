#!/usr/bin/env python3
"""Replace only frozen afterGeo seed covariance by the analytic coordinate map."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np

from prepare_charge_component_tasks import _jacobian_to_after
from validate_charge_components import input_rows


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def prepare(root):
    repo = Path(__file__).resolve().parents[1]
    chart = root/'common_prior_ckf_chart_v1'
    reference = chart/'width_scan_v1'
    baseline = json.loads((reference/'manifest.json').read_text())
    output = chart/'analytic_covariance_v1'
    output.mkdir(exist_ok=False)
    code = output/'run_code'
    for relative in ['configs/slope_tracking.py', 'scripts/run_slope_task.py',
                     'scripts/prepare_charge_component_tasks.py',
                     'scripts/prepare_charge_analytic_covariance.py',
                     'scripts/prepare_cluster_measurements.py',
                     'scripts/validate_charge_components.py']:
        target = code/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, target)
    points, tasks, checks = [], [], []
    for point in baseline['points']:
        series, dsid = point['series'], point['dsid']
        before_folder, old_folder = Path(point['beforeGeo']), Path(point['afterGeo'])
        assert all((p/'exit_code.txt').read_text().strip() == '0' for p in (before_folder, old_folder))
        before_task = json.loads((before_folder/'source_task.json').read_text())
        old_task = json.loads((old_folder/'source_task.json').read_text())
        before_path, old_path = Path(before_task['component_file']), Path(old_task['component_file'])
        before, after = input_rows(before_path), input_rows(old_path)
        assert before.keys() == after.keys()
        lines, changed, maximum = [], 0, 0.
        for line in old_path.read_text().splitlines():
            if not line.strip():
                continue
            fields = line.split()
            key = int(fields[0]), int(fields[2])
            a, b = before[key], after[key]
            assert a['seed'] == b['seed']
            # Existing ACTS construction canonicalizes negative theta without
            # rotating C. Frozen before inputs are diagonal, so equivalent
            # angle representations preserve their physical local covariance.
            assert np.array_equal(a['cov'], np.diag(np.diag(a['cov'])))
            jacobian = _jacobian_to_after(a['bound'], a['surface'], b['surface'], method='analytic')
            covariance = jacobian @ a['cov'] @ jacobian.T
            covariance = (covariance+covariance.T)/2.
            np.linalg.cholesky(covariance)
            changed += not np.array_equal(covariance, b['cov'])
            maximum = max(maximum, float(np.max(np.abs(covariance-b['cov']))))
            fields[9:45] = [format(v, '.17g') for v in covariance.reshape(-1)]
            lines.append(' '.join(fields)+'\n')
        component_file = output/'inputs'/series/dsid/'afterGeo_components.txt'
        component_file.parent.mkdir(parents=True)
        component_file.write_text(''.join(lines))
        new_rows = input_rows(component_file)
        for key, old in after.items():
            new = new_rows[key]
            assert old['seed'] == new['seed']
            for field in ['bound', 'surface']:
                np.testing.assert_array_equal(old[field], new[field])
        # Also check all non-C tokens, including event/seed/component identity.
        keep = lambda s: [f[:9]+f[45:] for f in map(str.split, s.splitlines()) if f]
        assert keep(component_file.read_text()) == keep(old_path.read_text())
        options = old_task['options'].copy()
        index = options.index('--seed-component-input-file')+1
        assert options[index] == str(old_path)
        options[index] = str(component_file)
        folder = output/'runs'/series/dsid/'afterGeo_analytic'
        task = dict(old_task, runner=str(code/'configs/slope_tracking.py'), options=options,
                    output=str(folder), component_file=str(component_file),
                    covariance_jacobian='analytic', reference_covariance_task=str(old_folder/'source_task.json'))
        task_path = output/'tasks'/(series+'_'+dsid+'_afterGeo.json')
        write(task_path, task)
        tasks.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                                                   '--task', str(task_path)], output=str(folder)))
        points.append(dict(series=series, dsid=dsid, beforeGeo=str(before_folder),
                           afterGeo=str(folder), afterGeo_reference=str(old_folder)))
        checks.append(dict(series=series, dsid=dsid, components=len(after),
            changed_covariances=changed, maximum_absolute_C_change_native_units=maximum,
            non_covariance_input_tokens_exact=True, before_input_reused=str(before_path),
            after_reference_input=str(old_path), after_analytic_input=str(component_file),
            before_C_diagonal_verified=True))
    write(output/'afterGeo_tasks.json', tasks)
    write(output/'manifest.json', dict(points=points, series=baseline['series'],
        samples=baseline['samples'], before_library=baseline['before_library'],
        after_library=baseline['after_library'], isolated_smoothing=True,
        refit_diagonal_scale=1., source_events=60, checks=checks,
        scope='Actual unchanged Calypso chain; only after seed covariance differential changes '
              'from the frozen central difference to analytic. All means, surfaces, before '
              'inputs, finite component families and reconstruction settings are fixed. '
              'Not a finite Gaussian distribution transform or a mixture posterior.'))
    print(json.dumps(dict(output=str(output), new_reconstructions=len(tasks),
                          components=sum(r['components'] for r in checks)), indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve())


if __name__ == '__main__':
    main()
