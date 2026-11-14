#!/usr/bin/env python3
"""Prepare finite common physical-prior CKF tasks for the actual Calypso.

Split the recorded beforeGeo Gaussian into a finite mixture whose mean and
total covariance equal the recorded input. Map every component mean exactly,
and its covariance at that component's own mean. Finite-width Gaussian chart
transforms remain an approximation; record their physical direction error so
that width and component-count convergence can be checked. All fitting runs
in the actual Calypso CKF.
"""
import argparse
import gzip
import json
from pathlib import Path
import shutil

import numpy as np

from prepare_cluster_measurements import dump


ROTATE_F_TO_A = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
SAMPLES = ('100044', '100116', '100117')
VARIANTS = ('beforeGeo', 'afterGeo')


def _read_seed_rows(path):
    rows = {}
    current = None
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            row = json.loads(line)
            if row['stage'] == 'event':
                current = int(row.get('local_event_index', row['event']))
            elif row['stage'] == 'seed':
                assert current is not None
                rows[(current, int(row['seed']))] = row
    return rows


def _direction(phi, theta):
    st = np.sin(theta)
    return np.array([np.cos(phi) * st, np.sin(phi) * st, np.cos(theta)])


def _angle_delta(value, reference):
    return np.arctan2(np.sin(value - reference), np.cos(value - reference))


def _bound_to_physical(mean, surface, frame):
    transform = np.asarray(surface, dtype=float).reshape(4, 4)
    local = np.array([mean[0], mean[1], 0., 1.])
    native_position = (transform @ local)[:3]
    native_direction = _direction(float(mean[2]), float(mean[3]))
    if frame == 'FASER':
        position_f, direction_f = native_position, native_direction
    else:
        position_f, direction_f = ROTATE_F_TO_A.T @ native_position, ROTATE_F_TO_A.T @ native_direction
    return position_f, direction_f, float(mean[4]), float(mean[5])


def _physical_to_bound(position_f, direction_f, qop, time, surface, frame):
    transform = np.asarray(surface, dtype=float).reshape(4, 4)
    if frame == 'FASER':
        native_position, native_direction = position_f, direction_f
    else:
        native_position, native_direction = ROTATE_F_TO_A @ position_f, ROTATE_F_TO_A @ direction_f
    local = np.linalg.solve(transform, np.r_[native_position, 1.])
    direction = np.asarray(native_direction, dtype=float)
    direction /= np.linalg.norm(direction)
    phi = np.arctan2(direction[1], direction[0])
    theta = np.arctan2(np.hypot(direction[0], direction[1]), direction[2])
    return np.array([local[0], local[1], phi, theta, qop, time], dtype=float)


def _physical_map(mean_f, surface_f, surface_a, qop=None):
    position, direction, source_qop, time = _bound_to_physical(mean_f, surface_f, 'FASER')
    return _physical_to_bound(position, direction, source_qop if qop is None else qop,
                              time, surface_a, 'A')


def _jacobian_to_after(mean_f, surface_f, surface_a):
    """Numerical exact physical F->A bound Jacobian for covariance transport."""
    base = _physical_map(mean_f, surface_f, surface_a)
    jacobian = np.zeros((6, 6), dtype=float)
    steps = np.array([1e-4, 1e-4, 1e-6, 1e-6, 1e-8, 1e-4])
    for column, step in enumerate(steps):
        plus, minus = mean_f.copy(), mean_f.copy()
        plus[column] += step
        minus[column] -= step
        delta = _physical_map(plus, surface_f, surface_a) - _physical_map(minus, surface_f, surface_a)
        delta[2] = _angle_delta(_physical_map(plus, surface_f, surface_a)[2],
                                _physical_map(minus, surface_f, surface_a)[2])
        jacobian[:, column] = delta / (2. * step)
    # The time coordinate is a passive state component and must remain exact.
    jacobian[5, :] = 0.
    jacobian[5, 5] = 1.
    assert np.allclose(base[4], mean_f[4], rtol=0, atol=1e-14)
    return jacobian


def _standard_cloud(n_components, seed):
    if n_components < 14 or n_components % 2:
        raise ValueError('use an even component count of at least 14 for a full-rank 6D cloud')
    rng = np.random.default_rng(seed)
    half = rng.standard_normal((n_components // 2, 6))
    values = np.concatenate([half, -half])
    eigenvalues, eigenvectors = np.linalg.eigh(values.T @ values / n_components)
    inverse_root = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.T
    values = values @ inverse_root
    np.testing.assert_allclose(values.mean(axis=0), 0., rtol=0, atol=1e-14)
    np.testing.assert_allclose(values.T @ values / n_components, np.eye(6), rtol=0, atol=1e-12)
    return values


def _component_rows(before, after, n_components, seed, within_fraction):
    if not 0. < within_fraction <= 1.:
        raise ValueError('within-component covariance fraction must be in (0, 1]')
    mean_f = np.asarray(before['bound'], dtype=float)
    covariance_f = np.asarray(before['cov'], dtype=float)
    surface_f = np.asarray(before['surface'], dtype=float)
    surface_a = np.asarray(after['surface'], dtype=float)
    # Check the source row and its exact rotated counterpart before sampling.
    position_f, direction_f, qop, time = _bound_to_physical(mean_f, surface_f, 'FASER')
    np.testing.assert_allclose(position_f, np.asarray(before['position_native_mm']), rtol=0, atol=1e-8)
    np.testing.assert_allclose(direction_f, np.asarray(before['momentum_native_GeV']) /
                               np.linalg.norm(before['momentum_native_GeV']), rtol=0, atol=1e-10)
    mean_a = _physical_map(mean_f, surface_f, surface_a)
    np.testing.assert_allclose(mean_a[:2], [after['bound'][0], after['bound'][1]], rtol=0, atol=1e-7)
    standard = _standard_cloud(n_components, seed)
    chol = np.linalg.cholesky(covariance_f)
    within_f = within_fraction * covariance_f
    centers = mean_f + np.sqrt(1. - within_fraction) * (standard @ chol.T)
    centered = centers - mean_f
    total = centered.T @ centered / n_components + within_f
    scaled_error = np.max(np.abs(total - covariance_f) /
                          np.sqrt(np.outer(np.diag(covariance_f), np.diag(covariance_f))))
    assert scaled_error < 1e-11
    inner_offsets = np.sqrt(6. * within_fraction) * np.concatenate([chol.T, -chol.T])
    rows = []
    physical_error = 0.
    for sampled_f in centers:
        position, direction, sampled_qop, sampled_time = _bound_to_physical(
            sampled_f, surface_f, 'FASER')
        sampled_a = _physical_to_bound(position, direction, sampled_qop, sampled_time,
                                       surface_a, 'A')
        jacobian = _jacobian_to_after(sampled_f, surface_f, surface_a)
        covariance_a = jacobian @ within_f @ jacobian.T
        covariance_a = (covariance_a + covariance_a.T) / 2.
        np.linalg.cholesky(covariance_a)
        for delta in inner_offsets:
            exact_direction = _bound_to_physical(sampled_f + delta, surface_f, 'FASER')[1]
            linear_direction = _bound_to_physical(sampled_a + jacobian @ delta, surface_a, 'A')[1]
            physical_error = max(physical_error, float(np.linalg.norm(exact_direction - linear_direction)))
        assert np.isfinite(sampled_f).all() and np.isfinite(sampled_a).all()
        rows.append((sampled_f, sampled_a, within_f, covariance_a))
    return rows, dict(total_F_covariance_scaled_error=float(scaled_error),
        maximum_inner_sigma_direction_error=physical_error,
        inner_sigma_points_per_component=len(inner_offsets),
        finite_width_chart_mapping='first order at each component mean; convergence required')


def _line(entry, seed, component, mean, covariance, surface):
    values = [entry, seed, component, *mean.tolist(),
              *covariance.reshape(-1).tolist(), *surface.reshape(-1).tolist()]
    return ' '.join(format(v, '.17g') if isinstance(v, (float, np.floating)) else str(int(v))
                    for v in values) + '\n'


def prepare(root, output, n_components, within_fraction):
    output.mkdir(exist_ok=False)
    code = output / 'code'
    (code / 'scripts').mkdir(parents=True)
    (code / 'configs').mkdir()
    repo = Path(__file__).resolve().parents[1]
    for name in ['run_slope_task.py', 'prepare_cluster_measurements.py']:
        shutil.copy2(repo / 'scripts' / name, code / 'scripts' / name)
    shutil.copy2(repo / 'configs' / 'slope_tracking.py', code / 'configs' / 'slope_tracking.py')
    tasks = {variant: [] for variant in VARIANTS}
    maps = json.loads((root / 'physical_prior_material_bins_v4/manifest.json').read_text())
    manifest = dict(n_components=n_components, within_component_covariance_fraction=within_fraction,
                    samples={}, rotation=ROTATE_F_TO_A.tolist(),
                    material_maps={v:maps[k] for v,k in [('beforeGeo','source_map'), ('afterGeo','corrected_map')]},
                    scope='Finite moment-matched split of the recorded beforeGeo Gaussian. '
                          'Exact component means, per-component first-order covariance mapping, '
                          'native Calypso CKF/update/isolated smoothing. '
                          'Finite-width convergence and posterior likelihood aggregation not assumed.')
    for dsid in SAMPLES:
        source_rows = {variant: _read_seed_rows(
            root / 'common_prior_ckf_v1' / 'smoke' / variant / dsid / 'audit' /
            'prior_audit.jsonl.gz') for variant in VARIANTS}
        common = sorted(set(source_rows['beforeGeo']) & set(source_rows['afterGeo']))
        assert common
        lines = {variant: [] for variant in VARIANTS}
        records = []
        input_checks = []
        components_per_event = {}
        for entry, seed in common:
            before, after = source_rows['beforeGeo'][(entry, seed)], source_rows['afterGeo'][(entry, seed)]
            assert before['hits'] == after['hits'] and before['bound'][4] == after['bound'][4]
            rows, checks = _component_rows(before, after, n_components,
                2026111400 + int(dsid) * 10000 + entry * 101 + seed, within_fraction)
            input_checks.append(dict(entry=entry, seed=seed, **checks))
            component_base = components_per_event.get(entry, 0)
            for draw, (mean_f, mean_after, cov_before, cov_after) in enumerate(rows):
                component = component_base + draw
                lines['beforeGeo'].append(_line(entry, seed, component, mean_f,
                                                cov_before,
                                                np.asarray(before['surface'], dtype=float)))
                lines['afterGeo'].append(_line(entry, seed, component, mean_after, cov_after,
                                               np.asarray(after['surface'], dtype=float)))
                records.append(dict(entry=entry, seed=seed, component=component,
                                    mean_F=mean_f.tolist(), mean_A=mean_after.tolist()))
            components_per_event[entry] = component_base + len(rows)
        group = output / dsid
        group.mkdir()
        source_tasks = {}
        for variant in VARIANTS:
            source_task = json.loads((root / 'common_prior_ckf_v1' / 'smoke' / variant / dsid /
                                      'audit' / 'source_task.json').read_text())
            source_tasks[variant] = source_task
            component_file = group / f'{variant}_components.txt'
            component_file.write_text(''.join(lines[variant]))
            options = ['--no-seed-writer', '--actual-covariance', '--isolate-smoothing',
                       '--prior-component-audit', '--seed-component-input-file', str(component_file)]
            task = dict(source_task, runner=str(code / 'configs' / 'slope_tracking.py'),
                        options=options, output=str(group / variant), component_mode='finite_common_prior',
                        component_file=str(component_file), material=manifest['material_maps'][variant])
            task_path = group / f'{variant}_source.json'
            dump(task_path, task)
            tasks[variant].append(dict(kind='analysis', command=['python3', str(code / 'scripts' /
                                   'run_slope_task.py'), '--task', str(task_path)], output=task['output']))
        manifest['samples'][dsid] = dict(common_seed_rows=len(common), source_tasks=source_tasks,
                                         input_checks=input_checks,
                                         component_rows=len(records), component_files={
                                             v: str(output / dsid / f'{v}_components.txt') for v in VARIANTS})
    for variant in VARIANTS:
        dump(output / f'{variant}_tasks.json', tasks[variant])
    dump(output / 'manifest.json', manifest)


def validate(root, output):
    manifest = json.loads((output / 'manifest.json').read_text())
    checks = []
    for dsid, sample in manifest['samples'].items():
        for variant in VARIANTS:
            task = json.loads((output / dsid / f'{variant}_source.json').read_text())
            folder = Path(task['output'])
            assert (folder / 'exit_code.txt').read_text().strip() == '0', folder
            audit = folder / 'prior_audit.jsonl.gz'
            rows = []
            with gzip.open(audit, 'rt') as stream:
                for line in stream:
                    row = json.loads(line)
                    if row.get('stage') in ('seed', 'found', 'selected') and 'input_component' in row:
                        rows.append(row)
            input_count = sum(1 for line in Path(task['component_file']).read_text().splitlines() if line.strip())
            assert input_count == sample['component_rows']
            assert rows, (dsid, variant)
            checks.append(dict(dsid=dsid, variant=variant, input_components=input_count,
                               observed_component_rows=len(rows), output=str(folder)))
    dump(output / 'validation.json', dict(all_exit_zero=True, checks=checks))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['prepare', 'validate'], required=True)
    parser.add_argument('--n-components', type=int, default=32)
    parser.add_argument('--within-fraction', type=float, default=1./16.,
                        help='Fraction of original covariance within each component; remaining variance is between means')
    parser.add_argument('--component-directory', default='common_prior_components_v4')
    parser.add_argument('--output-name', default='smoke_v1',
                        help='Subdirectory below the selected component directory')
    args = parser.parse_args()
    root = args.root.resolve()
    output = root / args.component_directory / args.output_name
    if args.phase == 'prepare':
        if args.n_components < 14 or args.n_components % 2 or not 0. < args.within_fraction <= 1.:
            parser.error('use an even component count >=14 and within-fraction in (0,1]')
        prepare(root, output, args.n_components, args.within_fraction)
    else:
        validate(root, output)
    print(output)


if __name__ == '__main__':
    main()
