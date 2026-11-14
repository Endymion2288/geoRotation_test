#!/usr/bin/env python3
"""Sample a common finite physical prior before invoking native Calypso."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np


FRAME = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    output = root/'physical_prior_transport_v1'
    output.mkdir(exist_ok=False)
    code = output/'code'
    (code/'configs').mkdir(parents=True)
    (code/'scripts').mkdir()
    shutil.copy2(repo/'configs/charge_prior_transport.py', code/'configs/charge_prior_transport.py')
    shutil.copy2(repo/'scripts/run_slope_transport.py', code/'scripts/run_slope_transport.py')
    # Reuse the Calypso probe whose component registration and both native
    # geometry/field versions were already validated in S05.
    libraries = {v: str(root.parent.parent/'stage05_slopes/diagnostics_v1'/('software_'+v)/'lib')
                 for v in ['beforeGeo', 'afterGeo']}
    materials = {v: str(root.parent.parent/'stage02_closure'/name) for v, name in [
        ('beforeGeo', 'material-current-source.json'), ('afterGeo', 'material-current-rotated.json')]}
    tasks = {v: [] for v in libraries}
    samples = {}
    for dsid in ['100044', '100116', '100117']:
        source = root/'prior_smoke_v1'/dsid/'matched_audit'
        raw = np.atleast_1d(np.genfromtxt(source/'seed_covariance.csv', delimiter=',', names=True))
        before_angles = np.array([[np.arctan2(r['directionF_1'], r['directionF_0']),
            np.arctan2(np.hypot(r['directionF_0'], r['directionF_1']), r['directionF_2'])] for r in raw])
        selected = np.unique(np.argsort(before_angles[:, 1])[np.linspace(0, len(raw)-1, 5).astype(int)])
        rng = np.random.default_rng(2026111400+int(dsid))
        lines = {v: [] for v in libraries}
        clouds = []
        seed_definitions = []
        truth = np.load(source/'truth.npz')
        for ordinal, i in enumerate(selected):
            row = raw[i]
            basis_f, basis_a = [np.array([[row[f'{name}_{j}_{k}'] for k in range(3)] for j in range(3)])
                                for name in ['basisF', 'basisA']]
            mean_a = np.array([row[f'meanA_{j}'] for j in range(6)])
            mean = mean_a.copy()
            # Define the physical prior on the common plane z_F=0 mm. The
            # recorded transverse mean and all covariance entries are retained;
            # this is an explicit input plane, not the original seed surface.
            position_f = FRAME.T@basis_a[:, :2]@mean_a[:2]
            mean[:2] = basis_f[:, :2].T@position_f
            mean[2:4] = before_angles[i]
            covariance = np.array([[row[f'covF_{j}_{k}'] for k in range(6)] for j in range(6)])
            draws = rng.multivariate_normal(mean, covariance, size=512)
            st = np.sin(draws[:, 3])
            direction_f = np.column_stack([np.cos(draws[:, 2])*st, np.sin(draws[:, 2])*st, np.cos(draws[:, 3])])
            momentum_f = direction_f/np.abs(draws[:, 4, None])
            positions_f = draws[:, :2]@basis_f[:, :2].T
            np.testing.assert_allclose(positions_f[:, 2], 0., rtol=0, atol=1e-12)
            local_event = int(row['local_event_index'])
            primary = (truth['source_entry'] == local_event) & (truth['status'] == 1) & \
                (np.abs(truth['pdg']) == 13) & (truth['barcode'] > 0) & (truth['barcode'] < 200000)
            assert primary.sum() == 1
            barcode = int(truth['barcode'][primary][0])
            seed = int(row['seed'])
            assert np.isfinite(momentum_f).all() and np.all(draws[:, 4] != 0)
            for variant in lines:
                rotation = FRAME if variant == 'afterGeo' else np.eye(3)
                positions, momenta = positions_f@rotation.T, momentum_f@rotation.T
                for draw in range(len(draws)):
                    fields = [local_event, seed, ordinal*512+draw, barcode, 0, -1, 0, 0,
                              int(np.sign(draws[draw, 4]))]
                    lines[variant].append(' '.join(map(str, fields))+' '+
                        ' '.join(format(v, '.17g') for v in np.r_[positions[draw], momenta[draw]])+'\n')
            clouds.append(dict(source_entry=np.full(512, local_event, dtype=np.int64),
                candidate=np.arange(ordinal*512, (ordinal+1)*512, dtype=np.int64),
                position_F_mm=positions_f, momentum_F_GeV=momentum_f,
                qop_GeVinv=draws[:, 4]))
            seed_definitions.append(dict(source_seed_row=int(i), source_entry=local_event, seed=seed,
                mean_bound_F_at_common_plane=mean.tolist(), full_recorded_covariance_F=covariance.tolist(),
                angular_pole_crossing_draws=int(np.count_nonzero(draws[:, 3] < 0)),
                opposite_to_seed_charge_draws=int(np.count_nonzero(draws[:, 4]*mean[4] < 0))))
        truth.close()
        group = output/dsid
        group.mkdir()
        np.savez_compressed(group/'physical_cloud.npz', **{key: np.concatenate([r[key] for r in clouds]) for key in clouds[0]})
        original_task = json.loads((source/'source_task.json').read_text())
        for variant in libraries:
            states = group/(variant+'_states.txt')
            states.write_text(''.join(lines[variant]))
            task = dict(dsid=dsid, input=original_task['input'], variant=variant, material=materials[variant],
                events=20, skip=0, states=str(states), runner=str(code/'configs/charge_prior_transport.py'),
                output=str(group/variant))
            path = group/(variant+'_source.json')
            path.write_text(json.dumps(task, indent=2)+'\n')
            tasks[variant].append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_transport.py'),
                '--task', str(path)], output=task['output']))
        samples[dsid] = dict(sample_points=sum(len(c['candidate']) for c in clouds),
            source_prior=str(source/'seed_covariance.csv'), seeds=seed_definitions)
    for variant, jobs in tasks.items():
        (output/(variant+'_tasks.json')).write_text(json.dumps(jobs, indent=2)+'\n')
    manifest = dict(samples=samples, libraries=libraries, materials=materials,
        source='Full recorded F-frame seed Gaussian, sampled once per representative; identical physical points in both native Calypso versions',
        input_plane='Common z_F=0 mm; recorded transverse mean, direction, q/p, time and full covariance',
        targets=['z_F=0 mm', 'VetoNu z_F=-3112 mm', 'first station-1 SiHit plane of the primary in the source RDO event'],
        interpretation='Tests the full finite physical ensemble through actual geometry, field and material transport. Per-point covariance is absent: this does not assert complete Gaussian CKF/refit equivalence or replace tracking.',
        samples_per_seed=512, material_and_field_native=True, independent_tracking=False)
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(output)


if __name__ == '__main__':
    main()
