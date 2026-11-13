#!/usr/bin/env python3
"""Validate full seed Jacobians recorded by the actual matched-prior Calypso."""
import argparse
import json
from pathlib import Path

import numpy as np

from validate_charge_audit import compare_tables


FRAME = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])


def direction(phi, theta):
    return np.stack([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)], axis=-1)


def derivative(phi, theta):
    return np.array([[-np.sin(phi)*np.sin(theta), np.cos(phi)*np.cos(theta)],
                     [np.cos(phi)*np.sin(theta), np.sin(phi)*np.cos(theta)], [0., -np.sin(theta)]])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    summary = {}
    rng = np.random.default_rng(7541)
    for dsid in ['100044', '100116', '100117']:
        root = args.root/'prior_smoke_v1'/dsid
        for mode in ['default', 'matched', 'matched_audit']:
            assert (root/mode/'exit_code.txt').read_text().strip() == '0'
        default = compare_tables(root/'default', args.root/'charge_smoke/afterGeo'/dsid/'on')
        audit = compare_tables(root/'matched_audit', root/'matched')
        with np.load(root/'default/charge_lineage.npz') as a, np.load(root/'matched/charge_lineage.npz') as b:
            sa, sb = a['stage'] == 0, b['stage'] == 0
            for key in ['source_entry', 'seed', 'qop_GeVinv', 'phi', 'theta', 'loc0', 'loc1']:
                np.testing.assert_array_equal(a[key][sa], b[key][sb], err_msg='Same native seed mean: '+key)
        raw = np.atleast_1d(np.genfromtxt(root/'matched_audit/seed_covariance.csv', delimiter=',', names=True))
        maximum = dict(jacobian_finite_difference_scaled=0., covariance_mapping_relative=0., physical_covariance_relative=0.)
        directions, covariances = [], []
        for row in raw:
            d = np.array([row['directionF_'+str(i)] for i in range(3)])
            phi, theta = np.arctan2(d[1], d[0]), np.arctan2(np.hypot(d[0], d[1]), d[2])
            matrices = {name: np.array([[row[f'{name}_{i}_{j}'] for j in range(n)] for i in range(n)])
                        for name, n in [('basisF', 3), ('basisA', 3), ('covF', 6), ('jacobian', 6), ('covA', 6)]}
            cf, j, ca = [matrices[k] for k in ['covF', 'jacobian', 'covA']]
            np.testing.assert_allclose(np.diag(cf)[:5], [81., .0361, .1089, 8.41e-6, 1e-5], rtol=1e-14)
            expected = np.eye(6)
            expected[:2, :2] = matrices['basisA'][:, :2].T@FRAME@matrices['basisF'][:, :2]
            for column, step in [(2, 1e-5), (3, 1e-7)]:
                angles = np.array([phi, theta])
                high, low = angles.copy(), angles.copy()
                high[column-2] += step
                low[column-2] -= step
                ah, al = FRAME@direction(*high), FRAME@direction(*low)
                ph, pl = np.arctan2(ah[1], ah[0]), np.arctan2(al[1], al[0])
                th, tl = np.arctan2(np.hypot(ah[0], ah[1]), ah[2]), np.arctan2(np.hypot(al[0], al[1]), al[2])
                expected[2, column] = np.arctan2(np.sin(ph-pl), np.cos(ph-pl))/(2*step)
                expected[3, column] = (th-tl)/(2*step)
            error = float(np.max(np.abs(expected-j)/np.maximum(1., np.abs(j))))
            maximum['jacobian_finite_difference_scaled'] = max(maximum['jacobian_finite_difference_scaled'], error)
            assert error < 1e-7
            mapping = j@cf@j.T
            relative = float(np.linalg.norm(ca-mapping)/np.linalg.norm(mapping))
            maximum['covariance_mapping_relative'] = max(maximum['covariance_mapping_relative'], relative)
            assert relative < 1e-12
            a = FRAME@d
            pa, ta = np.arctan2(a[1], a[0]), np.arctan2(np.hypot(a[0], a[1]), a[2])
            before = FRAME@derivative(phi, theta)@cf[2:4, 2:4]@derivative(phi, theta).T@FRAME.T
            after = derivative(pa, ta)@ca[2:4, 2:4]@derivative(pa, ta).T
            relative = float(np.linalg.norm(after-before)/np.linalg.norm(before))
            maximum['physical_covariance_relative'] = max(maximum['physical_covariance_relative'], relative)
            assert relative < 1e-11
            directions.append([phi, theta, pa, ta])
            covariances.append((cf, ca))
        nonlinear = []
        directions = np.array(directions)
        representative = np.unique(np.argsort(directions[:, 1])[np.linspace(0, len(raw)-1, 5).astype(int)])
        for i in representative:
            phi, theta, pa, ta = directions[i]
            cf, ca = covariances[i]
            f = rng.multivariate_normal([phi, theta], cf[2:4, 2:4], size=100000)
            a = rng.multivariate_normal([pa, ta], ca[2:4, 2:4], size=100000)
            physical_f = direction(f[:, 0], f[:, 1])@FRAME.T
            physical_a = direction(a[:, 0], a[:, 1])
            cvf, cva = np.cov(physical_f, rowvar=False), np.cov(physical_a, rowvar=False)
            nonlinear.append(dict(seed_row=int(i), theta_F_rad=float(theta),
                sampled_physical_covariance_relative_difference=float(np.linalg.norm(cva-cvf)/np.linalg.norm(cvf)),
                mean_direction_difference=float(np.linalg.norm(physical_a.mean(axis=0)-physical_f.mean(axis=0)))))
        summary[dsid] = dict(default_fields_exact=len(default), audit_fields_exact=len(audit),
            default_reconstruction_unchanged=True, audit_read_only=True, seed_means_identical=True,
            actual_seed_rows=len(raw), maximum=maximum, nonlinear_prior_sampling=nonlinear)
    result = dict(samples=summary, jacobian_scope='Full local and angular bound Jacobian at each actual native seed mean',
        nonlinear_scope='Algebraic sampling of recorded seed priors only; not an independent tracking/simulation',
        interpretation='First-order physical seed covariance closes; finite-width spherical Gaussian priors need not define identical nonlinear physical distributions')
    (args.output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
