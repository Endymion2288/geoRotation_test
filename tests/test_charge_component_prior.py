import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_charge_component_tasks import (
    ROTATE_F_TO_A, _bound_to_physical, _component_rows, _jacobian_to_after,
    _physical_map, _standard_cloud,
)
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from charge_refit_input import physical


class CommonPhysicalPrior(unittest.TestCase):
    def seeds(self):
        surface_f = np.eye(4)
        surface_a = np.eye(4)
        surface_a[:3, :3] = ROTATE_F_TO_A
        mean = np.array([2., 3., -1.2, .003, .002, 0.])
        covariance = np.diag([81., .0361, .1089, 8.41e-6, 1e-5, 1.])
        position, direction, _, _ = _bound_to_physical(mean, surface_f, 'FASER')
        before = dict(bound=mean, cov=covariance, surface=surface_f,
                      position_native_mm=position, momentum_native_GeV=direction/abs(mean[4]))
        after = dict(surface=surface_a, bound=_physical_map(mean, surface_f, surface_a))
        return before, after

    def test_total_prior_covariance_is_preserved_instead_of_doubled(self):
        before, after = self.seeds()
        rows, checks = _component_rows(before, after, 32, 152, 1./16.)
        means = np.array([row[0] for row in rows])
        within = np.mean([row[2] for row in rows], axis=0)
        centered = means - before['bound']
        total = centered.T @ centered/len(rows) + within
        np.testing.assert_allclose(means.mean(axis=0), before['bound'], rtol=0, atol=1e-14)
        np.testing.assert_allclose(total, before['cov'], rtol=1e-12, atol=1e-13)
        self.assertLess(checks['total_F_covariance_scaled_error'], 1e-11)

    def test_rotated_component_means_are_the_same_physical_input(self):
        before, after = self.seeds()
        rows, _ = _component_rows(before, after, 32, 152, 1./16.)
        for mean_f, mean_a, _, covariance_a in rows:
            f = _bound_to_physical(mean_f, before['surface'], 'FASER')
            a = _bound_to_physical(mean_a, after['surface'], 'A')
            for i in [0, 1]:
                np.testing.assert_allclose(f[i], a[i], rtol=0, atol=1e-14)
            self.assertEqual(f[2:], a[2:])
            self.assertGreater(np.linalg.eigvalsh(covariance_a).min(), 0.)

    def test_finite_width_error_is_measured_and_shrinks(self):
        before, after = self.seeds()
        _, broad = _component_rows(before, after, 32, 152, 1./16.)
        _, narrow = _component_rows(before, after, 32, 152, 1./256.)
        self.assertGreater(broad['maximum_inner_sigma_direction_error'], 0.)
        self.assertLess(narrow['maximum_inner_sigma_direction_error'],
                        broad['maximum_inner_sigma_direction_error']/5.)

    def test_invalid_antithetic_cloud_sizes_are_rejected(self):
        for size in [0, 8, 15, 33]:
            with self.assertRaises(ValueError):
                _standard_cloud(size, 1)

    def test_analytic_map_preserves_correlated_physical_covariance(self):
        """Cross terms and near-pole/equivalent angles survive the coordinate map."""
        surface_f = np.eye(4)
        angle = .23
        surface_f[:3, :3] = [[np.cos(angle), -np.sin(angle), 0.],
                             [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]]
        surface_f[:3, 3] = [10., -7., 2050.]
        rotate = np.eye(4)
        rotate[:3, :3] = ROTATE_F_TO_A
        surface_a = rotate @ surface_f
        # An independent local basis and local origin on the same plane.
        basis = np.eye(4)
        basis[:2, :2] = [[np.cos(.37), -np.sin(.37)], [np.sin(.37), np.cos(.37)]]
        basis[:2, 3] = [4., -3.]
        surface_a = surface_a @ basis
        rng = np.random.default_rng(124)
        factor = rng.normal(size=(6, 6)) + 3.*np.eye(6)
        scale = np.diag([2., 3., .02, .0003, .0002, 1.])
        covariance = scale @ factor @ factor.T @ scale
        for phi, theta in [(1.2, .4), (np.pi-1e-9, .003), (-np.pi+1e-9, .002),
                           (1.1, 1e-8), (-.8, -.004), (.2, np.pi+.002)]:
            with self.subTest(phi=phi, theta=theta):
                mean = np.array([2., 3., phi, theta, -.002, .5])
                mapped = _physical_map(mean, surface_f, surface_a)
                jacobian = _jacobian_to_after(mean, surface_f, surface_a, method='analytic')
                def recorded(bound, cov, surface, frame):
                    position = (surface @ np.r_[bound[:2], 0., 1.])[:3]
                    direction = _bound_to_physical(bound, surface, frame)[1]
                    if frame == 'A':
                        direction = ROTATE_F_TO_A @ direction
                    suffix = 'native' if frame == 'FASER' else 'A'
                    return dict(bound=bound, cov=cov, surface=surface,
                                **{'position_'+suffix+'_mm': position,
                                   'momentum_'+suffix+'_GeV': direction/abs(bound[4])})
                mf, cf = physical(recorded(mean, covariance, surface_f, 'FASER'), 'beforeGeo')
                ma, ca = physical(recorded(mapped, jacobian @ covariance @ jacobian.T,
                                           surface_a, 'A'), 'afterGeo')
                np.testing.assert_allclose(mf, ma, rtol=0, atol=1e-12)
                diagonal = np.maximum(np.abs(np.diag(cf)), 1e-30)
                # Unit-direction C has a constrained radial mode. Near a pole
                # its individual diagonal can vanish; use the vector variance
                # scale so the check measures C accuracy rather than a ratio
                # to a floating-point representation of that zero mode.
                diagonal[3:6] = np.trace(cf[3:6, 3:6])
                error = np.max(np.abs(cf-ca)/np.sqrt(np.outer(diagonal, diagonal)))
                self.assertLess(error, 2e-12)
                np.testing.assert_array_equal(jacobian[4:], np.eye(6)[4:])
                numerical = _jacobian_to_after(mean, surface_f, surface_a)
                np.testing.assert_allclose(jacobian, numerical, rtol=1e-6, atol=3e-9)

    def test_analytic_component_option_changes_only_after_covariance(self):
        before, after = self.seeds()
        legacy, _ = _component_rows(before, after, 32, 152, 1./256.)
        explicit, _ = _component_rows(before, after, 32, 152, 1./256., 'finite_difference')
        analytic, _ = _component_rows(before, after, 32, 152, 1./256., 'analytic')
        changed = False
        for old, same, new in zip(legacy, explicit, analytic):
            for a, b in zip(old, same):
                np.testing.assert_array_equal(a, b)
            for index in [0, 1, 2]:
                np.testing.assert_array_equal(old[index], new[index])
            np.linalg.cholesky(new[3])
            changed |= not np.array_equal(old[3], new[3])
        self.assertTrue(changed)


if __name__ == '__main__':
    unittest.main()
