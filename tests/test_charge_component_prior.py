import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_charge_component_tasks import (
    ROTATE_F_TO_A, _bound_to_physical, _component_rows, _physical_map, _standard_cloud,
)


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


if __name__ == '__main__':
    unittest.main()
