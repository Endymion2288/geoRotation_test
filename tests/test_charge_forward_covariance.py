import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from charge_forward_covariance import chart, covariance_balance, forward_states
from charge_refit_input import ROTATE_F_TO_A


class ForwardCovarianceDiagnostics(unittest.TestCase):
    def test_common_chart_retains_direction_and_correlated_covariance(self):
        rng = np.random.default_rng(52)
        factor = rng.normal(size=(6, 6))
        covariance = factor @ factor.T
        for phi, theta in [(1.1, .006), (-2.1, 1e-8), (np.pi-1e-9, .003)]:
            bound = np.array([2., -3., phi, theta, -.002, 1.])
            mean, jacobian = chart(bound, 'beforeGeo')
            direction = ROTATE_F_TO_A @ np.array([np.cos(phi)*np.sin(theta),
                np.sin(phi)*np.sin(theta), np.cos(theta)])
            regular_direction = np.array([np.cos(mean[2])*np.cos(mean[3]),
                np.sin(mean[2])*np.cos(mean[3]), np.sin(mean[3])])
            np.testing.assert_allclose(direction, regular_direction, rtol=0, atol=3e-16)
            # Independently differentiate Cartesian directions on both ends.
            native = ROTATE_F_TO_A @ np.array([[-np.sin(phi)*np.sin(theta), np.cos(phi)*np.cos(theta)],
                [np.cos(phi)*np.sin(theta), np.sin(phi)*np.cos(theta)], [0., -np.sin(theta)]])
            regular = np.array([[-np.sin(mean[2])*np.cos(mean[3]), -np.cos(mean[2])*np.sin(mean[3])],
                [np.cos(mean[2])*np.cos(mean[3]), -np.sin(mean[2])*np.sin(mean[3])],
                [0., np.cos(mean[3])]])
            np.testing.assert_allclose(native, regular @ jacobian[2:4, 2:4], rtol=1e-13, atol=2e-16)
            c = covariance[2:4, 2:4]
            np.testing.assert_allclose(native @ c @ native.T,
                regular @ jacobian[2:4, 2:4] @ c @ jacobian[2:4, 2:4].T @ regular.T,
                rtol=2e-13, atol=1e-20)

    def test_after_polar_to_latitude_inverts_covariance_cross_terms(self):
        mean = np.array([1., 2., .01, np.pi/2.+.002, -.03, 1.])
        mapped, jacobian = chart(mean, 'afterGeo')
        self.assertAlmostEqual(mapped[3], -.002)
        c = np.eye(6)
        c[0, 3] = c[3, 0] = .2
        self.assertEqual((jacobian @ c @ jacobian.T)[0, 3], -.2)
        self.assertEqual((jacobian @ c @ jacobian.T)[3, 3], 1.)

    def test_covariance_balance_distinguishes_inheritance_and_added_noise(self):
        ca = np.eye(6)
        cb = ca.copy()
        cb[2, 2] += .3
        ja = jb = np.eye(6)
        noise_a, noise_b = np.zeros((6, 6)), np.zeros((6, 6))
        noise_a[3, 3], noise_b[3, 3] = .2, .7
        result = covariance_balance(ca, cb, ja, jb, ca+noise_a, cb+noise_b)
        np.testing.assert_allclose(result['inherited_C']['direction_matrix_rad2'], [[.3, 0.], [0., 0.]])
        np.testing.assert_allclose(result['saved_J']['direction_matrix_rad2'], 0.)
        np.testing.assert_allclose(result['remaining_balance']['direction_matrix_rad2'], [[0., 0.], [0., .5]])
        self.assertLess(result['decomposition_scaled_residual'], 1e-15)

    def test_forward_order_uses_linked_indices_and_keeps_passive_states(self):
        states = [dict(acts_state_index=0, acts_previous_index=1, path_length_mm=100., measurement_flag=True),
                  dict(acts_state_index=1, acts_previous_index=2, path_length_mm=90., measurement_flag=False),
                  dict(acts_state_index=2, acts_previous_index=4294967295, path_length_mm=80., measurement_flag=True)]
        ordered = forward_states(dict(states={3: states}), 3)
        self.assertEqual([r['path_length_mm'] for r in ordered], [80., 90., 100.])
        self.assertFalse(ordered[1]['measurement_flag'])

    def test_signed_path_length_does_not_replace_native_link_order(self):
        states = [dict(acts_state_index=0, acts_previous_index=1, path_length_mm=-10.),
                  dict(acts_state_index=1, acts_previous_index=4294967295, path_length_mm=0.)]
        self.assertEqual([r['acts_state_index'] for r in forward_states(dict(states={0: states}), 0)], [1, 0])


if __name__ == '__main__':
    unittest.main()
