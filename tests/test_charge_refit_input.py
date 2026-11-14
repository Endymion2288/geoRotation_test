"""Coordinate checks for the recorded native refit-prior comparison."""
from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'analysis'))
sys.path.insert(0, str(ROOT/'scripts'))
from charge_refit_input import physical, ROTATE_F_TO_A


def pair():
    phi, theta = .7, .02
    direction = np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
    transformed = ROTATE_F_TO_A @ direction
    pa, ta = np.arctan2(transformed[1], transformed[0]), np.arccos(transformed[2])
    derivative = np.column_stack(([-np.sin(phi)*np.sin(theta), np.cos(phi)*np.sin(theta), 0.],
                                  [np.cos(phi)*np.cos(theta), np.sin(phi)*np.cos(theta), -np.sin(theta)]))
    angular_jacobian = np.array([[-np.sin(pa)/np.sin(ta), np.cos(pa)/np.sin(ta), 0.],
                                  [np.cos(pa)*np.cos(ta), np.sin(pa)*np.cos(ta), -np.sin(ta)]])
    jacobian = np.eye(6)
    jacobian[2:4, 2:4] = angular_jacobian @ ROTATE_F_TO_A @ derivative
    surface = np.eye(4)
    surface[:3, :3] = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
    surface[:3, 3] = [0, 0, 1200]
    after_surface = surface.copy()
    after_surface[:3] = ROTATE_F_TO_A @ surface[:3]
    mean = np.array([20., -30., phi, theta, .0004, 0.])
    after_mean = mean.copy()
    after_mean[2:4] = pa, ta
    scales = np.array([1., .1, .003, .0002, 1e-5, 1.])
    lower = np.eye(6)
    lower[3, 2], lower[4, 3], lower[1, 0] = .5, -.2, .3
    covariance = scales[:, None]*(lower @ lower.T)*scales[None, :]
    rows = []
    for variant, bound, transform, cov, vec in [
            ('beforeGeo', mean, surface, covariance, direction),
            ('afterGeo', after_mean, after_surface, jacobian @ covariance @ jacobian.T, transformed)]:
        suffix = 'native' if variant == 'beforeGeo' else 'A'
        row = dict(bound=bound.tolist(), cov=cov.tolist(), surface=transform.tolist())
        row['position_'+suffix+'_mm'] = (transform[:3, :2] @ bound[:2]+transform[:3, 3]).tolist()
        row['momentum_'+suffix+'_GeV'] = (vec/abs(bound[4])).tolist()
        rows.append(row)
    return rows


class RefitInputCoordinates(unittest.TestCase):
    def test_physical_mean_and_complete_covariance_rotation(self):
        before, after = pair()
        ma, ca = physical(before, 'beforeGeo')
        mb, cb = physical(after, 'afterGeo')
        np.testing.assert_allclose(ma, mb, rtol=0, atol=1e-12)
        np.testing.assert_allclose(ca, cb, rtol=2e-12, atol=1e-20)

    def test_uniform_inflation_is_covariant(self):
        before, after = pair()
        for row in (before, after):
            row['cov'] = (10*np.asarray(row['cov'])).tolist()
        np.testing.assert_allclose(physical(before, 'beforeGeo')[1], physical(after, 'afterGeo')[1],
                                   rtol=2e-12, atol=1e-20)

    def test_diagonal_inflation_changes_the_physical_prior(self):
        before, after = pair()
        for row in (before, after):
            covariance = np.asarray(row['cov'])
            row['cov'] = (covariance+9*np.diag(np.diag(covariance))).tolist()
        ca, cb = physical(before, 'beforeGeo')[1], physical(after, 'afterGeo')[1]
        self.assertGreater(np.linalg.norm(ca[3:6, 3:6]-cb[3:6, 3:6])/np.linalg.norm(ca[3:6, 3:6]), .1)


if __name__ == '__main__':
    unittest.main()
