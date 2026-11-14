import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from bound_surface_coordinates import direction_f, rotate_tangent, surface_position, ROTATE_F_TO_A


class NativeSurfaceCoordinates(unittest.TestCase):
    def test_disc_radius_phi_are_not_plane_xy(self):
        bound = np.array([100., np.pi/3., .1, .003, -.02, 0.])
        transform = np.eye(4)
        transform[:3, 3] = [2., -3., 137.401]
        position, jacobian = surface_position(bound, transform, 2)
        np.testing.assert_allclose(position, [52., 50.*np.sqrt(3.)-3., 137.401], atol=2e-14)
        np.testing.assert_allclose(jacobian, [[.5, -50.*np.sqrt(3.)], [np.sqrt(3.)/2., 50.], [0., 0.]], atol=2e-14)
        self.assertGreater(np.linalg.norm(position-surface_position(bound, transform, 4)[0]), 90.)

    def test_global_rotation_preserves_disc_local_coordinates_and_full_tangent_covariance(self):
        bound = np.array([101.3, .62, 2.1, .005, -.002, 1.])
        rng = np.random.default_rng(93)
        factor = rng.normal(size=(6, 6))
        covariance = factor @ factor.T
        before, rotation = np.eye(4), np.eye(4)
        before[:3, 3] = [0., 0., 137.401]
        rotation[:3, :3] = ROTATE_F_TO_A
        after = rotation @ before
        mapped, mapped_c = rotate_tangent(bound, covariance, before, after)
        np.testing.assert_array_equal(mapped[:2], bound[:2])
        pa, ja = surface_position(bound, before, 2)
        pb, jb = surface_position(mapped, after, 2)
        np.testing.assert_allclose(pa, ROTATE_F_TO_A.T @ pb, atol=1e-14)
        da, ga = direction_f(bound, 'beforeGeo')
        db, gb = direction_f(mapped, 'afterGeo')
        np.testing.assert_allclose(da, db, atol=3e-16)
        # The mixed position/direction C checks off-diagonal covariance too.
        fa, fb = np.zeros((6, 6)), np.zeros((6, 6))
        fa[:3, :2], fb[:3, :2] = ja, ROTATE_F_TO_A.T @ jb
        fa[3:], fb[3:] = ga, gb
        np.testing.assert_allclose(fa @ covariance @ fa.T, fb @ mapped_c @ fb.T, rtol=1e-12, atol=2e-11)

    def test_unknown_native_surface_type_is_not_assumed_planar(self):
        with self.assertRaisesRegex(ValueError, 'surface type'):
            surface_position(np.zeros(6), np.eye(4), 1)


if __name__ == '__main__':
    unittest.main()
