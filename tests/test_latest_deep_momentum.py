import subprocess
from pathlib import Path
import unittest

import numpy as np

from analysis import latest_deep_reduce as ldr
from analysis import latest_deep_study as lds


class LatestDeepMomentumTests(unittest.TestCase):
    def test_rotated_seed_plane_axes_and_cov_loc_mapping(self):
        # Verify CircleFitTrackSeedTool.cxx:391-408 algebra:
        # rot maps FASER (x_F, y_F, z_F) -> ACTS A (z_F, y_F, -x_F).
        rot = np.array([[0., 0., 1.],
                        [0., 1., 0.],
                        [-1., 0., 0.]])
        plane_rot = np.column_stack([
            [0., 1., 0.],    # col(0) = eBoundLoc0 axis in ACTS A
            [0., 0., -1.],   # col(1) = eBoundLoc1 axis in ACTS A
            [-1., 0., 0.],   # col(2) = normal in ACTS A
        ])
        pos_f = np.array([12.5, -3.2, -5.0])
        pos_a = rot @ pos_f
        center_a = np.array([pos_f[2], 0., 0.])
        local = plane_rot.T @ (pos_a - center_a)
        self.assertAlmostEqual(local[0], pos_f[1], places=12)  # loc0 == +y_F (precise 0.019 mm strip axis)
        self.assertAlmostEqual(local[1], pos_f[0], places=12)  # loc1 == +x_F (coarse 0.9 mm stereo axis)
        self.assertAlmostEqual(local[2], 0.0, places=12)

    def test_symmetric_wide_curvature_folds_to_left_momentum_spike(self):
        rng = np.random.default_rng(42)
        p_true = 500.0
        q_true = -1.0 / p_true
        # Symmetric zero-mean curvature error with std = 8.5 * |q_true|
        eps = rng.normal(0.0, 8.5, size=200000)
        q_reco = q_true * (1.0 + eps)
        self.assertAlmostEqual(float(np.median(q_reco / q_true - 1.0)), 0.0, delta=0.05)
        res_p = lds.p_residual(q_reco, p_true)
        desc = lds.describe_p(res_p)
        # Over 80% of symmetric wide curvature states map to negative p residual (< -0.1)
        self.assertGreater(desc['below_minus01'], 0.75)
        self.assertLess(desc['quantiles']['p50'], -0.65)
        self.assertLess(desc['modal_bin'][1], -0.80)

    def test_native_formula_division_by_zero_matches_ttreeformula(self):
        q = np.array([0.0, -0.002, 0.002])
        tq = np.array([-0.002, 0.0, 0.002])
        res = lds.native_formula_residual(q, tq)
        # q=0 -> p=0 -> (0 - 500)/500 = -1.0; tq=0 -> truth=0 -> divide by 0 yields 0.0
        np.testing.assert_allclose(res, [-1.0, 0.0, 0.0], atol=1e-12)

    def test_good_seed_gate_excludes_phi_and_separates_sign(self):
        # SEED_FIELDS: seedIndex(0), maj(1), loc0(2), loc1(3), phi(4), theta(5), qop(6), t(7),
        #              t_loc0(8), t_loc1(9), t_phi(10), t_theta(11), t_qop(12), t_t(13)
        good_large_dphi = [0, 1, 0., 0., 0.25, 1.5708, -0.00205, 0., 0., 0., 0.0, 1.57082, -0.002, 0.]
        self.assertTrue(ldr.seed_good(good_large_dphi))
        wrong_sign = list(good_large_dphi)
        wrong_sign[6] = +0.00205
        self.assertFalse(ldr.seed_good(wrong_sign))
        bad_theta = list(good_large_dphi)
        bad_theta[5] = 1.5708 + 2e-4
        self.assertFalse(ldr.seed_good(bad_theta))

    def test_subset_files_strictly_excludes_calibration_file0(self):
        fake = dict(
            values=np.array([[0, 1.], [1, 2.], [2, 3.]]),
            fields=['file', 'x'],
            geometry=np.array([10, 20, 30], dtype=np.uint64),
            hits=np.array([1, 2, 3], dtype=np.uint64),
            tracks=np.array([[0, 0.], [1, 1.], [2, 2.]]),
            track_fields=['file', 't'],
            seeds=np.array([[0, 0.], [1, 1.]]),
            forensic_seeds=np.array([[0, 5.], [1, 6.], [2, 7.]]),
            forensic_seed_fields=['file', 's'],
            forensic_events=np.array([[0, 1.], [1, 1.], [2, 1.]]),
            forensic_event_fields=['file', 'e'],
            paths=[Path('file00'), Path('file01'), Path('file02')],
            coverage=[dict(file='f0', file_index=0, skip=0, events=1000),
                      dict(file='f1', file_index=1, skip=0, events=1000),
                      dict(file='f2', file_index=2, skip=0, events=1000)],
        )
        sub = lds.subset_files(fake, {1, 2})
        self.assertEqual({c['file_index'] for c in sub['coverage']}, {1, 2})
        self.assertNotIn(0, sub['values'][:, 0])
        self.assertNotIn(0, sub['forensic_seeds'][:, 0])
        self.assertEqual([p.name for p in sub['paths']], ['file01', 'file02'])

    def test_patch_0049_reverse_applies_cleanly_if_worktree_present(self):
        wt = Path('/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_momentum_v2')
        patch = Path(__file__).resolve().parents[1] / 'patches/0049_latest_seed_and_cov_audit.patch'
        self.assertTrue(patch.is_file())
        if wt.exists():
            diff = subprocess.check_output(['git', '-C', str(wt), 'diff', 'HEAD'], text=True)
            self.assertEqual(diff, patch.read_text())


if __name__ == '__main__':
    unittest.main()
