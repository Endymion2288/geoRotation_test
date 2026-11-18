import math
import unittest
import numpy as np


class LatestValidationFixesTests(unittest.TestCase):
    def test_rotated_plane_cov_loc_and_original_inflation(self):
        infl = [100.0, 100.0, 100.0, 100.0, 1000.0]
        sigma_loc0 = 0.019  # fine y_F
        sigma_loc1 = 0.900  # coarse stereo x_F
        sigma_phi = 0.033
        sigma_theta = 0.00029
        sigma_qop = 0.0001
        cov_loc0 = infl[0] * sigma_loc0 ** 2
        cov_loc1 = infl[1] * sigma_loc1 ** 2
        cov_phi = infl[2] * sigma_phi ** 2
        cov_theta = infl[3] * sigma_theta ** 2
        cov_qop = infl[4] * sigma_qop ** 2
        self.assertAlmostEqual(cov_loc0, 0.0361, places=10)
        self.assertAlmostEqual(cov_loc1, 81.0, places=10)
        self.assertAlmostEqual(cov_phi, 0.1089, places=10)
        self.assertAlmostEqual(cov_theta, 8.41e-6, places=12)
        self.assertAlmostEqual(cov_qop, 1.0e-5, places=12)

    def test_gety_nan_guard_and_strict_weak_ordering(self):
        cx, cy, r = 1200.0, 5.0, 400.0
        for z in [0.0, 799.9, 1200.0, 1600.1, 2440.0]:
            raw_rad = r * r - (z - cx) * (z - cx)
            sqt = math.sqrt(raw_rad) if raw_rad > 0.0 else 0.0
            y = cy - sqt if abs(cy - sqt) < abs(cy + sqt) else cy + sqt
            self.assertTrue(math.isfinite(y))
            if abs(z - cx) > r:
                dz_out = abs(z - cx) - r
                self.assertGreater(dz_out, 0.0)

        scores = [float('nan'), 12.5, 1.1, float('inf'), 3.4, float('nan')]
        def sort_key(val):
            return (0, val) if math.isfinite(val) else (1, 0.0)
        ordered = sorted(scores, key=sort_key)
        self.assertEqual(ordered[:3], [1.1, 3.4, 12.5])
        self.assertFalse(math.isfinite(ordered[-1]))

    def test_local_tangent_eliminates_chord_phi_bias_for_both_charges(self):
        # Circular arc in (z_F, y_F) with radius R = 1e6 mm (~142 GeV) starting at z=0, y=0, dy/dz(0)=0
        R = 1.0e6
        z0, z1, z2 = -10.0, 40.0, 1240.0
        for charge, cy in [(+1, -R), (-1, +R)]:
            cx = 0.0
            def y_true(z):
                sqt = math.sqrt(R * R - (z - cx) ** 2)
                return cy + sqt if cy < 0 else cy - sqt

            def dydz_true(z):
                y = y_true(z)
                return -(z - cx) / (y - cy)

            # Chord S1 -> S2 vs local circle tangent at origin z0
            dydz_chord = (y_true(z2) - y_true(z1)) / (z2 - z1)
            sqt0 = math.sqrt(R * R - (z0 - cx) ** 2)
            y_cand1, y_cand2 = cy - sqt0, cy + sqt0
            y0_est = y_cand1 if abs(y_cand1 - y_true(z1)) < abs(y_cand2 - y_true(z1)) else y_cand2
            dydz_tan = -(z0 - cx) / (y0_est - cy)

            # In rotated ACTS frame (x_A, y_A, z_A) = (z_F, y_F, -x_F) with x_F=0: phi_A = atan2(dydz, 1)
            phi_true = math.atan2(dydz_true(z0), 1.0)
            phi_chord = math.atan2(dydz_chord, 1.0)
            phi_tan = math.atan2(dydz_tan, 1.0)

            self.assertAlmostEqual(phi_tan - phi_true, 0.0, places=12)
            # Chord has ~0.65 mrad opposite-sign error!
            self.assertGreater(abs(phi_chord - phi_true) * 1e3, 0.5)
            if charge > 0:
                self.assertLess(phi_chord - phi_true, 0.0)
            else:
                self.assertGreater(phi_chord - phi_true, 0.0)

    def test_pdg_charge_covers_leptons_and_hadrons(self):
        def pdg_charge(pdg):
            a = abs(pdg)
            if a in (11, 13, 15):
                return -1.0 if pdg > 0 else +1.0
            if a in (211, 321, 2212):
                return +1.0 if pdg > 0 else -1.0
            return 0.0

        self.assertEqual(pdg_charge(11), -1.0)    # e-
        self.assertEqual(pdg_charge(-11), +1.0)   # e+
        self.assertEqual(pdg_charge(13), -1.0)    # mu-
        self.assertEqual(pdg_charge(-13), +1.0)   # mu+
        self.assertEqual(pdg_charge(211), +1.0)   # pi+
        self.assertEqual(pdg_charge(-211), -1.0)  # pi-
        self.assertEqual(pdg_charge(2212), +1.0)  # proton

    def test_two_dipole_strip_field_fit_recovers_exact_five_parameters(self):
        # Synthetic 18 wafers (3 stations x 3 layers x 2 stereo sides +/-20 mrad) in ACTS coords
        # Dipoles in [120, 1120] and [1320, 2320] mm with -B_zA = 1.7e-4 GeV/mm (~0.57 T)
        xs = np.linspace(-50.0, 2550.0, 2601)
        dx = 1.0
        ay = np.where(((xs >= 120.0) & (xs <= 1120.0)) | ((xs >= 1320.0) & (xs <= 2320.0)), 1.7e-4, 0.0)
        izero = int(round((0.0 - xs[0]) / dx))
        jy = np.zeros_like(xs)
        iy = np.zeros_like(xs)
        for i in range(izero + 1, len(xs)):
            jy[i] = jy[i - 1] + 0.5 * (ay[i - 1] + ay[i]) * dx
            iy[i] = iy[i - 1] + 0.5 * (jy[i - 1] + jy[i]) * dx
        for i in range(izero, 0, -1):
            jy[i - 1] = jy[i] - 0.5 * (ay[i - 1] + ay[i]) * dx
            iy[i - 1] = iy[i] - 0.5 * (jy[i - 1] + jy[i]) * dx

        theta_true = np.array([3.2, -0.0015, -7.4, 0.0022, -1.0 / 250.0])  # y0, ty, z0, tz, q/p
        wafers_x = [0, 3, 40, 43, 80, 83, 1200, 1203, 1240, 1243, 1280, 1283, 2400, 2403, 2440, 2443, 2480, 2483]
        H_rows, b_rows = [], []
        for idx, xw in enumerate(wafers_x):
            alpha = 0.02 if (idx % 2 == 0) else -0.02
            e0 = np.array([0.0, math.cos(alpha), math.sin(alpha)])
            iy_w = float(np.interp(xw, xs, iy))
            y_true = theta_true[0] + theta_true[1] * xw + theta_true[4] * iy_w
            z_true = theta_true[2] + theta_true[3] * xw
            u_meas = e0[1] * y_true + e0[2] * z_true
            H_rows.append([e0[1], e0[1] * xw, e0[2], e0[2] * xw, e0[1] * iy_w])
            b_rows.append(u_meas)
        H = np.asarray(H_rows)
        b = np.asarray(b_rows)
        theta_fit = np.linalg.solve(H.T @ H, H.T @ b)
        np.testing.assert_allclose(theta_fit, theta_true, rtol=1e-10, atol=1e-12)

    def test_union_and_s23_dedup_preserves_merged_s1_pair_and_kills_latin_square_fake(self):
        # 5 segments of 6 clusters each (30 clusters total):
        # S1_merged [0..5], S2_em [6..11], S2_ep [12..17], S3_em [18..23], S3_ep [24..29], plus S3_ghost [30..35]
        def mask(rng):
            s = set(rng)
            return s

        s1_m = mask(range(0, 6))
        s2_em, s2_ep = mask(range(6, 12)), mask(range(12, 18))
        s3_em, s3_ep, s3_gh = mask(range(18, 24)), mask(range(24, 30)), mask(range(30, 36))

        candidates = [
            ('real_em', s1_m | s2_em | s3_em, s2_em | s3_em),
            ('real_ep_merged_s1', s1_m | s2_ep | s3_ep, s2_ep | s3_ep),
            ('fake_cross_s23', s1_m | s2_em | s3_ep, s2_em | s3_ep),
            ('fake_latin_ghost', s1_m | s2_ep | s3_gh, s2_ep | s3_gh),
        ]
        selected = []
        used_union = set()
        for name, cset, s23 in candidates:
            kill = False
            if len(cset & used_union) > 6:
                kill = True
            for _, sel_c, sel_s23 in selected:
                if len(cset & sel_c) > 6 or len(s23 & sel_s23) >= 4:
                    kill = True
            if not kill:
                selected.append((name, cset, s23))
                used_union |= cset

        kept_names = [x[0] for x in selected]
        self.assertEqual(kept_names, ['real_em', 'real_ep_merged_s1'])

    def test_neutral_primary_geant4_decay_barcode_resolution(self):
        # K_S (PDG 310, bc=10001 -> cont 1010001) decaying in Geant4 to pi+ (200035) + pi- (200067)
        # with Geant4 continuation vertex 1200035 on pi+
        particles = [
            {'bc': 10001, 'pdg': 310, 'status': 1, 'in_bcs': [], 'z': -1750.0, 'p': 150.0},
            {'bc': 1010001, 'pdg': 310, 'status': 1, 'in_bcs': [10001], 'z': -1200.0, 'p': 150.0},
            {'bc': 200035, 'pdg': 211, 'status': 1, 'in_bcs': [1010001], 'z': -1200.0, 'p': 71.9},
            {'bc': 200067, 'pdg': -211, 'status': 1, 'in_bcs': [1010001], 'z': -1200.0, 'p': 78.1},
            {'bc': 1200035, 'pdg': 211, 'status': 1, 'in_bcs': [200035], 'z': 15.0, 'p': 71.8},
            {'bc': 200999, 'pdg': 11, 'status': 1, 'in_bcs': [200035], 'z': 40.0, 'p': 0.002},  # delta-ray
        ]
        pdg_map = {p['bc']: p['pdg'] for p in particles}

        def pdg_q(pdg):
            if abs(pdg) in (11, 13, 15):
                return -1.0 if pdg > 0 else 1.0
            if abs(pdg) in (211, 321, 2212):
                return 1.0 if pdg > 0 else -1.0
            return 0.0

        has_gen_charged = any(p['bc'] < 200000 and p['status'] == 1 and pdg_q(p['pdg']) != 0.0 for p in particles)
        self.assertFalse(has_gen_charged)

        signals = {}
        cont = {}
        for p in particles:
            bc = p['bc']
            if pdg_q(p['pdg']) == 0.0:
                continue
            from_neutral_primary = any((ibc % 1000000) < 200000 and pdg_q(pdg_map[ibc]) == 0.0 for ibc in p['in_bcs'])
            if bc < 1000000:
                if from_neutral_primary and p['z'] < 50.0:
                    signals[bc] = p['pdg']
                    cont[bc] = [(p['z'], p['p'])]
            else:
                base = bc % 1000000
                if base in signals and signals[base] == p['pdg']:
                    cont[base].append((p['z'], p['p']))

        self.assertEqual(sorted(signals.keys()), [200035, 200067])
        p_st1_piplus = min(p for z, p in cont[200035] if z <= 20.0)
        self.assertAlmostEqual(p_st1_piplus, 71.8, places=6)


if __name__ == '__main__':
    unittest.main()
