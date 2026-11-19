#!/usr/bin/env python3
"""Regression unit tests for Calypso NtupleDumper waveform & reduced VetoNu charge logic.

Covers:
  - Hit selection ((status & ~512) == 0) and +16.0 ns correction for TIMING_SHIFTED (512)
  - Automatic 2022-2026 data and MC timing lookup + explicit override precedence
  - 2025-2026 4-PMT Preshower Dig1/Dig2 alignment across partial-PMT subsets
  - VetoNu channel-to-index mapping with secondary hits (size >= 3) and swapped order
  - Empty raw waveform (wave_vector) safe fallback to raw_charge (never zeroing)
  - Source-level invariants in geoRotation_calypso_latest_fix_v1/PhysicsAnalysis/NtupleDumper
"""

from __future__ import annotations

import math
from pathlib import Path
import unittest

STATUS_THRESHOLD_FAILED = 1 << 0
STATUS_SECONDARY = 1 << 1
STATUS_OVERFLOW = 1 << 2
STATUS_GFIT_FAILED = 1 << 4
STATUS_TIMING_SHIFTED = 1 << 9  # 512
INT_TO_CHARGE = 2.0 / 50.0


def resolve_timing(
    is_mc: bool,
    run_number: int,
    mapped_veto2: list[int],
    override_veto: float = -999.0,
    override_ps: float = -999.0,
    override_ps_dig2: float = -999.0,
) -> tuple[float, float, float]:
    """Python mirror of NtupleDumperAlg::resolveReducedChargeTiming."""
    if is_mc:
        dv, dps, ddig2 = 27.21, 51.41, 0.0
    else:
        veto2_dig2 = (15845 <= run_number < 16928) or (run_number == 18193) or (22 in mapped_veto2) or (23 in mapped_veto2)
        if veto2_dig2:
            dv = 35.783
        elif run_number < 13847:
            dv = 33.055
        else:
            dv = 33.165

        if run_number < 13847:
            dps, ddig2 = 53.490, 0.0
        elif run_number < 15845:
            dps, ddig2 = 53.550, 0.0
        elif run_number < 16928:
            dps, ddig2 = 53.512, 0.0
        elif run_number < 18193:
            dps, ddig2 = 53.550, 0.0
        elif run_number == 18193:
            dps, ddig2 = 25.280, 3.970
        else:
            dps, ddig2 = 13.280, 3.970

    return (
        override_veto if override_veto > -900.0 else dv,
        override_ps if override_ps > -900.0 else dps,
        override_ps_dig2 if override_ps_dig2 > -900.0 else ddig2,
    )


def mean_scintillator_time(
    channels: list[int],
    hits: list[dict],
    dig2_extra_offset: float = 0.0,
) -> float:
    """Python mirror of NtupleDumperAlg::get_mean_scintillator_time."""
    total = 0.0
    n = 0
    for h in hits:
        if h["channel"] not in channels:
            continue
        st = h["status"]
        if (st & ~STATUS_TIMING_SHIFTED) != 0:
            continue
        t = float(h["localtime"])
        if st & STATUS_TIMING_SHIFTED:
            t += 16.0
        if h["channel"] >= 16:
            t -= dig2_extra_offset
        total += t
        n += 1
    return (total / n) if n > 0 else float("nan")


def compute_reduced_charge(
    is_mc: bool,
    run_number: int,
    mapped_veto2: list[int],
    veto_hits: list[dict],
    ps_hits: list[dict],
    vn_hits: list[dict],
    raw_charge_by_ch: dict[int, float],
    window: float = 35.0,
) -> tuple[list[float], list[float], int]:
    """Python mirror of NtupleDumperAlg reduced VetoNu charge calculation."""
    dv, dps, ddig2 = resolve_timing(is_mc, run_number, mapped_veto2)
    veto2_ch = list(dict.fromkeys(list(mapped_veto2) + [6, 7, 22, 23]))
    ps_ch = [12, 13, 24, 25]
    vn_ch = [4, 5]

    t_v2 = mean_scintillator_time(veto2_ch, veto_hits, 0.0)
    t_ps = mean_scintillator_time(ps_ch, ps_hits, ddig2)

    if not math.isnan(t_v2):
        fb = 1
        exp_t = t_v2 - dv
    elif not math.isnan(t_ps):
        fb = 2
        exp_t = t_ps - dps
    else:
        fb = 0
        exp_t = float("nan")

    q_red = [raw_charge_by_ch.get(vn_ch[0], 0.0), raw_charge_by_ch.get(vn_ch[1], 0.0)]
    q_val = [0.0, 0.0]
    if fb == 0 or math.isnan(exp_t):
        return q_red, q_val, fb

    good = [False, False]
    for h in vn_hits:
        ch = h["channel"]
        if ch == vn_ch[0]:
            idx = 0
        elif ch == vn_ch[1]:
            idx = 1
        else:
            continue
        st = h["status"]
        if ((st & STATUS_SECONDARY) == 0) or (not good[idx]):
            tv = h["time_vector"]
            wv = h["wave_vector"]
            if not tv:
                continue
            good[idx] = True
            is_valid = (st & ~STATUS_TIMING_SHIFTED) == 0
            t_shift = 16.0 if (st & STATUS_TIMING_SHIFTED) else 0.0
            s = sum(w for t, w in zip(tv, wv) if (t + t_shift - exp_t) < window)
            q = s * INT_TO_CHARGE
            q_red[idx] = q
            q_val[idx] = q if is_valid else 0.0
    return q_red, q_val, fb


class TestReducedChargeLogic(unittest.TestCase):
    def test_status512_selection_and_16ns_correction(self):
        hits = [
            {"channel": 22, "status": 0, "localtime": 182.0},
            {"channel": 23, "status": STATUS_TIMING_SHIFTED, "localtime": 166.0},
            {"channel": 22, "status": STATUS_TIMING_SHIFTED | STATUS_OVERFLOW, "localtime": 100.0},
            {"channel": 20, "status": 0, "localtime": 500.0},  # Veto1, must be ignored
        ]
        t = mean_scintillator_time([6, 7, 22, 23], hits)
        self.assertAlmostEqual(t, 182.0, places=6)

    def test_period_timing_boundaries_2022_to_2026_and_mc(self):
        self.assertEqual(resolve_timing(True, 100044, [6, 7]), (27.21, 51.41, 0.0))
        self.assertEqual(resolve_timing(False, 8090, [6, 7]), (33.055, 53.490, 0.0))
        self.assertEqual(resolve_timing(False, 11703, [6, 7]), (33.055, 53.490, 0.0))
        self.assertEqual(resolve_timing(False, 15387, [6, 7]), (33.165, 53.550, 0.0))
        # Early CaloNu (15821-15844): Veto still on Dig1 ch 6,7
        self.assertEqual(resolve_timing(False, 15838, [6, 7]), (33.165, 53.550, 0.0))
        # Main CaloNu (15845-16927): Veto on Dig2 ch 22,23
        self.assertEqual(resolve_timing(False, 15849, [22, 23]), (35.783, 53.512, 0.0))
        self.assertEqual(resolve_timing(False, 16923, [22, 23]), (35.783, 53.512, 0.0))
        # Post-CaloNu
        self.assertEqual(resolve_timing(False, 16929, [6, 7]), (33.165, 53.550, 0.0))
        # 2025-2026 4-PMT Preshower
        self.assertEqual(resolve_timing(False, 20012, [6, 7]), (33.165, 13.280, 3.970))
        self.assertEqual(resolve_timing(False, 22404, [6, 7]), (33.165, 13.280, 3.970))
        # Explicit overrides win over auto
        self.assertEqual(resolve_timing(False, 22404, [6, 7], override_veto=34.0, override_ps=14.5), (34.0, 14.5, 3.970))

    def test_2025_2026_preshower_dig1_dig2_alignment(self):
        # True t_VN = 150.0 ns -> ch12/13 at 163.28 ns, ch24/25 at 167.25 ns (or 151.25 ns with status=512)
        _, dps, ddig2 = resolve_timing(False, 22404, [6, 7])
        h_top = [{"channel": 12, "status": 0, "localtime": 163.28}, {"channel": 13, "status": 0, "localtime": 163.28}]
        h_bot = [
            {"channel": 24, "status": 0, "localtime": 167.25},
            {"channel": 25, "status": STATUS_TIMING_SHIFTED, "localtime": 151.25},
        ]
        # Top only, bottom only (like 2025 runs 19480-20080 when PS0,1 failed GFIT), and all 4 must give same exp_t = 150.0
        for subset in (h_top, h_bot, h_top + h_bot):
            t_ps = mean_scintillator_time([12, 13, 24, 25], subset, ddig2)
            self.assertAlmostEqual(t_ps - dps, 150.0, places=3)

    def test_vetonu_channel_indexing_secondary_and_swapped_order(self):
        # ch 5 appears before ch 4, plus two secondary hits (total size=4, would overflow old [2] array!)
        tv = [140.0, 150.0, 170.0, 200.0]
        vn_hits = [
            {"channel": 5, "status": 0, "time_vector": tv, "wave_vector": [10.0, 20.0, 30.0, 400.0]},
            {"channel": 5, "status": STATUS_SECONDARY, "time_vector": tv, "wave_vector": [999.0, 999.0, 999.0, 999.0]},
            {"channel": 4, "status": STATUS_SECONDARY, "time_vector": tv, "wave_vector": [888.0, 888.0, 888.0, 888.0]},
            {"channel": 4, "status": 0, "time_vector": tv, "wave_vector": [100.0, 200.0, 300.0, 500.0]},
        ]
        veto_hits = [{"channel": 22, "status": STATUS_TIMING_SHIFTED, "localtime": 150.0 + 35.783 - 16.0}]
        q_red, q_val, fb = compute_reduced_charge(
            is_mc=False,
            run_number=15849,
            mapped_veto2=[22, 23],
            veto_hits=veto_hits,
            ps_hits=[],
            vn_hits=vn_hits,
            raw_charge_by_ch={4: 99.0, 5: 88.0},
        )
        self.assertEqual(fb, 1)
        # exp_t = 150.0, cutoff = 185.0 -> includes 140, 150, 170 (excludes 200)
        self.assertAlmostEqual(q_red[0], (100.0 + 200.0 + 300.0) * INT_TO_CHARGE, places=5)
        self.assertAlmostEqual(q_red[1], (10.0 + 20.0 + 30.0) * INT_TO_CHARGE, places=5)
        self.assertAlmostEqual(q_val[0], q_red[0], places=5)
        self.assertAlmostEqual(q_val[1], q_red[1], places=5)

    def test_empty_wave_vector_preserves_raw_charge_fallback(self):
        veto_hits = [{"channel": 6, "status": 0, "localtime": 183.165}]
        vn_hits = [
            {"channel": 4, "status": 0, "time_vector": [], "wave_vector": []},
            {"channel": 5, "status": 0, "time_vector": [150.0], "wave_vector": [1000.0]},
        ]
        q_red, q_val, fb = compute_reduced_charge(
            is_mc=False,
            run_number=16929,
            mapped_veto2=[6, 7],
            veto_hits=veto_hits,
            ps_hits=[],
            vn_hits=vn_hits,
            raw_charge_by_ch={4: 123.4, 5: 56.7},
        )
        self.assertEqual(fb, 1)
        self.assertAlmostEqual(q_red[0], 123.4, places=4)
        self.assertAlmostEqual(q_val[0], 0.0, places=4)
        self.assertAlmostEqual(q_red[1], 1000.0 * INT_TO_CHARGE, places=4)

    def test_calypso_fix_v1_source_invariants(self):
        pkg = Path("/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/PhysicsAnalysis/NtupleDumper")
        if not pkg.is_dir():
            self.skipTest("fix_v1 worktree not present")
        cxx = (pkg / "src/NtupleDumperAlg.cxx").read_text()
        hdr = (pkg / "src/NtupleDumperAlg.h").read_text()
        maker = (pkg / "scripts/faser_ntuple_maker.py").read_text()

        self.assertIn("WaveformDefs.h", cxx)
        self.assertIn("TIMING_SHIFTED", cxx)
        self.assertIn("t += 16.0;", cxx)
        self.assertIn("(status & ~kTimingShiftedBit) != 0U", cxx)
        self.assertIn("resolveReducedChargeTiming", cxx)
        self.assertIn("VetoNu_reduced_charge_fallback", cxx)
        self.assertNotIn("m_vetonu_reduced_charge[idx] = 0;", cxx)
        self.assertIn("RotatedPhysicalAngles", hdr)
        self.assertIn("--vetoVetoNuTimeDiff", maker)
        self.assertNotIn("VetoVetoNuTimeDifference=33.165", maker)


if __name__ == "__main__":
    unittest.main()
