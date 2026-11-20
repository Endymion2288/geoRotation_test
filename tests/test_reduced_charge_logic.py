#!/usr/bin/env python3
"""Regression unit tests for Calypso NtupleDumper waveform & reduced VetoNu charge logic.

Covers:
  - Hit selection ((status & ~512) == 0) and +16.0 ns correction for TIMING_SHIFTED (512)
  - Python-side timing presets (NtupleDumperConfig.REDUCED_CHARGE_TIMING_PRESETS) & overrides
  - 2025-2026 4-PMT Preshower Dig1/Dig2 alignment across partial-PMT subsets
  - VetoNu channel-to-index mapping with secondary hits (size >= 3) and swapped order
  - Empty raw waveform (wave_vector) safe fallback to raw_charge (never zeroing)
  - Source-level invariants: ZERO hardcoded run numbers in C++ NtupleDumperAlg.{h,cxx}
"""

from __future__ import annotations

import ast
import math
from pathlib import Path
import unittest

STATUS_THRESHOLD_FAILED = 1 << 0
STATUS_SECONDARY = 1 << 1
STATUS_OVERFLOW = 1 << 2
STATUS_GFIT_FAILED = 1 << 4
STATUS_TIMING_SHIFTED = 1 << 9  # 512
INT_TO_CHARGE = 2.0 / 50.0

NTUPLE_PKG = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1/PhysicsAnalysis/NtupleDumper"
)

# Mirror of NtupleDumperConfig.REDUCED_CHARGE_TIMING_PRESETS
EXPECTED_PRESETS = {
    "data_default": {
        "VetoVetoNuTimeDifference": 33.165,
        "PreshowerVetoNuTimeDifference": 53.625,
        "PreshowerDig2ExtraOffset": 0.0,
    },
    "data_caloNu": {
        "VetoVetoNuTimeDifference": 35.783,
        "PreshowerVetoNuTimeDifference": 53.512,
        "PreshowerDig2ExtraOffset": 0.0,
    },
    "data_2025": {
        "VetoVetoNuTimeDifference": 33.165,
        "PreshowerVetoNuTimeDifference": 13.27,
        "PreshowerDig2ExtraOffset": 3.97,
    },
    "data_2026": {
        "VetoVetoNuTimeDifference": 33.165,
        "PreshowerVetoNuTimeDifference": 13.27,
        "PreshowerDig2ExtraOffset": 3.64,
    },
    "mc": {
        "VetoVetoNuTimeDifference": 27.21,
        "PreshowerVetoNuTimeDifference": 51.41,
        "PreshowerDig2ExtraOffset": 0.0,
    },
}


def load_presets_from_config(cfg_path: Path) -> dict[str, dict[str, float]]:
    """Parse REDUCED_CHARGE_TIMING_PRESETS from NtupleDumperConfig.py without importing Athena."""
    tree = ast.parse(cfg_path.read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "REDUCED_CHARGE_TIMING_PRESETS":
                    return ast.literal_eval(node.value)
    raise AssertionError("REDUCED_CHARGE_TIMING_PRESETS not found in NtupleDumperConfig.py")


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
        if dig2_extra_offset != 0.0 and h["channel"] >= 16:
            t -= dig2_extra_offset
        total += t
        n += 1
    return (total / n) if n > 0 else float("nan")


def compute_reduced_charge(
    timing_cfg: dict[str, float],
    mapped_veto2: list[int],
    mapped_preshower: list[int],
    veto_hits: list[dict],
    ps_hits: list[dict],
    vn_hits: list[dict],
    raw_charge_by_ch: dict[int, float],
    window: float = 35.0,
    mapped_vetonu: list[int] | None = None,
) -> tuple[list[float], list[float], int]:
    """Python mirror of refactored NtupleDumperAlg reduced VetoNu charge calculation."""
    dv = timing_cfg["VetoVetoNuTimeDifference"]
    dps = timing_cfg["PreshowerVetoNuTimeDifference"]
    ddig2 = timing_cfg.get("PreshowerDig2ExtraOffset", 0.0)
    vn_ch = mapped_vetonu if mapped_vetonu and len(mapped_vetonu) >= 2 else [4, 5]

    t_v2 = mean_scintillator_time(mapped_veto2, veto_hits, 0.0) if mapped_veto2 else float("nan")
    t_ps_dig1 = float("nan")
    t_ps_dig2 = float("nan")
    if math.isnan(t_v2) and mapped_preshower:
        ps_dig1 = [ch for ch in mapped_preshower if ch < 16]
        ps_dig2 = [ch for ch in mapped_preshower if ch >= 16]
        if ps_dig1:
            t_ps_dig1 = mean_scintillator_time(ps_dig1, ps_hits, 0.0)
        if math.isnan(t_ps_dig1) and ps_dig2:
            t_ps_dig2 = mean_scintillator_time(ps_dig2, ps_hits, ddig2)

    if not math.isnan(t_v2):
        fb = 1
        exp_t = t_v2 - dv
    elif not math.isnan(t_ps_dig1):
        fb = 2
        exp_t = t_ps_dig1 - dps
    elif not math.isnan(t_ps_dig2):
        fb = 3
        exp_t = t_ps_dig2 - dps
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
            {"channel": 20, "status": 0, "localtime": 500.0},  # Veto1, excluded by m_wave_map["veto2"]
        ]
        t = mean_scintillator_time([22, 23], hits)
        self.assertAlmostEqual(t, 182.0, places=6)

    def test_python_timing_presets_cover_all_periods_and_mc(self):
        cfg_path = NTUPLE_PKG / "python/NtupleDumperConfig.py"
        if not cfg_path.is_file():
            self.skipTest("fix_v1 worktree not present")
        presets = load_presets_from_config(cfg_path)
        self.assertEqual(set(presets.keys()), set(EXPECTED_PRESETS.keys()))
        for name, expected in EXPECTED_PRESETS.items():
            self.assertEqual(presets[name], expected)

    def test_2025_2026_preshower_two_tier_dig1_prefer_dig2_fallback(self):
        cfg = EXPECTED_PRESETS["data_2025"]
        tv = [140.0, 150.0, 180.0, 200.0]
        vn_hits = [
            {"channel": 4, "status": 0, "time_vector": tv, "wave_vector": [100.0, 200.0, 300.0, 500.0]},
            {"channel": 5, "status": 0, "time_vector": tv, "wave_vector": [50.0, 100.0, 150.0, 400.0]},
        ]
        # Saturated Veto2 forces Preshower fallback
        veto_ovf = [
            {"channel": 6, "status": STATUS_OVERFLOW, "localtime": 183.0},
            {"channel": 7, "status": STATUS_OVERFLOW, "localtime": 183.0},
        ]

        # Case 1: Both Dig1 (ch 12,13) and Dig2 (ch 24,25) valid -> MUST prefer Dig1 (fb=2) and ignore Dig2 jitter!
        ps_both = [
            {"channel": 12, "status": 0, "localtime": 163.27},
            {"channel": 13, "status": 0, "localtime": 163.27},
            {"channel": 24, "status": 0, "localtime": 175.0},  # intentionally jittered Dig2
            {"channel": 25, "status": STATUS_TIMING_SHIFTED, "localtime": 159.0},
        ]
        q_red, q_val, fb = compute_reduced_charge(
            timing_cfg=cfg,
            mapped_veto2=[6, 7],
            mapped_preshower=[12, 13, 24, 25],
            veto_hits=veto_ovf,
            ps_hits=ps_both,
            vn_hits=vn_hits,
            raw_charge_by_ch={4: 99.0, 5: 88.0},
        )
        self.assertEqual(fb, 2)
        # exp_t = 163.27 - 13.27 = 150.0 -> includes 140, 150, 180 (< 185), excludes 200
        self.assertAlmostEqual(q_red[0], (100.0 + 200.0 + 300.0) * INT_TO_CHARGE, places=5)
        self.assertAlmostEqual(q_val[0], q_red[0], places=5)

        # Case 2: Dig1 fails (GFIT_FAILED like 2025 p0013 or shower saturation), Dig2 valid (incl. status==512) -> fb=3 with alignment!
        ps_dig1_fail = [
            {"channel": 12, "status": STATUS_GFIT_FAILED, "localtime": 0.0},
            {"channel": 13, "status": STATUS_OVERFLOW, "localtime": 180.0},
            {"channel": 24, "status": 0, "localtime": 167.24},
            {"channel": 25, "status": STATUS_TIMING_SHIFTED, "localtime": 151.24},
        ]
        q_red2, q_val2, fb2 = compute_reduced_charge(
            timing_cfg=cfg,
            mapped_veto2=[6, 7],
            mapped_preshower=[12, 13, 24, 25],
            veto_hits=veto_ovf,
            ps_hits=ps_dig1_fail,
            vn_hits=vn_hits,
            raw_charge_by_ch={4: 99.0, 5: 88.0},
        )
        self.assertEqual(fb2, 3)
        # Dig2 aligned: (167.24 - 3.97) == 163.27 and (151.24 + 16.0 - 3.97) == 163.27 -> exp_t = 150.0
        self.assertAlmostEqual(q_red2[0], (100.0 + 200.0 + 300.0) * INT_TO_CHARGE, places=5)
        self.assertAlmostEqual(q_val2[0], q_red2[0], places=5)

        # Case 3: Both Dig1 and Dig2 fail -> fb=0, preserve raw_charge fallback
        ps_all_fail = [
            {"channel": 12, "status": STATUS_GFIT_FAILED, "localtime": 0.0},
            {"channel": 13, "status": STATUS_GFIT_FAILED, "localtime": 0.0},
            {"channel": 24, "status": STATUS_THRESHOLD_FAILED, "localtime": 0.0},
            {"channel": 25, "status": STATUS_OVERFLOW, "localtime": 0.0},
        ]
        q_red3, q_val3, fb3 = compute_reduced_charge(
            timing_cfg=cfg,
            mapped_veto2=[6, 7],
            mapped_preshower=[12, 13, 24, 25],
            veto_hits=veto_ovf,
            ps_hits=ps_all_fail,
            vn_hits=vn_hits,
            raw_charge_by_ch={4: 99.0, 5: 88.0},
        )
        self.assertEqual(fb3, 0)
        self.assertAlmostEqual(q_red3[0], 99.0, places=5)
        self.assertAlmostEqual(q_val3[0], 0.0, places=5)

    def test_vetonu_channel_indexing_secondary_and_swapped_order(self):
        tv = [140.0, 150.0, 170.0, 200.0]
        vn_hits = [
            {"channel": 5, "status": 0, "time_vector": tv, "wave_vector": [10.0, 20.0, 30.0, 400.0]},
            {"channel": 5, "status": STATUS_SECONDARY, "time_vector": tv, "wave_vector": [999.0] * 4},
            {"channel": 4, "status": STATUS_SECONDARY, "time_vector": tv, "wave_vector": [888.0] * 4},
            {"channel": 4, "status": 0, "time_vector": tv, "wave_vector": [100.0, 200.0, 300.0, 500.0]},
        ]
        veto_hits = [{"channel": 22, "status": STATUS_TIMING_SHIFTED, "localtime": 150.0 + 35.783 - 16.0}]
        q_red, q_val, fb = compute_reduced_charge(
            timing_cfg=EXPECTED_PRESETS["data_caloNu"],
            mapped_veto2=[22, 23],
            mapped_preshower=[12, 13],
            veto_hits=veto_hits,
            ps_hits=[],
            vn_hits=vn_hits,
            raw_charge_by_ch={4: 99.0, 5: 88.0},
        )
        self.assertEqual(fb, 1)
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
            timing_cfg=EXPECTED_PRESETS["data_default"],
            mapped_veto2=[6, 7],
            mapped_preshower=[12, 13],
            veto_hits=veto_hits,
            ps_hits=[],
            vn_hits=vn_hits,
            raw_charge_by_ch={4: 123.4, 5: 56.7},
        )
        self.assertEqual(fb, 1)
        self.assertAlmostEqual(q_red[0], 123.4, places=4)
        self.assertAlmostEqual(q_val[0], 0.0, places=4)
        self.assertAlmostEqual(q_red[1], 1000.0 * INT_TO_CHARGE, places=4)

    def test_calypso_fix_v1_zero_hardcoded_runs_and_source_invariants(self):
        if not NTUPLE_PKG.is_dir():
            self.skipTest("fix_v1 worktree not present")
        cxx = (NTUPLE_PKG / "src/NtupleDumperAlg.cxx").read_text()
        hdr = (NTUPLE_PKG / "src/NtupleDumperAlg.h").read_text()
        cfg = (NTUPLE_PKG / "python/NtupleDumperConfig.py").read_text()
        maker = (NTUPLE_PKG / "scripts/faser_ntuple_maker.py").read_text()

        # Zero hardcoded run-number logic in C++
        self.assertNotIn("resolveReducedChargeTiming", cxx)
        self.assertNotIn("resolveReducedChargeTiming", hdr)
        for forbidden_run in ("15845", "16928", "13847", "18193", "18194"):
            self.assertNotIn(forbidden_run, cxx, f"Hardcoded run {forbidden_run} found in NtupleDumperAlg.cxx")
            self.assertNotIn(forbidden_run, hdr, f"Hardcoded run {forbidden_run} found in NtupleDumperAlg.h")

        # Clean status/cable-map driven C++ invariants
        self.assertIn("WaveformDefs.h", cxx)
        self.assertIn("TIMING_SHIFTED", cxx)
        self.assertIn("t += 16.0;", cxx)
        self.assertIn("(status & ~kTimingShiftedBit) != 0U", cxx)
        self.assertIn("VetoNu_reduced_charge_fallback", cxx)
        self.assertIn("kDig2ChannelStart", cxx)
        self.assertIn("preshowerDig1Channels", cxx)
        self.assertIn("preshowerDig2Channels", cxx)
        self.assertIn("m_vetonu_reduced_charge_fallback = 2;", cxx)
        self.assertIn("m_vetonu_reduced_charge_fallback = 3;", cxx)
        self.assertNotIn("m_vetonu_reduced_charge[idx] = 0;", cxx)
        self.assertIn("RotatedPhysicalAngles", hdr)
        self.assertIn("PreshowerDig2ExtraOffset", hdr)

        # Python-side configurable timing presets and CLI flags
        self.assertIn("REDUCED_CHARGE_TIMING_PRESETS", cfg)
        self.assertIn("getReducedChargeTimingKwargs", cfg)
        self.assertIn("--timingPreset", maker)
        self.assertIn("--vetoVetoNuTimeDiff", maker)
        self.assertIn("--preshowerVetoNuTimeDiff", maker)
        self.assertIn("--preshowerDig2ExtraOffset", maker)


if __name__ == "__main__":
    unittest.main()
