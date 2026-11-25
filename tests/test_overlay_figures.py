#!/usr/bin/env python3
"""Unit tests for Stage 08 overlay plotting pipeline (both --mode baseline and --mode full)."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.overlay_figures import FIGURE_BUILDERS, generate_all_figures, slice_1d_range
from analysis.overlay_plot_style import SIGNAL_LABEL, choose_delta_ylim, choose_eff_ylim, choose_ratio_ylim, safe_asymmetric_yerr
from analysis.overlay_tier2_figures import ARM_ORDER, DELTA_PAIRS, generate_tier2_figures


def _make_eff_curve(edges: np.ndarray, base_eff: float = 0.935, dip_zero: bool = True) -> dict:
    centers = 0.5 * (edges[:-1] + edges[1:])
    hw = 0.5 * (edges[1:] - edges[:-1])
    den = np.full(len(centers), 1000, dtype=int)
    eff = np.full(len(centers), base_eff, dtype=float)
    if dip_zero:
        eff = eff - 0.06 * np.exp(-np.abs(centers) / 4.0)
    eff = np.clip(eff, 0.72, 0.995)
    passed = np.round(eff * den).astype(int)
    err = np.sqrt(eff * (1.0 - eff) / den)
    return {
        "name": "synthetic",
        "edges": edges.tolist(),
        "centers": centers.tolist(),
        "half_widths": hw.tolist(),
        "underflow_den": 0,
        "overflow_den": 12,
        "invalid_den": 0,
        "den": den.tolist(),
        "pass": passed.tolist(),
        "eff": eff.tolist(),
        "err_lo": err.tolist(),
        "err_hi": err.tolist(),
    }


def _make_ratio_curve(edges: np.ndarray, val: float = 1.01) -> dict:
    centers = 0.5 * (edges[:-1] + edges[1:])
    hw = 0.5 * (edges[1:] - edges[:-1])
    r = np.full(len(centers), val, dtype=float)
    err = np.full(len(centers), 0.005, dtype=float)
    return {
        "edges": edges.tolist(),
        "centers": centers.tolist(),
        "half_widths": hw.tolist(),
        "ratio": r.tolist(),
        "err_lo": err.tolist(),
        "err_hi": err.tolist(),
    }


def _make_paired_curve(edges: np.ndarray, delta: float = 0.002) -> dict:
    centers = 0.5 * (edges[:-1] + edges[1:])
    hw = 0.5 * (edges[1:] - edges[:-1])
    n = len(centers)
    return {
        "edges": edges.tolist(),
        "centers": centers.tolist(),
        "half_widths": hw.tolist(),
        "n_den": [1000] * n,
        "n_both_pass": [930] * n,
        "n_gain": [4] * n,
        "n_loss": [2] * n,
        "n_both_fail": [64] * n,
        "delta": [delta] * n,
        "delta_err": [0.0024] * n,
    }


def _make_synthetic_inputs(include_latest: bool = True) -> tuple[dict, dict]:
    curve_specs = {
        "dR1_wide": np.linspace(0.0, 100.0, 201),
        "dR1_zoom": np.linspace(0.0, 20.0, 41),
        "dR0_wide": np.linspace(0.0, 100.0, 201),
        "dR0_zoom": np.linspace(0.0, 20.0, 41),
        "abs_dY0_wide": np.linspace(0.0, 100.0, 101),
        "abs_dY0_zoom": np.linspace(0.0, 20.0, 41),
        "abs_dY1_wide": np.linspace(0.0, 100.0, 101),
        "abs_dY1_zoom": np.linspace(0.0, 20.0, 41),
        "abs_dX0_wide": np.linspace(0.0, 100.0, 101),
        "dX0_signed": np.linspace(-100.0, 100.0, 201),
        "abs_dX1_wide": np.linspace(0.0, 100.0, 101),
        "dX1_signed": np.linspace(-100.0, 100.0, 201),
        "dY1_signed": np.linspace(-100.0, 100.0, 201),
        "dY0_signed": np.linspace(-100.0, 100.0, 201),
    }
    versions = ["baseline_may30", "baseline_jul8"] + (["latest_georot"] if include_latest else [])
    datasets = ["All3Data", "10004X", "008749", "010698", "Pos_008749-Neg_010698", "008749+010698"]
    pairs = ["baseline_jul8_minus_baseline_may30"] + (
        ["latest_georot_minus_baseline_jul8", "latest_georot_minus_baseline_may30"] if include_latest else []
    )

    summary = {
        "schema_version": "1.0",
        "mode": "full" if include_latest else "baseline_only",
        "primary_version": "latest_georot" if include_latest else "baseline_jul8",
        "reference_version": "baseline_jul8",
        "versions": versions,
        "datasets": datasets,
        "provenance": {
            "frame": "FASER",
            "cond_tag": "OFLCOND-FASER-05",
            "geom_tag": "TI12MC04 / FASERNU-04",
            "matched_bin_pairs": {"008749": 81, "010698": 80, "Pos_008749-Neg_010698": 81, "10004X": 80, "total": 322},
            "og_identity_checks_passed": 6339599,
            "og_identity_checks_failed": 0,
        },
        "global_summary": {},
        "paired_global_deltas": {},
        "transition_matrices": {},
        "curves": {},
        "mc_data_ratios": {},
        "paired_curves": {},
        "multiplicity_curves": {},
        "momentum": {"mc_data_pz0_ratio": {}},
        "stage_of_loss": {},
        "large_dX_diagnostics": {},
        "momentum_reweighting": {},
    }

    pz_edges = np.geomspace(1.0, 2.0e7, 101)
    pz_cen = np.sqrt(pz_edges[:-1] * pz_edges[1:])
    res_edges = np.linspace(-0.5, 0.5, 101)
    res_cen = 0.5 * (res_edges[:-1] + res_edges[1:])
    p_rw_edges = np.geomspace(50.0, 5000.0, 41)
    p_rw_cen = np.sqrt(p_rw_edges[:-1] * p_rw_edges[1:])
    f_da = np.exp(-np.log(p_rw_cen / 150.0) ** 2 / 2.0)
    f_da /= f_da.sum()
    f_mc = np.exp(-np.log(p_rw_cen / 900.0) ** 2 / 1.5)
    f_mc /= f_mc.sum()
    summary["momentum_reweighting"]["og_p0_bins_gev"] = {
        "edges": p_rw_edges.tolist(),
        "centers": p_rw_cen.tolist(),
        "data_frac": f_da.tolist(),
        "mc_frac": f_mc.tolist(),
        "weights": (f_da / f_mc).tolist(),
    }

    for v_idx, ver in enumerate(versions):
        v_shift = 0.002 * v_idx
        summary["global_summary"][ver] = {}
        summary["curves"][ver] = {}
        summary["mc_data_ratios"][ver] = {}
        summary["multiplicity_curves"][ver] = {}
        summary["momentum"][ver] = {}
        summary["stage_of_loss"][ver] = {}
        summary["large_dX_diagnostics"][ver] = {}

        for ds in datasets:
            base = (0.940 if ds == "10004X" else 0.930) + v_shift
            eff_rec = {"pass": 93000, "den": 100000, "eff": base, "err_lo": 0.0008, "err_hi": 0.0008}
            summary["global_summary"][ver][ds] = {
                "n_files": 80 if ds == "10004X" else 242,
                "n_events": 2852041 if ds == "10004X" else 3487558,
                "n_og2": 2852041 if ds == "10004X" else 3487558,
                "lt_counts": {"0": 100, "1": 6000, "2": 93500, "3": 350, "4": 45, "5plus": 5},
                "lt2_fail_c2R": 500,
                "eff_pawan_c2R": eff_rec,
                "eff_lt2": eff_rec,
                "eff_ge2": eff_rec,
                "eff_matched": eff_rec,
                "regimes": {
                    "small_dR1_lt5": {"eff_pawan_c2R": eff_rec, "eff_lt2": eff_rec, "eff_ge2": eff_rec, "eff_matched": eff_rec},
                    "plateau_dR1_30_70": {"eff_pawan_c2R": eff_rec, "eff_lt2": eff_rec, "eff_ge2": eff_rec, "eff_matched": eff_rec},
                    "large_abs_dX0_70_100": {"eff_pawan_c2R": eff_rec, "eff_lt2": eff_rec, "eff_ge2": eff_rec, "eff_matched": eff_rec},
                },
            }

        for ds in ["All3Data", "10004X"]:
            base = (0.940 if ds == "10004X" else 0.930) + v_shift
            summary["curves"][ver][ds] = {}
            for ck, ed in curve_specs.items():
                summary["curves"][ver][ds][ck] = {
                    "eff_pawan_c2R": _make_eff_curve(ed, base),
                    "eff_lt2": _make_eff_curve(ed, base + 0.004),
                    "eff_ge2": _make_eff_curve(ed, base + 0.012),
                    "eff_matched": _make_eff_curve(ed, base - 0.002),
                }

            summary["multiplicity_curves"][ver][ds] = {}
            for ck in ["dR1_wide", "dR1_zoom", "abs_dX0_wide"]:
                ed = curve_specs[ck]
                nb = len(ed) - 1
                summary["multiplicity_curves"][ver][ds][ck] = {
                    "edges": ed.tolist(),
                    "centers": (0.5 * (ed[:-1] + ed[1:])).tolist(),
                    "half_widths": (0.5 * (ed[1:] - ed[:-1])).tolist(),
                    "den": [1000] * nb,
                    "p_lt0": [0.002] * nb,
                    "p_lt1": [0.055] * nb,
                    "p_lt2_pass_c2R": [0.930] * nb,
                    "p_lt2_fail_c2R": [0.005] * nb,
                    "p_lt_ge3": [0.008] * nb,
                }

            # Momentum spectra
            peak = 1.0e6 if ds == "10004X" else 1.5e5
            dist = np.exp(-0.5 * ((np.log10(pz_cen) - np.log10(peak)) / 0.45) ** 2)
            dist = np.where(pz_cen >= 4.0e4, dist, 1e-5)
            dist /= dist.sum()
            counts = np.round(dist * 500000).astype(int)
            res_dist = np.exp(-0.5 * (res_cen / 0.035) ** 2) + 0.01 * np.exp(-0.5 * (res_cen / 0.2) ** 2)
            res_dist /= res_dist.sum()
            summary["momentum"][ver][ds] = {
                "track_pz0": {
                    "edges": pz_edges.tolist(), "centers": pz_cen.tolist(),
                    "counts": counts.tolist(), "norm": dist.tolist(), "norm_err": (dist * 0.01).tolist(),
                    "total_tracks": int(counts.sum()), "underflow": 10, "overflow": 25,
                },
                "og_track_pz0": {
                    "edges": pz_edges.tolist(), "centers": pz_cen.tolist(),
                    "counts": counts.tolist(), "norm": dist.tolist(), "norm_err": (dist * 0.01).tolist(),
                    "total_tracks": int(counts.sum()), "underflow": 0, "overflow": 0,
                },
                "pz0_rel_res": {
                    "edges": res_edges.tolist(), "centers": res_cen.tolist(),
                    "counts": np.round(res_dist * 400000).astype(int).tolist(), "norm": res_dist.tolist(),
                    "total_matched": 400000, "mean": 0.0004, "std": 0.042, "median": 0.0001, "iqr68": 0.031,
                    "underflow": 50, "overflow": 60,
                },
            }

            # Stage of loss
            sl_counts = {
                "n_den": 100000, "pass_c2R": 93500, "lt2_fail_c2R": 450, "lt_ge3": 350, "lt0": 150,
                "lt1_total": 5550, "lt1_miss_seg_st1": 800, "lt1_miss_seg_st2": 900, "lt1_miss_seg_st3": 850,
                "lt1_miss_seg_multi": 600, "lt1_killed_by_gb": 1100, "lt1_all3gb_lost_in_ckf": 1200, "lt1_unmatched_kept": 100,
            }
            dx_e = np.linspace(0.0, 100.0, 21)
            summary["stage_of_loss"][ver][ds] = {
                "regimes": {
                    "global": sl_counts,
                    "small_dR1_lt5": sl_counts,
                    "plateau_dR1_30_70": sl_counts,
                    "large_abs_dX0_70_100": sl_counts,
                },
                "plateau_vs_nSeg": _make_eff_curve(np.arange(0, 14, dtype=float), base, dip_zero=False),
                "plateau_vs_nGBSeg": _make_eff_curve(np.arange(0, 14, dtype=float), base, dip_zero=False),
                "plateau_vs_min_og_p0": _make_eff_curve(np.geomspace(50.0, 5000.0, 21), base, dip_zero=False),
                "plateau_vs_max_og_R0": _make_eff_curve(np.linspace(0.0, 95.0, 20), base, dip_zero=False),
                "abs_dX0_stage_fractions": {
                    "edges": dx_e.tolist(),
                    "centers": (0.5 * (dx_e[:-1] + dx_e[1:])).tolist(),
                    "den": [5000] * 20,
                    "frac_pass_c2R": np.linspace(0.945, 0.900, 20).tolist(),
                    "frac_lt1_miss_seg": np.linspace(0.025, 0.055, 20).tolist(),
                    "frac_lt1_killed_gb": np.linspace(0.010, 0.015, 20).tolist(),
                    "frac_lt1_lost_ckf": np.linspace(0.010, 0.018, 20).tolist(),
                    "frac_lt0": [0.002] * 20,
                    "frac_lt_ge3": [0.004] * 20,
                    "frac_lt2_fail_c2R": [0.004] * 20,
                },
            }

            # Large |dX| diagnostics
            e2d = np.linspace(-100.0, 100.0, 41)
            xx, yy = np.meshgrid(0.5 * (e2d[:-1] + e2d[1:]), 0.5 * (e2d[:-1] + e2d[1:]))
            rr = np.hypot(xx, yy)
            den2d = np.where(rr < 95.0, 400, 0).astype(int)
            pl2d = np.where(rr < 95.0, 0.025 + 0.06 * (np.abs(xx) / 100.0) ** 2, 0.0)
            lost2d = np.round(den2d * pl2d).astype(int)
            summary["large_dX_diagnostics"][ver][ds] = {
                "eff_abs_dX0_by_max_abs_og_x0": {
                    "0_33": _make_eff_curve(dx_e, base + 0.008, dip_zero=False),
                    "33_66": _make_eff_curve(dx_e, base, dip_zero=False),
                    "66_95": _make_eff_curve(dx_e, base - 0.025, dip_zero=False),
                },
                "eff_abs_dX0_by_max_og_R0": {
                    "0_60": _make_eff_curve(dx_e, base + 0.005, dip_zero=False),
                    "60_80": _make_eff_curve(dx_e, base, dip_zero=False),
                    "80_95": _make_eff_curve(dx_e, base - 0.020, dip_zero=False),
                },
                "loss_map_st1": {"x_edges": e2d.tolist(), "y_edges": e2d.tolist(), "den": den2d.tolist(), "lost": lost2d.tolist(), "p_loss": pl2d.tolist()},
                "loss_map_st3": {"x_edges": e2d.tolist(), "y_edges": e2d.tolist(), "den": den2d.tolist(), "lost": lost2d.tolist(), "p_loss": pl2d.tolist()},
            }

        for ck, ed in curve_specs.items():
            summary["mc_data_ratios"][ver][ck] = {
                "eff_pawan_c2R": _make_ratio_curve(ed, 1.010),
                "eff_lt2": _make_ratio_curve(ed, 1.008),
            }

        mc_n = np.array(summary["momentum"][ver]["10004X"]["track_pz0"]["norm"])
        da_n = np.array(summary["momentum"][ver]["All3Data"]["track_pz0"]["norm"])
        summary["momentum"]["mc_data_pz0_ratio"][ver] = {
            "edges": pz_edges.tolist(),
            "centers": pz_cen.tolist(),
            "ratio": (mc_n / np.maximum(da_n, 1e-12)).tolist(),
            "err": (0.05 * np.ones_like(mc_n)).tolist(),
        }

        # Momentum reweighting per version
        ed_dr1 = curve_specs["dR1_wide"]
        nb = len(ed_dr1) - 1
        summary["momentum_reweighting"][ver] = {
            "n_events": 2852041,
            "sumw": 2852041.0,
            "sumw2": 8500000.0,
            "n_eff": 956957.0,
            "closure_mc_reweighted_frac": f_da.tolist(),
            "curves": {
                "dR1_wide": {
                    "edges": ed_dr1.tolist(),
                    "centers": (0.5 * (ed_dr1[:-1] + ed_dr1[1:])).tolist(),
                    "half_widths": (0.5 * (ed_dr1[1:] - ed_dr1[:-1])).tolist(),
                    "sumw_den": [1000.0] * nb,
                    "sumw2_den": [2500.0] * nb,
                    "sumw_pass": [934.0] * nb,
                    "sumw2_pass": [2300.0] * nb,
                    "sumw2_fail": [200.0] * nb,
                    "eff": [0.934] * nb,
                    "err": [0.004] * nb,
                }
            },
            "mc_data_ratios_reweighted": {
                "dR1_wide": _make_ratio_curve(ed_dr1, 1.004),
            },
        }

    for pk in pairs:
        summary["paired_global_deltas"][pk] = {}
        summary["transition_matrices"][pk] = {}
        summary["paired_curves"][pk] = {}
        for ds in datasets:
            drec = {"n_den": 100000, "n_both_pass": 92900, "n_gain": 350, "n_loss": 150, "n_both_fail": 6600, "delta": 0.002, "delta_err": 0.00022}
            summary["paired_global_deltas"][pk][ds] = {
                "global": {"eff_pawan_c2R": drec, "eff_lt2": drec, "eff_ge2": drec, "eff_matched": drec},
                "small_dR1_lt5": {"eff_pawan_c2R": drec, "eff_lt2": drec, "eff_ge2": drec},
                "plateau_dR1_30_70": {"eff_pawan_c2R": drec, "eff_lt2": drec, "eff_ge2": drec},
                "large_abs_dX0_70_100": {"eff_pawan_c2R": drec, "eff_lt2": drec, "eff_ge2": drec},
            }
        for ds in ["All3Data", "10004X"]:
            mat = [[80, 15, 5, 0, 0], [10, 5600, 350, 30, 10], [5, 150, 93200, 120, 25], [0, 10, 90, 240, 10], [0, 5, 10, 10, 25]]
            summary["transition_matrices"][pk][ds] = {
                "labels": ["0", "1", "2", "3", "4+"],
                "global": mat,
                "small_dR1_lt5": mat,
                "plateau_dR1_30_70": mat,
                "large_abs_dX0_70_100": mat,
            }
            summary["paired_curves"][pk][ds] = {
                ck: {"eff_pawan_c2R": _make_paired_curve(ed, 0.002), "eff_lt2": _make_paired_curve(ed, 0.002)}
                for ck, ed in curve_specs.items()
            }

    # Signal 110038 JSON
    sig_curves = {}
    for ck, ed in curve_specs.items():
        cen = 0.5 * (ed[:-1] + ed[1:])
        if "dY" in ck and "signed" in ck:
            w = np.where(cen >= 0, np.exp(-cen / 6.0), 0.0)
        else:
            w = np.exp(-np.abs(cen) / 3.0)
        w = w / max(w.sum(), 1e-12)
        sig_curves[ck] = {
            "edges": ed.tolist(),
            "centers": cen.tolist(),
            "half_widths": (0.5 * (ed[1:] - ed[:-1])).tolist(),
            "counts": np.round(w * 10000).astype(int).tolist(),
            "norm": w.tolist(),
            "norm_err": (w * 0.05).tolist(),
            "total_in_range": 10000,
            "underflow": 0,
            "overflow": 0,
        }
    signal_json = {
        "schema_version": "1.0",
        "sample_label": "A' MC (m = 100 MeV, eps = 1e-5)",
        "dsid": 110038,
        "sources": {
            "master_dev": {
                "input_dir": "/eos/experiment/faser/data0/sim/mc24/foresee/110038/phy/master_dev",
                "n_files": 10,
                "n_events_total": 20000,
                "n_events_reco_2trk_os_r95": 10000,
                "two_track_eff_total": 0.50,
                "two_track_eff_err": 0.0035,
                "curves": sig_curves,
            },
            "latest_georot": (
                {
                    "input_dir": "/tmp/latest_110038",
                    "n_files": 10,
                    "n_events_total": 20000,
                    "n_events_reco_2trk_os_r95": 10100,
                    "two_track_eff_total": 0.505,
                    "two_track_eff_err": 0.0035,
                    "curves": sig_curves,
                }
                if include_latest
                else None
            ),
        },
    }
    return summary, signal_json


class TestOverlayFigures(unittest.TestCase):
    def test_signal_label_uses_100_mev_not_gev(self) -> None:
        self.assertIn("100", SIGNAL_LABEL)
        self.assertIn("MeV", SIGNAL_LABEL)
        self.assertNotIn("GeV", SIGNAL_LABEL)

    def test_style_helpers(self) -> None:
        yerr = safe_asymmetric_yerr(np.array([0.9, 0.8]), err_lo=np.array([-1e-16, 0.02]), err_hi=np.array([0.01, np.nan]))
        self.assertEqual(yerr.shape, (2, 2))
        self.assertTrue(np.all(yerr >= 0.0))

        self.assertEqual(choose_eff_ylim([(np.array([0.85, 0.94]), None)]), (0.80, 1.00))
        lo, hi = choose_eff_ylim([(np.array([0.74, 0.94]), None)])
        self.assertLessEqual(lo, 0.75)
        self.assertEqual(hi, 1.00)

        dlo, dhi, _, _ = choose_delta_ylim([(np.array([-0.2, 0.3]), None)])
        self.assertEqual((dlo, dhi), (-1.0, 1.0))

        # Deep asymmetric crater (-25.6 pp) must never be clipped by percentiles or 15 pp cap
        crater = np.full(100, 1.2)
        crater[16] = -25.6
        dlo2, dhi2, _, _ = choose_delta_ylim([(crater, None)], err_arrays=[np.full(100, 0.2)])
        self.assertLess(dlo2, -25.6)
        self.assertGreaterEqual(dhi2, 6.0)

        # Ratio y-limits keep [0.95, 1.05] normally, expand cleanly for deep craters
        rlo, rhi, _, _ = choose_ratio_ylim([(np.array([0.99, 1.01]), np.array([0.005, 0.005]), None)])
        self.assertEqual((rlo, rhi), (0.95, 1.05))
        rlo2, rhi2, _, _ = choose_ratio_ylim([(np.array([0.915, 1.01]), np.array([0.005, 0.005]), None)])
        self.assertLess(rlo2, 0.915)
        self.assertGreater(rhi2, 1.08)

    def test_slice_1d_range(self) -> None:
        edges = np.linspace(-100.0, 100.0, 201)
        vals = np.arange(200, dtype=float)
        se, [sv] = slice_1d_range(edges, [vals], -20.0, 20.0)
        self.assertEqual(len(sv), 40)
        self.assertAlmostEqual(se[0], -20.0)
        self.assertAlmostEqual(se[-1], 20.0)

    def test_baseline_mode_generates_13_figures(self) -> None:
        summary, sig = _make_synthetic_inputs(include_latest=False)
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            sp = tdp / "overlay_summary.json"
            sgp = tdp / "signal_110038.json"
            out_dir = tdp / "figs_base"
            sp.write_text(json.dumps(summary), encoding="utf-8")
            sgp.write_text(json.dumps(sig), encoding="utf-8")

            idx = generate_all_figures(sp, sgp, out_dir, mode="baseline", workers=2)
            self.assertEqual(idx["page_count"], 13)
            self.assertEqual(len(idx["skipped_in_mode"]), 7)
            self.assertTrue((out_dir / "stage08_overlay_2track_report.pdf").is_file())
            self.assertTrue((out_dir / "figure_index.json").is_file())
            for rec in idx["figures"]:
                self.assertTrue(Path(rec["png"]).is_file())
                self.assertTrue(Path(rec["pdf"]).is_file())

    def test_full_mode_generates_all_20_figures_and_pdf_pages(self) -> None:
        self.assertEqual(len(FIGURE_BUILDERS), 20)
        summary, sig = _make_synthetic_inputs(include_latest=True)
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            sp = tdp / "overlay_summary.json"
            sgp = tdp / "signal_110038.json"
            out_dir = tdp / "figs_full"
            sp.write_text(json.dumps(summary), encoding="utf-8")
            sgp.write_text(json.dumps(sig), encoding="utf-8")

            idx = generate_all_figures(sp, sgp, out_dir, mode="full", workers=4)
            self.assertEqual(idx["page_count"], 20)
            self.assertEqual(idx["skipped_in_mode"], [])
            pdf_path = out_dir / "stage08_overlay_2track_report.pdf"
            self.assertTrue(pdf_path.is_file())
            self.assertGreater(pdf_path.stat().st_size, 50_000)

            # Verify exact PDF page count via pdfinfo if available
            res = subprocess.run(["pdfinfo", str(pdf_path)], capture_output=True, text=True, check=False)
            if res.returncode == 0:
                self.assertIn("Pages:", res.stdout)
                for line in res.stdout.splitlines():
                    if line.startswith("Pages:"):
                        self.assertEqual(int(line.split(":")[1].strip()), 20)

    def test_tier2_figures_and_22_page_append(self) -> None:
        # Build minimal synthetic Tier-2 summary matching tier2_schema.md
        curve_specs = {
            "abs_dY0_zoom": np.linspace(0.0, 20.0, 41),
            "dR1_zoom": np.linspace(0.0, 20.0, 41),
            "dR1_wide": np.linspace(0.0, 100.0, 201),
            "abs_dX0_wide": np.linspace(0.0, 100.0, 101),
            "dY1_zoom_signed": np.linspace(-25.0, 25.0, 51),
        }
        regimes = [
            "global", "crater_abs_dY0_12p5_18p5", "small_abs_dY0_lt2",
            "large_abs_dY0_ge19p5", "small_dR1_lt5", "plateau_dR1_30_70", "large_abs_dX0_70_100",
        ]
        pairs = [pk for _, pk in DELTA_PAIRS] + ["t2a_covphi_minus_latest_218813c2"]
        t2 = {
            "schema_version": "1.0",
            "arms": ARM_ORDER,
            "pairs": pairs,
            "datasets": ["008749", "008749_diag", "008749_offdiag", "10004X", "10004X_diag", "10004X_offdiag"],
            "regimes": regimes,
            "provenance": {"frame": "FASER", "cond_tag": "OFLCOND-FASER-05", "matched_quadruples": {"008749": 18, "10004X": 18, "total": 36}},
            "global_summary": {},
            "paired_deltas": {},
            "stage_of_loss": {},
            "curves": {},
            "paired_curves": {},
            "mc_data_ratios": {},
            "momentum_residual": {},
        }
        for arm in ARM_ORDER:
            t2["global_summary"][arm] = {}
            t2["stage_of_loss"][arm] = {}
            t2["curves"][arm] = {}
            t2["mc_data_ratios"][arm] = {ck: {"eff_pawan_c2R": _make_ratio_curve(ed, 1.005)} for ck, ed in curve_specs.items()}
            for ds in ["008749", "10004X"]:
                eff_rec = {"pass": 93000, "den": 100000, "eff": 0.93, "err_lo": 0.0008, "err_hi": 0.0008}
                t2["global_summary"][arm][ds] = {
                    "n_files": 18, "n_events": 260000 if ds == "008749" else 690000,
                    "lt_counts": {"0": 200, "1": 5000, "2": 94000, "3": 700, "4": 90, "5plus": 10},
                    "lt2_fail_c2R": 300,
                    "regimes": {r: {"eff_pawan_c2R": eff_rec, "eff_lt2": eff_rec, "eff_ge2": eff_rec, "eff_matched": eff_rec} for r in regimes},
                }
                sl = {
                    "n_den": 100000, "pass_c2R": 93500, "lt2_fail_c2R": 200, "lt_ge3": 300,
                    "lt0": 1800 if arm == "latest_218813c2" else 200,
                    "lt1_total": 4200, "lt1_miss_seg_st1": 400, "lt1_miss_seg_st2": 1000, "lt1_miss_seg_st3": 1200,
                    "lt1_miss_seg_multi": 200, "lt1_killed_by_gb": 50,
                    "lt1_all3gb_lost_in_ckf": 1300 if arm == "latest_218813c2" else 500,
                    "lt1_unmatched_kept": 50,
                }
                t2["stage_of_loss"][arm][ds] = {r: sl for r in regimes}
                t2["curves"][arm][ds] = {}
                for ck, ed in curve_specs.items():
                    c = _make_eff_curve(ed, 0.935 if arm == "t2a_covphi" else 0.925)
                    if arm == "latest_218813c2" and ck == "abs_dY0_zoom":
                        # Inject realistic -25 pp crater at 15-17 mm
                        for idx in range(30, 35):
                            c["eff"][idx] = 0.67
                    t2["curves"][arm][ds][ck] = {"eff_pawan_c2R": c, "eff_ge2": c}

        for pk in pairs:
            t2["paired_deltas"][pk] = {}
            t2["paired_curves"][pk] = {}
            for ds in ["008749", "10004X"]:
                t2["paired_deltas"][pk][ds] = {}
                for r in regimes:
                    dval = -0.22 if ("latest_218813c2" in pk and "crater" in r) else 0.015
                    drec = {"n_den": 50000, "n_both_pass": 45000, "n_gain": 900, "n_loss": 150, "n_both_fail": 3950, "delta": dval, "delta_err": 0.001}
                    t2["paired_deltas"][pk][ds][r] = {"eff_pawan_c2R": drec, "eff_lt2": drec, "eff_ge2": drec, "eff_matched": drec}
                t2["paired_curves"][pk][ds] = {}
                for ck, ed in curve_specs.items():
                    pc = _make_paired_curve(ed, 0.012 if "t2a_covphi" in pk else -0.010)
                    if "latest_218813c2" in pk and ck == "abs_dY0_zoom":
                        for idx in range(30, 35):
                            pc["delta"][idx] = -0.256
                    t2["paired_curves"][pk][ds][ck] = {"eff_pawan_c2R": pc, "eff_ge2": pc}

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            t2p = tdp / "tier2_attribution_summary.json"
            out_dir = tdp / "figs_t2"
            t2p.write_text(json.dumps(t2), encoding="utf-8")
            res = generate_tier2_figures(t2p, out_dir, append_report=False)
            self.assertEqual(len(res["figures"]), 2)
            self.assertTrue((out_dir / "png" / "fig21_tier2_covphi_crater_attribution.png").is_file())
            self.assertTrue((out_dir / "pdf" / "fig21_tier2_covphi_crater_attribution.pdf").is_file())
            self.assertTrue((out_dir / "png" / "fig22_tier2_covphi_regimes_and_loss.png").is_file())
            self.assertTrue((out_dir / "pdf" / "fig22_tier2_covphi_regimes_and_loss.pdf").is_file())
            self.assertTrue((out_dir / "stage08_tier2_attribution_figures.pdf").is_file())


if __name__ == "__main__":
    unittest.main()
