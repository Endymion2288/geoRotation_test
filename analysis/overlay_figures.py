#!/usr/bin/env python3
"""Stage 08 Two-Track Overlay Tracking Performance Figure Suite (20 figures).

Generates publication-grade English figures mirroring Pawan's FASER Tracking Meeting report
(Slides 7-16), adding full paired Baseline-vs-Latest GeoRotation comparisons and deep
physics/tracking diagnostics (multiplicity, stage-of-loss funnel, large-|dX| / 2D loss maps,
momentum reweighting closure, and May30-vs-Jul8 baseline validation).

Outputs:
- <out_dir>/png/figNN_*.png   (300 DPI)
- <out_dir>/pdf/figNN_*.pdf   (vector PDF)
- <out_dir>/stage08_overlay_2track_report.pdf  (unified multi-page vector PDF)
- <out_dir>/figure_index.json (manifest of figures, paths, references, and printed stats)
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Circle
from matplotlib.ticker import AutoMinorLocator, LogLocator, MultipleLocator
import numpy as np

# Ensure repository root is importable when executed as a standalone script
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from analysis.overlay_plot_style import (
    COLORS,
    DIAG_PALETTE,
    SIGNAL_LABEL,
    SIGNAL_LATEST_LABEL,
    add_figure_banners,
    apply_style,
    choose_delta_ylim,
    choose_eff_ylim,
    choose_ratio_ylim,
    combine_legend,
    configure_efficiency_axis,
    configure_ratio_axis,
    draw_signal_overlay,
    format_count_compact,
    format_faser_header,
    make_single_ratio_figure,
    make_two_column_ratio_figure,
    plot_binned_errorbar,
    safe_asymmetric_yerr,
    save_figure_outputs,
)


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------
VERSION_LABELS: Dict[str, str] = {
    "latest_georot": "Latest GeoRotation (fix_v1, 218813c2)",
    "baseline_jul8": "Baseline-Jul8 (OverlaidNTPs, 811b81b8)",
    "baseline_may30": "Baseline-May30 (OverlaidNTPs_bug, 54a9026b)",
}

VERSION_SHORT: Dict[str, str] = {
    "latest_georot": "Latest GeoRot",
    "baseline_jul8": "Baseline-Jul8",
    "baseline_may30": "Baseline-May30",
}

METRIC_FOOTNOTE = (
    "Metric: eff_pawan_c2R (longTracks == 2, chi2/nDoF <= 25, R0,R1 < 95 mm) | "
    "Frame: FASER [mm, MeV] | Cond: OFLCOND-FASER-05 (TI12MC04 / FASERNU-04)"
)


def to_float_array(seq: Optional[Sequence[Any]]) -> np.ndarray:
    """Convert JSON array (possibly containing None) to float64 numpy array."""
    if seq is None:
        return np.array([], dtype=float)
    return np.array([np.nan if v is None else float(v) for v in seq], dtype=float)


def to_int_array(seq: Optional[Sequence[Any]]) -> np.ndarray:
    if seq is None:
        return np.array([], dtype=int)
    return np.array([0 if v is None else int(v) for v in seq], dtype=int)


def slice_1d_range(
    edges: np.ndarray,
    arrays: Sequence[np.ndarray],
    xmin: float,
    xmax: float,
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Slice bin edges and parallel bin-value arrays to [xmin, xmax]."""
    edges = np.asarray(edges, dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    sel = (centers >= xmin) & (centers <= xmax)
    idx = np.where(sel)[0]
    if idx.size == 0:
        return edges, [np.asarray(a, dtype=float) for a in arrays]
    i0, i1 = int(idx[0]), int(idx[-1]) + 1
    sub_edges = edges[i0 : i1 + 1]
    sub_arrays = [np.asarray(a, dtype=float)[i0:i1] for a in arrays]
    return sub_edges, sub_arrays


def get_signal_curve(
    signal_data: Optional[Dict[str, Any]],
    curve_key: str,
    source: str = "master_dev",
) -> Optional[Dict[str, Any]]:
    if not signal_data:
        return None
    src = (signal_data.get("sources") or {}).get(source)
    if not src:
        return None
    return (src.get("curves") or {}).get(curve_key)


def get_scope_line(summary: Dict[str, Any], primary_ver: str, extra: str = "") -> str:
    gs = summary.get("global_summary", {}).get(primary_ver, {})
    n_data = int(gs.get("All3Data", {}).get("n_events", 0))
    n_mc = int(gs.get("10004X", {}).get("n_events", 0))
    prov = summary.get("provenance", {})
    pairs = prov.get("matched_bin_pairs", {})
    n_pairs = int(pairs.get("total", 322))
    ver_lbl = VERSION_LABELS.get(primary_ver, primary_ver)
    parts = [
        ver_lbl,
        f"All3Data: {format_count_compact(n_data)} evts",
        f"10004X MC: {format_count_compact(n_mc)} evts",
        f"{n_pairs} matched pairs",
    ]
    if extra:
        parts.append(extra)
    return "  |  ".join(parts)


def _bin_eff_stats(curve: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    edges = to_float_array(curve["edges"])
    eff = to_float_array(curve["eff"])
    den = to_int_array(curve["den"])
    yerr = safe_asymmetric_yerr(eff, err_lo=to_float_array(curve["err_lo"]), err_hi=to_float_array(curve["err_hi"]))
    return edges, eff, yerr, den > 0


def _ratio_stats(rcurve: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    edges = to_float_array(rcurve["edges"])
    ratio = to_float_array(rcurve["ratio"])
    yerr = safe_asymmetric_yerr(ratio, err_lo=to_float_array(rcurve["err_lo"]), err_hi=to_float_array(rcurve["err_hi"]))
    return edges, ratio, yerr, np.isfinite(ratio) & (ratio > 0)


def _paired_delta_pp(pcurve: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    edges = to_float_array(pcurve["edges"])
    delta = to_float_array(pcurve["delta"]) * 100.0
    derr = to_float_array(pcurve["delta_err"]) * 100.0
    den = to_int_array(pcurve["n_den"])
    return edges, delta, derr, (den > 0) & np.isfinite(delta)


def _draw_standard_suite_a_panel(
    ax_top: Axes,
    ax_bot: Axes,
    summary: Dict[str, Any],
    signal_data: Optional[Dict[str, Any]],
    version: str,
    curve_key: str,
    xlabel: str,
    xlim: Tuple[float, float],
    legend_loc: str = "lower right",
    panel_tag: str = "",
) -> Dict[str, Any]:
    """Render one (top efficiency + twin A' signal, bottom MC/Data ratio) column for Suite A."""
    mc_c = summary["curves"][version]["10004X"][curve_key]["eff_pawan_c2R"]
    da_c = summary["curves"][version]["All3Data"][curve_key]["eff_pawan_c2R"]
    ra_c = summary["mc_data_ratios"][version][curve_key]["eff_pawan_c2R"]

    edges_m, eff_m, yerr_m, mask_m = _bin_eff_stats(mc_c)
    edges_d, eff_d, yerr_d, mask_d = _bin_eff_stats(da_c)
    edges_r, rat, yerr_r, mask_r = _ratio_stats(ra_c)

    n_bins = len(edges_m) - 1
    ms = 4.0 if n_bins <= 60 else (3.1 if n_bins <= 120 else 2.5)
    ew = 0.95 if n_bins <= 60 else 0.78

    ylim = choose_eff_ylim([(eff_m, mask_m), (eff_d, mask_d)], default=(0.80, 1.00), min_floor=0.60)
    configure_efficiency_axis(ax_top, xlim=xlim, ylim=ylim, has_right_twin=(signal_data is not None))
    format_faser_header(ax_top, subtitle="Overlay Efficiency", right_text=panel_tag)

    extra_handles = []
    sig_c = get_signal_curve(signal_data, curve_key, source="master_dev")
    if sig_c is not None:
        _, h_sig = draw_signal_overlay(
            ax_top,
            to_float_array(sig_c["edges"]),
            to_float_array(sig_c["norm"]),
            label=SIGNAL_LABEL,
            color=COLORS["signal"],
        )
        extra_handles.append(h_sig)

    plot_binned_errorbar(
        ax_top, edges_m, eff_m, yerr_m, mask=mask_m,
        label="Muon MC Overlay (10004X)", color=COLORS["mc"], marker="o",
        markersize=ms, elinewidth=ew, zorder=6,
    )
    plot_binned_errorbar(
        ax_top, edges_d, eff_d, yerr_d, mask=mask_d,
        label="Data Overlay (All3Data)", color=COLORS["data"], marker="o",
        markersize=ms, elinewidth=ew, zorder=5,
    )
    combine_legend(ax_top, extra_handles=extra_handles, loc=legend_loc, fontsize=8.1)

    r_lo, r_hi, r_maj, r_min = choose_ratio_ylim([(rat, yerr_r, mask_r)], default=(0.95, 1.05))
    configure_ratio_axis(
        ax_bot, xlim=xlim, xlabel=xlabel, ylabel="MC / Data",
        ylim=(r_lo, r_hi), ref_value=1.0, major_ytick=r_maj, minor_ytick=r_min,
    )
    plot_binned_errorbar(
        ax_bot, edges_r, rat, yerr_r, mask=mask_r,
        label="", color=COLORS["ratio"], marker="o",
        markersize=max(2.2, ms - 0.4), elinewidth=ew, zorder=5,
    )

    return {
        "curve_key": curve_key,
        "mc_in_range": int(np.sum(to_int_array(mc_c["den"]))),
        "data_in_range": int(np.sum(to_int_array(da_c["den"]))),
        "mc_overflow": int(mc_c.get("overflow_den", 0)),
        "data_overflow": int(da_c.get("overflow_den", 0)),
        "ylim": list(ylim),
    }


# ---------------------------------------------------------------------------
# Suite A: Primary Report Figures (fig01 - fig08)
# ---------------------------------------------------------------------------
def _make_suite_a_pair(
    summary: Dict[str, Any],
    signal_data: Optional[Dict[str, Any]],
    out_dir: Path,
    version: str,
    fig_id: str,
    stem: str,
    title: str,
    pawan_ref: str,
    left_spec: Tuple[str, str, Tuple[float, float], str, str],
    right_spec: Tuple[str, str, Tuple[float, float], str, str],
) -> Dict[str, Any]:
    apply_style()
    fig, cols = make_two_column_ratio_figure()
    s_left = _draw_standard_suite_a_panel(
        cols[0][0], cols[0][1], summary, signal_data, version,
        curve_key=left_spec[0], xlabel=left_spec[1], xlim=left_spec[2],
        legend_loc=left_spec[3], panel_tag=left_spec[4],
    )
    s_right = _draw_standard_suite_a_panel(
        cols[1][0], cols[1][1], summary, signal_data, version,
        curve_key=right_spec[0], xlabel=right_spec[1], xlim=right_spec[2],
        legend_loc=right_spec[3], panel_tag=right_spec[4],
    )
    mode_tag = "" if version == "latest_georot" else " [Baseline-Jul8 Preview]"
    add_figure_banners(
        fig,
        title=f"{fig_id.upper()}: {title}{mode_tag}",
        scope_line=get_scope_line(summary, version),
        footer_line=f"{METRIC_FOOTNOTE} | Ref: Pawan {pawan_ref} (A' unit corrected to 100 MeV)",
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, stem)
    gs = summary["global_summary"][version]
    return {
        "fig_id": fig_id,
        "suite": "A",
        "title": f"{title}{mode_tag}",
        "pawan_ref": pawan_ref,
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": {
            "version": version,
            "data_global_eff": gs["All3Data"]["eff_pawan_c2R"]["eff"],
            "mc_global_eff": gs["10004X"]["eff_pawan_c2R"]["eff"],
            "left": s_left,
            "right": s_right,
        },
    }


def make_fig01(summary, signal_data, out_dir, version, **_):
    return _make_suite_a_pair(
        summary, signal_data, out_dir, version,
        fig_id="fig01",
        stem="fig01_dR1_efficiency",
        title="Two-Track Overlay Efficiency vs Station 3 Radial Separation (Delta R1)",
        pawan_ref="Slide 7",
        left_spec=("dR1_wide", r"Track Separation ($\Delta R_1$) [mm]", (0.0, 100.0), "lower right", "Wide [0, 100] mm (0.5 mm bins)"),
        right_spec=("dR1_zoom", r"Track Separation ($\Delta R_1$) [mm]", (0.0, 20.0), "lower right", "Zoom [0, 20] mm (0.5 mm bins)"),
    )


def make_fig02(summary, signal_data, out_dir, version, **_):
    return _make_suite_a_pair(
        summary, signal_data, out_dir, version,
        fig_id="fig02",
        stem="fig02_dR0_efficiency",
        title="Two-Track Overlay Efficiency vs Station 1 Radial Separation (Delta R0)",
        pawan_ref="Slide 8",
        left_spec=("dR0_wide", r"Track Separation ($\Delta R_0$) [mm]", (0.0, 100.0), "lower right", "Wide [0, 100] mm (0.5 mm bins)"),
        right_spec=("dR0_zoom", r"Track Separation ($\Delta R_0$) [mm]", (0.0, 20.0), "lower right", "Zoom [0, 20] mm (0.5 mm bins)"),
    )


def make_fig03(summary, signal_data, out_dir, version, **_):
    return _make_suite_a_pair(
        summary, signal_data, out_dir, version,
        fig_id="fig03",
        stem="fig03_abs_dY0_efficiency",
        title="Two-Track Overlay Efficiency vs Station 1 Vertical Separation (|Delta Y0|)",
        pawan_ref="Slide 9",
        left_spec=("abs_dY0_wide", r"Track Separation ($|\Delta Y_0|$) [mm]", (0.0, 100.0), "lower right", "Wide [0, 100] mm (1.0 mm bins)"),
        right_spec=("abs_dY0_zoom", r"Track Separation ($|\Delta Y_0|$) [mm]", (0.0, 20.0), "lower left", "Zoom [0, 20] mm (0.5 mm bins)"),
    )


def make_fig04(summary, signal_data, out_dir, version, **_):
    return _make_suite_a_pair(
        summary, signal_data, out_dir, version,
        fig_id="fig04",
        stem="fig04_abs_dY1_efficiency",
        title="Two-Track Overlay Efficiency vs Station 3 Vertical Separation (|Delta Y1|)",
        pawan_ref="Slide 10",
        left_spec=("abs_dY1_wide", r"Track Separation ($|\Delta Y_1|$) [mm]", (0.0, 100.0), "lower right", "Wide [0, 100] mm (1.0 mm bins)"),
        right_spec=("abs_dY1_zoom", r"Track Separation ($|\Delta Y_1|$) [mm]", (0.0, 20.0), "lower right", "Zoom [0, 20] mm (0.5 mm bins)"),
    )


def make_fig05(summary, signal_data, out_dir, version, **_):
    return _make_suite_a_pair(
        summary, signal_data, out_dir, version,
        fig_id="fig05",
        stem="fig05_dX0_efficiency",
        title="Two-Track Overlay Efficiency vs Station 1 Horizontal Separation (|Delta X0|, Delta X0)",
        pawan_ref="Slide 11",
        left_spec=("abs_dX0_wide", r"Track Separation ($|\Delta X_0|$) [mm]", (0.0, 100.0), "lower left", "Unsigned [0, 100] mm (1.0 mm bins)"),
        right_spec=("dX0_signed", r"Track Separation ($\Delta X_0 = x_0^+ - x_0^-$) [mm]", (-100.0, 100.0), "lower right", "Signed [-100, 100] mm (1.0 mm bins)"),
    )


def make_fig06(summary, signal_data, out_dir, version, **_):
    return _make_suite_a_pair(
        summary, signal_data, out_dir, version,
        fig_id="fig06",
        stem="fig06_dX1_efficiency",
        title="Two-Track Overlay Efficiency vs Station 3 Horizontal Separation (|Delta X1|, Delta X1)",
        pawan_ref="Slide 12",
        left_spec=("abs_dX1_wide", r"Track Separation ($|\Delta X_1|$) [mm]", (0.0, 100.0), "lower left", "Unsigned [0, 100] mm (1.0 mm bins)"),
        right_spec=("dX1_signed", r"Track Separation ($\Delta X_1 = x_1^+ - x_1^-$) [mm]", (-100.0, 100.0), "lower right", "Signed [-100, 100] mm (1.0 mm bins)"),
    )


def make_fig07(summary, signal_data, out_dir, version, **_):
    return _make_suite_a_pair(
        summary, signal_data, out_dir, version,
        fig_id="fig07",
        stem="fig07_signed_dY_efficiency",
        title="Two-Track Overlay Efficiency vs Signed Vertical Separation (Delta Y1, Delta Y0)",
        pawan_ref="Slide 13",
        left_spec=("dY1_signed", r"Track Separation ($\Delta Y_1 = y_1^+ - y_1^-$) [mm]", (-100.0, 100.0), "lower right", "Station 3 [-100, 100] mm (1.0 mm bins)"),
        right_spec=("dY0_signed", r"Track Separation ($\Delta Y_0 = y_0^+ - y_0^-$) [mm]", (-100.0, 100.0), "lower right", "Station 1 [-100, 100] mm (1.0 mm bins)"),
    )


def make_fig08(summary, signal_data, out_dir, version, **_):
    """Mirrors Slide 16: Reconstructed Track Pz0 [MeV] normalized distribution + MC/Data ratio."""
    apply_style()
    fig, (ax_top, ax_bot) = make_single_ratio_figure(figsize=(9.6, 6.4), left=0.10, right=0.95)
    mom = summary["momentum"][version]
    mc_p = mom["10004X"]["track_pz0"]
    da_p = mom["All3Data"]["track_pz0"]
    rat_p = summary["momentum"]["mc_data_pz0_ratio"][version]

    edges = to_float_array(mc_p["edges"])
    mc_n = to_float_array(mc_p["norm"])
    da_n = to_float_array(da_p["norm"])
    r_val = to_float_array(rat_p["ratio"])
    r_err = to_float_array(rat_p["err"])

    # Shaded step histograms matching Slide 16
    ax_top.fill_between(edges[:-1], 0.0, mc_n, step="post", color=COLORS["mc"], alpha=0.28, zorder=2)
    ax_top.step(edges, np.concatenate([mc_n, mc_n[-1:]]), where="post", color=COLORS["mc"], linewidth=1.25,
                label=f"Muon MC Overlay ({format_count_compact(mc_p['total_tracks'])} tracks)", zorder=4)
    ax_top.fill_between(edges[:-1], 0.0, da_n, step="post", color=COLORS["data"], alpha=0.28, zorder=3)
    ax_top.step(edges, np.concatenate([da_n, da_n[-1:]]), where="post", color=COLORS["data"], linewidth=1.25,
                label=f"Data Overlay ({format_count_compact(da_p['total_tracks'])} tracks)", zorder=5)

    ymax = max(0.08, float(np.nanmax(np.concatenate([mc_n, da_n]))) * 1.16)
    ax_top.set_xscale("log")
    ax_top.set_xlim(1.0, 2.0e7)
    ax_top.set_ylim(0.0, ymax)
    ax_top.set_ylabel("Normalized Tracks / Bin (Unit Area)")
    ax_top.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax_top.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax_top.grid(axis="y", color=COLORS["grid"], linestyle=":", linewidth=0.6, alpha=0.65)
    format_faser_header(ax_top, subtitle="Reconstructed Track Momentum Spectrum", right_text="100 log10 bins [1, 2e7] MeV")
    combine_legend(ax_top, loc="upper left", fontsize=8.8)

    # Compact flow & turn-on annotation
    ax_top.text(
        0.025, 0.72,
        f"SelectRDOAlg cut: 50 GeV <= og_p <= 5 TeV\n"
        f"MC under/overflow: {mc_p.get('underflow', 0):,} / {mc_p.get('overflow', 0):,}\n"
        f"Data under/overflow: {da_p.get('underflow', 0):,} / {da_p.get('overflow', 0):,}",
        transform=ax_top.transAxes, fontsize=8.0, color="#333333",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#dddddd", alpha=0.90),
    )

    # Bottom MC / Data ratio
    valid_r = np.isfinite(r_val) & (r_val > 0.0) & (to_int_array(da_p["counts"]) >= 5)
    rmax = float(np.nanpercentile(r_val[valid_r], 99)) if np.any(valid_r) else 9.0
    r_ylim_hi = max(9.0, float(np.ceil(min(rmax * 1.15, 16.0))))
    ax_bot.set_xscale("log")
    ax_bot.set_xlim(1.0, 2.0e7)
    ax_bot.set_ylim(-0.5, r_ylim_hi)
    ax_bot.set_xlabel(r"Reconstructed Track $p_{z,0}$ [MeV]")
    ax_bot.set_ylabel("MC / Data")
    ax_bot.yaxis.set_major_locator(MultipleLocator(2.0))
    ax_bot.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax_bot.axhline(1.0, color=COLORS["ref_line"], linestyle="--", linewidth=0.9, zorder=1)
    ax_bot.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax_bot.grid(axis="y", color=COLORS["grid"], linestyle=":", linewidth=0.55, alpha=0.55)
    centers_log = np.sqrt(edges[:-1] * edges[1:])
    ax_bot.errorbar(
        centers_log[valid_r], r_val[valid_r], yerr=np.where(np.isfinite(r_err[valid_r]), r_err[valid_r], 0.0),
        fmt="o", color=COLORS["ratio"], markersize=3.0, elinewidth=0.8, capsize=0.0, zorder=5,
    )

    mode_tag = "" if version == "latest_georot" else " [Baseline-Jul8 Preview]"
    add_figure_banners(
        fig,
        title=f"FIG08: Reconstructed Track Momentum Spectrum (Track Pz0){mode_tag}",
        scope_line=get_scope_line(summary, version),
        footer_line="All reco long tracks | Frame: FASER [MeV] | Cond: OFLCOND-FASER-05 (TI12MC04) | Ref: Pawan Slide 16",
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, "fig08_track_pz0_spectrum")
    return {
        "fig_id": "fig08",
        "suite": "A",
        "title": f"Reconstructed Track Pz0 Spectrum{mode_tag}",
        "pawan_ref": "Slide 16",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": {
            "version": version,
            "mc_total_tracks": int(mc_p["total_tracks"]),
            "data_total_tracks": int(da_p["total_tracks"]),
            "peak_mc_data_ratio": float(np.nanmax(r_val[valid_r])) if np.any(valid_r) else None,
        },
    }


# ---------------------------------------------------------------------------
# Suite B: Direct Paired Comparison — Baseline-Jul8 vs Latest GeoRotation (fig09 - fig15)
# ---------------------------------------------------------------------------
def _draw_comparison_column(
    ax_top: Axes,
    ax_delta: Axes,
    summary: Dict[str, Any],
    signal_data: Optional[Dict[str, Any]],
    curve_key: str,
    xlabel: str,
    xlim: Tuple[float, float],
    panel_tag: str = "",
    legend_loc: str = "lower right",
    ax_ratio: Optional[Axes] = None,
    include_may30: bool = False,
    slice_range: Optional[Tuple[float, float]] = None,
    include_latest_signal: bool = False,
) -> Dict[str, Any]:
    """Draw a paired Baseline-Jul8 vs Latest GeoRotation comparison column."""
    ver_new = "latest_georot"
    ver_old = "baseline_jul8"
    pair_key = "latest_georot_minus_baseline_jul8"

    def _extract(ver: str, ds: str):
        c = summary["curves"][ver][ds][curve_key]["eff_pawan_c2R"]
        e, v, ye, m = _bin_eff_stats(c)
        if slice_range is not None:
            e, [v, ye0, ye1, mf] = slice_1d_range(e, [v, ye[0], ye[1], m.astype(float)], *slice_range)
            ye = np.vstack([ye0, ye1])
            m = mf > 0.5
        return e, v, ye, m

    e_nm, v_nm, ye_nm, m_nm = _extract(ver_new, "10004X")
    e_nd, v_nd, ye_nd, m_nd = _extract(ver_new, "All3Data")
    e_om, v_om, ye_om, m_om = _extract(ver_old, "10004X")
    e_od, v_od, ye_od, m_od = _extract(ver_old, "All3Data")

    n_bins = len(e_nm) - 1
    ms = 3.8 if n_bins <= 60 else (2.9 if n_bins <= 120 else 2.3)
    ew = 0.85 if n_bins <= 60 else 0.72

    series_for_ylim = [(v_nm, m_nm), (v_nd, m_nd), (v_om, m_om), (v_od, m_od)]
    ylim = choose_eff_ylim(series_for_ylim, default=(0.80, 1.00), min_floor=0.60)
    configure_efficiency_axis(ax_top, xlim=xlim, ylim=ylim, has_right_twin=(signal_data is not None))
    format_faser_header(ax_top, subtitle="Jul8 vs Latest GeoRot", right_text=panel_tag, fontsize=9.6, right_fontsize=8.0)

    extra_handles = []
    sig_c = get_signal_curve(signal_data, curve_key, source="master_dev")
    ax_twin = None
    if sig_c is not None:
        se = to_float_array(sig_c["edges"])
        sn = to_float_array(sig_c["norm"])
        if slice_range is not None:
            se, [sn] = slice_1d_range(se, [sn], *slice_range)
        sig_lbl = SIGNAL_LABEL if not include_latest_signal else r"A' MC ($m = 100\ \mathrm{MeV},\ \varepsilon = 10^{-5}$, truth fid.)"
        ax_twin, h_sig = draw_signal_overlay(ax_top, se, sn, label=sig_lbl, color=COLORS["signal"])
        extra_handles.append(h_sig)
        if include_latest_signal:
            # Plot Latest GeoRotation reco_curves only if distinct from truth-fiducial curve
            src_l = (signal_data.get("sources") or {}).get("latest_georot") or {}
            sig_new = (src_l.get("reco_curves") or {}).get(curve_key)
            if sig_new is not None:
                se2 = to_float_array(sig_new["edges"])
                sn2 = to_float_array(sig_new["norm"])
                if slice_range is not None:
                    se2, [sn2] = slice_1d_range(se2, [sn2], *slice_range)
                if sn2.shape != sn.shape or not np.allclose(sn2, sn, atol=1e-6, equal_nan=True):
                    _, h_sig2 = draw_signal_overlay(
                        ax_top, se2, sn2, label=r"A' MC Latest GeoRot (reco 2-trk)",
                        color=COLORS["signal_latest"], linestyle="--", existing_twin=ax_twin,
                    )
                    extra_handles.append(h_sig2)

    if include_may30 and "baseline_may30" in summary.get("curves", {}):
        e_m30, v_m30, _, m_m30 = _extract("baseline_may30", "All3Data")
        centers_m30 = 0.5 * (e_m30[:-1] + e_m30[1:])
        ax_top.plot(
            centers_m30[m_m30], v_m30[m_m30],
            linestyle=":", color=COLORS["slate"], linewidth=1.0, alpha=0.8,
            label="Baseline-May30 Data (ref)", zorder=3,
        )

    # Baseline-Jul8: open triangles; Latest GeoRotation: solid circles
    plot_binned_errorbar(
        ax_top, e_om, v_om, ye_om, mask=m_om,
        label="Baseline-Jul8 MC (10004X)", color=COLORS["mc_base"], marker="^",
        markerfacecolor="white", markersize=ms + 0.3, elinewidth=ew, alpha=0.85, zorder=4,
    )
    plot_binned_errorbar(
        ax_top, e_od, v_od, ye_od, mask=m_od,
        label="Baseline-Jul8 Data (All3Data)", color=COLORS["data_base"], marker="^",
        markerfacecolor="white", markersize=ms + 0.3, elinewidth=ew, alpha=0.85, zorder=4,
    )
    plot_binned_errorbar(
        ax_top, e_nm, v_nm, ye_nm, mask=m_nm,
        label="Latest GeoRot MC (10004X)", color=COLORS["mc"], marker="o",
        markersize=ms, elinewidth=ew, zorder=6,
    )
    plot_binned_errorbar(
        ax_top, e_nd, v_nd, ye_nd, mask=m_nd,
        label="Latest GeoRot Data (All3Data)", color=COLORS["data"], marker="o",
        markersize=ms, elinewidth=ew, zorder=5,
    )
    combine_legend(ax_top, extra_handles=extra_handles, loc=legend_loc, fontsize=7.6)

    # Paired McNemar delta panel (percentage points)
    pm = summary["paired_curves"][pair_key]["10004X"][curve_key]["eff_pawan_c2R"]
    pd = summary["paired_curves"][pair_key]["All3Data"][curve_key]["eff_pawan_c2R"]
    e_pm, dm, dem, mm = _paired_delta_pp(pm)
    e_pd, dd, ded, md = _paired_delta_pp(pd)
    if slice_range is not None:
        e_pm, [dm, dem, mmf] = slice_1d_range(e_pm, [dm, dem, mm.astype(float)], *slice_range)
        mm = mmf > 0.5
        e_pd, [dd, ded, mdf] = slice_1d_range(e_pd, [dd, ded, md.astype(float)], *slice_range)
        md = mdf > 0.5

    d_lo, d_hi, d_maj, d_min = choose_delta_ylim(
        [(dm, mm), (dd, md)], min_half_span=1.0, max_half_span=35.0, err_arrays=[dem, ded]
    )
    delta_xlabel = "" if ax_ratio is not None else xlabel
    # Pad delta ylim slightly (8%) so extreme tick labels never collide with adjacent stacked pads
    pad = 0.08 * (d_hi - d_lo)
    configure_ratio_axis(
        ax_delta, xlim=xlim, xlabel=delta_xlabel,
        ylabel=r"$\Delta\varepsilon$ (Latest$-$Jul8) [pp]",
        ylim=(d_lo - pad, d_hi + pad), ref_value=0.0, major_ytick=d_maj, minor_ytick=d_min,
    )
    plot_binned_errorbar(
        ax_delta, e_pm, dm, dem, mask=mm,
        label="MC (10004X)", color=COLORS["mc"], marker="o",
        markersize=max(2.2, ms - 0.5), elinewidth=ew, zorder=6,
    )
    plot_binned_errorbar(
        ax_delta, e_pd, dd, ded, mask=md,
        label="Data (All3Data)", color=COLORS["data"], marker="s",
        markersize=max(2.1, ms - 0.6), elinewidth=ew, zorder=5,
    )
    combine_legend(ax_delta, loc="upper right", ncol=2, fontsize=7.3)

    # Optional third panel: MC/Data ratio comparison
    if ax_ratio is not None:
        r_new = summary["mc_data_ratios"][ver_new][curve_key]["eff_pawan_c2R"]
        r_old = summary["mc_data_ratios"][ver_old][curve_key]["eff_pawan_c2R"]
        e_rn, rn, yern, mrn = _ratio_stats(r_new)
        e_ro, ro, yero, mro = _ratio_stats(r_old)
        if slice_range is not None:
            e_rn, [rn, yern0, yern1, mrnf] = slice_1d_range(e_rn, [rn, yern[0], yern[1], mrn.astype(float)], *slice_range)
            yern, mrn = np.vstack([yern0, yern1]), mrnf > 0.5
            e_ro, [ro, yero0, yero1, mrof] = slice_1d_range(e_ro, [ro, yero[0], yero[1], mro.astype(float)], *slice_range)
            yero, mro = np.vstack([yero0, yero1]), mrof > 0.5
        configure_ratio_axis(ax_ratio, xlim=xlim, xlabel=xlabel, ylabel="MC / Data", ylim=(0.95, 1.05), ref_value=1.0)
        plot_binned_errorbar(
            ax_ratio, e_ro, ro, yero, mask=mro,
            label="Baseline-Jul8", color="#666666", marker="^", markerfacecolor="white",
            markersize=max(2.3, ms - 0.3), elinewidth=ew, alpha=0.85, zorder=4,
        )
        plot_binned_errorbar(
            ax_ratio, e_rn, rn, yern, mask=mrn,
            label="Latest GeoRot", color=COLORS["ratio"], marker="o",
            markersize=max(2.2, ms - 0.4), elinewidth=ew, zorder=6,
        )
        combine_legend(ax_ratio, loc="lower right", ncol=2, fontsize=7.3)

    return {
        "curve_key": curve_key,
        "mc_delta_pp_range": [float(np.nanmin(dm[mm])), float(np.nanmax(dm[mm]))] if np.any(mm) else None,
        "data_delta_pp_range": [float(np.nanmin(dd[md])), float(np.nanmax(dd[md]))] if np.any(md) else None,
    }


def make_fig09(summary, signal_data, out_dir, version, **_):
    """Mirrors Slide 14: Primary Baseline-Jul8 vs Latest GeoRotation comparison on Delta R1 (3-row)."""
    apply_style()
    fig, cols = make_two_column_ratio_figure(
        figsize=(14.6, 7.5), height_ratios=(3.2, 1.18, 1.08),
        top=0.87, bottom=0.095, wspace=0.34, hspace=0.085,
    )
    s_l = _draw_comparison_column(
        cols[0][0], cols[0][1], summary, signal_data,
        curve_key="dR1_wide", xlabel=r"Track Separation ($\Delta R_1$) [mm]", xlim=(0.0, 100.0),
        panel_tag="[0, 100] mm, 0.5 mm bins", legend_loc="lower right",
        ax_ratio=cols[0][2], include_may30=True, include_latest_signal=True,
    )
    s_r = _draw_comparison_column(
        cols[1][0], cols[1][1], summary, signal_data,
        curve_key="dR1_zoom", xlabel=r"Track Separation ($\Delta R_1$) [mm]", xlim=(0.0, 20.0),
        panel_tag="[0, 20] mm, 0.5 mm bins", legend_loc="lower right",
        ax_ratio=cols[1][2], include_may30=True, include_latest_signal=True,
    )

    # Optional A' 110038 2-track efficiency comparison note if latest_georot signal exists
    sig_note = ""
    if signal_data:
        src_m = (signal_data.get("sources") or {}).get("master_dev")
        src_l = (signal_data.get("sources") or {}).get("latest_georot")
        if src_m and src_l:
            if "two_track_eff_vs_truth_fid_c2r" in src_m and "two_track_eff_vs_truth_fid_c2r" in src_l:
                sig_note = (
                    f" | A' 110038 eff (fid R<95, c2R): master_dev {src_m['two_track_eff_vs_truth_fid_c2r']*100:.2f}% "
                    f"vs Latest {src_l['two_track_eff_vs_truth_fid_c2r']*100:.2f}%"
                )
            else:
                sig_note = (
                    f" | A' 110038 2-trk eff (all): master_dev {src_m.get('two_track_eff_total', 0)*100:.2f}% "
                    f"vs Latest {src_l.get('two_track_eff_total', 0)*100:.2f}%"
                )

    pd_g = summary["paired_global_deltas"]["latest_georot_minus_baseline_jul8"]
    da_d = pd_g["All3Data"]["global"]["eff_pawan_c2R"]
    mc_d = pd_g["10004X"]["global"]["eff_pawan_c2R"]
    extra = (
        f"Global Δε (Latest-Jul8): Data {da_d['delta']*100:+.3f}±{da_d['delta_err']*100:.3f} pp, "
        f"MC {mc_d['delta']*100:+.3f}±{mc_d['delta_err']*100:.3f} pp"
    )
    add_figure_banners(
        fig,
        title="FIG09: Baseline-Jul8 vs Latest GeoRotation Overlay Efficiency (Delta R1, Station 3)",
        scope_line=get_scope_line(summary, "latest_georot", extra=extra),
        footer_line=f"Metric: eff_pawan_c2R (lt==2, chi2/nDoF<=25, R0,R1<95mm) | Frame: FASER | Cond: OFLCOND-FASER-05 | Ref: Slide 14{sig_note}",
        footer_fontsize=7.0,
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, "fig09_cmp_dR1_baseline_vs_latest")
    return {
        "fig_id": "fig09",
        "suite": "B",
        "title": "Baseline-Jul8 vs Latest GeoRotation: Delta R1",
        "pawan_ref": "Slide 14",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": {"left": s_l, "right": s_r, "global_delta_data_pp": da_d["delta"] * 100.0, "global_delta_mc_pp": mc_d["delta"] * 100.0},
    }


def _make_suite_b_pair(
    summary, signal_data, out_dir, fig_id, stem, title, pawan_ref, left_kw, right_kw
) -> Dict[str, Any]:
    apply_style()
    fig, cols = make_two_column_ratio_figure(figsize=(14.4, 6.2), height_ratios=(3.0, 1.2))
    s_l = _draw_comparison_column(cols[0][0], cols[0][1], summary, signal_data, **left_kw)
    s_r = _draw_comparison_column(cols[1][0], cols[1][1], summary, signal_data, **right_kw)
    add_figure_banners(
        fig,
        title=f"{fig_id.upper()}: {title}",
        scope_line=get_scope_line(summary, "latest_georot"),
        footer_line=f"{METRIC_FOOTNOTE} | Paired McNemar delta (Latest GeoRot - Baseline-Jul8) | Ref: Pawan {pawan_ref}",
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, stem)
    return {
        "fig_id": fig_id,
        "suite": "B",
        "title": title,
        "pawan_ref": pawan_ref,
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": {"left": s_l, "right": s_r},
    }


def make_fig10(summary, signal_data, out_dir, version, **_):
    return _make_suite_b_pair(
        summary, signal_data, out_dir,
        fig_id="fig10",
        stem="fig10_cmp_dR0_baseline_vs_latest",
        title="Baseline-Jul8 vs Latest GeoRotation: Station 1 Radial Separation (Delta R0)",
        pawan_ref="Slide 8 cmp",
        left_kw=dict(curve_key="dR0_wide", xlabel=r"Track Separation ($\Delta R_0$) [mm]", xlim=(0.0, 100.0), panel_tag="Wide [0, 100] mm", legend_loc="lower right"),
        right_kw=dict(curve_key="dR0_zoom", xlabel=r"Track Separation ($\Delta R_0$) [mm]", xlim=(0.0, 20.0), panel_tag="Zoom [0, 20] mm", legend_loc="lower right"),
    )


def make_fig11(summary, signal_data, out_dir, version, **_):
    return _make_suite_b_pair(
        summary, signal_data, out_dir,
        fig_id="fig11",
        stem="fig11_cmp_abs_dX_baseline_vs_latest",
        title="Baseline-Jul8 vs Latest GeoRotation: Horizontal Separation (|Delta X0|, |Delta X1|)",
        pawan_ref="Slides 11-12 cmp",
        left_kw=dict(curve_key="abs_dX0_wide", xlabel=r"Track Separation ($|\Delta X_0|$) [mm]", xlim=(0.0, 100.0), panel_tag="Station 1 [0, 100] mm", legend_loc="lower left"),
        right_kw=dict(curve_key="abs_dX1_wide", xlabel=r"Track Separation ($|\Delta X_1|$) [mm]", xlim=(0.0, 100.0), panel_tag="Station 3 [0, 100] mm", legend_loc="lower left"),
    )


def make_fig12(summary, signal_data, out_dir, version, **_):
    return _make_suite_b_pair(
        summary, signal_data, out_dir,
        fig_id="fig12",
        stem="fig12_cmp_abs_dY_baseline_vs_latest",
        title="Baseline-Jul8 vs Latest GeoRotation: Vertical Separation (|Delta Y0|, |Delta Y1|)",
        pawan_ref="Slides 9-10 cmp",
        left_kw=dict(curve_key="abs_dY0_wide", xlabel=r"Track Separation ($|\Delta Y_0|$) [mm]", xlim=(0.0, 100.0), panel_tag="Station 1 [0, 100] mm", legend_loc="lower right"),
        right_kw=dict(curve_key="abs_dY1_wide", xlabel=r"Track Separation ($|\Delta Y_1|$) [mm]", xlim=(0.0, 100.0), panel_tag="Station 3 [0, 100] mm", legend_loc="lower right"),
    )


def make_fig13(summary, signal_data, out_dir, version, **_):
    return _make_suite_b_pair(
        summary, signal_data, out_dir,
        fig_id="fig13",
        stem="fig13_cmp_signed_dY1_baseline_vs_latest",
        title="Baseline-Jul8 vs Latest GeoRotation: Signed Station 3 Vertical Separation (Delta Y1)",
        pawan_ref="Slide 13 cmp",
        left_kw=dict(
            curve_key="dY1_signed", xlabel=r"Track Separation ($\Delta Y_1 = y_1^+ - y_1^-$) [mm]",
            xlim=(-20.0, 20.0), slice_range=(-20.0, 20.0), panel_tag="Zoom [-20, +20] mm (1.0 mm bins)", legend_loc="lower right",
        ),
        right_kw=dict(
            curve_key="dY1_signed", xlabel=r"Track Separation ($\Delta Y_1 = y_1^+ - y_1^-$) [mm]",
            xlim=(-100.0, 100.0), panel_tag="Wide [-100, +100] mm (1.0 mm bins)", legend_loc="lower right",
        ),
    )


def make_fig14(summary, signal_data, out_dir, version, **_):
    """Baseline-Jul8 vs Latest GeoRotation: Track Pz0 spectrum + matched relative momentum residual."""
    apply_style()
    fig = plt.figure(figsize=(14.4, 6.2))
    outer = GridSpec(1, 2, figure=fig, left=0.065, right=0.96, top=0.855, bottom=0.115, wspace=0.26)
    inner_l = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0, 0], height_ratios=(2.8, 1.15), hspace=0.06)
    ax_tl = fig.add_subplot(inner_l[0, 0])
    ax_bl = fig.add_subplot(inner_l[1, 0], sharex=ax_tl)
    plt.setp(ax_tl.get_xticklabels(), visible=False)
    ax_r = fig.add_subplot(outer[0, 1])

    mom = summary["momentum"]
    ver_new, ver_old = "latest_georot", "baseline_jul8"
    edges = to_float_array(mom[ver_new]["10004X"]["track_pz0"]["edges"])
    centers_log = np.sqrt(edges[:-1] * edges[1:])

    mc_new = to_float_array(mom[ver_new]["10004X"]["track_pz0"]["norm"])
    da_new = to_float_array(mom[ver_new]["All3Data"]["track_pz0"]["norm"])
    mc_old = to_float_array(mom[ver_old]["10004X"]["track_pz0"]["norm"])
    da_old = to_float_array(mom[ver_old]["All3Data"]["track_pz0"]["norm"])

    ax_tl.step(edges, np.concatenate([mc_new, mc_new[-1:]]), where="post", color=COLORS["mc"], lw=1.4, label="Latest GeoRot MC (10004X)")
    ax_tl.step(edges, np.concatenate([mc_old, mc_old[-1:]]), where="post", color=COLORS["mc_base"], lw=1.1, ls="--", label="Baseline-Jul8 MC (10004X)")
    ax_tl.step(edges, np.concatenate([da_new, da_new[-1:]]), where="post", color=COLORS["data"], lw=1.4, label="Latest GeoRot Data (All3Data)")
    ax_tl.step(edges, np.concatenate([da_old, da_old[-1:]]), where="post", color=COLORS["data_base"], lw=1.1, ls="--", label="Baseline-Jul8 Data (All3Data)")
    ymax = max(0.08, float(np.nanmax(np.concatenate([mc_new, da_new, mc_old, da_old]))) * 1.18)
    ax_tl.set_xscale("log")
    ax_tl.set_xlim(1.0, 2.0e7)
    ax_tl.set_ylim(0.0, ymax)
    ax_tl.set_ylabel("Normalized Tracks / Bin")
    ax_tl.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax_tl.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.6, alpha=0.6)
    format_faser_header(ax_tl, subtitle="Track Pz0 Spectrum", right_text="Latest vs Baseline-Jul8")
    combine_legend(ax_tl, loc="upper left", fontsize=8.0)

    # Latest / Baseline-Jul8 ratio
    c_mc_o = to_int_array(mom[ver_old]["10004X"]["track_pz0"]["counts"])
    c_da_o = to_int_array(mom[ver_old]["All3Data"]["track_pz0"]["counts"])
    m_mc = (mc_old > 0) & (c_mc_o >= 25)
    m_da = (da_old > 0) & (c_da_o >= 25)
    r_mc = np.where(m_mc, mc_new / np.where(mc_old > 0, mc_old, 1.0), np.nan)
    r_da = np.where(m_da, da_new / np.where(da_old > 0, da_old, 1.0), np.nan)

    ax_bl.set_xscale("log")
    ax_bl.set_xlim(1.0, 2.0e7)
    ax_bl.set_ylim(0.65, 1.55)
    ax_bl.set_xlabel(r"Reconstructed Track $p_{z,0}$ [MeV]")
    ax_bl.set_ylabel("Latest / Jul8")
    ax_bl.yaxis.set_major_locator(MultipleLocator(0.20))
    ax_bl.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax_bl.axhline(1.0, color=COLORS["ref_line"], ls="--", lw=0.9)
    ax_bl.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax_bl.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.55, alpha=0.55)
    ax_bl.plot(centers_log[m_mc], r_mc[m_mc], "o-", color=COLORS["mc"], ms=2.6, lw=0.9, label="MC")
    ax_bl.plot(centers_log[m_da], r_da[m_da], "s-", color=COLORS["data"], ms=2.5, lw=0.9, label="Data")
    combine_legend(ax_bl, loc="lower left", ncol=2, fontsize=7.6)

    # Right panel: matched relative momentum residual (Track_pz0 / og_pz0 - 1) on log Y
    res_stats = {}
    for ver, ls, lbl_v in [(ver_new, "-", "Latest"), (ver_old, "--", "Jul8")]:
        for ds, col, lbl_d in [("10004X", COLORS["mc"], "MC"), ("All3Data", COLORS["data"], "Data")]:
            h = mom[ver][ds]["pz0_rel_res"]
            re = to_float_array(h["edges"])
            rn = to_float_array(h["norm"])
            safe = np.where(rn > 0, rn, 1e-7)
            ax_r.step(
                re, np.concatenate([safe, safe[-1:]]), where="post",
                color=col, linestyle=ls, linewidth=1.3 if ls == "-" else 1.05,
                label=f"{lbl_v} {lbl_d} (med={h.get('median', 0):+.4f}, IQR68={h.get('iqr68', 0):.4f})",
            )
            res_stats[f"{ver}_{ds}"] = {
                "total_matched": int(h.get("total_matched", 0)),
                "mean": float(h.get("mean", 0.0)),
                "median": float(h.get("median", 0.0)),
                "iqr68": float(h.get("iqr68", 0.0)),
                "std": float(h.get("std", 0.0)),
            }

    ax_r.set_xlim(-0.5, 0.5)
    ax_r.set_yscale("log")
    ax_r.set_ylim(2e-5, 1.0)
    ax_r.set_xlabel(r"Matched Relative Residual $(p_{z,0}^{\mathrm{reco}} / p_{z,0}^{\mathrm{og}} - 1)$")
    ax_r.set_ylabel("Normalized Matched Tracks / Bin")
    ax_r.axvline(0.0, color=COLORS["ref_line"], ls="--", lw=0.85)
    ax_r.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_r.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax_r.grid(axis="y", which="major", color=COLORS["grid"], ls=":", lw=0.6, alpha=0.6)
    format_faser_header(ax_r, subtitle="Matched Momentum Residual", right_text="[-0.5, +0.5], 0.01 bins")
    combine_legend(ax_r, loc="upper right", fontsize=7.6)

    add_figure_banners(
        fig,
        title="FIG14: Track Momentum Spectrum & Matched Relative Residual (Baseline-Jul8 vs Latest GeoRotation)",
        scope_line=get_scope_line(summary, ver_new),
        footer_line=f"Matched tracks within < 5 mm at St.1 & St.3 | {METRIC_FOOTNOTE} | Ref: Pawan Slide 16 cmp",
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, "fig14_cmp_momentum_and_residual")
    return {
        "fig_id": "fig14",
        "suite": "B",
        "title": "Momentum Spectrum & Relative Residual: Baseline-Jul8 vs Latest GeoRotation",
        "pawan_ref": "Slide 16 cmp",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": res_stats,
    }


def make_fig15(summary, signal_data, out_dir, version, **_):
    """Per-dataset efficiency summary + 5x5 multiplicity transition matrices (Data & MC)."""
    apply_style()
    fig = plt.figure(figsize=(15.2, 6.4))
    gs = GridSpec(1, 3, figure=fig, width_ratios=(1.55, 1.0, 1.0), left=0.115, right=0.975, top=0.845, bottom=0.125, wspace=0.30)
    ax_f = fig.add_subplot(gs[0, 0])
    ax_hd = fig.add_subplot(gs[0, 1])
    ax_hm = fig.add_subplot(gs[0, 2])

    datasets = ["All3Data", "008749", "010698", "Pos_008749-Neg_010698", "008749+010698", "10004X"]
    ds_labels = ["All3Data\n(242 pairs)", "008749\n(2022)", "010698\n(2023)", "Pos8749-\nNeg10698", "008749+\n010698", "10004X\n(MC)"]
    y_pos = np.arange(len(datasets))[::-1].astype(float)

    ver_specs = [
        ("baseline_may30", -0.22, "s", COLORS["slate"], "white", "Baseline-May30"),
        ("baseline_jul8", 0.0, "^", COLORS["data_base"], "white", "Baseline-Jul8"),
        ("latest_georot", +0.22, "o", COLORS["mc"], COLORS["mc"], "Latest GeoRotation"),
    ]
    all_effs = []
    summary_stats = {}
    for ver, dy, mk, col, mfc, lbl in ver_specs:
        if ver not in summary.get("global_summary", {}):
            continue
        vals, elo, ehi = [], [], []
        for ds in datasets:
            rec = summary["global_summary"][ver][ds]["eff_pawan_c2R"]
            vals.append(rec["eff"] * 100.0)
            elo.append(rec["err_lo"] * 100.0)
            ehi.append(rec["err_hi"] * 100.0)
            summary_stats[f"{ver}_{ds}"] = rec["eff"]
        all_effs.extend(vals)
        ax_f.errorbar(
            vals, y_pos + dy, xerr=np.vstack([elo, ehi]),
            fmt=mk, color=col, markerfacecolor=mfc, markeredgecolor=col,
            markersize=6.2, elinewidth=1.2, capsize=2.5, label=lbl, zorder=5,
        )

    span = (max(all_effs) - min(all_effs)) if all_effs else 2.0
    xmin = max(80.0, float(np.floor(min(all_effs) - max(0.6, 0.25 * span)))) if all_effs else 88.0
    xmax = min(100.5, float(np.ceil(max(all_effs) + max(1.8, 0.55 * span)))) if all_effs else 96.5
    ax_f.set_xlim(xmin, xmax)
    ax_f.set_ylim(-1.35, len(datasets) - 0.35)
    ax_f.set_yticks(y_pos)
    ax_f.set_yticklabels(ds_labels, fontsize=8.8)
    ax_f.set_xlabel("Global Overlay Efficiency eff_pawan_c2R [%]")
    ax_f.grid(axis="x", color=COLORS["grid"], ls=":", lw=0.7, alpha=0.7)
    ax_f.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    format_faser_header(ax_f, subtitle="Per-Dataset Eff.", right_text="322 pairs", fontsize=9.2, right_fontsize=7.8)

    # Annotate paired Delta (Latest - Jul8) on right edge of each row
    pair_key = "latest_georot_minus_baseline_jul8"
    p_glob = summary.get("paired_global_deltas", {}).get(pair_key, {})
    for idx, ds in enumerate(datasets):
        if ds in p_glob:
            drec = p_glob[ds]["global"]["eff_pawan_c2R"]
            txt = f"{drec['delta']*100:+.2f}±{drec['delta_err']*100:.2f} pp"
            ax_f.text(xmax - 0.08, y_pos[idx], txt, ha="right", va="center", fontsize=7.6, color="#222222",
                      bbox=dict(boxstyle="round,pad=0.18", facecolor="#f8f9fa", edgecolor="#dddddd", alpha=0.9))
    combine_legend(ax_f, loc="lower left", fontsize=7.8)

    # Heatmaps: row-normalized 5x5 migration matrices M[lt_jul8, lt_latest]
    tm = summary.get("transition_matrices", {}).get(pair_key, {})
    for ax_h, ds, sub_t in [(ax_hd, "All3Data", "All3Data (Row %)"), (ax_hm, "10004X", "10004X MC (Row %)")]:
        mat_raw = np.asarray(tm[ds]["global"], dtype=float)
        labels = tm[ds].get("labels", ["0", "1", "2", "3", "4+"])
        row_sums = np.maximum(mat_raw.sum(axis=1, keepdims=True), 1.0)
        mat_pct = 100.0 * mat_raw / row_sums
        ax_h.imshow(mat_pct, cmap="Blues", vmin=0.0, vmax=100.0, aspect="equal")
        ax_h.set_xticks(range(len(labels)))
        ax_h.set_yticks(range(len(labels)))
        ax_h.set_xticklabels(labels)
        ax_h.set_yticklabels(labels)
        ax_h.set_xlabel("Latest GeoRot longTracks")
        ax_h.set_ylabel("Baseline-Jul8 longTracks")
        format_faser_header(ax_h, subtitle=sub_t, right_text="", fontsize=9.0)
        for iy in range(mat_raw.shape[0]):
            for ix in range(mat_raw.shape[1]):
                pct = mat_pct[iy, ix]
                cnt = int(mat_raw[iy, ix])
                txt_col = "white" if pct > 55.0 else "#111111"
                cnt_str = f"{cnt/1e3:.1f}k" if cnt >= 10000 else f"{cnt:,}"
                ax_h.text(ix, iy, f"{pct:.1f}%\n({cnt_str})", ha="center", va="center", fontsize=6.8, color=txt_col)

    add_figure_banners(
        fig,
        title="FIG15: Per-Dataset Efficiency Breakdown & Multiplicity Migration (Baseline-Jul8 -> Latest GeoRotation)",
        scope_line=get_scope_line(summary, "latest_georot"),
        footer_line=f"{METRIC_FOOTNOTE} | Right labels in panel 1 show paired Δε = Latest - Jul8 [pp]",
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, "fig15_dataset_summary_and_migration")
    return {
        "fig_id": "fig15",
        "suite": "B",
        "title": "Per-Dataset Summary & 5x5 Multiplicity Migration",
        "pawan_ref": "Summary",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": summary_stats,
    }


# ---------------------------------------------------------------------------
# Suite C: Physics & Tracking Diagnostics (fig16 - fig20)
# ---------------------------------------------------------------------------
def make_fig16(summary, signal_data, out_dir, version, **_):
    """Multiplicity & efficiency-definition breakdown vs Delta R1 (small-dR vs plateau)."""
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(14.4, 7.6))
    fig.subplots_adjust(left=0.068, right=0.97, top=0.865, bottom=0.105, wspace=0.20, hspace=0.30)

    metric_specs = [
        ("eff_pawan_c2R", "eff_pawan_c2R (lt==2 & c2R)", DIAG_PALETTE[0], "o"),
        ("eff_lt2", "eff_lt2 (longTracks == 2)", DIAG_PALETTE[1], "s"),
        ("eff_ge2", "eff_ge2 (longTracks >= 2)", DIAG_PALETTE[2], "^"),
        ("eff_matched", "eff_matched (1-to-1 <5mm & c2R)", DIAG_PALETTE[3], "D"),
    ]

    # Top row: 4 efficiency metrics vs dR1 [0, 20] mm (Data left, MC right)
    for ax, ds, ds_lbl in [(axes[0, 0], "All3Data", "Data Overlay (All3Data)"), (axes[0, 1], "10004X", "Muon MC Overlay (10004X)")]:
        curves_m = summary["curves"][version][ds]["dR1_zoom"]
        arrs = []
        for mkey, mlbl, col, mk in metric_specs:
            e, v, ye, m = _bin_eff_stats(curves_m[mkey])
            arrs.append((v, m))
            plot_binned_errorbar(
                ax, e, v, ye, mask=m, label=mlbl, color=col, marker=mk,
                markersize=3.4, elinewidth=0.85, show_xerr=True,
            )
        ylim = choose_eff_ylim(arrs, default=(0.80, 1.00), min_floor=0.65)
        configure_efficiency_axis(ax, xlim=(0.0, 20.0), ylim=ylim, has_right_twin=False)
        ax.set_xlabel(r"Track Separation ($\Delta R_1$) [mm]")
        format_faser_header(ax, subtitle=f"Metrics: {ds}", right_text="[0, 20] mm", fontsize=9.4, right_fontsize=8.0)
        combine_legend(ax, loc="lower right", fontsize=7.8)

    # Bottom row: Multiplicity failure fractions P(lt=0, 1, >=3, lt=2 fail c2R) vs dR1 [0, 20] and [0, 100]
    mult_specs = [
        ("p_lt1", "lt == 1 (1 track lost)", DIAG_PALETTE[0]),
        ("p_lt_ge3", "lt >= 3 (split / ghost)", DIAG_PALETTE[1]),
        ("p_lt2_fail_c2R", "lt == 2 & fail chi2/R", DIAG_PALETTE[4]),
        ("p_lt0", "lt == 0 (both lost)", DIAG_PALETTE[6]),
    ]
    for ax, ckey, xlim, rtag in [
        (axes[1, 0], "dR1_zoom", (0.0, 20.0), "Zoom [0, 20] mm"),
        (axes[1, 1], "dR1_wide", (0.0, 100.0), "Wide [0, 100] mm"),
    ]:
        md_da = summary["multiplicity_curves"][version]["All3Data"][ckey]
        md_mc = summary["multiplicity_curves"][version]["10004X"][ckey]
        e = to_float_array(md_da["edges"])
        c = 0.5 * (e[:-1] + e[1:])
        den_d = to_int_array(md_da["den"]) > 0
        den_m = to_int_array(md_mc["den"]) > 0
        ymax = 0.08
        for pkey, plbl, col in mult_specs:
            vd = to_float_array(md_da[pkey]) * 100.0
            vm = to_float_array(md_mc[pkey]) * 100.0
            if np.any(den_d):
                ymax = max(ymax, float(np.nanpercentile(vd[den_d], 99)) * 1.18)
                ax.plot(c[den_d], vd[den_d], color=col, lw=1.45, ls="-", label=f"Data: {plbl}")
            if np.any(den_m):
                ymax = max(ymax, float(np.nanpercentile(vm[den_m], 99)) * 1.18)
                ax.plot(c[den_m], vm[den_m], color=col, lw=1.15, ls="--", alpha=0.85, label=f"MC: {plbl}")
        ax.set_xlim(*xlim)
        ax.set_ylim(0.0, min(35.0, max(12.5, ymax * 1.32)))
        ax.set_xlabel(r"Track Separation ($\Delta R_1$) [mm]")
        ax.set_ylabel("Failure Mode Fraction [%]")
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)
        ax.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.6, alpha=0.65)
        format_faser_header(ax, subtitle="Failure Modes (Solid=Data, Dash=MC)", right_text=rtag, fontsize=9.2, right_fontsize=7.8)
        combine_legend(ax, loc="upper right", ncol=2, fontsize=6.9)

    mode_tag = "" if version == "latest_georot" else " [Baseline-Jul8 Preview]"
    add_figure_banners(
        fig,
        title=f"FIG16: Small-Separation Metric Comparison & Multiplicity Failure Breakdown{mode_tag}",
        scope_line=get_scope_line(summary, version),
        footer_line=f"{METRIC_FOOTNOTE} | Isolates split/ghost (lt>=3) at small Delta R1 vs single-track loss (lt==1) at plateau",
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, "fig16_multiplicity_and_metric_breakdown")
    gs = summary["global_summary"][version]
    return {
        "fig_id": "fig16",
        "suite": "C",
        "title": f"Multiplicity & Metric Breakdown{mode_tag}",
        "pawan_ref": "Small-dR diagnostic",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": {
            "version": version,
            "data_small_dR1_lt5": gs["All3Data"]["regimes"]["small_dR1_lt5"],
            "mc_small_dR1_lt5": gs["10004X"]["regimes"]["small_dR1_lt5"],
        },
    }


def make_fig17(summary, signal_data, out_dir, version, **_):
    """Mirrors Slide 15 diagnostic: Plateau (dR1 in [30, 70] mm) stage-of-loss funnel & dependencies."""
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(14.6, 7.8))
    fig.subplots_adjust(left=0.095, right=0.97, top=0.865, bottom=0.105, wspace=0.22, hspace=0.31)

    sol = summary["stage_of_loss"][version]

    # Top-left: Stage-of-loss inefficiency breakdown in plateau dR1 in [30, 70] mm
    ax0 = axes[0, 0]
    cats = [
        ("lt1_miss_seg_st1", "lt=1: Miss Seg St1"),
        ("lt1_miss_seg_st2", "lt=1: Miss Seg St2"),
        ("lt1_miss_seg_st3", "lt=1: Miss Seg St3"),
        ("lt1_miss_seg_multi", "lt=1: Miss >=2 Stns"),
        ("lt1_killed_by_gb", "lt=1: Killed by GB"),
        ("lt1_all3gb_lost_in_ckf", "lt=1: Lost in CKF"),
        ("lt1_unmatched_kept", "lt=1: Other/Unmatch"),
        ("lt0", "lt == 0 (Both lost)"),
        ("lt_ge3", "lt >= 3 (Ghost/Split)"),
        ("lt2_fail_c2R", "lt == 2 Fail chi2/R"),
    ]
    da_p = sol["All3Data"]["regimes"]["plateau_dR1_30_70"]
    mc_p = sol["10004X"]["regimes"]["plateau_dR1_30_70"]
    da_den = max(1, int(da_p["n_den"]))
    mc_den = max(1, int(mc_p["n_den"]))
    da_pct = [100.0 * int(da_p.get(k, 0)) / da_den for k, _ in cats]
    mc_pct = [100.0 * int(mc_p.get(k, 0)) / mc_den for k, _ in cats]

    y = np.arange(len(cats))[::-1].astype(float)
    h = 0.36
    ax0.barh(y + h / 2, da_pct, height=h, color=COLORS["data"], alpha=0.88, label=f"Data (tot loss {100.0*(1-da_p['pass_c2R']/da_den):.2f}%)")
    ax0.barh(y - h / 2, mc_pct, height=h, color=COLORS["mc"], alpha=0.88, label=f"MC (tot loss {100.0*(1-mc_p['pass_c2R']/mc_den):.2f}%)")
    xmax0 = max(2.0, max(da_pct + mc_pct) * 1.32)
    ax0.set_xlim(0.0, xmax0)
    ax0.set_yticks(y)
    ax0.set_yticklabels([lbl for _, lbl in cats], fontsize=8.0)
    ax0.set_xlabel("Fraction of Plateau Events [%]")
    ax0.grid(axis="x", color=COLORS["grid"], ls=":", lw=0.6, alpha=0.65)
    ax0.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    for i, (dv, mv) in enumerate(zip(da_pct, mc_pct)):
        if dv >= 0.05 or mv >= 0.05:
            ax0.text(max(dv, mv) + xmax0 * 0.015, y[i], f"{dv:.2f}% / {mv:.2f}%", va="center", ha="left", fontsize=6.9, color="#222222")
    format_faser_header(ax0, subtitle="Plateau Stage-of-Loss Funnel", right_text=r"$\Delta R_1 \in [30, 70]$ mm")
    combine_legend(ax0, loc="lower right", fontsize=7.6)

    # Top-right: Plateau eff vs nSeg and nGBSeg (integer bin index = edges[:-1])
    ax1 = axes[0, 1]
    for ds, col, lbl_d in [("All3Data", COLORS["data"], "Data"), ("10004X", COLORS["mc"], "MC")]:
        e_s, v_s, ye_s, m_s = _bin_eff_stats(sol[ds]["plateau_vs_nSeg"])
        e_g, v_g, ye_g, m_g = _bin_eff_stats(sol[ds]["plateau_vs_nGBSeg"])
        # Only plot bins with >= 50 events to avoid noisy 0-count edges; shift -0.5 so integer N sits at N
        m_s = m_s & (to_int_array(sol[ds]["plateau_vs_nSeg"]["den"]) >= 50)
        m_g = m_g & (to_int_array(sol[ds]["plateau_vs_nGBSeg"]["den"]) >= 50)
        plot_binned_errorbar(ax1, e_s - 0.5, v_s, ye_s, mask=m_s, label=f"{lbl_d}: nSeg (SegFit)", color=col, marker="o", markersize=4.6, show_xerr=False, x_offset=-0.12)
        plot_binned_errorbar(ax1, e_g - 0.5, v_g, ye_g, mask=m_g, label=f"{lbl_d}: nGBSeg (GhostBust)", color=col, marker="^", markerfacecolor="white", markersize=4.9, show_xerr=False, x_offset=+0.12)
    ax1.set_xlim(4.4, 12.6)
    ax1.set_ylim(-0.04, 1.06)
    ax1.set_xlabel("Total 3-Station Segments in Event (nSeg or nGBSeg)")
    ax1.set_ylabel("Plateau Overlay Efficiency")
    ax1.xaxis.set_major_locator(MultipleLocator(1.0))
    ax1.yaxis.set_major_locator(MultipleLocator(0.2))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax1.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax1.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.6, alpha=0.65)
    format_faser_header(ax1, subtitle="Plateau Eff. vs Segment Count", right_text="bins with N >= 50", fontsize=9.4, right_fontsize=8.0)
    combine_legend(ax1, loc="lower right", ncol=2, fontsize=7.4)

    # Bottom-left: Plateau eff vs min(og_p0[0], og_p0[1]) [GeV]
    ax2 = axes[1, 0]
    p_arrs = []
    for ds, col, lbl_d in [("10004X", COLORS["mc"], "MC (10004X)"), ("All3Data", COLORS["data"], "Data (All3Data)")]:
        c_p = sol[ds]["plateau_vs_min_og_p0"]
        e_p, v_p, ye_p, m_p = _bin_eff_stats(c_p)
        m_p = m_p & (to_int_array(c_p["den"]) >= 30)
        p_arrs.append((v_p, m_p))
        centers_log = np.sqrt(e_p[:-1] * e_p[1:])
        ax2.errorbar(
            centers_log[m_p], v_p[m_p], yerr=ye_p[:, m_p],
            fmt="o-", color=col, markersize=4.0, lw=1.0, elinewidth=0.9, label=lbl_d,
        )
    ax2.set_xscale("log")
    ax2.set_xlim(50.0, 5000.0)
    ax2.set_ylim(*choose_eff_ylim(p_arrs, default=(0.80, 1.00), min_floor=0.65))
    ax2.set_xlabel(r"$\min(p_{0}^{\mathrm{og},+}, p_{0}^{\mathrm{og},-})$ [GeV]")
    ax2.set_ylabel("Plateau Overlay Efficiency")
    ax2.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax2.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.6, alpha=0.65)
    format_faser_header(ax2, subtitle="Plateau Efficiency vs Softer Track Momentum", right_text="[50, 5000] GeV")
    combine_legend(ax2, loc="lower right", fontsize=8.0)

    # Bottom-right: Plateau eff vs max(og_R0) [mm]
    ax3 = axes[1, 1]
    r_arrs = []
    for ds, col, lbl_d in [("10004X", COLORS["mc"], "MC (10004X)"), ("All3Data", COLORS["data"], "Data (All3Data)")]:
        c_r = sol[ds]["plateau_vs_max_og_R0"]
        e_r, v_r, ye_r, m_r = _bin_eff_stats(c_r)
        m_r = m_r & (to_int_array(c_r["den"]) >= 30)
        r_arrs.append((v_r, m_r))
        plot_binned_errorbar(ax3, e_r, v_r, ye_r, mask=m_r, label=lbl_d, color=col, marker="o", markersize=4.0)
    ax3.set_xlim(0.0, 95.0)
    ax3.set_ylim(*choose_eff_ylim(r_arrs, default=(0.80, 1.00), min_floor=0.65))
    ax3.set_xlabel(r"$\max(R_0^{\mathrm{og},+}, R_0^{\mathrm{og},-})$ [mm]")
    ax3.set_ylabel("Plateau Overlay Efficiency")
    ax3.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax3.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.6, alpha=0.65)
    format_faser_header(ax3, subtitle="Plateau Efficiency vs Outer Track Radius", right_text="[0, 95] mm")
    combine_legend(ax3, loc="lower left", fontsize=8.0)

    mode_tag = "" if version == "latest_georot" else " [Baseline-Jul8 Preview]"
    add_figure_banners(
        fig,
        title=f"FIG17: Plateau Efficiency Gap Diagnostic (SegmentFit -> GhostBusters -> CKF Funnel){mode_tag}",
        scope_line=get_scope_line(summary, version),
        footer_line=f"{METRIC_FOOTNOTE} | Plateau: 30 <= Delta R1 <= 70 mm | Ref: Pawan Slide 15",
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, "fig17_plateau_stage_of_loss")
    return {
        "fig_id": "fig17",
        "suite": "C",
        "title": f"Plateau Stage-of-Loss Funnel & Dependencies{mode_tag}",
        "pawan_ref": "Slide 15",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": {"version": version, "plateau_data": da_p, "plateau_mc": mc_p},
    }


def make_fig18(summary, signal_data, out_dir, version, **_):
    """Mirrors Slides 11-12 diagnostic: Large-|dX0| stratification + stage fractions + 2D P_loss maps."""
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(14.6, 8.2))
    fig.subplots_adjust(left=0.068, right=0.955, top=0.87, bottom=0.095, wspace=0.22, hspace=0.29)

    ldx = summary["large_dX_diagnostics"][version]
    sol = summary["stage_of_loss"][version]

    # Top-left: eff(|dX0|) stratified by max(|og_x0|)
    ax0 = axes[0, 0]
    strata = [
        ("0_33", r"$\max(|x_0^{\mathrm{og}}|) \in [0, 33]$ mm", DIAG_PALETTE[0]),
        ("33_66", r"$\max(|x_0^{\mathrm{og}}|) \in [33, 66]$ mm", DIAG_PALETTE[2]),
        ("66_95", r"$\max(|x_0^{\mathrm{og}}|) \in [66, 95]$ mm", DIAG_PALETTE[1]),
    ]
    s_arrs = []
    for skey, slbl, col in strata:
        for ds, mk, mfc, tag in [("All3Data", "o", col, "Data"), ("10004X", "^", "white", "MC")]:
            c = ldx[ds]["eff_abs_dX0_by_max_abs_og_x0"][skey]
            e, v, ye, m = _bin_eff_stats(c)
            m = m & (to_int_array(c["den"]) >= 50)
            s_arrs.append((v, m))
            plot_binned_errorbar(
                ax0, e, v, ye, mask=m, label=f"{tag}: {slbl}", color=col,
                marker=mk, markerfacecolor=mfc, markersize=4.2, elinewidth=0.9,
            )
    ax0.set_xlim(0.0, 100.0)
    ylim_dx = choose_eff_ylim(s_arrs, default=(0.75, 1.00), min_floor=0.65)
    ax0.set_ylim(min(ylim_dx[0], 0.74), 1.00)
    ax0.set_xlabel(r"Track Separation ($|\Delta X_0|$) [mm]")
    ax0.set_ylabel("Overlay Efficiency")
    ax0.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax0.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.6, alpha=0.65)
    format_faser_header(ax0, subtitle=r"Eff. vs $|\Delta X_0|$ by $\max(|x_0^{\mathrm{og}}|)$", right_text="5 mm bins", fontsize=9.4, right_fontsize=8.0)
    combine_legend(ax0, loc="lower left", ncol=2, fontsize=7.0)

    # Top-right: Loss stage fractions vs |dX0|
    ax1 = axes[0, 1]
    sf_d = sol["All3Data"]["abs_dX0_stage_fractions"]
    sf_m = sol["10004X"]["abs_dX0_stage_fractions"]
    e1 = to_float_array(sf_d["edges"])
    c1 = 0.5 * (e1[:-1] + e1[1:])
    m1 = to_int_array(sf_d["den"]) >= 50
    stage_lines = [
        ("frac_lt1_miss_seg", "lt=1: Miss SegFit", DIAG_PALETTE[0]),
        ("frac_lt1_killed_gb", "lt=1: Killed by GB", DIAG_PALETTE[1]),
        ("frac_lt1_lost_ckf", "lt=1: Lost in CKF", DIAG_PALETTE[2]),
        ("frac_lt_ge3", "lt >= 3 (Ghost/Split)", DIAG_PALETTE[3]),
        ("frac_lt2_fail_c2R", "lt == 2 Fail chi2/R", DIAG_PALETTE[4]),
        ("frac_lt0", "lt == 0", DIAG_PALETTE[6]),
    ]
    ymax1 = 5.0
    for fkey, flbl, col in stage_lines:
        yd = to_float_array(sf_d[fkey]) * 100.0
        if np.any(m1):
            ymax1 = max(ymax1, float(np.nanmax(yd[m1])) * 1.25)
            ax1.plot(c1[m1], yd[m1], "o-", color=col, ms=3.4, lw=1.35, label=f"Data: {flbl}")
    # Total loss comparison (Data vs MC)
    tot_d = (1.0 - to_float_array(sf_d["frac_pass_c2R"])) * 100.0
    tot_m = (1.0 - to_float_array(sf_m["frac_pass_c2R"])) * 100.0
    m1m = to_int_array(sf_m["den"]) >= 50
    if np.any(m1):
        ymax1 = max(ymax1, float(np.nanmax(tot_d[m1])) * 1.18)
        ax1.plot(c1[m1], tot_d[m1], color="#111111", lw=1.6, ls="-", label="Data: Total loss")
    if np.any(m1m):
        ymax1 = max(ymax1, float(np.nanmax(tot_m[m1m])) * 1.18)
        ax1.plot(c1[m1m], tot_m[m1m], color=COLORS["mc"], lw=1.4, ls="--", label="MC: Total loss")
    ax1.set_xlim(0.0, 100.0)
    ax1.set_ylim(0.0, min(25.0, max(10.0, ymax1)))
    ax1.set_xlabel(r"Track Separation ($|\Delta X_0|$) [mm]")
    ax1.set_ylabel("Inefficiency Contribution [%]")
    ax1.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax1.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.6, alpha=0.65)
    format_faser_header(ax1, subtitle=r"Failure Stage Contributions vs $|\Delta X_0|$", right_text="5 mm bins")
    combine_legend(ax1, loc="upper left", ncol=2, fontsize=6.9)

    # Bottom row: 2D single-track loss probability maps P_loss(og_x0, og_y0) at Station 1 (Data & MC)
    maps = [
        (axes[1, 0], "All3Data", "Data Overlay (All3Data)"),
        (axes[1, 1], "10004X", "Muon MC Overlay (10004X)"),
    ]
    # Shared color scale across Data & MC
    vmax_2d = 0.08
    for _, ds, _ in maps:
        lm = ldx[ds]["loss_map_st1"]
        den2 = np.asarray(lm["den"], dtype=float)
        pl2 = np.asarray(lm["p_loss"], dtype=float)
        good = (den2 >= 50) & np.isfinite(pl2)
        if np.any(good):
            vmax_2d = max(vmax_2d, float(np.nanpercentile(pl2[good], 97.5)))
    vmax_2d = min(0.25, max(0.08, round(vmax_2d * 100.0) / 100.0))

    cmap = plt.get_cmap("inferno_r").copy()
    cmap.set_bad("#f0f0f0")
    map_stats = {}
    for ax_m, ds, ds_lbl in maps:
        lm = ldx[ds]["loss_map_st1"]
        xe = to_float_array(lm["x_edges"])
        ye = to_float_array(lm["y_edges"])
        den2 = np.asarray(lm["den"], dtype=float)
        pl2 = np.asarray(lm["p_loss"], dtype=float)
        masked = np.ma.masked_where((den2 < 20) | (~np.isfinite(pl2)), pl2 * 100.0)
        pcm = ax_m.pcolormesh(xe, ye, masked, cmap=cmap, vmin=0.0, vmax=vmax_2d * 100.0, shading="flat")
        cb = fig.colorbar(pcm, ax=ax_m, pad=0.02, fraction=0.046)
        cb.set_label("Single-Track Loss Rate [%]", fontsize=8.5)
        ax_m.add_patch(Circle((0.0, 0.0), 95.0, fill=False, edgecolor="#111111", lw=1.0, ls="--"))
        for xv in (-66.0, -33.0, 33.0, 66.0):
            ax_m.axvline(xv, color="#444444", lw=0.65, ls=":", alpha=0.65)
        ax_m.set_xlim(-100.0, 100.0)
        ax_m.set_ylim(-100.0, 100.0)
        ax_m.set_aspect("equal", adjustable="box")
        ax_m.set_xlabel(r"Original Track $x_0^{\mathrm{og}}$ [mm] (Station 1)")
        ax_m.set_ylabel(r"Original Track $y_0^{\mathrm{og}}$ [mm]")
        ax_m.tick_params(axis="both", which="both", direction="in", top=True, right=True)
        format_faser_header(ax_m, subtitle=f"St.1 Loss: {ds} (40x40)", right_text="", fontsize=9.0)
        good = (den2 >= 20) & np.isfinite(pl2)
        map_stats[ds] = {
            "median_loss_pct": float(np.nanmedian(pl2[good]) * 100.0) if np.any(good) else None,
            "p95_loss_pct": float(np.nanpercentile(pl2[good], 95.0) * 100.0) if np.any(good) else None,
        }

    mode_tag = "" if version == "latest_georot" else " [Baseline-Jul8 Preview]"
    add_figure_banners(
        fig,
        title=f"FIG18: Large-|Delta X0| Acceptance Stratification & 2D Single-Track Loss Maps{mode_tag}",
        scope_line=get_scope_line(summary, version),
        footer_line=f"{METRIC_FOOTNOTE} | Dashed circle: R0 = 95 mm; dotted lines: |x0| = 33, 66 mm | Ref: Pawan Slides 11-12",
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, "fig18_large_dX_and_2d_loss_maps")
    return {
        "fig_id": "fig18",
        "suite": "C",
        "title": f"Large-|dX| Diagnostic & 2D Loss Maps{mode_tag}",
        "pawan_ref": "Slides 11-12",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": {"version": version, "loss_map_st1": map_stats},
    }


def make_fig19(summary, signal_data, out_dir, version, **_):
    """Mirrors Slide 16 reweighting test: og_p0 closure + eff(dR1) & MC/Data before/after reweighting."""
    apply_style()
    fig, cols = make_two_column_ratio_figure(figsize=(14.4, 6.3), height_ratios=(2.8, 1.25), wspace=0.28)
    ax_tl, ax_bl = cols[0]
    ax_tr, ax_br = cols[1]

    rw = summary["momentum_reweighting"]
    bins_info = rw["og_p0_bins_gev"]
    rw_v = rw[version]

    edges_p = to_float_array(bins_info["edges"])
    centers_p = np.sqrt(edges_p[:-1] * edges_p[1:])
    f_da = to_float_array(bins_info["data_frac"])
    f_mc = to_float_array(bins_info["mc_frac"])
    w_p = to_float_array(bins_info["weights"])
    f_rw = to_float_array(rw_v["closure_mc_reweighted_frac"])

    # Left top: og_p0 closure
    ax_tl.step(edges_p, np.concatenate([f_da, f_da[-1:]]), where="post", color=COLORS["data"], lw=1.6, label="All3Data og_p0")
    ax_tl.step(edges_p, np.concatenate([f_mc, f_mc[-1:]]), where="post", color=COLORS["mc"], lw=1.3, ls="--", label="10004X MC Unweighted")
    ax_tl.plot(centers_p, f_rw, "o", color=COLORS["signal"], ms=3.5, label="10004X MC Reweighted (Closure)")
    ymax_p = max(0.05, float(np.nanmax(np.concatenate([f_da, f_mc, f_rw]))) * 1.22)
    ax_tl.set_xscale("log")
    ax_tl.set_xlim(50.0, 5000.0)
    ax_tl.set_ylim(0.0, ymax_p)
    ax_tl.set_ylabel("Normalized Single-Track Density")
    ax_tl.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax_tl.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.6, alpha=0.65)
    format_faser_header(ax_tl, subtitle="Momentum Reweighting Closure", right_text="40 log bins [50, 5000] GeV")
    combine_legend(ax_tl, loc="upper right", fontsize=8.0)

    n_ev = int(rw_v.get("n_events", 0))
    n_eff = float(rw_v.get("n_eff", 0.0))
    eff_frac = 100.0 * n_eff / max(1, n_ev)
    ax_tl.text(
        0.03, 0.86,
        f"W_evt = w(p0+) * w(p0-)\nN_MC = {n_ev:,}\nN_eff = {n_eff:,.0f} ({eff_frac:.1f}%)",
        transform=ax_tl.transAxes, fontsize=7.8, va="top", color="#222222",
        bbox=dict(boxstyle="round,pad=0.32", facecolor="white", edgecolor="#dddddd", alpha=0.92),
    )

    # Left bottom: weight function w(p)
    ax_bl.set_xscale("log")
    ax_bl.set_yscale("log")
    ax_bl.set_xlim(50.0, 5000.0)
    valid_w = np.isfinite(w_p) & (w_p > 0)
    wmin = max(0.05, float(np.nanmin(w_p[valid_w])) * 0.7) if np.any(valid_w) else 0.1
    wmax = min(50.0, float(np.nanmax(w_p[valid_w])) * 1.4) if np.any(valid_w) else 10.0
    ax_bl.set_ylim(wmin, wmax)
    ax_bl.set_xlabel(r"Original Single-Track Momentum $p_0^{\mathrm{og}}$ [GeV]")
    ax_bl.set_ylabel(r"Weight $w(p) = f_{\mathrm{Data}}/f_{\mathrm{MC}}$")
    ax_bl.axhline(1.0, color=COLORS["ref_line"], ls="--", lw=0.9)
    ax_bl.plot(centers_p[valid_w], w_p[valid_w], "o-", color="#222222", ms=3.2, lw=1.1)
    ax_bl.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax_bl.grid(axis="y", which="major", color=COLORS["grid"], ls=":", lw=0.55, alpha=0.6)

    # Right top: eff(dR1) before vs after MC->Data momentum reweighting
    mc_raw = summary["curves"][version]["10004X"]["dR1_wide"]["eff_pawan_c2R"]
    da_raw = summary["curves"][version]["All3Data"]["dR1_wide"]["eff_pawan_c2R"]
    mc_rw = rw_v["curves"]["dR1_wide"]

    e_d, v_d, ye_d, m_d = _bin_eff_stats(da_raw)
    e_m, v_m, ye_m, m_m = _bin_eff_stats(mc_raw)
    e_rw = to_float_array(mc_rw["edges"])
    v_rw = to_float_array(mc_rw["eff"])
    ye_rw = to_float_array(mc_rw["err"])
    m_rw = (to_float_array(mc_rw["sumw_den"]) > 0) & np.isfinite(v_rw)

    ylim_r = choose_eff_ylim([(v_d, m_d), (v_m, m_m), (v_rw, m_rw)], default=(0.80, 1.00), min_floor=0.65)
    configure_efficiency_axis(ax_tr, xlim=(0.0, 100.0), ylim=ylim_r, has_right_twin=False)
    format_faser_header(ax_tr, subtitle="MC->Data Reweighting Impact", right_text=r"$\Delta R_1$ [0, 100] mm", fontsize=9.5, right_fontsize=8.0)

    plot_binned_errorbar(ax_tr, e_m, v_m, ye_m, mask=m_m, label="10004X MC (Unweighted)", color=COLORS["mc"], marker="^", markerfacecolor="white", markersize=2.8, elinewidth=0.72, alpha=0.82, zorder=4)
    plot_binned_errorbar(ax_tr, e_rw, v_rw, ye_rw, mask=m_rw, label="10004X MC (p0-Reweighted)", color=DIAG_PALETTE[2], marker="o", markersize=2.6, elinewidth=0.78, zorder=6)
    plot_binned_errorbar(ax_tr, e_d, v_d, ye_d, mask=m_d, label="All3Data Overlay", color=COLORS["data"], marker="o", markersize=2.6, elinewidth=0.75, zorder=5)
    combine_legend(ax_tr, loc="lower right", fontsize=7.9)

    # Right bottom: MC/Data ratio before vs after reweighting
    r_before = summary["mc_data_ratios"][version]["dR1_wide"]["eff_pawan_c2R"]
    r_after = rw_v["mc_data_ratios_reweighted"]["dR1_wide"]
    e_rb, rb, yerb, mrb = _ratio_stats(r_before)
    e_ra, ra, yera, mra = _ratio_stats(r_after)

    configure_ratio_axis(ax_br, xlim=(0.0, 100.0), xlabel=r"Track Separation ($\Delta R_1$) [mm]", ylabel="MC / Data", ylim=(0.95, 1.05), ref_value=1.0)
    plot_binned_errorbar(ax_br, e_rb, rb, yerb, mask=mrb, label="Unweighted MC / Data", color="#777777", marker="^", markerfacecolor="white", markersize=2.5, elinewidth=0.7, alpha=0.85, zorder=4)
    plot_binned_errorbar(ax_br, e_ra, ra, yera, mask=mra, label="Reweighted MC / Data", color=DIAG_PALETTE[2], marker="o", markersize=2.5, elinewidth=0.78, zorder=6)
    combine_legend(ax_br, loc="lower right", ncol=2, fontsize=7.4)

    # Plateau mean MC/Data before vs after
    cen_r = 0.5 * (e_rb[:-1] + e_rb[1:])
    plat = (cen_r >= 30.0) & (cen_r <= 70.0)
    mean_rb = float(np.nanmean(rb[plat & mrb])) if np.any(plat & mrb) else np.nan
    mean_ra = float(np.nanmean(ra[plat & mra])) if np.any(plat & mra) else np.nan
    ax_br.text(
        0.02, 0.87,
        f"Plateau [30,70] mm mean MC/Data: {mean_rb:.4f} -> {mean_ra:.4f}",
        transform=ax_br.transAxes, fontsize=7.4, color="#222222",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#dddddd", alpha=0.90),
    )

    mode_tag = "" if version == "latest_georot" else " [Baseline-Jul8 Preview]"
    add_figure_banners(
        fig,
        title=f"FIG19: Two-Pass MC->Data Momentum Reweighting Closure & MC/Data Impact{mode_tag}",
        scope_line=get_scope_line(summary, version),
        footer_line=f"{METRIC_FOOTNOTE} | Weighted sumw/sumw2 variance | Ref: Pawan Slide 16",
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, "fig19_momentum_reweighting")
    return {
        "fig_id": "fig19",
        "suite": "C",
        "title": f"Momentum Reweighting Closure & MC/Data Impact{mode_tag}",
        "pawan_ref": "Slide 16",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": {
            "version": version,
            "n_events": n_ev,
            "n_eff": n_eff,
            "plateau_mc_data_before": mean_rb,
            "plateau_mc_data_after": mean_ra,
        },
    }


def make_fig20(summary, signal_data, out_dir, version, **_):
    """Direct validation: Baseline-May30 (OverlaidNTPs_bug) vs Baseline-Jul8 (OverlaidNTPs)."""
    apply_style()
    fig = plt.figure(figsize=(15.4, 6.2))
    outer = GridSpec(1, 3, figure=fig, left=0.055, right=0.975, top=0.85, bottom=0.115, wspace=0.27)
    cols = []
    for cidx in range(3):
        inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0, cidx], height_ratios=(2.8, 1.2), hspace=0.06)
        ax_t = fig.add_subplot(inner[0, 0])
        ax_b = fig.add_subplot(inner[1, 0], sharex=ax_t)
        plt.setp(ax_t.get_xticklabels(), visible=False)
        cols.append((ax_t, ax_b))

    pair_key = "baseline_jul8_minus_baseline_may30"
    v_a, v_b = "baseline_may30", "baseline_jul8"

    # Columns 0 & 1: dR1_zoom [0, 20] and dY1_signed [-20, 20]
    for col_idx, (ckey, xl, xlim, sr, tag, lloc) in enumerate([
        ("dR1_zoom", r"Track Separation ($\Delta R_1$) [mm]", (0.0, 20.0), None, "Slide 7 zoom", "lower right"),
        ("dY1_signed", r"Track Separation ($\Delta Y_1$) [mm]", (-20.0, 20.0), (-20.0, 20.0), "Slide 13 zoom", "lower left"),
    ]):
        ax_t, ax_b = cols[col_idx]
        arrs = []
        for ver, mk, mfc, vlbl in [(v_a, "^", "white", "May30"), (v_b, "o", None, "Jul8")]:
            for ds, col, dlbl in [("10004X", COLORS["mc"], "MC"), ("All3Data", COLORS["data"], "Data")]:
                e, v, ye, m = _bin_eff_stats(summary["curves"][ver][ds][ckey]["eff_pawan_c2R"])
                if sr is not None:
                    e, [v, ye0, ye1, mf] = slice_1d_range(e, [v, ye[0], ye[1], m.astype(float)], *sr)
                    ye, m = np.vstack([ye0, ye1]), mf > 0.5
                arrs.append((v, m))
                plot_binned_errorbar(
                    ax_t, e, v, ye, mask=m, label=f"{vlbl} {dlbl}",
                    color=col, marker=mk, markerfacecolor=mfc, markersize=3.5, elinewidth=0.8,
                )
        ylim = choose_eff_ylim(arrs, default=(0.80, 1.00), min_floor=0.65)
        configure_efficiency_axis(ax_t, xlim=xlim, ylim=ylim, has_right_twin=False)
        format_faser_header(ax_t, subtitle="May30 vs Jul8", right_text=tag, fontsize=8.8, right_fontsize=7.6)
        combine_legend(ax_t, loc=lloc, ncol=2, fontsize=7.2)

        # Bottom paired McNemar delta (Jul8 - May30) [pp]
        e_pm, dm, dem, mm = _paired_delta_pp(summary["paired_curves"][pair_key]["10004X"][ckey]["eff_pawan_c2R"])
        e_pd, dd, ded, md = _paired_delta_pp(summary["paired_curves"][pair_key]["All3Data"][ckey]["eff_pawan_c2R"])
        if sr is not None:
            e_pm, [dm, dem, mmf] = slice_1d_range(e_pm, [dm, dem, mm.astype(float)], *sr)
            mm = mmf > 0.5
            e_pd, [dd, ded, mdf] = slice_1d_range(e_pd, [dd, ded, md.astype(float)], *sr)
            md = mdf > 0.5
        d_lo, d_hi, d_maj, d_min = choose_delta_ylim([(dm, mm), (dd, md)], min_half_span=0.5, max_half_span=5.0)
        pad20 = 0.12 * (d_hi - d_lo)
        configure_ratio_axis(ax_b, xlim=xlim, xlabel=xl, ylabel=r"$\Delta\varepsilon$ (Jul8$-$May30) [pp]", ylim=(d_lo - pad20, d_hi + pad20), ref_value=0.0, major_ytick=d_maj, minor_ytick=d_min)
        plot_binned_errorbar(ax_b, e_pm, dm, dem, mask=mm, label="MC", color=COLORS["mc"], marker="o", markersize=2.8, elinewidth=0.8)
        plot_binned_errorbar(ax_b, e_pd, dd, ded, mask=md, label="Data", color=COLORS["data"], marker="s", markersize=2.7, elinewidth=0.8)
        combine_legend(ax_b, loc="lower right", ncol=2, fontsize=7.2)

    # Column 2: Track Pz0 May30 vs Jul8
    ax_t2, ax_b2 = cols[2]
    mom = summary["momentum"]
    edges_p = to_float_array(mom[v_b]["10004X"]["track_pz0"]["edges"])
    c_log = np.sqrt(edges_p[:-1] * edges_p[1:])
    mc_b = to_float_array(mom[v_b]["10004X"]["track_pz0"]["norm"])
    mc_a = to_float_array(mom[v_a]["10004X"]["track_pz0"]["norm"])
    da_b = to_float_array(mom[v_b]["All3Data"]["track_pz0"]["norm"])
    da_a = to_float_array(mom[v_a]["All3Data"]["track_pz0"]["norm"])

    ax_t2.step(edges_p, np.concatenate([mc_b, mc_b[-1:]]), where="post", color=COLORS["mc"], lw=1.3, label="Jul8 MC")
    ax_t2.step(edges_p, np.concatenate([mc_a, mc_a[-1:]]), where="post", color=COLORS["mc_base"], lw=1.0, ls="--", label="May30 MC")
    ax_t2.step(edges_p, np.concatenate([da_b, da_b[-1:]]), where="post", color=COLORS["data"], lw=1.3, label="Jul8 Data")
    ax_t2.step(edges_p, np.concatenate([da_a, da_a[-1:]]), where="post", color=COLORS["data_base"], lw=1.0, ls="--", label="May30 Data")
    ax_t2.set_xscale("log")
    ax_t2.set_xlim(1.0, 2.0e7)
    ax_t2.set_ylim(0.0, max(0.08, float(np.nanmax(np.concatenate([mc_b, da_b]))) * 1.18))
    ax_t2.set_ylabel("Normalized Tracks / Bin")
    ax_t2.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax_t2.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.6, alpha=0.6)
    format_faser_header(ax_t2, subtitle="Track Pz0", right_text="Slide 16", fontsize=8.8, right_fontsize=7.6)
    combine_legend(ax_t2, loc="upper left", fontsize=7.4)

    c_ma = to_int_array(mom[v_a]["10004X"]["track_pz0"]["counts"])
    c_da = to_int_array(mom[v_a]["All3Data"]["track_pz0"]["counts"])
    m_mc = (mc_a > 0) & (c_ma >= 25)
    m_da = (da_a > 0) & (c_da >= 25)
    ax_b2.set_xscale("log")
    ax_b2.set_xlim(1.0, 2.0e7)
    ax_b2.set_ylim(0.96, 1.04)
    ax_b2.set_xlabel(r"Reconstructed Track $p_{z,0}$ [MeV]")
    ax_b2.set_ylabel("Jul8 / May30")
    ax_b2.yaxis.set_major_locator(MultipleLocator(0.02))
    ax_b2.axhline(1.0, color=COLORS["ref_line"], ls="--", lw=0.9)
    ax_b2.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax_b2.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.55, alpha=0.55)
    ax_b2.plot(c_log[m_mc], mc_b[m_mc] / mc_a[m_mc], "o-", color=COLORS["mc"], ms=2.4, lw=0.9, label="MC")
    ax_b2.plot(c_log[m_da], da_b[m_da] / da_a[m_da], "s-", color=COLORS["data"], ms=2.3, lw=0.9, label="Data")
    combine_legend(ax_b2, loc="lower left", ncol=2, fontsize=7.2)

    pd_g = summary["paired_global_deltas"][pair_key]
    da_g = pd_g["All3Data"]["global"]["eff_pawan_c2R"]
    mc_g = pd_g["10004X"]["global"]["eff_pawan_c2R"]
    extra = (
        f"Global Δε (Jul8-May30): Data {da_g['delta']*100:+.3f}±{da_g['delta_err']*100:.3f} pp, "
        f"MC {mc_g['delta']*100:+.3f}±{mc_g['delta_err']*100:.3f} pp"
    )
    add_figure_banners(
        fig,
        title="FIG20: Baseline Validation — May30 (June 9 Slides) vs Jul8 (3-Station checkIsLongTrack Fix)",
        scope_line=get_scope_line(summary, "baseline_jul8", extra=extra),
        footer_line=f"{METRIC_FOOTNOTE} | Exact matched 322 pairs | Verifies Pawan baseline reproduction & Jul8 ntupler fix shift",
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, "fig20_baseline_may30_vs_jul8_validation")
    return {
        "fig_id": "fig20",
        "suite": "C",
        "title": "Baseline Validation: May30 vs Jul8",
        "pawan_ref": "Validation",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": {
            "global_delta_data_pp": da_g["delta"] * 100.0,
            "global_delta_mc_pp": mc_g["delta"] * 100.0,
        },
    }


# ---------------------------------------------------------------------------
# Registry, parallel worker dispatch, and multi-page PDF assembly
# ---------------------------------------------------------------------------
FIGURE_BUILDERS: List[Tuple[str, str, Callable[..., Dict[str, Any]]]] = [
    ("fig01", "A", make_fig01),
    ("fig02", "A", make_fig02),
    ("fig03", "A", make_fig03),
    ("fig04", "A", make_fig04),
    ("fig05", "A", make_fig05),
    ("fig06", "A", make_fig06),
    ("fig07", "A", make_fig07),
    ("fig08", "A", make_fig08),
    ("fig09", "B", make_fig09),
    ("fig10", "B", make_fig10),
    ("fig11", "B", make_fig11),
    ("fig12", "B", make_fig12),
    ("fig13", "B", make_fig13),
    ("fig14", "B", make_fig14),
    ("fig15", "B", make_fig15),
    ("fig16", "C", make_fig16),
    ("fig17", "C", make_fig17),
    ("fig18", "C", make_fig18),
    ("fig19", "C", make_fig19),
    ("fig20", "C", make_fig20),
]

_BUILDER_MAP: Dict[str, Callable[..., Dict[str, Any]]] = {fid: fn for fid, _, fn in FIGURE_BUILDERS}

# Worker process cache to avoid re-deserializing summary JSON for every task
_WORKER_CACHE: Dict[str, Any] = {}


def _init_worker(summary_path: str, signal_path: Optional[str]) -> None:
    _WORKER_CACHE["summary"] = json.loads(Path(summary_path).read_text(encoding="utf-8"))
    if signal_path and Path(signal_path).is_file():
        _WORKER_CACHE["signal"] = json.loads(Path(signal_path).read_text(encoding="utf-8"))
    else:
        _WORKER_CACHE["signal"] = None


def _run_worker_task(task: Tuple[str, str, str]) -> Dict[str, Any]:
    fig_id, out_dir_str, primary_version = task
    fn = _BUILDER_MAP[fig_id]
    return fn(
        summary=_WORKER_CACHE["summary"],
        signal_data=_WORKER_CACHE["signal"],
        out_dir=Path(out_dir_str),
        version=primary_version,
    )


def merge_pdfs(pdf_paths: Sequence[Path], output_pdf: Path) -> str:
    """Assemble ordered per-figure vector PDFs into one unified multi-page PDF."""
    existing = [str(p) for p in pdf_paths if p.is_file()]
    if not existing:
        raise RuntimeError("No per-figure PDF files found to assemble.")
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    tmp_pdf = output_pdf.with_suffix(".pdf.tmp")

    pdfunite = shutil.which("pdfunite")
    if pdfunite:
        cmd = [pdfunite, *existing, str(tmp_pdf)]
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if res.returncode == 0 and tmp_pdf.is_file() and tmp_pdf.stat().st_size > 0:
            tmp_pdf.replace(output_pdf)
            return "pdfunite"

    gs = shutil.which("gs")
    if gs:
        cmd = [
            gs, "-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite",
            f"-sOutputFile={tmp_pdf}", *existing,
        ]
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if res.returncode == 0 and tmp_pdf.is_file() and tmp_pdf.stat().st_size > 0:
            tmp_pdf.replace(output_pdf)
            return "ghostscript"

    raise RuntimeError("Neither pdfunite nor gs succeeded in merging vector PDFs.")


def generate_all_figures(
    summary_path: Path,
    signal_path: Optional[Path],
    out_dir: Path,
    mode: str = "full",
    selected_figs: Optional[Sequence[str]] = None,
    workers: int = 0,
) -> Dict[str, Any]:
    """Generate selected Stage 08 overlay figures, assemble report PDF, and write figure_index.json."""
    t0 = time.monotonic()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    signal_data = json.loads(signal_path.read_text(encoding="utf-8")) if (signal_path and signal_path.is_file()) else None

    versions = set(summary.get("versions", []))
    if mode == "full":
        if "latest_georot" not in versions:
            raise ValueError(
                f"--mode full requires 'latest_georot' in summary['versions'], found {sorted(versions)}. "
                "Use --mode baseline for baseline-only summaries."
            )
        primary_version = "latest_georot"
        allowed_suites = {"A", "B", "C"}
    elif mode == "baseline":
        primary_version = "baseline_jul8" if "baseline_jul8" in versions else sorted(versions)[0]
        allowed_suites = {"A", "C"}
    else:
        raise ValueError(f"Unsupported mode: {mode!r} (expected 'full' or 'baseline')")

    selected_set = set(selected_figs) if selected_figs else None
    to_run: List[str] = []
    skipped: List[str] = []
    for fid, suite, _ in FIGURE_BUILDERS:
        if suite not in allowed_suites:
            skipped.append(fid)
            continue
        if selected_set is not None and fid not in selected_set:
            continue
        to_run.append(fid)

    out_dir.mkdir(parents=True, exist_ok=True)
    n_workers = workers if workers > 0 else min(8, os.cpu_count() or 1, max(1, len(to_run)))

    results: List[Dict[str, Any]] = []
    if n_workers > 1 and len(to_run) > 1:
        tasks = [(fid, str(out_dir), primary_version) for fid in to_run]
        with ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_init_worker,
            initargs=(str(summary_path), str(signal_path) if signal_path else None),
        ) as ex:
            results = list(ex.map(_run_worker_task, tasks))
    else:
        for fid in to_run:
            results.append(_BUILDER_MAP[fid](summary=summary, signal_data=signal_data, out_dir=out_dir, version=primary_version))

    # Assemble unified multi-page PDF in canonical figure order
    report_pdf = out_dir / "stage08_overlay_2track_report.pdf"
    pdf_paths = [Path(r["pdf"]) for r in results]
    merge_tool = merge_pdfs(pdf_paths, report_pdf)

    elapsed = time.monotonic() - t0
    index_payload = {
        "schema_version": "1.0",
        "mode": mode,
        "primary_version": primary_version,
        "summary_input": str(summary_path),
        "signal_input": str(signal_path) if signal_path else None,
        "report_pdf": str(report_pdf),
        "merge_tool": merge_tool,
        "page_count": len(results),
        "skipped_in_mode": skipped,
        "elapsed_seconds": round(elapsed, 3),
        "figures": results,
    }
    index_path = out_dir / "figure_index.json"
    index_path.write_text(json.dumps(index_payload, indent=2), encoding="utf-8")
    return index_payload


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Stage 08 Two-Track Overlay figures (PNG + vector PDF + report PDF).")
    parser.add_argument("--summary", type=Path, required=True, help="Path to overlay_summary.json")
    parser.add_argument("--signal", type=Path, default=None, help="Path to signal_110038.json (optional)")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory for png/, pdf/, report PDF, figure_index.json")
    parser.add_argument("--mode", choices=["full", "baseline"], default="full", help="full (20 figs) or baseline (13 figs Suite A+C preview)")
    parser.add_argument("--figs", type=str, default=None, help="Comma-separated subset of figure IDs (e.g. fig01,fig09,fig17)")
    parser.add_argument("--workers", type=int, default=0, help="Parallel worker processes (0 = auto, 1 = sequential)")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    selected = [s.strip() for s in args.figs.split(",") if s.strip()] if args.figs else None
    payload = generate_all_figures(
        summary_path=args.summary,
        signal_path=args.signal,
        out_dir=args.out_dir,
        mode=args.mode,
        selected_figs=selected,
        workers=args.workers,
    )
    print(
        f"[overlay_figures] mode={payload['mode']} primary={payload['primary_version']} "
        f"generated={payload['page_count']} pages in {payload['elapsed_seconds']:.2f}s "
        f"-> {payload['report_pdf']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
