#!/usr/bin/env python3
"""Stage 08 Tier-2 Attribution Control Figure Suite (fig21 & fig22).

Visualizes the 4-way paired attribution ladder across the 36 matched quadruples
(18 shards 008749 Data + 18 shards 10004X MC):
  1. Baseline-Jul8 (Unrotated, 811b81b8)
  2. Pre-fix GeoRot (7fd40df8, legacySwappedCovLoc)
  3. Latest GeoRot (218813c2, default covPhi = 0.1089 rad^2)
  4. Latest GeoRot + covPhi = 8.41e-6 rad^2 (1-line CKF2Config fix)

Produces:
- <out_dir>/png/fig21_tier2_covphi_crater_attribution.png  (300 DPI)
- <out_dir>/pdf/fig21_tier2_covphi_crater_attribution.pdf  (vector PDF)
- <out_dir>/png/fig22_tier2_covphi_regimes_and_loss.png    (300 DPI)
- <out_dir>/pdf/fig22_tier2_covphi_regimes_and_loss.pdf    (vector PDF)
- Re-assembles <out_dir>/stage08_overlay_2track_report.pdf (22 pages if fig01..fig20 exist)
- Updates <out_dir>/figure_index.json (and writes <out_dir>/tier2_figure_index.json).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from analysis.overlay_figures import (
    _bin_eff_stats,
    _paired_delta_pp,
    _ratio_stats,
    merge_pdfs,
    to_float_array,
    to_int_array,
)
from analysis.overlay_plot_style import (
    COLORS,
    add_figure_banners,
    apply_style,
    choose_delta_ylim,
    choose_eff_ylim,
    choose_ratio_ylim,
    combine_legend,
    configure_efficiency_axis,
    configure_ratio_axis,
    format_count_compact,
    format_faser_header,
    plot_binned_errorbar,
    save_figure_outputs,
)


# ---------------------------------------------------------------------------
# 4-Arm visual styling (colorblind-safe, shape-distinct)
# ---------------------------------------------------------------------------
ARM_ORDER: List[str] = ["baseline_jul8", "t2b_7fd40df8", "latest_218813c2", "t2a_covphi"]

ARM_STYLE: Dict[str, Dict[str, Any]] = {
    "baseline_jul8": {
        "label": "Baseline-Jul8 (Unrotated)",
        "short": "Baseline-Jul8",
        "color": "#222222",
        "marker": "^",
        "mfc": "white",
        "zorder": 4,
    },
    "t2b_7fd40df8": {
        "label": "Pre-fix GeoRot (7fd40df8)",
        "short": "Pre-fix 7fd40df8",
        "color": "#D97706",
        "marker": "s",
        "mfc": "white",
        "zorder": 5,
    },
    "latest_218813c2": {
        "label": "Latest GeoRot (218813c2, covPhi=0.1089)",
        "short": "Latest 218813c2",
        "color": "#D62728",
        "marker": "o",
        "mfc": "#D62728",
        "zorder": 6,
    },
    "t2a_covphi": {
        "label": "Latest GeoRot + covPhi=8.41e-6",
        "short": "Latest + covPhi=8.41e-6",
        "color": "#008B64",
        "marker": "D",
        "mfc": "#008B64",
        "zorder": 7,
    },
}

DELTA_PAIRS: List[Tuple[str, str]] = [
    ("t2b_7fd40df8", "t2b_7fd40df8_minus_baseline_jul8"),
    ("latest_218813c2", "latest_218813c2_minus_baseline_jul8"),
    ("t2a_covphi", "t2a_covphi_minus_baseline_jul8"),
]

TIER2_FOOTNOTE = (
    "Tier-2 36-shard matched quadruples (18 008749 + 18 10004X) | "
    "Metric: eff_pawan_c2R (lt==2, chi2/nDoF<=25, R0,R1<95 mm) | Frame: FASER | Cond: OFLCOND-FASER-05"
)


def _tier2_scope_line(t2: Dict[str, Any]) -> str:
    gs = t2["global_summary"]["baseline_jul8"]
    n_da = int(gs["008749"]["n_events"])
    n_mc = int(gs["10004X"]["n_events"])
    prov = t2.get("provenance", {})
    nq = int((prov.get("matched_quadruples") or {}).get("total", 36))
    return (
        f"4-Way Paired Attribution ({nq} matched quadruples)  |  "
        f"008749 Data: {format_count_compact(n_da)} evts  |  "
        f"10004X MC: {format_count_compact(n_mc)} evts  |  "
        f"CKF2Config.py:67 covPhi = 0.1089 -> 8.41e-6 rad^2"
    )


def _draw_4way_eff_and_delta_panel(
    ax_top: Axes,
    ax_bot: Axes,
    t2: Dict[str, Any],
    dataset: str,
    curve_key: str,
    xlabel: str,
    xlim: Tuple[float, float],
    panel_subtitle: str,
    right_tag: str,
    legend_loc: str = "lower left",
    shade_crater: bool = False,
) -> Dict[str, Any]:
    """Render one 4-arm efficiency + 3-arm paired delta vs Baseline-Jul8 stacked subplot."""
    eff_series = []
    for arm in ARM_ORDER:
        c = t2["curves"][arm][dataset][curve_key]["eff_pawan_c2R"]
        e, v, ye, m = _bin_eff_stats(c)
        eff_series.append((arm, e, v, ye, m))

    ylim = choose_eff_ylim([(v, m) for _, _, v, _, m in eff_series], default=(0.80, 1.00), min_floor=0.50)
    configure_efficiency_axis(ax_top, xlim=xlim, ylim=ylim, has_right_twin=False)
    format_faser_header(ax_top, subtitle=panel_subtitle, right_text=right_tag, fontsize=9.3, right_fontsize=7.9)

    if shade_crater:
        ax_top.axvspan(12.5, 18.5, color="#f8d7da", alpha=0.32, zorder=0)
        ax_bot.axvspan(12.5, 18.5, color="#f8d7da", alpha=0.32, zorder=0)
        ax_top.text(
            15.5, ylim[1] - 0.012, "Crater\n[12.5, 18.5]",
            ha="center", va="top", fontsize=7.0, color="#8b1e24",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="#e2b4b8", alpha=0.85),
            zorder=9,
        )

    n_bins = len(eff_series[0][1]) - 1
    ms = 3.7 if n_bins <= 60 else 2.5
    ew = 0.85 if n_bins <= 60 else 0.72
    for arm, e, v, ye, m in eff_series:
        st = ARM_STYLE[arm]
        plot_binned_errorbar(
            ax_top, e, v, ye, mask=m,
            label=st["label"], color=st["color"], marker=st["marker"],
            markerfacecolor=st["mfc"], markersize=ms, elinewidth=ew,
            zorder=st["zorder"],
        )
    combine_legend(ax_top, loc=legend_loc, fontsize=7.1)

    # Bottom panel: 3 rotated arms minus Baseline-Jul8 [pp]
    delta_list = []
    err_list = []
    delta_data = []
    for arm, pkey in DELTA_PAIRS:
        pc = t2["paired_curves"][pkey][dataset][curve_key]["eff_pawan_c2R"]
        ep, dp, dep, mp = _paired_delta_pp(pc)
        delta_list.append((dp, mp))
        err_list.append(dep)
        delta_data.append((arm, ep, dp, dep, mp))

    d_lo, d_hi, d_maj, d_min = choose_delta_ylim(delta_list, min_half_span=2.0, max_half_span=45.0, err_arrays=err_list)
    pad = 0.10 * (d_hi - d_lo)
    configure_ratio_axis(
        ax_bot, xlim=xlim, xlabel=xlabel,
        ylabel=r"$\Delta\varepsilon$ vs Jul8 [pp]",
        ylim=(d_lo - pad, d_hi + pad), ref_value=0.0, major_ytick=d_maj, minor_ytick=d_min,
    )
    for arm, ep, dp, dep, mp in delta_data:
        st = ARM_STYLE[arm]
        plot_binned_errorbar(
            ax_bot, ep, dp, dep, mask=mp,
            label=st["short"], color=st["color"], marker=st["marker"],
            markerfacecolor=st["mfc"], markersize=max(2.3, ms - 0.4), elinewidth=ew,
            zorder=st["zorder"],
        )
    combine_legend(ax_bot, loc="lower left" if not shade_crater else "upper left", ncol=3, fontsize=6.6)

    # Extract min/max delta of Latest vs T2-A for manifest stats
    def _minmax(arm_target: str):
        for arm, _, dp, _, mp in delta_data:
            if arm == arm_target and np.any(mp):
                return [float(np.nanmin(dp[mp])), float(np.nanmax(dp[mp]))]
        return None

    return {
        "dataset": dataset,
        "curve_key": curve_key,
        "ylim_eff": list(ylim),
        "latest_218813c2_delta_pp": _minmax("latest_218813c2"),
        "t2a_covphi_delta_pp": _minmax("t2a_covphi"),
    }


def make_fig21(t2: Dict[str, Any], out_dir: Path) -> Dict[str, Any]:
    """FIG21: 4-way attribution on |Delta Y0| [0, 20] mm (crater erased) & Delta R1 [0, 20] mm (Data + MC)."""
    apply_style()
    fig = plt.figure(figsize=(15.0, 9.2))
    outer = GridSpec(2, 2, figure=fig, left=0.062, right=0.975, top=0.885, bottom=0.085, wspace=0.19, hspace=0.28)

    specs = [
        (0, 0, "008749", "abs_dY0_zoom", r"Station 1 Vertical Separation ($|\Delta Y_0|$) [mm]", (0.0, 20.0), "008749 Data: |Delta Y0| Crater", "0.5 mm bins", "lower left", True),
        (0, 1, "10004X", "abs_dY0_zoom", r"Station 1 Vertical Separation ($|\Delta Y_0|$) [mm]", (0.0, 20.0), "10004X MC: |Delta Y0| Crater", "0.5 mm bins", "lower left", True),
        (1, 0, "008749", "dR1_zoom", r"Station 3 Radial Separation ($\Delta R_1$) [mm]", (0.0, 20.0), "008749 Data: Small Delta R1", "0.5 mm bins", "lower right", False),
        (1, 1, "10004X", "dR1_zoom", r"Station 3 Radial Separation ($\Delta R_1$) [mm]", (0.0, 20.0), "10004X MC: Small Delta R1", "0.5 mm bins", "lower right", False),
    ]
    panel_stats = {}
    for row, col, ds, ck, xl, xlim, sub, rtag, lloc, shade in specs:
        inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[row, col], height_ratios=(2.45, 1.15), hspace=0.10)
        ax_t = fig.add_subplot(inner[0, 0])
        ax_b = fig.add_subplot(inner[1, 0], sharex=ax_t)
        plt.setp(ax_t.get_xticklabels(), visible=False)
        panel_stats[f"{ds}_{ck}"] = _draw_4way_eff_and_delta_panel(
            ax_t, ax_b, t2, ds, ck, xl, xlim, sub, rtag, legend_loc=lloc, shade_crater=shade
        )

    add_figure_banners(
        fig,
        title="FIG21: Tier-2 Attribution — covPhi=8.41e-6 Eliminates |Delta Y0| Crater & Restores Small-Delta R1",
        scope_line=_tier2_scope_line(t2),
        footer_line=f"{TIER2_FOOTNOTE} | Paired McNemar errors vs Baseline-Jul8",
        title_y=0.980,
        scope_y=0.946,
        footer_y=0.010,
        footer_fontsize=7.3,
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, "fig21_tier2_covphi_crater_attribution")
    return {
        "fig_id": "fig21",
        "suite": "Tier2",
        "title": "Tier-2 Attribution: |Delta Y0| Crater Elimination & Small-Delta R1 Recovery",
        "pawan_ref": "Tier-2 Control (Slides 7, 9)",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": panel_stats,
    }


def make_fig22(t2: Dict[str, Any], out_dir: Path) -> Dict[str, Any]:
    """FIG22: 4-way Delta R1 [0,100] & MC/Data closure, |Delta X0| [0,100], regime summary & stage-of-loss."""
    apply_style()
    fig = plt.figure(figsize=(15.4, 9.2))
    outer = GridSpec(2, 2, figure=fig, left=0.068, right=0.975, top=0.885, bottom=0.085, wspace=0.22, hspace=0.30)

    # ---- Top-Left: Delta R1 [0, 100] mm — Top: paired Delta (10004X & 008749); Bottom: 4-arm MC/Data ratio ----
    inner_tl = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0, 0], height_ratios=(1.8, 1.3), hspace=0.08)
    ax_tl_t = fig.add_subplot(inner_tl[0, 0])
    ax_tl_b = fig.add_subplot(inner_tl[1, 0], sharex=ax_tl_t)
    plt.setp(ax_tl_t.get_xticklabels(), visible=False)

    # Rebin dR1_wide (200 x 0.5 mm -> 50 x 2.0 mm) visually? No, keep native 0.5 mm paired points cleanly!
    tl_deltas, tl_errs = [], []
    tl_plot_items = []
    for ds, mk_override, ls_lbl in [("008749", "o", "Data"), ("10004X", "^", "MC")]:
        for arm, pkey in DELTA_PAIRS:
            pc = t2["paired_curves"][pkey][ds]["dR1_wide"]["eff_pawan_c2R"]
            ep, dp, dep, mp = _paired_delta_pp(pc)
            # Filter tiny-stat bins (< 30) only for clean scale estimation
            good = mp & (to_int_array(pc["n_den"]) >= 30)
            tl_deltas.append((dp, good))
            tl_errs.append(dep)
            tl_plot_items.append((ds, arm, ep, dp, dep, good))

    dlo, dhi, dmaj, dmin = choose_delta_ylim(tl_deltas, min_half_span=3.0, max_half_span=15.0, err_arrays=tl_errs)
    pad = 0.08 * (dhi - dlo)
    configure_ratio_axis(
        ax_tl_t, xlim=(0.0, 100.0), xlabel="", ylabel=r"$\Delta\varepsilon$ vs Jul8 [pp]",
        ylim=(dlo - pad, dhi + pad), ref_value=0.0, major_ytick=dmaj, minor_ytick=dmin,
    )
    format_faser_header(ax_tl_t, subtitle="Delta R1 Plateau Recovery & MC/Data Closure", right_text="[0, 100] mm", fontsize=9.2, right_fontsize=7.8)
    for ds, arm, ep, dp, dep, good in tl_plot_items:
        st = ARM_STYLE[arm]
        cen = 0.5 * (ep[:-1] + ep[1:])
        if ds == "008749":
            ax_tl_t.plot(cen[good], dp[good], "-", color=st["color"], lw=1.35, alpha=0.95, label=f"Data: {st['short']}", zorder=st["zorder"])
        else:
            ax_tl_t.plot(cen[good], dp[good], "--", color=st["color"], lw=1.15, alpha=0.85, label=f"MC: {st['short']}", zorder=st["zorder"] - 1)
    combine_legend(ax_tl_t, loc="lower right", ncol=2, fontsize=6.6)

    # Bottom of Top-Left: MC/Data ratio for all 4 arms
    r_items = []
    for arm in ARM_ORDER:
        rc = t2["mc_data_ratios"][arm]["dR1_wide"]["eff_pawan_c2R"]
        er, rv, yer, mr = _ratio_stats(rc)
        r_items.append((arm, er, rv, yer, mr))
    rlo, rhi, rmaj, rmin = choose_ratio_ylim([(rv, yer, mr) for _, _, rv, yer, mr in r_items], default=(0.95, 1.05), max_half_span=0.10)
    configure_ratio_axis(
        ax_tl_b, xlim=(0.0, 100.0), xlabel=r"Track Separation ($\Delta R_1$) [mm]", ylabel="MC / Data",
        ylim=(rlo, rhi), ref_value=1.0, major_ytick=rmaj, minor_ytick=rmin,
    )
    for arm, er, rv, yer, mr in r_items:
        st = ARM_STYLE[arm]
        cen = 0.5 * (er[:-1] + er[1:])
        ax_tl_b.plot(
            cen[mr], rv[mr], color=st["color"],
            ls="--" if arm == "baseline_jul8" else "-",
            lw=1.3 if arm in ("t2a_covphi", "latest_218813c2") else 1.05,
            alpha=0.92, label=st["short"], zorder=st["zorder"],
        )
    combine_legend(ax_tl_b, loc="lower right", ncol=2, fontsize=6.6)

    # ---- Top-Right: |Delta X0| [0, 100] mm — Top: 008749 Data paired delta; Bottom: 10004X MC paired delta ----
    inner_tr = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0, 1], height_ratios=(1.0, 1.0), hspace=0.14)
    ax_tr_t = fig.add_subplot(inner_tr[0, 0])
    ax_tr_b = fig.add_subplot(inner_tr[1, 0], sharex=ax_tr_t)
    plt.setp(ax_tr_t.get_xticklabels(), visible=False)

    for ax_x, ds, ds_lbl, show_xl in [
        (ax_tr_t, "008749", "008749 Data: Delta vs Jul8", False),
        (ax_tr_b, "10004X", "10004X MC: Delta vs Jul8", True),
    ]:
        xd, xe, xitems = [], [], []
        for arm, pkey in DELTA_PAIRS:
            pc = t2["paired_curves"][pkey][ds]["abs_dX0_wide"]["eff_pawan_c2R"]
            ep, dp, dep, mp = _paired_delta_pp(pc)
            good = mp & (to_int_array(pc["n_den"]) >= 30)
            xd.append((dp, good))
            xe.append(dep)
            xitems.append((arm, ep, dp, dep, good))
        xlo, xhi, xmaj, xmin_t = choose_delta_ylim(xd, min_half_span=2.0, max_half_span=10.0, err_arrays=xe)
        xpad = 0.10 * (xhi - xlo)
        configure_ratio_axis(
            ax_x, xlim=(0.0, 100.0),
            xlabel=r"Station 1 Horizontal Separation ($|\Delta X_0|$) [mm]" if show_xl else "",
            ylabel=r"$\Delta\varepsilon$ [pp]",
            ylim=(xlo - xpad, xhi + xpad), ref_value=0.0, major_ytick=xmaj, minor_ytick=xmin_t,
        )
        for arm, ep, dp, dep, good in xitems:
            st = ARM_STYLE[arm]
            plot_binned_errorbar(
                ax_x, ep, dp, dep, mask=good, label=st["short"],
                color=st["color"], marker=st["marker"], markerfacecolor=st["mfc"],
                markersize=2.6, elinewidth=0.72, zorder=st["zorder"],
            )
        if not show_xl:
            format_faser_header(ax_x, subtitle="Horizontal Separation |Delta X0| Shift", right_text="Top: Data | Bot: MC", fontsize=9.2, right_fontsize=7.8)
        ax_x.text(
            0.02, 0.88, ds_lbl, transform=ax_x.transAxes, fontsize=7.8, fontweight="bold", color="#222222",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#dddddd", alpha=0.88),
        )
        combine_legend(ax_x, loc="lower right" if show_xl else "upper right", ncol=3, fontsize=6.5)

    # ---- Bottom-Left: 7-Regime Paired Delta Summary (Data & MC, 3 GeoRot arms vs Baseline-Jul8) ----
    ax_bl = fig.add_subplot(outer[1, 0])
    regime_specs = [
        ("global", "Global (36 sh)"),
        ("crater_abs_dY0_12p5_18p5", "Crater |dY0|\n[12.5, 18.5]"),
        ("small_abs_dY0_lt2", "Small |dY0|\n< 2 mm"),
        ("large_abs_dY0_ge19p5", "Large |dY0|\n>= 19.5 mm"),
        ("small_dR1_lt5", "Small dR1\n< 5 mm"),
        ("plateau_dR1_30_70", "Plateau dR1\n[30, 70] mm"),
        ("large_abs_dX0_70_100", "Large |dX0|\n[70, 100] mm"),
    ]
    y_idx = np.arange(len(regime_specs))[::-1].astype(float)
    offsets = [
        ("008749", "t2b_7fd40df8", "t2b_7fd40df8_minus_baseline_jul8", +0.26, "s", "white", "Data: 7fd40df8"),
        ("008749", "latest_218813c2", "latest_218813c2_minus_baseline_jul8", +0.10, "o", None, "Data: 218813c2"),
        ("008749", "t2a_covphi", "t2a_covphi_minus_baseline_jul8", -0.06, "D", None, "Data: +covPhi=8.41e-6"),
        ("10004X", "latest_218813c2", "latest_218813c2_minus_baseline_jul8", -0.20, "o", "white", "MC: 218813c2"),
        ("10004X", "t2a_covphi", "t2a_covphi_minus_baseline_jul8", -0.34, "D", "white", "MC: +covPhi=8.41e-6"),
    ]
    all_dx = []
    regime_stats: Dict[str, Any] = {}
    for ds, arm, pkey, dy, mk, mfc_ov, lbl in offsets:
        xs, xerrs = [], []
        for rkey, _ in regime_specs:
            drec = t2["paired_deltas"][pkey][ds][rkey]["eff_pawan_c2R"]
            val = float(drec["delta"]) * 100.0
            err = float(drec["delta_err"]) * 100.0
            xs.append(val)
            xerrs.append(err)
            regime_stats[f"{ds}_{arm}_{rkey}_pp"] = round(val, 4)
        all_dx.extend(xs)
        col = ARM_STYLE[arm]["color"]
        ax_bl.errorbar(
            xs, y_idx + dy, xerr=xerrs, fmt=mk, color=col,
            markerfacecolor=col if mfc_ov is None else mfc_ov, markeredgecolor=col,
            markersize=5.2, elinewidth=1.0, capsize=2.0, label=lbl, zorder=5,
        )
    ax_bl.axvline(0.0, color=COLORS["ref_line"], ls="--", lw=1.0, zorder=1)
    x_min = min(-5.0, float(np.floor(min(all_dx) - 2.0)))
    x_max = max(8.0, float(np.ceil(max(all_dx) + 2.5)))
    ax_bl.set_xlim(x_min, x_max)
    ax_bl.set_ylim(-0.65, len(regime_specs) - 0.35)
    ax_bl.set_yticks(y_idx)
    ax_bl.set_yticklabels([lbl for _, lbl in regime_specs], fontsize=8.0)
    ax_bl.set_xlabel(r"Paired $\Delta\varepsilon$ (Arm $-$ Baseline-Jul8) [pp]")
    ax_bl.grid(axis="x", color=COLORS["grid"], ls=":", lw=0.65, alpha=0.65)
    ax_bl.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    format_faser_header(ax_bl, subtitle="Regime Paired-Delta Summary", right_text="Solid=Data | Open=MC", fontsize=9.2, right_fontsize=7.6)
    combine_legend(ax_bl, loc="lower left", ncol=1, fontsize=7.0)

    # ---- Bottom-Right: Global Stage-of-Loss & lt=0 Crash Elimination Across All 4 Arms ----
    ax_br = fig.add_subplot(outer[1, 1])
    loss_groups = [
        ("008749", "lt0", "Data\nlt==0"),
        ("008749", "ckf", "Data\nCKF Lost"),
        ("008749", "miss_seg", "Data\nMiss Seg"),
        ("10004X", "lt0", "MC\nlt==0"),
        ("10004X", "ckf", "MC\nCKF Lost"),
        ("10004X", "miss_seg", "MC\nMiss Seg"),
    ]

    def _loss_pct(arm: str, ds: str, kind: str) -> float:
        sl = t2["stage_of_loss"][arm][ds]["global"]
        den = max(1, int(sl["n_den"]))
        if kind == "lt0":
            num = int(sl["lt0"])
        elif kind == "ckf":
            num = int(sl["lt1_all3gb_lost_in_ckf"])
        elif kind == "miss_seg":
            num = int(sl["lt1_miss_seg_st1"]) + int(sl["lt1_miss_seg_st2"]) + int(sl["lt1_miss_seg_st3"]) + int(sl["lt1_miss_seg_multi"])
        else:
            num = 0
        return 100.0 * num / den

    x_g = np.arange(len(loss_groups), dtype=float)
    w = 0.19
    max_bar = 0.0
    loss_stats: Dict[str, Any] = {}
    for a_idx, arm in enumerate(ARM_ORDER):
        vals = [_loss_pct(arm, ds, kind) for ds, kind, _ in loss_groups]
        max_bar = max(max_bar, max(vals))
        st = ARM_STYLE[arm]
        offset = (a_idx - 1.5) * w
        bars = ax_br.bar(
            x_g + offset, vals, width=w * 0.92,
            color=st["color"], alpha=0.88 if arm != "baseline_jul8" else 0.72,
            edgecolor="#222222", linewidth=0.5, label=st["short"], zorder=4,
        )
        for (ds, kind, _), b, v in zip(loss_groups, bars, vals):
            loss_stats[f"{arm}_{ds}_{kind}_pct"] = round(v, 3)
            if v >= 0.01:
                ax_br.text(b.get_x() + b.get_width() / 2.0, v + 0.08, f"{v:.2f}", ha="center", va="bottom", fontsize=5.8, rotation=90, color="#222222")

    ax_br.set_xticks(x_g)
    ax_br.set_xticklabels([lbl for _, _, lbl in loss_groups], fontsize=8.2)
    ax_br.set_ylim(0.0, max_bar * 1.38)
    ax_br.set_ylabel("Global Event Fraction [%]")
    ax_br.axvline(2.5, color="#888888", ls="--", lw=0.85)
    ax_br.grid(axis="y", color=COLORS["grid"], ls=":", lw=0.65, alpha=0.65)
    ax_br.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    format_faser_header(ax_br, subtitle="Global Failure Funnel (lt=0 & CKF Eliminated)", right_text="36 matched shards", fontsize=9.1, right_fontsize=7.6)
    combine_legend(ax_br, loc="upper left", ncol=2, fontsize=7.0)

    add_figure_banners(
        fig,
        title="FIG22: Tier-2 Attribution — Plateau/Horizontal Recovery, Regime Summary & Failure Funnel",
        scope_line=_tier2_scope_line(t2),
        footer_line=f"{TIER2_FOOTNOTE} | covPhi=8.41e-6 slashes lt==0 & CKF loss across Data & MC",
        title_y=0.980,
        scope_y=0.946,
        footer_y=0.010,
        footer_fontsize=7.3,
    )
    png_path, pdf_path = save_figure_outputs(fig, out_dir, "fig22_tier2_covphi_regimes_and_loss")
    return {
        "fig_id": "fig22",
        "suite": "Tier2",
        "title": "Tier-2 Attribution: Plateau, Horizontal Separation, Regime Summary & Stage-of-Loss",
        "pawan_ref": "Tier-2 Control (Slides 7, 11, 15)",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "stats": {"regimes": regime_stats, "loss_funnel_pct": loss_stats},
    }


def generate_tier2_figures(
    tier2_summary_path: Path,
    out_dir: Path,
    append_report: bool = True,
) -> Dict[str, Any]:
    """Generate fig21 & fig22, optionally re-assemble unified 22-page PDF and update figure_index.json."""
    t0 = time.monotonic()
    t2 = json.loads(tier2_summary_path.read_text(encoding="utf-8"))
    out_dir.mkdir(parents=True, exist_ok=True)

    r21 = make_fig21(t2, out_dir)
    r22 = make_fig22(t2, out_dir)
    t2_figs = [r21, r22]

    # Always write standalone Tier-2 figure index
    t2_pdf = out_dir / "stage08_tier2_attribution_figures.pdf"
    merge_tool = merge_pdfs([Path(r21["pdf"]), Path(r22["pdf"])], t2_pdf)

    report_pdf = out_dir / "stage08_overlay_2track_report.pdf"
    index_path = out_dir / "figure_index.json"
    combined_pages = len(t2_figs)

    if append_report and index_path.is_file():
        existing_idx = json.loads(index_path.read_text(encoding="utf-8"))
        # Remove any stale fig21/fig22 entries before appending fresh ones
        base_figs = [f for f in existing_idx.get("figures", []) if f.get("fig_id") not in ("fig21", "fig22")]
        merged_figs = base_figs + t2_figs
        pdf_paths = [Path(f["pdf"]) for f in merged_figs]
        merge_tool = merge_pdfs(pdf_paths, report_pdf)
        combined_pages = len(merged_figs)

        existing_idx["figures"] = merged_figs
        existing_idx["page_count"] = combined_pages
        existing_idx["tier2_summary_input"] = str(tier2_summary_path)
        existing_idx["merge_tool"] = merge_tool
        index_path.write_text(json.dumps(existing_idx, indent=2), encoding="utf-8")

    elapsed = time.monotonic() - t0
    t2_payload = {
        "schema_version": "1.0",
        "tier2_summary_input": str(tier2_summary_path),
        "tier2_pdf": str(t2_pdf),
        "combined_report_pdf": str(report_pdf) if report_pdf.is_file() else None,
        "combined_page_count": combined_pages,
        "merge_tool": merge_tool,
        "elapsed_seconds": round(elapsed, 3),
        "figures": t2_figs,
    }
    (out_dir / "tier2_figure_index.json").write_text(json.dumps(t2_payload, indent=2), encoding="utf-8")
    return t2_payload


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Stage 08 Tier-2 Attribution Control figures (fig21 & fig22).")
    parser.add_argument("--tier2-summary", type=Path, required=True, help="Path to tier2_attribution_summary.json")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output figures directory (e.g. run_20261125_v1/figures)")
    parser.add_argument("--no-append-report", action="store_true", help="Do not rebuild unified 22-page stage08_overlay_2track_report.pdf")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    payload = generate_tier2_figures(
        tier2_summary_path=args.tier2_summary,
        out_dir=args.out_dir,
        append_report=(not args.no_append_report),
    )
    print(
        f"[overlay_tier2_figures] generated fig21-fig22 in {payload['elapsed_seconds']:.2f}s | "
        f"combined_pages={payload['combined_page_count']} -> {payload['combined_report_pdf'] or payload['tier2_pdf']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
