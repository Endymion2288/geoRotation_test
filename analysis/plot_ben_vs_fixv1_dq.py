#!/usr/bin/env python3
"""Plot FASER DQ-style distributions comparing geoRotation_calypso_latest_fix_v1
full reconstruction against Ben's baseline (Official r0021/r0022 PHYS +
El9CalypsoForWaveForm reduced charge) across all 5 periods (2022-2024 post-CaloNu).

Visual conventions mirror FASER_DQ_Analyis (notebooks/plotting_utils.py,
internal_note/finalize.py, and Electronic_2D_measurement_internal_note-6.pdf Figs. 5-6):
  - plt.style.use('default'), mplhep.histplot(..., histtype='errorbar', density=True)
  - Period color palette: 2022 red, 2023 blue, 2024_pre orange, 2024_caloNu purple, 2024_post green
  - Frameless legends, large clean axis labels, plus attached fix_v1 / Ben ratio strips.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import boost_histogram as bh
import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mplhep as hep
import numpy as np

DEFAULT_RUN_ROOT = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1"
)

PERIOD_ORDER = ["2022", "2023", "2024_pre", "2024_caloNu", "2024_post"]
PERIOD_META = {
    "2022":        {"label": "2022 (007961)",         "short": "2022",       "color": "#d62728"},  # red
    "2023":        {"label": "2023 (011523)",         "short": "2023",       "color": "#1f77b4"},  # blue
    "2024_pre":    {"label": "2024 pre (015387)",     "short": "2024 pre",   "color": "#ff7f0e"},  # orange
    "2024_caloNu": {"label": "2024 CaloNu (015849)",  "short": "2024 CaloNu","color": "#9467bd"},  # purple
    "2024_post":   {"label": "2024 post (016963)",    "short": "2024 post",  "color": "#2ca02c"},  # green
}


def to_bh(hdict: dict[str, Any]) -> bh.Histogram:
    edges = np.asarray(hdict["edges"], dtype=np.float64)
    counts = np.asarray(hdict["counts"], dtype=np.float64)
    h = bh.Histogram(bh.axis.Variable(edges), storage=bh.storage.Weight())
    h.view(flow=False).value[:] = counts
    h.view(flow=False).variance[:] = counts
    return h


def plot_dq_errorbar(
    ax,
    hdict: dict[str, Any],
    color: str,
    label: str | None = None,
    markersize: float = 2.6,
    density: bool = True,
):
    """FASER DQ-style errorbar histogram skipping empty (0-count) bins so log-y axes stay crisp."""
    edges = np.asarray(hdict["edges"], dtype=np.float64)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = edges[1:] - edges[:-1]
    counts = np.asarray(hdict["counts"], dtype=np.float64)
    tot = counts.sum()
    if tot <= 0:
        return
    m = counts > 0
    if density:
        norm = tot * widths
        y = counts[m] / norm[m]
        ye = np.sqrt(counts[m]) / norm[m]
    else:
        y = counts[m]
        ye = np.sqrt(counts[m])
    ax.errorbar(
        centers[m], y, xerr=0.5 * widths[m], yerr=ye,
        fmt="o", markersize=markersize, elinewidth=0.9, capsize=0,
        color=color, label=label,
    )


def density_ratio(h_num: dict[str, Any], h_den: dict[str, Any]):
    """Compute density-normalized ratio (fix_v1 / Ben) and statistical error per bin."""
    edges = np.asarray(h_num["edges"], dtype=np.float64)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = edges[1:] - edges[:-1]
    c1 = np.asarray(h_num["counts"], dtype=np.float64)
    c2 = np.asarray(h_den["counts"], dtype=np.float64)
    s1, s2 = c1.sum(), c2.sum()
    if s1 <= 0 or s2 <= 0:
        return centers, np.full_like(centers, np.nan), np.full_like(centers, np.nan)
    d1 = c1 / (s1 * widths)
    d2 = c2 / (s2 * widths)
    valid = (c1 > 0) & (c2 > 0)
    r = np.full_like(centers, np.nan)
    re = np.full_like(centers, np.nan)
    r[valid] = d1[valid] / d2[valid]
    re[valid] = r[valid] * np.sqrt(1.0 / c1[valid] + 1.0 / c2[valid])
    return centers, r, re


def plot_dq_overlay_with_ratio(
    ax_top,
    ax_bot,
    data: dict[str, Any],
    hkey: str,
    xlabel: str,
    title: str,
    logy: bool = False,
    cut_line: float | None = None,
    cut_label: str | None = None,
    ratio_ylim: tuple[float, float] = (0.75, 1.25),
    show_legend: bool = False,
):
    """Draw FASER DQ-style multi-period density overlay (fix_v1 solid vs Ben hollow/step) + ratio."""
    for p in PERIOD_ORDER:
        col = PERIOD_META[p]["color"]
        hm = data[p]["hists_fixv1"][hkey]
        hb = data[p]["hists_ben"][hkey]
        if sum(hm["counts"]) == 0 or sum(hb["counts"]) == 0:
            continue
        # Ben baseline: subtle step outline in same period color (alpha=0.65)
        hep.histplot(
            to_bh(hb), histtype="step", yerr=False, density=True, ax=ax_top,
            color=col, lw=1.15, ls="--", alpha=0.65,
        )
        # fix_v1: FASER DQ errorbar points (non-zero bins only)
        plot_dq_errorbar(ax_top, hm, color=col, label=PERIOD_META[p]["short"], markersize=2.6, density=True)
        # Ratio strip: fix_v1 / Ben
        xc, r, re = density_ratio(hm, hb)
        ok = np.isfinite(r) & (re < 0.35)
        ax_bot.errorbar(
            xc[ok], r[ok], yerr=re[ok], fmt=".", ms=2.6, elinewidth=0.75,
            color=col, alpha=0.85,
        )

    if cut_line is not None:
        ax_top.axvline(cut_line, color="#222222", ls=":", lw=1.2)
        ax_bot.axvline(cut_line, color="#222222", ls=":", lw=1.2)
        if cut_label:
            ax_top.text(
                cut_line, 0.94, f" {cut_label}", transform=ax_top.get_xaxis_transform(),
                fontsize=7.6, color="#222222", va="top", ha="left",
            )

    if logy:
        ax_top.set_yscale("log")
    else:
        ax_top.set_ylim(bottom=0.0)
        ymax = ax_top.get_ylim()[1]
        ax_top.set_ylim(0.0, ymax * 1.12)

    edges = data[PERIOD_ORDER[0]]["hists_fixv1"][hkey]["edges"]
    ax_top.set_xlim(edges[0], edges[-1])
    ax_bot.set_xlim(edges[0], edges[-1])
    ax_top.set_title(title, fontsize=10.5, pad=4)
    ax_top.set_ylabel("Prob. density", fontsize=9.0)
    ax_top.tick_params(labelbottom=False, labelsize=8.5, direction="in", top=True, right=True)

    ax_bot.axhline(1.0, color="#444444", ls="--", lw=0.9)
    ax_bot.set_ylim(0.72, 1.28)
    ax_bot.set_yticks([0.8, 1.0, 1.2])
    ax_bot.set_xlabel(xlabel, fontsize=9.8)
    ax_bot.set_ylabel("fix / Ben", fontsize=8.4)
    ax_bot.tick_params(labelsize=8.3, direction="in", top=True, right=True)
    ax_bot.grid(True, axis="y", color="#e6e6e6", lw=0.6)

    if show_legend:
        # Add style explanation proxies
        ax_top.plot([], [], color="#333333", marker="o", ms=3.0, ls="", label="fix_v1 (points)")
        ax_top.plot([], [], color="#333333", ls="--", lw=1.2, label="Ben (dashed)")
        ax_top.legend(frameon=False, fontsize=7.6, loc="upper right", ncol=1, labelspacing=0.22)


def make_grid_with_ratios(nrows: int, ncols: int, figsize: tuple[float, float], suptitle: str):
    fig = plt.figure(figsize=figsize)
    outer = gridspec.GridSpec(
        nrows, ncols, figure=fig, hspace=0.34, wspace=0.24,
        left=0.065, right=0.985, top=0.90, bottom=0.07,
    )
    pairs = []
    for r in range(nrows):
        for c in range(ncols):
            inner = gridspec.GridSpecFromSubplotSpec(
                2, 1, subplot_spec=outer[r, c], height_ratios=[3.2, 1.05], hspace=0.0,
            )
            ax_t = fig.add_subplot(inner[0, 0])
            ax_b = fig.add_subplot(inner[1, 0], sharex=ax_t)
            pairs.append((ax_t, ax_b))
    fig.suptitle(suptitle, fontsize=12.6, fontweight="bold", y=0.968)
    return fig, pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    args = parser.parse_args()

    plt.style.use("default")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "mathtext.default": "regular",
        "axes.linewidth": 0.9,
        "figure.dpi": 150,
    })

    data = json.loads((args.run_root / "analysis/ben_vs_fixv1_full.json").read_text())
    fig_dir = args.run_root / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    report_pdf = fig_dir / "ben_vs_fixv1_dq_report.pdf"
    index = []

    def save(fig, stem: str, desc: str, pdf: PdfPages):
        png = fig_dir / f"{stem}.png"
        pf = fig_dir / f"{stem}.pdf"
        fig.savefig(png, dpi=300)
        fig.savefig(pf)
        pdf.savefig(fig)
        plt.close(fig)
        index.append({"stem": stem, "description": desc, "png": str(png), "pdf": str(pf)})

    with PdfPages(report_pdf) as pdf:
        # ---- FIG 1: VetoNu0 & VetoNu1 Reduced Charge across ALL 5 periods (DQ selection, 2x5 with ratio) ----
        fig, pairs = make_grid_with_ratios(
            2, 5, figsize=(16.5, 7.6),
            suptitle=r"FASER DQ Style — VetoNu Reduced Charge: geoRotation_calypso_latest_fix_v1 vs Ben (2022–2024)",
        )
        for row_idx, (hkey, layer_lbl) in enumerate([
            ("VetoNu0_reduced_charge", r"FASER$\nu$ veto L0 reduced charge [pC]"),
            ("VetoNu1_reduced_charge", r"FASER$\nu$ veto L1 reduced charge [pC]"),
        ]):
            for col_idx, p in enumerate(PERIOD_ORDER):
                ax_t, ax_b = pairs[row_idx * 5 + col_idx]
                col = PERIOD_META[p]["color"]
                hm = data[p]["hists_fixv1"][hkey]
                hb = data[p]["hists_ben"][hkey]
                hep.histplot(
                    to_bh(hb), histtype="step", yerr=False, density=True, ax=ax_t,
                    color="#111111", lw=1.35, ls="--", label="Ben (El9WF)",
                )
                plot_dq_errorbar(ax_t, hm, color=col, label="fix_v1 (reco)", markersize=2.8, density=True)
                ax_t.axvline(30.0, color="#666666", ls=":", lw=1.1)
                ax_b.axvline(30.0, color="#666666", ls=":", lw=1.1)
                xc, r, re = density_ratio(hm, hb)
                ok = np.isfinite(r) & (re < 0.30)
                ax_b.errorbar(xc[ok], r[ok], yerr=re[ok], fmt="o", ms=2.4, elinewidth=0.8, color=col)
                ax_b.axhline(1.0, color="#444444", ls="--", lw=0.9)
                ax_t.set_xlim(0, 300)
                ax_b.set_xlim(0, 300)
                ax_t.set_ylim(0, ax_t.get_ylim()[1] * 1.18)
                ax_b.set_ylim(0.78, 1.22)
                ax_b.set_yticks([0.85, 1.0, 1.15])
                if row_idx == 0:
                    eb = data[p]["redq_ebye"]["all_vn0"]
                    ax_t.set_title(
                        f"{PERIOD_META[p]['label']}\n"
                        rf"$|\Delta Q_{{0}}|<0.01$: {100*eb['frac_lt_0p01']:.1f}%",
                        fontsize=9.4, pad=4,
                    )
                else:
                    eb = data[p]["redq_ebye"]["all_vn1"]
                    ax_t.set_title(rf"$|\Delta Q_{{1}}|<0.01$: {100*eb['frac_lt_0p01']:.1f}%", fontsize=8.8, pad=3)
                ax_t.tick_params(labelbottom=False, labelsize=8.2, direction="in", top=True, right=True)
                ax_b.tick_params(labelsize=8.0, direction="in", top=True, right=True)
                ax_b.grid(True, axis="y", color="#e6e6e6", lw=0.6)
                if col_idx == 0:
                    ax_t.set_ylabel("Prob. density [pC$^{-1}$]", fontsize=9.0)
                    ax_b.set_ylabel("fix / Ben", fontsize=8.5)
                if row_idx == 1:
                    ax_b.set_xlabel(layer_lbl, fontsize=8.8)
                else:
                    ax_b.set_xlabel("VetoNu L0 red. Q [pC]", fontsize=8.5)
                if col_idx == 4 and row_idx == 0:
                    ax_t.legend(frameon=False, fontsize=8.0, loc="upper right")
        save(
            fig, "dq_01_vetonu_reduced_charge_all_periods",
            "FASER DQ-style VetoNu0/1 reduced charge overlays and ratios across all 5 periods (fix_v1 vs Ben)", pdf,
        )

        # ---- FIG 2: nu_mu Pre-VetoNu Reduced Charge Zoom [0, 60] pC + Event-by-Event ΔQ_red ----
        fig = plt.figure(figsize=(14.5, 8.2))
        gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.20, left=0.07, right=0.98, top=0.90, bottom=0.08)
        fig.suptitle(
            r"Reduced Charge Closure Around the 30 pC $\nu_\mu$ Cut & Event-by-Event $\Delta Q_{\mathrm{red}} = Q_{\mathrm{fix\_v1}} - Q_{\mathrm{Ben}}$",
            fontsize=12.2, fontweight="bold", y=0.965,
        )

        # (a) VetoNu0 reduced charge zoom [0, 60] pC under nu_mu preselection (log y)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        for ax, hkey, ttl in [
            (ax1, "VetoNu0_reduced_charge_numu_zoom", r"(a) $\nu_\mu$-preselected VetoNu0 reduced charge ($<60$ pC)"),
            (ax2, "VetoNu1_reduced_charge_numu_zoom", r"(b) $\nu_\mu$-preselected VetoNu1 reduced charge ($<60$ pC)"),
        ]:
            for p in PERIOD_ORDER:
                col = PERIOD_META[p]["color"]
                hm = data[p]["hists_fixv1"][hkey]
                hb = data[p]["hists_ben"][hkey]
                if sum(hm["counts"]) == 0 or sum(hb["counts"]) == 0:
                    continue
                hep.histplot(to_bh(hb), histtype="step", yerr=False, density=True, ax=ax, color=col, ls="--", lw=1.2, alpha=0.65)
                plot_dq_errorbar(ax, hm, color=col, label=PERIOD_META[p]["short"], markersize=2.8, density=True)
            ax.axvline(30.0, color="#222222", ls=":", lw=1.4, label=r"$\nu_\mu$ cut ($<30$ pC)")
            ax.set_yscale("log")
            ax.set_xlim(0, 60)
            ax.set_ylim(bottom=2e-4)
            ax.set_xlabel("Reduced charge [pC]", fontsize=10.5)
            ax.set_ylabel("Probability density [pC$^{-1}$]", fontsize=10.0)
            ax.set_title(ttl, fontsize=10.4)
            ax.tick_params(direction="in", top=True, right=True)
            ax.grid(True, color="#ececec", lw=0.6)
        ax1.plot([], [], color="#333333", marker="o", ms=3.0, ls="", label="fix_v1 (points)")
        ax1.plot([], [], color="#333333", ls="--", lw=1.2, label="Ben (dashed)")
        ax1.legend(frameon=False, fontsize=8.0, loc="upper left", ncol=2)

        # (c) Event-by-event delta_redQ0 across all 5 periods (log-y)
        ax3 = fig.add_subplot(gs[1, 0])
        for p in PERIOD_ORDER:
            col = PERIOD_META[p]["color"]
            hd = data[p]["hists_ebye"]["delta_redQ0"]
            eb = data[p]["redq_ebye"]["all_vn0"]
            hep.histplot(
                to_bh(hd), histtype="step", yerr=False, density=True, ax=ax3, color=col, lw=1.5,
                label=f"{PERIOD_META[p]['short']} (|Δ|<0.01: {100*eb['frac_lt_0p01']:.1f}%, <1: {100*eb['frac_lt_1p0']:.1f}%)",
            )
        ax3.set_yscale("log")
        ax3.set_xlim(-20, 20)
        ax3.set_ylim(bottom=8e-6)
        ax3.set_xlabel(r"$\Delta Q_{\mathrm{red}}(\mathrm{VetoNu0}) = Q_{\mathrm{fix\_v1}} - Q_{\mathrm{Ben}}$ [pC]", fontsize=10.2)
        ax3.set_ylabel("Probability density [pC$^{-1}$]", fontsize=10.0)
        ax3.set_title(r"(c) Event-by-event VetoNu0 reduced-charge difference (DQ events)", fontsize=10.4)
        ax3.tick_params(direction="in", top=True, right=True)
        ax3.grid(True, color="#ececec", lw=0.6)
        ax3.legend(frameon=False, fontsize=7.8, loc="upper right")

        # (d) 2024 CaloNu breakdown: status == 0 vs status == 512 (+16 ns)
        ax4 = fig.add_subplot(gs[1, 1])
        cn = data["2024_caloNu"]
        eb0 = cn["redq_ebye"]["st0_vn0"]
        eb5 = cn["redq_ebye"]["st512_vn0"]
        hep.histplot(
            to_bh(cn["hists_ebye"]["delta_redQ0_st0"]), histtype="step", yerr=False, density=True, ax=ax4,
            color="#1f77b4", lw=1.8,
            label=f"2024 CaloNu status=0 (n={eb0['n']:,}, |Δ|<0.01: {100*eb0['frac_lt_0p01']:.1f}%)",
        )
        hep.histplot(
            to_bh(cn["hists_ebye"]["delta_redQ0_st512"]), histtype="step", yerr=False, density=True, ax=ax4,
            color="#d62728", lw=1.8, ls="--",
            label=f"2024 CaloNu status=512 (n={eb5['n']:,}, |Δ|<0.01: {100*eb5['frac_lt_0p01']:.1f}%)",
        )
        ax4.set_yscale("log")
        ax4.set_xlim(-20, 20)
        ax4.set_ylim(bottom=8e-6)
        ax4.set_xlabel(r"$\Delta Q_{\mathrm{red}}(\mathrm{VetoNu0}) = Q_{\mathrm{fix\_v1}} - Q_{\mathrm{Ben}}$ [pC]", fontsize=10.2)
        ax4.set_ylabel("Probability density [pC$^{-1}$]", fontsize=10.0)
        ax4.set_title("(d) 2024 CaloNu (run 015849): Dig2 status=0 vs status=512 (+16 ns)", fontsize=10.4)
        ax4.tick_params(direction="in", top=True, right=True)
        ax4.grid(True, color="#ececec", lw=0.6)
        ax4.legend(frameon=False, fontsize=8.2, loc="upper right")
        save(
            fig, "dq_02_reduced_charge_zoom_and_ebye_diff",
            "Zoomed nu_mu VetoNu reduced charge (<60 pC) and event-by-event fix_v1 - Ben differences including CaloNu status=512", pdf,
        )

        # ---- FIG 3: FASER DQ Fig. 5 Style — Scintillator & Calorimeter Charges (2x4 with ratio) ----
        fig, pairs = make_grid_with_ratios(
            2, 4, figsize=(15.5, 8.2),
            suptitle="FASER DQ Fig. 5 Style — Scintillator & Calorimeter Response: fix_v1 vs Official PHYS (2022–2024)",
        )
        scint_panels = [
            ("VetoNu0_charge",      "Charge [pC]",      r"FASER$\nu$ veto, layer 1",    False, (0.85, 1.15)),
            ("VetoNu1_charge",      "Charge [pC]",      r"FASER$\nu$ veto, layer 2",    False, (0.85, 1.15)),
            ("Veto20_charge",       "Charge [pC]",      "Veto station 2, layer 1",      False, (0.85, 1.15)),
            ("Veto21_charge",       "Charge [pC]",      "Veto station 2, layer 2",      False, (0.85, 1.15)),
            ("Timing_charge_total", "Total Charge [pC]","Timing station (sum 4 PMTs)",  False, (0.85, 1.15)),
            ("Preshower0_charge",   "Charge [pC]",      "Pre-shower station, layer 0",  False, (0.85, 1.15)),
            ("Preshower1_charge",   "Charge [pC]",      "Pre-shower station, layer 1",  False, (0.85, 1.15)),
            ("Calo_total_E_EM",     r"$E_{\mathrm{EM}}$ [GeV]", "Calorimeter total $E_{\\mathrm{EM}}$", True, (0.75, 1.25)),
        ]
        for idx, (hkey, xl, ttl, logy, ryl) in enumerate(scint_panels):
            ax_t, ax_b = pairs[idx]
            plot_dq_overlay_with_ratio(
                ax_t, ax_b, data, hkey, xl, ttl, logy=logy, ratio_ylim=ryl, show_legend=(idx == 3),
            )
        save(
            fig, "dq_03_scintillator_calo_dq_fig5_style",
            "FASER DQ Fig. 5-style scintillator and calorimeter distributions comparing fix_v1 vs Official PHYS across all 5 periods", pdf,
        )

        # ---- FIG 4: FASER DQ Fig. 6(a-g) Style — Track Kinematics (Common Physical FASER Frame, 2x4) ----
        fig, pairs = make_grid_with_ratios(
            2, 4, figsize=(15.5, 8.2),
            suptitle="FASER DQ Fig. 6 Style — Track Kinematics (Common Physical FASER Frame): fix_v1 vs Official PHYS",
        )
        trk_panels = [
            ("Track_Chi2",         r"$\chi^2$",              r"Track $\chi^2$",        False, None, None, (0.70, 1.30)),
            ("Track_nDoF",         "nDoF",                   "Track nDoF",             False, None, None, (0.70, 1.30)),
            ("Track_chi2_per_dof", r"$\chi^2$/NDF",          r"Track $\chi^2$/NDF",    False, 15.0, "15", (0.70, 1.30)),
            ("Track_pz_charge0",   r"$p_z \cdot q$ [GeV]",   r"Signed $p_z \cdot q$",  False, None, None, (0.75, 1.25)),
            ("Track_x0",           r"$x_0$ [mm]",            r"Upstream $x_0$",        False, None, None, (0.80, 1.20)),
            ("Track_y0",           r"$y_0$ [mm]",            r"Upstream $y_0$",        False, None, None, (0.80, 1.20)),
            ("Track_theta_x0",     r"$\theta_x$ [mrad]",     r"Upstream $\theta_x$",   False, None, None, (0.75, 1.25)),
            ("Track_theta_y0",     r"$\theta_y$ [mrad]",     r"Upstream $\theta_y$",   False, None, None, (0.75, 1.25)),
        ]
        for idx, (hkey, xl, ttl, logy, cl, clbl, ryl) in enumerate(trk_panels):
            ax_t, ax_b = pairs[idx]
            plot_dq_overlay_with_ratio(
                ax_t, ax_b, data, hkey, xl, ttl, logy=logy, cut_line=cl, cut_label=clbl,
                ratio_ylim=ryl, show_legend=(idx == 3),
            )
        save(
            fig, "dq_04_track_kinematics_dq_fig6_style",
            "FASER DQ Fig. 6-style track kinematics in common physical FASER frame (fix_v1 vs Official PHYS)", pdf,
        )

        # ---- FIG 5: FASER DQ Fig. 6(h-k) Style — Charge-Split Angles (2x2 with ratio) ----
        fig, pairs = make_grid_with_ratios(
            2, 2, figsize=(11.8, 8.2),
            suptitle=r"FASER DQ Fig. 6(h–k) Style — Charge-Separated Track Angles (Physical FASER Frame)",
        )
        ang_panels = [
            ("Track_theta_x0_pos", r"$\theta_{x, q>0}$ [mrad]", r"$\theta_x$ ($q > 0$)"),
            ("Track_theta_x0_neg", r"$\theta_{x, q<0}$ [mrad]", r"$\theta_x$ ($q < 0$)"),
            ("Track_theta_y0_pos", r"$\theta_{y, q>0}$ [mrad]", r"$\theta_y$ ($q > 0$)"),
            ("Track_theta_y0_neg", r"$\theta_{y, q<0}$ [mrad]", r"$\theta_y$ ($q < 0$)"),
        ]
        for idx, (hkey, xl, ttl) in enumerate(ang_panels):
            ax_t, ax_b = pairs[idx]
            plot_dq_overlay_with_ratio(
                ax_t, ax_b, data, hkey, xl, ttl, logy=False, ratio_ylim=(0.75, 1.25), show_legend=(idx == 1),
            )
        save(
            fig, "dq_05_track_angles_by_charge_dq_fig6hk_style",
            "FASER DQ Fig. 6(h-k)-style charge-split track angles (fix_v1 vs Official PHYS)", pdf,
        )

        # ---- FIG 6: nu_mu Analysis Selection Observables (LeadTrack + Multiplicity, 2x4 with ratio) ----
        fig, pairs = make_grid_with_ratios(
            2, 4, figsize=(15.5, 8.2),
            suptitle=r"$\nu_\mu$ Analysis Selection Observables (Scintillator-Preselected, Physical FASER Frame): fix_v1 vs Ben",
        )
        numu_panels = [
            ("LeadTrack_pz0",          r"Lead $p_z$ [GeV]",       r"Lead track $p_z$",            True,  100.0, ">100", (0.75, 1.25)),
            ("LeadTrack_Chi2_NDF",     r"Lead $\chi^2$/NDF",      r"Lead track $\chi^2$/NDF",     False, 15.0,  "<15",  (0.70, 1.30)),
            ("LeadTrack_nLayers",      "Lead nLayers",            "Lead track nLayers",           False, 6.5,   ">=7",  (0.70, 1.30)),
            ("LeadTrack_Theta",        r"Lead $\theta$ [mrad]",   r"Lead track polar $\theta$",   True,  25.0,  "<25",  (0.70, 1.30)),
            ("LeadTrack_r_atMaxRadius",r"$r_{\max}$ [mm]",        r"Track $r$ at max radius",     False, 95.0,  "<95",  (0.80, 1.20)),
            ("LeadTrack_rIFT",         r"$r_{\mathrm{IFT}}$ [mm]",r"Track $r$ at Veto2/IFT",      False, 95.0,  "<95",  (0.80, 1.20)),
            ("LeadTrack_rVetoNu",      r"$r_{\mathrm{VetoNu}}$ [mm]", r"Track $r$ at VetoNu",     False, 120.0, "<120", (0.80, 1.20)),
            ("longTracks",             "longTracks",              "Long-track multiplicity",      True,  0.5,   ">0",   (0.70, 1.30)),
        ]
        for idx, (hkey, xl, ttl, logy, cl, clbl, ryl) in enumerate(numu_panels):
            ax_t, ax_b = pairs[idx]
            plot_dq_overlay_with_ratio(
                ax_t, ax_b, data, hkey, xl, ttl, logy=logy, cut_line=cl, cut_label=clbl,
                ratio_ylim=ryl, show_legend=(idx == 0),
            )
        save(
            fig, "dq_06_numu_selection_observables",
            "nu_mu analysis selection observables comparing fix_v1 vs Ben across all 5 periods", pdf,
        )

        # ---- FIG 7: nu_mu Cutflow Ratio & Gained vs Lost Track Quality ----
        fig = plt.figure(figsize=(14.6, 8.0))
        gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.22, left=0.07, right=0.98, top=0.90, bottom=0.11)
        fig.suptitle(
            r"$\nu_\mu$ Cutflow Comparison & Gained (only fix_v1) vs Lost (only Official) Track Quality",
            fontsize=12.4, fontweight="bold", y=0.965,
        )

        # (a) Cutflow ratio fix_v1 / Ben across all pre-VetoNu steps
        ax_cf = fig.add_subplot(gs[0, :])
        cf_steps = [
            ("GoodTimes_notExcluded", "GRL+Colliding+Trig"),
            ("longTracks_gt0", "longTracks>0"),
            ("CaloNu_status_clean", "CaloNuStatus"),
            ("Veto2_gt40pC", "Veto2>40 pC"),
            ("Timing_gt20pC", "Timing>20 pC"),
            ("Preshower_gt2.5pC", "PS>2.5 pC"),
            ("LeadTrack_pz_gt100GeV", "pz>100 GeV"),
            ("LeadTrack_quality", "TrkQual"),
            ("LeadTrack_fiducial_angle", "Fid+θ<25"),
            ("RedQ_available", "RedQAvail"),
        ]
        x_cf = np.arange(len(cf_steps))
        for p in PERIOD_ORDER:
            col = PERIOD_META[p]["color"]
            ratios = []
            for sk, _ in cf_steps:
                n1 = data[p]["cutflow_fixv1"][sk]
                n2 = data[p]["cutflow_ben"][sk]
                ratios.append(n1 / n2 if n2 > 0 else np.nan)
            ax_cf.plot(x_cf, ratios, marker="o", ms=5.0, lw=1.7, color=col, label=PERIOD_META[p]["label"])
        ax_cf.axhline(1.0, color="#333333", ls="--", lw=1.0)
        ax_cf.set_xticks(x_cf, [lbl for _, lbl in cf_steps], fontsize=9.5)
        ax_cf.set_ylabel(r"Yield ratio $N(\mathrm{fix\_v1}) / N(\mathrm{Ben})$", fontsize=10.5)
        ax_cf.set_title(r"(a) Step-by-step $\nu_\mu$ pre-VetoNu cutflow yield ratio (common physical FASER frame)", fontsize=10.8)
        ax_cf.set_ylim(0.95, 1.15)
        ax_cf.grid(True, color="#e6e6e6", lw=0.6)
        ax_cf.tick_params(direction="in", top=True, right=True)
        ax_cf.legend(frameon=False, fontsize=8.8, ncol=5, loc="upper center")

        # (b) Gained vs Lost tracks: Chi2/nDoF & pz0 (summed across all 5 periods)
        ax_g1 = fig.add_subplot(gs[1, 0])
        ax_g2 = fig.add_subplot(gs[1, 1])
        for ax, gk, lk, xl, ttl in [
            (ax_g1, "gain_chi2_ndf", "loss_chi2_ndf", r"Lead track $\chi^2$/NDF", "(b) Gained vs Lost tracks: fit quality"),
            (ax_g2, "gain_pz0",      "loss_pz0",      r"Lead track $p_z$ [GeV]", "(c) Gained vs Lost tracks: momentum"),
        ]:
            edges = np.asarray(data[PERIOD_ORDER[0]]["hists_ebye"][gk]["edges"], dtype=np.float64)
            cg = sum(np.asarray(data[p]["hists_ebye"][gk]["counts"], dtype=np.float64) for p in PERIOD_ORDER)
            cl = sum(np.asarray(data[p]["hists_ebye"][lk]["counts"], dtype=np.float64) for p in PERIOD_ORDER)
            hg = bh.Histogram(bh.axis.Variable(edges), storage=bh.storage.Weight())
            hg.view(flow=False).value[:] = cg; hg.view(flow=False).variance[:] = cg
            hl = bh.Histogram(bh.axis.Variable(edges), storage=bh.storage.Weight())
            hl.view(flow=False).value[:] = cl; hl.view(flow=False).variance[:] = cl
            hep.histplot(hg, histtype="errorbar", ax=ax, color="#2ca02c", marker="o", markersize=3.5, xerr=True,
                         label=f"Only fix_v1 (gained, N={int(cg.sum()):,})")
            hep.histplot(hl, histtype="errorbar", ax=ax, color="#d62728", marker="s", markersize=3.5, xerr=True,
                         label=f"Only Official (lost, N={int(cl.sum()):,})")
            ax.set_xlabel(xl, fontsize=10.2)
            ax.set_ylabel("Events / bin (all 5 periods)", fontsize=9.8)
            ax.set_title(ttl, fontsize=10.5)
            ax.set_xlim(edges[0], edges[-1])
            ax.tick_params(direction="in", top=True, right=True)
            ax.grid(True, color="#e6e6e6", lw=0.6)
            ax.legend(frameon=False, fontsize=8.8, loc="upper right")
        ax_g2.set_yscale("log")
        save(
            fig, "dq_07_numu_cutflow_and_gained_tracks",
            "Step-by-step nu_mu cutflow yield ratio and gained-vs-lost track quality comparison", pdf,
        )

        # ---- FIG 8: Event-by-Event Tracking Residuals & Raw NtupleDumper Frame Audit ----
        fig = plt.figure(figsize=(14.8, 8.2))
        gs = fig.add_gridspec(2, 2, hspace=0.33, wspace=0.22, left=0.07, right=0.98, top=0.90, bottom=0.08)
        fig.suptitle(
            "Event-by-Event Tracking Closure (Physical FASER Frame) & Raw NtupleDumper Coordinate-Frame Audit",
            fontsize=12.2, fontweight="bold", y=0.965,
        )

        # (a) Event-by-event relative q/p shift
        ax1 = fig.add_subplot(gs[0, 0])
        for p in PERIOD_ORDER:
            col = PERIOD_META[p]["color"]
            hep.histplot(
                to_bh(data[p]["hists_ebye"]["delta_qoverp_rel"]), histtype="step", yerr=False, density=True, ax=ax1,
                color=col, lw=1.5,
                label=f"{PERIOD_META[p]['short']} (med={data[p]['trk_ebye']['d_qp_rel_med']:+.4f})",
            )
        ax1.axvline(0.0, color="#333333", ls="--", lw=1.0)
        ax1.set_yscale("log")
        ax1.set_xlim(-0.5, 0.5)
        ax1.set_xlabel(r"$(q/p_{\mathrm{fix\_v1}} - q/p_{\mathrm{off}}) \cdot p_{\mathrm{off}}$", fontsize=10.2)
        ax1.set_ylabel("Probability density", fontsize=9.8)
        ax1.set_title(r"(a) Event-by-event relative curvature shift $\Delta(q/p) \cdot p_{\mathrm{off}}$", fontsize=10.3)
        ax1.tick_params(direction="in", top=True, right=True)
        ax1.grid(True, color="#ececec", lw=0.6)
        ax1.legend(frameon=False, fontsize=7.8, loc="upper right")

        # (b) Charge-sign agreement vs official momentum bin
        ax2 = fig.add_subplot(gs[0, 1])
        x_pb = np.arange(5)
        p_lbls = ["<50", "50–150", "150–400", "400–1000", ">1000"]
        for p in PERIOD_ORDER:
            col = PERIOD_META[p]["color"]
            bdata = data[p]["trk_ebye"]["by_p_bin"]
            ys = [100.0 * b["q_agree_frac"] for b in bdata]
            ax2.plot(x_pb[:len(ys)], ys, marker="o", ms=5.0, lw=1.6, color=col, label=PERIOD_META[p]["short"])
        ax2.set_xticks(x_pb, p_lbls, fontsize=9.5)
        ax2.set_ylim(55, 102)
        ax2.set_xlabel(r"Official track momentum $p_{\mathrm{off}}$ [GeV]", fontsize=10.2)
        ax2.set_ylabel(r"Charge-sign agreement $P(q_{\mathrm{fix\_v1}} = q_{\mathrm{off}})$ [%]", fontsize=9.8)
        ax2.set_title("(b) Charge-sign agreement vs momentum (high-p straight-track turnover)", fontsize=10.3)
        ax2.tick_params(direction="in", top=True, right=True)
        ax2.grid(True, color="#ececec", lw=0.6)
        ax2.legend(frameon=False, fontsize=8.3, loc="lower left")

        # (c) Event-by-event delta_chi2
        ax3 = fig.add_subplot(gs[1, 0])
        for p in PERIOD_ORDER:
            col = PERIOD_META[p]["color"]
            hep.histplot(
                to_bh(data[p]["hists_ebye"]["delta_chi2"]), histtype="step", yerr=False, density=True, ax=ax3,
                color=col, lw=1.5,
                label=f"{PERIOD_META[p]['short']} (med Δχ²={data[p]['trk_ebye']['d_chi2_med']:+.2f})",
            )
        ax3.axvline(0.0, color="#333333", ls="--", lw=1.0)
        ax3.set_yscale("log")
        ax3.set_xlim(-20, 20)
        ax3.set_xlabel(r"$\Delta \chi^2 = \chi^2_{\mathrm{fix\_v1}} - \chi^2_{\mathrm{off}}$", fontsize=10.2)
        ax3.set_ylabel("Probability density", fontsize=9.8)
        ax3.set_title(r"(c) Event-by-event $\Delta\chi^2$ on matched tracks", fontsize=10.3)
        ax3.tick_params(direction="in", top=True, right=True)
        ax3.grid(True, color="#ececec", lw=0.6)
        ax3.legend(frameon=False, fontsize=7.9, loc="upper right")

        # (d) Raw NtupleDumper coordinate-frame audit (2024 CaloNu example)
        ax4 = fig.add_subplot(gs[1, 1])
        cn = data["2024_caloNu"]
        hep.histplot(
            to_bh(cn["hists_fixv1"]["raw_Track_pz0"]), histtype="step", yerr=False, density=True, ax=ax4,
            color="#d62728", lw=1.8, label=r"fix_v1 RAW `Track_pz0` (= $-p_x^{\mathrm{FASER}}$, ~0 GeV!)",
        )
        hep.histplot(
            to_bh(cn["hists_ben"]["raw_Track_pz0"]), histtype="step", yerr=False, density=True, ax=ax4,
            color="#111111", lw=1.5, ls="--", label=r"Official `Track_pz0` (= $+p_z^{\mathrm{FASER}}$, >20 GeV)",
        )
        ax4.set_xlim(-25, 50)
        ax4.set_ylim(bottom=5e-5)
        ax4.set_yscale("log")
        ax4.set_xlabel("Raw branch `Track_pz0` [GeV] (before $R^T$ un-rotation)", fontsize=10.0)
        ax4.set_ylabel("Probability density [GeV$^{-1}$]", fontsize=9.8)
        ax4.set_title("(d) Raw NtupleDumper frame audit: rotated ACTS frame written verbatim", fontsize=10.2)
        ax4.tick_params(direction="in", top=True, right=True)
        ax4.grid(True, color="#ececec", lw=0.6)
        ax4.legend(frameon=False, fontsize=8.0, loc="upper right")
        save(
            fig, "dq_08_ebye_tracking_and_frame_audit",
            "Event-by-event tracking residuals (q/p, charge agreement, chi2) and raw NtupleDumper coordinate-frame audit", pdf,
        )

    (fig_dir / "figure_index.json").write_text(json.dumps({
        "style_reference": "FASER_DQ_Analyis (notebooks/plotting_utils.py, internal_note/finalize.py, Figs. 5-6)",
        "report_pdf": str(report_pdf),
        "figures": index,
    }, indent=2) + "\n")
    print(f"Generated {len(index)} FASER DQ-style figures + {report_pdf}")


if __name__ == "__main__":
    main()
