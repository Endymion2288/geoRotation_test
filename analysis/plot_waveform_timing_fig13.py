#!/usr/bin/env python3
"""Generate Fig. 13-style English figures for the 2022-2026 & MC24 timing study.

Reads timing_summary.json and timing_histograms.json produced by
analysis/waveform_timing_study.py and writes PNG, per-figure PDF,
a combined multi-page PDF report, and figure_index.json.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

DEFAULT_DIR = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/"
    "stage04_waveform_reduced_charge/timing_study_fig13"
)

PERIOD_STYLES = {
    "2022": ("2022 data", "#1f77b4", "-", 1.8),
    "2023": ("2023 data", "#ff7f0e", "-", 1.8),
    "2024_preCaloNu": ("2024 pre-CaloNu", "#2ca02c", "-", 1.8),
    "2024_earlyCaloNu": ("2024 early-CaloNu (15821-15844)", "#17becf", "--", 1.5),
    "2024_duringCaloNu": ("2024 during-CaloNu (15845-16924)", "#d62728", "-", 2.2),
    "2024_postCaloNu": ("2024 post-CaloNu", "#9467bd", "-", 1.8),
    "2025": ("2025 data", "#8c564b", "-", 1.9),
    "2026": ("2026 data", "#e377c2", "-", 1.9),
    "MC24_Muon_FaserNu": ("MC24 Muon (FaserNu, 100124)", "#222222", "--", 1.9),
    "MC24_Muon_CaloNu": ("MC24 Muon (CaloNu, 100050)", "#666666", ":", 1.8),
    "MC24_Genie_22_23": ("MC24 GENIE 22/23 default", "#444444", "-.", 1.5),
}


def setup_style():
    plt.rcParams.update({
        "font.size": 10.5,
        "axes.labelsize": 11.5,
        "axes.titlesize": 12.0,
        "legend.fontsize": 8.8,
        "xtick.labelsize": 10.0,
        "ytick.labelsize": 10.0,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linestyle": "--",
        "figure.dpi": 150,
    })


def hist_xy_normalized(counts: List[float], spec: dict, rebin: int = 1):
    arr = np.asarray(counts, dtype=np.float64)
    nb = len(arr)
    if rebin > 1 and nb % rebin == 0:
        arr = arr.reshape(-1, rebin).sum(axis=1)
        nb = len(arr)
    edges = np.linspace(spec["xmin"], spec["xmax"], nb + 1)
    total = arr.sum()
    if total > 0:
        # Normalize to unit area (density * bin_width -> sum(y) = 1, matching Fig. 13 style)
        y = arr / total
    else:
        y = arr
    return edges, y, int(round(total))


def step_plot(ax, edges, y, label, color, ls="-", lw=1.8, alpha=0.95):
    ax.stairs(y, edges, label=label, color=color, linestyle=ls, linewidth=lw, alpha=alpha)


def save_fig(fig, out_dir: Path, stem: str, pdf_pages: PdfPages, index_entries: list, title: str):
    png_path = out_dir / f"{stem}.png"
    pdf_path = out_dir / f"{stem}.pdf"
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    pdf_pages.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    index_entries.append({
        "stem": stem,
        "title": title,
        "png": str(png_path),
        "pdf": str(pdf_path),
    })


def main():
    parser = argparse.ArgumentParser(description="Plot Fig. 13-style timing-difference figures")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DIR)
    args = parser.parse_args()

    setup_style()
    summary = json.loads((args.data_dir / "timing_summary.json").read_text())
    hists_doc = json.loads((args.data_dir / "timing_histograms.json").read_text())
    specs = hists_doc["hist_specs"]
    hists = hists_doc["histograms"]
    periods = summary["periods"]

    fig_dir = args.data_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    report_pdf = fig_dir / "timing_study_fig13_report.pdf"
    index_entries = []

    selection_note = (
        "Selection: FASER_DQ_RDF GoodTracks (nLayers>=7, chi2/ndof<25, nHits>=12, pz0>20 GeV) "
        "+ annual GRL + colliding BCID + TAP&4 + no timing saturation; (status & ~512)==0"
    )

    with PdfPages(report_pdf) as pdf:
        # ------------------------------------------------------------------
        # Figure 1: Direct Fig. 13 reproduction & extension (Veto2 - VetoNu)
        # ------------------------------------------------------------------
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 5.4))
        plot_order = [
            "2022", "2023", "2024_preCaloNu", "2024_earlyCaloNu",
            "2024_duringCaloNu", "2024_postCaloNu", "2025", "2026",
            "MC24_Muon_FaserNu",
        ]
        for p in plot_order:
            if p not in hists or "veto2_vetonu_raw" not in hists[p]:
                continue
            label_base, col, ls, lw = PERIOD_STYLES[p]
            st_raw = periods[p]["stats"]["veto2_vetonu_raw"]
            st_corr = periods[p]["stats"]["veto2_vetonu_corr"]
            e1, y1, n1 = hist_xy_normalized(hists[p]["veto2_vetonu_raw"], specs["veto2_vetonu_raw"], rebin=2)
            e2, y2, n2 = hist_xy_normalized(hists[p]["veto2_vetonu_corr"], specs["veto2_vetonu_corr"], rebin=1)
            step_plot(ax1, e1, y1, f"{label_base} (N={n1:,})", col, ls, lw)
            step_plot(
                ax2, e2, y2,
                f"{label_base}: mode={st_corr['mode_ns']:.2f} ns",
                col, ls, lw,
            )

        ax1.set_xlim(10.0, 50.0)
        ymax1 = ax1.get_ylim()[1]
        ax1.set_ylim(0.0, ymax1 * 1.12)
        ax1.set_xlabel(r"$t_{\mathrm{Veto2}}^{\mathrm{raw}} - t_{\mathrm{VetoNu}}$ [ns]")
        ax1.set_ylabel("Normalized (unit area / 0.4 ns bin)")
        ax1.set_title("(a) Raw localtime (Fig. 13 style: shows -16 ns bit-512 split & Dig2 shift)")
        ax1.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=8.0)
        ax1.annotate(
            "CaloNu status==512\n(-16 ns Dig2 trigger shift)",
            xy=(19.8, ymax1 * 0.13),
            xytext=(10.8, ymax1 * 0.42),
            arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.3),
            fontsize=8.8, color="#d62728", fontweight="bold",
        )

        ax2.set_xlim(21.5, 42.5)
        ymax2 = ax2.get_ylim()[1]
        ax2.set_ylim(0.0, ymax2 * 1.15)
        ax2.set_xlabel(r"$t_{\mathrm{Veto2}}^{\mathrm{+16\,ns\,corr}} - t_{\mathrm{VetoNu}}$ [ns]")
        ax2.set_ylabel("Normalized (unit area / 0.2 ns bin)")
        ax2.set_title("(b) After +16 ns correction: residual +2.7 ns CaloNu Dig2 & MC offsets")
        ax2.axvline(33.165, color="#1f77b4", ls=":", lw=1.1, alpha=0.7)
        ax2.axvline(35.783, color="#d62728", ls=":", lw=1.1, alpha=0.7)
        ax2.axvline(27.21, color="#222222", ls=":", lw=1.1, alpha=0.7)
        ax2.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=8.0)

        fig.suptitle(
            "Throughgoing-Muon Timing Difference: Veto St2 - VetoNu Across 2022-2026 Data & MC24",
            fontsize=13.0, fontweight="bold", y=1.01,
        )
        fig.text(0.5, -0.03, selection_note, ha="center", fontsize=8.8, color="#333333")
        fig.tight_layout()
        save_fig(fig, fig_dir, "fig01_veto2_vetonu_fig13_all_periods", pdf, index_entries,
                 "Veto2 - VetoNu timing distributions (raw and +16 ns corrected) across 2022-2026 & MC24")

        # ------------------------------------------------------------------
        # Figure 2: 2024 CaloNu status == 0 vs status == 512 decomposition
        # ------------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(9.5, 5.2))
        if "2024_preCaloNu" in hists:
            e, y, n = hist_xy_normalized(hists["2024_preCaloNu"]["veto2_vetonu_corr"], specs["veto2_vetonu_corr"], 1)
            m = periods["2024_preCaloNu"]["stats"]["veto2_vetonu_corr"]["mode_ns"]
            step_plot(ax, e, y, f"2024 pre-CaloNu (Dig1 ch 6,7): mode={m:.2f} ns (N={n:,})", "#2ca02c", "-", 2.0)
        if "2024_duringCaloNu" in hists:
            p = "2024_duringCaloNu"
            for hk, lbl, col, ls, lw in [
                ("veto2_vetonu_s0", "CaloNu Dig2 ch 22,23 (status == 0)", "#d62728", "-", 2.1),
                ("veto2_vetonu_s512_raw", "CaloNu Dig2 ch 22,23 (status == 512, raw)", "#ff7f0e", "--", 2.0),
                ("veto2_vetonu_s512_corr", "CaloNu Dig2 ch 22,23 (status == 512, +16 ns)", "#1f77b4", "-.", 2.1),
            ]:
                e, y, n = hist_xy_normalized(hists[p][hk], specs[hk], 1)
                m = periods[p]["stats"][hk]["mode_ns"]
                step_plot(ax, e, y, f"{lbl}: mode={m:.2f} ns (N={n:,})", col, ls, lw)
        ax.set_xlim(12.0, 44.0)
        ax.set_xlabel(r"$t_{\mathrm{Veto2}} - t_{\mathrm{VetoNu}}$ [ns]")
        ax.set_ylabel("Normalized (unit area / 0.2 ns bin)")
        ax.set_title("2024 CaloNu Digitizer-2 Decomposition: Discrete -16 ns Shift (bit 512) + Fixed +2.6 ns Board/Cable Offset")
        ax.legend(loc="upper left", frameon=True, framealpha=0.92)
        fig.text(0.5, -0.03, selection_note, ha="center", fontsize=8.8, color="#333333")
        fig.tight_layout()
        save_fig(fig, fig_dir, "fig02_caloNu_status512_decomposition", pdf, index_entries,
                 "2024 CaloNu Veto2 status==0 vs status==512 decomposition")

        # ------------------------------------------------------------------
        # Figure 3: Preshower (Dig1: ch 12,13) - VetoNu across 2022-2026 & MC24
        # ------------------------------------------------------------------
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 5.4), gridspec_kw={"width_ratios": [1.35, 1.0]})
        for p in plot_order:
            if p not in hists or "ps_dig1_vetonu" not in hists[p]:
                continue
            label_base, col, ls, lw = PERIOD_STYLES[p]
            st = periods[p]["stats"]["ps_dig1_vetonu"]
            e, y, n = hist_xy_normalized(hists[p]["ps_dig1_vetonu"], specs["ps_dig1_vetonu"], rebin=2)
            step_plot(ax1, e, y, f"{label_base}: mode={st['mode_ns']:.2f} ns (N={n:,})", col, ls, lw)
            if p in ("2022", "2023", "2024_preCaloNu", "2024_duringCaloNu", "2024_postCaloNu", "MC24_Muon_FaserNu"):
                e_z, y_z, _ = hist_xy_normalized(hists[p]["ps_dig1_vetonu"], specs["ps_dig1_vetonu"], rebin=1)
                step_plot(ax2, e_z, y_z, f"{label_base}: {st['mode_ns']:.2f} ns", col, ls, lw)

        ax1.set_xlim(0.0, 66.0)
        ym1 = ax1.get_ylim()[1]
        ax1.set_ylim(0.0, ym1 * 1.12)
        ax1.set_xlabel(r"$t_{\mathrm{Preshower\,Dig1\,(ch\,12,13)}} - t_{\mathrm{VetoNu}}$ [ns]")
        ax1.set_ylabel("Normalized (unit area / 0.4 ns bin)")
        ax1.set_title("(a) Full range: -40.3 ns shift after 2025 Preshower upgrade")
        ax1.legend(loc="upper left", frameon=True, framealpha=0.92, fontsize=8.0)
        ax1.annotate(
            "-40.3 ns hardware shift\n(new 2025-2026 PMTs/cables;\nCOOL TriggerOffset: +18 -> -40 ns)",
            xy=(17.0, ym1 * 0.16),
            xytext=(24.5, ym1 * 0.25),
            arrowprops=dict(arrowstyle="->", color="#8c564b", lw=1.4),
            fontsize=8.6, color="#8c564b", fontweight="bold",
        )

        ax2.set_xlim(46.0, 64.0)
        ym2 = ax2.get_ylim()[1]
        ax2.set_ylim(0.0, ym2 * 1.15)
        ax2.set_xlabel(r"$t_{\mathrm{Preshower\,Dig1\,(ch\,12,13)}} - t_{\mathrm{VetoNu}}$ [ns]")
        ax2.set_ylabel("Normalized (unit area / 0.2 ns bin)")
        ax2.set_title("(b) Zoom 2022-2024 (~53.5 ns) vs MC24 (~51.4 ns)")
        ax2.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=8.2)

        fig.suptitle(
            "Throughgoing-Muon Timing Difference: Preshower (Dig1: ch 12,13) - VetoNu Across 2022-2026 & MC24",
            fontsize=13.0, fontweight="bold", y=1.01,
        )
        fig.text(0.5, -0.03, selection_note, ha="center", fontsize=8.8, color="#333333")
        fig.tight_layout()
        save_fig(fig, fig_dir, "fig03_preshower_dig1_vetonu_all_periods", pdf, index_entries,
                 "Preshower Dig1 (ch 12,13) - VetoNu timing across 2022-2026 & MC24")

        # ------------------------------------------------------------------
        # Figure 4: 2025-2026 Preshower Dig2 (ch 24,25) extra offset study
        # ------------------------------------------------------------------
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16.5, 5.1))
        for p, col_base, ls in [("2025", "#8c564b", "-"), ("2026", "#e377c2", "--")]:
            if p not in hists or "ps_dig2_vetonu_corr" not in hists[p]:
                continue
            # Panel 1: Dig1 vs Dig2 raw vs Dig2 +16ns corr
            for hk, lbl, c, sty in [
                ("ps_dig1_vetonu", f"{p} Dig1 (ch 12,13)", "#1f77b4" if p == "2025" else "#17becf", ls),
                ("ps_dig2_vetonu_raw", f"{p} Dig2 (ch 24,25) raw", "#ff7f0e" if p == "2025" else "#bcbd22", ls),
                ("ps_dig2_vetonu_corr", f"{p} Dig2 (ch 24,25) +16 ns", "#d62728" if p == "2025" else "#9467bd", ls),
            ]:
                e, y, n = hist_xy_normalized(hists[p][hk], specs[hk], 2)
                m = periods[p]["stats"][hk]["mode_ns"]
                step_plot(ax1, e, y, f"{lbl}: {m:.2f} ns", c, sty, 1.8)

            # Panel 2: Direct event-by-event Dig2 - Dig1
            e, y, n = hist_xy_normalized(hists[p]["ps_dig2_minus_dig1_corr"], specs["ps_dig2_minus_dig1_corr"], 1)
            st_d = periods[p]["stats"]["ps_dig2_minus_dig1_corr"]
            step_plot(
                ax2, e, y,
                f"{p}: mode={st_d['mode_ns']:+.2f} ns, mean={st_d['core_mean_ns']:+.2f}+/-{st_d['core_std_ns']:.2f} (N={n:,})",
                col_base, ls, 2.0,
            )

            # Panel 3: Combined 4-PMT Preshower mean (raw vs +16ns only vs +16ns + Dig2 aligned)
            if p == "2025":
                for hk, lbl, c, sty in [
                    ("ps_all_vetonu_raw", "2025 4-PMT raw (uncorrected)", "#ff7f0e", ":"),
                    ("ps_all_vetonu_corr16", "2025 4-PMT +16 ns only (no Dig2 offset)", "#d62728", "--"),
                    ("ps_all_vetonu_aligned", "2025 4-PMT +16 ns & -4.0 ns Dig2 align", "#2ca02c", "-"),
                ]:
                    e, y, n = hist_xy_normalized(hists[p][hk], specs[hk], 1)
                    st_a = periods[p]["stats"][hk]
                    step_plot(
                        ax3, e, y,
                        f"{lbl}\n  mode={st_a['mode_ns']:.2f}, core_std={st_a['core_std_ns']:.2f} ns",
                        c, sty, 2.0,
                    )
            elif p == "2026":
                e, y, n = hist_xy_normalized(hists[p]["ps_all_vetonu_aligned"], specs["ps_all_vetonu_aligned"], 1)
                st_a = periods[p]["stats"]["ps_all_vetonu_aligned"]
                step_plot(
                    ax3, e, y,
                    f"2026 4-PMT +16 ns & -4.0 ns Dig2 align\n  mode={st_a['mode_ns']:.2f}, core_std={st_a['core_std_ns']:.2f} ns",
                    "#1f77b4", "-.", 1.9,
                )

        ax1.set_xlim(-6.0, 26.0)
        ax1.set_ylim(0.0, ax1.get_ylim()[1] * 1.28)
        ax1.set_xlabel(r"$t_{\mathrm{Preshower}} - t_{\mathrm{VetoNu}}$ [ns]")
        ax1.set_ylabel("Normalized (unit area / 0.4 ns bin)")
        ax1.set_title("(a) Dig1 (ch 12,13) vs Dig2 (ch 24,25)")
        ax1.legend(loc="upper left", fontsize=7.6, frameon=True, framealpha=0.92)

        ax2.set_xlim(-4.0, 12.0)
        ax2.set_ylim(0.0, ax2.get_ylim()[1] * 1.24)
        ax2.axvline(4.0, color="#222222", ls=":", lw=1.2)
        ax2.set_xlabel(r"$t_{\mathrm{PS\,Dig2}}^{\mathrm{+16\,ns}} - t_{\mathrm{PS\,Dig1}}$ [ns]")
        ax2.set_ylabel("Normalized (unit area / 0.2 ns bin)")
        ax2.set_title("(b) Event-by-event Dig2 - Dig1 extra offset")
        ax2.legend(loc="upper center", fontsize=7.8, frameon=True, framealpha=0.92)

        ax3.set_xlim(4.0, 25.0)
        ax3.set_ylim(0.0, ax3.get_ylim()[1] * 1.32)
        ax3.set_xlabel(r"$t_{\mathrm{Preshower\,(all\,4\,PMTs)}} - t_{\mathrm{VetoNu}}$ [ns]")
        ax3.set_ylabel("Normalized (unit area / 0.2 ns bin)")
        ax3.set_title("(c) Effect of Dig2 alignment on 4-PMT mean")
        ax3.legend(loc="upper right", fontsize=7.6, frameon=True, framealpha=0.92)

        fig.suptitle(
            "2025-2026 Upgraded 4-PMT Preshower: Why PreshowerDig2ExtraOffset (~+4.0 ns) Is Required",
            fontsize=13.0, fontweight="bold", y=1.02,
        )
        fig.text(0.5, -0.03, selection_note, ha="center", fontsize=8.8, color="#333333")
        fig.tight_layout()
        save_fig(fig, fig_dir, "fig04_preshower_dig2_extra_offset_2025_2026", pdf, index_entries,
                 "2025-2026 Preshower Dig2 vs Dig1 extra offset and 4-PMT alignment")

        # ------------------------------------------------------------------
        # Figure 5: All stations localtime & triggertime summary across periods
        # ------------------------------------------------------------------
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 5.4))
        p_list = [p for p in [
            "2022", "2023", "2024_preCaloNu", "2024_earlyCaloNu",
            "2024_duringCaloNu", "2024_postCaloNu", "2025", "2026",
            "MC24_Muon_FaserNu", "MC24_Muon_CaloNu",
        ] if p in periods]
        x = np.arange(len(p_list))
        short_map = {
            "2022": "2022",
            "2023": "2023",
            "2024_preCaloNu": "24\npre",
            "2024_earlyCaloNu": "24\nearly",
            "2024_duringCaloNu": "24\nCaloNu",
            "2024_postCaloNu": "24\npost",
            "2025": "2025",
            "2026": "2026",
            "MC24_Muon_FaserNu": "MC24\nFaserNu",
            "MC24_Muon_CaloNu": "MC24\nCaloNu",
        }
        short_lbls = [short_map.get(p, p) for p in p_list]

        for hk, lbl, col, mk in [
            ("veto1_vetonu", "Veto1 - VetoNu (localtime)", "#17becf", "v"),
            ("veto2_vetonu_corr", "Veto2 - VetoNu (localtime, +16 ns)", "#d62728", "o"),
            ("timing_vetonu", "Timing - VetoNu (localtime)", "#2ca02c", "s"),
            ("ps_dig1_vetonu", "Preshower Dig1 - VetoNu (localtime)", "#1f77b4", "^"),
            ("ps_dig2_vetonu_corr", "Preshower Dig2 - VetoNu (localtime, +16 ns)", "#9467bd", "D"),
        ]:
            ys, errs, xs = [], [], []
            for i, p in enumerate(p_list):
                st = periods[p]["stats"].get(hk, {})
                if st.get("entries", 0) > 100:
                    xs.append(i)
                    ys.append(st["mode_ns"])
                    errs.append(st["core_std_ns"])
            ax1.errorbar(xs, ys, yerr=errs, label=lbl, color=col, marker=mk, lw=1.6, ms=6, capsize=3)

        ax1.set_xticks(x)
        ax1.set_xticklabels(short_lbls, fontsize=8.5)
        ax1.set_ylabel(r"Peak $\Delta t_{\mathrm{localtime}}$ [ns] (bars: core $\sigma$)")
        ax1.set_title("(a) Station localtime differences across 2022-2026 & MC24")
        ax1.legend(loc="lower left", fontsize=8.3, frameon=True, framealpha=0.92)

        # Panel 2: Triggertime differences (data only) showing COOL TriggerOffset residuals
        data_p = [p for p in p_list if not p.startswith("MC24")]
        xd = np.arange(len(data_p))
        d_lbls = [short_map.get(p, p) for p in data_p]
        for hk, lbl, col, mk in [
            ("veto2_vetonu_trig", "Veto2 - VetoNu (triggertime)", "#d62728", "o"),
            ("ps_dig1_vetonu_trig", "Preshower Dig1 - VetoNu (triggertime)", "#1f77b4", "^"),
            ("ps_dig2_vetonu_trig", "Preshower Dig2 - VetoNu (triggertime)", "#9467bd", "D"),
        ]:
            ys, errs, xs = [], [], []
            for i, p in enumerate(data_p):
                st = periods[p]["stats"].get(hk, {})
                if st.get("entries", 0) > 100:
                    xs.append(i)
                    ys.append(st["mode_ns"])
                    errs.append(st["core_std_ns"])
            ax2.errorbar(xs, ys, yerr=errs, label=lbl, color=col, marker=mk, lw=1.6, ms=6, capsize=3)

        ax2.axhline(0.0, color="#444444", ls=":", lw=1.0)
        ax2.set_xticks(xd)
        ax2.set_xticklabels(d_lbls, fontsize=8.5)
        ax2.set_ylabel(r"Peak $\Delta t_{\mathrm{triggertime}}$ [ns] (bars: core $\sigma$)")
        ax2.set_title("(b) Triggertime differences (COOL TriggerOffset is coarse, not zeroed!)")
        ax2.legend(loc="upper left", fontsize=8.5, frameon=True, framealpha=0.92)

        fig.suptitle(
            "All Scintillator Stations: Localtime vs Triggertime Evolution (2022-2026 & MC24)",
            fontsize=13.0, fontweight="bold", y=1.01,
        )
        fig.text(0.5, -0.03, selection_note, ha="center", fontsize=8.8, color="#333333")
        fig.tight_layout()
        save_fig(fig, fig_dir, "fig05_station_localtime_and_triggertime_evolution", pdf, index_entries,
                 "All scintillator stations localtime and triggertime evolution")

        # ------------------------------------------------------------------
        # Figure 6: Per-run stability timeline across 2022-2026
        # ------------------------------------------------------------------
        data_runs = [r for r in summary["per_run"] if not r["is_mc"] and r["goodtracks"] > 500]
        data_runs.sort(key=lambda r: r["run"])
        if data_runs:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14.5, 9.5), sharex=True)
            r_idx = np.arange(len(data_runs))
            r_nums = [r["run"] for r in data_runs]

            # Period shading
            period_colors = {
                "2022": "#e8f2fa",
                "2023": "#fff3e6",
                "2024_preCaloNu": "#eaf7ea",
                "2024_earlyCaloNu": "#e0f7fa",
                "2024_duringCaloNu": "#fdeaea",
                "2024_postCaloNu": "#f3ebf9",
                "2025": "#f5eee9",
                "2026": "#fcecf5",
            }
            p_spans = []
            cur_p, start_i = None, 0
            for i, r in enumerate(data_runs):
                if r["period"] != cur_p:
                    if cur_p is not None:
                        p_spans.append((cur_p, start_i, i - 1))
                    cur_p, start_i = r["period"], i
            if cur_p is not None:
                p_spans.append((cur_p, start_i, len(data_runs) - 1))

            for ax in (ax1, ax2, ax3):
                for p, i0, i1 in p_spans:
                    ax.axvspan(i0 - 0.5, i1 + 0.5, color=period_colors.get(p, "#f5f5f5"), alpha=0.85, zorder=0)

            # Panel 1: Veto2 - VetoNu
            v2_x, v2_y = [], []
            for i, r in enumerate(data_runs):
                d = r.get("veto2_vetonu_corr", {})
                if d.get("entries", 0) > 200:
                    v2_x.append(i)
                    v2_y.append(d["mode_ns"])
            ax1.plot(v2_x, v2_y, "o-", color="#d62728", ms=4.2, lw=1.2, label="Veto2 - VetoNu (+16 ns corr)")
            ax1.axhline(33.165, color="#1f77b4", ls=":", lw=1.1, label="Standard ref (33.165 ns)")
            ax1.axhline(35.783, color="#d62728", ls="--", lw=1.1, label="CaloNu ref (35.783 ns)")
            ax1.set_ylim(31.2, 37.9)
            ax1.set_ylabel(r"$t_{\mathrm{Veto2}} - t_{\mathrm{VN}}$ [ns]")
            ax1.legend(loc="lower left", ncol=3, fontsize=8.5, framealpha=0.92)
            for p, i0, i1 in p_spans:
                ax1.text(0.5 * (i0 + i1), 37.65, p.replace("2024_", "24\n"),
                         ha="center", va="top", fontsize=7.6, fontweight="bold", color="#333333")

            # Panel 2: Preshower Dig1 & Dig2 - VetoNu
            p1_x, p1_y, p2_x, p2_y = [], [], [], []
            for i, r in enumerate(data_runs):
                d1 = r.get("ps_dig1_vetonu", {})
                if d1.get("entries", 0) > 200:
                    p1_x.append(i)
                    p1_y.append(d1["mode_ns"])
                d2 = r.get("ps_dig2_vetonu_corr", {})
                if d2.get("entries", 0) > 200:
                    p2_x.append(i)
                    p2_y.append(d2["mode_ns"])
            ax2.plot(p1_x, p1_y, "^-", color="#1f77b4", ms=4.5, lw=1.2, label="Preshower Dig1 (ch 12,13) - VetoNu")
            ax2.plot(p2_x, p2_y, "D-", color="#9467bd", ms=4.2, lw=1.2, label="Preshower Dig2 (ch 24,25) - VetoNu (+16 ns)")
            ax2.set_ylim(5.0, 60.0)
            ax2.set_ylabel(r"$t_{\mathrm{PS}} - t_{\mathrm{VN}}$ [ns]")
            ax2.legend(loc="lower left", ncol=2, fontsize=8.5)

            # Panel 3: Waveform status fractions (Veto20 bit 512, PS0 GFIT_FAILED, PS2 valid)
            s512_f, ps0_gf, ps2_vf = [], [], []
            for r in data_runs:
                sc = r["status_counts"]
                n = max(1, sc.get("n_goodtrk", 1))
                s512_f.append(100.0 * sc.get("v20_s512", 0) / n)
                ps0_gf.append(100.0 * sc.get("ps0_gfit", 0) / n)
                ps2_vf.append(100.0 * (sc.get("ps2_s0", 0) + sc.get("ps2_s512", 0)) / n)
            ax3.plot(r_idx, s512_f, "o-", color="#d62728", ms=3.8, lw=1.2, label="Veto20 status==512 [%]")
            ax3.plot(r_idx, ps0_gf, "s-", color="#ff7f0e", ms=3.8, lw=1.2, label="Preshower0 GFIT_FAILED (bit 16) [%]")
            ax3.plot(r_idx, ps2_vf, "^-", color="#2ca02c", ms=3.8, lw=1.2, label="Preshower2 valid (status==0 or 512) [%]")
            ax3.set_ylim(-5.0, 108.0)
            ax3.set_ylabel("Fraction of GoodTracks [%]")
            ax3.legend(loc="upper left", ncol=3, fontsize=8.5)

            step_tick = max(1, len(data_runs) // 25)
            ticks = list(range(0, len(data_runs), step_tick))
            ax3.set_xticks(ticks)
            ax3.set_xticklabels([str(r_nums[i]) for i in ticks], rotation=45, ha="right", fontsize=8.5)
            ax3.set_xlabel("Run Number (chronological GRL runs, 2022 - 2026)")

            fig.suptitle(
                "Per-Run Stability of Station Timing Offsets and Waveform Status Fractions (2022-2026)",
                fontsize=13.0, fontweight="bold", y=0.995,
            )
            fig.text(0.5, -0.01, selection_note, ha="center", fontsize=8.8, color="#333333")
            fig.tight_layout()
            save_fig(fig, fig_dir, "fig06_per_run_stability_timeline_2022_2026", pdf, index_entries,
                     "Per-run stability timeline of Veto2/Preshower timing and waveform status fractions")

    index_path = fig_dir / "figure_index.json"
    index_path.write_text(json.dumps({
        "report_pdf": str(report_pdf),
        "figures": index_entries,
    }, indent=2))
    print(f"Generated {len(index_entries)} figures + {report_pdf}")


if __name__ == "__main__":
    main()
