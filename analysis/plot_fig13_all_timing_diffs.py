#!/usr/bin/env python3
"""Produce exact Fig. 13-style timing-difference plots (modeled on
Electronic_2D_measurement_internal_note-6.pdf, p. 22) across all years
(2022, 2023, 2024 pre/early/during/post-CaloNu, 2025, 2026, and MC24)
for every configured reduced-charge timing offset.

Matches Fig. 13 visual conventions:
  - Single-panel clean FASER Internal canvas with official FASER logo above top-left
  - Title centered above top spine: "All years — <observable>"
  - Solid box spines, inward major + minor ticks on all 4 sides, light solid grey grid
  - Area-normalized probability density (1.0 ns bins on [-60, 60] ns), y-label "Normalized"
    at top-left, x-label right-aligned at bottom-right
  - Clean step histograms with lower-left boxed legend using exact Fig. 13 period naming
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np

DEFAULT_DATA_DIR = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/"
    "stage04_waveform_reduced_charge/timing_study_fig13"
)
LOGO_PATH = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/"
    "stage04_waveform_reduced_charge/ref_pages/faser_logo_600dpi.png"
)

# Exact Fig. 13 color sequence for 2022-2024 + natural tab10 extension for 2025, 2026, MC24
PERIOD_CFG = {
    "2022": ("2022", "#1f77b4", "-", 1.6),
    "2023": ("2023", "#ff7f0e", "-", 1.6),
    "2024_preCaloNu": ("2024_preCaloNu", "#2ca02c", "-", 1.6),
    "2024_earlyCaloNu": ("2024_earlyCaloNu", "#17becf", "--", 1.5),
    "2024_duringCaloNu": ("2024_duringCaloNu", "#d62728", "-", 1.8),
    "2024_postCaloNu": ("2024_postCaloNu", "#9467bd", "-", 1.6),
    "2025": ("2025", "#8c564b", "-", 1.7),
    "2026": ("2026", "#e377c2", "-", 1.7),
    "MC24_Muon_FaserNu": ("MC24", "#111111", "--", 1.6),
}

ALL_PERIODS = [
    "2022",
    "2023",
    "2024_preCaloNu",
    "2024_earlyCaloNu",
    "2024_duringCaloNu",
    "2024_postCaloNu",
    "2025",
    "2026",
    "MC24_Muon_FaserNu",
]


def configure_rc():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "mathtext.default": "regular",
        "axes.linewidth": 0.95,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 6.5,
        "ytick.major.size": 6.5,
        "xtick.minor.size": 3.5,
        "ytick.minor.size": 3.5,
        "xtick.major.width": 0.85,
        "ytick.major.width": 0.85,
        "xtick.minor.width": 0.65,
        "ytick.minor.width": 0.65,
        "xtick.labelsize": 11.5,
        "ytick.labelsize": 11.5,
        "figure.dpi": 150,
    })


def hist_density(counts: List[float], spec: dict, rebin: int = 5) -> Tuple[np.ndarray, np.ndarray]:
    """Rebin counts (default x5: 0.2 ns -> 1.0 ns) and return area-normalized probability density."""
    arr = np.asarray(counts, dtype=np.float64)
    nb = len(arr)
    if rebin > 1 and nb % rebin == 0:
        arr = arr.reshape(-1, rebin).sum(axis=1)
        nb = len(arr)
    edges = np.linspace(spec["xmin"], spec["xmax"], nb + 1)
    bw = edges[1] - edges[0]
    total = arr.sum()
    density = arr / (total * bw) if total > 0 else arr
    return edges, density


def make_fig13_canvas(
    title: str,
    xlabel: str,
    xlim: Tuple[float, float] = (-60.0, 60.0),
    major_x: float = 20.0,
    minor_x_div: int = 4,
    logo_img=None,
):
    """Create a single-panel figure matching Electronic_2D_measurement_internal_note-6.pdf Fig. 13."""
    fig = plt.figure(figsize=(7.2, 4.65))
    ax = fig.add_axes([0.14, 0.145, 0.82, 0.72])

    # Crisp FASER logo above top-left corner, matching Fig. 13 placement
    if logo_img is not None:
        lax = fig.add_axes([0.070, 0.870, 0.165, 0.100])
        lax.imshow(logo_img)
        lax.axis("off")

    ax.set_title(title, fontsize=13.0, pad=8.0)
    ax.set_xlim(xlim[0], xlim[1])
    ax.xaxis.set_major_locator(MultipleLocator(major_x))
    ax.xaxis.set_minor_locator(AutoMinorLocator(minor_x_div))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))

    ax.set_xlabel(xlabel, fontsize=12.5, loc="right", labelpad=6.0)
    ax.set_ylabel("Normalized", fontsize=12.5, loc="top", labelpad=8.0)

    ax.grid(True, which="major", color="#d8d8d8", linestyle="-", linewidth=0.75, zorder=0)
    ax.set_axisbelow(True)

    return fig, ax


def finalize_and_save(
    fig,
    ax,
    max_y: float,
    out_dir: Path,
    stem: str,
    pdf: PdfPages,
    index_entries: list,
    desc: str,
    legend_loc: str = "lower left",
    legend_fontsize: float = 9.8,
    legend_ncol: int = 1,
    badge_x: float = 0.045,
    badge_y: float = 0.905,
):
    ymax = max(0.25, max_y * 1.22)
    ymin = -0.15 * ymax
    ax.set_ylim(ymin, ymax)

    # Set clean y major ticks (0.1 or 0.2 or 0.5) so 0.0 is always shown cleanly
    if ymax <= 0.45:
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
    elif ymax <= 0.95:
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
    elif ymax <= 1.8:
        ax.yaxis.set_major_locator(MultipleLocator(0.4))
    else:
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    # Hide negative y tick labels (matching Fig. 13 where ticks start at 0.0)
    fig.canvas.draw()
    yticks = [t for t in ax.get_yticks() if t >= -1e-9 and t <= ymax + 1e-9]
    ax.set_yticks(yticks)

    ax.text(
        badge_x, badge_y, "FASER Internal",
        transform=ax.transAxes,
        fontsize=12.2, fontweight="bold",
        va="top", ha="left",
    )

    leg = ax.legend(
        loc=legend_loc,
        frameon=True,
        fancybox=False,
        edgecolor="#c8c8c8",
        facecolor="white",
        framealpha=0.95,
        fontsize=legend_fontsize,
        ncol=legend_ncol,
        handlelength=1.8,
        borderpad=0.55,
        labelspacing=0.35,
    )
    leg.get_frame().set_linewidth(0.8)

    png_path = out_dir / f"{stem}.png"
    pdf_path = out_dir / f"{stem}.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    pdf.savefig(fig)
    plt.close(fig)
    index_entries.append({
        "stem": stem,
        "description": desc,
        "png": str(png_path),
        "pdf": str(pdf_path),
    })


def main():
    ap = argparse.ArgumentParser(description="Generate exact Fig. 13-style timing-difference figures")
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    args = ap.parse_args()

    configure_rc()
    logo_img = mpimg.imread(str(LOGO_PATH)) if LOGO_PATH.is_file() else None

    summary = json.loads((args.data_dir / "timing_summary.json").read_text())
    hdoc = json.loads((args.data_dir / "timing_histograms.json").read_text())
    specs = hdoc["hist_specs"]
    hists = hdoc["histograms"]

    out_dir = args.data_dir / "fig13_exact"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_pdf = out_dir / "fig13_all_timing_diffs_report.pdf"
    index_entries = []

    def plot_all_years(
        stem: str,
        hist_key: str,
        title: str,
        xlabel: str,
        desc: str,
        periods: List[str] = ALL_PERIODS,
        xlim: Tuple[float, float] = (-60.0, 60.0),
        major_x: float = 20.0,
        minor_x_div: int = 4,
        rebin: int = 5,
        legend_loc: str = "lower left",
        legend_fontsize: float = 9.5,
    ):
        fig, ax = make_fig13_canvas(title, xlabel, xlim=xlim, major_x=major_x,
                                    minor_x_div=minor_x_div, logo_img=logo_img)
        max_y = 0.0
        for p in periods:
            if p not in hists or hist_key not in hists[p]:
                continue
            lbl, col, ls, lw = PERIOD_CFG[p]
            edges, dens = hist_density(hists[p][hist_key], specs[hist_key], rebin=rebin)
            if dens.max() <= 0:
                continue
            # Restrict max_y to visible xlim
            centers = 0.5 * (edges[:-1] + edges[1:])
            vis = (centers >= xlim[0]) & (centers <= xlim[1])
            if np.any(vis):
                max_y = max(max_y, float(dens[vis].max()))
            ax.stairs(dens, edges, label=lbl, color=col, linestyle=ls, linewidth=lw)
        finalize_and_save(fig, ax, max_y, out_dir, stem, pdf, index_entries, desc,
                          legend_loc=legend_loc, legend_fontsize=legend_fontsize)

    with PdfPages(report_pdf) as pdf:
        # 1. Exact Fig. 13 reproduction & extension: Veto St2 - VetoNu (+16 ns corrected), ALL years + MC24
        plot_all_years(
            stem="fig13_01_veto2_vetonu_all_years",
            hist_key="veto2_vetonu_corr",
            title=r"All years — $t_{\mathrm{Veto}} - t_{\nu\mathrm{veto}}$ (+16 ns corrected)",
            xlabel=r"$t_{\mathrm{Veto}} - t_{\nu\mathrm{veto}}$ [ns]",
            desc="Veto St2 - VetoNu (+16 ns corrected) across 2022-2026 & MC24 (justifies VetoVetoNuTimeDifference: 33.165 / 35.783 / 27.21 ns)",
        )

        # 2. Veto St2 - VetoNu RAW (before +16 ns status==512 fix) -> justifies bit-512 +16 ns correction
        plot_all_years(
            stem="fig13_02_veto2_vetonu_raw_uncorrected",
            hist_key="veto2_vetonu_raw",
            title=r"All years — $t_{\mathrm{Veto}}^{\mathrm{raw}} - t_{\nu\mathrm{veto}}$ (before +16 ns fix)",
            xlabel=r"$t_{\mathrm{Veto}}^{\mathrm{raw}} - t_{\nu\mathrm{veto}}$ [ns]",
            desc="Veto St2 - VetoNu raw localtime across 2022-2026 & MC24 (shows CaloNu status==512 -16 ns satellite peak at 19.9 ns)",
        )

        # 3. Veto St2 - VetoNu zoomed [15, 45] ns (0.4 ns bins) -> clear separation of +2.7 ns CaloNu & MC24 shifts
        plot_all_years(
            stem="fig13_03_veto2_vetonu_zoom",
            hist_key="veto2_vetonu_corr",
            title=r"All years — $t_{\mathrm{Veto}} - t_{\nu\mathrm{veto}}$ (zoomed)",
            xlabel=r"$t_{\mathrm{Veto}} - t_{\nu\mathrm{veto}}$ [ns]",
            desc="Zoomed Veto St2 - VetoNu (+16 ns corrected) highlighting 33.1-33.2 ns (Dig1), 35.9 ns (CaloNu Dig2), and 27.3 ns (MC24)",
            xlim=(10.0, 44.0),
            major_x=5.0,
            minor_x_div=5,
            rebin=2,
            legend_loc="lower left",
            legend_fontsize=9.0,
        )

        # 4. Preshower (Dig1: ch 12,13) - VetoNu across ALL years + MC24 -> justifies PreshowerVetoNuTimeDifference
        plot_all_years(
            stem="fig13_04_preshower_dig1_vetonu_all_years",
            hist_key="ps_dig1_vetonu",
            title=r"All years — $t_{\mathrm{Preshower\,(Dig1)}} - t_{\nu\mathrm{veto}}$",
            xlabel=r"$t_{\mathrm{Preshower\,(Dig1)}} - t_{\nu\mathrm{veto}}$ [ns]",
            desc="Preshower Dig1 (ch 12,13) - VetoNu across 2022-2026 & MC24 (justifies PreshowerVetoNuTimeDifference: 53.625 / 53.512 / 13.27 / 51.41 ns)",
        )

        # 5. Preshower Dig2 (ch 24,25) vs Dig1 (ch 12,13) - VetoNu in 2025-2026 (raw vs +16 ns)
        fig, ax = make_fig13_canvas(
            r"2025–2026 — $t_{\mathrm{Preshower\,(Dig1\,vs\,Dig2)}} - t_{\nu\mathrm{veto}}$",
            r"$t_{\mathrm{Preshower}} - t_{\nu\mathrm{veto}}$ [ns]",
            xlim=(-60.0, 60.0),
            major_x=20.0,
            minor_x_div=4,
            logo_img=logo_img,
        )
        max_y = 0.0
        curves_ps2 = [
            ("2025", "ps_dig1_vetonu", "2025 Dig1 (ch 12,13)", "#1f77b4", "-", 1.7),
            ("2026", "ps_dig1_vetonu", "2026 Dig1 (ch 12,13)", "#17becf", "--", 1.6),
            ("2025", "ps_dig2_vetonu_raw", "2025 Dig2 (ch 24,25) raw", "#ff7f0e", ":", 1.6),
            ("2026", "ps_dig2_vetonu_raw", "2026 Dig2 (ch 24,25) raw", "#bcbd22", ":", 1.5),
            ("2025", "ps_dig2_vetonu_corr", "2025 Dig2 (ch 24,25) +16 ns", "#d62728", "-", 1.8),
            ("2026", "ps_dig2_vetonu_corr", "2026 Dig2 (ch 24,25) +16 ns", "#9467bd", "--", 1.7),
        ]
        for p, hk, lbl, col, ls, lw in curves_ps2:
            edges, dens = hist_density(hists[p][hk], specs[hk], rebin=5)
            max_y = max(max_y, float(dens.max()))
            ax.stairs(dens, edges, label=lbl, color=col, linestyle=ls, linewidth=lw)
        finalize_and_save(
            fig, ax, max_y, out_dir, "fig13_05_preshower_dig2_vs_dig1_vetonu_2025_2026", pdf, index_entries,
            "2025-2026 Preshower Dig1 vs Dig2 (raw and +16 ns corrected) relative to VetoNu",
            legend_loc="lower left", legend_fontsize=9.2,
        )

        # 6. Direct event-by-event Preshower Dig2 - Dig1 (raw and +16 ns) -> justifies PreshowerDig2ExtraOffset
        fig, ax = make_fig13_canvas(
            r"2025–2026 — $t_{\mathrm{Preshower\,(Dig2)}} - t_{\mathrm{Preshower\,(Dig1)}}$",
            r"$t_{\mathrm{Preshower\,(Dig2)}} - t_{\mathrm{Preshower\,(Dig1)}}$ [ns]",
            xlim=(-60.0, 60.0),
            major_x=20.0,
            minor_x_div=4,
            logo_img=logo_img,
        )
        max_y = 0.0
        for p, hk, lbl, col, ls, lw in [
            ("2025", "ps_dig2_minus_dig1_raw", "2025 Dig2 − Dig1 (raw)", "#ff7f0e", ":", 1.6),
            ("2026", "ps_dig2_minus_dig1_raw", "2026 Dig2 − Dig1 (raw)", "#17becf", ":", 1.5),
            ("2025", "ps_dig2_minus_dig1_corr", "2025 Dig2 − Dig1 (+16 ns)", "#8c564b", "-", 1.8),
            ("2026", "ps_dig2_minus_dig1_corr", "2026 Dig2 − Dig1 (+16 ns)", "#e377c2", "--", 1.8),
        ]:
            edges, dens = hist_density(hists[p][hk], specs[hk], rebin=5)
            max_y = max(max_y, float(dens.max()))
            ax.stairs(dens, edges, label=lbl, color=col, linestyle=ls, linewidth=lw)
        finalize_and_save(
            fig, ax, max_y, out_dir, "fig13_06_preshower_dig2_minus_dig1_2025_2026", pdf, index_entries,
            "Direct event-by-event Preshower Dig2 - Dig1 offset (justifies PreshowerDig2ExtraOffset: 3.97 ns in 2025, 3.64 ns in 2026)",
            legend_loc="upper right", legend_fontsize=9.4,
        )

        # 7. Combined 4-PMT Preshower (all PMTs) - VetoNu in 2025-2026: before vs after PreshowerDig2ExtraOffset (zoomed)
        fig, ax = make_fig13_canvas(
            r"2025–2026 — $t_{\mathrm{Preshower\,(all\,4\,PMTs)}} - t_{\nu\mathrm{veto}}$",
            r"$t_{\mathrm{Preshower\,(all\,4\,PMTs)}} - t_{\nu\mathrm{veto}}$ [ns]",
            xlim=(-10.0, 40.0),
            major_x=10.0,
            minor_x_div=5,
            logo_img=logo_img,
        )
        max_y = 0.0
        for p, hk, lbl, col, ls, lw in [
            ("2025", "ps_all_vetonu_raw", "2025 4-PMT (raw)", "#ff7f0e", ":", 1.6),
            ("2025", "ps_all_vetonu_corr16", "2025 4-PMT (+16 ns only)", "#d62728", "--", 1.7),
            ("2025", "ps_all_vetonu_aligned", "2025 4-PMT (+16 ns, Dig2 −3.97 ns)", "#2ca02c", "-", 1.9),
            ("2026", "ps_all_vetonu_corr16", "2026 4-PMT (+16 ns only)", "#9467bd", "-.", 1.6),
            ("2026", "ps_all_vetonu_aligned", "2026 4-PMT (+16 ns, Dig2 −3.64 ns)", "#1f77b4", "-", 1.8),
        ]:
            edges, dens = hist_density(hists[p][hk], specs[hk], rebin=2)
            centers = 0.5 * (edges[:-1] + edges[1:])
            vis = (centers >= -10.0) & (centers <= 40.0)
            max_y = max(max_y, float(dens[vis].max()))
            ax.stairs(dens, edges, label=lbl, color=col, linestyle=ls, linewidth=lw)
        finalize_and_save(
            fig, ax, max_y, out_dir, "fig13_07_preshower_4pmt_dig2_alignment_effect", pdf, index_entries,
            "Effect of PreshowerDig2ExtraOffset on combined 4-PMT Preshower mean in 2025 and 2026",
            legend_loc="upper right", legend_fontsize=8.5,
        )

        # 8. Timing station - VetoNu across ALL years + MC24 (control proving VetoNu never shifted!)
        plot_all_years(
            stem="fig13_08_timing_vetonu_all_years",
            hist_key="timing_vetonu",
            title=r"All years — $t_{\mathrm{Timing}} - t_{\nu\mathrm{veto}}$",
            xlabel=r"$t_{\mathrm{Timing}} - t_{\nu\mathrm{veto}}$ [ns]",
            desc="Timing - VetoNu across 2022-2026 & MC24 (flat at 34.7-35.0 ns in all data periods, proving VetoNu did not shift)",
        )

        # 9. Veto St1 - VetoNu across ALL years + MC24 (independent confirmation of CaloNu Dig2 +2.8 ns shift!)
        plot_all_years(
            stem="fig13_09_veto1_vetonu_all_years",
            hist_key="veto1_vetonu",
            title=r"All years — $t_{\mathrm{Veto1}} - t_{\nu\mathrm{veto}}$ (+16 ns corrected)",
            xlabel=r"$t_{\mathrm{Veto1}} - t_{\nu\mathrm{veto}}$ [ns]",
            desc="Veto St1 - VetoNu (+16 ns corrected) across 2022-2026 & MC24 (shows identical +2.8 ns CaloNu Dig2 shift)",
        )

        # 10. Closure check: Veto St2 - VetoNu AFTER applying Python timing presets (all years collapse to 0 ns!)
        plot_all_years(
            stem="fig13_10_veto2_after_preset_alignment",
            hist_key="veto2_vetonu_preset_resid",
            title=r"All years — $(t_{\mathrm{Veto}} - t_{\nu\mathrm{veto}}) - \Delta t_{\mathrm{preset}}$",
            xlabel=r"$(t_{\mathrm{Veto}} - t_{\nu\mathrm{veto}}) - \Delta t_{\mathrm{preset}}$ [ns]",
            desc="Post-alignment closure for Veto St2 across all 2022-2026 periods and MC24 (all peaks aligned to 0.0 ns)",
        )

        # 11. Closure check: Preshower (all PMTs) - VetoNu AFTER applying Python timing presets (all years collapse to 0 ns!)
        plot_all_years(
            stem="fig13_11_preshower_after_preset_alignment",
            hist_key="ps_all_vetonu_preset_resid",
            title=r"All years — $(t_{\mathrm{Preshower}} - t_{\nu\mathrm{veto}}) - \Delta t_{\mathrm{preset}}$",
            xlabel=r"$(t_{\mathrm{Preshower}} - t_{\nu\mathrm{veto}}) - \Delta t_{\mathrm{preset}}$ [ns]",
            desc="Post-alignment closure for Preshower (incl. 4-PMT Dig2 alignment) across all 2022-2026 periods and MC24 (all peaks aligned to 0.0 ns)",
        )

    (out_dir / "figure_index.json").write_text(json.dumps({
        "style_reference": "/afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/docs/Electronic_2D_measurement_internal_note-6.pdf (Fig. 13, p. 22)",
        "report_pdf": str(report_pdf),
        "figures": index_entries,
    }, indent=2))
    print(f"Generated {len(index_entries)} exact Fig. 13-style figures + {report_pdf}")


if __name__ == "__main__":
    main()
