#!/usr/bin/env python3
"""Generate English validation figures from validation_summary_v2.json.

Produces PNG, per-plot PDF, combined multi-page PDF report, and figure_index.json for:
  1. 2022-2026 real-data validation: CaloNu status==512 recovery, 2025-2026 Preshower fallback
     recovery, Ben Wilson 5-folder parity, and elimination of OOB clockPhase corruption.
  2. All 9 official default-flux MC24 GENIE neutrino runs (FASER_FaserMC.pdf):
     2022/23 Crossing angle, 2024 FaserNu geometry, and 2024 CaloNu geometry.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

DEFAULT_DIR = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/"
    "stage04_waveform_reduced_charge/validation_v2"
)


def setup_style():
    plt.rcParams.update({
        "font.size": 10.5,
        "axes.labelsize": 11.2,
        "axes.titlesize": 11.8,
        "legend.fontsize": 8.6,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linestyle": "--",
        "figure.dpi": 150,
    })


def save_fig(fig, fig_dir: Path, stem: str, pdf: PdfPages, index_entries: list, title: str):
    png = fig_dir / f"{stem}.png"
    pdf_p = fig_dir / f"{stem}.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf_p, bbox_inches="tight")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    index_entries.append({"stem": stem, "title": title, "png": str(png), "pdf": str(pdf_p)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DIR)
    args = ap.parse_args()

    setup_style()
    doc = json.loads((args.data_dir / "validation_summary_v2.json").read_text())
    samples = doc["samples"]
    fig_dir = args.data_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    report_pdf = fig_dir / "validation_v2_report.pdf"
    index_entries = []
    edges_q = np.linspace(0.0, 250.0, 51)

    with PdfPages(report_pdf) as pdf:
        # ------------------------------------------------------------------
        # Fig V1: 2024 CaloNu status==512 & 2025/2026 Preshower fallback recovery
        # ------------------------------------------------------------------
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16.2, 5.1))

        # Panel 1: 2024 CaloNu Veto2 status==512
        s_cn = samples.get("2024_CaloNu_015849", {}).get("fix_vs_base", {})
        h_cn = s_cn.get("hists", {})
        if "v2_512_fix" in h_cn:
            for key, lbl, col, ls in [
                ("v2_512_raw", "Raw VetoNu0 charge", "#555555", ":"),
                ("v2_512_base", "Baseline NtupleDumper (broken)", "#d62728", "--"),
                ("v2_512_fix", "Refactored (--timingPreset data_caloNu)", "#2ca02c", "-"),
            ]:
                y = np.asarray(h_cn[key], dtype=float)
                if y.sum() > 0:
                    y = y / y.sum()
                ax1.stairs(y, edges_q, label=lbl, color=col, ls=ls, lw=2.0)
            st = s_cn["mu_v2_status512"]
            ax1.set_title(
                f"(a) 2024 CaloNu (15849): Veto2 status==512\n"
                f"N={st['n']}, med: fix={st['fix_med']:.1f} vs base={st['base_med']:.1f} pC"
            )
        ax1.axvline(30.0, color="#000000", ls=":", lw=1.1, label="30 pC veto threshold")
        ax1.set_xlim(0.0, 220.0)
        ax1.set_ylim(0.0, ax1.get_ylim()[1] * 1.25)
        ax1.set_xlabel("VetoNu0 charge [pC]")
        ax1.set_ylabel("Normalized (unit area / 5 pC bin)")
        ax1.legend(loc="upper right", fontsize=8.0, framealpha=0.92)

        # Panel 2: 2025 Preshower fallback (sum 20012 + 20481)
        y_raw25 = np.zeros(50)
        y_fix25 = np.zeros(50)
        y_bas25 = np.zeros(50)
        n25 = 0
        for nm in ("2025_20012", "2025_20481"):
            fb = samples.get(nm, {}).get("fix_vs_base", {})
            h = fb.get("hists", {})
            if "ps_fb_fix" in h:
                y_raw25 += np.asarray(h["ps_fb_raw"])
                y_fix25 += np.asarray(h["ps_fb_fix"])
                y_bas25 += np.asarray(h["ps_fb_base"])
                n25 += fb["mu_preshower_fallback"]["n"]
        if n25 > 0:
            for y, lbl, col, ls in [
                (y_raw25, "Raw VetoNu0 charge", "#555555", ":"),
                (y_bas25, "Baseline (53.625 ns -> truncated!)", "#d62728", "--"),
                (y_fix25, "Refactored (--timingPreset data_2025)", "#1f77b4", "-"),
            ]:
                ax2.stairs(y / max(1.0, y.sum()), edges_q, label=lbl, color=col, ls=ls, lw=2.0)
            ax2.set_title(f"(b) 2025 Data (20012+20481): Preshower fallback\nN={n25} muons (Veto2 saturated/invalid)")
        ax2.axvline(30.0, color="#000000", ls=":", lw=1.1, label="30 pC veto threshold")
        ax2.set_xlim(0.0, 220.0)
        ax2.set_ylim(0.0, ax2.get_ylim()[1] * 1.25)
        ax2.set_xlabel("VetoNu0 charge [pC]")
        ax2.set_ylabel("Normalized (unit area / 5 pC bin)")
        ax2.legend(loc="upper right", fontsize=8.0, framealpha=0.92)

        # Panel 3: 2026 Preshower fallback (sum 22404 + 23176)
        y_raw26 = np.zeros(50)
        y_fix26 = np.zeros(50)
        y_bas26 = np.zeros(50)
        n26 = 0
        for nm in ("2026_22404", "2026_23176"):
            fb = samples.get(nm, {}).get("fix_vs_base", {})
            h = fb.get("hists", {})
            if "ps_fb_fix" in h:
                y_raw26 += np.asarray(h["ps_fb_raw"])
                y_fix26 += np.asarray(h["ps_fb_fix"])
                y_bas26 += np.asarray(h["ps_fb_base"])
                n26 += fb["mu_preshower_fallback"]["n"]
        if n26 > 0:
            for y, lbl, col, ls in [
                (y_raw26, "Raw VetoNu0 charge", "#555555", ":"),
                (y_bas26, "Baseline (53.625 ns -> truncated!)", "#d62728", "--"),
                (y_fix26, "Refactored (--timingPreset data_2026)", "#9467bd", "-"),
            ]:
                ax3.stairs(y / max(1.0, y.sum()), edges_q, label=lbl, color=col, ls=ls, lw=2.0)
            ax3.set_title(f"(c) 2026 Data (22404+23176): Preshower fallback\nN={n26} muons (Veto2 saturated/invalid)")
        ax3.axvline(30.0, color="#000000", ls=":", lw=1.1, label="30 pC veto threshold")
        ax3.set_xlim(0.0, 220.0)
        ax3.set_ylim(0.0, ax3.get_ylim()[1] * 1.25)
        ax3.set_xlabel("VetoNu0 charge [pC]")
        ax3.set_ylabel("Normalized (unit area / 5 pC bin)")
        ax3.legend(loc="upper right", fontsize=8.0, framealpha=0.92)

        fig.suptitle(
            "End-to-End Compiled NtupleDumper Validation: 2024 CaloNu & 2025-2026 Preshower Fallback Recovery",
            fontsize=13.0, fontweight="bold", y=1.02,
        )
        fig.tight_layout()
        save_fig(fig, fig_dir, "figV1_data_caloNu_and_2025_2026_recovery", pdf, index_entries,
                 "2024 CaloNu status==512 and 2025-2026 Preshower fallback recovery")

        # ------------------------------------------------------------------
        # Fig V2: All 9 Official Default-Flux MC24 GENIE Neutrino Runs
        # ------------------------------------------------------------------
        groups = [
            ("2022/23 Crossing Angle\n(FaserNu, s0014)", "MC24_Genie_22_23",
             ["MC24_Genie2223_200082_epos_pi", "MC24_Genie2223_200083_epos_k", "MC24_Genie2223_200089_charm"],
             "#1f77b4"),
            ("2024 Crossing Angle\n(FaserNu geom, s0014)", "MC24_Genie_24_FaserNu",
             ["MC24_Genie24FN_200093_epos_pi", "MC24_Genie24FN_200094_epos_k", "MC24_Genie24FN_200100_charm"],
             "#2ca02c"),
            ("2024 Crossing Angle\n(CaloNu geom, s0015)", "MC24_Genie_24_CaloNu",
             ["MC24_Genie24CN_200139_epos_pi", "MC24_Genie24CN_200140_epos_k", "MC24_Genie24CN_200146_charm"],
             "#d62728"),
        ]

        fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))
        for ax, (gtitle, gkey, rnames, col) in zip(axes, groups):
            h_raw = np.zeros(50)
            h_red = np.zeros(50)
            n_sel = 0
            n_raw_fail = 0
            n_rescued = 0
            fb_tot = {0: 0, 1: 0, 2: 0}
            for rn in rnames:
                fb = samples.get(rn, {}).get("fix_vs_base", {})
                gs = fb.get("genie_stats")
                if not gs:
                    continue
                h_raw += np.asarray(fb["hists"]["nu_raw_max"])
                h_red += np.asarray(fb["hists"]["nu_fix_red_max"])
                n_sel += gs["n_nu_downstream_sel"]
                n_raw_fail += gs["n_raw_ge40_backsplash"]
                n_rescued += gs["n_rescued_by_reduced_lt30"]
                for k, v in gs["fb_counts_nu_sel"].items():
                    fb_tot[int(k)] = fb_tot.get(int(k), 0) + int(v)
            if h_raw.sum() > 0:
                ax.stairs(h_raw / h_raw.sum(), edges_q, label="max(VetoNu0,1) raw charge", color="#555555", ls="--", lw=1.9)
                ax.stairs(h_red / h_red.sum(), edges_q, label="max(VetoNu0,1) reduced charge", color=col, ls="-", lw=2.1)
            ax.axvline(30.0, color="#000000", ls=":", lw=1.1, label="30 pC reduced threshold")
            ax.axvline(40.0, color="#888888", ls="-.", lw=1.0, label="40 pC raw threshold")
            ax.set_yscale("log")
            ax.set_xlim(0.0, 240.0)
            ax.set_ylim(1e-3, 1.8)
            ax.set_xlabel(r"$\max(Q_{\mathrm{VetoNu0}}, Q_{\mathrm{VetoNu1}})$ [pC]")
            ax.set_ylabel("Normalized (log scale / 5 pC bin)")
            rescue_pct = 100.0 * n_rescued / max(1, n_raw_fail)
            fb2_pct = 100.0 * fb_tot[2] / max(1, n_sel)
            ax.set_title(
                f"{gtitle}\n"
                f"N_sel={n_sel:,}, Veto2->PS fallback={fb2_pct:.1f}%\n"
                f"Backsplash rescued: {n_rescued}/{n_raw_fail} ({rescue_pct:.1f}%)"
            )
            ax.legend(loc="upper right", fontsize=7.8, framealpha=0.92)

        fig.suptitle(
            "Official Default-Flux MC24 GENIE Neutrino Validation (9 Runs from FASER_FaserMC.pdf, --timingPreset mc)",
            fontsize=13.0, fontweight="bold", y=1.03,
        )
        fig.tight_layout()
        save_fig(fig, fig_dir, "figV2_genie_default_flux_9runs_validation", pdf, index_entries,
                 "Official default-flux MC24 GENIE 9-run reduced-charge validation")

        # ------------------------------------------------------------------
        # Fig V3: Per-run summary bar chart (9 GENIE runs fallback + backsplash rescue + OOB bugs fixed)
        # ------------------------------------------------------------------
        genie_keys = [k for k, v in samples.items() if v.get("is_genie") and "fix_vs_base" in v]
        if genie_keys:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 5.2))
            x = np.arange(len(genie_keys))
            short_g = [k.replace("MC24_Genie", "").replace("_", "\n") for k in genie_keys]
            fb1_pcts, fb2_pcts, rescue_pcts, v2ovf_pcts = [], [], [], []
            for k in genie_keys:
                fb = samples[k]["fix_vs_base"]
                gs = fb["genie_stats"]
                n_s = max(1, gs["n_nu_downstream_sel"])
                fb1_pcts.append(100.0 * gs["fb_counts_nu_sel"].get("1", gs["fb_counts_nu_sel"].get(1, 0)) / n_s)
                fb2_pcts.append(100.0 * gs["fb_counts_nu_sel"].get("2", gs["fb_counts_nu_sel"].get(2, 0)) / n_s)
                rescue_pcts.append(100.0 * (gs["rescue_frac"] or 0.0))
                v2ovf_pcts.append(100.0 * (gs["v2_overflow_frac"] or 0.0))

            ax1.bar(x, fb1_pcts, color="#1f77b4", label="Fallback 1: Veto St2 timing [%]")
            ax1.bar(x, fb2_pcts, bottom=fb1_pcts, color="#ff7f0e", label="Fallback 2: Preshower timing (Veto2 sat.) [%]")
            ax1.set_xticks(x)
            ax1.set_xticklabels(short_g, fontsize=7.8)
            ax1.set_ylim(0, 122)
            ax1.set_ylabel("Share of selected neutrino events [%]")
            ax1.set_title("(a) Timing Reference Source on Selected GENIE Events")
            ax1.legend(loc="upper right", fontsize=8.2)

            ax2.bar(x - 0.19, v2ovf_pcts, width=0.38, color="#d62728", alpha=0.85, label="Veto2 overflow rate (status & 4) [%]")
            ax2.bar(x + 0.19, rescue_pcts, width=0.38, color="#2ca02c", alpha=0.90, label="Backsplash rescued (raw>=40 -> red<30) [%]")
            ax2.set_xticks(x)
            ax2.set_xticklabels(short_g, fontsize=7.8)
            ax2.set_ylim(0, 92)
            ax2.set_ylabel("Fraction [%]")
            ax2.set_title("(b) Veto2 Saturation Rate vs Backsplash Rescue Rate")
            ax2.legend(loc="upper left", fontsize=8.2)

            fig.suptitle(
                "Per-Run Validation Across All 9 Official Default-Flux MC24 GENIE Runs (3,000 events/run)",
                fontsize=13.0, fontweight="bold", y=1.01,
            )
            fig.tight_layout()
            save_fig(fig, fig_dir, "figV3_genie_9runs_fallback_and_oob_summary", pdf, index_entries,
                     "Per-run GENIE 9-run fallback mode, Veto2 overflow, and backsplash rescue summary")

    (fig_dir / "figure_index.json").write_text(json.dumps({
        "report_pdf": str(report_pdf),
        "figures": index_entries,
    }, indent=2))
    print(f"Generated {len(index_entries)} validation figures + {report_pdf}")


if __name__ == "__main__":
    main()
