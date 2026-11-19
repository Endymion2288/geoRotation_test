#!/usr/bin/env python3
"""Generate diagnostic figures for the 2025 Preshower GFIT_FAILED (p0013) FASER_DQ_Analyis audit."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np

AUDIT_DIR = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/"
    "stage04_waveform_reduced_charge/dq_audit_2025_p0013"
)
TIMING_HISTS = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/"
    "stage04_waveform_reduced_charge/timing_study_fig13/timing_histograms.json"
)
LOGO_PATH = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/"
    "stage04_waveform_reduced_charge/ref_pages/faser_logo_600dpi.png"
)

FOCAL_P13 = ["19415_p0013", "19524_p0013", "19866_p0013", "19923_p0013", "19970_p0013", "20081_p0013"]
OVERLAP_P13 = ["19878_p0013", "20012_p0013"]
OVERLAP_P14 = ["19336_p0014", "19878_p0014", "20012_p0014"]
REF_P14 = ["20170_p0014", "20279_p0014", "20481_p0014", "20591_p0014"]


def setup_rc():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "mathtext.default": "regular",
        "axes.linewidth": 0.95,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 5.5,
        "ytick.major.size": 5.5,
        "xtick.minor.size": 3.0,
        "ytick.minor.size": 3.0,
        "xtick.labelsize": 10.0,
        "ytick.labelsize": 10.0,
        "figure.dpi": 150,
    })


def sum_hists(fa: dict, keys: list[str], hkey: str) -> tuple[np.ndarray, np.ndarray]:
    counts = None
    xmin = xmax = None
    for k in keys:
        h = fa[k]["histograms"][hkey]
        arr = np.asarray(h["counts"], dtype=np.float64)
        counts = arr.copy() if counts is None else counts + arr
        xmin, xmax = h["xmin"], h["xmax"]
    edges = np.linspace(xmin, xmax, len(counts) + 1)
    return edges, counts


def to_density(edges: np.ndarray, counts: np.ndarray) -> np.ndarray:
    bw = edges[1] - edges[0]
    tot = counts.sum()
    return counts / (tot * bw) if tot > 0 else counts


def main():
    setup_rc()
    out_dir = AUDIT_DIR / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    full = json.loads((AUDIT_DIR / "dq_audit_2025_full.json").read_text())
    fa = full["full_audit"]
    surv = full["survey_2025"]
    th25 = full["thresholds_2025"]
    wf = json.loads((AUDIT_DIR / "preshower0_mean_waveform_20012.json").read_text())
    thdoc = json.loads(TIMING_HISTS.read_text())
    logo_img = mpimg.imread(str(LOGO_PATH)) if LOGO_PATH.is_file() else None

    report_pdf = out_dir / "dq_audit_2025_p0013_report.pdf"
    index = []

    def save(fig, stem: str, desc: str, pdf: PdfPages):
        png = out_dir / f"{stem}.png"
        pf = out_dir / f"{stem}.pdf"
        fig.savefig(png, dpi=300)
        fig.savefig(pf)
        pdf.savefig(fig)
        plt.close(fig)
        index.append({"stem": stem, "description": desc, "png": str(png), "pdf": str(pf)})

    with PdfPages(report_pdf) as pdf:
        # ------------------------------------------------------------------
        # FIG DQ1: ps_dig1_vetonu in 2025 (Linear vs Log-y) — p0013 bump proof
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.7))
        fig.subplots_adjust(left=0.075, right=0.98, top=0.84, bottom=0.14, wspace=0.22)
        if logo_img is not None:
            lax = fig.add_axes([0.045, 0.875, 0.115, 0.095])
            lax.imshow(logo_img)
            lax.axis("off")
        fig.suptitle(
            r"2025 $t_{\mathrm{Preshower\,(Dig1)}} - t_{\nu\mathrm{veto}}$: Impact of p0013 GFIT-failed runs (Linear vs Log)",
            fontsize=12.8, y=0.955,
        )

        # 2025 combined from fig13_04 (rebinned 0.2 -> 1.0 ns)
        c25_raw = np.asarray(thdoc["histograms"]["2025"]["ps_dig1_vetonu"], dtype=np.float64).reshape(-1, 5).sum(axis=1)
        e25 = np.linspace(-60.0, 60.0, len(c25_raw) + 1)
        d25_comb = to_density(e25, c25_raw)

        # Full-run p0014 only vs p0013 only (6 focal runs)
        e14, c14 = sum_hists(fa, OVERLAP_P14 + REF_P14, "dt_ps1_vn")
        e13, c13 = sum_hists(fa, FOCAL_P13, "dt_ps1_vn")
        d14 = to_density(e14, c14)
        d13 = to_density(e13, c13)
        # Also show p0013 normalized to total 2025 fig13_04 denominator (weighted contribution!)
        # Weight by actual fraction of 6 focal p0013 runs in timing_study_fig13 (95,161 / 6,837,434 = 1.39%)
        d13_weighted = d13 * (95161.0 / 6837434.0)

        for idx_ax, (ax, is_log) in enumerate(zip(axes, [False, True])):
            ax.grid(True, which="major", color="#dcdcdc", linewidth=0.75)
            ax.set_axisbelow(True)
            ax.stairs(d25_comb, e25, label="2025 combined (Fig. 13-04, N=6.84M)", color="#8c564b", lw=2.0)
            ax.stairs(d14, e14, label="2025 p0014 only (7 runs, N=39.3M)", color="#1f77b4", ls="--", lw=1.7)
            if not is_log:
                ax.stairs(d13, e13, label="2025 p0013 only (6 runs, unit-norm, N=0.70M)", color="#d62728", ls="-", lw=1.8)
                ax.set_ylim(-0.04, 0.54)
                ax.set_title("Linear scale (unit-area normalized)", fontsize=11.0, pad=6)
                ax.text(0.04, 0.92, "FASER Internal", transform=ax.transAxes, fontsize=11.0, fontweight="bold", va="top")
                ax.annotate(
                    "p0013 surviving ~4%\nclipped at t >= 178 ns\npeak at 30–34 ns",
                    xy=(32.5, 0.088), xytext=(18.0, 0.24),
                    arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.2),
                    fontsize=8.8, color="#b01e1e",
                )
                ax.legend(loc="lower left", fontsize=8.3, frameon=True, fancybox=False, edgecolor="#cccccc", framealpha=0.95)
            else:
                ax.stairs(d13_weighted, e13, label="p0013 weighted share in Fig. 13-04 (1.39%)", color="#d62728", ls="-", lw=1.8)
                ax.set_yscale("log")
                ax.set_ylim(2e-6, 1.2)
                ax.set_title("Log scale (exposes 1.30% p0013 bump at 22–42 ns)", fontsize=11.0, pad=6)
                ax.text(0.04, 0.92, "FASER Internal", transform=ax.transAxes, fontsize=11.0, fontweight="bold", va="top")
                ax.annotate(
                    "88,895 events (1.30%)\nin [22, 42] ns from\n6 p0013 runs\n(density ~10⁻³)",
                    xy=(32.0, 1.0e-3), xytext=(22.0, 1.8e-2),
                    arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.2),
                    fontsize=8.6, color="#b01e1e",
                )
                ax.legend(loc="lower left", fontsize=8.3, frameon=True, fancybox=False, edgecolor="#cccccc", framealpha=0.95)
            ax.set_xlim(-60, 60)
            ax.xaxis.set_major_locator(MultipleLocator(20))
            ax.xaxis.set_minor_locator(AutoMinorLocator(4))
            ax.set_xlabel(r"$t_{\mathrm{Preshower\,(Dig1)}} - t_{\nu\mathrm{veto}}$ [ns]", fontsize=11.0, loc="right")
            ax.set_ylabel("Normalized", fontsize=11.0, loc="top")
        save(fig, "figDQ1_ps_dig1_vetonu_p0013_contribution",
             "2025 Preshower(Dig1) - VetoNu linear vs log showing 1.30% secondary bump from p0013 runs at 30-34 ns", pdf)

        # ------------------------------------------------------------------
        # FIG DQ2: Root-cause mechanism (2x2)
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(2, 2, figsize=(12.6, 8.8))
        fig.subplots_adjust(left=0.075, right=0.975, top=0.90, bottom=0.075, wspace=0.23, hspace=0.30)
        fig.suptitle(
            "Root Cause of 2025 p0013 Preshower0/1 GFIT_FAILED: Stale COOL TimingOffset (+18 ns -> Fit Window [178, 378] ns)",
            fontsize=12.4, fontweight="bold", y=0.965,
        )

        # (a) Real 20012 p0014 Preshower0 waveform + p0013 vs p0014 windows
        ax = axes[0, 0]
        ax.grid(True, color="#e0e0e0", lw=0.7)
        t_wf = np.asarray(wf["time_ns"])
        w_wf = np.asarray(wf["mean_mv"])
        ax.axvspan(120, 178, color="#1f77b4", alpha=0.12, label="Missed by p0013 ([120, 178] ns)")
        ax.axvspan(178, 260, color="#d62728", alpha=0.10, label="p0013 fit window starts at 178 ns")
        ax.axvline(120, color="#1f77b4", ls="--", lw=1.4, label="p0014 lo_edge = 120 ns (Offset = −40 ns)")
        ax.axvline(178, color="#d62728", ls="-", lw=1.8, label="p0013 lo_edge = 178 ns (Offset = +18 ns)")
        ax.step(t_wf, w_wf, where="mid", color="#111111", lw=1.8, label="Real 20012 p0014 ch=12 waveform")
        ax.axvspan(158, 169, color="#2ca02c", alpha=0.18, label="Typical 2025 PS0 localtime (158–169 ns)")
        ax.set_xlim(110, 250)
        ax.set_ylim(-1.0, 15.5)
        ax.set_xlabel("Digitizer sample time [ns]", fontsize=10.5)
        ax.set_ylabel("Baseline-subtracted amplitude [mV]", fontsize=10.5)
        ax.set_title("(a) Fit-window truncation on real 2025 Preshower0 pulse", fontsize=10.8)
        ax.legend(loc="upper right", fontsize=7.8, frameon=True, framealpha=0.94)

        # (b) Preshower0_charge (exact FASER_DQ_Analyis Fig. 5j)
        ax = axes[0, 1]
        ax.grid(True, color="#e0e0e0", lw=0.7)
        for k, col, ls in [
            ("19415_p0013", "#d62728", "-"),
            ("19524_p0013", "#ff7f0e", "-"),
            ("19923_p0013", "#e377c2", "--"),
            ("20012_p0013", "#8c564b", ":"),
            ("20012_p0014", "#1f77b4", "-"),
            ("19878_p0014", "#2ca02c", "--"),
            ("20170_p0014", "#17becf", ":"),
        ]:
            h = fa[k]["histograms"]["Preshower0_charge"]
            ed = np.linspace(h["xmin"], h["xmax"], h["nbins"] + 1)
            dn = to_density(ed, np.asarray(h["counts"], dtype=float))
            ax.stairs(dn, ed, label=k, color=col, ls=ls, lw=1.6)
        ax.set_xlim(0, 15)
        ax.set_xlabel("Preshower0_charge [pC] (FASER_DQ_Analyis Fig. 5j)", fontsize=10.5)
        ax.set_ylabel("Probability density", fontsize=10.5)
        ax.set_title("(b) Preshower0 charge: p0013 collapsed vs p0014 MIP peak", fontsize=10.8)
        ax.legend(loc="upper right", fontsize=8.0, frameon=True, ncol=2)

        # (c) Preshower0_localtime_wide (showing 0 ns GFIT_FAILED pileup & >=178 ns clipped edge)
        ax = axes[1, 0]
        ax.grid(True, color="#e0e0e0", lw=0.7)
        for k, col, ls in [
            ("20012_p0013", "#d62728", "-"),
            ("19415_p0013", "#ff7f0e", "--"),
            ("20012_p0014", "#1f77b4", "-"),
            ("20170_p0014", "#2ca02c", "--"),
        ]:
            h = fa[k]["histograms"]["Preshower0_localtime_wide"]
            ed = np.linspace(h["xmin"], h["xmax"], h["nbins"] + 1)
            dn = to_density(ed, np.asarray(h["counts"], dtype=float))
            ax.stairs(dn, ed, label=k, color=col, ls=ls, lw=1.7)
        ax.axvline(178, color="#d62728", ls=":", lw=1.3, label="p0013 window edge (178 ns)")
        ax.set_yscale("log")
        ax.set_xlim(-10, 270)
        ax.set_ylim(1e-5, 1.5)
        ax.set_xlabel("Preshower0_localtime [ns] (raw_charge > 2.5 pC)", fontsize=10.5)
        ax.set_ylabel("Probability density (log)", fontsize=10.5)
        ax.set_title("(c) Preshower0 localtime: GFIT fail (0 ns) & 178 ns clip", fontsize=10.8)
        ax.legend(loc="upper left", bbox_to_anchor=(0.06, 0.98), fontsize=8.0, frameon=True, framealpha=0.95)

        # (d) All-2025 survey: Preshower0 GFIT_FAILED % vs run number
        ax = axes[1, 1]
        ax.grid(True, color="#e0e0e0", lw=0.7)
        p13_s = [r for r in surv if r["status"] == "ok" and r["tag"] == "p0013" and r["n_goodtrk"] > 0]
        p14_s = [r for r in surv if r["status"] == "ok" and r["tag"] == "p0014" and r["n_goodtrk"] > 0]
        ax.scatter([r["run"] for r in p13_s], [100 * r["ps0_gfit_frac"] for r in p13_s],
                   s=22, color="#d62728", label=f"p0013 ({len(p13_s)}/67 runs, Offset=+18 ns)")
        ax.scatter([r["run"] for r in p14_s], [100 * r["ps0_gfit_frac"] for r in p14_s],
                   s=22, color="#1f77b4", marker="s", label=f"p0014 ({len(p14_s)}/94 runs, Offset=−40 ns)")
        # Highlight 6 focal runs
        focal_s = [r for r in p13_s if r["run"] in {19415, 19524, 19866, 19923, 19970, 20081}]
        ax.scatter([r["run"] for r in focal_s], [100 * r["ps0_gfit_frac"] for r in focal_s],
                   s=65, facecolors="none", edgecolors="#111111", lw=1.4, label="6 sampled focal p0013 runs")
        ax.set_ylim(-5, 105)
        ax.set_xlabel("2025 Run number", fontsize=10.5)
        ax.set_ylabel("Preshower0 GFIT_FAILED [%]", fontsize=10.5)
        ax.set_title("(d) All-2025 survey: 100% of p0013 runs affected (67/67)", fontsize=10.8)
        ax.legend(loc="center right", fontsize=8.4, frameon=True)
        save(fig, "figDQ2_p0013_preshower_window_bug_and_all2025",
             "Root cause of p0013 Preshower0/1 GFIT_FAILED: fit window truncation, charge collapse, localtime clipping, and all-2025 survey", pdf)

        # ------------------------------------------------------------------
        # FIG DQ3: Full FASER_DQ_Analyis multi-subsystem audit (2x2)
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.0))
        fig.subplots_adjust(left=0.075, right=0.975, top=0.90, bottom=0.125, wspace=0.22, hspace=0.35)
        fig.suptitle(
            "Complete FASER_DQ_Analyis Multi-Subsystem Audit: Focal p0013 Runs vs Overlap A/B (19878, 20012) & p0014 Reference",
            fontsize=12.4, fontweight="bold", y=0.965,
        )

        # (a) Status OK% across all 16 channels
        ax = axes[0, 0]
        ax.grid(True, axis="y", color="#e0e0e0", lw=0.7)
        chs = ["VetoNu0", "VetoNu1", "Veto10", "Veto11", "Veto20", "Veto21",
               "Timing0", "Timing1", "Timing2", "Timing3",
               "Preshower0", "Preshower1", "Preshower2", "Preshower3", "CaloLo0", "CaloHi0"]
        short = ["VN0", "VN1", "V10", "V11", "V20", "V21", "Tm0", "Tm1", "Tm2", "Tm3", "PS0", "PS1", "PS2", "PS3", "CL0", "CH0"]
        x = np.arange(len(chs))
        y_foc = [np.mean([fa[k]["status_summary"][c]["ok_pct"] for k in FOCAL_P13]) for c in chs]
        y_20012_13 = [fa["20012_p0013"]["status_summary"][c]["ok_pct"] for c in chs]
        y_20012_14 = [fa["20012_p0014"]["status_summary"][c]["ok_pct"] for c in chs]
        y_ref = [np.mean([fa[k]["status_summary"][c]["ok_pct"] for k in REF_P14]) for c in chs]
        w = 0.20
        ax.bar(x - 1.5 * w, y_foc, width=w, color="#d62728", label="6 focal p0013 mean")
        ax.bar(x - 0.5 * w, y_20012_13, width=w, color="#ff7f0e", label="20012 p0013")
        ax.bar(x + 0.5 * w, y_20012_14, width=w, color="#1f77b4", label="20012 p0014")
        ax.bar(x + 1.5 * w, y_ref, width=w, color="#2ca02c", label="4 ref p0014 mean")
        ax.set_xticks(x)
        ax.set_xticklabels(short, fontsize=8.5, rotation=30)
        ax.set_ylim(0, 124)
        ax.set_ylabel("Good hit fraction (status == 0/512) [%]", fontsize=10.0)
        ax.set_title("(a) All 16 channels status OK%: ONLY PS0/PS1 broken", fontsize=10.6)
        ax.legend(loc="upper right", ncol=2, fontsize=7.8, frameon=True)

        # (b) Veto20_charge (Fig. 5d) — proving non-PS scintillators identical
        ax = axes[0, 1]
        ax.grid(True, color="#e0e0e0", lw=0.7)
        for k, col, ls in [
            ("19415_p0013", "#d62728", "-"),
            ("19524_p0013", "#ff7f0e", "-"),
            ("19923_p0013", "#e377c2", "--"),
            ("20012_p0013", "#8c564b", ":"),
            ("20012_p0014", "#1f77b4", "-"),
            ("20170_p0014", "#2ca02c", "--"),
        ]:
            h = fa[k]["histograms"]["Veto20_charge"]
            ed = np.linspace(h["xmin"], h["xmax"], h["nbins"] + 1)
            dn = to_density(ed, np.asarray(h["counts"], dtype=float))
            ax.stairs(dn, ed, label=k, color=col, ls=ls, lw=1.6)
        ax.set_xlim(0, 300)
        ax.set_xlabel("Veto20_charge [pC] (FASER_DQ_Analyis Fig. 5d)", fontsize=10.2)
        ax.set_ylabel("Probability density", fontsize=10.2)
        ax.set_title("(b) Veto20 charge: 100% consistent between p0013 & p0014", fontsize=10.6)
        ax.legend(loc="upper right", fontsize=8.0, frameon=True, ncol=2)

        # (c) Track_Chi2 (Fig. 6a) — showing TRACKER-ALIGN-04 -> 05 shift!
        ax = axes[1, 0]
        ax.grid(True, color="#e0e0e0", lw=0.7)
        for k, col, ls in [
            ("19415_p0013", "#d62728", "-"),
            ("19923_p0013", "#ff7f0e", "--"),
            ("19878_p0013", "#8c564b", "-."),
            ("20012_p0013", "#e377c2", ":"),
            ("19878_p0014", "#1f77b4", "-"),
            ("20012_p0014", "#2ca02c", "--"),
            ("20170_p0014", "#17becf", ":"),
        ]:
            h = fa[k]["histograms"]["Track_Chi2"]
            ed = np.linspace(h["xmin"], h["xmax"], h["nbins"] + 1)
            dn = to_density(ed, np.asarray(h["counts"], dtype=float))
            ax.stairs(dn, ed, label=f"{k} (μ={h['mean']:.2f})", color=col, ls=ls, lw=1.6)
        ax.set_xlim(0, 50)
        ax.set_xlabel(r"Track $\chi^2$ (FASER_DQ_Analyis Fig. 6a)", fontsize=10.2)
        ax.set_ylabel("Probability density", fontsize=10.2)
        ax.set_title(r"(c) Track $\chi^2$: subtle TRACKER-ALIGN-04 $\to$ 05 (−0.74) shift", fontsize=10.6)
        ax.legend(loc="upper right", fontsize=7.7, frameon=True)

        # (d) GoodTracks/pb^-1 vs FASER_DQ_Analyis thresholds.json
        ax = axes[1, 1]
        ax.grid(True, axis="y", color="#e0e0e0", lw=0.7)
        ds_order = FOCAL_P13 + ["19336_p0013", "19878_p0013", "20012_p0013"] + OVERLAP_P14 + REF_P14
        vals = [fa[k]["rates_per_pb"]["good_tracks"] for k in ds_order]
        cols_bar = []
        for k in ds_order:
            fl = fa[k]["dq_yield_flags"]["track_flag"]
            if fl != "normal":
                cols_bar.append("#ff7f0e")
            elif k.endswith("p0013"):
                cols_bar.append("#d62728")
            else:
                cols_bar.append("#1f77b4")
        xi = np.arange(len(ds_order))
        ax.axhspan(th25["tracks"]["low"], th25["tracks"]["high"], color="#2ca02c", alpha=0.12,
                   label=f"DQ normal band [4641, 10330] (m±5s)")
        ax.axhline(th25["tracks"]["median"], color="#2ca02c", ls="--", lw=1.3, label=f"2025 median ({th25['tracks']['median']:.0f})")
        ax.bar(xi, vals, color=cols_bar, width=0.68, edgecolor="#333333", lw=0.5)
        ax.set_xticks(xi)
        ax.set_xticklabels([k.replace("_p00", " p") for k in ds_order], fontsize=7.8, rotation=42, ha="right")
        ax.set_ylim(0, 11800)
        ax.set_ylabel(r"GoodTracks / pb$^{-1}$", fontsize=10.2)
        ax.set_title("(d) Track yields vs DQ m±5s band (flags 19970 & 19336 p13)", fontsize=10.6)
        ax.legend(loc="upper right", fontsize=8.0, frameon=True)
        save(fig, "figDQ3_faser_dq_all_subsystems_check",
             "Full FASER_DQ_Analyis multi-subsystem check across status OK%, Veto20 charge, Track Chi2 alignment shift, and GoodTracks/pb^-1", pdf)

    (out_dir / "figure_index.json").write_text(json.dumps({
        "report_pdf": str(report_pdf),
        "figures": index,
    }, indent=2))
    print(f"Generated {len(index)} DQ audit figures + {report_pdf}")


if __name__ == "__main__":
    main()
