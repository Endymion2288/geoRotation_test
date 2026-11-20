#!/usr/bin/env python3
"""Plot exact Fig. 12-style multi-station waveforms across all 2022-2026 data periods
and MC24 (modeled on Electronic_2D_measurement_internal_note-6.pdf, p. 21 Fig. 12),
plus quantitative integration-window optimization & pulse-shape comparison figures.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator, FixedLocator
import numpy as np

STUDY_DIR = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/"
    "stage04_waveform_reduced_charge/fig12_waveform_window_study"
)


def fmt_ns(val: float) -> str:
    return f"{val:.3f}".rstrip("0").rstrip(".")


def _clean_upper_yticks(ax, ymax: float):
    """Keep 0 and clean major y-ticks well below ymax so stacked hspace=0 panels never collide."""
    raw = [float(t) for t in ax.get_yticks() if 0.0 <= t <= 0.87 * ymax]
    if len(raw) >= 2:
        ax.yaxis.set_major_locator(FixedLocator(raw))
        step = raw[1] - raw[0]
        if step > 0:
            minors = [t + 0.5 * step for t in raw if t + 0.5 * step < ymax]
            ax.yaxis.set_minor_locator(FixedLocator(minors))
    else:
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def configure_rc():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "mathtext.default": "regular",
        "axes.linewidth": 0.9,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 4.5,
        "ytick.major.size": 4.5,
        "xtick.minor.size": 2.5,
        "ytick.minor.size": 2.5,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "figure.dpi": 150,
    })


def draw_fig12_column(
    axes: list,
    ev: dict[str, Any],
    dv: float,
    dps: float,
    ddig2: float,
    col_title: str,
    show_512_fix: bool = False,
    is_neutrino: bool = False,
    compact: bool = False,
):
    """Draw one 3-panel stacked Fig. 12 waveform column (Top: VetoNu, Mid: Veto2, Bot: Preshower)."""
    ax_vn, ax_v2, ax_ps = axes
    vn = ev["vn"]
    v2 = ev["v2"]
    ps1 = ev["ps1"]
    ps2 = ev.get("ps2")

    tv_vn = np.asarray(vn["tv"], dtype=float)
    wv_vn = np.asarray(vn["wv"], dtype=float)
    tv_v2 = np.asarray(v2["tv"], dtype=float)
    wv_v2 = np.asarray(v2["wv"], dtype=float)
    tv_p1 = np.asarray(ps1["tv"], dtype=float)
    wv_p1 = np.asarray(ps1["wv"], dtype=float)

    t_vn = vn["localtime"]
    v2_st512 = (v2["status"] & 512) != 0
    t_v2_raw = v2["localtime"]
    t_v2 = t_v2_raw + (16.0 if v2_st512 else 0.0)
    t_p1 = ps1["localtime"]
    exp_v2 = t_v2 - dv
    exp_ps = t_p1 - dps
    t_max_int = exp_v2 + 35.0

    # Shared x-limits centered cleanly on the waveforms (200 ns span matching Fig. 12)
    x_lo = float(tv_vn[0])
    x_hi = float(tv_vn[-1])

    # Compute exact reduced charge at 35 ns
    q_red35 = float(wv_vn[(tv_vn - exp_v2) < 35.0].sum()) * (2.0 / 50.0)

    ms = 2.0 if compact else 2.8
    lw_wf = 1.0 if compact else 1.25
    lw_v = 1.15 if compact else 1.5
    fs_leg = 6.9 if compact else 8.4
    fs_lbl = 8.2 if compact else 9.8

    # ---- TOP PANEL: FASERnu Veto (VetoNu) ----
    ax_vn.plot(tv_vn, wv_vn, color="#1f77b4", lw=lw_wf, zorder=3)
    ax_vn.plot(tv_vn, wv_vn, "k.", ms=ms, zorder=4)
    ymax_vn = max(20.0, float(np.max(wv_vn)) * 1.24)
    ax_vn.set_ylim(0.0, ymax_vn)

    # Shaded integration region: start of VetoNu waveform up to exp_v2 + 35 ns
    ax_vn.axvspan(x_lo, t_max_int, color="#c8c8c8", alpha=0.65, zorder=1,
                  label=f"Max. Int. time (+35 ns): {t_max_int:.1f} ns")
    if show_512_fix and v2_st512:
        # Also show where the wrong (uncorrected -16 ns) integration cutoff would have landed!
        ax_vn.axvline(t_max_int - 16.0, color="#7f7f7f", ls=":", lw=1.2, zorder=2,
                      label=f"No +16 ns cutoff: {t_max_int-16.0:.1f} ns")

    ax_vn.axvline(t_v2, color="#ff7f0e", ls="-.", lw=lw_v, zorder=5,
                  label=f"Veto time: {t_v2:.1f} ns" + (" (+16 ns)" if v2_st512 else ""))
    ax_vn.axvline(t_vn, color="#d62728", ls="--", lw=lw_v, zorder=6,
                  label=f"FASERν time: {t_vn:.1f} ns")
    ax_vn.axvline(exp_v2, color="#d4b400", ls=":", lw=lw_v + 0.3, zorder=7,
                  label=f"Exp. time (V2): {exp_v2:.1f} ns")

    # Yellow horizontal arrow for Delta t_Veto (matching Fig. 12!)
    arrow_y = ymax_vn * 0.90
    ax_vn.annotate(
        "", xy=(exp_v2, arrow_y), xytext=(t_v2, arrow_y),
        arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.25", color="#d4b400", lw=1.8),
        zorder=8,
    )
    # Dummy legend patch for Delta t_Veto
    ax_vn.plot([], [], color="#d4b400", lw=4.0, label=rf"$\Delta\bar{{t}}_{{\mathrm{{Veto}}}} = {fmt_ns(dv)}$ ns")

    ax_vn.set_title(col_title, fontsize=9.6 if compact else 10.3, fontweight="bold", pad=6)
    ax_vn.set_ylabel("FASERν veto [mV]" if compact else "Amplitude FASERν veto [mV]", fontsize=fs_lbl)
    leg1 = ax_vn.legend(
        loc="upper right", fontsize=fs_leg, frameon=True, fancybox=False,
        facecolor="white", edgecolor="#cccccc", framealpha=0.92, handlelength=1.8, labelspacing=0.22,
    )
    leg1.get_frame().set_linewidth(0.7)
    # Compact charge annotation at lower right of top panel
    pct = 100.0 * q_red35 / max(vn["raw_charge"], 1e-3)
    q_txt = (
        f"Raw Q = {vn['raw_charge']:.1f} pC\nRed Q(35 ns) = {q_red35:.1f} pC"
        + (f" ({pct:.0f}%, rescued <30)" if (is_neutrino and q_red35 < 30.0) else f" ({pct:.0f}%)")
    )
    ax_vn.text(
        0.98, 0.34 if is_neutrino else 0.08, q_txt, transform=ax_vn.transAxes, ha="right", va="bottom",
        fontsize=fs_leg, zorder=9, bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="#cccccc", alpha=0.92),
    )

    # ---- MIDDLE PANEL: Veto St2 ----
    if show_512_fix and v2_st512:
        # Show raw (-16 ns) waveform faded AND shifted (+16 ns) waveform in blue!
        ax_v2.plot(tv_v2, wv_v2, color="#b0b0b0", ls="--", lw=lw_wf, label=f"Raw Veto2 (status=512, t={t_v2_raw:.1f} ns)")
        ax_v2.plot(tv_v2 + 16.0, wv_v2, color="#1f77b4", lw=lw_wf, zorder=3, label="Shifted +16 ns")
        ax_v2.plot(tv_v2 + 16.0, wv_v2, "k.", ms=ms, zorder=4)
        ax_v2.axvline(t_v2_raw, color="#999999", ls=":", lw=1.1, zorder=2)
    else:
        ax_v2.plot(tv_v2, wv_v2, color="#1f77b4", lw=lw_wf, zorder=3)
        ax_v2.plot(tv_v2, wv_v2, "k.", ms=ms, zorder=4)
    ymax_v2 = max(20.0, float(np.max(wv_v2)) * 1.18)
    ax_v2.set_ylim(0.0, ymax_v2)
    ax_v2.axvline(t_v2, color="#ff7f0e", ls="-.", lw=lw_v, zorder=5)
    ax_v2.set_ylabel("Veto St2 [mV]" if compact else "Amplitude Veto [mV]", fontsize=fs_lbl)
    v2_info = f"Veto2 ch {v2['ch']} (st={v2['status']})\nQ = {v2['raw_charge']:.1f} pC"
    if show_512_fix and v2_st512:
        ax_v2.legend(loc="upper right", fontsize=fs_leg - 0.3, frameon=True, fancybox=False, edgecolor="#cccccc", framealpha=0.92)
    else:
        ax_v2.text(0.98, 0.88, v2_info, transform=ax_v2.transAxes, ha="right", va="top",
                   fontsize=fs_leg, bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="#cccccc", alpha=0.88))

    # ---- BOTTOM PANEL: Preshower (Dig1, plus Dig2 overlay in 2025-2026) ----
    ax_ps.plot(tv_p1, wv_p1, color="#1f77b4", lw=lw_wf, zorder=3, label=f"PS Dig1 ch {ps1['ch']}: {t_p1:.1f} ns")
    ax_ps.plot(tv_p1, wv_p1, "k.", ms=ms, zorder=4)
    ymax_ps = max(8.0, float(np.max(wv_p1)) * 1.28)
    if ps2 is not None and ps2.get("has_wf"):
        tv_p2 = np.asarray(ps2["tv"], dtype=float)
        wv_p2 = np.asarray(ps2["wv"], dtype=float)
        p2_512 = (ps2["status"] & 512) != 0
        t_p2_corr = ps2["localtime"] + (16.0 if p2_512 else 0.0)
        ax_ps.plot(tv_p2 + (16.0 if p2_512 else 0.0), wv_p2, color="#9467bd", ls="-", lw=lw_wf, zorder=3,
                   label=f"PS Dig2 ch {ps2['ch']}: {t_p2_corr:.1f} ns (extra −{ddig2:.2f})")
        ax_ps.axvline(t_p2_corr, color="#9467bd", ls=":", lw=lw_v - 0.2, zorder=5)
        ymax_ps = max(ymax_ps, float(np.max(wv_p2)) * 1.35)

    ax_ps.set_ylim(0.0, ymax_ps)
    ax_ps.axvline(t_p1, color="#2ca02c", ls="-.", lw=lw_v, zorder=5)
    ax_ps.axvline(exp_ps, color="#d4b400", ls=":", lw=lw_v + 0.2, zorder=6,
                  label=f"Exp. time (PS): {exp_ps:.1f} ns")
    # Green/gold arrow for Delta t_PS from Preshower Dig1 back to Exp. time!
    arrow_y_ps = ymax_ps * 0.86
    ax_ps.annotate(
        "", xy=(exp_ps, arrow_y_ps), xytext=(t_p1, arrow_y_ps),
        arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.25", color="#2ca02c", lw=1.7),
        zorder=8,
    )
    ax_ps.plot([], [], color="#2ca02c", lw=3.5, label=rf"$\Delta\bar{{t}}_{{\mathrm{{PS}}}} = {fmt_ns(dps)}$ ns")

    ax_ps.set_ylabel("Preshower [mV]" if compact else "Amplitude Preshower [mV]", fontsize=fs_lbl)
    ax_ps.set_xlabel("Time [ns]", fontsize=fs_lbl)
    leg3 = ax_ps.legend(
        loc="upper right", fontsize=fs_leg - 0.2, frameon=True, fancybox=False,
        facecolor="white", edgecolor="#cccccc", framealpha=0.92, handlelength=1.6, labelspacing=0.20,
    )
    leg3.get_frame().set_linewidth(0.7)

    for ax in axes:
        ax.set_xlim(x_lo, x_hi)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.grid(True, which="major", color="#e6e6e6", lw=0.6, zorder=0)
    ax_vn.yaxis.set_minor_locator(AutoMinorLocator(2))
    _clean_upper_yticks(ax_v2, ymax_v2)
    _clean_upper_yticks(ax_ps, ymax_ps)
    ax_vn.tick_params(labelbottom=False)
    ax_v2.tick_params(labelbottom=False)


def make_pair_fig(
    ev_left: dict, cfg_left: dict, title_left: str,
    ev_right: dict, cfg_right: dict, title_right: str,
    suptitle: str,
    show_512_right: bool = False,
    is_nu_left: bool = False,
    is_nu_right: bool = False,
):
    fig = plt.figure(figsize=(13.2, 7.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.25, 1.0, 1.0], hspace=0.0, wspace=0.17,
                          left=0.075, right=0.98, top=0.88, bottom=0.08)
    ax_l = [fig.add_subplot(gs[r, 0]) for r in range(3)]
    ax_r = [fig.add_subplot(gs[r, 1]) for r in range(3)]
    draw_fig12_column(ax_l, ev_left, cfg_left["dv"], cfg_left["dps"], cfg_left["ddig2"], title_left,
                      is_neutrino=is_nu_left)
    draw_fig12_column(ax_r, ev_right, cfg_right["dv"], cfg_right["dps"], cfg_right["ddig2"], title_right,
                      show_512_fix=show_512_right, is_neutrino=is_nu_right)
    fig.suptitle(suptitle, fontsize=12.4, fontweight="bold", y=0.955)
    return fig


def main():
    configure_rc()
    data = json.loads((STUDY_DIR / "fig12_waveform_study.json").read_text())
    out_dir = STUDY_DIR / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_pdf = out_dir / "fig12_all_periods_report.pdf"
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
        # 1. 2022 vs 2023 (showing 820 ns -> 220 ns NominalTriggerTime shift, identical Delta-t = 33.165 / 53.625 ns!)
        fig = make_pair_fig(
            data["2022"]["rep_muon_evt"], data["2022"],
            "2022 Data (run 8090) — data_default (Nom. trig 820 ns)",
            data["2023"]["rep_muon_evt"], data["2023"],
            "2023 Data (run 11703) — data_default (Nom. trig 220 ns)",
            "Fig. 12 Style — 2022 vs 2023 Real Muons (VetoVetoNu = 33.165 ns, PreshowerVetoNu = 53.625 ns)",
        )
        save(fig, "fig12_01_2022_vs_2023_waveforms",
             "Fig. 12-style 3-station waveforms for 2022 and 2023 real muons (NominalTriggerTime 820 -> 220 ns; identical Δt)", pdf)

        # 2. 2024 pre-CaloNu vs 2024 early-CaloNu (runs 15821-15844: still Dig1!)
        fig = make_pair_fig(
            data["2024_preCaloNu"]["rep_muon_evt"], data["2024_preCaloNu"],
            "2024 pre-CaloNu (run 15387) — data_default (Veto2 Dig1)",
            data["2024_earlyCaloNu"]["rep_muon_evt"], data["2024_earlyCaloNu"],
            "2024 early-CaloNu (run 15838) — data_default (Veto2 Dig1)",
            "Fig. 12 Style — 2024 pre-CaloNu vs early-CaloNu Real Muons (both Veto2 Dig1: 33.165 / 53.625 ns)",
        )
        save(fig, "fig12_02_2024_pre_vs_early_caloNu_waveforms",
             "Fig. 12-style waveforms for 2024 pre-CaloNu and early-CaloNu (both use data_default 33.165 / 53.625 ns)", pdf)

        # 3. 2024 during-CaloNu status==0 vs status==512 (+16 ns!)
        fig = make_pair_fig(
            data["2024_duringCaloNu"]["rep_muon_evt"], data["2024_duringCaloNu"],
            "2024 CaloNu (run 15849, Veto2 status=0) — data_caloNu",
            data["2024_duringCaloNu"]["rep_512_evt"], data["2024_duringCaloNu"],
            "2024 CaloNu (run 15849, Veto2 status=512: +16 ns fix)",
            "Fig. 12 Style — 2024 CaloNu (Veto2 Dig2 ch 22,23): status=0 vs status=512 (35.783 / 53.512 ns)",
            show_512_right=True,
        )
        save(fig, "fig12_03_2024_caloNu_status0_vs_status512",
             "Fig. 12-style waveforms for 2024 CaloNu showing Veto2 Dig2 +2.6 ns shift and status==512 +16 ns restoration", pdf)

        # 4. 2024 post-CaloNu vs 2025 (showing dramatic -40.3 ns Preshower upgrade jump!)
        fig = make_pair_fig(
            data["2024_postCaloNu"]["rep_muon_evt"], data["2024_postCaloNu"],
            "2024 post-CaloNu (run 16929) — data_default (PS = 53.625 ns)",
            data["2025"]["rep_muon_evt"], data["2025"],
            "2025 Data (run 20481, p0014) — data_2025 (PS = 13.27 ns)",
            "Fig. 12 Style — 2024 post-CaloNu vs 2025 Real Muons (Preshower jumps −40.3 ns ahead of Veto2!)",
        )
        save(fig, "fig12_04_2024_post_vs_2025_preshower_shift",
             "Fig. 12-style waveforms comparing 2024 post-CaloNu (PS=53.625 ns) vs 2025 upgraded 4-PMT Preshower (PS=13.27 ns, Dig2 +3.97 ns)", pdf)

        # 5. 2025 vs 2026 (4-PMT Preshower Dig1 + Dig2 consistency)
        fig = make_pair_fig(
            data["2025"]["rep_muon_evt"], data["2025"],
            "2025 Data (run 20481, p0014) — data_2025 (Dig2 +3.97 ns)",
            data["2026"]["rep_muon_evt"], data["2026"],
            "2026 Data (run 22404, p0015) — data_2026 (Dig2 +3.64 ns)",
            "Fig. 12 Style — 2025 vs 2026 Real Muons (Upgraded 4-PMT Preshower: Dig1 preferred, Dig2 fallback)",
        )
        save(fig, "fig12_05_2025_vs_2026_4pmt_waveforms",
             "Fig. 12-style waveforms for 2025 and 2026 showing Dig1 (13.27 ns) and Dig2 (+3.97 / +3.64 ns) on shared axis", pdf)

        # 6. MC24 Muon vs MC24 GENIE FaserNu Neutrino (exact Fig. 12 counterpart!)
        fig = make_pair_fig(
            data["MC24_Muon"]["rep_muon_evt"], data["MC24_Muon"],
            "Simulated Muon (MC24 100044) — mc preset",
            data["MC24_Genie_FaserNu"]["rep_genie_evt"], data["MC24_Genie_FaserNu"],
            "Simulated Neutrino (MC24 GENIE 200093, FaserNu)",
            "Fig. 12 Style — MC24 Simulated Muon vs Simulated Neutrino Interaction (27.21 / 51.41 ns)",
            is_nu_right=True,
        )
        save(fig, "fig12_06_mc24_muon_vs_genie_fasernu",
             "Exact Fig. 12 counterpart on MC24 Muon vs GENIE FaserNu delayed backsplash (raw 56.3 pC -> reduced <30 pC rescued)", pdf)

        # 7. MC24 GENIE FaserNu vs CaloNu Neutrino backsplash waveforms
        fig = make_pair_fig(
            data["MC24_Genie_FaserNu"]["rep_genie_evt"], data["MC24_Genie_FaserNu"],
            "Simulated Neutrino (GENIE 200093, FaserNu)",
            data["MC24_Genie_CaloNu"]["rep_genie_evt"], data["MC24_Genie_CaloNu"],
            "Simulated Neutrino (GENIE 200139, CaloNu)",
            "Fig. 12 Style — MC24 GENIE Neutrino Backsplash: FaserNu vs CaloNu Geometry (27.21 / 51.41 ns)",
            is_nu_left=True,
            is_nu_right=True,
        )
        save(fig, "fig12_07_mc24_genie_fasernu_vs_calonu",
             "Fig. 12-style waveforms for MC24 GENIE FaserNu and CaloNu neutrino backsplash events", pdf)

        # 8. VetoNu ensemble-averaged pulse shapes & cumulative charge capture across ALL 9 periods
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 4.9))
        fig.subplots_adjust(left=0.068, right=0.98, top=0.86, bottom=0.13, wspace=0.20)
        fig.suptitle(
            "VetoNu Waveform Shape & Cumulative Charge Capture Across All Periods (2022–2026 & MC24)",
            fontsize=12.4, fontweight="bold", y=0.955,
        )
        p_styles = [
            ("2022", "2022 (pre=42 ns, FWHM=20 ns)", "#1f77b4", "-", 1.6),
            ("2023", "2023 (pre=42 ns, FWHM=20 ns)", "#ff7f0e", "-", 1.6),
            ("2024_preCaloNu", "2024 pre (pre=23 ns, FWHM=20 ns)", "#2ca02c", "-", 1.5),
            ("2024_duringCaloNu", "2024 CaloNu (pre=24 ns, FWHM=20 ns)", "#d62728", "-", 1.7),
            ("2024_postCaloNu", "2024 post (pre=24 ns, FWHM=20 ns)", "#9467bd", "--", 1.5),
            ("2025", "2025 (pre=26 ns, FWHM=20 ns)", "#8c564b", "-", 1.6),
            ("2026", "2026 (pre=26 ns, FWHM=20 ns)", "#e377c2", "--", 1.6),
            ("MC24_Muon", "MC24 Muon (pre=12 ns, FWHM=18 ns)", "#111111", "--", 2.0),
        ]

        ax1.grid(True, color="#e0e0e0", lw=0.7)
        ax1.axvspan(-50, 35.0, color="#d8d8d8", alpha=0.40, label="Integration window (< +35 ns)")
        ax1.axvline(0.0, color="#d4b400", ls=":", lw=1.6, label=r"$t = t_{\mathrm{exp}}$ (40% rising edge)")
        ax1.axvline(35.0, color="#333333", ls="-.", lw=1.4, label=r"Cutoff $\Delta T = 35$ ns")
        for pk, lbl, col, ls, lw in p_styles:
            aa = data[pk]["aligned_avg"]
            x = np.asarray(aa["rel_t_ns"], dtype=float)
            y = np.asarray([np.nan if v is None else v for v in aa["vn_norm"]], dtype=float)
            ax1.plot(x, y, label=lbl, color=col, ls=ls, lw=lw)
        ax1.set_xlim(-46, 115)
        ax1.set_ylim(-0.04, 1.12)
        ax1.set_xlabel(r"$t - t_{\mathrm{exp}}$ [ns]", fontsize=10.8)
        ax1.set_ylabel("Peak-normalized VetoNu amplitude", fontsize=10.8)
        ax1.set_title("(a) Aligned average VetoNu pulse shape (2022–2026 identical; MC24 narrower)", fontsize=10.3)
        ax1.legend(loc="upper right", fontsize=7.3, frameon=True, framealpha=0.94, labelspacing=0.22)

        ax2.grid(True, color="#e0e0e0", lw=0.7)
        dts_num = np.array([15, 20, 25, 30, 35, 40, 45, 50, 60, 80], dtype=float)
        short_lbl = {
            "2022": "2022",
            "2023": "2023",
            "2024_preCaloNu": "2024 pre-CaloNu",
            "2024_duringCaloNu": "2024 CaloNu",
            "2024_postCaloNu": "2024 post-CaloNu",
            "2025": "2025",
            "2026": "2026",
            "MC24_Muon": "MC24 Muon",
        }
        ax2.axvline(35.0, color="#333333", ls="-.", lw=1.5, label=r"Nominal $\Delta T = 35$ ns (83.7–85.6%)")
        for pk, lbl, col, ls, lw in p_styles:
            sc = data[pk]["window_scan_v2"]
            caps = [100.0 * sc[f"{int(dt)}"]["cap_med"] for dt in dts_num]
            ax2.plot(dts_num, caps, marker="o", ms=3.5, label=short_lbl[pk], color=col, ls=ls, lw=lw)
        ax2.set_xlim(13, 82)
        ax2.set_ylim(44, 101)
        ax2.set_xlabel(r"Integration window $\Delta T = t_{\max} - t_{\mathrm{exp}}$ [ns]", fontsize=10.8)
        ax2.set_ylabel(r"Median charge capture $Q_{\mathrm{red}}(\Delta T) / Q_{\mathrm{raw}}$ [%]", fontsize=10.8)
        ax2.set_title("(b) Cumulative VetoNu charge capture vs integration cutoff", fontsize=10.3)
        ax2.legend(loc="lower right", fontsize=7.9, ncol=2, frameon=True, framealpha=0.94)
        save(fig, "fig12_08_vetonu_pulse_shape_and_capture_all_periods",
             "Aligned average VetoNu pulse shapes and cumulative charge capture curves across all 2022-2026 periods and MC24", pdf)

        # 9. Integration window optimization: Muon leakage (<30 pC) vs GENIE backsplash rescue
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 4.9))
        fig.subplots_adjust(left=0.068, right=0.98, top=0.86, bottom=0.13, wspace=0.22)
        fig.suptitle(
            r"Integration Window Optimization: Muon Leakage ($Q_{\mathrm{red}} < 30$ pC) vs GENIE Backsplash Rescue",
            fontsize=12.4, fontweight="bold", y=0.955,
        )

        ax1.grid(True, color="#e0e0e0", lw=0.7)
        ax1.axvline(35.0, color="#333333", ls="-.", lw=1.5, label=r"Nominal $\Delta T = 35$ ns")
        for pk, lbl, col, ls, mk in [
            ("2022", "2022 Veto2 (fb=1)", "#1f77b4", "-", "o"),
            ("2023", "2023 Veto2 (fb=1)", "#ff7f0e", "-", "s"),
            ("2024_duringCaloNu", "2024 CaloNu Veto2 (fb=1)", "#d62728", "-", "^"),
            ("2025", "2025 Veto2 (fb=1)", "#8c564b", "--", "v"),
            ("2026", "2026 Veto2 (fb=1)", "#e377c2", "--", "d"),
        ]:
            sc = data[pk]["window_scan_v2"]
            lk = [100.0 * sc[f"{int(dt)}"]["leak_lt30_frac"] for dt in dts_num]
            ax1.plot(dts_num, lk, marker=mk, ms=4.2, label=lbl, color=col, ls=ls, lw=1.6)
        # Overlay 2025/2026 Preshower Dig1 (fb=2) vs Dig2 (fb=3)!
        for pk, col_d1, col_d2 in [("2025", "#2ca02c", "#9467bd"), ("2026", "#17becf", "#111111")]:
            lk1 = [100.0 * data[pk]["window_scan_ps1"][f"{int(dt)}"]["leak_lt30_frac"] for dt in dts_num]
            lk2 = [100.0 * data[pk]["window_scan_ps2"][f"{int(dt)}"]["leak_lt30_frac"] for dt in dts_num]
            ax1.plot(dts_num, lk1, marker="o", ms=3.5, ls=":", lw=1.5, color=col_d1, label=f"{pk} PS Dig1 (fb=2: 0.00% at 35 ns)")
            ax1.plot(dts_num, lk2, marker="x", ms=5.0, ls="-", lw=1.8, color=col_d2, label=f"{pk} PS Dig2 (fb=3: far-tail outliers)")
        ax1.set_xlim(13, 82)
        ax1.set_ylim(-0.15, 3.2)
        ax1.set_xlabel(r"Integration window $\Delta T$ [ns]", fontsize=10.8)
        ax1.set_ylabel(r"Muon leakage $P(Q_{\mathrm{red}} < 30\ \mathrm{pC}\mid Q_{\mathrm{raw}} > 40\ \mathrm{pC})$ [%]", fontsize=10.2)
        ax1.set_title("(a) Muon leakage: Veto2 & PS Dig1 reach 0% by 30–35 ns", fontsize=10.3)
        ax1.legend(loc="upper right", fontsize=7.2, frameon=True, framealpha=0.94, labelspacing=0.20)

        ax2.grid(True, color="#e0e0e0", lw=0.7)
        ax2.axvline(35.0, color="#333333", ls="-.", lw=1.5, label=r"Nominal $\Delta T = 35$ ns (Pareto optimum)")
        g_fn = [100.0 * data["MC24_Genie_FaserNu"]["genie_summary"]["scan"][f"{int(dt)}"]["rescue_lt30_frac"] for dt in dts_num]
        g_fn2 = [100.0 * data["MC24_Genie_FaserNu"]["genie_summary"]["scan"][f"{int(dt)}"]["rescue_lt30_twosided_frac"] for dt in dts_num]
        g_cn = [100.0 * data["MC24_Genie_CaloNu"]["genie_summary"]["scan"][f"{int(dt)}"]["rescue_lt30_frac"] for dt in dts_num]
        ax2.plot(dts_num, g_fn, marker="o", ms=4.5, color="#1f77b4", lw=2.0,
                 label=r"GENIE FaserNu (200093): 1-sided ($t < t_{\mathrm{exp}}+\Delta T$)")
        ax2.plot(dts_num, g_fn2, marker="x", ms=5.0, color="#17becf", ls=":", lw=1.6,
                 label=r"GENIE FaserNu: 2-sided ($[-20, \Delta T]$ — identical)")
        ax2.plot(dts_num, g_cn, marker="s", ms=4.5, color="#d62728", lw=2.0,
                 label=r"GENIE CaloNu (200139): 1-sided ($t < t_{\mathrm{exp}}+\Delta T$)")
        ax2.axvspan(13, 30, color="#d62728", alpha=0.08, label=r"Muon leakage risk zone ($\Delta T < 30$ ns)")
        ax2.set_xlim(13, 82)
        ax2.set_ylim(0, 78)
        ax2.set_xlabel(r"Integration window $\Delta T$ [ns]", fontsize=10.8)
        ax2.set_ylabel(r"Backsplash rescue $P(Q_{\mathrm{red}} < 30\ \mathrm{pC}\mid Q_{\mathrm{raw}} \geq 40\ \mathrm{pC})$ [%]", fontsize=10.2)
        ax2.set_title("(b) Neutrino backsplash rescue drops steeply for ΔT > 35–40 ns", fontsize=10.3)
        ax2.legend(loc="upper right", fontsize=7.6, frameon=True, framealpha=0.94, labelspacing=0.24)
        save(fig, "fig12_09_window_optimization_leakage_vs_rescue",
             "Quantitative integration-window optimization comparing muon leakage across all periods vs GENIE backsplash rescue efficiency", pdf)

    (out_dir / "figure_index.json").write_text(json.dumps({
        "style_reference": "/afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis/docs/Electronic_2D_measurement_internal_note-6.pdf (Fig. 12, p. 21)",
        "report_pdf": str(report_pdf),
        "figures": index,
    }, indent=2))
    print(f"Generated {len(index)} Fig. 12-style waveform & window-study figures + {report_pdf}")


if __name__ == "__main__":
    main()
