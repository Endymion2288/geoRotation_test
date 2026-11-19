#!/usr/bin/env python3
"""Generate English validation figures (PNG, per-plot PDF, multi-page report PDF,
and bin-count JSON) for the 2022-2026 & MC NtupleDumper reduced-charge upgrade.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import uproot


def save_figure(fig: plt.Figure, out_dir: Path, stem: str, pdf_pages: PdfPages) -> dict[str, str]:
    png = out_dir / f"{stem}.png"
    pdf = out_dir / f"{stem}.pdf"
    fig.savefig(png, dpi=180, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    pdf_pages.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    return {"png": str(png), "pdf": str(pdf)}


def plot_fig01_timeline(
    phys_json: Path, xaod_json: Path, out_dir: Path, pdf_pages: PdfPages
) -> tuple[dict[str, str], dict[str, Any]]:
    phys = json.loads(phys_json.read_text())
    xaod = json.loads(xaod_json.read_text())

    # Combine run-level measurements (prefer PHYS when present, add xAOD for 2024 CaloNu & 2026 p0015)
    by_run: dict[int, dict[str, Any]] = {}
    for r in xaod:
        if "error" in r or r["period"].startswith("MC"):
            continue
        v2 = r.get("veto2_minus_vetonu", {}).get("peak")
        pst = r.get("ps_top_minus_vetonu", {}).get("peak")
        psb = r.get("ps_bot_minus_vetonu", {}).get("peak")
        if v2 is None:
            continue
        if int(r["run"]) >= 18194 and pst is not None and pst > 25.0:
            pst = None
        by_run[int(r["run"])] = {"run": int(r["run"]), "src": "xAOD", "v2": v2, "pst": pst, "psb": psb, "period": r["period"]}
    for r in phys:
        if "error" in r or r["period"].startswith("MC"):
            continue
        v2 = r.get("veto2_minus_vetonu", {}).get("peak")
        pst = r.get("ps_top_minus_vetonu", {}).get("peak")
        psb = r.get("ps_bot_minus_vetonu", {}).get("peak")
        if v2 is None:
            continue
        # Skip 2025 runs where PS0,1 failed GFIT (pst ~ 34 ns noise) when plotting clean PS top
        if int(r["run"]) >= 18194 and pst is not None and pst > 25.0:
            pst = None
        by_run[int(r["run"])] = {"run": int(r["run"]), "src": "PHYS", "v2": v2, "pst": pst, "psb": psb, "period": r["period"]}

    rows = [by_run[k] for k in sorted(by_run)]
    runs = np.array([r["run"] for r in rows])
    v2s = np.array([r["v2"] for r in rows], dtype=float)
    psts = np.array([np.nan if r["pst"] is None else r["pst"] for r in rows], dtype=float)
    psbs = np.array([np.nan if r["psb"] is None else r["psb"] for r in rows], dtype=float)

    # Also collect MC24 points
    mc_pts = []
    for r in phys:
        if r["period"].startswith("MC24") and "Genie" not in r["period"] and "error" not in r:
            v2 = r["veto2_minus_vetonu"]["peak"]
            ps = r["ps_all_minus_vetonu"]["peak"]
            if v2 is not None and ps is not None:
                mc_pts.append((r["period"], v2, ps))

    fig, axes = plt.subplots(2, 1, figsize=(12.5, 7.6), sharex=False, gridspec_kw={"hspace": 0.34})

    # Map runs onto ordinal x with clean year boundaries for readability
    x = np.arange(len(rows))
    ax = axes[0]
    ax.plot(x, v2s, "o", ms=5, color="#1f77b4", label="Measured throughgoing-muon peak (Data 2022-2026)")
    # Auto lookup curve
    auto_v2 = np.array([
        35.783 if (15845 <= r < 16928 or r == 18193) else (33.055 if r < 13847 else 33.165)
        for r in runs
    ])
    ax.step(x, auto_v2, where="mid", color="#d62728", lw=2.0, label="Auto lookup: Veto2 - VetoNu (33.055 / 33.165 / 35.783 ns)")
    ax.axhline(33.165, color="gray", ls="--", lw=1.0, alpha=0.6, label="Old hardcoded data default (33.165 ns)")
    # Highlight CaloNu Dig2 Veto runs
    cn_idx = np.where((runs >= 15845) & (runs < 16928))[0]
    if len(cn_idx):
        ax.axvspan(cn_idx[0] - 0.5, cn_idx[-1] + 0.5, color="#ff9896", alpha=0.22, label="2024 CaloNu (Veto2 on Dig2 ch 22,23)")
    ax.set_ylabel("Peak $\\Delta t = t_{\\mathrm{Veto2}} - t_{\\mathrm{VetoNu}}$ [ns]")
    ax.set_title("Throughgoing-Muon Scintillator Localtime Differences Across 2022-2026 Data (and MC24)")
    ax.set_ylim(31.8, 37.2)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.5, ncol=2)
    # Year dividers and x-ticks (evenly spaced without crowded final tick)
    tick_pos = list(range(0, len(rows), max(1, len(rows) // 14)))
    if (len(rows) - 1) - tick_pos[-1] >= 3:
        tick_pos.append(len(rows) - 1)
    else:
        tick_pos[-1] = len(rows) - 1
    ax.set_xticks(tick_pos)
    ax.set_xticklabels([str(runs[i]) for i in tick_pos], rotation=35, ha="right", fontsize=8)
    ax.set_xlabel("Run Number (2022 -> 2026 chronological index)")

    ax2 = axes[1]
    m_top = np.isfinite(psts)
    m_bot = np.isfinite(psbs)
    ax2.plot(x[m_top], psts[m_top], "o", ms=5, color="#2ca02c", label="Preshower Dig1 (ch 12,13) - VetoNu")
    ax2.plot(x[m_bot], psbs[m_bot], "s", ms=5, color="#ff7f0e", label="Preshower Dig2 (ch 24,25, +16 ns corr) - VetoNu")
    auto_ps_top = np.array([
        53.490 if r < 13847 else (53.550 if r < 15845 else (53.512 if r < 16928 else (53.550 if r < 18193 else 13.280)))
        for r in runs
    ])
    auto_ps_bot = np.where(runs >= 18194, 13.280 + 3.970, np.nan)
    ax2.step(x, auto_ps_top, where="mid", color="#1b7837", lw=2.0, label="Auto lookup: Dig1 (53.49-53.55 -> 13.28 ns)")
    m_ab = np.isfinite(auto_ps_bot)
    if np.any(m_ab):
        ax2.step(x[m_ab], auto_ps_bot[m_ab], where="mid", color="#d95f02", lw=2.0, label="Auto lookup: Dig2 (13.28 + 3.97 = 17.25 ns)")
    ax2.axhline(53.625, color="gray", ls="--", lw=1.0, alpha=0.6, label="Old hardcoded default (53.625 ns; -40 ns off in 2025-26)")
    if mc_pts:
        mc_v2 = np.mean([p[1] for p in mc_pts])
        mc_ps = np.mean([p[2] for p in mc_pts])
        ax2.text(
            0.985, 0.80,
            f"MC24 throughgoing muons ({len(mc_pts)} samples):\n"
            f"Veto2 - VetoNu = {mc_v2:.2f} ns (auto: 27.21 ns)\n"
            f"Preshower - VetoNu = {mc_ps:.2f} ns (auto: 51.41 ns)",
            transform=ax2.transAxes, fontsize=8.3, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f4ff", edgecolor="#4c72b0"),
        )
    ax2.set_ylabel("Peak $\\Delta t = t_{\\mathrm{Preshower}} - t_{\\mathrm{VetoNu}}$ [ns]")
    ax2.set_xlabel("Run Number (2022 -> 2026 chronological index)")
    ax2.set_ylim(5.0, 62.0)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="lower left", fontsize=8.0, ncol=1, framealpha=0.92)
    ax2.set_xticks(tick_pos)
    ax2.set_xticklabels([str(runs[i]) for i in tick_pos], rotation=35, ha="right", fontsize=8)

    paths = save_figure(fig, out_dir, "fig01_timing_diff_timeline_2022_2026_mc", pdf_pages)
    meta = {"n_runs_plotted": len(rows), "run_min": int(runs.min()), "run_max": int(runs.max()), "mc24_samples": len(mc_pts)}
    return paths, meta


def plot_fig02_status512(val_dir: Path, out_dir: Path, pdf_pages: PdfPages) -> tuple[dict[str, str], dict[str, Any]]:
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.6))
    meta: dict[str, Any] = {}

    # Panel A: 2024 CaloNu run 15849 Veto20 (ch 22, Dig2)
    fp_cn = val_dir / "ntuples/2024_CaloNu_015849/2024_CaloNu_015849_fix_auto.root"
    with uproot.open(fp_cn) as f:
        t = f["nt"]
        a = t.arrays(
            ["VetoNu0_localtime", "VetoNu0_status", "VetoNu1_localtime", "VetoNu1_status",
             "VetoNu0_raw_charge", "VetoNu1_raw_charge", "Veto20_localtime", "Veto20_status"],
            library="np",
        )
    mu = (a["VetoNu0_status"] == 0) & (a["VetoNu1_status"] == 0) & (a["VetoNu0_raw_charge"] > 40) & (a["VetoNu1_raw_charge"] > 40)
    tvn = 0.5 * (a["VetoNu0_localtime"] + a["VetoNu1_localtime"])
    dt0 = (a["Veto20_localtime"] - tvn)[mu & (a["Veto20_status"] == 0)]
    dt5r = (a["Veto20_localtime"] - tvn)[mu & (a["Veto20_status"] == 512)]
    dt5c = dt5r + 16.0

    bins_a = np.linspace(14.0, 42.0, 113)
    h0, _ = np.histogram(dt0, bins=bins_a)
    h5r, _ = np.histogram(dt5r, bins=bins_a)
    h5c, _ = np.histogram(dt5c, bins=bins_a)
    ax = axes[0]
    ax.step(bins_a[:-1], h0 / max(1, h0.max()), where="post", lw=2.0, color="#1f77b4", label=f"status == 0 (N={len(dt0)}, med={np.median(dt0):.2f} ns)")
    ax.step(bins_a[:-1], h5r / max(1, h5r.max()), where="post", lw=1.8, ls="--", color="#d62728", label=f"status == 512 raw (N={len(dt5r)}, med={np.median(dt5r):.2f} ns)")
    ax.step(bins_a[:-1], h5c / max(1, h5c.max()), where="post", lw=2.0, color="#2ca02c", label=f"status == 512 + 16 ns (med={np.median(dt5c):.2f} ns)")
    ax.set_xlabel("$t_{\\mathrm{Veto20}} - t_{\\mathrm{VetoNu}}$ [ns] (local time)")
    ax.set_ylabel("Peak-normalized counts")
    ax.set_title("2024 CaloNu (Run 15849): Veto2 Ch 22 (Dig2)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.2, loc="upper left")
    meta["calonu_15849_v20"] = {
        "n_st0": int(len(dt0)), "med_st0": round(float(np.median(dt0)), 3),
        "n_st512": int(len(dt5r)), "med_st512_raw": round(float(np.median(dt5r)), 3),
        "med_st512_corr": round(float(np.median(dt5c)), 3),
    }

    # Panel B: 2026 run 22404 Preshower2 (ch 24, Dig2)
    fp_26 = val_dir / "ntuples/2026_22404/2026_22404_fix_auto.root"
    with uproot.open(fp_26) as f:
        t = f["nt"]
        a = t.arrays(
            ["VetoNu0_localtime", "VetoNu0_status", "VetoNu1_localtime", "VetoNu1_status",
             "VetoNu0_raw_charge", "VetoNu1_raw_charge", "Preshower2_localtime", "Preshower2_status"],
            library="np",
        )
    mu = (a["VetoNu0_status"] == 0) & (a["VetoNu1_status"] == 0) & (a["VetoNu0_raw_charge"] > 40) & (a["VetoNu1_raw_charge"] > 40)
    tvn = 0.5 * (a["VetoNu0_localtime"] + a["VetoNu1_localtime"])
    dt0 = (a["Preshower2_localtime"] - tvn)[mu & (a["Preshower2_status"] == 0)]
    dt5r = (a["Preshower2_localtime"] - tvn)[mu & (a["Preshower2_status"] == 512)]
    dt5c = dt5r + 16.0

    bins_b = np.linspace(-6.0, 24.0, 121)
    h0, _ = np.histogram(dt0, bins=bins_b)
    h5r, _ = np.histogram(dt5r, bins=bins_b)
    h5c, _ = np.histogram(dt5c, bins=bins_b)
    ax = axes[1]
    ax.step(bins_b[:-1], h0 / max(1, h0.max()), where="post", lw=2.0, color="#1f77b4", label=f"status == 0 (N={len(dt0)}, med={np.median(dt0):.2f} ns)")
    ax.step(bins_b[:-1], h5r / max(1, h5r.max()), where="post", lw=1.8, ls="--", color="#d62728", label=f"status == 512 raw (N={len(dt5r)}, med={np.median(dt5r):.2f} ns)")
    ax.step(bins_b[:-1], h5c / max(1, h5c.max()), where="post", lw=2.0, color="#2ca02c", label=f"status == 512 + 16 ns (med={np.median(dt5c):.2f} ns)")
    ax.set_xlabel("$t_{\\mathrm{Preshower2}} - t_{\\mathrm{VetoNu}}$ [ns] (local time)")
    ax.set_ylabel("Peak-normalized counts")
    ax.set_title("2026 Data (Run 22404): Preshower Ch 24 (Dig2)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.2, loc="upper left")
    meta["data_2026_22404_ps2"] = {
        "n_st0": int(len(dt0)), "med_st0": round(float(np.median(dt0)), 3),
        "n_st512": int(len(dt5r)), "med_st512_raw": round(float(np.median(dt5r)), 3),
        "med_st512_corr": round(float(np.median(dt5c)), 3),
    }

    fig.suptitle("Verification of 2nd-Digitizer TIMING_SHIFTED (512) +16.0 ns Correction on Real Data", fontsize=12)
    paths = save_figure(fig, out_dir, "fig02_status512_16ns_alignment", pdf_pages)
    return paths, meta


def plot_fig03_ben_parity(val_dir: Path, summary: dict[str, Any], out_dir: Path, pdf_pages: PdfPages) -> tuple[dict[str, str], dict[str, Any]]:
    ben_samples = ["2022_008090", "2023_011703", "2024_pre_015387", "2024_CaloNu_015849", "2024_post_016929"]
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.7))
    meta: dict[str, Any] = {}

    # Panel A: Summary exact-parity bar chart across all 5 Ben folders
    labels = []
    exact_fracs = []
    auto_1pc_fracs = []
    n_acts = []
    for name in ben_samples:
        rec = summary.get(name, {})
        ov = rec.get("vs_ben_exact_override", {})
        au = rec.get("vs_ben_auto", {})
        if ov.get("n_active_timing", 0) > 0:
            labels.append(name.replace("_", "\n"))
            exact_fracs.append(100.0 * ov["dq0_exact_1e3_frac"])
            auto_1pc_fracs.append(100.0 * au.get("dq0_within_1pC_frac", 0.0))
            n_acts.append(ov["n_active_timing"])
            meta[name] = {"override": ov, "auto": au}

    x = np.arange(len(labels))
    w = 0.36
    ax = axes[0]
    b1 = ax.bar(x - w / 2, exact_fracs, width=w, color="#1f77b4", label="Exact override: $|\\Delta q_{35}| < 10^{-3}$ pC")
    b2 = ax.bar(x + w / 2, auto_1pc_fracs, width=w, color="#2ca02c", label="Default auto lookup: $|\\Delta q_{35}| < 1.0$ pC")
    for bar, val in zip(b1, exact_fracs):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.04, f"{val:.2f}%", ha="center", va="bottom", fontsize=7.4, rotation=0)
    for bar, val in zip(b2, auto_1pc_fracs):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.04, f"{val:.2f}%", ha="center", va="bottom", fontsize=7.4, rotation=0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylim(98.5, 100.35)
    ax.set_ylabel("Event-by-event agreement [%] (zoomed 98.5%-100%)")
    ax.set_title("Parity vs Ben Wilson Reference Files (5 Periods)")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.2)

    # Panel B: Event-by-event scatter on 2024 CaloNu 15849 (Fixed Auto vs Ben Reference)
    name = "2024_CaloNu_015849"
    nt_fp = val_dir / f"ntuples/{name}/{name}_fix_auto.root"
    ben_fp = summary[name]["vs_ben_auto"]["ben_file"]
    with uproot.open(nt_fp) as f:
        a = f["nt"].arrays(["eventID", "VetoNu0_reduced_charge", "VetoNu_reduced_charge_fallback"], library="np")
    with uproot.open(ben_fp) as f:
        b = f["tree"].arrays(["event", "myVetoNu0_rawCharge35"], library="np")
    idx_map = {int(e): i for i, e in enumerate(b["event"])}
    q_nt, q_ben = [], []
    for i, e in enumerate(a["eventID"]):
        if a["VetoNu_reduced_charge_fallback"][i] <= 0:
            continue
        j = idx_map.get(int(e))
        if j is not None and np.isfinite(b["myVetoNu0_rawCharge35"][j]):
            q_nt.append(a["VetoNu0_reduced_charge"][i])
            q_ben.append(b["myVetoNu0_rawCharge35"][j])
    q_nt_arr = np.asarray(q_nt)
    q_ben_arr = np.asarray(q_ben)
    ax2 = axes[1]
    ax2.scatter(q_ben_arr, q_nt_arr, s=4, alpha=0.35, color="#1f77b4", rasterized=True)
    lim = max(10.0, float(np.percentile(np.concatenate([q_ben_arr, q_nt_arr]), 99.5)))
    ax2.plot([-5, lim], [-5, lim], "--", color="#d62728", lw=1.4, label="$y = x$ identity")
    ax2.set_xlim(-5, lim)
    ax2.set_ylim(-5, lim)
    ax2.set_xlabel("Ben Wilson myVetoNu0_rawCharge35 [pC]")
    ax2.set_ylabel("Fixed NtupleDumper VetoNu0_reduced_charge [pC]")
    max_diff = float(np.max(np.abs(q_nt_arr - q_ben_arr))) if len(q_nt_arr) else 0.0
    ax2.set_title(f"2024 CaloNu Run 15849 Event-by-Event (N={len(q_nt_arr)}, max |$\\Delta$|={max_diff:.2e} pC)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper left", fontsize=8.5)

    paths = save_figure(fig, out_dir, "fig03_parity_vs_ben_wilson", pdf_pages)
    return paths, meta


def plot_fig04_smoking_guns(val_dir: Path, summary: dict[str, Any], out_dir: Path, pdf_pages: PdfPages) -> tuple[dict[str, str], dict[str, Any]]:
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.8))
    meta: dict[str, Any] = {}

    # (a) 2024 CaloNu 15849: Veto2 status==512 & no good PS (baseline fell back to RAW, fix recovered Veto2!)
    ax = axes[0, 0]
    with uproot.open(val_dir / "ntuples/2024_CaloNu_015849/2024_CaloNu_015849_fix_auto.root") as ff, \
         uproot.open(val_dir / "ntuples/2024_CaloNu_015849/2024_CaloNu_015849_baseline.root") as fb:
        af = ff["nt"].arrays(
            ["VetoNu0_raw_charge", "VetoNu1_raw_charge", "VetoNu0_reduced_charge",
             "Veto20_status", "Veto21_status", "Preshower0_status", "Preshower1_status", "VetoNu_reduced_charge_fallback"],
            library="np",
        )
        ab = fb["nt"].arrays(["VetoNu0_reduced_charge"], library="np")
    m512_nops = (
        ((af["Veto20_status"] == 512) | (af["Veto21_status"] == 512))
        & ((af["Veto20_status"] & ~512) == 0)
        & ((af["Veto21_status"] & ~512) == 0)
        & (af["Preshower0_status"] != 0)
        & (af["Preshower1_status"] != 0)
    )
    bins = np.linspace(0, 260, 53)
    ax.hist(np.clip(af["VetoNu0_raw_charge"][m512_nops], 0, 260), bins=bins, histtype="step", lw=1.6, ls=":", color="black",
            label=f"Raw charge (med={np.median(af['VetoNu0_raw_charge'][m512_nops]):.1f} pC)")
    ax.hist(np.clip(ab["VetoNu0_reduced_charge"][m512_nops], 0, 260), bins=bins, histtype="step", lw=2.0, color="#d62728",
            label=f"Baseline latest: fell back to raw (med={np.median(ab['VetoNu0_reduced_charge'][m512_nops]):.1f})")
    ax.hist(np.clip(af["VetoNu0_reduced_charge"][m512_nops], 0, 260), bins=bins, histtype="step", lw=2.0, color="#1f77b4",
            label=f"Fixed fix_v1: Veto2 +16 ns (med={np.median(af['VetoNu0_reduced_charge'][m512_nops]):.1f})")
    ax.set_xlabel("VetoNu0 Charge [pC] (clipped at 260 pC)")
    ax.set_ylabel("Events / 5 pC")
    ax.set_title(f"(a) 2024 CaloNu 15849: Veto2 status==512, No Good PS (N={int(m512_nops.sum())})")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7.8)

    # (b) 2025-2026 Preshower-fallback (fallback==2) muons: Fixed vs Baseline!
    ax = axes[0, 1]
    q_raw_l, q_fix_l, q_base_l = [], [], []
    for nm in ["2025_20012", "2025_20481", "2026_22404", "2026_23176"]:
        fp_f = val_dir / f"ntuples/{nm}/{nm}_fix_auto.root"
        fp_b = val_dir / f"ntuples/{nm}/{nm}_baseline.root"
        if not (fp_f.is_file() and fp_b.is_file()):
            continue
        with uproot.open(fp_f) as ff, uproot.open(fp_b) as fb:
            af_s = ff["nt"].arrays(["VetoNu0_raw_charge", "VetoNu1_raw_charge", "VetoNu0_reduced_charge", "VetoNu_reduced_charge_fallback"], library="np")
            ab_s = fb["nt"].arrays(["VetoNu0_reduced_charge"], library="np")
        m = (af_s["VetoNu0_raw_charge"] > 40) & (af_s["VetoNu1_raw_charge"] > 40) & (af_s["VetoNu_reduced_charge_fallback"] == 2)
        q_raw_l.append(af_s["VetoNu0_raw_charge"][m])
        q_fix_l.append(af_s["VetoNu0_reduced_charge"][m])
        q_base_l.append(ab_s["VetoNu0_reduced_charge"][m])
    qr = np.concatenate(q_raw_l) if q_raw_l else np.array([])
    qf = np.concatenate(q_fix_l) if q_fix_l else np.array([])
    qb = np.concatenate(q_base_l) if q_base_l else np.array([])
    bins_b = np.linspace(-5, 300, 62)
    if len(qf):
        ax.hist(np.clip(qr, -5, 300), bins=bins_b, histtype="step", lw=1.6, ls=":", color="black", label=f"Raw charge (med={np.median(qr):.1f} pC)")
        ax.hist(np.clip(qb, -5, 300), bins=bins_b, histtype="step", lw=2.0, color="#d62728", label=f"Baseline latest: <40 pC = {100*np.mean(qb<40):.1f}% (med={np.median(qb):.1f})")
        ax.hist(np.clip(qf, -5, 300), bins=bins_b, histtype="step", lw=2.0, color="#1f77b4", label=f"Fixed fix_v1: <40 pC = {100*np.mean(qf<40):.1f}% (med={np.median(qf):.1f})")
        ax.axvline(40.0, color="gray", ls="--", lw=1.0)
        meta["ps_fallback_2025_2026"] = {
            "n": int(len(qf)),
            "raw_med": round(float(np.median(qr)), 2),
            "base_med": round(float(np.median(qb)), 2),
            "fix_med": round(float(np.median(qf)), 2),
            "base_lt40_pct": round(100.0 * float(np.mean(qb < 40)), 2),
            "fix_lt40_pct": round(100.0 * float(np.mean(qf < 40)), 2),
        }
    ax.set_xlabel("VetoNu0 Charge [pC] (clipped at 300 pC)")
    ax.set_ylabel("Events / 5 pC")
    ax.set_title(f"(b) 2025-2026 Preshower-Fallback Muons (N={len(qf)})")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7.8)

    # (c) 2024 CaloNu 15849 V2_st0 muons: Baseline(33.165 ns) - Fixed(35.783 ns) timing-bias residual
    ax = axes[1, 0]
    m_v2_st0 = (af["Veto20_status"] == 0) & (af["Veto21_status"] == 0) & (af["VetoNu0_raw_charge"] > 40)
    dq_st0 = ab["VetoNu0_reduced_charge"][m_v2_st0] - af["VetoNu0_reduced_charge"][m_v2_st0]
    bins_c = np.linspace(-2.0, 20.0, 89)
    ax.hist(np.clip(dq_st0, -2.0, 20.0), bins=bins_c, histtype="stepfilled", alpha=0.35, color="#d62728", edgecolor="#d62728", lw=1.5,
            label=f"Baseline (33.165 ns) $-$ Fixed (35.783 ns)\nmed={np.median(dq_st0):+.2f} pC, >1 pC: {100*np.mean(dq_st0>1):.1f}%")
    ax.axvline(0.0, color="black", ls="--", lw=1.2)
    ax.set_xlabel("$q_{35}^{\\mathrm{baseline}} - q_{35}^{\\mathrm{fixed}}$ [pC] (clipped to [-2, 20] pC)")
    ax.set_ylabel("Events / 0.25 pC")
    ax.set_title(f"(c) 2024 CaloNu 15849: Veto2 status==0 Timing-Bias Residual (N={int(m_v2_st0.sum())})")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.2, loc="upper right")
    meta["calonu_st0_bias"] = {
        "n": int(m_v2_st0.sum()),
        "med_diff_pC": round(float(np.median(dq_st0)), 3),
        "p95_diff_pC": round(float(np.percentile(dq_st0, 95)), 3),
        "gt1pC_pct": round(100.0 * float(np.mean(dq_st0 > 1.0)), 2),
    }

    # (d) MC24 muon (100044) & Genie neutrino (200139): raw vs reduced charge
    ax = axes[1, 1]
    bins_d = np.linspace(0, 300, 61)
    for nm, col, lbl in [
        ("MC24_PG_100044", "#1f77b4", "MC24 Muon 100044"),
        ("MC24_Genie_200139", "#2ca02c", "MC24 Genie $\\nu$ 200139"),
    ]:
        fp_f = val_dir / f"ntuples/{nm}/{nm}_fix_auto.root"
        if not fp_f.is_file():
            continue
        with uproot.open(fp_f) as ff:
            af = ff["nt"].arrays(["VetoNu0_raw_charge", "VetoNu0_reduced_charge", "VetoNu_reduced_charge_fallback"], library="np")
        act = af["VetoNu_reduced_charge_fallback"] > 0
        qr_m = af["VetoNu0_raw_charge"][act]
        qf_m = af["VetoNu0_reduced_charge"][act]
        ax.hist(np.clip(qr_m, 0, 300), bins=bins_d, density=True, histtype="step", lw=1.5, ls="--", color=col, label=f"{lbl} raw (<40 pC: {100*np.mean(qr_m<40):.1f}%)")
        ax.hist(np.clip(qf_m, 0, 300), bins=bins_d, density=True, histtype="step", lw=2.0, color=col, label=f"{lbl} reduced (<40 pC: {100*np.mean(qf_m<40):.1f}%)")
        meta[nm] = {
            "n_active": int(act.sum()),
            "raw_lt40_pct": round(100.0 * float(np.mean(qr_m < 40)), 2),
            "red_lt40_pct": round(100.0 * float(np.mean(qf_m < 40)), 2),
            "raw_med": round(float(np.median(qr_m)), 2),
            "red_med": round(float(np.median(qf_m)), 2),
        }
    ax.axvline(40.0, color="gray", ls="--", lw=1.0)
    ax.set_xlabel("VetoNu0 Charge [pC] (clipped at 300 pC)")
    ax.set_ylabel("Normalized density")
    ax.set_title("(d) MC24 Muon Preservation vs Genie Neutrino Backsplash Rejection")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7.8)

    fig.tight_layout()
    paths = save_figure(fig, out_dir, "fig04_fix_vs_baseline_smoking_guns", pdf_pages)
    return paths, meta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phys-json", type=Path, required=True)
    ap.add_argument("--xaod-json", type=Path, required=True)
    ap.add_argument("--val-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = json.loads((args.val_dir / "validation_summary.json").read_text())

    report_pdf = args.out_dir / "waveform_reduced_charge_report.pdf"
    index: dict[str, Any] = {"report_pdf": str(report_pdf), "figures": {}}
    with PdfPages(report_pdf) as pdf_pages:
        p1, m1 = plot_fig01_timeline(args.phys_json, args.xaod_json, args.out_dir, pdf_pages)
        index["figures"]["fig01"] = {"paths": p1, "meta": m1}
        p2, m2 = plot_fig02_status512(args.val_dir, args.out_dir, pdf_pages)
        index["figures"]["fig02"] = {"paths": p2, "meta": m2}
        p3, m3 = plot_fig03_ben_parity(args.val_dir, summary, args.out_dir, pdf_pages)
        index["figures"]["fig03"] = {"paths": p3, "meta": m3}
        p4, m4 = plot_fig04_smoking_guns(args.val_dir, summary, args.out_dir, pdf_pages)
        index["figures"]["fig04"] = {"paths": p4, "meta": m4}

    idx_path = args.out_dir / "figures_index.json"
    idx_path.write_text(json.dumps(index, indent=2) + "\n")
    print("Generated all figures and report PDF:", report_pdf)


if __name__ == "__main__":
    main()
