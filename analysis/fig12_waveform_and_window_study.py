#!/usr/bin/env python3
"""Multi-core real-xAOD waveform & reduced-charge integration window study across
2022-2026 data and MC24 (Muon + GENIE FaserNu/CaloNu).

Extracts:
  1. Representative single-event Fig. 12-style waveforms (VetoNu, Veto St2, Preshower Dig1,
     and in 2025-2026 Preshower Dig2) with exact localtime, expected time, Delta-t, and
     integration shading for every data-taking period and MC24.
  2. Ensemble-averaged VetoNu pulse shapes aligned to (t - t_exp) across all periods.
  3. Comprehensive integration-window scan (Delta-T = 15..80 ns) & pre-pulse lower-bound
     audit (Q_pre for t < t_exp - 15 ns, undershoot, Veto2 vs Preshower Dig1 vs Dig2
     timing jitter impact, muon leakage rate vs GENIE backsplash rescue efficiency).
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np

DEFAULT_OUTDIR = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/"
    "stage04_waveform_reduced_charge/fig12_waveform_window_study"
)

STATUS_TIMING_SHIFTED = 1 << 9  # 512
STATUS_SECONDARY = 1 << 1       # 2
INT_TO_CHARGE = 2.0 / 50.0      # 2 ns / 50 Ohm -> pC

DT_SCAN = [15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 80.0]

# Exact period table matching NtupleDumperConfig.REDUCED_CHARGE_TIMING_PRESETS
PERIOD_SAMPLES: list[dict[str, Any]] = [
    {
        "period": "2022",
        "label": "2022 Data (run 8090)",
        "preset": "data_default",
        "dv": 33.165,
        "dps": 53.625,
        "ddig2": 0.0,
        "v2_ch": [6, 7],
        "ps_dig1_ch": [12, 13],
        "ps_dig2_ch": [],
        "vn_ch": [4, 5],
        "xaod": "/eos/experiment/faser/rec/2022/r0021/008090/Faser-Physics-008090-00010-xAOD.root",
        "is_mc": False,
        "is_genie": False,
        "max_events": 25000,
    },
    {
        "period": "2023",
        "label": "2023 Data (run 11703)",
        "preset": "data_default",
        "dv": 33.165,
        "dps": 53.625,
        "ddig2": 0.0,
        "v2_ch": [6, 7],
        "ps_dig1_ch": [12, 13],
        "ps_dig2_ch": [],
        "vn_ch": [4, 5],
        "xaod": "/eos/experiment/faser/rec/2023/r0021/011703/Faser-Physics-011703-00010-xAOD.root",
        "is_mc": False,
        "is_genie": False,
        "max_events": 25000,
    },
    {
        "period": "2024_preCaloNu",
        "label": "2024 pre-CaloNu (run 15387)",
        "preset": "data_default",
        "dv": 33.165,
        "dps": 53.625,
        "ddig2": 0.0,
        "v2_ch": [6, 7],
        "ps_dig1_ch": [12, 13],
        "ps_dig2_ch": [],
        "vn_ch": [4, 5],
        "xaod": "/eos/experiment/faser/rec/2024/r0022/015387/Faser-Physics-015387-00010-r0022-xAOD.root",
        "is_mc": False,
        "is_genie": False,
        "max_events": 25000,
    },
    {
        "period": "2024_earlyCaloNu",
        "label": "2024 early-CaloNu (run 15838)",
        "preset": "data_default",
        "dv": 33.165,
        "dps": 53.625,
        "ddig2": 0.0,
        "v2_ch": [6, 7],
        "ps_dig1_ch": [12, 13],
        "ps_dig2_ch": [],
        "vn_ch": [4, 5],
        "xaod": "/eos/experiment/faser/rec/2024/r0022/015838/Faser-Physics-015838-00010-r0022-xAOD.root",
        "is_mc": False,
        "is_genie": False,
        "max_events": 25000,
    },
    {
        "period": "2024_duringCaloNu",
        "label": "2024 CaloNu (run 15849)",
        "preset": "data_caloNu",
        "dv": 35.783,
        "dps": 53.512,
        "ddig2": 0.0,
        "v2_ch": [22, 23],
        "ps_dig1_ch": [12, 13],
        "ps_dig2_ch": [],
        "vn_ch": [4, 5],
        "xaod": "/eos/experiment/faser/rec/2024/r0022/015849/Faser-Physics-015849-00010-r0022-xAOD.root",
        "is_mc": False,
        "is_genie": False,
        "max_events": 30000,
    },
    {
        "period": "2024_postCaloNu",
        "label": "2024 post-CaloNu (run 16929)",
        "preset": "data_default",
        "dv": 33.165,
        "dps": 53.625,
        "ddig2": 0.0,
        "v2_ch": [6, 7],
        "ps_dig1_ch": [12, 13],
        "ps_dig2_ch": [],
        "vn_ch": [4, 5],
        "xaod": "/eos/experiment/faser/rec/2024/r0022/016929/Faser-Physics-016929-00010-r0022-xAOD.root",
        "is_mc": False,
        "is_genie": False,
        "max_events": 25000,
    },
    {
        "period": "2025",
        "label": "2025 Data (run 20481, p0014)",
        "preset": "data_2025",
        "dv": 33.165,
        "dps": 13.27,
        "ddig2": 3.97,
        "v2_ch": [6, 7],
        "ps_dig1_ch": [12, 13],
        "ps_dig2_ch": [24, 25],
        "vn_ch": [4, 5],
        "xaod": "/eos/experiment/faser/rec/2025/p0014/020481/Faser-Physics-020481-00276-p0014-xAOD.root",
        "is_mc": False,
        "is_genie": False,
        "max_events": 25000,
    },
    {
        "period": "2026",
        "label": "2026 Data (run 22404, p0015)",
        "preset": "data_2026",
        "dv": 33.165,
        "dps": 13.27,
        "ddig2": 3.64,
        "v2_ch": [6, 7],
        "ps_dig1_ch": [12, 13],
        "ps_dig2_ch": [24, 25],
        "vn_ch": [4, 5],
        "xaod": "/eos/experiment/faser/rec/2026/p0015/022404/Faser-Physics-022404-00010-p0015-xAOD.root",
        "is_mc": False,
        "is_genie": False,
        "max_events": 25000,
    },
    {
        "period": "MC24_Muon",
        "label": "MC24 Muon (100044)",
        "preset": "mc",
        "dv": 27.21,
        "dps": 51.41,
        "ddig2": 0.0,
        "v2_ch": [6, 7],
        "ps_dig1_ch": [12, 13],
        "ps_dig2_ch": [],
        "vn_ch": [4, 5],
        "xaod": "/eos/experiment/faser/data0/sim/mc24/particle_gun/100044/rec/s0013-r0023/FaserMC-MC24_PG_mupl_fasernu_5mrad_flukaE-100044-00100-00199-s0013-r0023-xAOD.root",
        "is_mc": True,
        "is_genie": False,
        "max_events": 6000,
    },
    {
        "period": "MC24_Genie_FaserNu",
        "label": "MC24 GENIE FaserNu (200093)",
        "preset": "mc",
        "dv": 27.21,
        "dps": 51.41,
        "ddig2": 0.0,
        "v2_ch": [6, 7],
        "ps_dig1_ch": [12, 13],
        "ps_dig2_ch": [],
        "vn_ch": [4, 5],
        "xaod": "/eos/experiment/fasernu-data0/faser/sim/mc24/genie/200093/rec/s0014/FaserMC-MC24_Genie_24_xing_angle_wo_caloNu_epos_pi-200093-00000-00009-s0014-xAOD.root",
        "is_mc": True,
        "is_genie": True,
        "max_events": 10000,
    },
    {
        "period": "MC24_Genie_CaloNu",
        "label": "MC24 GENIE CaloNu (200139)",
        "preset": "mc",
        "dv": 27.21,
        "dps": 51.41,
        "ddig2": 0.0,
        "v2_ch": [6, 7],
        "ps_dig1_ch": [12, 13],
        "ps_dig2_ch": [],
        "vn_ch": [4, 5],
        "xaod": "/eos/experiment/fasernu-data0/faser/sim/mc24/genie/200139/rec/s0015/FaserMC-MC24_Genie_24_xing_angle_caloNu_epos_pi-200139-00000-00009-s0015-xAOD.root",
        "is_mc": True,
        "is_genie": True,
        "max_events": 10000,
    },
]


def extract_primary_hits(container, aux, target_chs: set[int]) -> dict[int, dict[str, Any]]:
    """Extract primary (non-secondary) WaveformHit info by channel."""
    container.setStore(aux)
    out: dict[int, dict[str, Any]] = {}
    for h in container:
        ch = int(h.channel())
        if ch not in target_chs:
            continue
        st = int(h.hit_status())
        is_sec = (st & STATUS_SECONDARY) != 0
        if ch in out and is_sec:
            continue
        tv_cpp = h.time_vector()
        wv_cpp = h.wave_vector()
        n_wf = len(tv_cpp)
        rec = {
            "ch": ch,
            "status": st,
            "localtime": float(h.localtime()),
            "raw_charge": float(h.raw_integral()) / 50.0,
            "charge": float(h.integral()) / 50.0,
            "raw_peak": float(h.raw_peak()),
            "peak": float(h.peak()),
            "width": float(h.width()),
            "has_wf": n_wf > 10,
        }
        if n_wf > 10:
            rec["tv"] = [float(x) for x in tv_cpp]
            rec["wv"] = [float(x) for x in wv_cpp]
        out[ch] = rec
    return out


def station_mean_time(hits: dict[int, dict[str, Any]], chs: list[int], dig2_extra: float = 0.0) -> tuple[float, bool]:
    """Mirror NtupleDumperAlg::get_mean_scintillator_time; also returns whether any bit-512 hit was used."""
    s = 0.0
    n = 0
    used_512 = False
    for ch in chs:
        h = hits.get(ch)
        if h is None:
            continue
        st = h["status"]
        if (st & ~STATUS_TIMING_SHIFTED) != 0:
            continue
        t = h["localtime"]
        if (st & STATUS_TIMING_SHIFTED) != 0:
            t += 16.0
            used_512 = True
        if dig2_extra != 0.0 and ch >= 16:
            t -= dig2_extra
        s += t
        n += 1
    return (s / n if n > 0 else float("nan"), used_512)


def process_period(cfg: dict[str, Any]) -> dict[str, Any]:
    import ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    ROOT.gSystem.Load("libxAODFaserWaveform.so")

    t0 = time.perf_counter()
    period = cfg["period"]
    f = ROOT.TFile.Open(cfg["xaod"])
    tree = f.Get("CollectionTree")
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("VetoNuWaveformHits*", 1)
    tree.SetBranchStatus("VetoWaveformHits*", 1)
    tree.SetBranchStatus("PreshowerWaveformHits*", 1)

    n_entries = int(tree.GetEntries())
    max_ev = min(n_entries, cfg["max_events"])
    # Start slightly into data runs to skip initial empty/ramp events
    start_ev = 0 if cfg["is_mc"] else min(2000, max(0, n_entries - max_ev))
    end_ev = min(n_entries, start_ev + max_ev)

    vn_set = set(cfg["vn_ch"])
    v2_set = set(cfg["v2_ch"])
    ps1_set = set(cfg["ps_dig1_ch"])
    ps2_set = set(cfg["ps_dig2_ch"])
    ps_all_set = ps1_set | ps2_set

    dv = cfg["dv"]
    dps = cfg["dps"]
    ddig2 = cfg["ddig2"]

    # Aligned average VetoNu waveform grid: rel_t in [-50, +130] ns (step 2 ns -> 91 bins)
    rel_grid = np.arange(-50.0, 132.0, 2.0)
    avg_vn_sum = np.zeros_like(rel_grid)
    avg_vn_norm_sum = np.zeros_like(rel_grid)  # normalized to unit peak
    avg_vn_cnt = np.zeros_like(rel_grid)
    avg_v2_sum = np.zeros_like(rel_grid)
    avg_v2_cnt = np.zeros_like(rel_grid)
    avg_ps1_sum = np.zeros_like(rel_grid)
    avg_ps1_cnt = np.zeros_like(rel_grid)

    # Per-muon metrics for window study
    vn_t_start_list = []
    vn_t_end_list = []
    vn_pre_len_list = []     # t_exp - t_start
    vn_width_list = []       # CB width sigma
    vn_fwhm_list = []
    vn_tpeak_rel_list = []   # t_peak - t_vn_localtime
    vn_q_pre15_list = []     # integral for (t - t_exp) < -15 ns
    vn_q_pre20_list = []     # integral for (t - t_exp) < -20 ns
    vn_min_post_list = []    # min mV in [t_exp+20, t_exp+120] (undershoot check)
    vn_raw_q_list = []
    vn_st512_cnt = 0

    # Capture fraction & leakage arrays vs DT_SCAN for Veto2 timing, PS Dig1 timing, PS Dig2 timing
    q_red_v2_by_dt = {dt: [] for dt in DT_SCAN}
    q_red_v2_twosided_by_dt = {dt: [] for dt in DT_SCAN}  # -20 <= exp_time < dt
    q_red_ps1_by_dt = {dt: [] for dt in DT_SCAN}
    q_red_ps2_by_dt = {dt: [] for dt in DT_SCAN}
    q_raw_ps1_list = []
    q_raw_ps2_list = []
    dt_v2_meas = []
    dt_ps1_meas = []
    dt_ps2_meas = []

    # GENIE backsplash metrics (max of VetoNu0 & VetoNu1)
    genie_raw_max = []
    genie_red_max_by_dt = {dt: [] for dt in DT_SCAN}
    genie_red_twosided_by_dt = {dt: [] for dt in DT_SCAN}
    genie_dt_late = []

    # Best representative single-event Fig. 12 displays
    best_muon_evt = None
    best_muon_score = 1e9
    best_512_evt = None
    best_512_score = 1e9
    best_genie_evt = None
    best_genie_score = 1e9

    n_muons = 0
    n_muons_with_all_wf = 0
    n_genie_sel = 0

    for i in range(start_ev, end_ev):
        tree.GetEntry(i)
        vn_hits = extract_primary_hits(tree.VetoNuWaveformHits, getattr(tree, "VetoNuWaveformHitsAux."), vn_set)
        v2_hits = extract_primary_hits(tree.VetoWaveformHits, getattr(tree, "VetoWaveformHitsAux."), v2_set)
        ps_hits = extract_primary_hits(tree.PreshowerWaveformHits, getattr(tree, "PreshowerWaveformHitsAux."), ps_all_set)

        vn0 = vn_hits.get(cfg["vn_ch"][0])
        vn1 = vn_hits.get(cfg["vn_ch"][1])
        if vn0 is None or vn1 is None:
            continue
        if (vn0["status"] & STATUS_TIMING_SHIFTED) or (vn1["status"] & STATUS_TIMING_SHIFTED):
            vn_st512_cnt += 1

        t_v2, v2_used_512 = station_mean_time(v2_hits, cfg["v2_ch"], 0.0)
        t_ps1, _ = station_mean_time(ps_hits, cfg["ps_dig1_ch"], 0.0)
        t_ps2, ps2_used_512 = station_mean_time(ps_hits, cfg["ps_dig2_ch"], ddig2) if cfg["ps_dig2_ch"] else (float("nan"), False)

        # GENIE neutrino backsplash study
        if cfg["is_genie"]:
            # Require downstream activity (Veto2 raw > 40 or Preshower Dig1 raw > 2.5) and valid timing
            v20 = v2_hits.get(cfg["v2_ch"][0])
            v21 = v2_hits.get(cfg["v2_ch"][1])
            ps0 = ps_hits.get(cfg["ps_dig1_ch"][0])
            ps1 = ps_hits.get(cfg["ps_dig1_ch"][1])
            has_ds = (
                (v20 and v21 and v20["raw_charge"] > 40.0 and v21["raw_charge"] > 40.0)
                and (ps0 and ps1 and ps0["raw_charge"] > 2.5 and ps1["raw_charge"] > 2.5)
            )
            if not has_ds:
                continue
            if not math.isnan(t_v2):
                exp_t = t_v2 - dv
                fb_src = 1
            elif not math.isnan(t_ps1):
                exp_t = t_ps1 - dps
                fb_src = 2
            else:
                continue
            if not (vn0["has_wf"] and vn1["has_wf"]):
                continue
            n_genie_sel += 1
            qrmax = max(vn0["raw_charge"], vn1["raw_charge"])
            genie_raw_max.append(qrmax)
            tv0, wv0 = np.asarray(vn0["tv"]), np.asarray(vn0["wv"])
            tv1, wv1 = np.asarray(vn1["tv"]), np.asarray(vn1["wv"])
            sh0 = 16.0 if (vn0["status"] & STATUS_TIMING_SHIFTED) else 0.0
            sh1 = 16.0 if (vn1["status"] & STATUS_TIMING_SHIFTED) else 0.0
            rel0 = tv0 + sh0 - exp_t
            rel1 = tv1 + sh1 - exp_t
            for dt in DT_SCAN:
                q0 = float(wv0[rel0 < dt].sum()) * INT_TO_CHARGE
                q1 = float(wv1[rel1 < dt].sum()) * INT_TO_CHARGE
                genie_red_max_by_dt[dt].append(max(q0, q1))
                q0_2s = float(wv0[(rel0 >= -20.0) & (rel0 < dt)].sum()) * INT_TO_CHARGE
                q1_2s = float(wv1[(rel1 >= -20.0) & (rel1 < dt)].sum()) * INT_TO_CHARGE
                genie_red_twosided_by_dt[dt].append(max(q0_2s, q1_2s))
            if qrmax >= 40.0 and (vn0["status"] == 0 or vn1["status"] == 0):
                vn_dom = vn0 if vn0["raw_charge"] >= vn1["raw_charge"] else vn1
                if vn_dom["status"] == 0:
                    genie_dt_late.append(vn_dom["localtime"] - exp_t)
                # Candidate representative Fig. 12 neutrino backsplash event (with Veto2 + Preshower waveforms saved!)
                if (
                    fb_src == 1
                    and v20 and v20["has_wf"] and v20["status"] == 0
                    and ps0 and ps0["has_wf"] and ps0["status"] == 0
                    and vn_dom["status"] == 0
                    and 8.0 <= (vn_dom["localtime"] - exp_t) <= 60.0
                    and 50.0 <= vn_dom["raw_charge"] <= 400.0
                ):
                    q35 = float(np.asarray(vn_dom["wv"])[np.asarray(vn_dom["tv"]) - exp_t < 35.0].sum()) * INT_TO_CHARGE
                    # Score prefers clean delayed backsplash (q35 < 25 pC, delay ~15-35 ns)
                    sc = abs((vn_dom["localtime"] - exp_t) - 22.0) + (0.0 if q35 < 25.0 else 50.0)
                    if sc < best_genie_score:
                        best_genie_score = sc
                        best_genie_evt = {
                            "entry": i,
                            "exp_t_v2": exp_t,
                            "exp_t_ps1": t_ps1 - dps if not math.isnan(t_ps1) else None,
                            "v2_mean_t": t_v2,
                            "ps1_mean_t": t_ps1 if not math.isnan(t_ps1) else None,
                            "q_red35_vn": q35,
                            "vn": vn_dom,
                            "v2": v20,
                            "ps1": ps0,
                            "ps2": None,
                        }
            continue

        # Throughgoing muon selection (VetoNu0/1 > 40 pC, status == 0, valid Veto2 > 40 pC)
        v20 = v2_hits.get(cfg["v2_ch"][0])
        v21 = v2_hits.get(cfg["v2_ch"][1])
        ps0 = ps_hits.get(cfg["ps_dig1_ch"][0])
        ps1 = ps_hits.get(cfg["ps_dig1_ch"][1])
        ps2 = ps_hits.get(cfg["ps_dig2_ch"][0]) if cfg["ps_dig2_ch"] else None

        is_muon = (
            (vn0["status"] & ~STATUS_TIMING_SHIFTED) == 0
            and (vn1["status"] & ~STATUS_TIMING_SHIFTED) == 0
            and vn0["raw_charge"] > 40.0 and vn1["raw_charge"] > 40.0
            and not math.isnan(t_v2)
            and v20 is not None and v21 is not None
            and v20["raw_charge"] > 40.0 and v21["raw_charge"] > 40.0
            and vn0["has_wf"]
        )
        if not is_muon:
            continue

        n_muons += 1
        t_vn = 0.5 * (vn0["localtime"] + vn1["localtime"])
        exp_t_v2 = t_v2 - dv
        dt_v2_meas.append(t_v2 - t_vn)

        tv0 = np.asarray(vn0["tv"], dtype=np.float64)
        wv0 = np.asarray(vn0["wv"], dtype=np.float64)
        sh0 = 16.0 if (vn0["status"] & STATUS_TIMING_SHIFTED) else 0.0
        rel0 = tv0 + sh0 - exp_t_v2

        vn_t_start_list.append(float(tv0[0]))
        vn_t_end_list.append(float(tv0[-1]))
        vn_pre_len_list.append(float(exp_t_v2 - tv0[0]))
        vn_width_list.append(float(vn0["width"]))
        vn_raw_q_list.append(float(vn0["raw_charge"]))

        # Peak time & FWHM of VetoNu0
        ipk = int(np.argmax(wv0))
        pk_val = float(wv0[ipk])
        vn_tpeak_rel_list.append(float(tv0[ipk] - vn0["localtime"]))
        above_half = np.where(wv0 >= 0.5 * pk_val)[0]
        if above_half.size >= 2:
            vn_fwhm_list.append(float(tv0[above_half[-1]] - tv0[above_half[0]]))

        # Pre-pulse charge & post-pulse undershoot
        vn_q_pre15_list.append(float(wv0[rel0 < -15.0].sum()) * INT_TO_CHARGE)
        vn_q_pre20_list.append(float(wv0[rel0 < -20.0].sum()) * INT_TO_CHARGE)
        post_mask = (rel0 >= 20.0) & (rel0 <= 120.0)
        if np.any(post_mask):
            vn_min_post_list.append(float(wv0[post_mask].min()))

        # Integration window scan (Veto2 timing)
        for dt in DT_SCAN:
            q_red_v2_by_dt[dt].append(float(wv0[rel0 < dt].sum()) * INT_TO_CHARGE)
            q_red_v2_twosided_by_dt[dt].append(float(wv0[(rel0 >= -20.0) & (rel0 < dt)].sum()) * INT_TO_CHARGE)

        # Integration window scan (Preshower Dig1 timing)
        if not math.isnan(t_ps1) and ps0 is not None and ps0["raw_charge"] > 2.5:
            exp_t_ps1 = t_ps1 - dps
            dt_ps1_meas.append(t_ps1 - t_vn)
            rel_ps1 = tv0 + sh0 - exp_t_ps1
            q_raw_ps1_list.append(float(vn0["raw_charge"]))
            for dt in DT_SCAN:
                q_red_ps1_by_dt[dt].append(float(wv0[rel_ps1 < dt].sum()) * INT_TO_CHARGE)

        # Integration window scan (Preshower Dig2 aligned timing, 2025-2026)
        if not math.isnan(t_ps2):
            exp_t_ps2 = t_ps2 - dps
            dt_ps2_meas.append(t_ps2 - t_vn)
            rel_ps2 = tv0 + sh0 - exp_t_ps2
            q_raw_ps2_list.append(float(vn0["raw_charge"]))
            for dt in DT_SCAN:
                q_red_ps2_by_dt[dt].append(float(wv0[rel_ps2 < dt].sum()) * INT_TO_CHARGE)

        # Accumulate aligned average VetoNu waveform (for clean MIPs 50-200 pC)
        if 50.0 <= vn0["raw_charge"] <= 200.0 and pk_val > 20.0 and abs(t_vn - exp_t_v2) < 3.0:
            interp_w = np.interp(rel_grid, rel0, wv0, left=np.nan, right=np.nan)
            valid_g = np.isfinite(interp_w)
            avg_vn_sum[valid_g] += interp_w[valid_g]
            avg_vn_norm_sum[valid_g] += interp_w[valid_g] / pk_val
            avg_vn_cnt[valid_g] += 1.0

        # Check if ALL stations have waveforms saved (prescale 1/23 events!)
        has_all_wf = (
            vn0["has_wf"]
            and v20 is not None and v20["has_wf"]
            and ps0 is not None and ps0["has_wf"]
            and (not cfg["ps_dig2_ch"] or (ps2 is not None and ps2["has_wf"]))
        )
        if has_all_wf:
            n_muons_with_all_wf += 1
            # Accumulate aligned Veto2 & PS1 waveforms relative to exp_t_v2!
            if 50.0 <= v20["raw_charge"] <= 200.0 and v20["status"] == 0:
                tv_v2 = np.asarray(v20["tv"], dtype=np.float64)
                wv_v2 = np.asarray(v20["wv"], dtype=np.float64)
                iw = np.interp(rel_grid, tv_v2 - exp_t_v2, wv_v2, left=np.nan, right=np.nan)
                vg = np.isfinite(iw)
                avg_v2_sum[vg] += iw[vg]
                avg_v2_cnt[vg] += 1.0
            if ps0["status"] == 0 and 2.5 <= ps0["raw_charge"] <= 15.0:
                tv_p1 = np.asarray(ps0["tv"], dtype=np.float64)
                wv_p1 = np.asarray(ps0["wv"], dtype=np.float64)
                iw = np.interp(rel_grid, tv_p1 - exp_t_v2, wv_p1, left=np.nan, right=np.nan)
                vg = np.isfinite(iw)
                avg_ps1_sum[vg] += iw[vg]
                avg_ps1_cnt[vg] += 1.0

            # Pick cleanest representative Fig. 12 event (all status==0)
            if (
                vn0["status"] == 0 and v20["status"] == 0 and ps0["status"] == 0
                and 65.0 <= vn0["raw_charge"] <= 140.0
                and 65.0 <= v20["raw_charge"] <= 140.0
                and 3.5 <= ps0["raw_charge"] <= 10.0
                and (not cfg["ps_dig2_ch"] or (ps2["status"] == 0 and 1.2 <= ps2["raw_charge"] <= 6.0))
            ):
                err_v2 = abs((v20["localtime"] - vn0["localtime"]) - dv)
                err_ps = abs((ps0["localtime"] - vn0["localtime"]) - dps)
                err_ps2 = abs((ps2["localtime"] - ddig2 - vn0["localtime"]) - dps) if (cfg["ps_dig2_ch"] and ps2) else 0.0
                sc = err_v2 + 0.6 * err_ps + 0.4 * err_ps2 + 0.005 * abs(vn0["raw_charge"] - 95.0)
                if sc < best_muon_score:
                    best_muon_score = sc
                    best_muon_evt = {
                        "entry": i,
                        "vn": vn0,
                        "v2": v20,
                        "ps1": ps0,
                        "ps2": ps2,
                        "exp_t_v2": v20["localtime"] - dv,
                        "exp_t_ps1": ps0["localtime"] - dps,
                        "exp_t_ps2": (ps2["localtime"] - ddig2 - dps) if (cfg["ps_dig2_ch"] and ps2) else None,
                    }

            # Also pick a clean status==512 event in 2024_duringCaloNu!
            if (
                period == "2024_duringCaloNu"
                and vn0["status"] == 0 and v20["status"] == STATUS_TIMING_SHIFTED and ps0["status"] == 0
                and 65.0 <= vn0["raw_charge"] <= 140.0
                and 65.0 <= v20["raw_charge"] <= 140.0
                and 3.5 <= ps0["raw_charge"] <= 10.0
            ):
                err_v2 = abs((v20["localtime"] + 16.0 - vn0["localtime"]) - dv)
                err_ps = abs((ps0["localtime"] - vn0["localtime"]) - dps)
                sc = err_v2 + 0.6 * err_ps + 0.005 * abs(vn0["raw_charge"] - 95.0)
                if sc < best_512_score:
                    best_512_score = sc
                    best_512_evt = {
                        "entry": i,
                        "vn": vn0,
                        "v2": v20,
                        "ps1": ps0,
                        "ps2": None,
                        "exp_t_v2": v20["localtime"] + 16.0 - dv,
                        "exp_t_ps1": ps0["localtime"] - dps,
                    }

    f.Close()

    def summarize_scan(q_by_dt: dict[float, list[float]], raw_arr: np.ndarray) -> dict[str, dict[str, float]]:
        out = {}
        if raw_arr.size == 0:
            return out
        for dt in DT_SCAN:
            arr = np.asarray(q_by_dt[dt], dtype=np.float64)
            ratio = arr / np.maximum(raw_arr, 1e-3)
            out[f"{dt:.0f}"] = {
                "n": int(arr.size),
                "q_med": round(float(np.median(arr)), 2),
                "q_q05": round(float(np.quantile(arr, 0.05)), 2),
                "cap_med": round(float(np.median(ratio)), 4),
                "cap_q05": round(float(np.quantile(ratio, 0.05)), 4),
                "cap_mean": round(float(np.mean(ratio)), 4),
                "leak_lt30_frac": round(float(np.mean(arr < 30.0)), 5),
                "leak_lt40_frac": round(float(np.mean(arr < 40.0)), 5),
            }
        return out

    raw_v2_arr = np.asarray(vn_raw_q_list, dtype=np.float64)
    raw_ps1_arr = np.asarray(q_raw_ps1_list, dtype=np.float64)
    raw_ps2_arr = np.asarray(q_raw_ps2_list, dtype=np.float64)

    genie_summary = None
    if cfg["is_genie"] and genie_raw_max:
        graw = np.asarray(genie_raw_max, dtype=np.float64)
        bs = graw >= 40.0
        n_bs = int(bs.sum())
        scan_g = {}
        for dt in DT_SCAN:
            gred = np.asarray(genie_red_max_by_dt[dt], dtype=np.float64)
            gred2 = np.asarray(genie_red_twosided_by_dt[dt], dtype=np.float64)
            scan_g[f"{dt:.0f}"] = {
                "rescue_lt30_frac": round(float(np.mean(gred[bs] < 30.0)), 4) if n_bs else 0.0,
                "rescue_lt30_twosided_frac": round(float(np.mean(gred2[bs] < 30.0)), 4) if n_bs else 0.0,
                "n_rescued": int(np.sum(gred[bs] < 30.0)),
                "pass_total_frac": round(float(np.mean(gred < 30.0)), 4),
            }
        genie_summary = {
            "n_sel": n_genie_sel,
            "n_raw_ge40_backsplash": n_bs,
            "raw_lt40_frac": round(float(np.mean(graw < 40.0)), 4),
            "late_delay_med_ns": round(float(np.median(genie_dt_late)), 2) if genie_dt_late else None,
            "scan": scan_g,
        }

    # Build average aligned waveform curves
    vn_avg = np.where(avg_vn_cnt > 10, avg_vn_sum / np.maximum(avg_vn_cnt, 1), np.nan)
    vn_norm_avg = np.where(avg_vn_cnt > 10, avg_vn_norm_sum / np.maximum(avg_vn_cnt, 1), np.nan)
    v2_avg = np.where(avg_v2_cnt > 5, avg_v2_sum / np.maximum(avg_v2_cnt, 1), np.nan)
    ps1_avg = np.where(avg_ps1_cnt > 5, avg_ps1_sum / np.maximum(avg_ps1_cnt, 1), np.nan)

    def clean_list(arr):
        return [None if not np.isfinite(x) else round(float(x), 4) for x in arr]

    return {
        "period": period,
        "label": cfg["label"],
        "preset": cfg["preset"],
        "dv": dv,
        "dps": dps,
        "ddig2": ddig2,
        "v2_ch": cfg["v2_ch"],
        "ps_dig1_ch": cfg["ps_dig1_ch"],
        "ps_dig2_ch": cfg["ps_dig2_ch"],
        "vn_ch": cfg["vn_ch"],
        "xaod": cfg["xaod"],
        "wall_s": round(time.perf_counter() - t0, 2),
        "n_scanned": end_ev - start_ev,
        "n_muons": n_muons,
        "n_muons_with_all_wf": n_muons_with_all_wf,
        "vn_status512_count": vn_st512_cnt,
        "waveform_geometry": {
            "vn_t_start_med": round(float(np.median(vn_t_start_list)), 1) if vn_t_start_list else None,
            "vn_t_end_med": round(float(np.median(vn_t_end_list)), 1) if vn_t_end_list else None,
            "vn_pre_len_med_ns": round(float(np.median(vn_pre_len_list)), 2) if vn_pre_len_list else None,
            "vn_pre_len_q05_ns": round(float(np.quantile(vn_pre_len_list, 0.05)), 2) if vn_pre_len_list else None,
            "vn_width_med_ns": round(float(np.median(vn_width_list)), 2) if vn_width_list else None,
            "vn_fwhm_med_ns": round(float(np.median(vn_fwhm_list)), 2) if vn_fwhm_list else None,
            "vn_tpeak_rel_med_ns": round(float(np.median(vn_tpeak_rel_list)), 2) if vn_tpeak_rel_list else None,
            "vn_q_pre15_med_pC": round(float(np.median(vn_q_pre15_list)), 3) if vn_q_pre15_list else None,
            "vn_q_pre15_q99_pC": round(float(np.quantile(vn_q_pre15_list, 0.99)), 2) if vn_q_pre15_list else None,
            "vn_q_pre15_gt5_frac": round(float(np.mean(np.asarray(vn_q_pre15_list) > 5.0)), 4) if vn_q_pre15_list else None,
            "vn_q_pre20_gt5_frac": round(float(np.mean(np.asarray(vn_q_pre20_list) > 5.0)), 4) if vn_q_pre20_list else None,
            "vn_min_post_med_mV": round(float(np.median(vn_min_post_list)), 2) if vn_min_post_list else None,
            "vn_raw_q_med_pC": round(float(np.median(raw_v2_arr)), 2) if raw_v2_arr.size else None,
        },
        "dt_measured": {
            "v2_vn_med": round(float(np.median(dt_v2_meas)), 3) if dt_v2_meas else None,
            "v2_vn_std": round(float(np.std(dt_v2_meas)), 3) if dt_v2_meas else None,
            "ps1_vn_med": round(float(np.median(dt_ps1_meas)), 3) if dt_ps1_meas else None,
            "ps1_vn_std": round(float(np.std(dt_ps1_meas)), 3) if dt_ps1_meas else None,
            "ps2_vn_aligned_med": round(float(np.median(dt_ps2_meas)), 3) if dt_ps2_meas else None,
            "ps2_vn_aligned_std": round(float(np.std(dt_ps2_meas)), 3) if dt_ps2_meas else None,
        },
        "window_scan_v2": summarize_scan(q_red_v2_by_dt, raw_v2_arr),
        "window_scan_v2_twosided": summarize_scan(q_red_v2_twosided_by_dt, raw_v2_arr),
        "window_scan_ps1": summarize_scan(q_red_ps1_by_dt, raw_ps1_arr),
        "window_scan_ps2": summarize_scan(q_red_ps2_by_dt, raw_ps2_arr),
        "genie_summary": genie_summary,
        "aligned_avg": {
            "rel_t_ns": rel_grid.tolist(),
            "vn_mv": clean_list(vn_avg),
            "vn_norm": clean_list(vn_norm_avg),
            "v2_mv": clean_list(v2_avg),
            "ps1_mv": clean_list(ps1_avg),
            "n_vn_avg": int(np.max(avg_vn_cnt)) if avg_vn_cnt.size else 0,
        },
        "rep_muon_evt": best_muon_evt,
        "rep_512_evt": best_512_evt,
        "rep_genie_evt": best_genie_evt,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--workers", type=int, default=11)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    print(f"Scanning {len(PERIOD_SAMPLES)} period xAODs across {args.workers} workers...", flush=True)
    results = {}
    with cf.ProcessPoolExecutor(max_workers=args.workers) as ex:
        for res in ex.map(process_period, PERIOD_SAMPLES):
            results[res["period"]] = res
            wg = res["waveform_geometry"]
            s35 = res["window_scan_v2"].get("35", {})
            print(
                f"  done {res['period']:20s}: muons={res['n_muons']:5d} (all_wf={res['n_muons_with_all_wf']:4d}) "
                f"t_vn=[{wg['vn_t_start_med']}..{wg['vn_t_end_med']}] pre={wg['vn_pre_len_med_ns']}ns "
                f"fwhm={wg['vn_fwhm_med_ns']}ns cap35={s35.get('cap_med')} leak35={s35.get('leak_lt30_frac')} "
                f"rep_mu={'OK' if res['rep_muon_evt'] else 'NONE'} ({res['wall_s']}s)",
                flush=True,
            )
    print(f"All {len(results)} periods complete in {time.perf_counter() - t0:.1f}s", flush=True)
    (args.outdir / "fig12_waveform_study.json").write_text(json.dumps(results, indent=2))

    # Compact summary without raw waveform arrays
    compact = {}
    for k, v in results.items():
        compact[k] = {kk: vv for kk, vv in v.items() if kk not in ("aligned_avg", "rep_muon_evt", "rep_512_evt", "rep_genie_evt")}
    (args.outdir / "fig12_window_summary.json").write_text(json.dumps(compact, indent=2))
    print(f"Wrote {args.outdir / 'fig12_waveform_study.json'} and fig12_window_summary.json")


if __name__ == "__main__":
    main()
