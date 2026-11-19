#!/usr/bin/env python3
"""Comprehensive FASER_DQ_Analyis audit for the 2025 Preshower GFIT_FAILED (p0013) runs.

Performs two multi-core ROOT RDataFrame passes:
  1. All-2025 survey (all 199 GRL runs in 2025_filelist.txt, 1 valid file/run):
     identifies every run's production tag (p0013 vs p0014), Preshower0/1 GFIT_FAILED rate,
     inferred COOL TriggerOffset, VetoNu/Veto2/Timing/Calo status, and Colliding BCID status.
  2. Full-run FASER_DQ_Analyis DQ audit (all valid PHYS files per run) on:
       - 6 focal p0013 runs: 19415, 19524, 19866, 19923, 19970, 20081
       - 3 overlap runs in BOTH p0013 and p0014: 19336, 19878, 20012
       - 4 clean p0014 reference runs: 20170, 20279, 20481, 20591
     Runs the exact FASER_DQ_RDF.py cutflow, yields/pb^-1 (vs thresholds.json),
     all 7 YAML histogram suites (Track, VetoNu, Veto1/2, Timing, Preshower0-3, CaloHi, CaloLo),
     beam crossingAngle, status bit breakdowns, and p0013 vs p0014 ps_dig1_vetonu contribution.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import json
import os
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

DQ_REPO = Path("/afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis")
GRL_PREFIX = Path("/cvmfs/faser.cern.ch/repo/sw/runlist/v10/faser_runlist_2025_stable")
DQ_THRESHOLDS = Path("/eos/home-x/xcheng/FASER/DQ_internal_note_2022_2026_L10/yield_checks/thresholds.json")
DQ_CONFIG = DQ_REPO / "internal_note" / "config.json"
DEFAULT_OUTDIR = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/"
    "stage04_waveform_reduced_charge/dq_audit_2025_p0013"
)

FOCAL_P0013_RUNS = [19415, 19524, 19866, 19923, 19970, 20081]
OVERLAP_RUNS = [19336, 19878, 20012]
REF_P0014_RUNS = [20170, 20279, 20481, 20591]


def load_grl_2025() -> Tuple[Dict[int, dict], Dict[int, float]]:
    grl_json = json.loads(GRL_PREFIX.with_suffix(".json").read_text())
    lumis: Dict[int, float] = {}
    with GRL_PREFIX.with_suffix(".csv").open() as f:
        for row in csv.DictReader((l for l in f if not l.startswith("#")), skipinitialspace=True):
            lumis[int(row["Run"])] = float(row["LumiRec"])
    return {int(k): v for k, v in grl_json.items()}, lumis


def build_time_cuts(grl_entry: dict, run: int) -> Tuple[str, str]:
    stable = grl_entry.get("stable_list", [])
    good_cut = " || ".join(
        f"((eventTime >= {iv['start_utime']}) && (eventTime <= {iv['stop_utime']}) && (run == {run}))"
        for iv in stable
    )
    excl = grl_entry.get("excluded_list", [])
    excl_cut = " || ".join(
        f"((eventTime >= {iv['start_utime']}) && (eventTime <= {iv['stop_utime']}) && (run == {run}))"
        for iv in excl
    )
    return good_cut or "true", excl_cut


def filter_valid_files(files: List[str], unavailable: set) -> Tuple[List[str], List[str]]:
    valid, bad = [], []
    for uri in files:
        if uri in unavailable:
            bad.append(uri)
            continue
        p = Path(uri.replace("root://eospublic.cern.ch/", ""))
        if not p.is_file() or p.stat().st_size < 10000:
            bad.append(uri)
            continue
        valid.append(uri)
    return valid, bad


def list_tag_files(run: int, tag: str, unavailable: set) -> Tuple[List[str], List[str]]:
    d = Path(f"/eos/experiment/faser/phys/2025/{tag}/{run:06d}")
    if not d.is_dir():
        return [], []
    raw = [f"root://eospublic.cern.ch/{p}" for p in sorted(d.glob("*.root"))]
    return filter_valid_files(raw, unavailable)


def hist_stats(h) -> dict:
    nb = h.GetNbinsX()
    counts = np.array([h.GetBinContent(i) for i in range(1, nb + 1)], dtype=np.float64)
    edges = np.array([h.GetBinLowEdge(i) for i in range(1, nb + 2)], dtype=np.float64)
    centers = 0.5 * (edges[:-1] + edges[1:])
    tot = float(counts.sum())
    under = float(h.GetBinContent(0))
    over = float(h.GetBinContent(nb + 1))
    out = {
        "entries": int(round(h.GetEntries())),
        "in_range": int(round(tot)),
        "underflow": int(round(under)),
        "overflow": int(round(over)),
        "mean": float(h.GetMean()),
        "std": float(h.GetStdDev()),
        "xmin": float(edges[0]),
        "xmax": float(edges[-1]),
        "nbins": int(nb),
        "counts": counts.tolist(),
    }
    if tot > 0:
        cdf = np.cumsum(counts) / tot
        out["q16"] = float(centers[np.searchsorted(cdf, 0.16)])
        out["median"] = float(centers[np.searchsorted(cdf, 0.50)])
        out["q84"] = float(centers[np.searchsorted(cdf, 0.84)])
        imax = int(np.argmax(counts))
        out["mode"] = float(centers[imax])
    return out


# ---------------------------------------------------------------------------
# Worker 1: Fast all-2025 survey (1 file per run)
# ---------------------------------------------------------------------------
def survey_worker(task: dict) -> dict:
    import ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    ROOT.gInterpreter.Declare("using namespace ROOT::VecOps;")
    ROOT.gInterpreter.Declare(f'#include "{DQ_REPO / "RDFDefines.h"}"')

    run = task["run"]
    tag = task["tag"]
    fp = task["file"]
    good_cut = task["good_cut"]
    excl_cut = task["excl_cut"]

    try:
        tf = ROOT.TFile.Open(fp)
        if not tf or tf.IsZombie() or not tf.Get("nt"):
            return {"run": run, "tag": tag, "status": "bad_file", "file": fp}
        tf.Close()

        df = ROOT.RDataFrame("nt", fp)
        cols = set(str(c) for c in df.GetColumnNames())
        n_total = df.Count()
        df_g = df.Filter(good_cut)
        if excl_cut:
            df_g = df_g.Filter(f"!({excl_cut})")
        df_dq = (
            df_g.Filter("(Timing0_status & 4) == 0 && (Timing1_status & 4) == 0 && (Timing2_status & 4) == 0 && (Timing3_status & 4) == 0")
            .Filter("distanceToCollidingBCID == 0")
            .Filter("(TAP & 4) != 0")
        )
        df_dq = (
            df_dq.Define("Track_nHits", "Track_nDoF + 5")
            .Define("Track_chi2_per_dof", "Track_Chi2 / Track_nDoF")
            .Define("GoodTracks", "Track_nLayers >= 7 && Track_chi2_per_dof < 25 && Track_nHits >= 12 && Track_pz0 > 20000 && RemoveDuplicates(Track_p0)")
            .Define("NGoodTracks", "Track_nLayers[GoodTracks].size()")
        )
        n_dq = df_dq.Count()
        df_gt = df_dq.Filter("NGoodTracks >= 1")
        n_gt = df_gt.Count()

        # Check raw colliding & status on goodtracks (or on total if gt==0)
        n_coll = df.Filter("distanceToCollidingBCID == 0").Count()
        h_bcid = df.Histo1D(("bcid", "", 100, -4000, 4000), "distanceToCollidingBCID")

        ps0_s0 = df_gt.Filter("Preshower0_status == 0 || Preshower0_status == 512").Count()
        ps0_gfit = df_gt.Filter("(Preshower0_status & 16) != 0").Count()
        ps0_thr = df_gt.Filter("(Preshower0_status & 1) != 0").Count()
        ps1_s0 = df_gt.Filter("Preshower1_status == 0 || Preshower1_status == 512").Count()
        ps1_gfit = df_gt.Filter("(Preshower1_status & 16) != 0").Count()
        vn0_s0 = df_gt.Filter("(VetoNu0_status == 0 || VetoNu0_status == 512) && VetoNu0_raw_charge > 40").Count()
        v20_s0 = df_gt.Filter("(Veto20_status == 0 || Veto20_status == 512) && Veto20_raw_charge > 40").Count()
        tm0_s0 = df_gt.Filter("(Timing0_status == 0 || Timing0_status == 512) && Timing0_raw_charge > 10").Count()
        cl0_s0 = df_gt.Filter("CaloLo0_status == 0").Count()

        # Offset inferred from localtime - triggertime (always exact even when GFIT failed, since trigger_time = localtime - (200 + offset)!)
        df_gt2 = df_gt.Define("ps0_off_infer", "Preshower0_localtime - Preshower0_triggertime - 200.0")
        h_off = df_gt2.Histo1D(("off", "", 240, -120.0, 120.0), "ps0_off_infer")
        h_ps0_q = df_gt.Histo1D(("ps0_q", "", 150, 0.0, 15.0), "Preshower0_charge")
        h_ps0_rq = df_gt.Histo1D(("ps0_rq", "", 150, 0.0, 15.0), "Preshower0_raw_charge")
        h_ps0_lt = df_gt.Filter("(Preshower0_status == 0 || Preshower0_status == 512) && Preshower0_raw_charge > 2.5").Histo1D(("ps0_lt", "", 200, 100.0, 300.0), "Preshower0_localtime")

        has_ps2 = "Preshower2_status" in cols
        ps2_gfit = df_gt.Filter("(Preshower2_status & 16) != 0").Count() if has_ps2 else None
        ps2_ok = df_gt.Filter("Preshower2_status == 0 || Preshower2_status == 512").Count() if has_ps2 else None

        ntot = int(n_total.GetValue())
        ndq = int(n_dq.GetValue())
        ngt = int(n_gt.GetValue())
        denom = max(1, ngt)
        return {
            "run": run,
            "tag": tag,
            "status": "ok",
            "n_files_total": task["n_files_total"],
            "lumi_pb": task["lumi_pb"],
            "n_total": ntot,
            "n_coll_raw": int(n_coll.GetValue()),
            "bcid_mean": float(h_bcid.GetMean()),
            "n_dq": ndq,
            "n_goodtrk": ngt,
            "ps0_s0_frac": float(ps0_s0.GetValue()) / denom if ngt else 0.0,
            "ps0_gfit_frac": float(ps0_gfit.GetValue()) / denom if ngt else 0.0,
            "ps0_thr_frac": float(ps0_thr.GetValue()) / denom if ngt else 0.0,
            "ps1_s0_frac": float(ps1_s0.GetValue()) / denom if ngt else 0.0,
            "ps1_gfit_frac": float(ps1_gfit.GetValue()) / denom if ngt else 0.0,
            "ps2_ok_frac": (float(ps2_ok.GetValue()) / denom if (has_ps2 and ngt) else None),
            "ps2_gfit_frac": (float(ps2_gfit.GetValue()) / denom if (has_ps2 and ngt) else None),
            "vn0_ok_frac": float(vn0_s0.GetValue()) / denom if ngt else 0.0,
            "v20_ok_frac": float(v20_s0.GetValue()) / denom if ngt else 0.0,
            "tm0_ok_frac": float(tm0_s0.GetValue()) / denom if ngt else 0.0,
            "cl0_ok_frac": float(cl0_s0.GetValue()) / denom if ngt else 0.0,
            "ps0_trigger_offset_ns": float(h_off.GetMean()) if ngt else None,
            "ps0_charge_mean": float(h_ps0_q.GetMean()) if ngt else None,
            "ps0_raw_charge_mean": float(h_ps0_rq.GetMean()) if ngt else None,
            "ps0_lt_ok_mean": float(h_ps0_lt.GetMean()) if h_ps0_lt.GetEntries() > 0 else None,
        }
    except Exception as e:
        return {"run": run, "tag": tag, "status": f"error: {e}", "file": fp}


# ---------------------------------------------------------------------------
# Worker 2: Full-run FASER_DQ_Analyis complete DQ audit
# ---------------------------------------------------------------------------
def full_dq_worker(task: dict) -> dict:
    import ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    ROOT.ROOT.EnableImplicitMT(task["threads"])
    ROOT.gInterpreter.Declare("using namespace ROOT::VecOps;")
    ROOT.gInterpreter.Declare(f'#include "{DQ_REPO / "RDFDefines.h"}"')

    run = task["run"]
    tag = task["tag"]
    files = task["files"]
    lumi_pb = task["lumi_pb"]
    good_cut = task["good_cut"]
    excl_cut = task["excl_cut"]

    t0 = time.perf_counter()
    df = ROOT.RDataFrame("nt", files)
    cols = set(str(c) for c in df.GetColumnNames())

    # Exact FASER_DQ_RDF.py cuts in exact order
    df = df.Define("GoodTimes", good_cut).Filter("GoodTimes", "Good times")
    if excl_cut:
        df = df.Define("ExcludedTimes", excl_cut).Filter("!ExcludedTimes", "Excluded times")
    df = df.Filter(
        "(Timing0_status & 4) == 0 && (Timing1_status & 4) == 0 && (Timing2_status & 4) == 0 && (Timing3_status & 4) == 0",
        "No timing saturation",
    )
    df = df.Filter("distanceToCollidingBCID == 0", "Colliding")
    df = df.Filter("(TAP & 4) != 0", "Timing Trigger")

    # Exact FASER_DQ_RDF.py definitions
    df = (
        df.Define("NTracks", "Track_nLayers.size()")
        .Define("NPosTracks", "Track_nLayers[Track_charge > 0].size()")
        .Define("NNegTracks", "Track_nLayers[Track_charge < 0].size()")
        .Define("Track_nHits", "Track_nDoF + 5")
        .Define("Track_chi2_per_dof", "Track_Chi2 / Track_nDoF")
        .Define("GoodTracks", "Track_nLayers >= 7 && Track_chi2_per_dof < 25 && Track_nHits >= 12 && Track_pz0 > 20000 && RemoveDuplicates(Track_p0)")
        .Define("NGoodTracks", "Track_nLayers[GoodTracks].size()")
        .Define("NGoodPosTracks", "Track_nLayers[GoodTracks && Track_charge > 0].size()")
        .Define("NGoodNegTracks", "Track_nLayers[GoodTracks && Track_charge < 0].size()")
        .Define("Track_p0_GeV", "Track_p0 * 0.001")
        .Define("Track_pz_charge0_GeV", "Track_pz0 * Track_charge * 0.001")
        .Define("Track_theta_x0", "asin(Track_px0 / Track_p0)")
        .Define("Track_theta_y0", "asin(Track_py0 / Track_p0)")
        .Define("Track_theta_x0_pos", "Track_theta_x0[Track_charge > 0]")
        .Define("Track_theta_x0_neg", "Track_theta_x0[Track_charge < 0]")
        .Define("Track_theta_y0_pos", "Track_theta_y0[Track_charge > 0]")
        .Define("Track_theta_y0_neg", "Track_theta_y0[Track_charge < 0]")
        .Define("Track_x0_pos", "Track_x0[Track_charge > 0]")
        .Define("Track_x0_neg", "Track_x0[Track_charge < 0]")
        .Define("Track_y0_pos", "Track_y0[Track_charge > 0]")
        .Define("Track_y0_neg", "Track_y0[Track_charge < 0]")
        .Define("Track_X_atTrig_valid", "Track_X_atTrig[abs(Track_X_atTrig) < 500]")
        .Define("Track_Y_atTrig_valid", "Track_Y_atTrig[abs(Track_Y_atTrig) < 500]")
        .Define("CaloHiYield", "int(CaloHi_total_E_EM > 100000.0)")
        .Define("CaloLoYield", "int(CaloLo_total_E_EM > 10000.0 && CaloLo_total_E_EM < 100000.0)")
        .Define("CaloHi_total_E_EM_GeV", "CaloHi_total_E_EM * 0.001")
        .Define("CaloLo_total_E_EM_GeV", "CaloLo_total_E_EM * 0.001")
        .Define("CaloHi_total_raw_E_EM_GeV", "CaloHi_total_raw_E_EM * 0.001")
        .Define("CaloLo_total_raw_E_EM_GeV", "CaloLo_total_raw_E_EM * 0.001")
    )

    # Yields & cutflow
    cutflow_rep = df.Report()
    yields_ptr = {
        "selected_events": df.Count(),
        "TrkYield": df.Sum("NTracks"),
        "PosTrkYield": df.Sum("NPosTracks"),
        "NegTrkYield": df.Sum("NNegTracks"),
        "GoodTrkYield": df.Sum("NGoodTracks"),
        "GoodPosTrkYield": df.Sum("NGoodPosTracks"),
        "GoodNegTrkYield": df.Sum("NGoodNegTracks"),
        "CaloHiYield": df.Sum("CaloHiYield"),
        "CaloLoYield": df.Sum("CaloLoYield"),
        "NoTrackEvents": df.Filter("NTracks == 0").Count(),
        "GoodTrackEvents": df.Filter("NGoodTracks >= 1").Count(),
    }

    # DQ histograms (matching FASER_DQ_RDF.py YAML configs + Fig. 5/6 + timing_scan + beam_scan)
    hists = {}
    # 1. Track histograms (Fig. 6 / track_hists.yaml)
    hists["NTracks"] = df.Histo1D(("NTracks", "", 15, 0, 15), "NTracks")
    hists["NGoodTracks"] = df.Histo1D(("NGoodTracks", "", 15, 0, 15), "NGoodTracks")
    hists["Track_nLayers"] = df.Histo1D(("Track_nLayers", "", 10, 0, 10), "Track_nLayers")
    hists["Track_nHits"] = df.Histo1D(("Track_nHits", "", 20, 0, 20), "Track_nHits")
    hists["Track_Chi2"] = df.Histo1D(("Track_Chi2", "", 50, 0, 50), "Track_Chi2")
    hists["Track_nDoF"] = df.Histo1D(("Track_nDoF", "", 20, 0, 20), "Track_nDoF")
    hists["Track_chi2_per_dof"] = df.Histo1D(("Track_chi2_per_dof", "", 50, 0, 25), "Track_chi2_per_dof")
    hists["Track_p0_GeV"] = df.Histo1D(("Track_p0_GeV", "", 100, 0, 2000), "Track_p0_GeV")
    hists["Track_pz_charge0_GeV"] = df.Histo1D(("Track_pz_charge0_GeV", "", 100, -2000, 2000), "Track_pz_charge0_GeV")
    hists["Track_x0"] = df.Histo1D(("Track_x0", "", 50, -125, 125), "Track_x0")
    hists["Track_y0"] = df.Histo1D(("Track_y0", "", 50, -125, 125), "Track_y0")
    hists["Track_theta_x0"] = df.Histo1D(("Track_theta_x0", "", 50, -0.005, 0.005), "Track_theta_x0")
    hists["Track_theta_y0"] = df.Histo1D(("Track_theta_y0", "", 50, -0.005, 0.005), "Track_theta_y0")
    hists["Track_X_atTrig"] = df.Histo1D(("Track_X_atTrig", "", 100, -200, 200), "Track_X_atTrig_valid")
    hists["Track_Y_atTrig"] = df.Histo1D(("Track_Y_atTrig", "", 100, -200, 200), "Track_Y_atTrig_valid")
    if "crossingAngle" in cols:
        hists["crossingAngle"] = df.Histo1D(("crossingAngle", "", 100, -250, 250), "crossingAngle")

    # 2. Scintillator charges (exact Fig. 5 / YAML ranges) + raw charges + baselines
    for ch, xmax in [
        ("VetoNu0", 300.0), ("VetoNu1", 300.0),
        ("Veto10", 300.0), ("Veto11", 300.0),
        ("Veto20", 300.0), ("Veto21", 300.0),
        ("Timing0", 300.0), ("Timing1", 300.0), ("Timing2", 300.0), ("Timing3", 300.0),
        ("Preshower0", 15.0), ("Preshower1", 15.0),
    ]:
        hists[f"{ch}_charge"] = df.Histo1D((f"{ch}_charge", "", 50, 0.01, xmax), f"{ch}_charge")
        hists[f"{ch}_raw_charge"] = df.Histo1D((f"{ch}_raw_charge", "", 50, 0.01, xmax), f"{ch}_raw_charge")
        hists[f"{ch}_baseline"] = df.Histo1D((f"{ch}_baseline", "", 200, 1850.0, 1950.0), f"{ch}_baseline")
        hists[f"{ch}_baseline_rms"] = df.Histo1D((f"{ch}_baseline_rms", "", 100, 0.0, 5.0), f"{ch}_baseline_rms")

    # Preshower2/3 (Dig2, low-gain)
    for ch in ["Preshower2", "Preshower3"]:
        if f"{ch}_charge" in cols:
            hists[f"{ch}_charge"] = df.Histo1D((f"{ch}_charge", "", 50, 0.01, 15.0), f"{ch}_charge")
            hists[f"{ch}_raw_charge"] = df.Histo1D((f"{ch}_raw_charge", "", 50, 0.01, 15.0), f"{ch}_raw_charge")
            hists[f"{ch}_baseline"] = df.Histo1D((f"{ch}_baseline", "", 200, 1850.0, 1950.0), f"{ch}_baseline")

    # Timing peak/width/charge-windows (stage10/11)
    for ch in ["Timing0", "Timing1", "Timing2", "Timing3"]:
        hists[f"{ch}_peak"] = df.Histo1D((f"{ch}_peak", "", 100, 0.0, 1000.0), f"{ch}_peak")
        hists[f"{ch}_width"] = df.Histo1D((f"{ch}_width", "", 100, 0.0, 20.0), f"{ch}_width")
        hists[f"{ch}_charge_0_80"] = df.Histo1D((f"{ch}_charge_0_80", "", 80, 0.0, 80.0), f"{ch}_charge")

    # Preshower0..3 peak/raw_peak/width/localtime/triggertime/bcidtime (exact preshower_hists.yaml + full range)
    for ch in ["Preshower0", "Preshower1"]:
        df_ps = df.Filter(f"{ch}_raw_charge > 2.5")
        hists[f"{ch}_triggertime_yaml"] = df_ps.Histo1D((f"{ch}_tt", "", 100, -75.0, 50.0), f"{ch}_triggertime")
        hists[f"{ch}_localtime_yaml"] = df_ps.Histo1D((f"{ch}_lt", "", 100, 150.0, 225.0), f"{ch}_localtime")
        hists[f"{ch}_localtime_wide"] = df_ps.Histo1D((f"{ch}_ltw", "", 300, -10.0, 290.0), f"{ch}_localtime")
        hists[f"{ch}_triggertime_wide"] = df_ps.Histo1D((f"{ch}_ttw", "", 300, -250.0, 50.0), f"{ch}_triggertime")
        hists[f"{ch}_bcidtime_yaml"] = df_ps.Histo1D((f"{ch}_bt", "", 100, -25.0, 35.0), f"{ch}_bcidtime")
        hists[f"{ch}_raw_peak"] = df_ps.Histo1D((f"{ch}_rp", "", 100, 0.0, 100.0), f"{ch}_raw_peak")
        hists[f"{ch}_peak"] = df_ps.Histo1D((f"{ch}_pk", "", 100, 0.0, 100.0), f"{ch}_peak")
        hists[f"{ch}_width"] = df_ps.Histo1D((f"{ch}_wd", "", 100, 0.0, 20.0), f"{ch}_width")

    # 3. Calorimeter Hi & Lo (caloHi_hists.yaml, caloLo_hists.yaml)
    hists["CaloHi_total_E_EM_GeV"] = df.Histo1D(("CaloHi_total_E_EM_GeV", "", 100, 100.0, 2000.0), "CaloHi_total_E_EM_GeV")
    hists["CaloLo_total_E_EM_GeV"] = df.Histo1D(("CaloLo_total_E_EM_GeV", "", 100, 0.0, 200.0), "CaloLo_total_E_EM_GeV")
    hists["CaloHi_total_nMIP"] = df.Histo1D(("CaloHi_total_nMIP", "", 100, 10.0, 3000.0), "CaloHi_total_nMIP")
    hists["CaloLo_total_nMIP"] = df.Histo1D(("CaloLo_total_nMIP", "", 100, 0.0, 200.0), "CaloLo_total_nMIP")
    for ch in ["CaloHi0", "CaloHi1", "CaloHi2", "CaloHi3"]:
        hists[f"{ch}_charge"] = df.Histo1D((f"{ch}_charge", "", 100, 0.01, 4.0), f"{ch}_charge")
        hists[f"{ch}_localtime"] = df.Filter(f"{ch}_status == 0 && CaloLo{ch[-1]}_raw_charge > 10 * {ch}_raw_charge").Histo1D((f"{ch}_lt", "", 100, 150.0, 225.0), f"{ch}_localtime")
    for ch in ["CaloLo0", "CaloLo1", "CaloLo2", "CaloLo3"]:
        hists[f"{ch}_charge"] = df.Histo1D((f"{ch}_charge", "", 100, 0.01, 100.0), f"{ch}_charge")
        hists[f"{ch}_localtime"] = df.Filter(f"{ch}_status == 0").Histo1D((f"{ch}_lt", "", 100, 150.0, 225.0), f"{ch}_localtime")

    # 4. Status bit breakdown across EVERY scintillator & calorimeter channel (on GoodTracks >= 1)
    df_gt = df.Filter("NGoodTracks >= 1")
    status_channels = [
        "VetoNu0", "VetoNu1", "Veto10", "Veto11", "Veto20", "Veto21",
        "Timing0", "Timing1", "Timing2", "Timing3",
        "Preshower0", "Preshower1", "Preshower2", "Preshower3",
        "CaloLo0", "CaloLo1", "CaloLo2", "CaloLo3",
        "CaloHi0", "CaloHi1", "CaloHi2", "CaloHi3",
    ]
    status_ptrs = {}
    for ch in status_channels:
        if f"{ch}_status" not in cols:
            continue
        status_ptrs[ch] = {
            "ok": df_gt.Filter(f"{ch}_status == 0 || {ch}_status == 512").Count(),
            "thr": df_gt.Filter(f"({ch}_status & 1) != 0").Count(),
            "ovf": df_gt.Filter(f"({ch}_status & 4) != 0").Count(),
            "gfit": df_gt.Filter(f"({ch}_status & 16) != 0").Count(),
            "cfit": df_gt.Filter(f"({ch}_status & 32) != 0").Count(),
            "t512": df_gt.Filter(f"({ch}_status & 512) != 0").Count(),
        }

    # 5. Exact Fig. 13 timing differences on GoodTracks >= 1
    df_t = (
        df_gt
        .Define("vn0_ok", "VetoNu0_raw_charge > 40 && (VetoNu0_status == 0 || VetoNu0_status == 512)")
        .Define("vn1_ok", "VetoNu1_raw_charge > 40 && (VetoNu1_status == 0 || VetoNu1_status == 512)")
        .Define("vn_n", "int(vn0_ok) + int(vn1_ok)")
        .Define("vn_t", "vn_n > 0 ? ((vn0_ok ? VetoNu0_localtime : 0.0) + (vn1_ok ? VetoNu1_localtime : 0.0)) / double(vn_n) : -9999.0")
        .Define("v20_ok", "Veto20_raw_charge > 40 && (Veto20_status == 0 || Veto20_status == 512)")
        .Define("v21_ok", "Veto21_raw_charge > 40 && (Veto21_status == 0 || Veto21_status == 512)")
        .Define("v2_n", "int(v20_ok) + int(v21_ok)")
        .Define("v2_t", "v2_n > 0 ? ((v20_ok ? (Veto20_localtime + ((Veto20_status & 512) ? 16.0 : 0.0)) : 0.0) + (v21_ok ? (Veto21_localtime + ((Veto21_status & 512) ? 16.0 : 0.0)) : 0.0)) / double(v2_n) : -9999.0")
        .Define("ps0_ok", "Preshower0_raw_charge > 2.5 && (Preshower0_status == 0 || Preshower0_status == 512)")
        .Define("ps1_ok", "Preshower1_raw_charge > 2.5 && (Preshower1_status == 0 || Preshower1_status == 512)")
        .Define("ps12_n", "int(ps0_ok) + int(ps1_ok)")
        .Define("ps12_t", "ps12_n > 0 ? ((ps0_ok ? Preshower0_localtime : 0.0) + (ps1_ok ? Preshower1_localtime : 0.0)) / double(ps12_n) : -9999.0")
        .Define("ps2_ok", "Preshower2_status == 0 || Preshower2_status == 512")
        .Define("ps3_ok", "Preshower3_status == 0 || Preshower3_status == 512")
        .Define("ps34_n", "int(ps2_ok) + int(ps3_ok)")
        .Define("ps34_t", "ps34_n > 0 ? ((ps2_ok ? (Preshower2_localtime + ((Preshower2_status & 512) ? 16.0 : 0.0)) : 0.0) + (ps3_ok ? (Preshower3_localtime + ((Preshower3_status & 512) ? 16.0 : 0.0)) : 0.0)) / double(ps34_n) : -9999.0")
        .Define("dt_v2_vn", "v2_t - vn_t")
        .Define("dt_ps1_vn", "ps12_t - vn_t")
        .Define("dt_ps2_vn", "ps34_t - vn_t")
        .Define("dt_ps2_ps1", "ps34_t - ps12_t")
        .Define("ps0_off_infer", "Preshower0_localtime - Preshower0_triggertime - 200.0")
        .Define("ps2_off_infer", "Preshower2_localtime - Preshower2_triggertime - 200.0 + ((Preshower2_status & 512) ? 16.0 : 0.0)")
    )
    hists["dt_v2_vn"] = df_t.Filter("vn_n > 0 && v2_n > 0").Histo1D(("dt_v2_vn", "", 120, -60.0, 60.0), "dt_v2_vn")
    hists["dt_ps1_vn"] = df_t.Filter("vn_n > 0 && ps12_n > 0").Histo1D(("dt_ps1_vn", "", 120, -60.0, 60.0), "dt_ps1_vn")
    hists["dt_ps2_vn"] = df_t.Filter("vn_n > 0 && ps34_n > 0").Histo1D(("dt_ps2_vn", "", 120, -60.0, 60.0), "dt_ps2_vn")
    hists["dt_ps2_ps1"] = df_t.Filter("ps12_n > 0 && ps34_n > 0").Histo1D(("dt_ps2_ps1", "", 120, -60.0, 60.0), "dt_ps2_ps1")
    hists["ps0_off_infer"] = df_t.Histo1D(("ps0_off_infer", "", 240, -120.0, 120.0), "ps0_off_infer")
    hists["ps2_off_infer"] = df_t.Histo1D(("ps2_off_infer", "", 240, -120.0, 120.0), "ps2_off_infer")

    # Materialize everything in ONE event loop!
    cutflow = {}
    for c in cutflow_rep:
        cutflow[c.GetName()] = {"all": int(c.GetAll()), "passed": int(c.GetPass())}

    yields = {k: int(round(float(v.GetValue()))) for k, v in yields_ptr.items()}
    ngt = max(1, yields["GoodTrackEvents"])

    status_summary = {}
    for ch, d_ch in status_ptrs.items():
        status_summary[ch] = {k: int(v.GetValue()) for k, v in d_ch.items()}
        status_summary[ch]["ok_pct"] = round(100.0 * status_summary[ch]["ok"] / ngt, 3)
        status_summary[ch]["gfit_pct"] = round(100.0 * status_summary[ch]["gfit"] / ngt, 3)
        status_summary[ch]["thr_pct"] = round(100.0 * status_summary[ch]["thr"] / ngt, 3)
        status_summary[ch]["ovf_pct"] = round(100.0 * status_summary[ch]["ovf"] / ngt, 3)
        status_summary[ch]["t512_pct"] = round(100.0 * status_summary[ch]["t512"] / ngt, 3)

    hist_out = {k: hist_stats(h.GetValue()) for k, h in hists.items()}

    # Yield normalization & DQ threshold check
    ev_per_pb = yields["selected_events"] / lumi_pb if lumi_pb > 0 else 0.0
    trk_per_pb = yields["GoodTrkYield"] / lumi_pb if lumi_pb > 0 else 0.0
    calo_hi_per_pb = yields["CaloHiYield"] / lumi_pb if lumi_pb > 0 else 0.0
    calo_lo_per_pb = yields["CaloLoYield"] / lumi_pb if lumi_pb > 0 else 0.0
    pos_neg_ratio = yields["GoodPosTrkYield"] / max(1, yields["GoodNegTrkYield"])

    return {
        "key": f"{run}_{tag}",
        "run": run,
        "tag": tag,
        "n_files_valid": len(files),
        "n_files_bad": len(task["bad_files"]),
        "bad_files": task["bad_files"],
        "lumi_pb": lumi_pb,
        "wall_s": round(time.perf_counter() - t0, 2),
        "event_loops": int(df.GetNRuns()),
        "cutflow": cutflow,
        "yields": yields,
        "rates_per_pb": {
            "events": round(ev_per_pb, 2),
            "good_tracks": round(trk_per_pb, 2),
            "calo_hi_gt100GeV": round(calo_hi_per_pb, 3),
            "calo_lo_10_100GeV": round(calo_lo_per_pb, 3),
            "pos_neg_ratio": round(pos_neg_ratio, 4),
            "goodtrk_per_selected_evt": round(yields["GoodTrkYield"] / max(1, yields["selected_events"]), 4),
        },
        "status_summary": status_summary,
        "histograms": hist_out,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--audit-workers", type=int, default=8)
    ap.add_argument("--threads-per-audit", type=int, default=3)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    grl_json, lumis = load_grl_2025()
    dq_cfg = json.loads(DQ_CONFIG.read_text())
    unavailable = set(dq_cfg.get("unavailable_inputs", {}).keys())
    thresholds = json.loads(DQ_THRESHOLDS.read_text())["2025"]

    # Parse 2025_filelist.txt
    fl_lines = (DQ_REPO / "faser_filelists" / "2025_filelist.txt").read_text().splitlines()
    run_files: Dict[int, List[str]] = defaultdict(list)
    run_tag: Dict[int, str] = {}
    for line in fl_lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = re.search(r"/phys/2025/([^/]+)/(\d+)/", line)
        if m:
            r = int(m.group(2))
            run_tag[r] = m.group(1)
            run_files[r].append(line)

    # ------------------------------------------------------------------
    # PASS 1: All-2025 run-by-run survey (1 valid file/run)
    # ------------------------------------------------------------------
    t0 = time.perf_counter()
    survey_tasks = []
    for run in sorted(run_files):
        if run not in grl_json:
            continue
        valid, _ = filter_valid_files(run_files[run], unavailable)
        if not valid:
            continue
        good_cut, excl_cut = build_time_cuts(grl_json[run], run)
        survey_tasks.append({
            "run": run,
            "tag": run_tag[run],
            "file": valid[0],
            "n_files_total": len(run_files[run]),
            "lumi_pb": lumis.get(run, 0.0),
            "good_cut": good_cut,
            "excl_cut": excl_cut,
        })

    print(f"[Pass 1] Surveying ALL {len(survey_tasks)} 2025 GRL runs ({args.workers} workers)...", flush=True)
    survey_results = []
    with cf.ProcessPoolExecutor(max_workers=args.workers) as ex:
        for res in ex.map(survey_worker, survey_tasks):
            survey_results.append(res)
    survey_results.sort(key=lambda x: x["run"])
    print(f"[Pass 1] Done in {time.perf_counter() - t0:.1f}s", flush=True)

    # Print Pass 1 high-level summary
    p13_ok = [r for r in survey_results if r["status"] == "ok" and r["tag"] == "p0013" and r["n_goodtrk"] > 0]
    p14_ok = [r for r in survey_results if r["status"] == "ok" and r["tag"] == "p0014" and r["n_goodtrk"] > 0]
    p13_gfit_hi = [r for r in p13_ok if r["ps0_gfit_frac"] > 0.5]
    p14_gfit_hi = [r for r in p14_ok if r["ps0_gfit_frac"] > 0.5]
    zero_coll = [r for r in survey_results if r["status"] == "ok" and r["n_goodtrk"] == 0]
    vn_off = [r for r in survey_results if r["status"] == "ok" and r["n_goodtrk"] > 0 and r["vn0_ok_frac"] < 0.1]
    print(f"  p0013 non-zero runs: {len(p13_ok)} (GFIT_FAILED > 50%: {len(p13_gfit_hi)}, min={min((r['ps0_gfit_frac'] for r in p13_ok), default=0)*100:.1f}%, max={max((r['ps0_gfit_frac'] for r in p13_ok), default=0)*100:.1f}%)")
    print(f"  p0014 non-zero runs: {len(p14_ok)} (GFIT_FAILED > 50%: {len(p14_gfit_hi)}, max={max((r['ps0_gfit_frac'] for r in p14_ok), default=0)*100:.2f}%)")
    print(f"  zero-GoodTrack runs: {len(zero_coll)} -> {[r['run'] for r in zero_coll]}")
    print(f"  VetoNu0 < 10% ok runs: {[(r['run'], r['tag'], round(r['vn0_ok_frac']*100,2)) for r in vn_off]}")

    # ------------------------------------------------------------------
    # PASS 2: Full-run complete FASER_DQ_Analyis audit on focal + overlap + ref runs
    # ------------------------------------------------------------------
    t1 = time.perf_counter()
    audit_targets = (
        [(r, "p0013") for r in FOCAL_P0013_RUNS]
        + [(r, "p0013") for r in OVERLAP_RUNS]
        + [(r, "p0014") for r in OVERLAP_RUNS]
        + [(r, "p0014") for r in REF_P0014_RUNS]
    )
    audit_tasks = []
    for run, tag in audit_targets:
        valid, bad = list_tag_files(run, tag, unavailable)
        good_cut, excl_cut = build_time_cuts(grl_json[run], run)
        audit_tasks.append({
            "run": run,
            "tag": tag,
            "files": valid,
            "bad_files": bad,
            "lumi_pb": lumis.get(run, 0.0),
            "good_cut": good_cut,
            "excl_cut": excl_cut,
            "threads": args.threads_per_audit,
        })

    print(f"\n[Pass 2] Full-run FASER_DQ_Analyis audit on {len(audit_tasks)} (run, tag) datasets ({args.audit_workers} workers x {args.threads_per_audit} threads)...", flush=True)
    audit_results = {}
    with cf.ProcessPoolExecutor(max_workers=args.audit_workers) as ex:
        for res in ex.map(full_dq_worker, audit_tasks):
            # Apply FASER_DQ_Analyis yield threshold flags
            ev_r = res["rates_per_pb"]["events"]
            tr_r = res["rates_per_pb"]["good_tracks"]
            res["dq_yield_flags"] = {
                "event_flag": ("incomplete_input" if res["n_files_bad"] > 0 else ("low" if ev_r < thresholds["events"]["low"] else ("high" if ev_r > thresholds["events"]["high"] else "normal"))),
                "track_flag": ("incomplete_input" if res["n_files_bad"] > 0 else ("low" if tr_r < thresholds["tracks"]["low"] else ("high" if tr_r > thresholds["tracks"]["high"] else "normal"))),
            }
            audit_results[res["key"]] = res
            print(f"  done {res['key']:12s}: files={res['n_files_valid']:2d} (bad={res['n_files_bad']}) "
                  f"sel={res['yields']['selected_events']:9,d} gt={res['yields']['GoodTrkYield']:9,d} "
                  f"trk/pb={tr_r:7.1f} [{res['dq_yield_flags']['track_flag']}] "
                  f"PS0_gfit={res['status_summary']['Preshower0']['gfit_pct']:5.1f}% "
                  f"PS0_off={res['histograms']['ps0_off_infer']['mean']:+5.1f}ns ({res['wall_s']:.1f}s)", flush=True)

    print(f"[Pass 2] Done in {time.perf_counter() - t1:.1f}s", flush=True)

    # Save compact summary (without raw bin arrays) and full histograms separately
    full_out = {
        "thresholds_2025": thresholds,
        "survey_2025": survey_results,
        "full_audit": audit_results,
    }
    (args.outdir / "dq_audit_2025_full.json").write_text(json.dumps(full_out, indent=2))

    compact_audit = {}
    for k, v in audit_results.items():
        v2 = {kk: vv for kk, vv in v.items() if kk != "histograms"}
        v2["hist_summary"] = {
            hk: {sk: sv for sk, sv in hv.items() if sk != "counts"}
            for hk, hv in v["histograms"].items()
        }
        compact_audit[k] = v2
    (args.outdir / "dq_audit_2025_summary.json").write_text(json.dumps({
        "thresholds_2025": thresholds,
        "survey_counts": {
            "total_surveyed": len(survey_results),
            "p0013_nonzero": len(p13_ok),
            "p0013_gfit_gt50pct": len(p13_gfit_hi),
            "p0014_nonzero": len(p14_ok),
            "p0014_gfit_gt50pct": len(p14_gfit_hi),
            "zero_goodtracks_runs": [r["run"] for r in zero_coll],
            "vetonu_low_runs": [(r["run"], r["tag"], round(r["vn0_ok_frac"] * 100, 2)) for r in vn_off],
        },
        "full_audit": compact_audit,
    }, indent=2))
    print(f"Saved outputs to {args.outdir}")


if __name__ == "__main__":
    main()
