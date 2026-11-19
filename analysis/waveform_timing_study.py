#!/usr/bin/env python3
"""Real-data (2022-2026) and MC24 scintillator station timing-difference study.

Follows the exact GoodTracks and Data Quality selection from
FASER_DQ_Analyis/FASER_DQ_RDF.py and produces Fig. 13-style overlaid
normalized timing distributions (modeled on Electronic_2D_measurement_internal_note-6.pdf)
across all Run 3 periods to empirically explain:
  1. Veto2 - VetoNu time difference (default, CaloNu Dig2 +2.6 ns shift, -16 ns bit 512, MC)
  2. Preshower (Dig1: ch 12,13) - VetoNu time difference (2022-2024 vs 2025-2026 upgrade vs MC)
  3. Preshower (Dig2: ch 24,25) - VetoNu and Dig2 - Dig1 extra offset in 2025-2026
  4. Station-by-station localtime vs triggertime comparison and per-run stability timeline.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
import re
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

DQ_REPO = Path("/afs/cern.ch/user/x/xcheng/cernbox/FASER/FASER_DQ_Analyis")
DATASETS_JSON = DQ_REPO / "faser_filelists" / "datasets.json"
RDF_DEFINES_H = DQ_REPO / "RDFDefines.h"

DEFAULT_OUT_DIR = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/"
    "stage04_waveform_reduced_charge/timing_study_fig13"
)

# Histogram specifications: (nbins, xmin, xmax)
HIST_SPECS: Dict[str, Tuple[int, float, float]] = {
    # Main Fig. 13 style distributions (0.2 ns bin width)
    "veto2_vetonu_raw": (500, -20.0, 80.0),
    "veto2_vetonu_corr": (500, -20.0, 80.0),
    "veto2_vetonu_s0": (500, -20.0, 80.0),
    "veto2_vetonu_s512_raw": (500, -20.0, 80.0),
    "veto2_vetonu_s512_corr": (500, -20.0, 80.0),
    "ps_dig1_vetonu": (500, -20.0, 80.0),
    "ps_dig2_vetonu_raw": (500, -20.0, 80.0),
    "ps_dig2_vetonu_corr": (500, -20.0, 80.0),
    "ps_dig2_minus_dig1_raw": (400, -30.0, 50.0),
    "ps_dig2_minus_dig1_corr": (400, -30.0, 50.0),
    "ps_all_vetonu_raw": (500, -20.0, 80.0),
    "ps_all_vetonu_corr16": (500, -20.0, 80.0),
    "ps_all_vetonu_aligned": (500, -20.0, 80.0),
    "timing_vetonu": (500, -20.0, 80.0),
    "veto1_vetonu": (500, -20.0, 80.0),
    # Triggertime differences (to compare against COOL TriggerOffset)
    "veto2_vetonu_trig": (400, -40.0, 40.0),
    "ps_dig1_vetonu_trig": (400, -40.0, 40.0),
    "ps_dig2_vetonu_trig": (400, -40.0, 40.0),
    # Absolute localtimes (wide 100..900 ns covers 2022/MC ~800 ns and 2023-2026 ~150-240 ns)
    "vetonu_localtime": (800, 100.0, 900.0),
    "veto2_localtime_corr": (800, 100.0, 900.0),
    "ps_dig1_localtime": (800, 100.0, 900.0),
    "ps_dig2_localtime_corr": (800, 100.0, 900.0),
    # Fallback-relevant Preshower2/3 charge distribution
    "ps_dig2_raw_charge": (200, 0.0, 40.0),
}

CPP_HELPER_CODE = r"""
#ifndef WAVEFORM_TIMING_STUDY_HELPERS
#define WAVEFORM_TIMING_STUDY_HELPERS

#include <cmath>
#include <limits>

inline double mean2_local(float t0, unsigned int s0, float q0, double qmin0,
                          float t1, unsigned int s1, float q1, double qmin1,
                          bool apply_16ns, bool only_s0, bool only_s512) {
  int n = 0;
  double sum = 0.0;
  auto add = [&](float t, unsigned int s, float q, double qmin) {
    if (q <= qmin) return;
    if ((s & ~512U) != 0U) return;
    if (only_s0 && s != 0U) return;
    if (only_s512 && s != 512U) return;
    double tv = static_cast<double>(t);
    if (apply_16ns && (s & 512U)) tv += 16.0;
    sum += tv;
    ++n;
  };
  add(t0, s0, q0, qmin0);
  add(t1, s1, q1, qmin1);
  return (n > 0) ? (sum / n) : std::numeric_limits<double>::quiet_NaN();
}

inline double mean2_trig(float t0, unsigned int s0, float q0, double qmin0,
                         float t1, unsigned int s1, float q1, double qmin1) {
  int n = 0;
  double sum = 0.0;
  auto add = [&](float t, unsigned int s, float q, double qmin) {
    if (q <= qmin) return;
    if ((s & ~512U) != 0U) return;
    sum += static_cast<double>(t);
    ++n;
  };
  add(t0, s0, q0, qmin0);
  add(t1, s1, q1, qmin1);
  return (n > 0) ? (sum / n) : std::numeric_limits<double>::quiet_NaN();
}

inline double mean4_local(float t0, unsigned int s0, float q0, double qmin0,
                          float t1, unsigned int s1, float q1, double qmin1,
                          float t2, unsigned int s2, float q2, double qmin2,
                          float t3, unsigned int s3, float q3, double qmin3,
                          bool apply_16ns, double dig2_extra) {
  int n = 0;
  double sum = 0.0;
  auto add = [&](float t, unsigned int s, float q, double qmin, double extra) {
    if (q <= qmin) return;
    if ((s & ~512U) != 0U) return;
    double tv = static_cast<double>(t);
    if (apply_16ns && (s & 512U)) tv += 16.0;
    tv -= extra;
    sum += tv;
    ++n;
  };
  add(t0, s0, q0, qmin0, 0.0);
  add(t1, s1, q1, qmin1, 0.0);
  add(t2, s2, q2, qmin2, dig2_extra);
  add(t3, s3, q3, qmin3, dig2_extra);
  return (n > 0) ? (sum / n) : std::numeric_limits<double>::quiet_NaN();
}

inline double diff_valid(double a, double b) {
  if (std::isnan(a) || std::isnan(b)) return -9999.0;
  return a - b;
}

#endif
"""


def parse_year_filelists() -> Tuple[Dict[str, dict], Dict[str, Dict[int, List[str]]]]:
    datasets = json.loads(DATASETS_JSON.read_text())
    filelists: Dict[str, Dict[int, List[str]]] = {}
    for year in ["2022", "2023", "2024", "2025", "2026"]:
        fl_path = DQ_REPO / "faser_filelists" / f"{year}_filelist.txt"
        runs: Dict[int, List[str]] = {}
        for line in fl_path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            m = re.search(r"/0*(\d{4,6})/Faser-", line)
            if not m:
                continue
            run = int(m.group(1))
            runs.setdefault(run, []).append(line)
        filelists[year] = runs
    return datasets, filelists


def load_grl_intervals(grl_prefix: str) -> Dict[int, dict]:
    data = json.loads(Path(f"{grl_prefix}.json").read_text())
    return {int(k): v for k, v in data.items()}


def pick_evenly_spaced(runs: List[int], k: int, must_include: List[int] = None) -> List[int]:
    runs_s = sorted(runs)
    if len(runs_s) <= k:
        return runs_s
    chosen = set()
    if must_include:
        for r in must_include:
            if r in runs_s:
                chosen.add(r)
    idx = np.linspace(0, len(runs_s) - 1, k, dtype=int)
    for i in idx:
        chosen.add(runs_s[int(i)])
    return sorted(chosen)


def build_run_tasks(files_per_run: int = 2, runs_per_period: int = 12) -> List[dict]:
    datasets, filelists = parse_year_filelists()
    tasks: List[dict] = []

    def add_data_period(period: str, year: str, candidate_runs: List[int],
                        max_runs: int, must_include: List[int] = None):
        grl = load_grl_intervals(datasets[year]["grl"])
        valid_runs = [
            r for r in sorted(candidate_runs)
            if r in grl and len(grl[r].get("stable_list", [])) > 0
            and len(filelists[year].get(r, [])) > 0
        ]
        selected = pick_evenly_spaced(valid_runs, max_runs, must_include=must_include)
        for r in selected:
            all_files = sorted(filelists[year][r])
            # Pick evenly spaced files from this run (up to files_per_run)
            if len(all_files) <= files_per_run:
                chosen_files = all_files
            else:
                fidx = np.linspace(0, len(all_files) - 1, files_per_run, dtype=int)
                chosen_files = [all_files[int(i)] for i in sorted(set(fidx))]
            tasks.append({
                "period": period,
                "year": year,
                "run": r,
                "is_mc": False,
                "files": chosen_files,
                "stable_list": grl[r].get("stable_list", []),
                "excluded_list": grl[r].get("excluded_list", []),
            })

    add_data_period("2022", "2022", list(filelists["2022"].keys()), runs_per_period,
                    must_include=[7918, 8655, 9188])
    add_data_period("2023", "2023", list(filelists["2023"].keys()), runs_per_period,
                    must_include=[10555, 11066, 11703])

    runs_2024 = sorted(filelists["2024"].keys())
    add_data_period("2024_preCaloNu", "2024", [r for r in runs_2024 if r < 15821], runs_per_period,
                    must_include=[14593, 15270, 15820])
    # All available early-CaloNu runs (CaloNu installed, Veto2 still on Dig1 ch 6,7)
    add_data_period("2024_earlyCaloNu", "2024", [r for r in runs_2024 if 15821 <= r <= 15844], 10,
                    must_include=[15838, 15839, 15840, 15844])
    add_data_period("2024_duringCaloNu", "2024", [r for r in runs_2024 if 15845 <= r <= 16924], runs_per_period + 3,
                    must_include=[15849, 16000, 16517, 16923])
    add_data_period("2024_postCaloNu", "2024", [r for r in runs_2024 if r >= 16925], runs_per_period,
                    must_include=[16929, 17023, 17230])

    runs_2025 = sorted(filelists["2025"].keys())
    add_data_period("2025", "2025", runs_2025, runs_per_period + 6,
                    must_include=[19116, 19336, 19518, 19707, 19878, 20012, 20481, 20677])
    runs_2026 = sorted(filelists["2026"].keys())
    add_data_period("2026", "2026", runs_2026, runs_per_period + 3,
                    must_include=[22132, 22404, 22700, 23176, 23241])

    # Add MC24 throughgoing muon samples (FaserNu 2024 & CaloNu geometry)
    muon_mc_specs = [
        (
            "MC24_Muon_FaserNu",
            100124,
            "/eos/experiment/faser/data0/sim/mc24/particle_gun/100124/phy/s0013-r0022",
            4,
        ),
        (
            "MC24_Muon_CaloNu",
            100050,
            "/eos/experiment/faser/data0/sim/mc24/particle_gun/100050/bck/s0015",
            4,
        ),
    ]
    for period, rnum, mdir, maxf in muon_mc_specs:
        mfiles = sorted(str(p) for p in Path(mdir).glob("*PHYS.root"))[:maxf]
        if mfiles:
            tasks.append({
                "period": period,
                "year": "MC24",
                "run": rnum,
                "is_mc": True,
                "files": mfiles,
                "stable_list": [],
                "excluded_list": [],
            })

    # Add MC24 Genie default-flux samples (all 3 categories from FASER_FaserMC.pdf)
    genie_groups = [
        ("MC24_Genie_22_23", [200082, 200083, 200089], "s0014"),
        ("MC24_Genie_24_FaserNu", [200093, 200094, 200100], "s0014"),
        ("MC24_Genie_24_CaloNu", [200139, 200140, 200146], "s0015"),
    ]
    for period, rlist, tag in genie_groups:
        for r in rlist:
            gdir = Path(f"/eos/experiment/fasernu-data0/faser/sim/mc24/genie/{r}/phy/{tag}")
            gfiles = sorted(str(p) for p in gdir.glob("*PHYS.root"))[:6]
            if gfiles:
                tasks.append({
                    "period": period,
                    "year": "MC24",
                    "run": r,
                    "is_mc": True,
                    "files": gfiles,
                    "stable_list": [],
                    "excluded_list": [],
                })

    return tasks


def make_interval_cut(intervals: List[dict], run: int) -> str:
    if not intervals:
        return ""
    return " || ".join(
        f"((eventTime >= {iv['start_utime']}) && (eventTime <= {iv['stop_utime']}) && (run == {run}))"
        for iv in intervals
    )


def _worker_init():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gInterpreter.Declare("using namespace ROOT::VecOps;")
    ROOT.gInterpreter.Declare(f'#include "{RDF_DEFINES_H}"')
    ROOT.gInterpreter.Declare(CPP_HELPER_CODE)


def process_one_task(task: dict) -> dict:
    import ROOT

    period = task["period"]
    run = task["run"]
    is_mc = task["is_mc"]
    files = task["files"]

    # Validate readable files with non-empty nt tree
    valid_files = []
    for fp in files:
        tf = ROOT.TFile.Open(fp)
        if not tf or tf.IsZombie():
            continue
        t = tf.Get("nt")
        if t and t.GetEntries() > 0:
            valid_files.append(fp)
        tf.Close()

    if not valid_files:
        return {
            "period": period,
            "run": run,
            "is_mc": is_mc,
            "status": "no_valid_files",
            "files_used": 0,
            "cutflow": {},
            "hists": {},
            "status_counts": {},
        }

    df = ROOT.RDataFrame("nt", valid_files)
    cols = set(str(c) for c in df.GetColumnNames())
    # Compatibility alias for older MC/data ntuples using Veto0/1/2 instead of Veto20/21/10
    if "Veto20_localtime" not in cols and "Veto0_localtime" in cols:
        for suf in ("localtime", "triggertime", "raw_charge", "status"):
            if f"Veto0_{suf}" in cols:
                df = df.Alias(f"Veto20_{suf}", f"Veto0_{suf}")
            if f"Veto1_{suf}" in cols:
                df = df.Alias(f"Veto21_{suf}", f"Veto1_{suf}")
            if f"Veto2_{suf}" in cols:
                df = df.Alias(f"Veto10_{suf}", f"Veto2_{suf}")
        cols = set(str(c) for c in df.GetColumnNames())
    has_ps23 = ("Preshower2_localtime" in cols) and ("Preshower3_localtime" in cols)
    has_veto11 = ("Veto11_localtime" in cols)
    has_triggertime = ("VetoNu0_triggertime" in cols) and (not is_mc)

    c_total = df.Count()

    if not is_mc:
        good_cut = make_interval_cut(task["stable_list"], run)
        if good_cut:
            df = df.Filter(good_cut, "GoodTimes")
        excl_cut = make_interval_cut(task["excluded_list"], run)
        if excl_cut:
            df = df.Filter(f"!({excl_cut})", "NotExcludedTimes")
        c_grl = df.Count()
        df = df.Filter(
            "(Timing0_status & 4) == 0 && (Timing1_status & 4) == 0 && "
            "(Timing2_status & 4) == 0 && (Timing3_status & 4) == 0",
            "NoTimingSat",
        )
        df = df.Filter("distanceToCollidingBCID == 0", "Colliding")
        df = df.Filter("(TAP & 4) != 0", "TimingTrigger")
    else:
        c_grl = c_total
        df = df.Filter(
            "(Timing0_status & 4) == 0 && (Timing1_status & 4) == 0 && "
            "(Timing2_status & 4) == 0 && (Timing3_status & 4) == 0",
            "NoTimingSat",
        )
    c_dq = df.Count()

    # Exact FASER_DQ_RDF.py GoodTracks definition
    df = df.Define("Track_nHits", "Track_nDoF + 5")
    df = df.Define("Track_chi2_per_dof", "Track_Chi2 / Track_nDoF")
    df = df.Define(
        "GoodTracks",
        "Track_nLayers >= 7 && Track_chi2_per_dof < 25 && Track_nHits >= 12 "
        "&& Track_pz0 > 20000 && RemoveDuplicates(Track_p0)",
    )
    df = df.Define("NGoodTracks", "Track_nLayers[GoodTracks].size()")
    df = df.Filter("NGoodTracks >= 1", "NGoodTracks>=1")
    c_goodtrk = df.Count()

    # Status counts on GoodTracks events
    sc = {}
    sc["n_goodtrk"] = c_goodtrk
    sc["v20_s0"] = df.Filter("Veto20_status == 0").Count()
    sc["v20_s512"] = df.Filter("Veto20_status == 512").Count()
    sc["v20_ovf"] = df.Filter("(Veto20_status & 4) != 0").Count()
    sc["ps0_s0"] = df.Filter("Preshower0_status == 0").Count()
    sc["ps0_gfit"] = df.Filter("(Preshower0_status & 16) != 0").Count()
    sc["ps0_thr"] = df.Filter("(Preshower0_status & 1) != 0").Count()
    if has_ps23:
        sc["ps2_s0"] = df.Filter("Preshower2_status == 0").Count()
        sc["ps2_s512"] = df.Filter("Preshower2_status == 512").Count()
        sc["ps2_thr"] = df.Filter("(Preshower2_status & 1) != 0").Count()
        sc["ps3_s0"] = df.Filter("Preshower3_status == 0").Count()
        sc["ps3_s512"] = df.Filter("Preshower3_status == 512").Count()
        sc["ps3_thr"] = df.Filter("(Preshower3_status & 1) != 0").Count()

    # Define station mean times using FASER_DQ_RDF.py hit thresholds (VetoNu/Veto2 > 40 pC, PS0/1 > 2.5 pC)
    # and status & ~512 == 0 (matching MyWaveformTimeRangeAlg / NtupleDumperAlg)
    df = df.Define(
        "t_vn",
        "mean2_local(VetoNu0_localtime, VetoNu0_status, VetoNu0_raw_charge, 40.0, "
        "VetoNu1_localtime, VetoNu1_status, VetoNu1_raw_charge, 40.0, true, false, false)",
    )
    df = df.Define(
        "t_v2_raw",
        "mean2_local(Veto20_localtime, Veto20_status, Veto20_raw_charge, 40.0, "
        "Veto21_localtime, Veto21_status, Veto21_raw_charge, 40.0, false, false, false)",
    )
    df = df.Define(
        "t_v2_corr",
        "mean2_local(Veto20_localtime, Veto20_status, Veto20_raw_charge, 40.0, "
        "Veto21_localtime, Veto21_status, Veto21_raw_charge, 40.0, true, false, false)",
    )
    df = df.Define(
        "t_v2_s0",
        "mean2_local(Veto20_localtime, Veto20_status, Veto20_raw_charge, 40.0, "
        "Veto21_localtime, Veto21_status, Veto21_raw_charge, 40.0, false, true, false)",
    )
    df = df.Define(
        "t_v2_s512_raw",
        "mean2_local(Veto20_localtime, Veto20_status, Veto20_raw_charge, 40.0, "
        "Veto21_localtime, Veto21_status, Veto21_raw_charge, 40.0, false, false, true)",
    )
    df = df.Define(
        "t_v2_s512_corr",
        "mean2_local(Veto20_localtime, Veto20_status, Veto20_raw_charge, 40.0, "
        "Veto21_localtime, Veto21_status, Veto21_raw_charge, 40.0, true, false, true)",
    )
    df = df.Define(
        "t_ps_dig1",
        "mean2_local(Preshower0_localtime, Preshower0_status, Preshower0_raw_charge, 2.5, "
        "Preshower1_localtime, Preshower1_status, Preshower1_raw_charge, 2.5, true, false, false)",
    )
    df = df.Define(
        "t_timing",
        "mean4_local(Timing0_localtime, Timing0_status, Timing0_raw_charge, 10.0, "
        "Timing1_localtime, Timing1_status, Timing1_raw_charge, 10.0, "
        "Timing2_localtime, Timing2_status, Timing2_raw_charge, 10.0, "
        "Timing3_localtime, Timing3_status, Timing3_raw_charge, 10.0, true, 0.0)",
    )
    if has_veto11:
        df = df.Define(
            "t_v1",
            "mean2_local(Veto10_localtime, Veto10_status, Veto10_raw_charge, 40.0, "
            "Veto11_localtime, Veto11_status, Veto11_raw_charge, 40.0, true, false, false)",
        )
    else:
        df = df.Define(
            "t_v1",
            "mean2_local(Veto10_localtime, Veto10_status, Veto10_raw_charge, 40.0, "
            "0.0f, 1U, 0.0f, 40.0, true, false, false)",
        )

    df = df.Define("dt_v2_vn_raw", "diff_valid(t_v2_raw, t_vn)")
    df = df.Define("dt_v2_vn_corr", "diff_valid(t_v2_corr, t_vn)")
    df = df.Define("dt_v2_vn_s0", "diff_valid(t_v2_s0, t_vn)")
    df = df.Define("dt_v2_vn_s512_raw", "diff_valid(t_v2_s512_raw, t_vn)")
    df = df.Define("dt_v2_vn_s512_corr", "diff_valid(t_v2_s512_corr, t_vn)")
    df = df.Define("dt_ps1_vn", "diff_valid(t_ps_dig1, t_vn)")
    df = df.Define("dt_tim_vn", "diff_valid(t_timing, t_vn)")
    df = df.Define("dt_v1_vn", "diff_valid(t_v1, t_vn)")

    # Absolute localtimes (only when valid)
    df = df.Define("lt_vn", "std::isnan(t_vn) ? -9999.0 : t_vn")
    df = df.Define("lt_v2", "std::isnan(t_v2_corr) ? -9999.0 : t_v2_corr")
    df = df.Define("lt_ps1", "std::isnan(t_ps_dig1) ? -9999.0 : t_ps_dig1")

    if has_ps23:
        # Note: Preshower2/3 are lower-gain bottom PMTs; status & ~512 == 0 already requires above threshold
        df = df.Define(
            "t_ps_dig2_raw",
            "mean2_local(Preshower2_localtime, Preshower2_status, Preshower2_raw_charge, 0.0, "
            "Preshower3_localtime, Preshower3_status, Preshower3_raw_charge, 0.0, false, false, false)",
        )
        df = df.Define(
            "t_ps_dig2_corr",
            "mean2_local(Preshower2_localtime, Preshower2_status, Preshower2_raw_charge, 0.0, "
            "Preshower3_localtime, Preshower3_status, Preshower3_raw_charge, 0.0, true, false, false)",
        )
        df = df.Define("dt_ps2_vn_raw", "diff_valid(t_ps_dig2_raw, t_vn)")
        df = df.Define("dt_ps2_vn_corr", "diff_valid(t_ps_dig2_corr, t_vn)")
        df = df.Define("dt_ps2_ps1_raw", "diff_valid(t_ps_dig2_raw, t_ps_dig1)")
        df = df.Define("dt_ps2_ps1_corr", "diff_valid(t_ps_dig2_corr, t_ps_dig1)")
        df = df.Define("lt_ps2", "std::isnan(t_ps_dig2_corr) ? -9999.0 : t_ps_dig2_corr")
        # Combined 4-PMT Preshower mean (matching NtupleDumperAlg behavior: status & ~512 == 0)
        df = df.Define(
            "t_ps_all_raw",
            "mean4_local(Preshower0_localtime, Preshower0_status, Preshower0_raw_charge, 2.5, "
            "Preshower1_localtime, Preshower1_status, Preshower1_raw_charge, 2.5, "
            "Preshower2_localtime, Preshower2_status, Preshower2_raw_charge, 0.0, "
            "Preshower3_localtime, Preshower3_status, Preshower3_raw_charge, 0.0, false, 0.0)",
        )
        df = df.Define(
            "t_ps_all_corr16",
            "mean4_local(Preshower0_localtime, Preshower0_status, Preshower0_raw_charge, 2.5, "
            "Preshower1_localtime, Preshower1_status, Preshower1_raw_charge, 2.5, "
            "Preshower2_localtime, Preshower2_status, Preshower2_raw_charge, 0.0, "
            "Preshower3_localtime, Preshower3_status, Preshower3_raw_charge, 0.0, true, 0.0)",
        )
        df = df.Define(
            "t_ps_all_aligned",
            "mean4_local(Preshower0_localtime, Preshower0_status, Preshower0_raw_charge, 2.5, "
            "Preshower1_localtime, Preshower1_status, Preshower1_raw_charge, 2.5, "
            "Preshower2_localtime, Preshower2_status, Preshower2_raw_charge, 0.0, "
            "Preshower3_localtime, Preshower3_status, Preshower3_raw_charge, 0.0, true, 4.0)",
        )
        df = df.Define("dt_psall_vn_raw", "diff_valid(t_ps_all_raw, t_vn)")
        df = df.Define("dt_psall_vn_corr16", "diff_valid(t_ps_all_corr16, t_vn)")
        df = df.Define("dt_psall_vn_aligned", "diff_valid(t_ps_all_aligned, t_vn)")
        df = df.Define(
            "ps2_valid_q",
            "(Preshower2_status == 0U || Preshower2_status == 512U) ? Preshower2_raw_charge : -9999.0f",
        )

    if has_triggertime:
        df = df.Define(
            "tt_vn",
            "mean2_trig(VetoNu0_triggertime, VetoNu0_status, VetoNu0_raw_charge, 40.0, "
            "VetoNu1_triggertime, VetoNu1_status, VetoNu1_raw_charge, 40.0)",
        )
        df = df.Define(
            "tt_v2",
            "mean2_trig(Veto20_triggertime, Veto20_status, Veto20_raw_charge, 40.0, "
            "Veto21_triggertime, Veto21_status, Veto21_raw_charge, 40.0)",
        )
        df = df.Define(
            "tt_ps1",
            "mean2_trig(Preshower0_triggertime, Preshower0_status, Preshower0_raw_charge, 2.5, "
            "Preshower1_triggertime, Preshower1_status, Preshower1_raw_charge, 2.5)",
        )
        df = df.Define("dtt_v2_vn", "diff_valid(tt_v2, tt_vn)")
        df = df.Define("dtt_ps1_vn", "diff_valid(tt_ps1, tt_vn)")
        if has_ps23:
            df = df.Define(
                "tt_ps2",
                "mean2_trig(Preshower2_triggertime, Preshower2_status, Preshower2_raw_charge, 0.0, "
                "Preshower3_triggertime, Preshower3_status, Preshower3_raw_charge, 0.0)",
            )
            df = df.Define("dtt_ps2_vn", "diff_valid(tt_ps2, tt_vn)")

    booked = {}

    def book(key: str, col: str):
        nb, lo, hi = HIST_SPECS[key]
        booked[key] = df.Filter(f"{col} > -9000.0").Histo1D((f"h_{key}_{run}", key, nb, lo, hi), col)

    book("veto2_vetonu_raw", "dt_v2_vn_raw")
    book("veto2_vetonu_corr", "dt_v2_vn_corr")
    book("veto2_vetonu_s0", "dt_v2_vn_s0")
    book("veto2_vetonu_s512_raw", "dt_v2_vn_s512_raw")
    book("veto2_vetonu_s512_corr", "dt_v2_vn_s512_corr")
    book("ps_dig1_vetonu", "dt_ps1_vn")
    book("timing_vetonu", "dt_tim_vn")
    book("veto1_vetonu", "dt_v1_vn")
    book("vetonu_localtime", "lt_vn")
    book("veto2_localtime_corr", "lt_v2")
    book("ps_dig1_localtime", "lt_ps1")

    if has_ps23:
        book("ps_dig2_vetonu_raw", "dt_ps2_vn_raw")
        book("ps_dig2_vetonu_corr", "dt_ps2_vn_corr")
        book("ps_dig2_minus_dig1_raw", "dt_ps2_ps1_raw")
        book("ps_dig2_minus_dig1_corr", "dt_ps2_ps1_corr")
        book("ps_all_vetonu_raw", "dt_psall_vn_raw")
        book("ps_all_vetonu_corr16", "dt_psall_vn_corr16")
        book("ps_all_vetonu_aligned", "dt_psall_vn_aligned")
        book("ps_dig2_localtime_corr", "lt_ps2")
        book("ps_dig2_raw_charge", "ps2_valid_q")

    if has_triggertime:
        book("veto2_vetonu_trig", "dtt_v2_vn")
        book("ps_dig1_vetonu_trig", "dtt_ps1_vn")
        if has_ps23:
            book("ps_dig2_vetonu_trig", "dtt_ps2_vn")

    # Trigger single event loop
    out_hists = {}
    for key, hptr in booked.items():
        h = hptr.GetValue()
        nb = h.GetNbinsX()
        counts = [float(h.GetBinContent(i)) for i in range(1, nb + 1)]
        out_hists[key] = {
            "entries": int(h.GetEntries()),
            "underflow": float(h.GetBinContent(0)),
            "overflow": float(h.GetBinContent(nb + 1)),
            "counts": counts,
        }

    return {
        "period": period,
        "run": run,
        "is_mc": is_mc,
        "has_ps23": has_ps23,
        "status": "ok",
        "files_used": len(valid_files),
        "cutflow": {
            "total": int(c_total.GetValue()),
            "grl": int(c_grl.GetValue()),
            "dq": int(c_dq.GetValue()),
            "goodtracks": int(c_goodtrk.GetValue()),
        },
        "status_counts": {k: int(v.GetValue()) for k, v in sc.items()},
        "hists": out_hists,
    }


def summarize_hist(counts: np.ndarray, xmin: float, xmax: float) -> dict:
    total = float(counts.sum())
    if total <= 0:
        return {"entries": 0}
    nb = len(counts)
    edges = np.linspace(xmin, xmax, nb + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bw = (xmax - xmin) / nb

    imax = int(np.argmax(counts))
    mode = float(centers[imax])
    # Parabolic interpolation around mode if interior
    if 0 < imax < nb - 1:
        y0, y1, y2 = counts[imax - 1], counts[imax], counts[imax + 1]
        denom = (y0 - 2.0 * y1 + y2)
        if denom < 0:
            mode = float(centers[imax] + 0.5 * bw * (y0 - y2) / denom)

    # Core mean/std within +/- 4 ns of mode
    core_mask = np.abs(centers - mode) <= 4.0
    w_core = counts[core_mask]
    c_core = centers[core_mask]
    if w_core.sum() > 0:
        core_mean = float(np.sum(w_core * c_core) / w_core.sum())
        core_std = float(np.sqrt(np.sum(w_core * (c_core - core_mean) ** 2) / w_core.sum()))
        core_frac = float(w_core.sum() / total)
    else:
        core_mean = mode
        core_std = float("nan")
        core_frac = 0.0

    cdf = np.cumsum(counts) / total
    q16 = float(edges[min(nb, int(np.searchsorted(cdf, 0.16)))])
    q50 = float(edges[min(nb, int(np.searchsorted(cdf, 0.50)))])
    q84 = float(edges[min(nb, int(np.searchsorted(cdf, 0.84)))])
    mean = float(np.sum(counts * centers) / total)
    std = float(np.sqrt(np.sum(counts * (centers - mean) ** 2) / total))

    return {
        "entries": int(round(total)),
        "mode_ns": round(mode, 3),
        "core_mean_ns": round(core_mean, 3),
        "core_std_ns": round(core_std, 3),
        "core_frac_pm4ns": round(core_frac, 4),
        "median_ns": round(q50, 3),
        "q16_ns": round(q16, 3),
        "q84_ns": round(q84, 3),
        "mean_ns": round(mean, 3),
        "std_ns": round(std, 3),
    }


def main():
    parser = argparse.ArgumentParser(description="Real-data 2022-2026 & MC24 station timing-difference study")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--files-per-run", type=int, default=2)
    parser.add_argument("--runs-per-period", type=int, default=12)
    args = parser.parse_args()

    t0 = time.time()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    tasks = build_run_tasks(files_per_run=args.files_per_run, runs_per_period=args.runs_per_period)
    print(f"Prepared {len(tasks)} run tasks across periods; launching {args.workers} workers...")

    results = []
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.workers, initializer=_worker_init
    ) as pool:
        fut_to_task = {pool.submit(process_one_task, t): t for t in tasks}
        done = 0
        for fut in concurrent.futures.as_completed(fut_to_task):
            res = fut.result()
            results.append(res)
            done += 1
            if done % 15 == 0 or done == len(tasks):
                print(f"  [{done}/{len(tasks)}] finished {res['period']} run {res['run']} "
                      f"(goodtracks={res['cutflow'].get('goodtracks', 0)}, elapsed={time.time()-t0:.1f}s)",
                      flush=True)

    results.sort(key=lambda r: (r["period"], r["run"]))

    # Aggregate histograms per period & compute per-run summaries
    periods_order = [
        "2022",
        "2023",
        "2024_preCaloNu",
        "2024_earlyCaloNu",
        "2024_duringCaloNu",
        "2024_postCaloNu",
        "2025",
        "2026",
        "MC24_Muon_FaserNu",
        "MC24_Muon_CaloNu",
        "MC24_Genie_22_23",
        "MC24_Genie_24_FaserNu",
        "MC24_Genie_24_CaloNu",
    ]

    agg_hists: Dict[str, Dict[str, np.ndarray]] = {}
    agg_meta: Dict[str, dict] = {}
    per_run_rows: List[dict] = []

    for r in results:
        p = r["period"]
        if r["status"] != "ok":
            continue
        meta = agg_meta.setdefault(p, {
            "runs": [],
            "files_used": 0,
            "cutflow": {"total": 0, "grl": 0, "dq": 0, "goodtracks": 0},
            "status_counts": {},
        })
        meta["runs"].append(r["run"])
        meta["files_used"] += r["files_used"]
        for ck, cv in r["cutflow"].items():
            meta["cutflow"][ck] += cv
        for sk, sv in r["status_counts"].items():
            meta["status_counts"][sk] = meta["status_counts"].get(sk, 0) + sv

        ph = agg_hists.setdefault(p, {})
        run_row = {
            "period": p,
            "run": r["run"],
            "is_mc": r["is_mc"],
            "files_used": r["files_used"],
            "goodtracks": r["cutflow"]["goodtracks"],
            "status_counts": r["status_counts"],
        }
        for hk, hdata in r["hists"].items():
            arr = np.asarray(hdata["counts"], dtype=np.float64)
            if hk not in ph:
                ph[hk] = arr.copy()
            else:
                ph[hk] += arr
            if hk in (
                "veto2_vetonu_raw", "veto2_vetonu_corr",
                "ps_dig1_vetonu", "ps_dig2_vetonu_corr",
                "ps_dig2_minus_dig1_corr", "timing_vetonu",
                "vetonu_localtime", "veto2_localtime_corr",
                "ps_dig1_localtime", "ps_dig2_localtime_corr",
            ):
                nb, lo, hi = HIST_SPECS[hk]
                run_row[hk] = summarize_hist(arr, lo, hi)
        per_run_rows.append(run_row)

    summary = {
        "created_elapsed_s": round(time.time() - t0, 2),
        "selection": {
            "source": str(DQ_REPO / "FASER_DQ_RDF.py"),
            "dq_cuts": [
                "Annual GRL GoodTimes && !ExcludedTimes",
                "(Timing0_status & 4)==0 && (Timing1_status & 4)==0 && (Timing2_status & 4)==0 && (Timing3_status & 4)==0",
                "distanceToCollidingBCID == 0",
                "(TAP & 4) != 0",
            ],
            "good_tracks": "Track_nLayers >= 7 && (Track_Chi2/Track_nDoF) < 25 && (Track_nDoF+5) >= 12 && Track_pz0 > 20000 && RemoveDuplicates(Track_p0); NGoodTracks >= 1",
            "scintillator_cuts": {
                "VetoNu": "VetoNu{0,1}_raw_charge > 40 pC && (status & ~512) == 0",
                "Veto2": "Veto2{0,1}_raw_charge > 40 pC && (status & ~512) == 0",
                "Preshower_Dig1": "Preshower{0,1}_raw_charge > 2.5 pC && (status & ~512) == 0",
                "Preshower_Dig2": "(status & ~512) == 0 (low-gain bottom PMTs, threshold bit 1 already 0)",
                "Timing": "Timing{0..3}_raw_charge > 10 pC && (status & ~512) == 0",
            },
        },
        "hist_specs": {k: {"nbins": v[0], "xmin": v[1], "xmax": v[2]} for k, v in HIST_SPECS.items()},
        "periods": {},
        "per_run": per_run_rows,
    }

    hist_dump = {}
    for p in periods_order:
        if p not in agg_meta:
            continue
        p_stats = {}
        hist_dump[p] = {}
        for hk, arr in agg_hists[p].items():
            nb, lo, hi = HIST_SPECS[hk]
            p_stats[hk] = summarize_hist(arr, lo, hi)
            hist_dump[p][hk] = arr.tolist()
        summary["periods"][p] = {
            "n_runs": len(agg_meta[p]["runs"]),
            "runs": agg_meta[p]["runs"],
            "files_used": agg_meta[p]["files_used"],
            "cutflow": agg_meta[p]["cutflow"],
            "status_counts": agg_meta[p]["status_counts"],
            "stats": p_stats,
        }

    (args.out_dir / "timing_summary.json").write_text(json.dumps(summary, indent=2))
    (args.out_dir / "timing_histograms.json").write_text(json.dumps({
        "hist_specs": summary["hist_specs"],
        "histograms": hist_dump,
    }))

    print("\n=== SUMMARY TABLE (mode / core_mean +/- core_std [ns]) ===")
    for p in periods_order:
        if p not in summary["periods"]:
            continue
        sp = summary["periods"][p]
        st = sp["stats"]
        def fmt(k):
            d = st.get(k, {})
            if d.get("entries", 0) == 0:
                return "N/A"
            return f"{d['mode_ns']:6.2f} ({d['core_mean_ns']:6.2f}+/-{d['core_std_ns']:.2f}, N={d['entries']})"
        print(
            f"{p:22s} runs={sp['n_runs']:2d} goodTrk={sp['cutflow']['goodtracks']:9d} | "
            f"V2-VN(corr)={fmt('veto2_vetonu_corr')} | "
            f"PS1-VN={fmt('ps_dig1_vetonu')} | "
            f"PS2-VN(corr)={fmt('ps_dig2_vetonu_corr')} | "
            f"PS2-PS1(corr)={fmt('ps_dig2_minus_dig1_corr')}"
        )
    print(f"\nSaved results to {args.out_dir} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
