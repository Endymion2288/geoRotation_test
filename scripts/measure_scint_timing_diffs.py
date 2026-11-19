#!/usr/bin/env python3
"""Measure throughgoing-muon scintillator localtime differences across 2022-2026 data and MC.

Computes peak (Gaussian core fit around histogram mode) and median of:
  - t_Veto2 - t_VetoNu
  - t_Preshower(all) - t_VetoNu
  - t_Preshower(top: 0,1) - t_VetoNu
  - t_Preshower(bot: 2,3) - t_VetoNu (2025-2026 4-PMT Preshower)
  - Separate status==0 vs status==512 (+16.0 ns) validation on 2nd-digitizer channels.

All times use the exact hit-selection and 2nd-digitizer correction rule:
  good iff (status & ~512) == 0
  t_corr = localtime + (16.0 if (status & 512) else 0.0)
"""

from __future__ import annotations

import argparse
import concurrent.futures
import glob
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
import uproot

STATUS_TIMING_SHIFTED = 512  # 1 << xAOD::WaveformStatus::TIMING_SHIFTED


def pick_runs(base_dir: str, max_runs: int) -> list[tuple[int, str]]:
    """Return evenly spaced (run_number, first_nonempty_phys_file) tuples."""
    p = Path(base_dir)
    if not p.is_dir():
        return []
    run_dirs = sorted([d for d in p.iterdir() if d.is_dir() and d.name.isdigit()])
    if not run_dirs:
        return []
    if len(run_dirs) > max_runs:
        idxs = np.linspace(0, len(run_dirs) - 1, max_runs, dtype=int)
        run_dirs = [run_dirs[i] for i in sorted(set(idxs.tolist()))]
    out: list[tuple[int, str]] = []
    for rd in run_dirs:
        files = sorted(rd.glob("*.root"))
        if not files:
            continue
        # Pick a mid-run file (often stable beams rather than ramp/tail)
        # Fallback to largest file in the run directory.
        files_by_size = sorted(files, key=lambda x: x.stat().st_size, reverse=True)
        out.append((int(rd.name), str(files_by_size[0])))
    return out


def estimate_peak(values: np.ndarray, bin_width: float = 0.2, window: float = 2.0) -> dict[str, Any]:
    """Return peak (parabolic/Gaussian-core mean around mode), median, std_core, and count."""
    v = values[np.isfinite(values)]
    n = int(v.size)
    if n < 20:
        return {"n": n, "peak": None, "median": None, "core_mean": None, "core_std": None, "iqr": None}

    med = float(np.median(v))
    q16, q84 = np.percentile(v, [16.0, 84.0])
    iqr = float(0.5 * (q84 - q16))

    # Build histogram in [med - 25, med + 25]
    lo, hi = med - 25.0, med + 25.0
    in_rng = v[(v >= lo) & (v <= hi)]
    if in_rng.size < 20:
        lo, hi = float(np.percentile(v, 1.0)), float(np.percentile(v, 99.0))
        in_rng = v[(v >= lo) & (v <= hi)]
    nbins = max(10, int(round((hi - lo) / bin_width)))
    counts, edges = np.histogram(in_rng, bins=nbins, range=(lo, hi))
    imax = int(np.argmax(counts))
    mode_center = float(0.5 * (edges[imax] + edges[imax + 1]))

    # Weighted log-Gaussian / quadratic peak estimate around mode (+/- 3 bins)
    i0 = max(0, imax - 4)
    i1 = min(len(counts), imax + 5)
    centers = 0.5 * (edges[i0:i1] + edges[i0 + 1 : i1 + 1])
    c_win = counts[i0:i1].astype(float)
    pos = c_win > max(3.0, 0.25 * counts[imax])
    peak = mode_center
    if np.count_nonzero(pos) >= 3:
        try:
            # Quadratic fit to log(counts) -> Gaussian peak
            coeffs = np.polyfit(centers[pos], np.log(c_win[pos]), 2, w=np.sqrt(c_win[pos]))
            if coeffs[0] < 0:
                cand = -coeffs[1] / (2.0 * coeffs[0])
                if abs(cand - mode_center) <= 1.5:
                    peak = float(cand)
        except Exception:
            peak = mode_center

    core = v[np.abs(v - peak) <= window]
    core_mean = float(np.mean(core)) if core.size >= 10 else peak
    core_std = float(np.std(core)) if core.size >= 10 else iqr
    return {
        "n": n,
        "n_core": int(core.size),
        "peak": round(peak, 4),
        "mode": round(mode_center, 4),
        "median": round(med, 4),
        "core_mean": round(core_mean, 4),
        "core_std": round(core_std, 4),
        "iqr": round(iqr, 4),
    }


def station_mean_time(
    arr: dict[str, np.ndarray],
    prefixes: list[str],
    time_suffix: str,
    status_filter: str = "all",
    apply_512_shift: bool = True,
) -> np.ndarray:
    """Compute per-event mean corrected localtime across given channel prefixes."""
    n = len(next(iter(arr.values())))
    tsum = np.zeros(n, dtype=np.float64)
    cnt = np.zeros(n, dtype=np.int32)
    for pref in prefixes:
        tk = f"{pref}_{time_suffix}"
        sk = f"{pref}_status"
        if tk not in arr or sk not in arr:
            continue
        t = arr[tk].astype(np.float64)
        st = arr[sk].astype(np.int64)
        if status_filter == "all":
            good = (st & ~STATUS_TIMING_SHIFTED) == 0
        elif status_filter == "only0":
            good = st == 0
        elif status_filter == "only512":
            good = st == STATUS_TIMING_SHIFTED
        else:
            raise ValueError(status_filter)
        # Guard against zero-initialized / unphysical localtimes in empty hits
        good &= np.isfinite(t) & (t > 1.0)
        if apply_512_shift:
            t_corr = np.where((st & STATUS_TIMING_SHIFTED) != 0, t + 16.0, t)
        else:
            t_corr = t
        tsum += np.where(good, t_corr, 0.0)
        cnt += good.astype(np.int32)
    out = np.full(n, np.nan, dtype=np.float64)
    ok = cnt > 0
    out[ok] = tsum[ok] / cnt[ok]
    return out


def analyze_file(task: tuple[str, int, str, int]) -> dict[str, Any]:
    """Process one PHYS ROOT file and extract scintillator time-difference peaks."""
    period, run, fpath, max_entries = task
    try:
        with uproot.open(fpath) as f:
            if "nt" not in f:
                return {"period": period, "run": run, "file": fpath, "error": "no nt tree"}
            t = f["nt"]
            n_tot = t.num_entries
            if n_tot < 50:
                return {"period": period, "run": run, "file": fpath, "n_entries": n_tot, "error": "too_few_entries"}
            keys = set(t.keys())
            time_suffix = "localtime" if "VetoNu0_localtime" in keys else "time"

            # Detect channel prefixes
            vetonu_pref = [p for p in ("VetoNu0", "VetoNu1") if f"{p}_{time_suffix}" in keys]
            if "Veto20_" + time_suffix in keys:
                veto2_pref = [p for p in ("Veto20", "Veto21") if f"{p}_{time_suffix}" in keys]
            elif "VetoSt20_" + time_suffix in keys:
                veto2_pref = [p for p in ("VetoSt20", "VetoSt21") if f"{p}_{time_suffix}" in keys]
            else:
                veto2_pref = [p for p in ("Veto0", "Veto1") if f"{p}_{time_suffix}" in keys]

            ps_top_pref = [p for p in ("Preshower0", "Preshower1") if f"{p}_{time_suffix}" in keys]
            ps_bot_pref = [p for p in ("Preshower2", "Preshower3") if f"{p}_{time_suffix}" in keys]
            ps_all_pref = ps_top_pref + ps_bot_pref

            needed = set()
            for pref in vetonu_pref + veto2_pref + ps_all_pref:
                needed.add(f"{pref}_{time_suffix}")
                needed.add(f"{pref}_status")
                rc = f"{pref}_raw_charge"
                if rc in keys:
                    needed.add(rc)
            for meta in ("longTracks", "run", "VetoNu0_reduced_charge"):
                if meta in keys:
                    needed.add(meta)

            stop = min(n_tot, max_entries) if max_entries > 0 else n_tot
            arr = t.arrays(sorted(needed), entry_stop=stop, library="np")
    except Exception as exc:
        return {"period": period, "run": run, "file": fpath, "error": repr(exc)}

    n_read = len(next(iter(arr.values())))
    # Throughgoing muon selection: 1 long track (if branch present)
    if "longTracks" in arr:
        mu_mask = arr["longTracks"] == 1
    else:
        mu_mask = np.ones(n_read, dtype=bool)
    # Require real VetoNu signal (> 40 pC on at least one layer if available)
    if "VetoNu0_raw_charge" in arr and "VetoNu1_raw_charge" in arr:
        mu_mask &= (arr["VetoNu0_raw_charge"] > 40.0) & (arr["VetoNu1_raw_charge"] > 40.0)

    # Check status bit 512 counts
    st_counts: dict[str, dict[str, int]] = {}
    for label, prefs in [("vetonu", vetonu_pref), ("veto2", veto2_pref), ("ps_top", ps_top_pref), ("ps_bot", ps_bot_pref)]:
        if not prefs:
            continue
        n0 = n512 = nother = 0
        for p in prefs:
            st = arr[f"{p}_status"][mu_mask].astype(np.int64)
            n0 += int((st == 0).sum())
            n512 += int((st == STATUS_TIMING_SHIFTED).sum())
            nother += int(((st & ~STATUS_TIMING_SHIFTED) != 0).sum())
        st_counts[label] = {"st0": n0, "st512": n512, "other": nother}

    t_vn = station_mean_time(arr, vetonu_pref, time_suffix, "all", True)
    t_v2 = station_mean_time(arr, veto2_pref, time_suffix, "all", True)
    t_v2_0 = station_mean_time(arr, veto2_pref, time_suffix, "only0", True)
    t_v2_512 = station_mean_time(arr, veto2_pref, time_suffix, "only512", True)
    t_v2_512_noshift = station_mean_time(arr, veto2_pref, time_suffix, "only512", False)

    t_ps_all = station_mean_time(arr, ps_all_pref, time_suffix, "all", True)
    t_ps_top = station_mean_time(arr, ps_top_pref, time_suffix, "all", True)
    t_ps_bot = station_mean_time(arr, ps_bot_pref, time_suffix, "all", True) if ps_bot_pref else np.full(n_read, np.nan)
    t_ps_bot_0 = station_mean_time(arr, ps_bot_pref, time_suffix, "only0", True) if ps_bot_pref else np.full(n_read, np.nan)
    t_ps_bot_512 = station_mean_time(arr, ps_bot_pref, time_suffix, "only512", True) if ps_bot_pref else np.full(n_read, np.nan)
    t_ps_bot_512_noshift = station_mean_time(arr, ps_bot_pref, time_suffix, "only512", False) if ps_bot_pref else np.full(n_read, np.nan)

    # Also check individual channel offsets for 2025-2026 Preshower (0,1,2,3)
    per_ch_ps = {}
    for p in ps_all_pref:
        tp = station_mean_time(arr, [p], time_suffix, "all", True)
        per_ch_ps[p] = estimate_peak((tp - t_vn)[mu_mask])

    actual_run = int(arr["run"][0]) if "run" in arr and len(arr["run"]) > 0 else run

    return {
        "period": period,
        "run": actual_run,
        "file": fpath,
        "time_branch": time_suffix,
        "n_entries": n_tot,
        "n_read": n_read,
        "n_mu_sel": int(mu_mask.sum()),
        "ps_channels": ps_all_pref,
        "status_counts": st_counts,
        "veto2_minus_vetonu": estimate_peak((t_v2 - t_vn)[mu_mask]),
        "veto2_st0_minus_vetonu": estimate_peak((t_v2_0 - t_vn)[mu_mask]),
        "veto2_st512_corr_minus_vetonu": estimate_peak((t_v2_512 - t_vn)[mu_mask]),
        "veto2_st512_raw_minus_vetonu": estimate_peak((t_v2_512_noshift - t_vn)[mu_mask]),
        "ps_all_minus_vetonu": estimate_peak((t_ps_all - t_vn)[mu_mask]),
        "ps_top_minus_vetonu": estimate_peak((t_ps_top - t_vn)[mu_mask]),
        "ps_bot_minus_vetonu": estimate_peak((t_ps_bot - t_vn)[mu_mask]),
        "ps_bot_st0_minus_vetonu": estimate_peak((t_ps_bot_0 - t_vn)[mu_mask]),
        "ps_bot_st512_corr_minus_vetonu": estimate_peak((t_ps_bot_512 - t_vn)[mu_mask]),
        "ps_bot_st512_raw_minus_vetonu": estimate_peak((t_ps_bot_512_noshift - t_vn)[mu_mask]),
        "per_ch_ps": per_ch_ps,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--max-entries", type=int, default=200000)
    ap.add_argument("--workers", type=int, default=min(16, os.cpu_count() or 8))
    args = ap.parse_args()

    tasks: list[tuple[str, int, str, int]] = []

    # 2022
    for run, fp in pick_runs("/eos/experiment/faser/phys/2022/r0021", 12):
        tasks.append(("2022_r0021", run, fp, args.max_entries))
    for run, fp in pick_runs("/eos/experiment/faser/phys/2022/r0022", 6):
        tasks.append(("2022_r0022", run, fp, args.max_entries))
    # 2023
    for run, fp in pick_runs("/eos/experiment/faser/phys/2023/r0021", 15):
        tasks.append(("2023_r0021", run, fp, args.max_entries))
    # 2024
    for run, fp in pick_runs("/eos/experiment/faser/phys/2024/r0023", 5):
        tasks.append(("2024_pre_r0023", run, fp, args.max_entries))
    for run, fp in pick_runs("/eos/experiment/faser/phys/2024/dev", 10):
        tag = "2024_CaloNu_dev" if 15845 <= run <= 16927 else ("2024_post_dev" if run >= 16928 else "2024_pre_dev")
        tasks.append((tag, run, fp, args.max_entries))
    # 2025
    for run, fp in pick_runs("/eos/experiment/faser/phys/2025/p0013", 18):
        tasks.append(("2025_p0013", run, fp, args.max_entries))
    for run, fp in pick_runs("/eos/experiment/faser/phys/2025/p0014", 15):
        tasks.append(("2025_p0014", run, fp, args.max_entries))
    for run, fp in pick_runs("/eos/experiment/faser/phys/2025/r0023", 5):
        tasks.append(("2025_r0023", run, fp, args.max_entries))
    # 2026
    for run, fp in pick_runs("/eos/experiment/faser/phys/2026/dev", 10):
        tasks.append(("2026_dev", run, fp, args.max_entries))

    # MC samples
    mc_dirs = [
        ("MC22_fluka_210001", 210001, "/eos/experiment/faser/data0/sim/mc22/fluka/210001/phy/r0013"),
        ("MC24_fluka_210010", 210010, "/eos/experiment/faser/data0/sim/mc24/fluka/210010/phy/s0013-r0019"),
        ("MC24_PG_100012_r23", 100012, "/eos/experiment/faser/data0/sim/mc24/particle_gun/100012/phy/s0012-r0023"),
        ("MC24_PG_100043_r23", 100043, "/eos/experiment/faser/data0/sim/mc24/particle_gun/100043/phy/s0013-r0023"),
        ("MC24_PG_100044_r23", 100044, "/eos/experiment/faser/data0/sim/mc24/particle_gun/100044/phy/s0013-r0023"),
        ("MC24_PG_100049_s15", 100049, "/eos/experiment/faser/data0/sim/mc24/particle_gun/100049/phy/s0015"),
        ("MC24_PG_100050_s15", 100050, "/eos/experiment/faser/data0/sim/mc24/particle_gun/100050/phy/s0015"),
        ("MC24_Genie_200139", 200139, "/eos/experiment/faser/data0/sim/mc24/genie/200139/phy/caloNu"),
        ("MC24_FS_110051", 110051, "/eos/experiment/faser/data0/sim/mc24/foresee/110051/phy/s0008-dev"),
    ]
    for label, dsid, d in mc_dirs:
        fs = sorted(glob.glob(os.path.join(d, "*.root")))
        if fs:
            tasks.append((label, dsid, fs[0], args.max_entries))

    print(f"Scanning {len(tasks)} PHYS files with {args.workers} workers...")
    results: list[dict[str, Any]] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as ex:
        for res in ex.map(analyze_file, tasks):
            results.append(res)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(results, fh, indent=2)

    # Print concise table grouped by period
    print(f"\n{'Period':18s} {'Run':>6s} {'N_mu':>7s} | {'V2-VN peak':>10s} {'std':>5s} | {'PSall-VN':>10s} {'PStop-VN':>10s} {'PSbot-VN':>10s} | {'V2_512corr':>10s} {'PSb_512c':>9s}")
    print("-" * 118)
    for r in sorted(results, key=lambda x: (x["period"], x["run"])):
        if "error" in r:
            print(f"{r['period']:18s} {r['run']:6d} ERROR: {r['error']}")
            continue
        v2 = r["veto2_minus_vetonu"]
        psa = r["ps_all_minus_vetonu"]
        pst = r["ps_top_minus_vetonu"]
        psb = r["ps_bot_minus_vetonu"]
        v2_512 = r["veto2_st512_corr_minus_vetonu"]
        psb_512 = r["ps_bot_st512_corr_minus_vetonu"]
        fmt = lambda x: f"{x:10.3f}" if x is not None else f"{'-':>10s}"
        fmt5 = lambda x: f"{x:5.2f}" if x is not None else f"{'-':>5s}"
        fmt9 = lambda x: f"{x:9.3f}" if x is not None else f"{'-':>9s}"
        print(
            f"{r['period']:18s} {r['run']:6d} {r['n_mu_sel']:7d} | "
            f"{fmt(v2['peak'])} {fmt5(v2['core_std'])} | "
            f"{fmt(psa['peak'])} {fmt(pst['peak'])} {fmt(psb['peak'])} | "
            f"{fmt(v2_512['peak'])} {fmt9(psb_512['peak'])}"
        )


if __name__ == "__main__":
    main()
