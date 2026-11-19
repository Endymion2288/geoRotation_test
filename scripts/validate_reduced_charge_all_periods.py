#!/usr/bin/env python3
"""End-to-end multi-core validation of fixed NtupleDumper reduced VetoNu charge
across 2022-2026 data and MC, comparing:
  1. Fixed NtupleDumper (default auto timing) vs Baseline (unmodified latest)
  2. Event-by-event comparison against Ben Wilson's MyWaveformTimeRange reference files
  3. Smoking-gun checks:
     - 2024 CaloNu Veto2 ch 22,23 status==512 (+16 ns) recovery
     - 2025-2026 4-PMT Preshower fallback (Veto2 overflow / status!=0) recovery
     - Elimination of VetoNu secondary-hit out-of-bounds memory corruption (clockPhase / nClusters0)
     - MC24 throughgoing muon (100044) and Genie neutrino backsplash (200139) behavior
"""

from __future__ import annotations

import argparse
import concurrent.futures
import glob
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np
import uproot


def find_ben_part2(ben_dir: str) -> str | None:
    p = Path(ben_dir)
    if not p.is_dir():
        return None
    hits = sorted(p.glob("*.part2.*.root"))
    if hits:
        return str(hits[0])
    hits = sorted(p.glob("*.root"))
    return str(hits[min(2, len(hits) - 1)]) if hits else None


def find_xaod(run_dir: str, seg_hint: str = "00010") -> str:
    p = Path(run_dir)
    hits = sorted(p.glob(f"*-{seg_hint}-*xAOD.root"))
    if hits:
        return str(hits[0])
    all_fs = sorted(p.glob("*xAOD.root"))
    assert all_fs, f"No xAOD files in {run_dir}"
    return str(all_fs[min(10, len(all_fs) // 2)])


SAMPLES: list[dict[str, Any]] = [
    {
        "name": "2022_008090",
        "period": "2022",
        "xaod": "/eos/experiment/faser/rec/2022/r0021/008090/Faser-Physics-008090-00010-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/waveforms-2022/008090",
        "ben_diffs": (33.165, 53.625),
        "is_mc": False,
        "nevents": 5000,
    },
    {
        "name": "2023_011703",
        "period": "2023",
        "xaod": "/eos/experiment/faser/rec/2023/r0021/011703/Faser-Physics-011703-00010-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/Waveforms-2023/011703",
        "ben_diffs": (33.165, 53.625),
        "is_mc": False,
        "nevents": 5000,
    },
    {
        "name": "2024_pre_015387",
        "period": "2024_preCaloNu",
        "xaod": "/eos/experiment/faser/rec/2024/r0022/015387/Faser-Physics-015387-00010-r0022-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/waveforms-5FilesPerChunk/015387",
        "ben_diffs": (33.165, 53.625),
        "is_mc": False,
        "nevents": 5000,
    },
    {
        "name": "2024_cn_early_015838",
        "period": "2024_earlyCaloNu_Dig1Veto",
        "xaod": "/eos/experiment/faser/rec/2024/r0022/015838/Faser-Physics-015838-00010-r0022-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/WaveformsCaloNuPeriod/015838",
        "ben_diffs": (35.783, 53.512),
        "is_mc": False,
        "nevents": 5000,
    },
    {
        "name": "2024_CaloNu_015849",
        "period": "2024_CaloNu_Dig2Veto",
        "xaod": "/eos/experiment/faser/rec/2024/r0022/015849/Faser-Physics-015849-00010-r0022-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/WaveformsCaloNuPeriod/015849",
        "ben_diffs": (35.783, 53.512),
        "is_mc": False,
        "nevents": 8000,
    },
    {
        "name": "2024_post_016929",
        "period": "2024_postCaloNu",
        "xaod": "/eos/experiment/faser/rec/2024/r0022/016929/Faser-Physics-016929-00010-r0022-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/waveforms_postCaloNu/016929",
        "ben_diffs": (33.165, 53.625),
        "is_mc": False,
        "nevents": 5000,
    },
    {
        "name": "2025_20012",
        "period": "2025_PS01_gfit_fail",
        "xaod": "/eos/experiment/faser/rec/2025/r0023/020012/Faser-Physics-020012-00010-r0023-xAOD.root",
        "ben_dir": None,
        "ben_diffs": None,
        "is_mc": False,
        "nevents": 8000,
    },
    {
        "name": "2025_20481",
        "period": "2025_4PMT_PS",
        "xaod": "/eos/experiment/faser/rec/2025/p0014/020481/Faser-Physics-020481-00276-p0014-xAOD.root",
        "ben_dir": None,
        "ben_diffs": None,
        "is_mc": False,
        "nevents": 8000,
    },
    {
        "name": "2026_22404",
        "period": "2026_4PMT_PS",
        "xaod": "/eos/experiment/faser/rec/2026/p0015/022404/Faser-Physics-022404-00010-p0015-xAOD.root",
        "ben_dir": None,
        "ben_diffs": None,
        "is_mc": False,
        "nevents": 8000,
    },
    {
        "name": "2026_23176",
        "period": "2026_late_pp",
        "xaod": find_xaod("/eos/experiment/faser/rec/2026/p0015/023176", "00090"),
        "ben_dir": None,
        "ben_diffs": None,
        "is_mc": False,
        "nevents": 8000,
    },
    {
        "name": "MC24_PG_100044",
        "period": "MC24_muon",
        "xaod": "/eos/experiment/faser/data0/sim/mc24/particle_gun/100044/rec/s0013-r0023/FaserMC-MC24_PG_mupl_fasernu_5mrad_flukaE-100044-00100-00199-s0013-r0023-xAOD.root",
        "ben_dir": None,
        "ben_diffs": None,
        "is_mc": True,
        "nevents": 3000,
    },
    {
        "name": "MC24_Genie_200139",
        "period": "MC24_genie_nu",
        "xaod": "/eos/experiment/faser/data0/sim/mc24/genie/200139/rec/caloNu/FaserMC-MC24_Genie_24_xing_angle_caloNu_epos_pi-200139-00100-00109-caloNu-xAOD.root",
        "ben_dir": None,
        "ben_diffs": None,
        "is_mc": True,
        "nevents": 1000,
    },
]


def run_ntuple_job(job: dict[str, Any]) -> dict[str, Any]:
    """Run one faser_ntuple_maker.py job in an isolated workdir."""
    work_dir = Path(job["work_dir"])
    work_dir.mkdir(parents=True, exist_ok=True)
    out_root = work_dir / f"{job['tag']}.root"
    log_path = work_dir / f"{job['tag']}.log"

    env_script = job["env_script"]
    overlay = job.get("overlay")
    latest_run = job["latest_run"]
    maker = (
        f"{overlay}/bin/faser_ntuple_maker.py"
        if overlay
        else f"{latest_run}/bin/faser_ntuple_maker.py"
    )

    extra_flags = []
    if job["is_mc"]:
        extra_flags.append("--isMC")
    else:
        extra_flags += ["--unblind", "--NoTrackFilt", "--no_stable"]
    if job.get("override_diffs") is not None:
        v_diff, ps_diff = job["override_diffs"]
        extra_flags += [
            f"--vetoVetoNuTimeDiff={v_diff}",
            f"--preshowerVetoNuTimeDiff={ps_diff}",
        ]
    cmd_parts = [
        "python3",
        maker,
        job["xaod"],
        "--outfile",
        str(out_root),
        "-n",
        str(job["nevents"]),
        *extra_flags,
    ]
    overlay_exports = ""
    if overlay:
        overlay_exports = (
            f'export LD_LIBRARY_PATH="{overlay}/lib:$LD_LIBRARY_PATH"\n'
            f'export PYTHONPATH="{overlay}/python:$PYTHONPATH"\n'
            f'export GAUDI_PLUGIN_PATH="{overlay}/lib:${{GAUDI_PLUGIN_PATH:-}}"\n'
            f'export JOBOPTSEARCHPATH="{overlay}/python:${{JOBOPTSEARCHPATH:-}}"\n'
        )

    bash_cmd = (
        f"set -e\n"
        f"source {env_script} 2>/dev/null\n"
        f"{overlay_exports}"
        f"cd {work_dir}\n"
        + " ".join(cmd_parts)
    )
    proc = subprocess.run(["bash", "-c", bash_cmd], capture_output=True, text=True)
    log_path.write_text(proc.stdout + "\n=== STDERR ===\n" + proc.stderr)
    return {
        "name": job["name"],
        "mode": job["mode"],
        "rc": proc.returncode,
        "out_root": str(out_root),
        "log": str(log_path),
        "timing_line": next(
            (ln.strip() for ln in proc.stdout.splitlines() if "Reduced VetoNu charge timing config" in ln),
            "",
        ),
    }


def compare_with_ben(nt_path: str, ben_path: str) -> dict[str, Any]:
    with uproot.open(nt_path) as f:
        t = f["nt"]
        cols = [
            c
            for c in [
                "run",
                "eventID",
                "VetoNu0_reduced_charge",
                "VetoNu1_reduced_charge",
                "VetoNu0_reduced_charge_valid",
                "VetoNu1_reduced_charge_valid",
                "VetoNu_reduced_charge_fallback",
                "VetoNu0_raw_charge",
                "VetoNu1_raw_charge",
                "Veto20_status",
                "Veto21_status",
            ]
            if c in t
        ]
        a = t.arrays(cols, library="np")
    with uproot.open(ben_path) as f:
        bt = f["tree"]
        bcols = [c for c in ["run", "event", "fallback_mode", "myVetoNu0_rawCharge35", "myVetoNu1_rawCharge35", "myVetoNu0_rawCharge35_valid"] if c in bt]
        b = bt.arrays(bcols, library="np")

    # Map event -> index in Ben's tree
    ev_nt = a["eventID"].astype(np.int64)
    ev_ben = b["event"].astype(np.int64)
    idx_map = {int(e): i for i, e in enumerate(ev_ben)}
    matched_nt = []
    matched_ben = []
    for i, e in enumerate(ev_nt):
        j = idx_map.get(int(e))
        if j is not None:
            matched_nt.append(i)
            matched_ben.append(j)
    if not matched_nt:
        return {"n_matched": 0, "ben_file": ben_path}
    mi = np.array(matched_nt, dtype=np.int64)
    mj = np.array(matched_ben, dtype=np.int64)

    q0_nt = a["VetoNu0_reduced_charge"][mi].astype(np.float64)
    q1_nt = a["VetoNu1_reduced_charge"][mi].astype(np.float64)
    q0_ben = b["myVetoNu0_rawCharge35"][mj].astype(np.float64)
    q1_ben = b["myVetoNu1_rawCharge35"][mj].astype(np.float64)

    fb_nt = a["VetoNu_reduced_charge_fallback"][mi] if "VetoNu_reduced_charge_fallback" in a else np.ones_like(mi)
    has_fb_ben = "fallback_mode" in b
    if has_fb_ben:
        fb_ben = b["fallback_mode"][mj]
        fb_agree = float(np.mean(fb_nt == fb_ben))
        active = (fb_nt > 0) & (fb_ben > 0) & np.isfinite(q0_ben) & np.isfinite(q1_ben)
    else:
        fb_agree = None
        active = (fb_nt > 0) & np.isfinite(q0_ben) & np.isfinite(q1_ben)

    dq0 = np.abs(q0_nt[active] - q0_ben[active])
    dq1 = np.abs(q1_nt[active] - q1_ben[active])
    res: dict[str, Any] = {
        "ben_file": ben_path,
        "n_matched": int(len(mi)),
        "n_active_timing": int(active.sum()),
        "fallback_mode_agree": round(fb_agree, 6) if fb_agree is not None else None,
        "dq0_med": round(float(np.median(dq0)), 6) if dq0.size else None,
        "dq0_p95": round(float(np.percentile(dq0, 95)), 6) if dq0.size else None,
        "dq0_p99": round(float(np.percentile(dq0, 99)), 6) if dq0.size else None,
        "dq0_max": round(float(np.max(dq0)), 6) if dq0.size else None,
        "dq0_exact_1e3_frac": round(float(np.mean(dq0 < 1e-3)), 6) if dq0.size else None,
        "dq1_exact_1e3_frac": round(float(np.mean(dq1 < 1e-3)), 6) if dq1.size else None,
        "dq0_within_1pC_frac": round(float(np.mean(dq0 < 1.0)), 6) if dq0.size else None,
    }
    if "myVetoNu0_rawCharge35_valid" in b and "VetoNu0_reduced_charge_valid" in a:
        v0_nt = a["VetoNu0_reduced_charge_valid"][mi][active]
        v0_ben = b["myVetoNu0_rawCharge35_valid"][mj][active]
        dv0 = np.abs(v0_nt - v0_ben)
        res["dq0_valid_exact_1e3_frac"] = round(float(np.mean(dv0 < 1e-3)), 6) if dv0.size else None
    return res


def compare_fix_vs_base(fix_path: str, base_path: str) -> dict[str, Any]:
    with uproot.open(fix_path) as ff, uproot.open(base_path) as fb:
        tf, tb = ff["nt"], fb["nt"]
        common = [
            c
            for c in [
                "eventID",
                "longTracks",
                "clockPhase",
                "nClusters0",
                "VetoNu0_raw_charge",
                "VetoNu1_raw_charge",
                "VetoNu0_reduced_charge",
                "VetoNu1_reduced_charge",
                "VetoNu0_status",
                "VetoNu1_status",
                "Veto20_status",
                "Veto21_status",
                "Veto20_localtime",
                "Veto21_localtime",
                "Preshower0_status",
                "Preshower1_status",
                "Preshower2_status",
                "Preshower3_status",
            ]
            if c in tf and c in tb
        ]
        extra_f = [c for c in ["VetoNu_reduced_charge_fallback", "VetoNu_expected_localtime", "VetoNu0_reduced_charge_valid"] if c in tf]
        af = tf.arrays(common + extra_f, library="np")
        ab = tb.arrays(common, library="np")

    assert len(af["eventID"]) == len(ab["eventID"]) and np.array_equal(af["eventID"], ab["eventID"]), "Event mismatch fix vs base"
    n = len(af["eventID"])
    # Check clockPhase corruption from out-of-bounds m_vetonu_reduced_charge[2] write in baseline!
    cp_diff = np.abs(af["clockPhase"].astype(np.float64) - ab["clockPhase"].astype(np.float64))
    n_cp_corrupt = int(np.sum(cp_diff > 1e-4))
    nc0_diff = af["nClusters0"] != ab["nClusters0"]
    n_nc0_corrupt = int(np.sum(nc0_diff))

    fb_vals, fb_cnts = np.unique(af["VetoNu_reduced_charge_fallback"], return_counts=True)
    fb_hist = {int(k): int(v) for k, v in zip(fb_vals, fb_cnts)}

    # Throughgoing muons: VetoNu0_raw_charge > 40 & VetoNu1_raw_charge > 40
    mu = (af["VetoNu0_raw_charge"] > 40.0) & (af["VetoNu1_raw_charge"] > 40.0)
    # Subset: Veto2 has status == 512 (CaloNu Dig2 shifted!)
    v2_512 = np.zeros(n, dtype=bool)
    if "Veto20_status" in af and "Veto21_status" in af:
        v2_512 = ((af["Veto20_status"] == 512) | (af["Veto21_status"] == 512)) & ((af["Veto20_status"] & ~512) == 0) & ((af["Veto21_status"] & ~512) == 0)
    # Subset: Preshower fallback in fixed (fallback == 2, e.g., Veto2 overflow/invalid!)
    ps_fb = af["VetoNu_reduced_charge_fallback"] == 2

    def summarize_q(mask: np.ndarray) -> dict[str, Any]:
        m = mask & mu
        cnt = int(m.sum())
        if cnt == 0:
            return {"n": 0}
        qf = af["VetoNu0_reduced_charge"][m]
        qb = ab["VetoNu0_reduced_charge"][m]
        qr = af["VetoNu0_raw_charge"][m]
        return {
            "n": cnt,
            "raw_med": round(float(np.median(qr)), 2),
            "fix_med": round(float(np.median(qf)), 2),
            "base_med": round(float(np.median(qb)), 2),
            "fix_lt40_frac": round(float(np.mean(qf < 40.0)), 4),
            "base_lt40_frac": round(float(np.mean(qb < 40.0)), 4),
            "fix_over_raw_med": round(float(np.median(qf / np.maximum(qr, 1e-3))), 4),
            "base_over_raw_med": round(float(np.median(qb / np.maximum(qr, 1e-3))), 4),
        }

    return {
        "n_events": n,
        "fallback_counts": fb_hist,
        "n_baseline_clockPhase_corrupted": n_cp_corrupt,
        "n_baseline_nClusters0_corrupted": n_nc0_corrupt,
        "baseline_max_clockPhase_corrupt_val": round(float(np.max(np.abs(ab["clockPhase"][cp_diff > 1e-4]))), 3) if n_cp_corrupt > 0 else 0.0,
        "mu_all": summarize_q(np.ones(n, dtype=bool)),
        "mu_v2_status512": summarize_q(v2_512),
        "mu_preshower_fallback": summarize_q(ps_fb),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlay", required=True)
    ap.add_argument("--latest-run", default="/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/run")
    ap.add_argument("--env-script", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", nargs="*", default=None, help="Optional subset of sample names to run and merge into existing summary")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    selected = [s for s in SAMPLES if (not args.only) or (s["name"] in args.only)]

    jobs: list[dict[str, Any]] = []
    for s in selected:
        # 1. Fixed auto mode
        jobs.append({
            "name": s["name"],
            "mode": "fix_auto",
            "tag": f"{s['name']}_fix_auto",
            "xaod": s["xaod"],
            "is_mc": s["is_mc"],
            "nevents": s["nevents"],
            "overlay": args.overlay,
            "latest_run": args.latest_run,
            "env_script": args.env_script,
            "work_dir": str(outdir / "ntuples" / s["name"]),
            "override_diffs": None,
        })
        # 2. Baseline (unmodified latest)
        jobs.append({
            "name": s["name"],
            "mode": "baseline",
            "tag": f"{s['name']}_baseline",
            "xaod": s["xaod"],
            "is_mc": s["is_mc"],
            "nevents": s["nevents"],
            "overlay": None,
            "latest_run": args.latest_run,
            "env_script": args.env_script,
            "work_dir": str(outdir / "ntuples" / s["name"]),
            "override_diffs": None,
        })
        # 3. Fixed exact-Ben-override parity check (only when Ben diffs differ from auto)
        if s["ben_dir"] is not None:
            jobs.append({
                "name": s["name"],
                "mode": "fix_ben_parity",
                "tag": f"{s['name']}_fix_ben_parity",
                "xaod": s["xaod"],
                "is_mc": s["is_mc"],
                "nevents": s["nevents"],
                "overlay": args.overlay,
                "latest_run": args.latest_run,
                "env_script": args.env_script,
                "work_dir": str(outdir / "ntuples" / s["name"]),
                "override_diffs": s["ben_diffs"],
            })

    print(f"Launching {len(jobs)} NtupleDumper validation jobs with {args.workers} workers...", flush=True)
    job_results: dict[tuple[str, str], dict[str, Any]] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for r in ex.map(run_ntuple_job, jobs):
            job_results[(r["name"], r["mode"])] = r
            status = "OK" if r["rc"] == 0 else f"FAIL(rc={r['rc']})"
            print(f"  [{status}] {r['name']:22s} {r['mode']:15s} {r['timing_line']}", flush=True)

    out_json = outdir / "validation_summary.json"
    summary: dict[str, Any] = json.loads(out_json.read_text()) if (args.only and out_json.is_file()) else {}
    if args.only and "2025_20485" in summary and "2025_20481" in args.only:
        summary.pop("2025_20485", None)
    for s in selected:
        name = s["name"]
        rec: dict[str, Any] = {"period": s["period"], "xaod": s["xaod"], "is_mc": s["is_mc"]}
        r_fix = job_results[(name, "fix_auto")]
        r_base = job_results[(name, "baseline")]
        rec["fix_rc"] = r_fix["rc"]
        rec["base_rc"] = r_base["rc"]
        rec["timing_log"] = r_fix["timing_line"]
        if r_fix["rc"] == 0 and r_base["rc"] == 0:
            rec["fix_vs_base"] = compare_fix_vs_base(r_fix["out_root"], r_base["out_root"])
        if s["ben_dir"] is not None:
            ben_fp = find_ben_part2(s["ben_dir"])
            if ben_fp and r_fix["rc"] == 0:
                rec["vs_ben_auto"] = compare_with_ben(r_fix["out_root"], ben_fp)
            r_par = job_results.get((name, "fix_ben_parity"))
            if ben_fp and r_par and r_par["rc"] == 0:
                rec["vs_ben_exact_override"] = compare_with_ben(r_par["out_root"], ben_fp)
        summary[name] = rec

    out_json = outdir / "validation_summary.json"
    out_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nWrote validation summary to {out_json}")

    print("\n=== FIX vs BASELINE SUMMARY ===")
    for name, rec in summary.items():
        fb = rec.get("fix_vs_base")
        if not fb:
            print(f"{name:22s}: FAILED (fix_rc={rec['fix_rc']}, base_rc={rec['base_rc']})")
            continue
        mu_ps = fb["mu_preshower_fallback"]
        mu_512 = fb["mu_v2_status512"]
        print(
            f"{name:22s}: N={fb['n_events']:5d} cp_corrupt={fb['n_baseline_clockPhase_corrupted']:3d} "
            f"fb={fb['fallback_counts']} | "
            f"PS_fb(n={mu_ps.get('n',0):3d}, fix_med={mu_ps.get('fix_med')}, base_med={mu_ps.get('base_med')}, fix<40={mu_ps.get('fix_lt40_frac')}, base<40={mu_ps.get('base_lt40_frac')}) | "
            f"V2_512(n={mu_512.get('n',0):4d}, fix_med={mu_512.get('fix_med')}, base_med={mu_512.get('base_med')})"
        )

    print("\n=== EVENT-BY-EVENT PARITY vs BEN WILSON REFERENCE FILES ===")
    for name, rec in summary.items():
        if "vs_ben_exact_override" in rec:
            p = rec["vs_ben_exact_override"]
            a = rec.get("vs_ben_auto", {})
            print(
                f"{name:22s}: matched={p['n_matched']:5d} active={p['n_active_timing']:5d} "
                f"fb_agree={p['fallback_mode_agree']} | "
                f"override_exact1e-3={p['dq0_exact_1e3_frac']} max_dq0={p['dq0_max']} | "
                f"auto_med_dq0={a.get('dq0_med')} auto_within1pC={a.get('dq0_within_1pC_frac')}"
            )


if __name__ == "__main__":
    main()
