#!/usr/bin/env python3
"""End-to-end multi-core validation of refactored NtupleDumper reduced VetoNu charge
across 2022-2026 data and all 9 official default-flux MC24 GENIE neutrino runs
from FASER_FaserMC.pdf (plus throughgoing MC24 muon control 100044).

Validates:
  1. Zero hardcoded run numbers in C++: all period/MC timing differences supplied via
     Python --timingPreset (and optional explicit CLI overrides).
  2. Event-by-event parity against Ben Wilson's MyWaveformTimeRange reference files
     across all 5 reference folders (2022, 2023, 2024 pre-CaloNu, CaloNu, post-CaloNu).
  3. Smoking-gun comparisons vs unmodified baseline NtupleDumper:
     - 2024 CaloNu Veto2 ch 22,23 status==512 (+16 ns) recovery
     - 2025-2026 4-PMT Preshower fallback recovery (eliminating -40.3 ns shift & Dig2 bias)
     - Elimination of VetoNu secondary-hit out-of-bounds corruption (clockPhase / nClusters0)
  4. All 9 official default-flux GENIE neutrino MC runs from FASER_FaserMC.pdf:
     - 2022/23 Crossing angle distribution: 200082 (epos_pi), 200083 (epos_k), 200089 (charm_central)
     - 2024 Crossing angle, FaserNu geometry: 200093 (epos_pi), 200094 (epos_k), 200100 (charm_central)
     - 2024 Crossing angle, CaloNu geometry: 200139 (epos_pi), 200140 (epos_k), 200146 (charm_central)
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
from pathlib import Path
import shutil
import subprocess
import time
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
    # --- 2022-2026 Real Data ---
    {
        "name": "2022_008090",
        "period": "2022",
        "timing_preset": "data_default",
        "xaod": "/eos/experiment/faser/rec/2022/r0021/008090/Faser-Physics-008090-00010-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/waveforms-2022/008090",
        "is_mc": False,
        "nevents": 5000,
    },
    {
        "name": "2023_011703",
        "period": "2023",
        "timing_preset": "data_default",
        "xaod": "/eos/experiment/faser/rec/2023/r0021/011703/Faser-Physics-011703-00010-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/Waveforms-2023/011703",
        "is_mc": False,
        "nevents": 5000,
    },
    {
        "name": "2024_pre_015387",
        "period": "2024_preCaloNu",
        "timing_preset": "data_default",
        "xaod": "/eos/experiment/faser/rec/2024/r0022/015387/Faser-Physics-015387-00010-r0022-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/waveforms-5FilesPerChunk/015387",
        "is_mc": False,
        "nevents": 5000,
    },
    {
        "name": "2024_early_015838",
        "period": "2024_earlyCaloNu",
        "timing_preset": "data_default",
        "xaod": "/eos/experiment/faser/rec/2024/r0022/015838/Faser-Physics-015838-00010-r0022-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/WaveformsCaloNuPeriod/015838",
        "ben_override": (35.783, 53.512),  # Ben ran 15838 with CaloNu diffs even though Veto2 was still on Dig1
        "is_mc": False,
        "nevents": 5000,
    },
    {
        "name": "2024_CaloNu_015849",
        "period": "2024_duringCaloNu",
        "timing_preset": "data_caloNu",
        "xaod": "/eos/experiment/faser/rec/2024/r0022/015849/Faser-Physics-015849-00010-r0022-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/WaveformsCaloNuPeriod/015849",
        "is_mc": False,
        "nevents": 8000,
    },
    {
        "name": "2024_post_016929",
        "period": "2024_postCaloNu",
        "timing_preset": "data_default",
        "xaod": "/eos/experiment/faser/rec/2024/r0022/016929/Faser-Physics-016929-00010-r0022-xAOD.root",
        "ben_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/waveforms_postCaloNu/016929",
        "is_mc": False,
        "nevents": 5000,
    },
    {
        "name": "2025_20012",
        "period": "2025_4PMT_PS",
        "timing_preset": "data_2025",
        "xaod": "/eos/experiment/faser/rec/2025/r0023/020012/Faser-Physics-020012-00010-r0023-xAOD.root",
        "ben_dir": None,
        "is_mc": False,
        "nevents": 8000,
    },
    {
        "name": "2025_20481",
        "period": "2025_4PMT_PS",
        "timing_preset": "data_2025",
        "xaod": "/eos/experiment/faser/rec/2025/p0014/020481/Faser-Physics-020481-00276-p0014-xAOD.root",
        "ben_dir": None,
        "is_mc": False,
        "nevents": 8000,
    },
    {
        "name": "2026_22404",
        "period": "2026_4PMT_PS",
        "timing_preset": "data_2026",
        "xaod": "/eos/experiment/faser/rec/2026/p0015/022404/Faser-Physics-022404-00010-p0015-xAOD.root",
        "ben_dir": None,
        "is_mc": False,
        "nevents": 8000,
    },
    {
        "name": "2026_23176",
        "period": "2026_4PMT_PS",
        "timing_preset": "data_2026",
        "xaod": find_xaod("/eos/experiment/faser/rec/2026/p0015/023176", "00090"),
        "ben_dir": None,
        "is_mc": False,
        "nevents": 8000,
    },
    # --- MC24 Throughgoing Muon Control ---
    {
        "name": "MC24_Muon_100044",
        "period": "MC24_Muon",
        "timing_preset": "mc",
        "xaod": "/eos/experiment/faser/data0/sim/mc24/particle_gun/100044/rec/s0013-r0023/FaserMC-MC24_PG_mupl_fasernu_5mrad_flukaE-100044-00100-00199-s0013-r0023-xAOD.root",
        "ben_dir": None,
        "is_mc": True,
        "nevents": 3000,
    },
    # --- Official Default-Flux GENIE Neutrino MC (FASER_FaserMC.pdf) ---
    # 1. 2022/23 Crossing angle distribution (p. 11)
    {
        "name": "MC24_Genie2223_200082_epos_pi",
        "period": "MC24_Genie_22_23",
        "timing_preset": "mc",
        "xaod": "/eos/experiment/fasernu-data0/faser/sim/mc24/genie/200082/rec/s0014/FaserMC-MC24_Genie_epos_pi-200082-00000-00009-s0014-xAOD.root",
        "ben_dir": None,
        "is_mc": True,
        "is_genie": True,
        "nevents": 3000,
    },
    {
        "name": "MC24_Genie2223_200083_epos_k",
        "period": "MC24_Genie_22_23",
        "timing_preset": "mc",
        "xaod": "/eos/experiment/fasernu-data0/faser/sim/mc24/genie/200083/rec/s0014/FaserMC-MC24_Genie_epos_k-200083-00000-00009-s0014-xAOD.root",
        "ben_dir": None,
        "is_mc": True,
        "is_genie": True,
        "nevents": 3000,
    },
    {
        "name": "MC24_Genie2223_200089_charm",
        "period": "MC24_Genie_22_23",
        "timing_preset": "mc",
        "xaod": "/eos/experiment/fasernu-data0/faser/sim/mc24/genie/200089/rec/s0014/FaserMC-MC24_Genie_powheg_pythia_charm_central-200089-00000-00009-s0014-xAOD.root",
        "ben_dir": None,
        "is_mc": True,
        "is_genie": True,
        "nevents": 3000,
    },
    # 2. 2024 Crossing angle distribution, FaserNu geometry (p. 13-14)
    {
        "name": "MC24_Genie24FN_200093_epos_pi",
        "period": "MC24_Genie_24_FaserNu",
        "timing_preset": "mc",
        "xaod": "/eos/experiment/fasernu-data0/faser/sim/mc24/genie/200093/rec/s0014/FaserMC-MC24_Genie_24_xing_angle_wo_caloNu_epos_pi-200093-00000-00009-s0014-xAOD.root",
        "ben_dir": None,
        "is_mc": True,
        "is_genie": True,
        "nevents": 3000,
    },
    {
        "name": "MC24_Genie24FN_200094_epos_k",
        "period": "MC24_Genie_24_FaserNu",
        "timing_preset": "mc",
        "xaod": "/eos/experiment/fasernu-data0/faser/sim/mc24/genie/200094/rec/s0014/FaserMC-MC24_Genie_24_xing_angle_wo_caloNu_epos_k-200094-00000-00009-s0014-xAOD.root",
        "ben_dir": None,
        "is_mc": True,
        "is_genie": True,
        "nevents": 3000,
    },
    {
        "name": "MC24_Genie24FN_200100_charm",
        "period": "MC24_Genie_24_FaserNu",
        "timing_preset": "mc",
        "xaod": "/eos/experiment/fasernu-data0/faser/sim/mc24/genie/200100/rec/s0014/FaserMC-MC24_Genie_24_xing_angle_wo_caloNu_powheg_pythia_charm_central-200100-00000-00009-s0014-xAOD.root",
        "ben_dir": None,
        "is_mc": True,
        "is_genie": True,
        "nevents": 3000,
    },
    # 3. 2024 Crossing angle distribution, CaloNu geometry (p. 15-16)
    {
        "name": "MC24_Genie24CN_200139_epos_pi",
        "period": "MC24_Genie_24_CaloNu",
        "timing_preset": "mc",
        "xaod": "/eos/experiment/fasernu-data0/faser/sim/mc24/genie/200139/rec/s0015/FaserMC-MC24_Genie_24_xing_angle_caloNu_epos_pi-200139-00000-00009-s0015-xAOD.root",
        "ben_dir": None,
        "is_mc": True,
        "is_genie": True,
        "nevents": 3000,
    },
    {
        "name": "MC24_Genie24CN_200140_epos_k",
        "period": "MC24_Genie_24_CaloNu",
        "timing_preset": "mc",
        "xaod": "/eos/experiment/fasernu-data0/faser/sim/mc24/genie/200140/rec/s0015/FaserMC-MC24_Genie_24_xing_angle_caloNu_epos_k-200140-00000-00009-s0015-xAOD.root",
        "ben_dir": None,
        "is_mc": True,
        "is_genie": True,
        "nevents": 3000,
    },
    {
        "name": "MC24_Genie24CN_200146_charm",
        "period": "MC24_Genie_24_CaloNu",
        "timing_preset": "mc",
        "xaod": "/eos/experiment/fasernu-data0/faser/sim/mc24/genie/200146/rec/s0015/FaserMC-MC24_Genie_24_xing_angle_caloNu_powheg_pythia_charm_central-200146-00000-00009-s0015-xAOD.root",
        "ben_dir": None,
        "is_mc": True,
        "is_genie": True,
        "nevents": 3000,
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
    if overlay and job.get("timing_preset"):
        extra_flags += ["--timingPreset", job["timing_preset"]]
    if overlay and job.get("override_diffs") is not None:
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
            (ln.strip() for ln in proc.stdout.splitlines() if "Reduced VetoNu charge configuration:" in ln),
            "",
        ),
    }


def compare_with_ben(nt_path: str, ben_path: str) -> dict[str, Any]:
    with uproot.open(nt_path) as f:
        t = f["nt"]
        cols = [
            c for c in [
                "eventID",
                "VetoNu0_reduced_charge",
                "VetoNu1_reduced_charge",
                "VetoNu0_reduced_charge_valid",
                "VetoNu1_reduced_charge_valid",
                "VetoNu_reduced_charge_fallback",
            ] if c in t
        ]
        a = t.arrays(cols, library="np")
    with uproot.open(ben_path) as f:
        bt = f["tree"]
        bcols = [
            c for c in [
                "event", "fallback_mode",
                "myVetoNu0_rawCharge35", "myVetoNu1_rawCharge35",
                "myVetoNu0_rawCharge35_valid",
            ] if c in bt
        ]
        b = bt.arrays(bcols, library="np")

    ev_nt = a["eventID"].astype(np.int64)
    ev_ben = b["event"].astype(np.int64)
    idx_map = {int(e): i for i, e in enumerate(ev_ben)}
    matched_nt, matched_ben = [], []
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
    if "fallback_mode" in b:
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
        "dq0_max": round(float(np.max(dq0)), 6) if dq0.size else None,
        "dq0_exact_1e3_frac": round(float(np.mean(dq0 < 1e-3)), 6) if dq0.size else None,
        "dq1_exact_1e3_frac": round(float(np.mean(dq1 < 1e-3)), 6) if dq1.size else None,
        "dq0_within_1pC_frac": round(float(np.mean(dq0 < 1.0)), 6) if dq0.size else None,
    }
    if "myVetoNu0_rawCharge35_valid" in b and "VetoNu0_reduced_charge_valid" in a:
        dv0 = np.abs(a["VetoNu0_reduced_charge_valid"][mi][active] - b["myVetoNu0_rawCharge35_valid"][mj][active])
        res["dq0_valid_exact_1e3_frac"] = round(float(np.mean(dv0 < 1e-3)), 6) if dv0.size else None
    return res


def analyze_fix_vs_base(fix_path: str, base_path: str, is_genie: bool = False) -> dict[str, Any]:
    with uproot.open(fix_path) as ff, uproot.open(base_path) as fb:
        tf, tb = ff["nt"], fb["nt"]
        common = [
            c for c in [
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
                "Veto20_raw_charge",
                "Veto21_raw_charge",
                "Preshower0_raw_charge",
                "Preshower1_raw_charge",
            ] if c in tf and c in tb
        ]
        extra_f = [
            c for c in [
                "VetoNu_reduced_charge_fallback",
                "VetoNu_expected_localtime",
                "VetoNu0_reduced_charge_valid",
                "VetoNu1_reduced_charge_valid",
            ] if c in tf
        ]
        af = tf.arrays(common + extra_f, library="np")
        ab = tb.arrays(common, library="np")

    assert len(af["eventID"]) == len(ab["eventID"]) and np.array_equal(af["eventID"], ab["eventID"]), "Event mismatch"
    n = len(af["eventID"])

    cp_diff = np.abs(af["clockPhase"].astype(np.float64) - ab["clockPhase"].astype(np.float64))
    n_cp_corrupt = int(np.sum(cp_diff > 1e-4))
    n_nc0_corrupt = int(np.sum(af["nClusters0"] != ab["nClusters0"]))

    fb_vals, fb_cnts = np.unique(af["VetoNu_reduced_charge_fallback"], return_counts=True)
    fb_hist = {int(k): int(v) for k, v in zip(fb_vals, fb_cnts)}

    mu = (af["VetoNu0_raw_charge"] > 40.0) & (af["VetoNu1_raw_charge"] > 40.0)
    v2_512 = np.zeros(n, dtype=bool)
    if "Veto20_status" in af and "Veto21_status" in af:
        v2_512 = (
            ((af["Veto20_status"] == 512) | (af["Veto21_status"] == 512))
            & ((af["Veto20_status"] & ~512) == 0)
            & ((af["Veto21_status"] & ~512) == 0)
        )
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
            "fix_lt30_frac": round(float(np.mean(qf < 30.0)), 4),
            "base_lt30_frac": round(float(np.mean(qb < 30.0)), 4),
            "fix_lt40_frac": round(float(np.mean(qf < 40.0)), 4),
            "base_lt40_frac": round(float(np.mean(qb < 40.0)), 4),
            "fix_over_raw_med": round(float(np.median(qf / np.maximum(qr, 1e-3))), 4),
            "base_over_raw_med": round(float(np.median(qb / np.maximum(qr, 1e-3))), 4),
        }

    # Compact charge histograms (0..250 pC, 50 bins) for plotting
    bins_q = np.linspace(0.0, 250.0, 51)
    hists = {}
    if is_genie:
        # Neutrino backsplash study: events with >=1 long track & downstream veto2/preshower activity
        nu_sel = (af["longTracks"] > 0)
        if "Veto20_raw_charge" in af and "Veto21_raw_charge" in af:
            nu_sel = nu_sel & (af["Veto20_raw_charge"] > 40.0) & (af["Veto21_raw_charge"] > 40.0)
        if "Preshower0_raw_charge" in af and "Preshower1_raw_charge" in af:
            nu_sel = nu_sel & (af["Preshower0_raw_charge"] > 2.5) & (af["Preshower1_raw_charge"] > 2.5)
        n_nu = int(nu_sel.sum())
        q_raw_max = np.maximum(af["VetoNu0_raw_charge"][nu_sel], af["VetoNu1_raw_charge"][nu_sel])
        q_red_max = np.maximum(af["VetoNu0_reduced_charge"][nu_sel], af["VetoNu1_reduced_charge"][nu_sel])
        q_base_max = np.maximum(ab["VetoNu0_reduced_charge"][nu_sel], ab["VetoNu1_reduced_charge"][nu_sel])
        fb_nu = af["VetoNu_reduced_charge_fallback"][nu_sel]
        v2_ovf_nu = (
            ((af["Veto20_status"][nu_sel] & 4) != 0) | ((af["Veto21_status"][nu_sel] & 4) != 0)
        ) if ("Veto20_status" in af and "Veto21_status" in af) else np.zeros(n_nu, dtype=bool)
        raw_veto_fail = q_raw_max >= 40.0
        n_raw_fail = int(raw_veto_fail.sum())
        n_rescued_lt30 = int(np.sum(raw_veto_fail & (q_red_max < 30.0)))
        genie_stats = {
            "n_nu_downstream_sel": n_nu,
            "pass_raw_lt40": int(np.sum(q_raw_max < 40.0)),
            "pass_fix_red_lt30": int(np.sum(q_red_max < 30.0)),
            "pass_base_red_lt30": int(np.sum(q_base_max < 30.0)),
            "n_raw_ge40_backsplash": n_raw_fail,
            "n_rescued_by_reduced_lt30": n_rescued_lt30,
            "rescue_frac": round(n_rescued_lt30 / n_raw_fail, 4) if n_raw_fail > 0 else None,
            "v2_overflow_frac": round(float(np.mean(v2_ovf_nu)), 4) if n_nu > 0 else None,
            "fb_counts_nu_sel": {int(k): int(v) for k, v in zip(*np.unique(fb_nu, return_counts=True))} if n_nu > 0 else {},
            "q_max_diff_vs_base_max": round(float(np.max(np.abs(q_red_max - q_base_max))), 4) if n_nu > 0 else 0.0,
        }
        hists["nu_raw_max"] = np.histogram(np.clip(q_raw_max, 0, 249.9), bins=bins_q)[0].tolist()
        hists["nu_fix_red_max"] = np.histogram(np.clip(q_red_max, 0, 249.9), bins=bins_q)[0].tolist()
        hists["nu_base_red_max"] = np.histogram(np.clip(q_base_max, 0, 249.9), bins=bins_q)[0].tolist()
    else:
        genie_stats = None
        if int((mu & ps_fb).sum()) > 0:
            m = mu & ps_fb
            hists["ps_fb_raw"] = np.histogram(np.clip(af["VetoNu0_raw_charge"][m], 0, 249.9), bins=bins_q)[0].tolist()
            hists["ps_fb_fix"] = np.histogram(np.clip(af["VetoNu0_reduced_charge"][m], 0, 249.9), bins=bins_q)[0].tolist()
            hists["ps_fb_base"] = np.histogram(np.clip(ab["VetoNu0_reduced_charge"][m], 0, 249.9), bins=bins_q)[0].tolist()
        if int((mu & v2_512).sum()) > 0:
            m = mu & v2_512
            hists["v2_512_raw"] = np.histogram(np.clip(af["VetoNu0_raw_charge"][m], 0, 249.9), bins=bins_q)[0].tolist()
            hists["v2_512_fix"] = np.histogram(np.clip(af["VetoNu0_reduced_charge"][m], 0, 249.9), bins=bins_q)[0].tolist()
            hists["v2_512_base"] = np.histogram(np.clip(ab["VetoNu0_reduced_charge"][m], 0, 249.9), bins=bins_q)[0].tolist()

    return {
        "n_events": n,
        "fallback_counts": fb_hist,
        "n_baseline_clockPhase_corrupted": n_cp_corrupt,
        "n_baseline_nClusters0_corrupted": n_nc0_corrupt,
        "baseline_max_clockPhase_corrupt_val": round(float(np.max(np.abs(ab["clockPhase"][cp_diff > 1e-4]))), 3) if n_cp_corrupt > 0 else 0.0,
        "mu_all": summarize_q(np.ones(n, dtype=bool)),
        "mu_v2_status512": summarize_q(v2_512),
        "mu_preshower_fallback": summarize_q(ps_fb),
        "genie_stats": genie_stats,
        "hists": hists,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlay", required=True)
    ap.add_argument("--latest-run", default="/eos/home-x/xcheng/FASER/geoRotation_calypso_latest/run")
    ap.add_argument("--env-script", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--keep-root", action="store_true", help="Keep temporary PHYS ROOT outputs")
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()

    t0 = time.time()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    selected = [s for s in SAMPLES if (not args.only) or (s["name"] in args.only)]

    jobs: list[dict[str, Any]] = []
    for s in selected:
        # 1. Refactored NtupleDumper with explicit Python --timingPreset
        jobs.append({
            "name": s["name"],
            "mode": "fix_preset",
            "tag": f"{s['name']}_fix_preset",
            "xaod": s["xaod"],
            "is_mc": s["is_mc"],
            "timing_preset": s["timing_preset"],
            "nevents": s["nevents"],
            "overlay": args.overlay,
            "latest_run": args.latest_run,
            "env_script": args.env_script,
            "work_dir": str(outdir / "_tmp_ntuples" / s["name"]),
            "override_diffs": None,
        })
        # 2. Unmodified baseline (geoRotation_calypso_latest)
        jobs.append({
            "name": s["name"],
            "mode": "baseline",
            "tag": f"{s['name']}_baseline",
            "xaod": s["xaod"],
            "is_mc": s["is_mc"],
            "timing_preset": None,
            "nevents": s["nevents"],
            "overlay": None,
            "latest_run": args.latest_run,
            "env_script": args.env_script,
            "work_dir": str(outdir / "_tmp_ntuples" / s["name"]),
            "override_diffs": None,
        })
        # 3. Optional override parity job (e.g. 15838 where Ben used CaloNu diffs)
        if s.get("ben_override") is not None:
            jobs.append({
                "name": s["name"],
                "mode": "fix_ben_override",
                "tag": f"{s['name']}_fix_ben_override",
                "xaod": s["xaod"],
                "is_mc": s["is_mc"],
                "timing_preset": s["timing_preset"],
                "nevents": s["nevents"],
                "overlay": args.overlay,
                "latest_run": args.latest_run,
                "env_script": args.env_script,
                "work_dir": str(outdir / "_tmp_ntuples" / s["name"]),
                "override_diffs": s["ben_override"],
            })

    print(f"Launching {len(jobs)} NtupleDumper validation jobs across {args.workers} workers...", flush=True)
    job_results: dict[tuple[str, str], dict[str, Any]] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for r in ex.map(run_ntuple_job, jobs):
            job_results[(r["name"], r["mode"])] = r
            st = "OK" if r["rc"] == 0 else f"FAIL(rc={r['rc']})"
            print(f"  [{st}] {r['name']:32s} {r['mode']:16s} {r['timing_line']}", flush=True)

    summary: dict[str, Any] = {}
    for s in selected:
        name = s["name"]
        rec: dict[str, Any] = {
            "period": s["period"],
            "timing_preset": s["timing_preset"],
            "xaod": s["xaod"],
            "is_mc": s["is_mc"],
            "is_genie": bool(s.get("is_genie", False)),
        }
        r_fix = job_results[(name, "fix_preset")]
        r_base = job_results[(name, "baseline")]
        rec["fix_rc"] = r_fix["rc"]
        rec["base_rc"] = r_base["rc"]
        rec["timing_log"] = r_fix["timing_line"]
        if r_fix["rc"] == 0 and r_base["rc"] == 0:
            rec["fix_vs_base"] = analyze_fix_vs_base(
                r_fix["out_root"], r_base["out_root"], is_genie=bool(s.get("is_genie", False))
            )
        if s["ben_dir"] is not None:
            ben_fp = find_ben_part2(s["ben_dir"])
            if ben_fp and r_fix["rc"] == 0:
                rec["vs_ben_preset"] = compare_with_ben(r_fix["out_root"], ben_fp)
            r_ov = job_results.get((name, "fix_ben_override"))
            if ben_fp and r_ov and r_ov["rc"] == 0:
                rec["vs_ben_override"] = compare_with_ben(r_ov["out_root"], ben_fp)
        summary[name] = rec

    out_json = outdir / "validation_summary_v2.json"
    out_json.write_text(json.dumps({
        "elapsed_s": round(time.time() - t0, 2),
        "overlay": args.overlay,
        "samples": summary,
    }, indent=2) + "\n")

    if not args.keep_root:
        tmp_dir = outdir / "_tmp_ntuples"
        for root_f in tmp_dir.rglob("*.root"):
            try:
                root_f.unlink()
            except OSError:
                pass
        print("Cleaned temporary ROOT outputs (kept logs & compact JSON summary).")

    print(f"\n=== REAL DATA & MUON MC SUMMARY ({time.time()-t0:.1f}s) ===")
    for name, rec in summary.items():
        if rec["is_genie"]:
            continue
        fb = rec.get("fix_vs_base")
        if not fb:
            print(f"{name:24s}: FAILED fix_rc={rec['fix_rc']} base_rc={rec['base_rc']}")
            continue
        mu_ps = fb["mu_preshower_fallback"]
        mu_512 = fb["mu_v2_status512"]
        ben = rec.get("vs_ben_preset")
        ben_str = (
            f"ben_exact={ben['dq0_exact_1e3_frac']:.4f} fb_agree={ben['fallback_mode_agree']}"
            if ben and ben.get("dq0_exact_1e3_frac") is not None else "-"
        )
        print(
            f"{name:24s} [{rec['timing_preset']:12s}] N={fb['n_events']:4d} cp_oob={fb['n_baseline_clockPhase_corrupted']:3d} "
            f"PS_fb(n={mu_ps.get('n',0):3d}, fix={mu_ps.get('fix_med')}, base={mu_ps.get('base_med')}, <40:{mu_ps.get('fix_lt40_frac')}/{mu_ps.get('base_lt40_frac')}) "
            f"V2_512(n={mu_512.get('n',0):4d}, fix={mu_512.get('fix_med')}, base={mu_512.get('base_med')}) | {ben_str}"
        )

    print("\n=== OFFICIAL DEFAULT-FLUX GENIE NEUTRINO MC SUMMARY (9 RUNS, FASER_FaserMC.pdf) ===")
    for name, rec in summary.items():
        if not rec["is_genie"]:
            continue
        fb = rec.get("fix_vs_base")
        if not fb:
            print(f"{name:34s}: FAILED fix_rc={rec['fix_rc']} base_rc={rec['base_rc']}")
            continue
        gs = fb["genie_stats"]
        print(
            f"{name:34s} [{rec['period']:21s}] N={fb['n_events']:4d} cp_oob={fb['n_baseline_clockPhase_corrupted']:3d} "
            f"nu_sel={gs['n_nu_downstream_sel']:4d} v2_ovf={gs['v2_overflow_frac']} "
            f"fb={gs['fb_counts_nu_sel']} | raw>=40:{gs['n_raw_ge40_backsplash']:3d} -> red<30:{gs['n_rescued_by_reduced_lt30']:3d} "
            f"(rescue={gs['rescue_frac']})"
        )


if __name__ == "__main__":
    main()
