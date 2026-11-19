#!/usr/bin/env python3
"""Parallel xAOD waveform scanner driver invoking compiled scan_xaod_waveform."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
from pathlib import Path
import subprocess
from typing import Any

import numpy as np


def pick_xaod(run_dir: Path) -> str | None:
    fs = sorted(run_dir.glob("*xAOD.root"))
    if not fs:
        return None
    return str(fs[len(fs) // 2])


def run_one(args_tuple: tuple[str, str, int, str, int]) -> dict[str, Any]:
    binary, period, run, fpath, max_events = args_tuple
    proc = subprocess.run([binary, fpath, str(max_events)], capture_output=True, text=True)
    line = ""
    for ln in reversed(proc.stdout.splitlines()):
        ln = ln.strip()
        if ln.startswith("{") and ln.endswith("}"):
            line = ln
            break
    if not line:
        return {"period": period, "run": run, "file": fpath, "error": f"rc={proc.returncode} {proc.stderr[-200:]}"}
    rec = json.loads(line)
    rec["period"] = period
    rec["run"] = run
    rec["file"] = fpath
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--max-events", type=int, default=40000)
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()

    tasks: list[tuple[str, str, int, str, int]] = []
    add = lambda p, r, fp: tasks.append((args.binary, p, r, fp, args.max_events))

    # 2022 & 2023
    for label, rdir in [
        ("2022_r0021", Path("/eos/experiment/faser/rec/2022/r0021/007729")),
        ("2022_r0022", Path("/eos/experiment/faser/rec/2022/r0022/008090")),
        ("2023_r0021", Path("/eos/experiment/faser/rec/2023/r0021/010417")),
        ("2023_r0021", Path("/eos/experiment/faser/rec/2023/r0021/011703")),
    ]:
        fp = pick_xaod(rdir)
        if fp:
            add(label, int(rdir.name), fp)

    # 2024: pre-CaloNu, CaloNu-early (15821-15844), CaloNu-main (15845-16927), post-CaloNu
    rec24_22 = Path("/eos/experiment/faser/rec/2024/r0022")
    rec24_23 = Path("/eos/experiment/faser/rec/2024/r0023")
    for run in [14587, 15200, 15814, 15820, 15838, 15839, 15840, 15844, 15849, 15850, 16020, 16300, 16600, 16901, 16923, 16929, 16932, 17147, 17230]:
        for base, tag in [(rec24_22, "r0022"), (rec24_23, "r0023")]:
            rd = base / f"{run:06d}"
            if rd.is_dir():
                fp = pick_xaod(rd)
                if fp:
                    if run < 15821:
                        p = f"2024_pre_{tag}"
                    elif run < 15845:
                        p = f"2024_cn_early_{tag}"
                    elif run < 16928:
                        p = f"2024_CaloNu_{tag}"
                    else:
                        p = f"2024_post_{tag}"
                    add(p, run, fp)
                break
    for run in [13832, 15848, 17605]:
        rd = rec24_23 / f"{run:06d}"
        if rd.is_dir():
            fp = pick_xaod(rd)
            if fp:
                p = "2024_pre_r0023" if run < 15821 else ("2024_CaloNu_r0023" if run < 16928 else "2024_post_r0023")
                add(p, run, fp)

    # 2025
    for base, tag in [
        (Path("/eos/experiment/faser/rec/2025/p0014"), "2025_p0014"),
        (Path("/eos/experiment/faser/rec/2025/r0023"), "2025_r0023"),
    ]:
        if base.is_dir():
            rds = sorted([d for d in base.iterdir() if d.is_dir() and d.name.isdigit()])
            idxs = np.linspace(0, len(rds) - 1, min(6, len(rds)), dtype=int)
            for idx in sorted(set(idxs.tolist())):
                fp = pick_xaod(rds[idx])
                if fp:
                    add(tag, int(rds[idx].name), fp)

    # 2026 FULL YEAR (429 runs: 22126..24068)
    rec26 = Path("/eos/experiment/faser/rec/2026/p0015")
    if rec26.is_dir():
        rds = sorted([d for d in rec26.iterdir() if d.is_dir() and d.name.isdigit()])
        idxs = np.linspace(0, len(rds) - 1, min(20, len(rds)), dtype=int)
        for idx in sorted(set(idxs.tolist())):
            fp = pick_xaod(rds[idx])
            if fp:
                add("2026_p0015", int(rds[idx].name), fp)

    # MC
    mc_xaod_dirs = [
        ("MC24_PG_100044_r23", 100044, Path("/eos/experiment/faser/data0/sim/mc24/particle_gun/100044/rec/s0013-r0023")),
        ("MC24_PG_100049_s15", 100049, Path("/eos/experiment/faser/data0/sim/mc24/particle_gun/100049/rec/s0015")),
        ("MC24_fluka_210010",  210010, Path("/eos/experiment/faser/data0/sim/mc24/fluka/210010/rec/s0013-r0019")),
        ("MC24_genie_200139",  200139, Path("/eos/experiment/faser/data0/sim/mc24/genie/200139/rec/caloNu")),
        ("MC22_fluka_210001",  210001, Path("/eos/experiment/faser/data0/sim/mc22/fluka/210001/rec/r0013")),
    ]
    for label, dsid, d in mc_xaod_dirs:
        if d.is_dir():
            fp = pick_xaod(d)
            if fp:
                add(label, dsid, fp)

    print(f"Scanning {len(tasks)} xAOD files with {args.workers} workers...")
    results: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for res in ex.map(run_one, tasks):
            results.append(res)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(results, indent=2) + "\n")

    print(f"\n{'Period':18s} {'Run':>6s} {'Nev':>6s} {'VNw%':>5s} {'ge3':>5s} {'swp':>4s} | {'V2-VN':>8s} {'V2_0':>8s} {'V2_512c':>8s} {'V2_512r':>8s} | {'PSall':>8s} {'PStop':>8s} {'PSbot':>8s} {'PSb512c':>8s} | ch")
    print("-" * 145)
    for r in sorted(results, key=lambda x: (x["period"], x["run"])):
        if "error" in r:
            print(f"{r['period']:18s} {r['run']:6d} ERROR: {r['error']}")
            continue
        w_pct = 100.0 * r["n_vn_nonempty_wave"] / max(1, r["n_vn_hits"])
        f8 = lambda d: f"{d['peak']:8.3f}" if d and d.get("peak") is not None else f"{'-':>8s}"
        print(
            f"{r['period']:18s} {r['run']:6d} {r['n_events']:6d} {w_pct:5.1f} {r['n_vn_ge3_hits']:5d} {r['n_vn_ch5_before_ch4']:4d} | "
            f"{f8(r['veto2_minus_vetonu'])} {f8(r['veto2_st0_minus_vetonu'])} {f8(r['veto2_st512_corr_minus_vetonu'])} {f8(r['veto2_st512_raw_minus_vetonu'])} | "
            f"{f8(r['ps_all_minus_vetonu'])} {f8(r['ps_top_minus_vetonu'])} {f8(r['ps_bot_minus_vetonu'])} {f8(r['ps_bot_st512_corr_minus_vetonu'])} | "
            f"v2={r['v2_channels']} ps={r['ps_channels']}"
        )


if __name__ == "__main__":
    main()
