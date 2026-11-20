#!/usr/bin/env python3
"""Freeze the 50-file (5-period) Ben-vs-geoRotation_calypso_latest_fix_v1 reconstruction
task manifest, select the least-loaded compatible bigbird schedd, and optionally submit."""

from __future__ import annotations

import argparse
import datetime
import glob
import json
import os
from pathlib import Path
import shutil
import subprocess

PERIOD_SPECS = [
    {
        "period": "2022",
        "run": 7961,
        "segments": list(range(10, 20)),  # Ben part2 + part3 (all inside stable colliding beams)
        "ben_parts": [2, 3],
        "raw_dir": "/eos/experiment/faser/raw/2022/007961",
        "off_phys": "/eos/experiment/faser/data0/phys/2022_back/r0021/007961/Faser-Physics-007961-00000-00048-PHYS.root",
        "ben_wf_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/waveforms-2022/007961",
        "grl_json": "/cvmfs/faser.cern.ch/repo/sw/runlist/v10/faser_goodrunlist_2022_physics.json",
        "cond_tag": "OFLCOND-FASER-04",
        "timing_preset": "data_default",
    },
    {
        "period": "2023",
        "run": 11523,
        "segments": list(range(30, 40)),  # Ben part6 + part7 (stable colliding beams)
        "ben_parts": [6, 7],
        "raw_dir": "/eos/experiment/faser/staged/raw/2023/011523",
        "off_phys": "/eos/experiment/faser/data0/phys/2023_back/r0021/011523/Faser-Physics-011523-00000-00049-PHYS.root",
        "ben_wf_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/Waveforms-2023/011523",
        "grl_json": "/cvmfs/faser.cern.ch/repo/sw/runlist/v10/faser_goodrunlist_2023_physics.json",
        "cond_tag": "OFLCOND-FASER-04",
        "timing_preset": "data_default",
    },
    {
        "period": "2024_pre",
        "run": 15387,
        "segments": list(range(10, 20)),  # Ben part2 + part3
        "ben_parts": [2, 3],
        "raw_dir": "/eos/experiment/faser/raw/2024/015387",
        "off_phys": "/eos/experiment/faser/data0/phys/2024_back/r0022/015387/Faser-Physics-015387-00000-00049-r0022-PHYS.root",
        "ben_wf_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/waveforms-5FilesPerChunk/015387",
        "grl_json": "/cvmfs/faser.cern.ch/repo/sw/runlist/v10/faser_goodrunlist_2024_physics.json",
        "cond_tag": "OFLCOND-FASER-05",
        "timing_preset": "data_default",
    },
    {
        "period": "2024_caloNu",
        "run": 15849,
        "segments": list(range(10, 20)),  # Ben part2 + part3 (Dig2 Veto2 ch 22,23 + status==512)
        "ben_parts": [2, 3],
        "raw_dir": "/eos/experiment/faser/raw/2024/015849",
        "off_phys": "/eos/experiment/faser/data0/phys/2024_back/r0022/015849/Faser-Physics-015849-00000-00035-r0022-PHYS.root",
        "ben_wf_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/WaveformsCaloNuPeriod/015849",
        "grl_json": "/cvmfs/faser.cern.ch/repo/sw/runlist/v10/faser_goodrunlist_2024_physics.json",
        "cond_tag": "OFLCOND-FASER-05",
        "timing_preset": "data_caloNu",
    },
    {
        "period": "2024_post",
        "run": 16963,
        "segments": list(range(25, 35)),  # Ben part5 + part6 (stable colliding beams)
        "ben_parts": [5, 6],
        "raw_dir": "/eos/experiment/faser/raw/2024/016963",
        "off_phys": "/eos/experiment/faser/data0/phys/2024_back/r0022/016963/Faser-Physics-016963-00000-00049-r0022-PHYS.root",
        "ben_wf_dir": "/eos/user/b/bewilson/El9CalypsoForWaveForm/waveforms_postCaloNu/016963",
        "grl_json": "/cvmfs/faser.cern.ch/repo/sw/runlist/v10/faser_goodrunlist_2024_physics.json",
        "cond_tag": "OFLCOND-FASER-05",
        "timing_preset": "data_default",
    },
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default="run_20261120_v1")
    parser.add_argument("--flavour", default="workday")
    parser.add_argument("--memory-mb", type=int, default=3500)
    parser.add_argument("--disk-gb", type=int, default=8)
    parser.add_argument("--submit", action="store_true", help="Actually submit to HTCondor")
    parser.add_argument("--exclude-schedd", action="append", default=[])
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    out_root = Path(
        f"/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/{args.run_id}"
    )
    shards_root = out_root / "shards"
    logs_dir = out_root / "condor_logs"
    code_dir = out_root / "frozen_worker"
    out_root.mkdir(parents=True, exist_ok=True)
    shards_root.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    code_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    period_manifest = {}
    for spec in PERIOD_SPECS:
        rstr = f"{spec['run']:06d}"
        assert os.path.isfile(spec["off_phys"]), f"Missing official PHYS: {spec['off_phys']}"
        assert os.path.isfile(spec["grl_json"]), f"Missing GRL: {spec['grl_json']}"
        ben_files = []
        for pt in spec["ben_parts"]:
            matches = sorted(glob.glob(f"{spec['ben_wf_dir']}/*.part{pt}.{rstr}.root"))
            assert len(matches) == 1, f"Expected 1 Ben part{pt} in {spec['ben_wf_dir']}, found {matches}"
            ben_files.append(matches[0])
        period_manifest[spec["period"]] = {
            **spec,
            "ben_wf_files": ben_files,
        }
        for seg in spec["segments"]:
            raw_path = f"{spec['raw_dir']}/Faser-Physics-{rstr}-{seg:05d}.raw"
            assert os.path.isfile(raw_path), f"Missing RAW file: {raw_path}"
            task_id = f"{spec['period']}_{rstr}_{seg:05d}"
            shard_dir = shards_root / task_id
            shard_dir.mkdir(parents=True, exist_ok=True)
            tasks.append({
                "task_id": task_id,
                "period": spec["period"],
                "run": spec["run"],
                "segment": f"{seg:05d}",
                "raw_path": raw_path,
                "raw_url": f"root://eospublic.cern.ch/{raw_path}",
                "cond_tag": spec["cond_tag"],
                "timing_preset": spec["timing_preset"],
                "output_dir": str(shard_dir),
            })

    tasks_json = out_root / "tasks.json"
    tasks_json.write_text(json.dumps(tasks, indent=2) + "\n")
    (out_root / "periods.json").write_text(json.dumps(period_manifest, indent=2) + "\n")

    frozen_worker = code_dir / "ben_comparison_worker.sh"
    shutil.copy2(repo / "scripts/ben_comparison_worker.sh", frozen_worker)
    frozen_worker.chmod(0o755)

    # Query live schedd load and pick the least loaded ordinary bigbird schedd
    load = subprocess.check_output(
        ["condor_status", "-schedd", "-af", "Name", "TotalRunningJobs", "TotalIdleJobs", "MaxJobsRunning"],
        text=True,
    )
    (out_root / "schedd_load.txt").write_text(load)
    candidates = []
    for line in load.splitlines():
        parts = line.split()
        if len(parts) != 4:
            continue
        name, running, idle, cap = parts
        short = name.split(".")[0]
        if (
            short.startswith("bigbird")
            and len(short) == 9
            and short[7:9].isdigit()
            and int(short[7:9]) < 100
            and name not in args.exclude_schedd
        ):
            candidates.append((int(running) + int(idle), int(idle), int(running), int(cap), name))
    candidates.sort()
    assert candidates, "No compatible ordinary bigbird schedd found"
    schedd = candidates[0][4]
    queue_totals = subprocess.check_output(["condor_q", "-name", schedd, "-totals"], text=True)
    (out_root / "schedd_queue.txt").write_text(queue_totals)

    submit_dir = Path.home() / "condor/geoRotation_test" / f"stage04_ben_comparison_{args.run_id}"
    submit_dir.mkdir(parents=True, exist_ok=True)
    sub_file = submit_dir / "jobs.sub"
    sub_file.write_text(
        f"""universe = vanilla
executable = {frozen_worker}
arguments = {tasks_json} $(Process)
request_cpus = 1
request_memory = {args.memory_mb}MB
request_disk = {args.disk_gb}GB
+WantOS = "el9"
+JobFlavour = "{args.flavour}"
+GeoRotationStage = "S04_BenComparison"
getenv = False
should_transfer_files = YES
transfer_output_files = ""
output = {logs_dir}/job.$(Cluster).$(Process).out
error = {logs_dir}/job.$(Cluster).$(Process).err
log = {logs_dir}/condor_events.log
queue {len(tasks)}
"""
    )

    manifest = {
        "time_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "run_id": args.run_id,
        "n_tasks": len(tasks),
        "schedd": schedd,
        "top_schedds": [
            {"name": c[4], "running": c[2], "idle": c[1], "total": c[0], "capacity": c[3]}
            for c in candidates[:5]
        ],
        "selection_basis": "Minimum (TotalRunningJobs + TotalIdleJobs) among ordinary EL9/EOS/CVMFS bigbird schedds",
        "submit_file": str(sub_file),
        "tasks_json": str(tasks_json),
        "frozen_worker": str(frozen_worker),
        "flavour": args.flavour,
        "memory_mb": args.memory_mb,
        "disk_gb": args.disk_gb,
        "submitted": False,
        "submission_stdout": None,
    }
    if args.submit:
        res = subprocess.check_output(["condor_submit", "-name", schedd, str(sub_file)], text=True)
        manifest["submitted"] = True
        manifest["submission_stdout"] = res.strip()

    (out_root / "submission.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
