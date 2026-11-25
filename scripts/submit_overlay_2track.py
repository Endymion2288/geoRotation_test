#!/usr/bin/env python3
"""Freeze Stage 08 two-track overlay manifest, compile frame normalizer, select
the least-loaded compatible EL9 bigbird schedd, and prepare/submit HTCondor jobs.

Supports:
  --dry-run                Prepare frozen manifest, binaries, submit files, and schedd ranking without submitting
  --submit                 Actually submit prepared jobs via condor_submit -name <schedd>
  --resubmit-missing       Scan shards for missing/non-ok summary.json or empty *_nt.root and resubmit only those
  --include-signal-110038  Also prepare/submit the 10-shard A' (100 MeV, 1e-5) 110038 signal reconstruction
  --subset {all,tier2_36}  Run full 323 overlay shards or the 36-shard attribution control subset
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any, Iterable, Sequence

DEFAULT_OVERLAY_ROOT = Path("/eos/experiment/faser/data0/overlay")
DEFAULT_SIGNAL_110038_DIR = Path("/eos/experiment/faser/data0/sim/mc24/foresee/110038/rdo/master_dev")
DEFAULT_CALYPSO_SRC = Path("/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1")
DEFAULT_CALYPSO_RUN = Path("/eos/home-x/xcheng/FASER/geoRotation_test_results/calypso_sazhang_mr/run")
DEFAULT_CALYPSO_COMMIT = "218813c2bcfc656fb8c395a167affcc033250867"
DEFAULT_ROTATED_MATERIAL_MAP = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/calypso_master_sync_v1/material/material-maps-alma9-rotated-rowmajor.json"
)
DEFAULT_UNROTATED_MATERIAL_MAP = Path(
    "/cvmfs/faser.cern.ch/repo/sw/database/DBRelease/current/acts/material-maps-alma9.json"
)
DEFAULT_RESULTS_BASE = Path("/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track")

OVERLAY_DATASETS: Sequence[dict[str, Any]] = (
    {"name": "008749", "is_mc": False, "label": "2022 Data Run 008749"},
    {"name": "010698", "is_mc": False, "label": "2023 Data Run 010698"},
    {"name": "Pos_008749-Neg_010698", "is_mc": False, "label": "Cross-run Data (008749 Pos + 010698 Neg)"},
    {"name": "10004X", "is_mc": True, "label": "MC24 Particle Gun 100047/100048"},
)

# 18 pairs x 2 datasets (008749, 10004X) = 36 shards for optional Tier-2 control
TIER2_CONTROL_DATASETS = {"008749", "10004X"}
TIER2_CONTROL_PAIRS = {
    (0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (7, 7), (8, 8),
    (0, 3), (1, 2), (5, 6), (7, 8), (0, 8), (1, 7), (2, 8), (3, 7), (4, 7),
}

RDO_BIN_RE = re.compile(r"^myOverlay_PosBin(\d+)_NegBin(\d+)-RDO\.pool\.root$")


def acts_to_faser(vx: float, vy: float, vz: float) -> tuple[float, float, float]:
    """Exact inverse rotation R_{A->F}: (x_A, y_A, z_A) -> (-z_A, y_A, x_A)."""
    return (-vz, vy, vx)


def faser_to_acts(vx: float, vy: float, vz: float) -> tuple[float, float, float]:
    """Exact forward rotation R_{F->A}: (x_F, y_F, z_F) -> (z_F, y_F, -x_F)."""
    return (vz, vy, -vx)


def parse_bool(value: str) -> bool:
    norm = value.strip().lower()
    if norm in ("1", "true", "yes", "on"):
        return True
    if norm in ("0", "false", "no", "off"):
        return False
    raise argparse.ArgumentTypeError(f"Expected boolean string, got {value!r}")


def rank_bigbird_schedds(
    schedd_output: str,
    exclude_schedds: Iterable[str] = (),
) -> list[dict[str, Any]]:
    """Filter and rank compatible ordinary 2-digit EL9 bigbirdNN.cern.ch schedds."""
    excluded = set(exclude_schedds)
    candidates: list[tuple[int, int, int, int, str, str]] = []
    for raw_line in schedd_output.splitlines():
        parts = raw_line.split()
        if len(parts) not in (4, 5):
            continue
        if len(parts) == 5:
            name, machine, running_s, idle_s, cap_s = parts
        else:
            name, running_s, idle_s, cap_s = parts
            machine = name
        short = name.split(".")[0]
        if not (short.startswith("bigbird") and len(short) == 9 and short[7:9].isdigit()):
            continue
        if int(short[7:9]) >= 100 or name in excluded:
            continue
        try:
            running = int(running_s)
            idle = int(idle_s)
            cap = int(cap_s)
        except ValueError:
            continue
        if cap <= 0:
            continue
        total = running + idle
        candidates.append((total, idle, running, cap, name, machine))
    candidates.sort()
    return [
        {
            "name": c[4],
            "machine": c[5],
            "running": c[2],
            "idle": c[1],
            "total": c[0],
            "capacity": c[3],
        }
        for c in candidates
    ]


def build_overlay_manifest(
    overlay_root: Path,
    out_root: Path,
    *,
    cond_tag: str,
    geom: str,
    geometry_rotation: bool,
    material_map: Path,
    calypso_source: Path,
    calypso_run: Path,
    calypso_commit: str,
    overlay_reco_script: Path,
    acts_to_faser_bin: Path,
    subset: str = "all",
    nevents: int = -1,
    seed_cov_phi: float | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Discover all existing OverlaidRDOs and build the frozen shard manifest."""
    tasks: list[dict[str, Any]] = []
    per_dataset_summary: dict[str, Any] = {}

    for spec in OVERLAY_DATASETS:
        ds = spec["name"]
        rdo_dir = overlay_root / ds / "OverlaidRDOs"
        jul8_dir = overlay_root / ds / "OverlaidNTPs"
        may30_dir = overlay_root / ds / "OverlaidNTPs_bug"

        if not rdo_dir.is_dir():
            raise FileNotFoundError(f"Missing OverlaidRDOs directory: {rdo_dir}")

        found_pairs: dict[tuple[int, int], Path] = {}
        for entry in sorted(rdo_dir.iterdir()):
            if not entry.is_file():
                continue
            m = RDO_BIN_RE.match(entry.name)
            if m:
                found_pairs[(int(m.group(1)), int(m.group(2)))] = entry

        ds_total = 0
        ds_common = 0
        missing_rdo_pairs = [
            [i, j] for i in range(9) for j in range(9) if (i, j) not in found_pairs
        ]
        missing_baseline_pairs: list[list[int]] = []

        for (pos_bin, neg_bin), rdo_path in sorted(found_pairs.items()):
            if subset == "tier2_36":
                if ds not in TIER2_CONTROL_DATASETS or (pos_bin, neg_bin) not in TIER2_CONTROL_PAIRS:
                    continue

            ntp_name = f"myOverlay_PosBin{pos_bin}_NegBin{neg_bin}-RDO.pool-ntuple.root"
            jul8_path = jul8_dir / ntp_name
            may30_path = may30_dir / ntp_name
            has_jul8 = jul8_path.is_file()
            has_may30 = may30_path.is_file()
            in_common = has_jul8 and has_may30
            if not in_common:
                missing_baseline_pairs.append([pos_bin, neg_bin])

            task_id = f"{ds}_PosBin{pos_bin}_NegBin{neg_bin}"
            shard_dir = out_root / "shards" / ds / f"PosBin{pos_bin}_NegBin{neg_bin}"
            output_ntp = shard_dir / f"{task_id}_nt.root"

            tasks.append(
                {
                    "task_id": task_id,
                    "dataset": ds,
                    "is_mc": bool(spec["is_mc"]),
                    "pos_bin": pos_bin,
                    "neg_bin": neg_bin,
                    "input_rdo": str(rdo_path),
                    "baseline_jul8_ntp": str(jul8_path) if has_jul8 else None,
                    "baseline_may30_ntp": str(may30_path) if has_may30 else None,
                    "in_common_baseline_mask": in_common,
                    "output_dir": str(shard_dir),
                    "output_ntuple": str(output_ntp),
                    "cond_tag": cond_tag,
                    "geom": geom,
                    "geometry_rotation": geometry_rotation,
                    "material_map": str(material_map),
                    "calypso_source": str(calypso_source),
                    "calypso_run": str(calypso_run),
                    "calypso_commit": calypso_commit,
                    "overlay_reco_script": str(overlay_reco_script),
                    "acts_to_faser_bin": str(acts_to_faser_bin),
                    "nevents": nevents,
                    "seed_cov_phi": seed_cov_phi,
                }
            )
            ds_total += 1
            if in_common:
                ds_common += 1

        per_dataset_summary[ds] = {
            "is_mc": bool(spec["is_mc"]),
            "rdo_shards": ds_total,
            "common_baseline_shards": ds_common,
            "extra_non_baseline_shards": ds_total - ds_common,
            "missing_rdo_pairs": missing_rdo_pairs,
            "missing_baseline_pairs_among_rdos": missing_baseline_pairs,
        }

    meta = {
        "schema_version": 1,
        "subset": subset,
        "total_shards": len(tasks),
        "common_baseline_shards": sum(1 for t in tasks if t["in_common_baseline_mask"]),
        "extra_non_baseline_shards": sum(1 for t in tasks if not t["in_common_baseline_mask"]),
        "datasets": per_dataset_summary,
    }
    return tasks, meta


def build_signal_110038_tasks(
    signal_rdo_dir: Path,
    out_root: Path,
    *,
    cond_tag: str,
    geom: str,
    material_map: Path,
    calypso_source: Path,
    calypso_run: Path,
    calypso_commit: str,
    nevents: int = -1,
) -> list[dict[str, Any]]:
    if not signal_rdo_dir.is_dir():
        raise FileNotFoundError(f"Missing signal 110038 RDO dir: {signal_rdo_dir}")
    rdo_files = sorted(p for p in signal_rdo_dir.iterdir() if p.is_file() and p.name.endswith("-RDO.root"))
    if not rdo_files:
        raise FileNotFoundError(f"No *-RDO.root files found in {signal_rdo_dir}")

    tasks: list[dict[str, Any]] = []
    for idx, rdo_path in enumerate(rdo_files):
        task_id = f"signal_110038_{idx:02d}"
        shard_dir = out_root / "signal_110038" / task_id
        tasks.append(
            {
                "task_id": task_id,
                "shard_index": idx,
                "dataset": "signal_110038",
                "is_mc": True,
                "input_rdo": str(rdo_path),
                "output_dir": str(shard_dir),
                "output_phys": str(shard_dir / f"{task_id}_phys.root"),
                "cond_tag": cond_tag,
                "geom": geom,
                "material_map": str(material_map),
                "calypso_source": str(calypso_source),
                "calypso_run": str(calypso_run),
                "calypso_commit": calypso_commit,
                "nevents": nevents,
            }
        )
    return tasks


def shard_is_complete(task: dict[str, Any], output_key: str) -> bool:
    out_file = Path(task[output_key])
    summary_file = Path(task["output_dir"]) / "summary.json"
    if not (out_file.is_file() and out_file.stat().st_size > 0 and summary_file.is_file()):
        return False
    try:
        data = json.loads(summary_file.read_text())
    except Exception:
        return False
    return data.get("status") == "ok"


def compile_acts_to_faser(src_cpp: Path, out_bin: Path, calypso_source: Path) -> None:
    """Compile overlay_acts_to_faser inside Athena 24.0.41 CVMFS environment so
    the frozen binary matches the exact GCC 13 / ROOT ABI on every EL9 worker."""
    out_bin.parent.mkdir(parents=True, exist_ok=True)
    tmp_env = out_bin.parent / f".compile_env_{os.getpid()}"
    tmp_env.mkdir(parents=True, exist_ok=True)
    try:
        build_script = f"""
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
cd "{tmp_env}"
source "${{ATLAS_LOCAL_ROOT_BASE}}/user/atlasLocalSetup.sh" --quiet
source "${{AtlasSetup}}/scripts/asetup.sh" --input="{calypso_source}/asetup.faser" Athena,24.0.41
set -euo pipefail
g++ -O3 -std=c++17 -Wall -Wextra "{src_cpp}" -o "{out_bin}" $(root-config --cflags --libs)
"""
        subprocess.run(["bash", "-c", build_script], check=True)
        out_bin.chmod(0o755)
    finally:
        shutil.rmtree(tmp_env, ignore_errors=True)


def query_best_schedd(out_root: Path, exclude_schedds: Sequence[str], suffix: str = "") -> tuple[str, list[dict[str, Any]], str]:
    load_txt = subprocess.check_output(
        ["condor_status", "-schedd", "-af", "Name", "Machine", "TotalRunningJobs", "TotalIdleJobs", "MaxJobsRunning"],
        text=True,
    )
    (out_root / f"schedd_load{suffix}.txt").write_text(load_txt)
    ranked = rank_bigbird_schedds(load_txt, exclude_schedds=exclude_schedds)
    if not ranked:
        raise RuntimeError("No compatible ordinary 2-digit EL9 bigbird schedd found")
    best = ranked[0]["name"]
    queue_txt = subprocess.check_output(["condor_q", "-name", best, "-totals"], text=True)
    (out_root / f"schedd_queue{suffix}.txt").write_text(queue_txt)
    return best, ranked, queue_txt


def write_submit_file(
    sub_path: Path,
    *,
    executable: Path,
    tasks_json: Path,
    n_tasks: int,
    logs_dir: Path,
    log_prefix: str,
    stage_label: str,
    cpus: int,
    memory_mb: int,
    disk_gb: int,
    flavour: str,
) -> None:
    sub_path.parent.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    content = f"""universe = vanilla
executable = {executable}
arguments = {tasks_json} $(Process)
request_cpus = {cpus}
request_memory = {memory_mb}MB
request_disk = {disk_gb}GB
+WantOS = "el9"
+JobFlavour = "{flavour}"
+GeoRotationStage = "{stage_label}"
getenv = False
should_transfer_files = YES
transfer_output_files = ""
output = {logs_dir}/{log_prefix}.$(Cluster).$(Process).out
error = {logs_dir}/{log_prefix}.$(Cluster).$(Process).err
log = {logs_dir}/{log_prefix}.events.log
queue {n_tasks}
"""
    sub_path.write_text(content)


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare and submit Stage 08 two-track overlay HTCondor jobs.")
    parser.add_argument("--run-id", default="run_20261125_v1", help="Unique run identifier")
    parser.add_argument("--out-root", type=Path, default=None, help="Override EOS output root")
    parser.add_argument("--submit-dir", type=Path, default=None, help="Override AFS submit dir")
    parser.add_argument("--overlay-root", type=Path, default=DEFAULT_OVERLAY_ROOT)
    parser.add_argument("--signal-rdo-dir", type=Path, default=DEFAULT_SIGNAL_110038_DIR)
    parser.add_argument("--calypso-source", type=Path, default=DEFAULT_CALYPSO_SRC)
    parser.add_argument("--calypso-run", type=Path, default=DEFAULT_CALYPSO_RUN)
    parser.add_argument("--calypso-commit", default=DEFAULT_CALYPSO_COMMIT)
    parser.add_argument("--geom", default="TI12MC04")
    parser.add_argument("--cond-tag", default="OFLCOND-FASER-05")
    parser.add_argument("--geometry-rotation", type=parse_bool, default=True)
    parser.add_argument(
        "--seed-cov-phi",
        type=float,
        default=None,
        help="Optional override for CircleFitTrackSeedTool.covPhi [rad^2] (default: None -> unmodified CKF2Cfg default)",
    )
    parser.add_argument("--material-map", type=Path, default=None)
    parser.add_argument("--subset", choices=["all", "tier2_36"], default="all")
    parser.add_argument("--nevents", type=int, default=-1)
    parser.add_argument("--flavour", default="longlunch")
    parser.add_argument("--memory-mb", type=int, default=3000)
    parser.add_argument("--disk-gb", type=int, default=4)
    parser.add_argument("--include-signal-110038", action="store_true")
    parser.add_argument("--signal-flavour", default="workday")
    parser.add_argument("--signal-memory-mb", type=int, default=3500)
    parser.add_argument("--signal-disk-gb", type=int, default=12)
    parser.add_argument("--exclude-schedd", action="append", default=[])
    parser.add_argument("--dry-run", action="store_true", help="Prepare frozen artifacts & submit files without submitting")
    parser.add_argument("--submit", action="store_true", help="Submit prepared jobs to HTCondor")
    parser.add_argument("--resubmit-missing", action="store_true", help="Resubmit only incomplete/failed shards")
    args = parser.parse_args()

    if args.submit and args.dry_run:
        parser.error("Cannot specify both --submit and --dry-run")

    repo = Path(__file__).resolve().parents[1]
    out_root = args.out_root or (DEFAULT_RESULTS_BASE / args.run_id)
    submit_dir = args.submit_dir or (Path.home() / "condor" / "geoRotation_test" / f"stage08_overlay_2track_{args.run_id}")
    material_map = args.material_map or (
        DEFAULT_ROTATED_MATERIAL_MAP if args.geometry_rotation else DEFAULT_UNROTATED_MATERIAL_MAP
    )

    if not material_map.is_file():
        raise FileNotFoundError(f"Material map not found: {material_map}")
    if not (args.calypso_run / "setup.sh").is_file():
        raise FileNotFoundError(f"Missing Calypso run setup.sh: {args.calypso_run / 'setup.sh'}")

    frozen_code = out_root / "frozen_code"
    frozen_bin = out_root / "frozen_bin"
    logs_dir = out_root / "condor_logs"
    for d in (out_root, frozen_code, frozen_bin, logs_dir, submit_dir):
        d.mkdir(parents=True, exist_ok=True)

    frozen_reco_script = frozen_code / "overlay_reco_ntuple.py"
    frozen_worker = frozen_code / "overlay_worker.sh"
    frozen_signal_worker = frozen_code / "overlay_signal_worker.sh"
    frozen_acts_bin = frozen_bin / "overlay_acts_to_faser"

    # Freeze scripts & compile C++ frame converter unless resubmitting against existing frozen artifacts
    if not args.resubmit_missing or not frozen_acts_bin.is_file():
        shutil.copy2(repo / "configs" / "overlay_reco_ntuple.py", frozen_reco_script)
        frozen_reco_script.chmod(0o755)
        shutil.copy2(repo / "scripts" / "overlay_worker.sh", frozen_worker)
        frozen_worker.chmod(0o755)
        shutil.copy2(repo / "scripts" / "overlay_signal_worker.sh", frozen_signal_worker)
        frozen_signal_worker.chmod(0o755)
        compile_acts_to_faser(repo / "analysis" / "overlay_acts_to_faser.cpp", frozen_acts_bin, args.calypso_source)

    tasks_json = out_root / "tasks.json"
    manifest_json = out_root / "manifest.json"
    signal_tasks_json = out_root / "signal_110038_tasks.json"

    if args.resubmit_missing:
        if not tasks_json.is_file():
            raise FileNotFoundError(f"--resubmit-missing requires existing {tasks_json}")
        all_overlay_tasks = json.loads(tasks_json.read_text())
        overlay_to_run = [t for t in all_overlay_tasks if not shard_is_complete(t, "output_ntuple")]
        signal_to_run: list[dict[str, Any]] = []
        if args.include_signal_110038 or signal_tasks_json.is_file():
            if signal_tasks_json.is_file():
                all_sig = json.loads(signal_tasks_json.read_text())
                signal_to_run = [t for t in all_sig if not shard_is_complete(t, "output_phys")]

        retry_idx = 1
        while (out_root / f"submission_retry_{retry_idx}.json").exists():
            retry_idx += 1
        suffix = f"_retry_{retry_idx}"
        active_overlay_json = out_root / f"tasks{suffix}.json"
        active_signal_json = out_root / f"signal_110038_tasks{suffix}.json"
        overlay_sub = submit_dir / f"jobs_overlay{suffix}.sub"
        signal_sub = submit_dir / f"jobs_signal_110038{suffix}.sub"
        overlay_log_prefix = f"overlay{suffix}"
        signal_log_prefix = f"signal_110038{suffix}"
    else:
        all_overlay_tasks, manifest_meta = build_overlay_manifest(
            args.overlay_root,
            out_root,
            cond_tag=args.cond_tag,
            geom=args.geom,
            geometry_rotation=args.geometry_rotation,
            material_map=material_map,
            calypso_source=args.calypso_source,
            calypso_run=args.calypso_run,
            calypso_commit=args.calypso_commit,
            overlay_reco_script=frozen_reco_script,
            acts_to_faser_bin=frozen_acts_bin,
            subset=args.subset,
            nevents=args.nevents,
            seed_cov_phi=args.seed_cov_phi,
        )
        # Pre-create shard parent dirs once to prevent parallel EOS mkdir races
        for t in all_overlay_tasks:
            Path(t["output_dir"]).mkdir(parents=True, exist_ok=True)
        tasks_json.write_text(json.dumps(all_overlay_tasks, indent=2) + "\n")
        manifest_meta["run_id"] = args.run_id
        manifest_meta["created_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        manifest_meta["cond_tag"] = args.cond_tag
        manifest_meta["geom"] = args.geom
        manifest_meta["geometry_rotation"] = args.geometry_rotation
        manifest_meta["seed_cov_phi"] = args.seed_cov_phi
        manifest_meta["material_map"] = str(material_map)
        manifest_meta["calypso_source"] = str(args.calypso_source)
        manifest_meta["calypso_run"] = str(args.calypso_run)
        manifest_meta["calypso_commit"] = args.calypso_commit
        manifest_json.write_text(json.dumps(manifest_meta, indent=2) + "\n")

        overlay_to_run = all_overlay_tasks
        signal_to_run = []
        if args.include_signal_110038:
            signal_to_run = build_signal_110038_tasks(
                args.signal_rdo_dir,
                out_root,
                cond_tag=args.cond_tag,
                geom=args.geom,
                material_map=material_map,
                calypso_source=args.calypso_source,
                calypso_run=args.calypso_run,
                calypso_commit=args.calypso_commit,
                nevents=args.nevents,
            )
            for t in signal_to_run:
                Path(t["output_dir"]).mkdir(parents=True, exist_ok=True)
            signal_tasks_json.write_text(json.dumps(signal_to_run, indent=2) + "\n")

        suffix = ""
        active_overlay_json = tasks_json
        active_signal_json = signal_tasks_json
        overlay_sub = submit_dir / "jobs_overlay.sub"
        signal_sub = submit_dir / "jobs_signal_110038.sub"
        overlay_log_prefix = "overlay"
        signal_log_prefix = "signal_110038"

    if args.resubmit_missing:
        if overlay_to_run:
            active_overlay_json.write_text(json.dumps(overlay_to_run, indent=2) + "\n")
        if signal_to_run:
            active_signal_json.write_text(json.dumps(signal_to_run, indent=2) + "\n")

    schedd, top_schedds, _ = query_best_schedd(out_root, args.exclude_schedd, suffix=suffix)

    if overlay_to_run:
        write_submit_file(
            overlay_sub,
            executable=frozen_worker,
            tasks_json=active_overlay_json,
            n_tasks=len(overlay_to_run),
            logs_dir=logs_dir,
            log_prefix=overlay_log_prefix,
            stage_label="S08_Overlay2Track",
            cpus=1,
            memory_mb=args.memory_mb,
            disk_gb=args.disk_gb,
            flavour=args.flavour,
        )
    if signal_to_run:
        write_submit_file(
            signal_sub,
            executable=frozen_signal_worker,
            tasks_json=active_signal_json,
            n_tasks=len(signal_to_run),
            logs_dir=logs_dir,
            log_prefix=signal_log_prefix,
            stage_label="S08_Signal110038",
            cpus=1,
            memory_mb=args.signal_memory_mb,
            disk_gb=args.signal_disk_gb,
            flavour=args.signal_flavour,
        )

    submission_record: dict[str, Any] = {
        "time_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "run_id": args.run_id,
        "mode": "resubmit_missing" if args.resubmit_missing else ("submit" if args.submit else "dry_run"),
        "out_root": str(out_root),
        "submit_dir": str(submit_dir),
        "schedd": schedd,
        "top_schedds": top_schedds[:5],
        "excluded_schedds": list(args.exclude_schedd),
        "selection_basis": "Minimum (TotalRunningJobs + TotalIdleJobs) among ordinary 2-digit EL9 bigbirdNN schedds",
        "overlay_tasks": len(overlay_to_run),
        "overlay_tasks_json": str(active_overlay_json) if overlay_to_run else None,
        "overlay_submit_file": str(overlay_sub) if overlay_to_run else None,
        "overlay_submit_bytes": overlay_sub.stat().st_size if overlay_to_run and overlay_sub.exists() else 0,
        "signal_tasks": len(signal_to_run),
        "signal_tasks_json": str(active_signal_json) if signal_to_run else None,
        "signal_submit_file": str(signal_sub) if signal_to_run else None,
        "signal_submit_bytes": signal_sub.stat().st_size if signal_to_run and signal_sub.exists() else 0,
        "frozen_acts_to_faser_bin": str(frozen_acts_bin),
        "frozen_overlay_reco_script": str(frozen_reco_script),
        "submitted": False,
        "overlay_cluster_stdout": None,
        "signal_cluster_stdout": None,
    }

    if args.submit:
        if overlay_to_run:
            res_ov = subprocess.check_output(["condor_submit", "-name", schedd, str(overlay_sub)], text=True)
            submission_record["overlay_cluster_stdout"] = res_ov.strip()
        if signal_to_run:
            res_sg = subprocess.check_output(["condor_submit", "-name", schedd, str(signal_sub)], text=True)
            submission_record["signal_cluster_stdout"] = res_sg.strip()
        submission_record["submitted"] = True

    sub_json = out_root / f"submission{suffix}.json"
    sub_json.write_text(json.dumps(submission_record, indent=2) + "\n")
    print(json.dumps(submission_record, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
