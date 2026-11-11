#!/usr/bin/env python3
"""Prepare immutable event ranges and small Condor descriptors; does not submit."""
import argparse
import datetime
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--submit-dir", type=Path, required=True)
    parser.add_argument("--events-per-sample", type=int, default=10000)
    parser.add_argument("--shard-events", type=int, default=1000)
    parser.add_argument("--mode", choices=["current", "historical", "consistent"], default="current")
    parser.add_argument("--dsids",nargs="+",help="Restrict samples for worker checks or controls")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    args.output.mkdir(parents=True, exist_ok=True)
    args.submit_dir.mkdir(parents=True, exist_ok=True)
    tasks = []
    for sample in json.loads((repo / "configs/samples.json").read_text())["samples"]:
        if args.dsids and sample["dsid"] not in args.dsids:continue
        remaining = min(args.events_per_sample, sample["total_entries"])
        for file_index, source in enumerate(sample["files"]):
            used = min(remaining, source["entries"])
            for skip in range(0, used, args.shard_events):
                events = min(args.shard_events, used - skip)
                shard = f"{sample['dsid']}_f{file_index:02d}_e{skip:06d}_{events:05d}_{args.mode}"
                tasks.append(dict(dsid=sample["dsid"], file_index=file_index, skip=skip, events=events, mode=args.mode, shard=shard))
            remaining -= used
            if not remaining:
                break
    (args.output / "tasks.json").write_text(json.dumps(tasks, indent=2) + "\n")
    task_file = args.submit_dir / "tasks.txt"
    task_file.write_text("".join(f"{t['dsid']} {t['file_index']} {t['skip']} {t['events']} {t['mode']} {t['shard']}\n" for t in tasks))
    submit = f"""universe = vanilla
executable = {repo}/scripts/condor_worker.sh
arguments = {repo} {args.output} $(dsid) $(file_index) $(skip) $(events) $(mode) $(shard)
request_cpus = 1
request_memory = 3500MB
request_disk = 4GB
+JobFlavour = \"microcentury\"
+WantOS = \"el9\"
+GeoRotationStage = \"S03\"
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = \"\"
getenv = False
output = {args.output}/$(shard).out
error = {args.output}/$(shard).err
log = {args.output}/condor.events.log
queue dsid, file_index, skip, events, mode, shard from {task_file}
"""
    (args.submit_dir / "jobs.sub").write_text(submit)
    os.chmod(repo / "scripts/condor_worker.sh", 0o755)
    print(f"Prepared {len(tasks)} jobs, {sum(t['events'] for t in tasks)} events")


if __name__ == "__main__":
    main()
