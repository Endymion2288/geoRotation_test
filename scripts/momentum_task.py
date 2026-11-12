#!/usr/bin/env python3
"""Execute a frozen S04 task and retain compact audit plus requested states."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tasks", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    args = parser.parse_args()
    task = json.loads(args.tasks.read_text())[args.index]
    destination = Path(task["output"])
    destination.mkdir(parents=True, exist_ok=False)
    kind=task.get('kind')
    if kind == 'full_sample':
        from full_momentum_task import run
        return run(task, destination)
    if kind=='analysis':
        command=task['command']
    elif kind=='toy':
        command=[task['binary'],str(task['replicas']),'toy.csv']
    elif kind in ['diagnostic_simulation','mode_matrix']:
        driver=task['simulation_driver' if kind=='diagnostic_simulation' else 'modes_runner']
        command=['python3',driver,'--tasks',str(args.tasks),'--index',str(args.index)]
    else:
        command = ["python3", task["runner"], "--dsid", task["dsid"],
                   "--file-index", str(task["file_index"]), "--skip", str(task["skip"]),
                   "--events", str(task["events"]), "--audit", "--compact-only",
                   "--material", task["material"], *task.get("options", [])]
        if task.get("input"):command.extend(["--input",task["input"]])
    with Path("reconstruction.log").open("w") as log:
        result = subprocess.run(["/usr/bin/time", "-v", "-o", "resource.txt", *command],
                                stdout=log, stderr=subprocess.STDOUT)
    (destination / "task.json").write_text(json.dumps(task, indent=2) + "\n")
    for name in ["run.json", "resource.txt", "toy.csv", "ntuple_check.json"]:
        if Path(name).exists():shutil.copy2(name, destination / name)
    with Path("reconstruction.log").open("rb") as source, gzip.open(destination / "reconstruction.log.gz", "wb") as output:
        shutil.copyfileobj(source, output)
    if result.returncode:
        (destination / "exit_code.txt").write_text(str(result.returncode) + "\n")
        return result.returncode
    if Path("stages.jsonl").exists():
        keep = {(c["particle"]["final"]["seed"], c["particle"]["final"]["candidate"])
                for c in task.get("cases", [])}
        counts = dict(kept=0, omitted_states=0)
        def records():
            if task.get('compact_traces'):
                from audit_reduce import reduce_event
                event=[]
                with Path('stages.jsonl').open() as source:
                    for line in source:
                        row=json.loads(line)
                        if row['stage']=='event' and event:
                            yield from reduce_event(event);event=[]
                        event.append(row)
                if event:yield from reduce_event(event)
            else:
                with Path('stages.jsonl').open() as source:
                    for line in source:yield json.loads(line)
        with gzip.open(destination / "stages.jsonl.gz", "wt") as output:
            for row in records():
                if keep and row["stage"]=="found_state" and (row["seed"], row["candidate"]) not in keep:
                    counts["omitted_states"] += 1
                    continue
                output.write(json.dumps(row,separators=(',',':'),allow_nan=False)+'\n')
                counts["kept"] += 1
        (destination / "retained_records.json").write_text(json.dumps(counts, indent=2) + "\n")
    (destination / "exit_code.txt").write_text(str(result.returncode) + "\n")
    return result.returncode


if __name__ == "__main__":raise SystemExit(main())
