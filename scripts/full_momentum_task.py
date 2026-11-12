#!/usr/bin/env python3
"""Run paired Calypso reconstructions; publish histograms, never bulk states."""
import gzip
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE/"analysis" if (CODE/"analysis").exists() else CODE.parent/"analysis"))
from momentum_full_sample import consume, dump, empty, event_groups, finalize, pack_event, read_records


def packed_events(path, dsid, limit):
    for records in itertools.islice(event_groups(read_records(path)), limit):
        yield pack_event(records, dsid)


def run(task, destination):
    started = time.monotonic()
    work = Path.cwd()
    result = empty(task["dsid"])
    resources = dict(dsid=task["dsid"], file_index=task["file_index"], skip=task["skip"], events=task["events"], variants={})
    packed = work/"current_packets.jsonl.gz"
    for variant in ["current", "strict"]:
        folder = work/variant
        folder.mkdir(exist_ok=False)
        if "existing" in task:
            path = Path(task["existing"][variant])/"stages.jsonl.gz"
        else:
            command = ["python3", task["runner"], "--dsid", task["dsid"], "--input", task["input"],
                "--file-index", str(task["file_index"]), "--skip", str(task["skip"]), "--events", str(task["events"]),
                "--audit", "--compact-only", "--material", task["materials"][variant], *task["options"]]
            begin = time.monotonic()
            with (folder/"reconstruction.log").open("w") as log:
                status = subprocess.run(["/usr/bin/time", "-f", "%U %S %e %M", "-o", str(folder/"resource.txt"), *command],
                    cwd=folder, stdout=log, stderr=subprocess.STDOUT).returncode
            with (folder/"reconstruction.log").open("rb") as source, gzip.open(destination/(variant+".log.gz"), "wb") as target:
                shutil.copyfileobj(source, target)
            for name in ["run.json", "resource.txt"]:
                if (folder/name).exists():
                    shutil.copy2(folder/name, destination/(variant+"_"+name))
            if status:
                dump(destination/"resources.json", resources)
                (destination/"exit_code.txt").write_text(str(status)+"\n")
                return status
            values = (folder/"resource.txt").read_text().split()
            resources["variants"][variant] = dict(user_s=float(values[0]), system_s=float(values[1]),
                elapsed_s=float(values[2]), maxrss_kB=int(values[3]), invocation_wall_s=time.monotonic()-begin)
            path = folder/"stages.jsonl"
            resources["variants"][variant]["raw_audit_bytes"] = path.stat().st_size
        count = 0
        if variant == "current":
            with gzip.open(packed, "wt", compresslevel=1) as stream:
                for pack in packed_events(path, task["dsid"], task["events"] if "existing" in task else None):
                    stream.write(json.dumps(pack, separators=(",", ":"), allow_nan=False)+"\n")
                    count += 1
        else:
            with gzip.open(packed, "rt") as stream:
                previous = (json.loads(line) for line in stream)
                for pack in packed_events(path, task["dsid"], task["events"] if "existing" in task else None):
                    consume(result, dict(current=next(previous), strict=pack))
                    count += 1
                assert next(previous, None) is None
        assert count == task["events"], (variant, count, task["events"])
        shutil.rmtree(folder)
        print(variant, count, "events reduced", flush=True)
    packed.unlink()
    resources["worker_elapsed_s"] = time.monotonic()-started
    dump(destination/"histograms.json.gz", finalize(result))
    dump(destination/"resources.json", resources)
    dump(destination/"task.json", task)
    (destination/"exit_code.txt").write_text("0\n")
    return 0
