#!/usr/bin/env python3
"""Combine gzip members without decoding the large audit; preserve replay offsets."""
import argparse
import json
from pathlib import Path
import shutil
import tarfile


def main():
    parser=argparse.ArgumentParser()
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input",type=Path)
    source.add_argument("--restore-from",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--prune-shards",action="store_true")
    args=parser.parse_args()
    if args.restore_from:
        args.output.mkdir(parents=True,exist_ok=False)
        index=json.loads((args.restore_from/"index.json").read_text())
        with tarfile.open(args.restore_from/"job_records.tar.gz","r:gz") as source:
            source.extractall(args.output)
        for task in index["tasks"]:
            output=args.output/task["shard"]/"stages.jsonl.gz"
            output.parent.mkdir(parents=True,exist_ok=True)
            with (args.restore_from/task["audit_archive"]).open("rb") as source, output.open("wb") as destination:
                source.seek(task["compressed_offset"])
                remaining=task["compressed_bytes"]
                while remaining:
                    data=source.read(min(remaining,8*1024*1024))
                    if not data:raise EOFError(task["audit_archive"])
                    destination.write(data);remaining-=len(data)
        print(f"Restored {len(index['tasks'])} shards in {args.output}")
        return
    tasks=json.loads((args.input/"tasks.json").read_text())
    assert all((args.input/t["shard"]/"exit_code.txt").read_text().strip()=="0" for t in tasks)
    args.output.mkdir(parents=True,exist_ok=False)
    index={"input":str(args.input),"audit_format":"Concatenated gzip members; indexed offsets recover each original shard",
           "tasks":[],"archives":{}}
    source_files=[]
    for dsid in sorted({t["dsid"] for t in tasks}):
        archive=args.output/f"{dsid}_stages.jsonl.gz"
        with archive.open("wb") as output:
            for task in (t for t in tasks if t["dsid"]==dsid):
                folder=args.input/task["shard"]
                source=folder/"stages.jsonl.gz"
                offset=output.tell()
                with source.open("rb") as stream:shutil.copyfileobj(stream,output,8*1024*1024)
                entry=dict(task,audit_archive=archive.name,compressed_offset=offset,
                    compressed_bytes=source.stat().st_size,run=json.loads((folder/"run.json").read_text()),
                    resource=(folder/"resource.txt").read_text())
                index["tasks"].append(entry)
                source_files.append(source)
        expected=sum(t["compressed_bytes"] for t in index["tasks"] if t["dsid"]==dsid)
        assert archive.stat().st_size==expected
        index["archives"][archive.name]=dict(bytes=expected)
    # Small manifests and logs stay available as one archive, including Condor
    # event history and submission provenance. Large audit members are above.
    with tarfile.open(args.output/"job_records.tar.gz","w:gz") as output:
        for source in sorted(args.input.rglob("*")):
            if source.is_file() and source.name!="stages.jsonl.gz":
                output.add(source,arcname=str(source.relative_to(args.input)),recursive=False)
                source_files.append(source)
    (args.output/"index.json").write_text(json.dumps(index,indent=2)+"\n")
    if args.prune_shards:
        for source in source_files:
            if source.parent!=args.input:source.unlink()
        for task in tasks:(args.input/task["shard"]).rmdir()
        # Per-shard Condor stdout/stderr are also in job_records.tar.gz.
        for source in args.input.iterdir():
            if source.suffix in (".out",".err"):source.unlink()
    print(json.dumps({"tasks":len(tasks),"audit_bytes":sum(a["bytes"] for a in index["archives"].values()),
                      "output":str(args.output),"pruned":args.prune_shards},indent=2))


if __name__=="__main__":main()
