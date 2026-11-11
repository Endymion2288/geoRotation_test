#!/usr/bin/env python3
"""Freeze the small worker code, tracking plugin and exact external source patch."""
import argparse
import datetime
import json
from pathlib import Path
import shutil
import subprocess


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args=parser.parse_args()
    repo=Path(__file__).resolve().parents[1]
    calypso=Path("/eos/home-x/xcheng/FASER/geoRotation_calypso")
    args.output.mkdir(parents=True, exist_ok=False)
    for directory in ["configs","scripts"]:
        shutil.copytree(repo/directory,args.output/"code"/directory,ignore=shutil.ignore_patterns("__pycache__"))
    (args.output/"lib").mkdir()
    library=calypso/"run/lib/libFaserActsKalmanFilter.so"
    shutil.copy2(library,args.output/"lib"/library.name)
    patch=git(calypso,"diff","--","Tracking/Acts/FaserActsKalmanFilter")
    for name in ["TrackingAudit.h","TrackingClosure.h"]:
        relative=Path("Tracking/Acts/FaserActsKalmanFilter/src")/name
        lines=(calypso/relative).read_text().splitlines(keepends=True)
        patch+=f"diff --git a/{relative} b/{relative}\nnew file mode 100644\n--- /dev/null\n+++ b/{relative}\n@@ -0,0 +1,{len(lines)} @@\n"
        patch+="".join("+"+line for line in lines)
    (args.output/"tracking.patch").write_text(patch)
    manifest=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        repository=str(repo),repository_commit=git(repo,"rev-parse","HEAD").strip(),
        calypso=str(calypso),calypso_commit=git(calypso,"rev-parse","HEAD").strip(),
        athena="24.0.41",platform="x86_64-el9-gcc13-opt",audit_schema=2,
        library_source=str(library),library_size=library.stat().st_size,
        library_mtime_ns=library.stat().st_mtime_ns,
        tracking_patch="tracking.patch",code="code",library="lib/libFaserActsKalmanFilter.so",
        baseline_options=dict(correct_interfaces=False,reset_event_maps=False),
        files=[dict(path=str(p.relative_to(args.output)),bytes=p.stat().st_size)
               for p in sorted(args.output.rglob("*")) if p.is_file()])
    (args.output/"software.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(f"Frozen {len(manifest['files'])} files in {args.output}")


if __name__=="__main__":main()
