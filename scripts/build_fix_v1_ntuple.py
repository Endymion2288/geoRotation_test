#!/usr/bin/env python3
"""Compile PhysicsAnalysis/NtupleDumper from geoRotation_calypso_latest_fix_v1
against geoRotation_calypso_latest/run installed libraries and headers.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1"),
    )
    parser.add_argument(
        "--latest",
        type=Path,
        default=Path("/eos/home-x/xcheng/FASER/geoRotation_calypso_latest"),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()

    old = Path("/eos/home-x/xcheng/FASER/geoRotation_calypso")
    latest = args.latest.resolve()
    source = args.source.resolve()
    output = args.output.resolve()
    package = Path("PhysicsAnalysis/NtupleDumper")

    recipe_file = old / "build/compile_commands.json"
    records = [
        r
        for r in json.loads(recipe_file.read_text())
        if Path(r["directory"]).name == "NtupleDumper" and "NtupleDumper.dir/" in r["output"]
    ]
    assert len(records) == 2, f"Expected 2 compile units for NtupleDumper, found {len(records)}"

    if output.exists():
        shutil.rmtree(output)
    scratch = output / "_build_tmp"
    lib = output / "lib"
    bin_dir = output / "bin"
    scratch.mkdir(parents=True, exist_ok=True)
    lib.mkdir(parents=True, exist_ok=True)
    bin_dir.mkdir(parents=True, exist_ok=True)

    commands = []
    objects = []
    for i, rec in enumerate(records):
        src_rel = Path(rec["file"]).relative_to(old / package)
        src_file = source / package / src_rel
        assert src_file.is_file(), f"Missing source file: {src_file}"
        obj = scratch / f"unit_{i}_{src_file.stem}.o"
        objects.append(str(obj))
        tokens = shlex.split(rec["command"])
        cleaned = []
        skip_next = False
        for t in tokens:
            if skip_next:
                skip_next = False
                continue
            if t == "-g":
                continue
            if "libchecker_gccplugins.so" in t:
                continue
            cleaned.append(t.replace(str(old), str(latest)))
        cleaned[cleaned.index("-o") + 1] = str(obj)
        cleaned[-1] = str(src_file)
        cleaned.insert(1, "-I" + str(latest / "run/include"))
        cleaned.insert(1, "-I" + str(source / package / "src"))
        cleaned.insert(1, "-I" + str(source / package))
        commands.append(cleaned)

    (output / "compile_recipes.json").write_text(json.dumps(commands, indent=2) + "\n")

    def compile_one(cmd: list[str]) -> None:
        print("Compiling", cmd[-1], flush=True)
        subprocess.run(cmd, check=True)

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(compile_one, commands))

    original = shlex.split((old / "build" / package / "CMakeFiles/NtupleDumper.dir/link.txt").read_text())
    if original[0].endswith("atlas_build_run.sh"):
        original = original[1:]

    linked = []
    skip = False
    for token in original:
        if skip:
            skip = False
            continue
        if token == "-o":
            skip = True
            continue
        if token.endswith(".o") or "libchecker_gccplugins.so" in token:
            continue
        if not token.startswith("-") and token.endswith((".so", ".a")) and not Path(token).is_absolute():
            token = str((old / "build" / package / token).resolve())
        if token.startswith(str(old)) and token.endswith((".so", ".a")):
            cand = latest / "run/lib" / Path(token).name
            assert cand.is_file(), f"Missing library in latest/run/lib: {cand}"
            token = str(cand)
        token = token.replace(str(old / "build"), str(latest / "run"))
        linked.append(token)

    linked += ["-o", str(lib / "libNtupleDumper.so")]
    linked[1:1] = objects
    linked += ["-L" + str(latest / "run/lib"), "-lTrackingTruthUtilsLib"]
    (output / "link_recipe.json").write_text(json.dumps(linked, indent=2) + "\n")
    print("Linking libNtupleDumper.so...", flush=True)
    subprocess.run(linked, check=True)

    py_dir = output / "python/NtupleDumper"
    shutil.copytree(source / package / "python", py_dir, dirs_exist_ok=True)
    (py_dir / "__init__.py").touch()

    for script in (source / package / "scripts").glob("*.py"):
        dst = bin_dir / script.name
        shutil.copy2(script, dst)
        dst.chmod(0o755)

    tools = Path(
        "/cvmfs/atlas.cern.ch/repo/sw/software/24.0/AthenaExternals/24.0.41/InstallArea/x86_64-el9-gcc13-opt/bin"
    )
    env = dict(
        os.environ,
        LD_LIBRARY_PATH=f"{lib}:{latest / 'run/lib'}:{os.environ.get('LD_LIBRARY_PATH', '')}",
    )
    gen = scratch / "genconf"
    gen.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(tools / "listcomponents"), "--output", str(lib / "libNtupleDumper.components"), "libNtupleDumper.so"],
        env=env,
        check=True,
    )
    subprocess.run(
        [str(tools / "genconf"), "-o", str(gen), "-p", "NtupleDumper", "--no-init", "-i", "libNtupleDumper.so"],
        env=env,
        check=True,
    )
    shutil.copy2(gen / "NtupleDumper.confdb", lib / "NtupleDumper.confdb")
    inputs = gen / "inputs.txt"
    inputs.write_text(str(gen / "NtupleDumper.confdb2_part") + "\n")
    merge = Path(
        "/cvmfs/atlas.cern.ch/repo/sw/software/24.0/Athena/24.0.41/InstallArea/x86_64-el9-gcc13-opt/cmake/modules/scripts/mergeConfdb2.py"
    )
    subprocess.run([str(merge), str(lib / "NtupleDumper.confdb2"), str(inputs)], env=env, check=True)
    shutil.copy2(gen / "NtupleDumperConf.py", py_dir / "NtupleDumperConf.py")

    manifest = {
        "source": str(source),
        "latest": str(latest),
        "package": str(package),
        "commit": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "dirty_diff": subprocess.check_output(["git", "-C", str(source), "diff", "--", str(package)], text=True),
        "shared_dependencies": str(latest / "run"),
    }
    (output / "software.json").write_text(json.dumps(manifest, indent=2) + "\n")
    shutil.rmtree(scratch, ignore_errors=True)
    print("Built NtupleDumper overlay cleanly at:", output)


if __name__ == "__main__":
    main()
