#!/usr/bin/env python3
"""Assemble a native Calypso component CKF with the slope observer.

The component build contains the instrumented CKF/geometry libraries but not
the observer package.  The source prior runs used a complete overlay, so this
utility starts from that exact overlay and replaces only the Faser Acts
packages produced by the component build.  The resulting directory is passed
as GEOROTATION_LIBRARY_DIR to the real Calypso worker.
"""
import argparse
import json
from pathlib import Path
import shutil


PYTHON_PACKAGES = ("FaserActsGeometry", "FaserActsKalmanFilter")


def _copy_file(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def assemble(component, base, output, variant):
    component = component.resolve()
    base = base.resolve()
    output = output.resolve()
    if output.exists():
        raise FileExistsError(output)
    if not (component / "lib").is_dir() or not (component / "python").is_dir():
        raise FileNotFoundError(f"component build is incomplete: {component}")
    if not (base / "lib").is_dir() or not (base / "python").is_dir():
        raise FileNotFoundError(f"base overlay is incomplete: {base}")

    for folder in ("lib", "python"):
        (output / folder).mkdir(parents=True)
    # Retain every observer, clusterization and conditions package from the
    # validated source overlay.  Copying rather than symlinking leaves an
    # immutable, self-describing worker input directory.
    for folder in ("lib", "python"):
        for item in (base / folder).iterdir():
            destination = output / folder / item.name
            if item.is_dir():
                shutil.copytree(item, destination)
            else:
                shutil.copy2(item, destination)

    replaced = []
    # The geometry package does not provide a component/confdb file in every
    # build; copy exactly the files emitted by this component build rather than
    # assuming a fixed target list.
    for source in sorted((component / "lib").iterdir()):
        if not source.is_file():
            continue
        _copy_file(source, output / "lib" / source.name)
        replaced.append(str((output / "lib" / source.name).relative_to(output)))
    for package in PYTHON_PACKAGES:
        source = component / "python" / package
        destination = output / "python" / package
        if not source.is_dir():
            raise FileNotFoundError(source)
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(source, destination)
        replaced.append(str(destination.relative_to(output)))

    manifest = {
        "variant": variant,
        "component_build": str(component),
        "base_overlay": str(base),
        "actual_calypso": True,
        "isolate_smoothing": True,
        "replaced": replaced,
        "retained_base_overlay": "observer, clusterization, conditions and all other packages",
    }
    (output / "software.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--component", type=Path, required=True)
    parser.add_argument("--base-overlay", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--variant", required=True)
    args = parser.parse_args()
    print(json.dumps(assemble(args.component, args.base_overlay, args.output, args.variant), indent=2))


if __name__ == "__main__":
    main()
