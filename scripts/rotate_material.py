#!/usr/bin/env python3
"""Rotate a material map rigidly, retaining every geometry ID and material bin."""
import argparse
import copy
import json
from pathlib import Path

import numpy as np

R = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])


def transform(value):
    rotation = value.get("rotation")
    return (np.eye(3) if rotation is None else np.asarray(rotation).reshape(3, 3, order="F"),
            np.asarray(value.get("translation", [0., 0., 0.])))


def rotate(source):
    result = copy.deepcopy(source)
    for section in ("Surfaces", "Volumes"):
        for entry in result[section]["entries"]:
            utility = entry["value"].get("material", {}).get("binUtility")
            if utility is None:
                continue
            old_rotation, translation = transform(utility.get("transform", {}))
            utility["transform"] = {"rotation": (R @ old_rotation).flatten(order="F").tolist(),
                                    "translation": (R @ translation).tolist()}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=Path("/cvmfs/faser.cern.ch/repo/sw/database/DBRelease/current/acts/material-maps-alma9.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.source.read_text())
    rotated = rotate(source)
    max_error = 0.
    rng = np.random.default_rng(81723)
    points = rng.uniform(-100., 100., (4096, 3))
    points[:, 2] *= 25.
    entry_counts = {}
    for section in source:
        entry_counts[section] = len(source[section]["entries"])
        for original, output in zip(source[section]["entries"], rotated[section]["entries"]):
            assert {k: v for k, v in original.items() if k != "value"} == {k: v for k, v in output.items() if k != "value"}
            material = original["value"].get("material", {})
            if "binUtility" not in material:
                continue
            assert material["data"] == output["value"]["material"]["data"]
            u, t = transform(material["binUtility"].get("transform", {}))
            v, s = transform(output["value"]["material"]["binUtility"]["transform"])
            local_f = (points - t) @ u
            local_a = (points @ R.T - s) @ v
            max_error = max(max_error, np.max(np.abs(local_f - local_a)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rotated, separators=(",", ":")) + "\n")
    report = dict(source=str(args.source.resolve()), output=str(args.output), entries=entry_counts,
                  points_per_surface=len(points), max_local_coordinate_difference_mm=max_error,
                  geometry_ids_identical=True, material_bin_data_identical=True,
                  rotation=R.tolist())
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
