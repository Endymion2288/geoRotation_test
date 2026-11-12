#!/usr/bin/env python3
"""Compare ACTS material JSON payloads and binning transforms; no track fitting."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np


def identity(entry):
    return tuple(sorted((key, value) for key, value in entry.items() if key != "value"))


def describe(value):
    material = value.get("material", {})
    grid = material.get("data", [])
    cells = [cell for row in grid for cell in row]
    utility = material.get("binUtility", {})
    return dict(bins=len(cells), rows=len(grid), columns=sorted({len(row) for row in grid}),
                vacuum_bins=sum(cell["material"] is None and cell["thickness"] == 0 for cell in cells),
                transform=utility.get("transform"), binningdata=utility.get("binningdata"),
                metadata={key: val for key, val in material.items() if key not in ["data", "binUtility"]})


def transform(utility):
    data = utility.get("transform", {})
    rotation = data.get("rotation")
    matrix = np.eye(3) if rotation is None else np.asarray(rotation).reshape(3, 3, order="F")
    return matrix, np.asarray(data.get("translation", [0., 0., 0.]))


def compare_values(first, second):
    a, b = first["material"], second["material"]
    u, v = a["binUtility"], b["binUtility"]
    ra, ta = transform(u)
    rb, tb = transform(v)
    result = dict(first=describe(first), second=describe(second),
        data_equal=a["data"] == b["data"],
        material_metadata_equal={k: x for k, x in a.items() if k not in ["data", "binUtility"]} == {
            k: x for k, x in b.items() if k not in ["data", "binUtility"]},
        binning_equal={k: x for k, x in u.items() if k != "transform"} == {
            k: x for k, x in v.items() if k != "transform"},
        rotation_equal=bool(np.array_equal(ra, rb)),
        translation_first_minus_second_mm=(ta-tb).tolist())
    first_normalized = dict(first)
    first_normalized["material"] = dict(a)
    first_normalized["material"]["binUtility"] = dict(u)
    first_normalized["material"]["binUtility"]["transform"] = dict(u["transform"], translation=tb.tolist())
    result["equal_after_replacing_translation"] = first_normalized == second
    if result["rotation_equal"]:
        offset = ra.T @ (tb-ta)
        result["local_lookup_first_minus_second_mm"] = offset.tolist()
        axes = [item["value"] for item in u["binningdata"]]
        result["only_local_xy_binning"] = all(axis in ["binX", "binY"] for axis in axes)
        result["lookup_xy_equal_for_any_global_point"] = bool(
            result["only_local_xy_binning"] and np.array_equal(offset[:2], [0., 0.]))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=Path, required=True)
    parser.add_argument("--second", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sources = []
    maps = []
    for path in [args.first, args.second]:
        stat = path.stat()
        sources.append(dict(path=str(path.resolve()), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns))
        maps.append(json.loads(path.read_text()))
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(), sources=sources,
                  operation="Read-only material configuration comparison; no reconstruction", sections={})
    for section in ["Surfaces", "Volumes"]:
        a, b = [data[section] for data in maps]
        left, right = [{identity(entry): entry["value"] for entry in data["entries"]} for data in [a, b]]
        assert len(left) == len(a["entries"]) and len(right) == len(b["entries"])
        common = [dict(identity=dict(key), **compare_values(left[key], right[key]))
                  for key in sorted(left.keys() & right.keys())]
        only_first = [dict(identity=dict(key), **describe(left[key])) for key in sorted(left.keys()-right.keys())]
        only_second = [dict(identity=dict(key), **describe(right[key])) for key in sorted(right.keys()-left.keys())]
        result["sections"][section] = dict(first_entries=len(left), second_entries=len(right),
            headers_equal={k: x for k, x in a.items() if k != "entries"} == {
                k: x for k, x in b.items() if k != "entries"},
            common=common, only_first=only_first, only_second=only_second,
            common_bins=sum(row["first"]["bins"] for row in common),
            common_data_equal=all(row["data_equal"] for row in common),
            common_equal_except_translation=all(row["equal_after_replacing_translation"] for row in common),
            common_lookup_xy_equal=all(row.get("lookup_xy_equal_for_any_global_point", False) for row in common),
            added_bins=sum(row["bins"] for row in only_second),
            added_nonvacuum_bins=sum(row["bins"]-row["vacuum_bins"] for row in only_second))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({section: {key: value for key, value in values.items()
                               if key not in ["common", "only_first", "only_second"]}
                      for section, values in result["sections"].items()}, indent=2))


if __name__ == "__main__":
    main()
