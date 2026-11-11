#!/usr/bin/env python3
"""Read POOL metadata and schema once per production; run in Calypso."""
import argparse
import json
from pathlib import Path

import uproot
from PyUtils.MetaReader import read_metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = {"schema_version": 1, "samples": []}
    for sample in json.loads(args.samples.read_text())["samples"]:
        path = sample["files"][0]["path"]
        metadata = read_metadata([path], mode="full")
        with uproot.open(path) as root_file:
            tree = root_file["CollectionTree"]
            schema = tree.typenames(recursive=False)
        row = {**sample, "metadata": metadata, "collection_schema": schema}
        result["samples"].append(row)
        print(sample["dsid"], "entries", sample["total_entries"], "branches", len(schema), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, default=str) + "\n")


if __name__ == "__main__":
    main()
