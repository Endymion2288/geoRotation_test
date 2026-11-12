#!/usr/bin/env python3
"""Freeze every file/entry for selected production samples (default: 110051)."""
import argparse
import datetime
import json
from pathlib import Path

import uproot


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--events-per-job", type=int, default=5000)
    parser.add_argument("--dsid", action="append", help="DSID to freeze (repeatable; default: 110051)")
    args = parser.parse_args()
    selected_dsids = set(args.dsid or ["110051"])
    repo = Path(__file__).resolve().parents[1]
    stage = Path("/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_momentum")
    pilot = json.loads((repo/"configs/samples.json").read_text())
    previous = json.loads((stage/"operations_tasks.json").read_text())
    current = {(t["dsid"], t["file_index"], t["skip"]): t for t in previous if t["variant"] == "current"}
    strict = {(t["dsid"], t["file_index"], t["skip"]): t for t in previous if t["variant"] == "strict"}
    template = previous[0]
    samples = []
    tasks = []
    for sample in pilot["samples"]:
        if sample["dsid"] not in selected_dsids:
            continue
        dsid = sample["dsid"]
        files = []
        for index, path in enumerate(sorted(Path(sample["files"][0]["path"]).parent.glob("*.root"))):
            with uproot.open(path) as root_file:
                entries = root_file["CollectionTree"].num_entries
            files.append(dict(path=str(path), entries=entries, size_bytes=path.stat().st_size))
            entry = 0
            while entry < entries:
                key = (dsid, index, entry)
                old = current.get(key)
                n = old["events"] if old else min(args.events_per_job, entries-entry)
                task = dict(kind="full_sample", dsid=dsid, file_index=index, skip=entry, events=n,
                    input=str(path), runner=template["runner"], options=template["options"],
                    materials=dict(current=template["material"], strict="/eos/home-x/xcheng/FASER/geoRotation_test_results/stage02_closure/material-current-rotated.json"))
                if old:
                    assert old["input"] == str(path) and old["events"] == strict[key]["events"]
                    task["existing"] = {v: source[key]["output"] for v, source in [("current", current), ("strict", strict)]}
                shard = f"{dsid}_f{index:02d}_e{entry:06d}_{n:05d}"
                task["output"] = str(args.output/"shards"/shard)
                tasks.append(task)
                entry += n
        samples.append(dict(dsid=dsid, species=sample["species"], files=files, total_entries=sum(f["entries"] for f in files),
                            previous_discovery_events=10000))
        print(dsid, len(files), samples[-1]["total_entries"], flush=True)
    # Small complete electron samples finish early; file indices stay lexicographic within each DSID.
    order = {d: i for i, d in enumerate(["110033", "110051", "100116", "100117", "100044"])}
    tasks.sort(key=lambda t: (order[t["dsid"]], t["file_index"], t["skip"]))
    provenance = json.loads((stage/"software_controls_v3/software.json").read_text())
    manifest = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), samples=samples,
        scope="All CollectionTree entries in every ROOT file in the selected production RDO directories",
        calypso=provenance["calypso"], calypso_commit=provenance["calypso_commit"],
        library=str(stage/"software_controls_v3/lib"), software_manifest=str(stage/"software_controls_v3/software.json"),
        options=template["options"], geometry="FASERNU-04", conditions="OFLCOND-FASER-05",
        source_style_commit="33ca2281b961935456dfd4dc9fc06235e24b551e",
        samples_are_for_full_input_descriptions=True, independent_holdout_claim=False)
    args.output.mkdir(exist_ok=True, parents=True)
    (args.output/"inputs.json").write_text(json.dumps(manifest, indent=2)+"\n")
    (args.output/"tasks.json").write_text(json.dumps(tasks, indent=2)+"\n")
    print(json.dumps(dict(tasks=len(tasks), reused=sum("existing" in t for t in tasks),
        source_events=sum(t["events"] for t in tasks), new_event_processings=2*sum(t["events"] for t in tasks if "existing" not in t))))


if __name__ == "__main__":
    main()
