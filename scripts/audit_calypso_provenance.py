#!/usr/bin/env python3
"""Check the Calypso provenance of existing S04 reconstruction jobs."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path


CALYPSO = Path("/eos/home-x/xcheng/FASER/geoRotation_calypso")
BASE_COMMIT = "88c457c1e4bf4154e67b96114ee8e7020d9b1a97"


def inspect_job(task):
    folder = Path(task["output"])
    assert task.get("kind") is None, (folder, "Expected a Calypso reconstruction task")
    assert (folder / "exit_code.txt").read_text().strip() == "0", folder
    run = json.loads((folder / "run.json").read_text())
    assert run["success"] and run["calypso_dir"] == str(CALYPSO / "run"), folder
    assert run["athena"] == "24.0.41", folder
    assert run["input"] == task["input"], folder
    for field in ["dsid", "file_index", "skip", "events"]:
        assert run["args"][field] == task[field], (folder, field)
    assert task["material"] == run["material_source"] == run["args"]["material"], folder
    for field in ["correct_interfaces", "correct_covariance", "correct_refit_surface",
                  "no_energy_loss", "no_scattering", "bypass_refit", "fixed_truth_hits",
                  "electron_hypothesis"]:
        assert not run["args"][field], (folder, field)
    assert run["args"]["seed_override"] is None, folder
    library = Path(run["library_override"]) / "libFaserActsKalmanFilter.so"
    assert library.is_file(), library
    with gzip.open(folder / "reconstruction.log.gz", "rt", errors="replace") as stream:
        log = stream.read()
    material_line = "Configured to use material input: " + task["material"]
    assert material_line in log, folder
    assert "FaserActsTrackingGeometrySvc is initializing" in log, folder
    assert "Acts version is: v32.0.2" in log, folder
    assert "CKF_Audit" in log and "CombinatorialKalmanFilterAlg::finalize()" in log, folder
    assert f"{task['events']} events processed." in log, folder
    assert "Application Manager Terminated successfully" in log, folder
    return dict(dsid=task["dsid"], variant=task["variant"], file_index=task["file_index"],
                skip=task["skip"], events=task["events"], input=task["input"], output=str(folder),
                runner=task["runner"], material=task["material"], library=str(library),
                options=task.get("options", []), hostname=run["hostname"],
                checks="Task, run manifest, Calypso initialization/finalization and material log agree")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(), calypso=str(CALYPSO),
                  calypso_base_commit=BASE_COMMIT, athena="24.0.41", acts="32.0.2",
                  standalone_reconstruction_used=False, historical_process_maps_recorded=False,
                  library_evidence="run.json library_override, frozen software.json and worker LD_LIBRARY_PATH; historical /proc maps were not saved",
                  groups={}, software={})
    for group in ["independent", "operations"]:
        task_file = args.root / (group + "_tasks.json")
        tasks = json.loads(task_file.read_text())
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            jobs = list(pool.map(inspect_job, tasks))
        pairs = {}
        for task in tasks:
            key = (task["dsid"], task["file_index"], task["skip"], task["events"])
            assert task["variant"] not in pairs.setdefault(key, {}), (group, key)
            pairs[key][task["variant"]] = task
        for pair in pairs.values():
            assert set(pair) == {"current", "strict"}
            left, right = pair["current"], pair["strict"]
            assert {k: v for k, v in left.items() if k not in ["variant", "material", "output"]} == {
                k: v for k, v in right.items() if k not in ["variant", "material", "output"]}
        events = Counter()
        for job in jobs:
            events[job["variant"]] += job["events"]
            library = Path(job["library"])
            if str(library) in result["software"]:
                continue
            software_file = library.parent.parent / "software.json"
            software = json.loads(software_file.read_text())
            assert software["calypso"] == str(CALYPSO)
            assert software["calypso_commit"] == BASE_COMMIT
            assert software["library_source"] == str(CALYPSO / "run/lib/libFaserActsKalmanFilter.so")
            assert library.stat().st_size == software["library_size"]
            patch = software_file.parent / software["tracking_patch"]
            assert patch.is_file()
            result["software"][str(library)] = dict(manifest=str(software_file),
                tracking_patch=str(patch), bytes=library.stat().st_size,
                original_library_source=software["library_source"], created_utc=software["created_utc"])
        result["groups"][group] = dict(task_file=str(task_file), jobs=len(jobs),
            events_by_variant=dict(events), paired_ranges=len(pairs), all_checks_passed=True, records=jobs)
    result["all_checks_passed"] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"groups": {g: {k: v for k, v in values.items() if k != "records"}
                                for g, values in result["groups"].items()},
                      "software": result["software"], "all_checks_passed": True}, indent=2))


if __name__ == "__main__":
    main()
