#!/usr/bin/env python3
"""Reduce actual Calypso states to mergeable reference plots and P/F/S counts."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from pilot import candidate_rank, event_groups, good_track, match, pnorm, read_records, truth_fiducial
from momentum_reference_figures import COMPONENTS, DSIDS, Histogram, NAMES, VARIANTS, layer_index

EDGES = np.linspace(-1.2, 1.2, 101)
SIGMA_EDGES = np.linspace(-5., 5., 121)
COHORTS = ("good95", "same_path_good_seed")


def dump(path, data):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as stream:
        json.dump(data, stream, separators=(",", ":"), allow_nan=False)
        stream.write("\n")


def load(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        return json.load(stream)


def finite(value):
    return value is not None and math.isfinite(value)


def qratio(bound, truth_p):
    q = bound[4] if bound else None
    return 1. / abs(q) / truth_p if finite(q) and q else None


def good_seed(bound, p_ratio, truth):
    pt = math.sqrt(sum(v*v for v in truth["momentum_F_MeV"]))
    theta = math.acos(max(-1., min(1., -truth["momentum_F_MeV"][0]/pt)))
    return bool(p_ratio is not None and abs(p_ratio-1.) < .1
                and finite(bound[4]) and bound[4]*(-truth["pdg"]) > 0
                and abs(bound[3]-theta) < 1e-4)


def pack_event(records, dsid):
    """Apply the frozen S04 selection; retain only final ancestry on scratch."""
    truth = {r["barcode"]: r for r in records if r["stage"] == "truth" and r["signal"]
             and r["status"] == 1 and abs(r["pdg"]) == (13 if dsid.startswith("100") else 11)}
    eligible = {b: t for b, t in truth.items() if truth_fiducial(t) and sum(t["hit_counts"][1:]) >= 18}
    seeds = {r["seed"]: r for r in records if r["stage"] == "seed"}
    selected = {(r["seed"], r["candidate"]): r for r in records if r["stage"] == "selected"}
    final_candidates = defaultdict(list)
    for r in records:
        if r["stage"] == "final" and match(r) in eligible and pnorm(r) is not None:
            final_candidates[match(r)].append(r)
    finals = {b: max(rs, key=candidate_rank) for b, rs in final_candidates.items()}
    finals = {b: r for b, r in finals.items() if good_track(r)}
    paths = {(r["seed"], r["candidate"]) for r in finals.values()}
    traces = {(r["seed"], r["candidate"]): r for r in records if r["stage"] == "selected_state_trace"
              and (r["seed"], r["candidate"]) in paths}
    raw = defaultdict(list)
    for r in records:
        if r["stage"] == "selected_state" and (r["seed"], r["candidate"]) in paths:
            raw[(r["seed"], r["candidate"])].append(r)
    if raw:
        from audit_reduce import trace_record
        traces.update({key: trace_record(rs) for key, rs in raw.items()})
    header = next(r for r in records if r["stage"] == "event")
    assert header["event_weight"] == 1
    result = dict(production_key=[header[k] for k in ["run", "lumi", "event"]],
                  eligible_truth=len(eligible), seed_candidates=[], particles=[])
    for r in seeds.values():
        t = eligible.get(match(r))
        if t is not None:
            pt = math.sqrt(sum(v*v for v in t["momentum_F_MeV"]))*1e-3
            p = pnorm(r)
            result["seed_candidates"].append([p/pt if p else None,
                abs(r["bound"][4])-1./pt if p else None])
    for b, r in finals.items():
        t = eligible[b]
        pt = math.sqrt(sum(v*v for v in t["momentum_F_MeV"]))*1e-3
        seed = seeds[r["seed"]]
        target = selected[(r["seed"], r["candidate"])]
        trace = traces[(r["seed"], r["candidate"])]
        seed_p = pnorm(seed)
        seed_ratio = seed_p/pt if seed_p else None
        result["particles"].append(dict(barcode=b, seed=r["seed"], candidate=r["candidate"],
            seed_bound=seed["bound"], seed_ratio=seed_ratio, good_seed=good_seed(seed["bound"], seed_ratio, t),
            selected_hit_ids=target["hits"], ckf_ratio=qratio(target["bound"], pt),
            final_ratio=pnorm(r)/pt, final_wrong=r["charge"]*(-t["pdg"]) < 0,
            refit=r["refit"], truth_p_GeV=pt, truth_charge=-1 if t["pdg"] > 0 else 1,
            steps=trace["steps"]))
    return result


def histogram(edges=EDGES):
    return Histogram(edges).result()


def add(h, value):
    if not finite(value):
        h["invalid"] += 1
        return
    h["n"] += 1
    h["sumw"] += 1
    h["sumw2"] += 1
    h["tail_half"] += value < -.5
    edges = h["bin_edges"]
    if value < edges[0]:
        h["underflow"] += 1
    elif value >= edges[-1]:
        h["overflow"] += 1
    else:
        # Uniform bins. Correct a floating-point edge rounding with one comparison.
        index = min(len(edges)-2, max(0, int((value-edges[0])/(edges[1]-edges[0]))))
        if index and value < edges[index]:
            index -= 1
        elif index+1 < len(edges)-1 and value >= edges[index+1]:
            index += 1
        h["counts"][index] += 1


def metrics():
    return dict(counts=Counter(), masks=[0]*8,
        sigma={c: histogram(SIGMA_EDGES) for c in COMPONENTS},
        qop_error={c: histogram(np.linspace(-10., 10., 201)) for c in COMPONENTS},
        max_update_error_over_truth=0., max_smoothing_error_over_truth=0.,
        max_sigma_S_over_F=0.)


def empty(dsid):
    def grid(factory):
        return [[factory() for _ in range(3)] for _ in range(3)]
    result = dict(schema=1, dsid=dsid, paired_counts=Counter(events=0), variants={})
    for v in VARIANTS:
        result["variants"][v] = dict(counts=Counter(events=0),
            histograms={k: histogram() for k in ["seed_candidates", "ckf_target", "final", "common_final",
                "same_path_good_seed_ckf", "same_path_good_seed_final"]},
            seed_abs_qop=histogram(np.linspace(-.03, .03, 101)),
            layers={cohort: grid(lambda: {c: histogram(np.linspace(-1.1, 1.1, 101) if cohort == "same_path_good_seed" else EDGES)
                                         for c in COMPONENTS}) for cohort in COHORTS},
            metrics={cohort: grid(metrics) for cohort in COHORTS},
            terminal=Counter(), terminal_max_F_S_over_truth=0.,
            material_X0=histogram(np.linspace(-5., 4., 109)))
    return result


def step_metrics(m, step, pt, charge):
    counts = m["counts"]
    counts["states"] += 1
    qs = {c: step.get(c, [None]*5)[4] for c in COMPONENTS}
    valid = all(finite(q) and q != 0 for q in qs.values())
    if valid:
        tails = {c: abs(q)*pt > 2. for c, q in qs.items()}
        m["masks"][sum(int(tails[c]) << i for i, c in enumerate(COMPONENTS))] += 1
        counts["all_three_finite"] += 1
        for a, b, name in [("predicted", "filtered", "P_F"), ("filtered", "smoothed", "F_S"),
                            ("predicted", "smoothed", "P_S")]:
            counts[name+"_repaired"] += tails[a] and not tails[b]
            counts[name+"_introduced"] += not tails[a] and tails[b]
            counts[name+"_both_tail"] += tails[a] and tails[b]
            counts[name+"_both_nontail"] += not tails[a] and not tails[b]
    else:
        counts["missing_or_invalid_parameter"] += 1
    for c, q in qs.items():
        if finite(q):
            add(m["qop_error"][c], q*pt*charge-1.)
        else:
            add(m["qop_error"][c], None)
        var = step.get(c+"_variance_qop")
        add(m["sigma"][c], math.log10(math.sqrt(var)*pt) if finite(var) and var > 0 else None)
    fvar, svar = [step.get(c+"_variance_qop") for c in ["filtered", "smoothed"]]
    if finite(fvar) and finite(svar) and fvar > 0 and svar >= 0:
        value = math.sqrt(svar/fvar)
        m["max_sigma_S_over_F"] = max(m["max_sigma_S_over_F"], value)
        counts["sigma_S_greater_F_over_1e-6"] += value > 1.+1e-6
    for field, target in [("update_identity_error", "max_update_error_over_truth"),
                          ("smoothing_identity_error", "max_smoothing_error_over_truth")]:
        value = step.get(field)
        if finite(value):
            counts[field+"_checked"] += 1
            m[target] = max(m[target], abs(value)*pt)
            counts[field+"_over_1e-6"] += abs(value)*pt > 1e-6
    counts["invalid_matrix_check"] += bool(step.get("invalid_matrix_check"))


def consume(result, packs):
    indexes = {v: {r["barcode"]: r for r in packs[v]["particles"]} for v in VARIANTS}
    common = indexes["current"].keys() & indexes["strict"].keys()
    same_path = set()
    for b in common:
        a, z = [indexes[v][b] for v in VARIANTS]
        if a["seed"] == z["seed"] and a["selected_hit_ids"] == z["selected_hit_ids"]:
            assert a["seed_bound"] == z["seed_bound"], "Paired seed means differ"
            assert a["good_seed"] == z["good_seed"]
            same_path.add(b)
    good_seeds = {b for b in same_path if indexes["current"][b]["good_seed"]}
    result["paired_counts"]["events"] += 1
    result["paired_counts"]["common_good95"] += len(common)
    result["paired_counts"]["same_seed_and_selected_hits"] += len(same_path)
    result["paired_counts"]["same_path_good_seed"] += len(good_seeds)
    assert packs["current"]["production_key"] == packs["strict"]["production_key"]
    assert packs["current"]["seed_candidates"] == packs["strict"]["seed_candidates"]
    for v, pack in packs.items():
        out = result["variants"][v]
        counts = out["counts"]
        counts["events"] += 1
        counts["fiducial_H18_truth"] += pack["eligible_truth"]
        for p, q in pack["seed_candidates"]:
            add(out["histograms"]["seed_candidates"], p-1. if p is not None else None)
            add(out["seed_abs_qop"], q)
        for row in pack["particles"]:
            b = row["barcode"]
            pt = row["truth_p_GeV"]
            counts["good95_particles"] += 1
            counts["wrong_charge"] += row["final_wrong"]
            counts["refit_fallback"] += row["refit"] == 0
            for key, value in [("ckf_target", row["ckf_ratio"]), ("final", row["final_ratio"])]:
                add(out["histograms"][key], value-1. if value is not None else None)
            if b in common:
                add(out["histograms"]["common_final"], row["final_ratio"]-1.)
            if b in good_seeds:
                counts["same_path_good_seed_particles"] += 1
                add(out["histograms"]["same_path_good_seed_ckf"], row["ckf_ratio"]-1. if row["ckf_ratio"] else None)
                add(out["histograms"]["same_path_good_seed_final"], row["final_ratio"]-1.)
            material = sum(s.get("material_full_X0", 0.)*s.get("material_path_correction", 1.) for s in row["steps"])
            add(out["material_X0"], math.log10(material) if material > 0 else None)
            measurement_steps = [s for s in row["steps"] if s["measurement_flag"]]
            for name, s in [("tip", row["steps"][-1]), ("last_measurement", measurement_steps[-1])]:
                f, z = s.get("filtered"), s.get("smoothed")
                if f and z and finite(f[4]) and finite(z[4]):
                    delta = abs(f[4]-z[4])*pt
                    out["terminal"][name+"_checked"] += 1
                    out["terminal"][name+"_exact_equal"] += f[4] == z[4]
                    out["terminal"][name+"_different_over_1e-9"] += delta > 1e-9
                    out["terminal_max_F_S_over_truth"] = max(out["terminal_max_F_S_over_truth"], delta)
            for s in measurement_steps:
                station, layer = layer_index(s["hit"])
                assert station in [1, 2, 3] and layer in [0, 1, 2]
                assert s["geometry_id"] >> 56 == {1: 25, 2: 27, 3: 29}[station]
                assert (s["geometry_id"] >> 36) & 0xfff == (layer+1)*2
                for cohort in COHORTS:
                    if cohort == "same_path_good_seed" and b not in good_seeds:
                        continue
                    cell = out["layers"][cohort][station-1][layer]
                    for c in COMPONENTS:
                        value = qratio(s.get(c), pt)
                        add(cell[c], value-1. if value is not None else None)
                    step_metrics(out["metrics"][cohort][station-1][layer], s, pt, row["truth_charge"])


def finalize(result):
    def visit(obj):
        if isinstance(obj, dict):
            if "bin_edges" in obj:
                assert obj["n"] == sum(obj["counts"])+obj["underflow"]+obj["overflow"]
                if "tail_half" in obj:
                    obj["tail_fraction"] = obj["tail_half"]/obj["n"] if obj["n"] else None
            else:
                for v in obj.values():
                    visit(v)
        elif isinstance(obj, list):
            for v in obj:
                visit(v)
    visit(result)
    assert result["variants"]["current"]["histograms"]["seed_candidates"] == result["variants"]["strict"]["histograms"]["seed_candidates"]
    for v in VARIANTS:
        data = result["variants"][v]
        assert data["counts"]["events"] == result["paired_counts"]["events"]
        # Only momentum residuals use the p/p(prod)<0.5 tail definition.
        auxiliary = [data["material_X0"], data["seed_abs_qop"]]
        for grid in data["metrics"].values():
            for station in grid:
                for m in station:
                    auxiliary.extend(h for group in ["sigma", "qop_error"] for h in m[group].values())
        for h in auxiliary:
            h.pop("tail_half", None)
            h.pop("tail_fraction", None)
    return result


def merge_values(a, b, key=""):
    if key in ["dsid", "schema", "bin_edges"]:
        assert a == b, key
        return a
    if isinstance(a, dict):
        return {k: merge_values(a.get(k, 0), b.get(k, 0), k) for k in a.keys() | b.keys()}
    if isinstance(a, list):
        assert len(a) == len(b)
        return [merge_values(x, y, key) for x, y in zip(a, b)]
    if key == "tail_fraction":
        return None
    if key.startswith("max_") or "_max_" in key:
        return max(a, b)
    return a+b


def verify_coverage(manifest, tasks):
    expected = {(s["dsid"], i): f for s in manifest["samples"] for i, f in enumerate(s["files"])}
    spans = defaultdict(list)
    for task in tasks:
        key = (task["dsid"], task["file_index"])
        assert key in expected and task["input"] == expected[key]["path"], (key, task["input"])
        spans[key].append((task["skip"], task["skip"]+task["events"]))
    assert spans.keys() == expected.keys()
    for key, f in expected.items():
        end = 0
        for start, stop in sorted(spans[key]):
            assert start == end and stop > start, (key, start, end)
            end = stop
        assert end == f["entries"], (f["path"], end)


def merge(manifest, tasks, output):
    verify_coverage(manifest, tasks)
    by_dsid = {}
    resources = []
    provenance_counts = Counter()
    for task in tasks:
        folder = Path(task["output"])
        assert (folder/"exit_code.txt").read_text().strip() == "0", folder
        part = load(folder/"histograms.json.gz")
        assert part["paired_counts"]["events"] == task["events"]
        dsid = task["dsid"]
        by_dsid[dsid] = merge_values(by_dsid[dsid], part) if dsid in by_dsid else part
        resources.append(load(folder/"resources.json"))
        for variant in VARIANTS:
            run_path = (Path(task["existing"][variant])/"run.json" if "existing" in task
                        else folder/(variant+"_run.json"))
            run = load(run_path)
            assert run["success"] and run["library_override"] == manifest["library"], run_path
            assert run["calypso_dir"] == manifest["calypso"]+"/run", run_path
            assert run["args"]["dsid"] == dsid and run["args"]["file_index"] == task["file_index"], run_path
            for key, expected in [("input", task["input"]), ("skip", task["skip"]), ("events", task["events"]),
                                  ("material", task["materials"][variant])]:
                assert run["args"][key] == expected, (run_path, key)
            for key in ["states", "selected_states_only", "reset_event_maps", "isolate_smoothing"]:
                assert run["args"][key], (run_path, key)
            provenance_counts["reused_runs" if "existing" in task else "new_runs"] += 1
    for sample in manifest["samples"]:
        assert by_dsid[sample["dsid"]]["paired_counts"]["events"] == sample["total_entries"]
    nominal_check = {}
    if "110051" in by_dsid:
        stage = Path(manifest["software_manifest"]).parent.parent
        for variant in VARIANTS:
            discovery = load(stage/f"operations_analysis/{variant}/summary.json")["samples"]["110051"]["counts"]
            independent = load(stage/f"independent_merged/{variant}/summary.json")["samples"]["110051"]["good95"]
            new = by_dsid["110051"]["variants"][variant]["histograms"]["final"]
            expected = dict(n=discovery["good95"]+independent["n"], tail_half=discovery["production_tail"]+independent["tail_half"])
            nominal_check[variant] = dict(previous_discovery_plus_independent=expected,
                new_full_sample=dict(n=new["n"], tail_half=new["tail_half"]),
                counts_identical=new["n"] == expected["n"] and new["tail_half"] == expected["tail_half"],
                scope="Aggregate final Good95 and tail counts, not an event-by-event observer equivalence claim")
    output.mkdir(parents=True, exist_ok=True)
    for d, result in by_dsid.items():
        dump(output/(d+"_histograms.json.gz"), finalize(result))
    table = []
    for dsid, result in by_dsid.items():
        for variant in VARIANTS:
            for cohort in COHORTS:
                data = result["variants"][variant]
                for station in range(3):
                    for layer in range(3):
                        h = data["layers"][cohort][station][layer]
                        m = data["metrics"][cohort][station][layer]
                        row = dict(dsid=dsid, variant=variant, cohort=cohort, station=station+1, layer=layer)
                        for c in COMPONENTS:
                            for field in ["n", "tail_half", "invalid", "underflow", "overflow"]:
                                row[c+"_"+field] = h[c][field]
                            interval = quantile_hist(m["sigma"][c])
                            row[c+"_relative_sigma_median_bin_low"] = 10**interval[0] if interval else None
                            row[c+"_relative_sigma_median_bin_high"] = 10**interval[1] if interval else None
                        for field in ["all_three_finite", "F_S_repaired", "F_S_introduced", "F_S_both_tail",
                                      "P_F_repaired", "P_F_introduced", "missing_or_invalid_parameter"]:
                            row[field] = m["counts"].get(field, 0)
                        row["max_update_error_over_truth"] = m["max_update_error_over_truth"]
                        row["max_smoothing_error_over_truth"] = m["max_smoothing_error_over_truth"]
                        table.append(row)
    with (output/"layer_metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    dump(output/"coverage.json", dict(complete=True, source_events=sum(s["total_entries"] for s in manifest["samples"]),
        event_processings=2*sum(s["total_entries"] for s in manifest["samples"]), files=sum(len(s["files"]) for s in manifest["samples"]),
        shards=len(tasks), provenance_checks=dict(provenance_counts), nominal_final_count_check=nominal_check,
        resources=resources, manifest=manifest))
    return by_dsid


def validate_discovery(root, output):
    """Compare the reused discovery shards with the frozen reference counts."""
    tasks = [t for t in load(root/"tasks.json") if "existing" in t]
    results = {}
    for task in tasks:
        folder = Path(task["output"])
        assert (folder/"exit_code.txt").read_text().strip() == "0", folder
        data = load(folder/"histograms.json.gz")
        dsid = task["dsid"]
        results[dsid] = merge_values(results[dsid], data) if dsid in results else data
    checks = {}
    output.mkdir(parents=True, exist_ok=True)
    for dsid, data in results.items():
        finalize(data)
        old = load(root.parent/"reference_figures/source_reductions"/dsid/"histograms.json")
        checks[dsid] = dict(events=data["paired_counts"]["events"])
        for variant in VARIANTS:
            a, b = data["variants"][variant], old["variants"][variant]
            assert a["counts"]["events"] == b["counts"]["events"] == 10000
            assert a["counts"]["good95_particles"] == b["counts"]["good95_particles"]
            for cohort in COHORTS:
                assert a["layers"][cohort] == b["layers"][cohort], (dsid, variant, cohort)
            assert a["seed_abs_qop"] == {k: val for k, val in b["seed"]["abs_qop_difference"].items()
                                         if k not in ["tail_half", "tail_fraction"]}
            for key in ["n", "invalid", "tail_half"]:
                assert a["histograms"]["seed_candidates"][key] == b["seed"]["relative_p"][key]
            checks[dsid][variant] = dict(good95=a["counts"]["good95_particles"], layer_bin_counts_identical=True,
                seed_abs_qop_bin_counts_identical=True)
        dump(output/(dsid+"_histograms.json.gz"), data)
    assert sorted(results) == sorted(DSIDS)
    dump(output/"checks.json", checks)


def examples(root, output):
    cases = []
    for task in load(root/"tasks.json"):
        if "existing" not in task:
            continue
        streams = [event_groups(read_records(Path(task["existing"][v])/"stages.jsonl.gz")) for v in VARIANTS]
        for ordinal, (a, b) in enumerate(zip(*streams)):
            packs = {v: pack_event(records, task["dsid"]) for v, records in zip(VARIANTS, [a, b])}
            indexes = {v: {r["barcode"]: r for r in pack["particles"]} for v, pack in packs.items()}
            for barcode in indexes["current"].keys() & indexes["strict"].keys():
                before, after = [indexes[v][barcode] for v in VARIANTS]
                if not (after["good_seed"] and before["seed_bound"] == after["seed_bound"]
                        and before["selected_hit_ids"] == after["selected_hit_ids"]):
                    continue
                pt = after["truth_p_GeV"]
                targets = [s for s in after["steps"] if s["measurement_flag"] and layer_index(s["hit"]) == (2, 1)]
                targets = [s for s in targets if (qratio(s.get("filtered"), pt) or 0) < .5
                           and (qratio(s.get("smoothed"), pt) or 0) > .8]
                if not targets:
                    continue
                case = dict(dsid=task["dsid"], file_index=task["file_index"], source_entry=task["skip"]+ordinal,
                    input=task["input"], barcode=barcode, truth_p_GeV=pt, truth_charge=after["truth_charge"],
                    seed_p_ratio=after["seed_ratio"], variants={})
                for v, row in [("current", before), ("strict", after)]:
                    steps = []
                    for step in row["steps"]:
                        if not step["measurement_flag"]:
                            continue
                        station, layer = layer_index(step["hit"])
                        item = dict(station=station, layer=layer, state=step["state"], hit=step["hit"], z_F_mm=step["z_F_mm"])
                        for c in COMPONENTS:
                            p = qratio(step.get(c), pt)
                            var = step.get(c+"_variance_qop")
                            item[c] = dict(qop_GeV_inverse=step.get(c, [None]*5)[4], p_GeV=p*pt if p else None,
                                           sigma_qop_GeV_inverse=math.sqrt(var) if finite(var) and var >= 0 else None)
                        for key in ["update_residual_mm", "update_gain_qop", "update_delta_qop", "update_identity_error", "smoothing_delta_qop"]:
                            item[key] = step.get(key)
                        steps.append(item)
                    case["variants"][v] = dict(ckf_p_ratio=row["ckf_ratio"], final_p_ratio=row["final_ratio"], steps=steps)
                cases.append(case)
                if len(cases) == 3:
                    output.mkdir(parents=True, exist_ok=True)
                    dump(output/"cases.json", cases)
                    return
    raise RuntimeError("Fewer than three qualifying same-path examples")


def quantile_hist(hist, p=.5):
    if not hist["n"]:
        return None
    cumulative = hist["underflow"]
    for i, n in enumerate(hist["counts"]):
        cumulative += n
        if cumulative >= p*hist["n"]:
            return [hist["bin_edges"][i], hist["bin_edges"][i+1]]
    return None


def plot(reductions, output):
    from momentum_reference_figures import COLORS, LABELS, draw_hist, layer_figure, source_seed_figure, stat_label
    output.mkdir(parents=True, exist_ok=True)
    coverage = load(reductions/"coverage.json")
    datasets = [load(reductions/(s["dsid"]+"_histograms.json.gz")) for s in coverage["manifest"]["samples"]]
    assert coverage["complete"]
    inventory = []
    plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42})
    def save(fig, name, pdf):
        fig.savefig(output/(name+".png"), dpi=180)
        fig.savefig(output/(name+".pdf"))
        pdf.savefig(fig)
        plt.close(fig)
        inventory.append(name)
    with PdfPages(output/"momentum_full_sample.pdf") as pdf:
        for data in datasets:
            d = data["dsid"]
            data["same_path_good_seed_particles"] = data["paired_counts"].get("same_path_good_seed", 0)
            file_count = next(len(s["files"]) for s in coverage["manifest"]["samples"] if s["dsid"] == d)
            data["scope_label"] = f"Full input: {file_count} files, {data['paired_counts']['events']:,} source events per configuration"
            reference = {"variants": {v: dict(layers=data["variants"][v]["layers"], seed={
                "relative_p": data["variants"][v]["histograms"]["seed_candidates"],
                "abs_qop_difference": data["variants"][v]["seed_abs_qop"]}) for v in VARIANTS}}
            for cohort, suffix in [("good95", "nine_layers"), ("same_path_good_seed", "good_seed_nine_layers")]:
                save(layer_figure(data, cohort, reference), d+"_"+suffix, pdf)
            seed_fig = source_seed_figure(data, reference)
            seed_fig.text(.5, .90, data["scope_label"], ha="center", fontsize=10)
            save(seed_fig, d+"_seed_views", pdf)
            fig, axes = plt.subplots(1, 3, figsize=(16, 5.8))
            for ax, key, title in zip(axes, ["seed_candidates", "ckf_target", "final"],
                                      ["All matched seed candidates", "Selected CKF target", "Final refit / fallback"]):
                ymax = 0
                for v in VARIANTS:
                    h = data["variants"][v]["histograms"][key]
                    ymax = max(ymax, draw_hist(ax, h, COLORS[v], LABELS[v]+"\n"+stat_label(h),
                                              "--" if v == "strict" else "-", raw=True))
                ax.set(xlabel=r"$(p-p_{\mathrm{prod}})/p_{\mathrm{prod}}$", ylabel="Entries / bin", title=title,
                       xlim=(-1.2, 1.2), ylim=(0, 1.4*ymax))
                ax.axvline(-.5, color=".4", ls="--", lw=.8)
                ax.legend(fontsize=8)
                ax.grid(axis="y", alpha=.18)
            fig.suptitle(NAMES[d]+" | "+data["scope_label"], y=.98)
            fig.text(.5, .035, "Rotated detector geometry in both runs; S04 fiducial+H18 / SDO matching / final Good95; production truth.\n"
                     "Seed candidates use truth acceptance and matching; CKF/final use final-track ancestry. Raw counts; overflow retained in statistics.", ha="center")
            fig.subplots_adjust(top=.84, bottom=.22, left=.06, right=.99, wspace=.28)
            save(fig, d+"_seed_ckf_final", pdf)
            fig, axes = plt.subplots(2, 2, figsize=(13, 9))
            for j, v in enumerate(VARIANTS):
                cells = [c for s in data["variants"][v]["layers"]["good95"] for c in s]
                x = np.arange(9)
                for c, color in COMPONENTS.items():
                    axes[0, j].plot(x, [100*cell[c]["tail_fraction"] if cell[c]["n"] else np.nan for cell in cells],
                                    "o-", color=color, label=c.title())
                ms = [c for s in data["variants"][v]["metrics"]["good95"] for c in s]
                n = np.array([m["counts"].get("all_three_finite", 0) for m in ms], dtype=float)
                for field, label, color in [("F_S_repaired", "Filtered tail corrected by smoothing", "#244de1"),
                                             ("F_S_introduced", "Tail introduced by smoothing", "#d62728")]:
                    axes[1, j].plot(x, [100*m["counts"].get(field, 0)/den if den else np.nan for m, den in zip(ms, n)],
                                    "o-", label=label, color=color)
                for ax in axes[:, j]:
                    ax.set_xticks(x)
                    ax.set_xticklabels([f"S{s}L{l}" for s in [1, 2, 3] for l in [0, 1, 2]])
                    ax.set_ylabel("Measurement states [%]")
                    ax.grid(alpha=.2)
                    ax.legend(fontsize=8)
                axes[0, j].set_title(LABELS[v])
                axes[1, j].legend(loc="lower center", bbox_to_anchor=(.5, 1.01), fontsize=8)
            for row in axes:
                ymax = max(ax.get_ylim()[1] for ax in row)
                for ax in row:
                    ax.set_ylim(0, ymax)
            fig.suptitle(NAMES[d]+" | Forward estimates and backward smoothing", y=.98)
            fig.text(.5, .035, data["scope_label"]+"; same-state comparisons; tail: p / p(prod) < 0.5.\n"
                     "Final Good95 ancestry; two strip sides pooled. Lower row uses states with all P/F/S finite.", ha="center")
            fig.subplots_adjust(top=.9, bottom=.13, hspace=.3, wspace=.25)
            save(fig, d+"_pfs_tail_transitions", pdf)
            fig, axes = plt.subplots(2, 2, figsize=(13, 9))
            for j, v in enumerate(VARIANTS):
                ms = [cell for station in data["variants"][v]["metrics"]["good95"] for cell in station]
                for c, color in COMPONENTS.items():
                    medians = [quantile_hist(m["sigma"][c]) for m in ms]
                    axes[0, j].plot(np.arange(9), [10**(.5*sum(q)) if q else np.nan for q in medians],
                                    "o-", color=color, label=c.title())
                    draw_hist(axes[1, j], ms[4]["qop_error"][c], color, c.title())
                axes[0, j].set(yscale="log", ylabel=r"Median $\sigma_{q/p}/|q/p_{\mathrm{prod}}|$", title=LABELS[v])
                axes[0, j].set_xticks(np.arange(9))
                axes[0, j].set_xticklabels([f"S{s}L{l}" for s in [1, 2, 3] for l in [0, 1, 2]])
                axes[1, j].set(xlabel=r"$(q/p)/(q/p)_{\mathrm{prod}}-1$", ylabel="Measurement states [% / bin]",
                               title="Station 2, layer 1: signed curvature error", xlim=(-10, 10))
                for limit in [-3., 1.]:
                    axes[1, j].axvline(limit, color=".4", ls="--", lw=.8)
                for ax in axes[:, j]:
                    ax.grid(alpha=.2)
                    ax.legend(fontsize=8)
            ymin = min(ax.get_ylim()[0] for ax in axes[0])
            ymax = max(ax.get_ylim()[1] for ax in axes[0])
            for ax in axes[0]:
                ax.set_ylim(ymin, ymax)
            ymax = max(ax.get_ylim()[1] for ax in axes[1])
            for ax in axes[1]:
                ax.set_ylim(0, ymax)
            fig.suptitle(NAMES[d]+" | Curvature uncertainty and interim momentum tails", y=.98)
            fig.text(.5, .035, data["scope_label"]+"; final Good95 ancestry. Medians estimated from fixed log-sigma bins.\n"
                     "For the lower axis, momentum tail p/p(prod)<0.5 lies below -3 or above +1; overflow remains in normalization.", ha="center")
            fig.subplots_adjust(top=.9, bottom=.13, hspace=.3, wspace=.28)
            save(fig, d+"_curvature_information", pdf)
        example_file = reductions.parent/"examples/cases.json"
        if example_file.exists():
            cases = load(example_file)
            fig, axes = plt.subplots(len(cases), 2, figsize=(14, 3.4*len(cases)), squeeze=False)
            for i, case in enumerate(cases):
                limits = []
                for j, v in enumerate(VARIANTS):
                    ax = axes[i, j]
                    steps = case["variants"][v]["steps"]
                    for c, color in COMPONENTS.items():
                        values = [s[c]["p_GeV"]/case["truth_p_GeV"] if s[c]["p_GeV"] else np.nan for s in steps]
                        limits.extend(x for x in values if math.isfinite(x) and x > 0)
                        ax.plot(np.arange(1, len(steps)+1), values, "o-", ms=3, color=color, label=c.title())
                    for k in range(1, len(steps)):
                        if steps[k]["station"] != steps[k-1]["station"]:
                            ax.axvline(k+.5, color=".6", ls=":", lw=1.)
                    ax.plot(0, case["seed_p_ratio"], "ks", label="Seed")
                    ax.axhline(1., color=".4", lw=.8)
                    ax.axhline(.5, color=".4", lw=.8, ls="--")
                    ax.set(yscale="log", xlabel="Forward measurement index (station boundaries dotted)",
                           ylabel=r"$p/p_{\mathrm{prod}}$", title=f"{LABELS[v]} | file {case['file_index']}, entry {case['source_entry']}, barcode {case['barcode']}")
                    ax.grid(axis="y", alpha=.18)
                for ax in axes[i]:
                    ax.set_ylim(min(limits)*.7, max(limits)*1.4)
            fig.suptitle("110051: same seed and selected hits | Examples of forward tails corrected by smoothing", y=.99)
            fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(.5, .965),
                       ncol=4, frameon=False, fontsize=10)
            fig.text(.5, .015, "Three selected examples from file 0 of the 20-file sample; production truth; both strip sides shown. These cases do not estimate population rates.", ha="center", fontsize=9)
            fig.tight_layout(rect=(0, .04, 1, .925))
            save(fig, "110051_selected_state_evolution", pdf)
    summary = dict(complete=True, source_events=coverage["source_events"], event_processings=coverage["event_processings"],
        figures=inventory, pages=len(inventory), samples={d["dsid"]: dict(events=d["paired_counts"]["events"],
            paired_counts=d["paired_counts"], variants={v: dict(counts=d["variants"][v]["counts"],
                final=d["variants"][v]["histograms"]["final"], terminal=d["variants"][v]["terminal"],
                terminal_max_F_S_over_truth=d["variants"][v]["terminal_max_F_S_over_truth"]) for v in VARIANTS}) for d in datasets})
    dump(output/"summary.json", summary)
    lines = ["# Full-input Calypso momentum figures", "", "[Complete PDF](momentum_full_sample.pdf)", "",
             f"{coverage['source_events']:,} source events; {coverage['event_processings']:,} event processings; both material configurations.", "",
             "Source style: zsm2024/FASER_tracking/test_pResolution, commit 33ca2281b961935456dfd4dc9fc06235e24b551e. S04 selection and production truth.", ""]
    lines += ["[Layer counts and P/F/S transitions (CSV)](../reductions/layer_metrics.csv)", ""]
    lines += [f"- [{stem} PNG]({stem}.png) / [PDF]({stem}.pdf)" for stem in inventory]
    (output/"figure_index.md").write_text("\n".join(lines)+"\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["merge", "plot", "validate-discovery", "examples"], required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.phase == "merge":
        merge(load(args.root/"inputs.json"), load(args.root/"tasks.json"), args.output)
    elif args.phase == "plot":
        plot(args.root, args.output)
    elif args.phase == "validate-discovery":
        validate_discovery(args.root, args.output)
    else:
        examples(args.root, args.output)


if __name__ == "__main__":
    main()
