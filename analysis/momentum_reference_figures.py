#!/usr/bin/env python3
"""Reproduce the seed and nine-layer residual views from problem.pdf page 3."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np

from pilot import event_groups, match, pnorm, read_records, truth_fiducial


DSIDS = ["100044", "100116", "100117", "110033", "110051"]
NAMES = {"100044": r"100044: $\mu^+$ (5 mrad FLUKA E)",
         "100116": r"100116: $\mu^+$ (2D FLUKA)",
         "100117": r"100117: $\mu^-$ (2D FLUKA)",
         "110033": r"110033: $A'\to e^+e^-$, $\epsilon=10^{-5}$",
         "110051": r"110051: $A'\to e^+e^-$, $\epsilon=10^{-4}$"}
VARIANTS = ["current", "strict"]
LABELS = {"current": "Before material-map rotation", "strict": "After material-map rotation"}
COLORS = {"current": "#c43c39", "strict": "#0072b2"}
COMPONENTS = {"predicted": "#2a9d55", "filtered": "#d62728", "smoothed": "#244de1"}
EDGES = np.linspace(-2., 2., 101)
XLABEL = r"$(p-p_{\mathrm{prod}})/p_{\mathrm{prod}}$"


class Histogram:
    def __init__(self, edges=EDGES):
        self.edges = edges
        self.counts = np.zeros(len(edges)-1, dtype=np.int64)
        self.n = self.invalid = self.underflow = self.overflow = self.tail = 0

    def add(self, value):
        if value is None or not math.isfinite(value):
            self.invalid += 1
            return
        self.n += 1
        self.tail += value < -.5
        if value < self.edges[0]:
            self.underflow += 1
        elif value >= self.edges[-1]:
            self.overflow += 1
        else:
            self.counts[int(np.searchsorted(self.edges, value, side="right"))-1] += 1

    def result(self):
        assert self.n == int(self.counts.sum()) + self.underflow + self.overflow
        return dict(counts=self.counts.tolist(), bin_edges=self.edges.tolist(), n=self.n, invalid=self.invalid,
                    underflow=self.underflow, overflow=self.overflow, tail_half=self.tail,
                    tail_fraction=self.tail/self.n if self.n else None,
                    sumw=self.n, sumw2=self.n)


def rows(path):
    with gzip.open(path, "rt") as stream:
        for line in stream:
            yield json.loads(line)


def residual_from_bound(bound, truth_p):
    qop = bound[4] if bound else None
    return 1./abs(qop)/truth_p-1. if qop and math.isfinite(qop) else None


def is_good_seed(seed, truth):
    p = pnorm(seed)
    p_truth = math.sqrt(sum(x*x for x in truth["momentum_F_MeV"]))*1e-3
    theta_truth_A = math.acos(-truth["momentum_F_MeV"][0]*1e-3/p_truth)
    return bool(p and abs(p/p_truth-1.) < .1 and seed["bound"][4]*(-truth["pdg"]) > 0
                and abs(seed["bound"][3]-theta_truth_A) < 1e-4)


def layer_index(identifier):
    return (identifier >> 59) & 3, (identifier >> 57) & 3


def reduce(root, dsid, output):
    selected = {}
    for variant in VARIANTS:
        selected[variant] = {}
        for row in rows(root/"operations_analysis"/variant/"particles.jsonl.gz"):
            if row["dsid"] != dsid:
                continue
            key = (row["file_index"], row["source_entry"], row["barcode"])
            assert key not in selected[variant]
            selected[variant][key] = {k: row[k] for k in ["seed", "candidate", "truth_p_GeV",
                                                        "final_p_ratio", "selected_hit_ids"]}
    common = selected["current"].keys() & selected["strict"].keys()
    same_path = {k for k in common if selected["current"][k]["seed"] == selected["strict"][k]["seed"]
                 and selected["current"][k]["selected_hit_ids"] == selected["strict"][k]["selected_hit_ids"]}
    seed_means = {}
    good_seeds = set()
    source_tasks = [t for t in json.loads((root/"operations_tasks.json").read_text()) if t["dsid"] == dsid]
    result = dict(dsid=dsid, residual="(positive momentum norm - production truth p) / production truth p",
                  bin_edges=EDGES.tolist(), bin_interval="left closed, right open",
                  normalization="Percent of all finite entries, including entries outside the plotted range",
                  common_good95=len(common), same_seed_and_selected_hits=len(same_path), variants={})
    for variant in VARIANTS:
        hist = {name: Histogram() for name in ["seed_candidates", "ckf_target", "final", "independent_final",
                                               "common_final", "same_path_good_seed_ckf", "same_path_good_seed_final"]}
        layers = {cohort: [[{c: Histogram() for c in COMPONENTS} for _ in range(3)] for _ in range(3)]
                  for cohort in ["good95", "same_path_good_seed"]}
        by_event = defaultdict(list)
        for key in selected[variant]:
            by_event[key[:2]].append(key)
        counts = Counter()
        for task in [t for t in source_tasks if t["variant"] == variant]:
            seen_events = 0
            folder = Path(task["output"])
            assert (folder/"exit_code.txt").read_text().strip() == "0", folder
            for ordinal, records in enumerate(event_groups(read_records(folder/"stages.jsonl.gz"))):
                seen_events += 1
                assert records[0]["event_weight"] == 1
                event_key = (task["file_index"], task["skip"]+ordinal)
                truth = {r["barcode"]: r for r in records if r["stage"] == "truth" and r["signal"]
                         and r["status"] == 1 and abs(r["pdg"]) == (13 if dsid.startswith("100") else 11)}
                seeds = {r["seed"]: r for r in records if r["stage"] == "seed"}
                lookup = {(r["stage"], r["seed"], r["candidate"]): r for r in records
                          if r["stage"] in ["selected", "selected_state_trace"]}
                for record in records:
                    if record["stage"] == "truth_hit":
                        assert layer_index(record["sensor"]) == (record["station"], record["layer"])
                        counts["validated_layer_decodes"] += 1
                for seed in seeds.values():
                    barcode = match(seed)
                    t = truth.get(barcode)
                    if t and truth_fiducial(t) and sum(t["hit_counts"][1:]) >= 18:
                        p = pnorm(seed)
                        pt = math.sqrt(sum(x*x for x in t["momentum_F_MeV"]))*1e-3
                        hist["seed_candidates"].add(p/pt-1. if p else None)
                for key in by_event.get(event_key, []):
                    item = selected[variant][key]
                    seed_id, candidate = item["seed"], item["candidate"]
                    seed = seeds[seed_id]
                    if key in same_path:
                        qualifies = is_good_seed(seed, truth[key[2]])
                        if variant == "current":
                            seed_means[key] = seed["bound"]
                            if qualifies:
                                good_seeds.add(key)
                        else:
                            assert seed_means[key] == seed["bound"], (key, "seed changed")
                            assert qualifies == (key in good_seeds)
                    target = lookup[("selected", seed_id, candidate)]
                    trace = lookup[("selected_state_trace", seed_id, candidate)]
                    pt = item["truth_p_GeV"]
                    ckf = residual_from_bound(target["bound"], pt)
                    hist["ckf_target"].add(ckf)
                    hist["final"].add(item["final_p_ratio"]-1.)
                    counts["good95_particles"] += 1
                    if key in common:
                        hist["common_final"].add(item["final_p_ratio"]-1.)
                    if key in good_seeds:
                        hist["same_path_good_seed_ckf"].add(ckf)
                        hist["same_path_good_seed_final"].add(item["final_p_ratio"]-1.)
                        counts["same_path_good_seed_particles"] += 1
                    for step in trace["steps"]:
                        if not step["measurement_flag"]:
                            continue
                        station, layer = layer_index(step["hit"])
                        assert station in [1, 2, 3] and layer in [0, 1, 2]
                        for component in COMPONENTS:
                            value = residual_from_bound(step.get(component), pt)
                            layers["good95"][station-1][layer][component].add(value)
                            if key in good_seeds:
                                layers["same_path_good_seed"][station-1][layer][component].add(value)
                counts["events"] += 1
            assert seen_events == task["events"], task["shard"]
            print("Reduced", task["shard"], flush=True)
        assert counts["good95_particles"] == len(selected[variant])
        independent_events = set()
        for row in rows(root/"independent_merged"/variant/"particles.jsonl.gz"):
            if row["dsid"] != dsid:
                continue
            independent_events.add((row["file_index"], row["source_entry"]))
            if row["fiducial"] and row["truth_H18"] and row["final"] and row["final"]["good95"]:
                assert row["event_weight"] == 1
                hist["independent_final"].add(row["final"]["p_ratio"]-1.)
        assert counts["events"] == len(independent_events) == 10000
        result["variants"][variant] = dict(counts=dict(counts), histograms={k: h.result() for k, h in hist.items()},
            layers={cohort: [[{c: h.result() for c, h in cell.items()} for cell in station] for station in grid]
                    for cohort, grid in layers.items()})
    result["same_path_good_seed_particles"] = len(good_seeds)
    a, b = [result["variants"][v]["histograms"]["seed_candidates"] for v in VARIANTS]
    assert a == b, "Material lookup must preserve seed candidate distributions"
    for variant in VARIANTS:
        assert result["variants"][variant]["counts"].get("same_path_good_seed_particles", 0) == len(good_seeds)
        expected = json.loads((root/"independent_merged"/variant/"summary.json").read_text())["samples"][dsid]["good95"]
        observed = result["variants"][variant]["histograms"]["independent_final"]
        assert (observed["n"], observed["tail_half"]) == (expected["n"], expected["tail_half"])
    result["source_tasks"] = source_tasks
    (output/"histograms.json").write_text(json.dumps(result, indent=2)+"\n")


def reduce_reference(root, dsid, output):
    """Use source-macro axes/counting while retaining the frozen S04 selection."""
    selected = {}
    for variant in VARIANTS:
        selected[variant] = {}
        for row in rows(root/"operations_analysis"/variant/"particles.jsonl.gz"):
            if row["dsid"] == dsid:
                key = (row["file_index"], row["source_entry"], row["barcode"])
                selected[variant][key] = {k: row[k] for k in ["seed", "candidate", "selected_hit_ids"]}
    common = selected["current"].keys() & selected["strict"].keys()
    same_path = {k for k in common if selected["current"][k]["seed"] == selected["strict"][k]["seed"]
                 and selected["current"][k]["selected_hit_ids"] == selected["strict"][k]["selected_hit_ids"]}
    result = dict(dsid=dsid, source_commit="33ca2281b961935456dfd4dc9fc06235e24b551e", variants={},
                  selection="S04 fiducial+H18, SDO>=0.7 and final Good95 for CKF; not the source macro quality cuts",
                  truth="Production truth; source seed writer instead uses propagated truth at the seed plane")
    tasks = [t for t in json.loads((root/"operations_tasks.json").read_text()) if t["dsid"] == dsid]
    for variant in VARIANTS:
        seeds_hist = {"relative_p": Histogram(), "abs_qop_difference": Histogram(np.linspace(-.03, .03, 101))}
        layers = {cohort: [[{c: Histogram(np.linspace(-limit, limit, 101)) for c in COMPONENTS}
                            for _ in range(3)] for _ in range(3)]
                  for cohort, limit in [("good95", 1.2), ("same_path_good_seed", 1.1)]}
        by_event = defaultdict(list)
        for key in selected[variant]:
            by_event[key[:2]].append(key)
        counts = Counter()
        for task in [t for t in tasks if t["variant"] == variant]:
            folder = Path(task["output"])
            for ordinal, records in enumerate(event_groups(read_records(folder/"stages.jsonl.gz"))):
                counts["events"] += 1
                truth = {r["barcode"]: r for r in records if r["stage"] == "truth" and r["signal"]
                         and r["status"] == 1 and abs(r["pdg"]) == (13 if dsid.startswith("100") else 11)}
                seeds = {r["seed"]: r for r in records if r["stage"] == "seed"}
                traces = {(r["seed"], r["candidate"]): r for r in records if r["stage"] == "selected_state_trace"}
                for seed in seeds.values():
                    t = truth.get(match(seed))
                    if t and truth_fiducial(t) and sum(t["hit_counts"][1:]) >= 18:
                        pt = math.sqrt(sum(x*x for x in t["momentum_F_MeV"]))*1e-3
                        p = pnorm(seed)
                        seeds_hist["relative_p"].add(p/pt-1. if p else None)
                        seeds_hist["abs_qop_difference"].add(abs(seed["bound"][4])-1./pt if p else None)
                for key in by_event.get((task["file_index"], task["skip"]+ordinal), []):
                    item = selected[variant][key]
                    seed = seeds[item["seed"]]
                    pt = math.sqrt(sum(x*x for x in truth[key[2]]["momentum_F_MeV"]))*1e-3
                    good_seed = key in same_path and is_good_seed(seed, truth[key[2]])
                    counts["good95_particles"] += 1
                    counts["same_path_good_seed_particles"] += good_seed
                    for step in traces[(item["seed"], item["candidate"])]["steps"]:
                        if not step["measurement_flag"]:
                            continue
                        station, layer = layer_index(step["hit"])
                        assert step["geometry_id"] >> 56 == {1: 25, 2: 27, 3: 29}[station]
                        assert (step["geometry_id"] >> 36) & 0xfff == (layer+1)*2
                        counts["validated_ACTS_volume_layer_mapping"] += 1
                        for component in COMPONENTS:
                            value = residual_from_bound(step.get(component), pt)
                            layers["good95"][station-1][layer][component].add(value)
                            if good_seed:
                                layers["same_path_good_seed"][station-1][layer][component].add(value)
            print("Source axes", task["shard"], flush=True)
        assert counts["events"] == 10000 and counts["good95_particles"] == len(selected[variant])
        result["variants"][variant] = dict(counts=dict(counts), seed={k: h.result() for k, h in seeds_hist.items()},
            layers={cohort: [[{c: h.result() for c, h in cell.items()} for cell in station] for station in grid]
                    for cohort, grid in layers.items()})
    assert result["variants"]["current"]["seed"] == result["variants"]["strict"]["seed"]
    (output/"histograms.json").write_text(json.dumps(result, indent=2)+"\n")


def draw_hist(ax, histogram, color, label, linestyle="-", raw=False):
    edges = np.asarray(histogram.get("bin_edges", EDGES))
    scale = 1. if raw else (100./histogram["n"] if histogram["n"] else 0.)
    values = np.asarray(histogram["counts"])*scale
    ax.step(edges, np.r_[values, values[-1]], where="post", color=color,
            linestyle=linestyle, lw=1.35, label=label)
    return max(values)


def style_axis(ax, states=False):
    ax.axvline(0, color="0.45", lw=.7, ls=":")
    ax.axvline(-.5, color="0.45", lw=.7, ls="--")
    ax.set_xlim(-2, 2)
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(("Measurement states" if states else "Entries")+" [% / 0.04]")
    ax.grid(axis="y", alpha=.18)


def stat_label(h):
    return f"N={h['n']:,}, tail={100*h['tail_fraction']:.2f}%" if h["n"] else "N=0"


def overview(data):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    for ax, name, title in zip(axes, ["seed_candidates", "ckf_target", "final"],
                               ["Initial seed candidates", "Selected CKF target", "Final track after refit/fallback"]):
        ymax = 0
        for variant in VARIANTS:
            h = data["variants"][variant]["histograms"][name]
            ymax = max(ymax, draw_hist(ax, h, COLORS[variant], LABELS[variant]+"\n"+stat_label(h),
                                      "--" if variant == "strict" else "-"))
        style_axis(ax)
        ax.set_ylim(0, ymax*1.36 if ymax else 1)
        ax.set_title(title)
        ax.legend(loc="upper left", fontsize=8, framealpha=.93)
    fig.suptitle(NAMES[data["dsid"]]+" | Momentum residuals before and after material-map rotation", fontsize=14, y=.99)
    fig.text(.5, .89, "Seed: all matched candidates in fiducial+H18 truth; CKF/final: each configuration's matched Good95 tracks",
             ha="center", fontsize=10)
    fig.text(.5, .025, "Rotated detector geometry in both runs; production-truth reference. Dashed vertical line: p/p(prod) = 0.5.\n"
             "Unit weights; normalization includes under/overflow. All initial seed candidates are identical between the two runs.",
             ha="center", fontsize=9)
    fig.subplots_adjust(top=.78, bottom=.20, left=.065, right=.985, wspace=.27)
    return fig


def layer_figure(data, cohort, reference=None):
    fig, axes = plt.subplots(3, 6, figsize=(22, 11.5))
    for station in range(3):
        for layer in range(3):
            ymax = 0
            for column, variant in enumerate(VARIANTS):
                ax = axes[station, column*3+layer]
                cell = (reference or data)["variants"][variant]["layers"][cohort][station][layer]
                for component, color in COMPONENTS.items():
                    ymax = max(ymax, draw_hist(ax, cell[component], color, component, raw=reference is not None))
                style_axis(ax, states=True)
                if reference:
                    edges = cell["smoothed"]["bin_edges"]
                    ax.set_xlim(edges[0], edges[-1])
                    ax.set_ylabel("Measurement states / bin")
                ax.set_title(f"Station {station+1}, layer {layer}", fontsize=11)
                ax.tick_params(labelsize=9)
                h = cell["smoothed"]
                text = f"N(states) = {h['n']:,}\nTail P/F/S = "+" / ".join(
                    f"{100*cell[c]['tail_fraction']:.1f}%" if cell[c]["n"] else "n/a" for c in COMPONENTS)
                ax.text(.025, .96, text, transform=ax.transAxes, va="top", fontsize=8,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=.9))
            for column in range(2):
                axes[station, column*3+layer].set_ylim(0, ymax*1.32 if ymax else 1)
    good_seed = cohort == "same_path_good_seed"
    subtitle = (f"Common matched Good95 truth, identical seeds/hit sequences; good-seed cuts; N(particles) = {data['same_path_good_seed_particles']:,} in each run"
                if good_seed else "Fiducial+H18 truth, matched Good95 final tracks in each configuration; two strip sides pooled per physical layer")
    prefix = "Source binning, S04 selection | " if reference else ""
    fig.suptitle(prefix+NAMES[data["dsid"]]+" | "+("Good-seed CKF momentum residuals" if good_seed else "CKF momentum residuals at nine tracker layers"),
                 fontsize=17, y=.994)
    fig.text(.5, .953, subtitle, ha="center", fontsize=11)
    for x, variant in [(.265, "current"), (.755, "strict")]:
        fig.text(x, .914, LABELS[variant], color=COLORS[variant], ha="center", fontsize=15, weight="bold")
    fig.legend(handles=[Line2D([0], [0], color=color, lw=2, label=c.title()) for c, color in COMPONENTS.items()],
               loc="upper center", bbox_to_anchor=(.5, .912), ncol=3, frameon=False, fontsize=11)
    footer = "Rotated detector geometry in both runs; production-truth reference; candidate-isolated smoothing; both strip sides pooled."
    footer += (r" Good seed: $|\delta p|<0.1$, correct charge, $|\Delta\theta_A|<10^{-4}$ rad." if good_seed else "")
    note = ("Raw counts and source binning; S04 selection and production truth are retained. Source macro quality cuts / propagated seed truth are not applied."
            if reference else "Each curve is normalized to all its finite measurement states, including residuals outside [-2, 2). Paired panels share y limits.")
    fig.text(.5, .028, footer+"\n"+note,
             ha="center", fontsize=10)
    fig.add_artist(Line2D([.507, .507], [.10, .85], transform=fig.transFigure, color=".6", ls=":", lw=1))
    fig.subplots_adjust(top=.85, bottom=.115, left=.045, right=.99, wspace=.39, hspace=.52)
    return fig


def source_seed_figure(data, reference):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, field, xlabel, title in zip(axes, ["relative_p", "abs_qop_difference"],
                                      [XLABEL, r"$|q/p_{\mathrm{seed}}|-|q/p_{\mathrm{prod}}|$ [GeV$^{-1}$]"],
                                      ["Relative momentum residual (PDF / commented formula)",
                                       "Absolute-curvature difference (active GitHub formula)"]):
        ymax = 0
        for variant in VARIANTS:
            h = reference["variants"][variant]["seed"][field]
            label = LABELS[variant]+f" (N={h['n']:,}, outside={h['underflow']+h['overflow']:,})"
            ymax = max(ymax, draw_hist(ax, h, COLORS[variant], label,
                                      "--" if variant == "strict" else "-", raw=True))
        ax.set(xlabel=xlabel, ylabel="Seed candidates / bin", title=title)
        edges = h["bin_edges"]
        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(0, ymax*1.3)
        ax.axvline(0, color=".5", ls=":", lw=.8)
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(axis="y", alpha=.18)
    fig.suptitle(NAMES[data["dsid"]]+" | Seed observables checked against test_pResolution", fontsize=14, y=.98)
    fig.text(.5, .045, "Raw counts; matched S04 fiducial+H18 seed candidates. Both runs use rotated detector geometry.\n"
             "Truth here is production momentum. The source seed ROOT writer uses truth propagated to the seed plane.", ha="center", fontsize=10)
    fig.subplots_adjust(top=.82, bottom=.21, left=.08, right=.99, wspace=.30)
    return fig


def independent_figure(datasets):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for ax, data in zip(axes.flat, datasets):
        ymax = 0
        for variant in VARIANTS:
            h = data["variants"][variant]["histograms"]["independent_final"]
            ymax = max(ymax, draw_hist(ax, h, COLORS[variant], LABELS[variant]+"\n"+stat_label(h)))
        style_axis(ax)
        ax.set_ylim(0, ymax*1.4)
        ax.set_title(NAMES[data["dsid"]], fontsize=11)
        ax.legend(fontsize=8, loc="upper left")
    axes[1, 2].axis("off")
    axes[1, 2].text(.02, .8, "Independent validation\n10,000 reserved source events per sample\nFiducial+H18 truth, matched Good95 tracks\n\nSame source events before/after\nEach configuration uses its own selected tracks\n\nResidual: positive momentum norm vs production truth\nUnit weights; finite-entry normalization\nUnder/overflow remain in the normalization\nVertical dashed line: p/p(prod) = 0.5", fontsize=11, va="top")
    fig.suptitle("Final momentum residuals | Independent validation of material-map rotation", fontsize=16)
    fig.subplots_adjust(top=.90, bottom=.10, left=.065, right=.99, wspace=.26, hspace=.45)
    return fig


def plot(root, reductions, output, source_reductions):
    datasets = [json.loads((reductions/dsid/"histograms.json").read_text()) for dsid in DSIDS]
    source_views = {d: json.loads((source_reductions/d/"histograms.json").read_text()) for d in DSIDS}
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "pdf.fonttype": 42})
    inventory = []
    def save(fig, stem, pdf):
        fig.savefig(output/(stem+".png"), dpi=180)
        fig.savefig(output/(stem+".pdf"))
        pdf.savefig(fig)
        inventory.append(stem)
        plt.close(fig)
    with PdfPages(output/"momentum_material_comparison.pdf") as pdf:
        save(independent_figure(datasets), "independent_final_residuals", pdf)
        for data in datasets:
            d = data["dsid"]
            save(overview(data), d+"_seed_ckf_final_residuals", pdf)
            save(layer_figure(data, "good95"), d+"_ckf_nine_layers", pdf)
            save(layer_figure(data, "same_path_good_seed"), d+"_good_seed_nine_layers", pdf)
            reference = source_views[d]
            for variant in VARIANTS:
                assert reference["variants"][variant]["counts"]["same_path_good_seed_particles"] == data["same_path_good_seed_particles"]
            save(source_seed_figure(data, reference), d+"_source_seed_views", pdf)
            save(layer_figure(data, "good95", reference), d+"_source_binning_nine_layers", pdf)
            save(layer_figure(data, "same_path_good_seed", reference), d+"_source_binning_good_seed_nine_layers", pdf)
    summary = dict(reference="docs/problem.pdf, PDF page 3", figures=inventory, pages=len(inventory),
                   material_comparison="Same rotated detector geometry; current22 vs strictly rotated current22 material map",
                   sample_events_per_dsid=10000, discovery_events_per_variant=50000, independent_events_per_variant=50000,
                   definition=datasets[0]["residual"], bin_edges=EDGES.tolist(),
                   selections={"seed_candidates": "Fiducial+H18 signal truth, SDO majority>=0.7, all matched seed candidates",
                               "good95": "Each variant's fiducial+H18 matched Good95 final tracks, final ancestry for CKF",
                               "same_path_good_seed": "Common Good95 truth, identical seed mean and selected measurement IDs; |delta p|<0.1, correct charge, |delta theta_A|<1e-4 rad"},
                   normalization=datasets[0]["normalization"],
                   source_binning_views=source_views,
                   samples={d["dsid"]: {k: v for k, v in d.items() if k != "source_tasks"} for d in datasets})
    expected = json.loads((root/"paired_independent.json").read_text())
    for d in datasets:
        for variant, reference in [("current", "baseline"), ("strict", "alternative")]:
            h = d["variants"][variant]["histograms"]["independent_final"]
            ref = expected[d["dsid"]]["distributions"][reference]
            assert h["n"] == ref["n"] and h["tail_half"] == ref["tails"]["0.5"]
    summary["checks"] = dict(independent_counts_match_frozen_S04=True, identical_seed_candidate_histograms=True,
                             histogram_conservation=True, same_path_seed_means_identical=True)
    (output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    index = ["# S04 momentum figure index", "", "Reference: problem.pdf page 3; GitHub test_pResolution commit 33ca2281b961935456dfd4dc9fc06235e24b551e.",
             "", "Same rotated detector geometry; only material-map coordinates differ. Source-binning plots retain the S04 sample selection and production truth.",
             "", "[Complete PDF](momentum_material_comparison.pdf)", ""]
    index.extend(f"- [{name} PNG]({name}.png) / [PDF]({name}.pdf)" for name in inventory)
    (output/"figure_index.md").write_text("\n".join(index)+"\n")
    print(json.dumps(dict(pages=len(inventory), checks=summary["checks"])), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["reduce", "reference-reduce", "plot"], required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dsid", choices=DSIDS)
    parser.add_argument("--reductions", type=Path)
    parser.add_argument("--source-reductions", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.phase == "reduce":
        assert args.dsid
        reduce(args.root, args.dsid, args.output)
    elif args.phase == "reference-reduce":
        assert args.dsid
        reduce_reference(args.root, args.dsid, args.output)
    else:
        assert args.reductions and args.source_reductions
        plot(args.root, args.reductions, args.output, args.source_reductions)


if __name__ == "__main__":
    main()
