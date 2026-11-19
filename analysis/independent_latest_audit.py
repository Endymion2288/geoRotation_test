#!/usr/bin/env python3
"""Recompute S04 audit quantities directly from the frozen worker tables."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np


DEEP_GROUPS = ("baseline", "fix_cov_loc_swap", "qop_1000_only", "qop_1000_loc_fixed",
               "seed_b_04737", "combined_phys_fix")
VALIDATION_GROUPS = ("C0_baseline", "C1_cov_loc", "C2_orig_var", "C3_seed_b", "C4_no_nan",
                     "C5_local_tan", "H1_electron_hypo", "S1_dedup_prefilter",
                     "S2_strip_field_fit", "S3_strip_only")
SAMPLES = ("110051", "110033", "100116", "100117", "100130")


def read_json(path):
    return json.loads(path.read_text())


def columns(data, values_key, names_key):
    return dict(zip(data[names_key].tolist(), data[values_key].T))


def residual(qop, p):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1.0 / (np.abs(qop) * p) - 1.0


def stats(values):
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    if not len(values):
        return {"n": 0}
    q16, median, q84 = np.quantile(values, [.16, .50, .84])
    return dict(n=len(values), median=float(median), sigma68=float((q84 - q16) / 2),
                below_m01=float(np.mean(values < -.1)), below_m05=float(np.mean(values < -.5)))


def paired_median(values, events):
    """Poisson bootstrap with one shared resampling weight per source event."""
    order = np.argsort(values, kind="stable")
    ordered = values[order]
    _, event_index = np.unique(np.asarray(events)[order], axis=0, return_inverse=True)
    rng = np.random.default_rng(110051)
    draws = []
    for _ in range(2000):
        cumulative = np.cumsum(rng.poisson(1, event_index.max() + 1)[event_index])
        lo = np.searchsorted(cumulative, (cumulative[-1] - 1) // 2, side="right")
        hi = np.searchsorted(cumulative, cumulative[-1] // 2, side="right")
        draws.append((ordered[lo] + ordered[hi]) / 2)
    return np.quantile(draws, [.025, .975]).tolist()


def load_chains(path):
    chains = defaultdict(list)
    for file_index in range(20):
        with np.load(path / f"file{file_index:02}" / "vertices.npz") as data:
            for entry, barcode, pdg, z, p, _ in data["values"]:
                chains[file_index, int(entry), int(barcode) % 1_000_000, int(pdg)].append((z, p))
    return chains


def seed_counts(data, file_index, chains):
    s = columns(data, "seeds", "seed_fields")
    e = columns(data, "events", "event_fields")
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = np.abs(s["writer_truth_qop"] / s["writer_reco_qop"]) - 1
    valid = np.isfinite(delta) & (s["writer_truth_qop"] != 0)
    pure = valid & (s["stations"] == 3) & (s["all_segs_same_barcode"] == 1)
    pure &= (s["has_secondary_seg"] == 0) & (s["min_seg_purity"] >= .8)
    low = np.zeros(len(valid), dtype=bool)
    hard = low.copy()
    missing = low.copy()
    for i in np.flatnonzero(pure):
        if np.isfinite(s["truth_pdg"][i]):
            key = (file_index, int(s["entry"][i]), int(s["writer_majority_barcode"][i]) % 1_000_000,
                   int(s["truth_pdg"][i]))
            history = sorted(chains.get(key, ()))
        else:
            history = []
        if history and history[0][1] > 0:
            ratio = min((p for z, p in history if z <= 2500), default=history[0][1]) / history[0][1]
            low[i], hard[i] = ratio >= .95, ratio < .95
        else:
            missing[i] = True
    left = valid & (delta < -.5)
    mix = valid & ~pure & (s["has_opp_primary_mix"] == 1)
    other = valid & ~pure & ~mix
    counts = dict(selected=len(valid), valid=int(valid.sum()), below_m05=int(left.sum()),
                  raw=int(e["raw_seed_combinations"].sum()),
                  raw_nonfinite=int(e["raw_nan_chi2_seeds"].sum()),
                  selected_nan=int(s["chi2_is_nan"].sum()),
                  negative_radicand=int((s["min_radicand"] < 0).sum()),
                  shared_six=int((s["max_shared_prev"] == 6).sum()))
    for name, mask in (("low", low), ("hard", hard), ("missing_chain", missing),
                       ("opposite", mix), ("other", other)):
        counts[name] = int(mask.sum())
        counts[f"tail_{name}"] = int((mask & left).sum())
    return counts, delta[low | missing]


def deep_study(root, chains):
    outputs, paths_by_group, hist_inputs = {}, {}, {}
    seed_total = defaultdict(int)
    clean_residuals = []
    input_files = {}
    for group in DEEP_GROUPS:
        same, paths, by_sensor = [], defaultdict(list), []
        for index in range(20):
            shard = root / "runs" / group / f"file{index:02}"
            task = read_json(shard / "source_task.json")
            run = read_json(shard / "run.json")
            assert task["file_index"] == index and task["skip"] == 0 and task["events"] == 1000
            assert input_files.setdefault(index, task["input"]) == task["input"]
            assert (shard / "exit_code.txt").read_text().strip() == "0"
            assert run["success"] and run["isolate_smoothing"]
            assert any("software_deep_v1/lib/libFaserActsKalmanFilter.so" in lib
                       for lib in run["loaded_acts_libraries"])
            if group == "baseline":
                with np.load(shard / "seed_forensics.npz") as source:
                    counts, clean = seed_counts(source, index, chains)
                for name, n in counts.items():
                    seed_total[name] += n
                clean_residuals.extend(clean)
            if index == 0 and group != "baseline":
                continue
            with np.load(shard / "states.npz") as source:
                s = columns(source, "values", "fields")
                t = columns(source, "tracks", "track_fields")
                track_lookup = {(int(e), int(c)): j for j, (e, c)
                                in enumerate(zip(t["entry"], t["candidate"]))}
                keep = ((s["measurement"] > 0) & (s["github"] > 0)
                        & (s["purity"] > .7) & (s["same_seed_truth"] > 0))
                selected = np.flatnonzero(keep)
                keys = np.column_stack((s["entry"][keep], s["candidate"][keep]))
                splits = np.flatnonzero(np.any(keys[1:] != keys[:-1], axis=1)) + 1
                hit_ids = source["hit_ids"]
                geometry_ids = source["geometry_ids"]
                for rows in np.split(selected, splits):
                    if not len(rows):
                        continue
                    station = s["station"][rows]
                    second = rows[station == 2]
                    if not np.any(station == 1) or not len(second):
                        continue
                    j = second[0]
                    p = s["p_production_GeV"][j]
                    truth_charge = s["truth_charge"][j]
                    pdg = int(-11 * np.sign(truth_charge)) if np.isfinite(truth_charge) else 0
                    entry, bc = int(s["entry"][j]), int(s["track_barcode"][j])
                    k = track_lookup[entry, int(s["candidate"][j])]
                    prt, flt, smt = [float(residual(s[f"q_{kind}"][j], p))
                                     for kind in ("predicted", "filtered", "smoothed")]
                    full_path = tuple(zip(geometry_ids[rows].tolist(), hit_ids[rows].tolist()))
                    key = (index, entry, bc, full_path)
                    good = s["good_seed"][j] > 0
                    paths[key].append((good, prt, flt, smt))
                    final_q = t["final_qop_GeV_inv"][k]
                    final_p = float(residual(final_q, p))
                    history = sorted(chains.get((index, entry, bc % 1_000_000, pdg), ()))
                    pre = minimum = np.nan
                    local = np.nan
                    if len(history) > 1 and history[0][1] > 0:
                        p0 = history[0][1]
                        pre = min((p for z, p in history if z <= 1200), default=p0) / p0
                        minimum = min((p for z, p in history if z <= 2500), default=p0) / p0
                        local = float(residual(s["q_smoothed"][j], p * pre))
                    same.append((index, entry, good, prt, flt, smt, final_p,
                                 t["refit_ok"][k], minimum, pre, local))
                    if index > 0 and good:
                        l0 = second[s["layer"][second] == 0]
                        for row in l0:
                            by_sensor.append((int(s["sensor_ordinal_in_layer"][row]),
                                              float(residual(s["q_predicted"][row], p))))
        data = np.asarray(same)
        held = data[:, 0] > 0
        good = data[:, 2] > 0
        selection = held & good
        outputs[group] = {
            "heldout_good": {name: stats(data[selection, i]) for name, i in
                             (("predicted", 3), ("filtered", 4), ("smoothed", 5))},
            "heldout_ungated_smoothed": stats(data[held, 5]),
        }
        if group == "baseline":
            for name, mask in (("ungated", np.ones(len(data), bool)), ("no_recorded_loss", data[:, 8] >= .995),
                               ("hard", data[:, 8] < .95), ("pre_s2_hard", data[:, 9] < .90)):
                outputs[group][name] = dict(production=stats(data[mask, 5]), pre_s2=stats(data[mask, 10]))
            with np.errstate(divide="ignore", invalid="ignore"):
                outputs[group]["good_refit_vs_smoothed"] = stats((1 + data[good, 6]) / (1 + data[good, 5]) - 1)
            outputs[group]["good_refit_success"] = float(np.mean(data[good, 7] == 1))
        paths_by_group[group] = paths
        hist_inputs[group] = data[selection, 3:6]
        if group in ("baseline", "qop_1000_only"):
            sensor = np.asarray(by_sensor)
            outputs[group]["s2_layer0_by_ordinal"] = {
                str(i): stats(sensor[sensor[:, 0] == i, 1]) for i in (0, 1, 2)}
    assert len(set(input_files.values())) == 20

    paired = {}
    baseline = paths_by_group["baseline"]
    for group in DEEP_GROUPS[1:]:
        control = paths_by_group[group]
        common = [(k, baseline[k][0], control[k][0]) for k in sorted(baseline.keys() & control.keys())
                  if k[0] > 0 and len(baseline[k]) == len(control[k]) == 1 and baseline[k][0][0]]
        events = np.asarray([k[:2] for k, _, _ in common])
        differences = np.asarray([b[1:] for _, _, b in common]) - np.asarray([a[1:] for _, a, _ in common])
        valid = np.isfinite(differences[:, 2])
        paired[group] = dict(predicted=stats(differences[:, 0]), filtered=stats(differences[:, 1]),
                             smoothed=stats(differences[:, 2]),
                             smoothed_ci95=paired_median(differences[valid, 2], events[valid]))
    return dict(jobs=120, reconstructed_events=120000, seeds=dict(seed_total),
                clean_seed=stats(clean_residuals), groups=outputs, paired=paired), hist_inputs


def validation_study(root):
    outputs, efficiency_inputs = {}, {}
    library_path = root / "software_fix_v1/lib/libFaserActsKalmanFilter.so"
    common_inputs = {}
    jobs, total_events = 0, 0
    for group in VALIDATION_GROUPS:
        for sample in SAMPLES:
            event_columns, signal_columns, seed_columns, track_columns = [], [], [], []
            group_inputs = set()
            for path in sorted((root / "runs" / group / sample).glob("*/shard_data.npz")):
                job = read_json(path.parent / "task.json")
                task = read_json(Path(job["command"][-1]))
                m = read_json(path.parent / "summary.json")
                args = m["run_meta_0"]["args"]
                assert (path.parent / "exit_code.txt").read_text().strip() == "0"
                assert str(library_path) in m["run_meta_0"]["loaded_acts_libraries"]
                assert m["run_meta_0"]["success"] and args["mask"] == "0,1,2,3,4,5"
                assert task["group"] == group and task["sample"] == sample
                ranges = {(r["input"], r["skip"], r["events"]) for r in task["sub_ranges"]}
                assert not group_inputs & ranges
                group_inputs |= ranges
                jobs += 1
                with np.load(path) as rows:
                    ev = columns(rows, "events", "event_cols")
                    assert len(ev["entry"]) == m["total_events"] == sum(r[2] for r in ranges)
                    assert len(np.unique(ev["entry"])) == len(ev["entry"])
                    total_events += len(ev["entry"])
                    if sample == "110051" and path.parent.name == "file00":
                        continue
                    event_columns.append(ev)
                    signal_columns.append(columns(rows, "signals", "signal_cols"))
                    seed_columns.append(columns(rows, "seeds", "seed_cols"))
                    track_columns.append(columns(rows, "tracks", "track_cols"))
            assert group_inputs == common_inputs.setdefault(sample, group_inputs)
            e, s, seed, t = [{key: np.concatenate([d[key] for d in blocks]) for key in blocks[0]}
                             for blocks in (event_columns, signal_columns, seed_columns, track_columns)]
            assert len(np.unique(e["entry"])) == len(e["entry"])
            assert len(e["entry"]) == (19000 if sample == "110051" else 10000)
            fid = s["fiducial"] > 0
            pass_ckf = s["n_ckf_match"][fid] > 0
            pair_identity, inv, multiplicity = np.unique(s["entry"][fid], return_inverse=True, return_counts=True)
            matched = np.bincount(inv, weights=pass_ckf)
            pair = multiplicity >= 2
            gh = (t["ndf_ckf"] > 9) & (t["chi2_ckf"] < 20 * t["ndf_ckf"]) & (t["rmax_ckf"] < 90)
            pur = gh & (t["trk_purity"] >= .7)
            same = pur & (t["seed_maj_bc"] == t["trk_maj_bc"]) & (t["trk_maj_bc"] > 0)
            seed_residual = np.full_like(seed["p_seed_GeV"], np.nan)
            np.divide(seed["p_seed_GeV"], seed["p_prod_GeV"], out=seed_residual,
                      where=seed["p_prod_GeV"] > 0)
            good = ((seed["seed_truth_valid"] > 0)
                    & (np.abs(seed_residual - 1) < .1)
                    & (seed["q_seed"] * seed["truth_q"] > 0)
                    & (np.abs(seed["theta_A"] - seed["t_theta_A"]) < 1e-4))
            good_seeds = set(zip(seed["entry"][good], seed["iseed"][good]))
            good_same = same & np.array([k in good_seeds for k in zip(t["entry"], t["iseed"])])
            refit_ok = t["refit_ok"] > 0
            good_w0 = refit_ok & np.isfinite(t["ref_w0_qop"])
            final_q = np.where(good_w0, t["ref_w0_qop"], t["ckf_w0_qop"])
            outputs[f"{group}/{sample}"] = dict(
                threshold=args["seed_max_reduced_chi2"], seeds=len(seed["entry"]),
                raw_nonfinite=int(e["raw_nan_chi2"].sum()), selected_nan=int(e["n_seeds_nan_chi2"].sum()),
                fiducial=int(fid.sum()), ckf_eff=float(np.mean(pass_ckf)),
                pair_eff=float(np.mean(matched[pair] >= 2)) if pair.any() else None,
                pair_events=int(pair.sum()), cross_track=float(np.mean(t["seed_maj_bc"][pur] != t["trk_maj_bc"][pur])),
                good_same=int(good_same.sum()), good_s2_smt=stats(residual(t["s2_first_smt_qop"][good_same], t["p_prod_GeV"][good_same])),
                final_w0=stats(residual(final_q[same], t["p_st1_GeV"][same])),
                final_vs_ckf_w0=stats(np.abs(t["ckf_w0_qop"][same & good_w0] / final_q[same & good_w0]) - 1),
                refit_w0_different_wafer=int((same & good_w0 & (t["ref_w0_geo"] != t["ckf_w0_geo"])).sum()),
                refit_w0_fallback=int((same & ~good_w0).sum()),
                s2_first_pull_unclipped_std=float(np.nanstd(t["s2_first_pull"][same])),
                s2_first_pull_clipped_std=float(np.std(t["s2_first_pull"][same & (np.abs(t["s2_first_pull"]) < 8)])),
            )
            efficiency_inputs[group, sample] = (pair_identity, multiplicity, matched)
    paired = {}
    for sample in SAMPLES:
        ids, count, base = efficiency_inputs["C5_local_tan", sample]
        for group in ("H1_electron_hypo", "S3_strip_only"):
            other_ids, other_count, other = efficiency_inputs[group, sample]
            assert np.array_equal(ids, other_ids) and np.array_equal(count, other_count)
            difference = other - base
            paired_events = count >= 2
            pair_diff = (other[paired_events] >= 2).astype(float) - (base[paired_events] >= 2)
            pair_rate = float(pair_diff.mean()) if len(pair_diff) else None
            pair_error = float(pair_diff.std() / np.sqrt(len(pair_diff))) if len(pair_diff) else None
            rate = difference.sum() / count.sum()
            error = np.linalg.norm(difference - rate * count) / count.sum()
            paired[f"{group}/{sample}"] = dict(
                diff=rate, ci95=[rate - 1.96 * error, rate + 1.96 * error],
                pair_diff=pair_rate,
                pair_ci95=[pair_rate - 1.96 * pair_error, pair_rate + 1.96 * pair_error]
                if pair_rate is not None else None,
            )
    return dict(jobs=jobs, reconstructed_events=total_events, groups=outputs, paired_efficiency=paired)


def field_closure(root):
    output = {}
    for sample in SAMPLES:
        rows = read_json(root / sample / "field_closure.json")["particles"]
        values = np.array([row["res_B0473662"] for row in rows])
        radii = np.array([row["r_max_mm"] for row in rows])
        edge = (radii >= 80) & (radii < 100)
        output[sample] = dict(
            entries=len(rows),
            total_quantiles=np.quantile(values, [.025, .16, .5, .84, .975]).tolist(),
            edge_entries=int(edge.sum()),
            edge_quantiles=np.quantile(values[edge], [.16, .5, .84]).tolist(),
        )
    return output


def plot_histograms(output, data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    edges = np.linspace(-1.2, 1.2, 101)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), sharey=True, constrained_layout=True)
    groups = (("baseline", "Baseline"), ("qop_1000_only", "Reduced q/p prior variance"),
              ("seed_b_04737", "Calibrated seed field"))
    counts = {}
    for i, (axis, stage) in enumerate(zip(axes, ("Predicted", "Filtered", "Smoothed"))):
        for group, label in groups:
            values = data[group][:, i]
            finite = values[np.isfinite(values)]
            n, _ = np.histogram(finite, edges)
            if np.any(finite == edges[-1]):
                n[-1] -= np.sum(finite == edges[-1])
            counts[f"{group}_{stage}"] = np.r_[(finite < edges[0]).sum(), n, (finite >= edges[-1]).sum()]
            axis.stairs(n, edges, label=f"{label} (N={len(finite):,})", linewidth=1.4)
        axis.set(title=f"First S2 measurement · {stage}", xlabel="(p reco − p production) / p production",
                 xlim=(-1.2, 1.2))
    axes[0].set_ylabel("Entries / bin")
    axes[-1].legend(fontsize=7, loc="upper left")
    fig.suptitle("DSID 110051 · 19 held-out files · 19,000 events · ACTS frame · good seed, same truth")
    fig.savefig(output / "station2_independent.png", dpi=170)
    fig.savefig(output / "station2_independent.pdf")
    np.savez_compressed(output / "station2_bin_counts.npz", edges=edges, **counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    deep_root = args.results / "stage04_latest_momentum"
    fix_root = args.results / "stage04_latest_validation/20261118_fix_v1"
    chains = load_chains(deep_root / "20261116_v1/truth_vertices")
    deep, hist_inputs = deep_study(deep_root / "20261118_v1", chains)
    validation = validation_study(fix_root)
    field = field_closure(fix_root / "field_closure")
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "independent_summary.json").write_text(
        json.dumps(dict(deep=deep, validation=validation, field_closure=field), indent=2, default=float) + "\n")
    plot_histograms(args.output, hist_inputs)


if __name__ == "__main__":
    main()
