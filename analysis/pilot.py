#!/usr/bin/env python3
"""Merge event-level audit shards and quantify the three baseline observables."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm


def wilson(k, n, level=.95):
    if not n:
        return [None, None]
    z = norm.ppf((1 + level) / 2)
    center = (k + z*z/2) / (n + z*z)
    half = z * math.sqrt(k*(n-k)/n + z*z/4) / (n + z*z)
    return [center-half, center+half]


def read_records(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        for line in stream:
            yield json.loads(line)


def pnorm(record):
    v = record.get("momentum_A_GeV")
    return float(np.linalg.norm(v)) if v and all(x is not None for x in v) else None


def good_track(record, radius=95.):
    ndof = record.get("ndof", 0)
    return (ndof > 0 and record.get("chi2", float("inf"))/ndof < 25.
            and record.get("layers", 0) >= 7 and ndof+5 >= 12
            and record.get("rmax_mm", float("inf")) < radius)


def truth_fiducial(record):
    positions = record["station_position_F_mm"]
    vertex = record.get("vertex_F_mm")
    return (vertex is not None and math.hypot(*vertex[:2]) < 100.
            and all(p is not None and math.hypot(*p[:2]) < 100. for p in positions[1:]))


def match(record):
    counts = record.get("matches", [])
    total = len(record.get("hits", []))
    if not counts or total == 0 or counts[0][1]/total < .7:
        return None
    return counts[0][0]


def process_event(key, records, dsid):
    truth = {r["barcode"]: r for r in records if r["stage"] == "truth" and r["signal"]
             and r["status"] == 1
             and abs(r["pdg"]) == (13 if dsid.startswith("100") else 11)}
    stages = defaultdict(list)
    for r in records:
        stages[r["stage"]].append(r)
    finals = stages["final"]
    weight = stages["event"][0]["event_weight"] if stages["event"] else 1.
    final_matches = Counter(match(r) for r in finals)
    rows = []
    for barcode, particle in truth.items():
        momentum = np.asarray(particle["momentum_F_MeV"])*.001
        truth_p = float(np.linalg.norm(momentum))
        truth_q = -1 if particle["pdg"] > 0 else 1
        fiducial = truth_fiducial(particle)
        strict_hits = sum(particle["hit_counts"][1:]) >= 18
        row = dict(dsid=dsid, run=key[0], lumi=key[1], event=key[2], barcode=barcode, event_weight=weight, truth_p_GeV=truth_p,
                   truth_charge=truth_q, truth_slope_y=momentum[1]/momentum[2] if momentum[2] else None,
                   fiducial=fiducial, truth_H18=strict_hits, n_final=len(finals),
                   n_good95=sum(good_track(r) for r in finals), n_good100=sum(good_track(r,100.) for r in finals),
                   n_matched_final=final_matches[barcode])
        for stage in ["seed", "found", "selected", "converted", "refit", "final"]:
            candidates = [r for r in stages[stage] if match(r) == barcode and pnorm(r) is not None]
            if not candidates:
                row[stage] = None
                continue
            candidate = max(candidates, key=lambda r: (r["matches"][0][1], len(r.get("hits", [])), -r.get("chi2", 0), -r["candidate"]))
            p = pnorm(candidate)
            px, py, pz = candidate["momentum_A_GeV"]
            q = candidate.get("charge", math.copysign(1., candidate.get("bound", [0]*5)[4]))
            row[stage] = dict(p_GeV=p, p_ratio=p/truth_p, charge=q, wrong=q*truth_q < 0,
                              slope_y=py/px if px else None, legacy_angle=math.atan(py/pz) if pz else None,
                              legacy_pz_GeV=pz, physical_pz_GeV=px,
                              seed=candidate["seed"], candidate=candidate["candidate"],
                              good95=good_track(candidate), good100=good_track(candidate,100.), refit=candidate.get("refit"),
                              majority_hits=candidate["matches"][0][1], measurements=candidate.get("measurements"))
        rows.append(row)
    return rows, dict(events=1, sumw=weight, sumw2=weight*weight, truth=len(truth),
                      zero_track_events=int(not finals),
                      unmatched_to_signal=sum(v for k,v in final_matches.items() if k not in truth),
                      duplicate_signal_tracks=sum(max(final_matches[k]-1,0) for k in truth),
                      invalid_final_momentum=sum(pnorm(r) is None or pnorm(r) <= 0 for r in finals),
                      seeds=len(stages["seed"]), found=len(stages["found"]),
                      selected=len(stages["selected"]), final=len(finals), finder_errors=len(stages["finder_error"]),
                      refit_fallback=sum(r.get("refit") == 0 for r in finals))


def summarize(rows):
    summary = {"truth_particles": len(rows), "fiducial_H18_truth": sum(r["fiducial"] and r["truth_H18"] for r in rows)}
    for selection in ["all_matched", "fiducial_H18", "good95", "one_good95", "alma_radius100_proxy"]:
        selected = [r for r in rows if r["final"] is not None]
        if selection not in ["all_matched","alma_radius100_proxy"]:
            selected = [r for r in selected if r["fiducial"] and r["truth_H18"]]
        if selection=="alma_radius100_proxy":
            selected=[r for r in selected if r["fiducial"] and r["final"]["good100"] and r["n_good100"]==1]
        if selection in ["good95", "one_good95"]:
            selected = [r for r in selected if r["final"]["good95"]]
        if selection == "one_good95":
            selected = [r for r in selected if r["n_good95"] == 1]
        n = len(selected)
        k = sum(r["final"]["wrong"] for r in selected)
        tails = sum(r["final"]["p_ratio"] < .5 for r in selected)
        values = [r["final"]["p_ratio"] for r in selected]
        summary[selection] = dict(n=n, wrong=k, misid=k/n if n else None, misid_wilson95=wilson(k,n),
                                  tail_half=tails, tail_fraction=tails/n if n else None, tail_wilson95=wilson(tails,n),
                                  p_ratio_quantiles=np.quantile(values,[.16,.5,.84]).tolist() if n else [])
        sumw=sum(r["event_weight"] for r in selected)
        summary[selection].update(sumw=sumw, sumw2=sum(r["event_weight"]**2 for r in selected),
                                  weighted_misid=sum(r["event_weight"] for r in selected if r["final"]["wrong"])/sumw if sumw else None,
                                  weighted_tail_fraction=sum(r["event_weight"] for r in selected if r["final"]["p_ratio"]<.5)/sumw if sumw else None)
        summary[selection]["p_ratio_overflow_ge2"]=sum(v>=2 for v in values)
        summary[selection]["p_ratio_nonpositive"]=sum(v<=0 for v in values)
        summary[selection]["wrong_and_tail_half"]=sum(r["final"]["wrong"] and r["final"]["p_ratio"]<.5 for r in selected)
        if selected:
            physical=[r["final"]["slope_y"] for r in selected if r["final"]["slope_y"] is not None]
            legacy=[r["final"]["legacy_angle"] for r in selected if r["final"]["legacy_angle"] is not None]
            residuals=[r["final"]["slope_y"]-r["truth_slope_y"] for r in selected if r["final"]["slope_y"] is not None and r["truth_slope_y"] is not None]
            summary[selection]["angles"]=dict(physical_outside_10mrad=sum(abs(math.atan(v))>=.01 for v in physical),
                                             legacy_outside_10mrad=sum(abs(v)>=.01 for v in legacy),
                                             physical_exact_zero=sum(v==0 for v in physical),
                                             physical_quantiles=np.quantile(physical,[.16,.5,.84]).tolist() if physical else [],
                                             production_slope_residual_quantiles=np.quantile(residuals,[.16,.5,.84]).tolist() if residuals else [])
            windows={}
            for label,field in [("reco","p_GeV"),("truth","truth_p_GeV")]:
                window=[r for r in selected if 50<(r["final"][field] if label=="reco" else r[field])<250]
                values=[r["final"]["slope_y"] for r in window if r["final"]["slope_y"] is not None]
                windows[label]=dict(n=len(window),slope_outside_8mrad=sum(abs(v)>=.008 for v in values),
                                   slope_quantiles=np.quantile(values,[.16,.5,.84]).tolist() if values else [])
            summary[selection]["momentum_windows_50_250_GeV"]=windows
        # Resample complete events, keeping both Aprime daughters together.
        groups = defaultdict(lambda: np.zeros(3, dtype=int))
        for r in selected:
            groups[(r.get("file_index",0),r["run"],r["lumi"],r["event"])]+= [1,int(r["final"]["wrong"]),int(r["final"]["p_ratio"]<.5)]
        if groups:
            counts=np.stack(list(groups.values()))
            rng=np.random.default_rng(19023)
            estimates=[]
            for batch in range(10):
                indices=rng.integers(0,len(counts),size=(100,len(counts)))
                totals=counts[indices].sum(axis=1)
                estimates.append(totals[:,1:]/totals[:,:1])
            intervals=np.quantile(np.concatenate(estimates),[.025,.975],axis=0)
            summary[selection]["misid_event_bootstrap95"]=intervals[:,0].tolist()
            summary[selection]["tail_event_bootstrap95"]=intervals[:,1].tolist()
    denominator = summary["fiducial_H18_truth"]
    numerator = summary["good95"]["n"]
    summary["good95_efficiency"] = dict(numerator=numerator, denominator=denominator,
                                        value=numerator/denominator if denominator else None,
                                        wilson95=wilson(numerator, denominator))
    truth_groups=defaultdict(lambda:np.zeros(3,dtype=int))
    for r in rows:
        if r["fiducial"] and r["truth_H18"]:
            passed=r["final"] is not None and r["final"]["good95"]
            truth_groups[(r.get("file_index",0),r["run"],r["lumi"],r["event"])]+=[1,int(passed),int(passed and not r["final"]["wrong"])]
    if truth_groups:
        counts=np.stack(list(truth_groups.values()));rng=np.random.default_rng(19024)
        estimates=[]
        for _ in range(10):
            totals=counts[rng.integers(0,len(counts),size=(100,len(counts)))].sum(axis=1)
            estimates.append(totals[:,1:]/totals[:,:1])
        limits=np.quantile(np.concatenate(estimates),[.025,.975],axis=0)
        summary["good95_efficiency"]["event_bootstrap95"]=limits[:,0].tolist()
        summary["truth_correct_efficiency"]=dict(numerator=int(counts[:,2].sum()),denominator=denominator,
            value=float(counts[:,2].sum()/denominator),event_bootstrap95=limits[:,1].tolist())
    summary["stages_fiducial_H18"]={}
    fiducial=[r for r in rows if r["fiducial"] and r["truth_H18"]]
    for stage in ["seed","found","selected","converted","refit","final"]:
        available=[r[stage] for r in fiducial if r[stage] is not None]
        n=len(available); tails=sum(r["p_ratio"]<.5 for r in available); wrong=sum(r["wrong"] for r in available)
        summary["stages_fiducial_H18"][stage]=dict(matched=n,tail_half=tails,wrong=wrong,
            tail_fraction=tails/n if n else None,misid=wrong/n if n else None,
            p_ratio_quantiles=np.quantile([r["p_ratio"] for r in available],[.16,.5,.84]).tolist() if n else [])
    summary["good95_truth_p_bins_GeV"]=[]
    for lo,hi in [(0,50),(50,100),(100,250),(250,500),(500,None)]:
        available=[r for r in fiducial if r["final"] and r["final"]["good95"] and
                   r["truth_p_GeV"]>=lo and (hi is None or r["truth_p_GeV"]<hi)]
        n=len(available);k=sum(r["final"]["wrong"] for r in available)
        summary["good95_truth_p_bins_GeV"].append(dict(low=lo,high=hi,n=n,wrong=k,
            misid=k/n if n else None,misid_wilson95=wilson(k,n)))
    return summary


def plots(rows, output):
    dsids = sorted(set(r["dsid"] for r in rows))
    for dsid in dsids:
        sample = [r for r in rows if r["dsid"] == dsid and r["fiducial"] and r["truth_H18"]]
        fig, axes = plt.subplots(1, 3, figsize=(16, 4))
        for stage in ["seed", "selected", "refit", "final"]:
            values = [r[stage]["p_ratio"] for r in sample if r[stage]]
            axes[0].hist(values, bins=np.linspace(0,2,101), histtype="step", label=f"{stage} (n={len(values)})")
        axes[0].set(xlabel="p(reco) / p(production truth)", ylabel="Particles", title=f"DSID {dsid}: momentum stages")
        axes[0].legend(fontsize=8)
        final = [r for r in sample if r["final"] and r["final"]["good95"]]
        for name, key in [("Physical angle", "slope_y"), ("Legacy component ratio", "legacy_angle")]:
            values = [math.atan(r["final"][key]) if key=="slope_y" else r["final"][key] for r in final if r["final"][key] is not None]
            outside=sum(abs(v)>=.01 for v in values)
            axes[1].hist(values,bins=np.linspace(-.01,.01,161),histtype="step",label=f"{name} (outside={outside})")
        values=[math.atan(r["truth_slope_y"]) for r in final if r["truth_slope_y"] is not None]
        axes[1].hist(values,bins=np.linspace(-.01,.01,161),histtype="step",linestyle="--",label="Production truth (selected)")
        axes[1].set(xlabel="theta_y [rad]",ylabel="Particles",title="Raw and physical definitions (Good95)")
        axes[1].legend(fontsize=8)
        low = [r for r in final if 50 < r["final"]["p_GeV"] < 250]
        true_low = [r for r in final if 50 < r["truth_p_GeV"] < 250]
        for name, group in [("Reco p: 50-250 GeV",low),("Truth p: 50-250 GeV",true_low)]:
            values = [r["final"]["slope_y"] for r in group if r["final"]["slope_y"] is not None]
            axes[2].hist(values,bins=np.linspace(-.008,.008,129),histtype="step",label=f"{name} (n={len(values)})")
        axes[2].set(xlabel="tan(theta_y) = py_A / px_A",ylabel="Particles",title="Momentum-window migration")
        axes[2].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(output / f"{dsid}_baseline.png",dpi=160)
        plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    tasks=json.loads((args.input/"tasks.json").read_text())
    all_rows=[]; cutflows=defaultdict(Counter); seen=set(); failures=[]; resources=defaultdict(list)
    for task in tasks:
        folder=args.input/task["shard"]
        record_file=folder/"stages.jsonl.gz"
        if not (folder/"exit_code.txt").exists() or (folder/"exit_code.txt").read_text().strip()!="0":
            failures.append(task["shard"]); continue
        records=defaultdict(list)
        for r in read_records(record_file):
            assert r["schema"]==2,("unexpected audit schema",task["shard"])
            records[(r["run"],r["lumi"],r["event"])].append(r)
        assert len(records)==task["events"],(task["shard"],len(records),task["events"])
        manifest=json.loads((folder/"run.json").read_text())
        assert manifest["success"] and manifest["args"]["events"]==task["events"]
        assert manifest["file_index"]==task["file_index"] and manifest["args"]["skip"]==task["skip"]
        usage={"wall_seconds":manifest["wall_seconds"],"events":task["events"],
               "compact_bytes":record_file.stat().st_size}
        for line in (folder/"resource.txt").read_text().splitlines():
            for field,prefix in [("user_seconds","User time (seconds):"),("system_seconds","System time (seconds):"),
                                 ("max_rss_kB","Maximum resident set size (kbytes):")]:
                if line.strip().startswith(prefix):usage[field]=float(line.split(":")[-1])
        resources[task["dsid"]].append(usage)
        for ordinal,(key,event) in enumerate(records.items()):
            uid=(task["dsid"],task["file_index"],*key)
            assert uid not in seen,("duplicate event",uid)
            assert sum(r["stage"]=="event" for r in event)==1,("repeated event key in shard",uid)
            seen.add(uid)
            rows,cutflow=process_event(key,event,task["dsid"])
            for r in rows:
                r.update(file_index=task["file_index"],source_entry=task["skip"]+ordinal)
            all_rows.extend(rows); cutflows[task["dsid"]].update(cutflow)
    summary={"schema_version":1,"missing_or_failed_shards":failures,"events":len(seen),"samples":{}}
    for dsid,counts in sorted(cutflows.items()):
        sample=[r for r in all_rows if r["dsid"]==dsid]
        summary["samples"][dsid]={"cutflow":dict(counts),**summarize(sample)}
        usage=resources[dsid]
        summary["samples"][dsid]["resources"]={"jobs":len(usage),
            "summed_wall_seconds":sum(r["wall_seconds"] for r in usage),
            "summed_user_seconds":sum(r.get("user_seconds",0) for r in usage),
            "summed_system_seconds":sum(r.get("system_seconds",0) for r in usage),
            "max_rss_kB":max(r.get("max_rss_kB",0) for r in usage),
            "compact_bytes":sum(r["compact_bytes"] for r in usage)}
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    with gzip.open(args.output/"particles.jsonl.gz","wt") as stream:
        for r in all_rows: stream.write(json.dumps(r,separators=(",",":"))+"\n")
    plots(all_rows,args.output)
    print(json.dumps(summary,indent=2))
    return 1 if failures else 0


if __name__=="__main__":
    raise SystemExit(main())
