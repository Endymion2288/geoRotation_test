#!/usr/bin/env python3
"""Bound local momentum with recorded G4 continuation vertices, not interpolation."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
import math
from pathlib import Path

from pilot import archive_records, event_groups, process_event


def envelope(records, barcode, pdg, z):
    chain=[r for r in records if r["stage"]=="truth" and r["barcode"]%1000000==barcode
           and r["pdg"]==pdg and r.get("vertex_F_mm") is not None]
    chain.sort(key=lambda r:r["vertex_F_mm"][2])
    if not chain:return None
    points=[dict(z_mm=r["vertex_F_mm"][2],p_GeV=math.sqrt(sum(v*v for v in r["momentum_F_MeV"]))*1e-3,
                 barcode=r["barcode"],status=r["status"]) for r in chain]
    monotonic=all(b["p_GeV"]<=a["p_GeV"]*(1+1e-6) and b["z_mm"]>=a["z_mm"] for a,b in zip(points,points[1:]))
    before=[r for r in points if r["z_mm"]<=z];after=[r for r in points if r["z_mm"]>z]
    if not before:return None
    return dict(reference_z_F_mm=z,upper_GeV=before[-1]["p_GeV"],lower_GeV=after[0]["p_GeV"] if after else 0.,
                before=before[-1],after=after[0] if after else None,chain_vertices=len(points),
                forward_monotonic=monotonic,exact_local_truth=False)


def category(p, bounds):
    if p is None:return "invalid_momentum"
    if not bounds or not bounds["forward_monotonic"]:return "invalid_chain"
    if p<.5*bounds["lower_GeV"]:return "guaranteed_local_tail"
    if p>=.5*bounds["upper_GeV"]:return "excluded_local_tail"
    return "unresolved_local_truth_interval"


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    tasks=json.loads((args.input/"index.json").read_text())["tasks"]
    summary=defaultdict(Counter);n=0
    with gzip.open(args.output/"local_bounds.jsonl.gz","wt") as output:
        for task in tasks:
            for ordinal,records in enumerate(event_groups(archive_records(args.input,task))):
                header=records[0];key=tuple(header[k] for k in ("run","lumi","event"))
                rows,_=process_event(key,records,task["dsid"])
                for row in rows:
                    final=row["final"]
                    if not (row["fiducial"] and row["truth_H18"] and final and final["good95"]):continue
                    pdg=(13 if task["dsid"].startswith('100') else 11)*(-row['truth_charge'])
                    bounds=envelope(records,row["barcode"],pdg,final["position_A_mm"][0])
                    seed_record=next(r for r in records if r["stage"]=="seed" and r["seed"]==final["seed"])
                    seed_bounds=envelope(records,row["barcode"],pdg,seed_record["position_A_mm"][0])
                    seed_vector=seed_record["momentum_A_GeV"]
                    seed_p=math.sqrt(sum(v*v for v in seed_vector)) if all(v is not None for v in seed_vector) else None
                    result=dict(dsid=task["dsid"],file_index=task["file_index"],source_entry=task["skip"]+ordinal,
                        barcode=row["barcode"],truth_p_GeV=row["truth_p_GeV"],final=final,bounds=bounds,
                        seed_p_GeV=seed_p,seed_bounds=seed_bounds,seed_local_class=category(seed_p,seed_bounds),
                        local_class=category(final["p_GeV"],bounds))
                    output.write(json.dumps(result,separators=(',',':'))+'\n');n+=1
                    counts=summary[task["dsid"]];counts['good95']+=1
                    counts[result['local_class']]+=1
                    if final['p_ratio']<.5:
                        counts['production_tail']+=1;counts['production_tail_'+result['local_class']]+=1
                        if bounds and bounds['upper_GeV']<.5*row['truth_p_GeV']:counts['production_tail_upstream_loss_below_half']+=1
                    if seed_p is None:counts['invalid_seed_momentum']+=1
                    elif seed_p<.5*row['truth_p_GeV']:
                        counts['seed_production_tail']+=1;counts['seed_production_tail_'+result['seed_local_class']]+=1
    result=dict(selected_particles=n,samples={k:dict(v) for k,v in summary.items()},
        interpretation='Monotone forward continuation p bounds; no continuous-loss interpolation; zero lower bound when no later continuation; exact sensor p unavailable')
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
