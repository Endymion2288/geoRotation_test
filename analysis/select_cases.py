#!/usr/bin/env python3
"""Select reproducible examples using the actual final candidate ancestry."""
import argparse
import json
from pathlib import Path

from pilot import read_records


def first_stage(row, observable):
    for stage,value in row["final_lineage"].items():
        if (value["p_ratio"]<.5 if observable=="tail" else value["wrong"]):return stage


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--particles",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    selected={}
    for row in read_records(args.particles):
        final=row["final"]
        if not (row["fiducial"] and row["truth_H18"] and final and final["good95"]
                and not final["matching_ambiguous"]):continue
        categories=[]
        if .9<final["p_ratio"]<1.1 and not final["wrong"] and 50<row["truth_p_GeV"]<250:categories.append("normal")
        if final["p_ratio"]<.5:categories.append(first_stage(row,"tail")+"_tail")
        if final["wrong"]:categories.append(first_stage(row,"charge")+"_wrong_charge")
        for category in categories:
            key=(row["dsid"],category)
            if key in selected:continue
            selected[key]=dict(category=category,particle=row,replay_args=["--dsid",row["dsid"],
                "--file-index",str(row["file_index"]),"--skip",str(row["source_entry"]),"--events","1",
                "--audit","--states","--material",
                "/eos/home-x/xcheng/FASER/geoRotation_test_results/stage02_closure/material-current-source.json"])
    output=dict(selection="fiducial + H18 + Good95 + unambiguous matching; first source entry per DSID/category; normal p_truth 50-250 GeV",
                note="Production truth momentum; observed-stage first occurrence is not a CKF update/root-cause diagnosis",
                cases=list(selected.values()))
    args.output.write_text(json.dumps(output,indent=2)+"\n")
    print(f"Selected {len(selected)} reproducible cases")


if __name__=="__main__":main()
