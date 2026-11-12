#!/usr/bin/env python3
"""Check actual ACTS update/smoother identities on the final candidate path."""
import argparse
from collections import Counter
import gzip
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def qop(vector):return float(vector[4])


def analyze(rows, case):
    particle = case["particle"]
    seed, candidate = (particle["final"][k] for k in ("seed", "candidate"))
    truth = particle["truth_charge"] / particle["truth_p_GeV"]
    result = dict(category=case["category"], dsid=particle["dsid"],
                  file_index=particle["file_index"], source_entry=particle["source_entry"],
                  barcode=particle["barcode"], seed=seed, candidate=candidate,
                  truth_qop_production=truth, stages={})
    for stage in ("selected_state", "refit_state"):
        states = list(reversed([r for r in rows if r["stage"]==stage and r["seed"]==seed and r["candidate"]==candidate]))
        steps = []
        for i, state in enumerate(states):
            step = {k:state[k] for k in ("state", "acts_state_index", "acts_previous_index", "geometry_id", "measurement_flag", "outlier_flag", "hole_flag")}
            step.update(x_A_mm=state["surface"][0][3], hit=state.get("hit"))
            for component in ("predicted", "filtered", "smoothed"):
                if component in state:
                    value = qop(state[component])
                    step[component] = dict(qop=value, p_GeV=abs(1/value),
                                           p_ratio=abs(truth/value), wrong=value*truth<0)
            if i and "filtered" in states[i-1] and "predicted" in state:
                step["prediction_delta_qop"] = qop(state["predicted"])-qop(states[i-1]["filtered"])
                step["prediction_delta_p_GeV"] = abs(1/qop(state["predicted"]))-abs(1/qop(states[i-1]["filtered"]))
            if state["measurement_flag"] and all(k in state for k in ("predicted", "filtered", "measurement", "projector")):
                predicted=np.asarray(state["predicted"]);filtered=np.asarray(state["filtered"])
                covariance=np.asarray(state["predicted_cov"]);h=np.asarray(state["projector"])
                residual=np.asarray(state["measurement"])-h@predicted
                innovation=h@covariance@h.T+np.asarray(state["measurement_cov"])
                gain=np.linalg.solve(innovation,h@covariance).T
                expected=gain@residual
                observed=filtered-predicted
                step["update"]=dict(residual_mm=residual.tolist(),innovation_mm2=innovation.tolist(),
                    gain_qop=gain[4].tolist(),delta_qop=float(observed[4]),expected_delta_qop=float(expected[4]),
                    qop_identity_error=float(observed[4]-expected[4]),
                    normalized_qop_identity_error=float(abs(observed[4]-expected[4])/max(abs(observed[4]),abs(truth),1e-12)),
                    chi2=float(residual@np.linalg.solve(innovation,residual)))
            if i+1<len(states) and "smoothed" in state and "filtered" in state:
                nxt=states[i+1]
                if all(k in nxt for k in ("predicted_cov", "jacobian", "smoothed", "predicted")):
                    cf=np.asarray(state["filtered_cov"]);cp=np.asarray(nxt["predicted_cov"]);j=np.asarray(nxt["jacobian"])
                    smooth_gain=np.linalg.solve(cp,j@cf).T
                    delta=np.asarray(nxt["smoothed"])-np.asarray(nxt["predicted"])
                    delta[2]=(delta[2]+math.pi)%(2*math.pi)-math.pi
                    expected=(smooth_gain@delta)[4]
                    observed=qop(state["smoothed"])-qop(state["filtered"])
                    step["smoothing"]=dict(delta_qop=observed,expected_delta_qop=float(expected),
                        normalized_qop_identity_error=float(abs(observed-expected)/max(abs(observed),abs(truth),1e-12)))
            steps.append(step)
        result["stages"][stage]=steps
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    results=[];errors=[]
    for folder in sorted(args.input.iterdir()):
        if not folder.is_dir():continue
        assert (folder/"exit_code.txt").read_text().strip()=="0"
        task=json.loads((folder/"task.json").read_text())
        rows=[json.loads(line) for line in gzip.open(folder/"stages.jsonl.gz","rt")]
        for case in task["cases"]:
            result=analyze(rows,case);results.append(result)
            final=next(r for r in rows if r["stage"]=="final" and r["seed"]==result["seed"] and r["candidate"]==result["candidate"])
            original=case["particle"]["final"]
            replay_p=math.sqrt(sum(v*v for v in final["momentum_A_GeV"]))
            assert math.isclose(replay_p,original["p_GeV"],rel_tol=1e-12)
            assert final["charge"]==original["charge"]
            for stage,steps in result["stages"].items():
                for step in steps:
                    for kind in ("update","smoothing"):
                        if kind in step:errors.append((kind,step[kind]["normalized_qop_identity_error"]))
    (args.output/"traces.json").write_text(json.dumps(results,indent=2)+"\n")
    summary=dict(cases=len(results),replay_matches_baseline=True,
                 maximum_normalized_qop_errors={kind:max(e for k,e in errors if k==kind) for kind in ("update","smoothing")},
                 categories=dict(Counter(r["category"] for r in results)))
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    for result in results:
        fig,axes=plt.subplots(1,2,figsize=(12,3.7),constrained_layout=True)
        for ax,(stage,steps) in zip(axes,result["stages"].items()):
            for component in ("predicted","filtered","smoothed"):
                usable=[s for s in steps if component in s]
                ax.plot([s["x_A_mm"] for s in usable],[s[component]["qop"] for s in usable],".-",label=component)
            ax.axhline(result["truth_qop_production"],color="black",linestyle="--",label="Production truth")
            ax.set(xlabel="z_F = x_A [mm]",ylabel="q/p [GeV^-1]",title=stage.replace("_state",""));ax.legend(fontsize=7)
        fig.suptitle(f"{result['dsid']}, entry {result['source_entry']}, barcode {result['barcode']}: {result['category']}")
        fig.savefig(args.output/f"{result['dsid']}_f{result['file_index']}_e{result['source_entry']}_{result['barcode']}_{result['category']}.png",dpi=130)
        plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=="__main__":main()
