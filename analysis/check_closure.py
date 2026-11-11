#!/usr/bin/env python3
"""Evaluate actual C++ audit closure records and paired reconstruction observers."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path

import numpy as np


def records(folder, name="stages.jsonl"):
    return [json.loads(line) for line in (folder/name).read_text().splitlines()]


def compare_tracks(left, right, subset=False, legacy_observer=False):
    def collect(folder):
        rows=records(folder,"events.jsonl")
        # S01's observer predates the lumi column. This explicit compatibility
        # mode is only for that fixed, single input file, not C++ audit records.
        key=lambda r: (r["run"],r["event"]) if legacy_observer else (r["run"],r["lumi"],r["event"])
        output={key(r):r["tracks"] for r in rows}
        assert len(output)==len(rows),"Repeated observer event key"
        return output
    a=collect(left); b=collect(right)
    if not subset: assert a.keys()==b.keys(),"Paired event sets differ"
    max_relative=0.; n=0
    for event in sorted(a.keys() & b.keys()):
        assert len(a[event])==len(b[event]),(event,"track count")
        for x,y in zip(a[event],b[event]):
            assert x.keys()==y.keys(),(event,"schema")
            for field in x:
                if isinstance(x[field],(int,float,bool,list)):
                    aa=np.asarray(x[field],dtype=float);bb=np.asarray(y[field],dtype=float)
                    max_relative=max(max_relative,float(np.max(np.abs(aa-bb)/(1.+np.abs(aa)))))
            n+=1
    assert a.keys() & b.keys(),"No paired events"
    return dict(events=len(a.keys() & b.keys()),tracks=n,max_normalized_difference=max_relative)


def evaluate(folder):
    rows=records(folder)
    grouped={s:[r for r in rows if r["stage"]==s] for s in {r["stage"] for r in rows}}
    keys=[]
    for r in rows:
        keys.append((r["run"],r["lumi"],r["event"],r["collection"],r["stage"],r["seed"],r["candidate"],r["state"],r.get("barcode")))
    assert len(set(keys))==len(keys),"Duplicate diagnostic lineage key"
    assert all(r["schema"]==2 for r in rows),"Unexpected diagnostic schema"
    events={(r["run"],r["lumi"],r["event"]):r for r in grouped["event"]}
    for event,e in events.items():
        erows=[r for r in rows if (r["run"],r["lumi"],r["event"])==event]
        seeds=[r for r in erows if r["stage"]=="seed"]
        assert len(seeds)==e["seeds"]
        seed_ids={r["seed"] for r in seeds}
        found={r["candidate"]:r for r in erows if r["stage"]=="found"}
        selected=[r for r in erows if r["stage"]=="selected"]
        final=[r for r in erows if r["stage"]=="final"]
        assert len(selected)==len(final)
        for r in selected+final:
            assert r["candidate"] in found and r["seed"] in seed_ids
            assert found[r["candidate"]]["seed"]==r["seed"]
    output={"events":len(events),"stage_counts":dict(Counter(r["stage"] for r in rows)),"unique_keys":len(keys)}
    for name,field in [("geometry","center_error_mm"),("geometry","basis_error"),("field","field_error_T"),
                       ("conversion","position_error_mm"),("conversion","conversion_chart_error"),
                       ("conversion","momentum_relative_error"),("conversion","covariance_relative_error")]:
        values=[r[field] for r in grouped.get("closure_"+name,[])]
        if values: output[f"max_{name}_{field}"]=max(values)
    prop=grouped.get("closure_propagation",[])
    output["propagation_success"]=sum(r["success"] for r in prop)
    output["propagation_trials"]=len(prop)
    if prop:
        output["max_p_relative_drift"]=max(abs(r.get("momentum_relative_drift",float("inf"))) for r in prop)
        errors=[]
        for r in prop:
            if r["model"]!=0 or not r["success"]:continue
            direction=np.array([math.sin(math.pi/2-.001)*math.cos(.002),math.sin(math.pi/2-.001)*math.sin(.002),math.cos(math.pi/2-.001)])
            errors.append(np.linalg.norm(np.asarray(r["position_A_mm"])-2000.*direction))
        output["zero_field_max_position_error_mm"]=max(errors,default=None)
        helix_errors=[]
        for r in prop:
            if r["model"]!=1 or not r["success"]:continue
            dx=math.sin(math.pi/2-.001)*math.cos(.002)
            dy=math.sin(math.pi/2-.001)*math.sin(.002)
            dz=math.cos(math.pi/2-.001)
            omega=r["charge"]*.000299792458*.55/150.
            angle=omega*2000.
            expected=[(dx*math.sin(angle)+dz*(math.cos(angle)-1))/omega,2000.*dy,
                      (dx*(1-math.cos(angle))+dz*math.sin(angle))/omega]
            helix_errors.append(np.linalg.norm(np.asarray(r["position_A_mm"])-expected))
        output["uniform_field_max_helix_error_mm"]=max(helix_errors,default=None)
    material=grouped.get("closure_material",[])
    output["material_trials"]=material
    if material:
        output["material_energy_accounting_error_GeV"]=max(abs(r.get("momentum_change_GeV",float("inf"))-r.get("summed_delta_p_GeV",0)) for r in material)
        # ACTS PointwiseMaterialInteraction stores both p and q/p as float.
        # Sum the half-ULP bounds from these two rounded representations over
        # the actual number of surface interactions, at the 150 GeV test p.
        per_surface=.5*float(np.spacing(np.float32(150.)))+.5*float(np.spacing(np.float32(1./150.)))*150.**2
        for r in material:
            r["float_accounting_bound_GeV"]=r["interactions"]*per_surface
            r["accounting_within_float_bound"]=abs(r["momentum_change_GeV"]-r["summed_delta_p_GeV"])<=r["float_accounting_bound_GeV"]
        output["material_accounting_within_float_bound"]=all(r["accounting_within_float_bound"] for r in material)
        mass=np.float32(.1056583755)
        energy=float(np.hypot(mass,np.float32(150.)))
        drift=math.sqrt(energy*energy-float(mass*mass))-150.
        no_loss=[r for r in material if not r["mode"]&1]
        output["float_no_loss_predicted_drift_GeV"]=drift
        output["float_no_loss_max_model_error_GeV"]=max((abs(r["momentum_change_GeV"]-drift) for r in no_loss),default=None)
    gradients=grouped.get("closure_field",[])
    output["max_gradient_relative_error"]=max((r["gradient_relative_error"] for r in gradients),default=None)
    output["max_gradient_absolute_error_T_per_mm"]=max((r["gradient_error_T_per_mm"] for r in gradients),default=None)
    if prop:
        convergence=[]
        for model in [0,1,2]:
            for q in [-1,1]:
                pair={r["max_step_mm"]:r for r in prop if r["model"]==model and r["charge"]==q and r["success"]}
                if 10. in pair and 100. in pair:
                    convergence.append(dict(model=model,charge=q,position_difference_mm=float(np.linalg.norm(np.array(pair[10.]["position_A_mm"])-pair[100.]["position_A_mm"]))))
        output["step_convergence"]=convergence
    return output


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--run",type=Path,required=True)
    parser.add_argument("--paired",type=Path)
    parser.add_argument("--legacy-observer",action="store_true",help="Pair S01's fixed-file observer lacking lumi")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    output=evaluate(args.run)
    if args.paired:output["paired"]=compare_tracks(args.paired,args.run,legacy_observer=args.legacy_observer)
    args.output.write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({k:v for k,v in output.items() if k!="material_trials"},indent=2))


if __name__=="__main__":main()
