#!/usr/bin/env python3
"""Trace the remaining local tails with fixed truth hits and true sensor p bounds."""
import argparse
from collections import Counter,defaultdict
import gzip
import json
from pathlib import Path

import numpy as np

from pilot import event_groups,read_records
from momentum_simulation import sensor_bounds,local_bound


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--tasks',type=Path,required=True)
    parser.add_argument('--particles',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    particles=defaultdict(dict)
    with gzip.open(args.particles,'rt') as stream:
        for row in map(json.loads,stream):
            if row['mode']=='T4' and row['final'] and row['final']['good95'] and row['fiducial'] and row['truth_H18']:
                particles[row['grid']][row['source_entry']]=row
    result={}
    for task in json.loads(args.tasks.read_text()):
        grid=task['shard'];sensors=sensor_bounds(Path(task['sensor_truth']))
        counts=Counter();examples=[];losses=defaultdict(list)
        for ordinal,records in enumerate(event_groups(read_records(Path(task['output'])/'T4/stages.jsonl.gz'))):
            row=particles[grid].get(ordinal)
            if row is None:continue
            event=records[0]['event'];truth=sensors[event]
            first=min(truth,key=lambda s:s['z_entry_mm']);last=max(truth,key=lambda s:s['z_entry_mm'])
            label=row['sensor_local_class'];loss=last['upper_GeV']/first['lower_GeV']
            losses[label].append(loss);counts['good95']+=1;counts[label]+=1
            if label!='guaranteed_local_tail':continue
            candidates=[];stage_finals=[]
            for trace in (r for r in records if r['stage'] in ['selected_state_trace','refit_state_trace']
                          and r['seed']==row['final']['seed'] and r['candidate']==row['final']['candidate']):
                steps=trace['steps'];stage=trace['trace_stage']
                for operation,ordered,a,b in [('update',steps,'predicted','filtered'),('smoothing',reversed(steps),'filtered','smoothed')]:
                    for step in ordered:
                        if not step['measurement_flag'] or a not in step or b not in step:continue
                        bound=local_bound(truth,step['z_F_mm'])
                        if not bound:continue
                        pa=1/abs(step[a][4]);pb=1/abs(step[b][4])
                        if pa>=.5*bound['upper_GeV'] and pb<.5*bound['lower_GeV']:
                            candidates.append(dict(stage=stage,operation=operation,z_F_mm=step['z_F_mm'],
                                p_before_GeV=pa,p_after_GeV=pb,truth_lower_GeV=bound['lower_GeV'],truth_upper_GeV=bound['upper_GeV'],
                                innovation=step.get('update_residual_mm'),gain_qop=step.get('update_gain_qop')))
                stage_finals.append(dict(stage=stage,first_sensor_smoothed_GeV=1/abs(steps[0]['smoothed'][4]) if steps and 'smoothed' in steps[0] else None))
            first_crossing=candidates[0] if candidates else None
            counts['first_local_crossing:'+ (first_crossing['stage']+':'+first_crossing['operation'] if first_crossing else 'not_bracketed_at_measurement')]+=1
            examples.append(dict(source_entry=ordinal,final_p_GeV=row['final']['p_GeV'],production_p_GeV=row['truth_p_GeV'],
                final_sensor_bound=row['sensor_bound'],first_sensor=first,last_sensor=last,
                last_over_first_truth_upper_over_lower=loss,first_local_crossing=first_crossing,
                measured_crossings=candidates,stage_first_sensor_parameters=stage_finals))
        result[grid]=dict(counts=dict(counts),loss_ratio_quantiles={k:np.quantile(v,[0,.16,.5,.84,1]).tolist() for k,v in losses.items()},local_tail_events=examples)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
