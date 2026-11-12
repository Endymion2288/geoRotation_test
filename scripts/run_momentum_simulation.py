#!/usr/bin/env python3
"""Generate/digitize once in scratch and run paired T0/T2/T3/T4 fits."""
import argparse
import csv
import gzip
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from audit_reduce import reduce_event


def groups(path):
    event=[]
    with path.open() as stream:
        for line in stream:
            row=json.loads(line)
            if row['stage']=='event' and event:yield event;event=[]
            event.append(row)
    if event:yield event


def archive(source,destination):
    compressed=source.with_name('archive_'+destination.name)
    with source.open('rb') as inp,gzip.open(compressed,'wb') as out:shutil.copyfileobj(inp,out)
    publish(compressed,destination)


def publish(source,destination):
    """Finish gzip locally, then commit and validate the EOS copy before success."""
    temporary=destination.with_name(destination.name+'.partial')
    with source.open('rb') as inp,temporary.open('wb') as out:
        shutil.copyfileobj(inp,out,8*1024*1024)
        out.flush();os.fsync(out.fileno())
    assert temporary.stat().st_size==source.stat().st_size
    with gzip.open(temporary,'rb') as check:
        while check.read(4*1024*1024):pass
    temporary.replace(destination)


def run(command,cwd,destination):
    cwd.mkdir();destination.mkdir(exist_ok=True)
    (destination/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    with (cwd/'process.log').open('w') as log:
        result=subprocess.run(['/usr/bin/time','-v','-o','resource.txt',*command],cwd=cwd,stdout=log,stderr=subprocess.STDOUT)
    archive(cwd/'process.log',destination/'process.log.gz')
    shutil.copy2(cwd/'resource.txt',destination/'resource.txt')
    (destination/'exit_code.txt').write_text(str(result.returncode)+'\n')
    if result.returncode:raise RuntimeError(f'Command failed in {cwd}: {result.returncode}')


def build_injections(path,task,output):
    lines=[];rows=[]
    for ordinal,event in enumerate(groups(path)):
        truth={r['barcode']:r for r in event if r['stage']=='truth' and r['signal'] and r['status']==1 and r['pdg']==task['pdg']}
        for seed in [r for r in event if r['stage']=='seed']:
            candidates=seed['matches'];barcode=candidates[0][0] if candidates else None
            row=dict(event_ordinal=ordinal,event=event[0]['event'],seed=seed['seed'],barcode=barcode,eligible=False)
            if barcode in truth:
                t=truth[barcode];surface=np.array(seed['surface']);v=np.array(t['vertex_F_mm']);p=np.array(t['momentum_F_MeV'])*.001
                rotation=np.array([[0,0,1],[0,1,0],[-1,0,0]])
                position=rotation@v;direction=rotation@p/np.linalg.norm(p)
                distance=float((position-surface[:3,3])@surface[:3,2])
                row.update(plane_distance_mm=distance,truth_p_GeV=float(np.linalg.norm(p)))
                if abs(distance)<1e-6:
                    local=surface[:3,:2].T@(position-surface[:3,3]);charge=-1 if t['pdg']>0 else 1
                    mean=[*local,math.atan2(direction[1],direction[0]),math.acos(direction[2]),charge/np.linalg.norm(p),0.]
                    original=seed['bound'];qop=original[4]
                    good=bool(qop and math.isfinite(qop) and abs(abs(mean[4]/qop)-1)<.1 and qop*mean[4]>0 and abs(original[3]-mean[3])<1e-4)
                    row.update(eligible=True,truth_mean=mean,original_mean=original,good_seed_T1=good,time_truth_ns=0.)
                    lines.append(f"{ordinal} {seed['seed']} "+' '.join(format(float(x),'.17g') for x in mean)+'\n')
            rows.append(row)
    override=output/'seed_override.txt';override.write_text(''.join(lines))
    local=Path('injections.jsonl.gz')
    with gzip.open(local,'wt') as out:
        for row in rows:out.write(json.dumps(row,separators=(',',':'))+'\n')
    publish(local,output/'injections.jsonl.gz')
    return override


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--tasks',type=Path,required=True);parser.add_argument('--index',type=int,required=True)
    args=parser.parse_args();task=json.loads(args.tasks.read_text())[args.index]
    output=Path(task['output']);scratch=Path.cwd();code=Path(task['runner']).parents[1]
    sim=code/'configs/momentum_simulation.py'
    command=['python3',str(sim),'simulate','--pdg',str(task['pdg']),'--p',str(task['p_GeV']),'--events',str(task['events']),'--seed',str(task['random_seed'])]
    run(command,scratch/'simulate',output/'simulate')
    archive(scratch/'simulate/sensor_truth.csv',output/'sensor_truth.csv.gz')
    shutil.copy2(scratch/'simulate/simulate_run.json',output/'simulate_run.json')
    if task.get('sidecar_only'):
        # Only event metadata changed. Compare every physical step to the
        # original mother simulation before reusing its reconstructed tracks.
        ignored={'event','source_entry','g4_event'}
        entries=set();events=set();count=0
        from itertools import zip_longest
        with gzip.open(task['original_sensor_truth'],'rt') as old, (scratch/'simulate/sensor_truth.csv').open() as new:
            for a,b in zip_longest(csv.DictReader(old),csv.DictReader(new)):
                assert a is not None and b is not None,'Sensor step count changed'
                assert {k:v for k,v in a.items() if k not in ignored}=={k:v for k,v in b.items() if k not in ignored},(count,a,b)
                entries.add(int(b['source_entry']));events.add(int(b['event']));count+=1
        assert entries==set(range(task['events'])) and len(events)==task['events'],(entries,events)
        (output/'sensor_validation.json').write_text(json.dumps(dict(physical_steps_identical=count,
            source_entries=len(entries),athena_events=len(events),original=task['original_sensor_truth']),indent=2)+'\n')
        return
    command=['python3',str(sim),'digitize','--input',str(scratch/'simulate/diagnostic.HITS.pool.root'),'--events',str(task['events'])]
    run(command,scratch/'digitize',output/'digitize')
    shutil.copy2(scratch/'digitize/digitize_run.json',output/'digitize_run.json')
    input_file=scratch/'digitize/diagnostic.RDO.pool.root'
    modes=task.get('modes',['T0','T0_current','T2','T3','T4','T3_current','T4_current'])
    extra_options={
        'T4_no_scatt_current':['--no-scattering'],
        'T4_no_eloss_current':['--no-energy-loss'],
        'T4_electron':['--electron-hypothesis'],
        'T4_seedC01':['--seed-covariance-scale','0.1'],
        'T4_interfaces':['--correct-interfaces'],
    }
    for mode in modes:
        material=task['material_current'] if mode.endswith('_current') else task['material']
        command=['python3',task['runner'],'--input',str(input_file),'--dsid',task['dsid'],'--events',str(task['events']),
            '--file-index',str(task['file_index']),'--audit','--states','--selected-states-only','--compact-only',
            '--reset-event-maps','--isolate-smoothing','--material',material]
        command+=extra_options.get(mode,[])
        if mode.startswith(('T2','T3','T4')):
            command+=['--seed-override',str(override),'--seed-override-mode','qop' if mode=='T2' else 'mean']
        if mode.startswith('T4'):command+=['--fixed-truth-hits']
        folder=scratch/mode;destination=output/mode
        run(command,folder,destination)
        if mode=='T0':override=build_injections(folder/'stages.jsonl',task,output)
        shutil.copy2(folder/'run.json',destination/'run.json')
        shutil.copy2(folder/'resource.txt',destination/'resource.txt')
        # The analysis entry has the same stable synthetic file index in all modes.
        record=dict(task,variant=mode,skip=0,shard=task['shard'],output=str(destination))
        (destination/'task.json').write_text(json.dumps(record,indent=2)+'\n')
        local=folder/'compact.jsonl.gz'
        with gzip.open(local,'wt') as out:
            for event in groups(folder/'stages.jsonl'):
                for row in reduce_event(event):out.write(json.dumps(row,separators=(',',':'),allow_nan=False)+'\n')
        publish(local,destination/'stages.jsonl.gz')
        (destination/'exit_code.txt').write_text('0\n')
    print(f"Completed {task['pdg']} {task['p_GeV']} GeV, {task['events']} events, {len(modes)} modes")


if __name__=='__main__':main()
