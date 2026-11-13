#!/usr/bin/env python3
"""Freeze Calypso worker scripts, compare compatible schedds, submit a task matrix."""
import argparse
import datetime
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--tasks",type=Path,required=True)
    parser.add_argument("--library",type=Path,required=True)
    parser.add_argument("--name",required=True)
    parser.add_argument("--stage",default="S04")
    parser.add_argument("--cpus",type=int,default=1)
    parser.add_argument("--memory-mb",type=int,default=3500)
    parser.add_argument("--disk-gb",type=int,default=8)
    parser.add_argument("--flavour",default="microcentury")
    parser.add_argument("--exclude-schedd",action='append',default=[],help="Schedd rejected by a recorded compatibility probe")
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    root=args.tasks.resolve().parent;code=root/(args.name+"_worker")
    code.mkdir(exist_ok=False)
    for name in ("momentum_worker.sh","momentum_task.py","audit_reduce.py","run_momentum_simulation.py","full_momentum_task.py","run_reference_task.py"):shutil.copy2(repo/"scripts"/name,code/name)
    shutil.copy2(repo/'configs/momentum_modes.py',code/'momentum_modes.py')
    (code/"momentum_worker.sh").chmod(0o755)
    load=subprocess.check_output(["condor_status","-schedd","-af","Name","TotalRunningJobs","TotalIdleJobs","MaxJobsRunning"],text=True)
    (code/"schedd_load.txt").write_text(load)
    ordinary=[]
    for line in load.splitlines():
        name,running,idle,capacity=line.split()
        if name.startswith("bigbird") and int(name[7:9])<100 and len(name.split('.')[0])==9 and name not in args.exclude_schedd:
            ordinary.append((int(running)+int(idle),int(idle),name))
    schedd=min(ordinary)[2]
    queue=subprocess.check_output(["condor_q","-name",schedd,"-totals"],text=True)
    (code/"schedd_queue.txt").write_text(queue)
    submit=Path('/afs/cern.ch/user/x/xcheng/condor/geoRotation_test')/args.name
    submit.mkdir(exist_ok=False,parents=True);(root/"condor_logs").mkdir(exist_ok=True)
    tasks=json.loads(args.tasks.read_text())
    # Create shared parents once; simultaneous worker mkdirs can race through EOS.
    # Each worker still creates its own fresh output directory.
    for task in tasks:
        Path(task['output']).parent.mkdir(parents=True, exist_ok=True)
    if any(task.get('kind') in ['analysis', 'full_sample', 'reference'] for task in tasks):
        shutil.copytree(repo/'analysis',code/'analysis',ignore=shutil.ignore_patterns('__pycache__'))
    (submit/"jobs.sub").write_text(f'''universe = vanilla
executable = {code}/momentum_worker.sh
arguments = {args.tasks.resolve()} $(Process) {args.library.resolve()} {code}/momentum_task.py
request_cpus = {args.cpus}
request_memory = {args.memory_mb}MB
request_disk = {args.disk_gb}GB
+WantOS = "el9"
+JobFlavour = "{args.flavour}"
+GeoRotationStage = "{args.stage}"
getenv = False
should_transfer_files = YES
transfer_output_files = ""
output = {root}/condor_logs/{args.name}.$(Cluster).$(Process).out
error = {root}/condor_logs/{args.name}.$(Cluster).$(Process).err
log = {root}/condor_logs/{args.name}.events.log
queue {len(tasks)}
''')
    result=subprocess.check_output(["condor_submit","-name",schedd,str(submit/"jobs.sub")],text=True)
    manifest=dict(time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),schedd=schedd,
        choice="Minimum running+idle among ordinary bigbird; same verified EL9/EOS/CVMFS requirements",
        task_file=str(args.tasks.resolve()),tasks=len(tasks),submission=result,library=str(args.library.resolve()))
    manifest['excluded_schedds']=args.exclude_schedd
    manifest['stage']=args.stage
    manifest['requested_resources']=dict(cpus=args.cpus,memory_mb=args.memory_mb,
        disk_gb=args.disk_gb,flavour=args.flavour)
    (code/"submission.json").write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
