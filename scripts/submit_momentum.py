#!/usr/bin/env python3
"""Freeze S04 worker scripts, compare compatible schedds, submit a task matrix."""
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
    parser.add_argument("--exclude-schedd",action='append',default=[],help="Schedd rejected by a recorded compatibility probe")
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    root=args.tasks.resolve().parent;code=root/(args.name+"_worker")
    code.mkdir(exist_ok=False)
    for name in ("momentum_worker.sh","momentum_task.py","audit_reduce.py","run_momentum_simulation.py"):shutil.copy2(repo/"scripts"/name,code/name)
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
    if any(task.get('kind')=='analysis' for task in tasks):
        shutil.copytree(repo/'analysis',code/'analysis',ignore=shutil.ignore_patterns('__pycache__'))
    (submit/"jobs.sub").write_text(f'''universe = vanilla
executable = {code}/momentum_worker.sh
arguments = {args.tasks.resolve()} $(Process) {args.library.resolve()} {code}/momentum_task.py
request_cpus = 1
request_memory = 3500MB
request_disk = 8GB
+WantOS = "el9"
+JobFlavour = "microcentury"
+GeoRotationStage = "S04"
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
    (code/"submission.json").write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
