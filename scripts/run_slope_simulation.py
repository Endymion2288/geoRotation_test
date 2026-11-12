#!/usr/bin/env python3
"""Full Calypso G4 -> digitization -> native reconstruction on worker scratch."""
import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def run(command,folder,output,env=None):
    folder.mkdir()
    output.mkdir(parents=True,exist_ok=True)
    (output/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    with (folder/'process.log').open('w') as log:
        status=subprocess.run(['/usr/bin/time','-v','-o','resource.txt',*command],cwd=folder,
                              stdout=log,stderr=subprocess.STDOUT,env=env).returncode
    with (folder/'process.log').open('rb') as source,gzip.open(output/'process.log.gz','wb') as target:
        shutil.copyfileobj(source,target)
    shutil.copy2(folder/'resource.txt',output/'resource.txt')
    (output/'exit_code.txt').write_text(str(status)+'\n')
    if status:raise RuntimeError(f'Calypso step failed: {output}')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args()
    task=json.loads(args.task.read_text())
    code=Path(task['code'])
    output=Path(task['output'])
    scratch=Path.cwd()
    simulation=code/'configs/momentum_simulation.py'
    common=['--events',str(task['events']),'--full-geometry','--gun-z',str(task.get('gun_z_F_mm',-3112.)),
            '--no-sidecar','--pdg',str(task['pdg']),'--p',str(task['p_GeV'])]
    command=[sys.executable,str(simulation),'simulate',*common,
             '--seed',str(task['seed']),'--uniform-angles','--gun-x',str(task.get('x_mm',0.)),
             '--gun-y',str(task.get('y_mm',0.)),'--xy-half-width',str(task.get('xy_half_width_mm',30.)),
             '--theta-half-width',str(task.get('theta_half_width_rad',.005))]
    run(command,scratch/'simulate',output/'simulate')
    command=[sys.executable,str(simulation),'digitize',*common,
             '--input',str(scratch/'simulate/diagnostic.HITS.pool.root')]
    run(command,scratch/'digitize',output/'digitize')
    modes=task.get('reconstruction_modes')
    if modes is None:
        modes=[dict(variant=v,cluster_mode=m,label=v if m=='native' else v+'_'+m)
               for v in task['variants'] for m in task.get('cluster_modes',['native'])]
    for mode in modes:
        variant,label=mode['variant'],mode['label']
        source=dict(input=str(scratch/'digitize/diagnostic.RDO.pool.root'),events=task['events'],skip=0,
            sample_entry=0,file_index=task.get('file_index',0),variant=variant,dsid=str(task['pdg']),
            material=task['materials'][variant],runner=str(code/'configs/slope_tracking.py'),
            output=str(output/label),options=['--actual-covariance','--no-seed-writer'])
        if mode.get('cluster_mode')=='truth_position':source['options'].append('--truth-cluster-position')
        source['options'].extend(mode.get('options',[]))
        sourcefile=scratch/(label+'_source.json')
        sourcefile.write_text(json.dumps(source,indent=2)+'\n')
        library=Path(task['libraries'][variant])
        env=dict(os.environ,GEOROTATION_LIBRARY_DIR=str(library),
                 LD_LIBRARY_PATH=str(library)+':'+os.environ.get('LD_LIBRARY_PATH',''),
                 PYTHONPATH=str(library.parent/'python')+':'+os.environ.get('PYTHONPATH',''))
        run([sys.executable,str(code/'scripts/run_slope_task.py'),'--task',str(sourcefile)],scratch/label,output/(label+'_wrapper'),env)
    for mode in ['simulate','digitize']:
        shutil.copy2(scratch/mode/(mode+'_run.json'),output/(mode+'_run.json'))
    (output/'source_task.json').write_text(json.dumps(task,indent=2)+'\n')
    (output/'truth_reference.json').write_text(json.dumps(dict(
        reference='Exact generated G4 launch momentum at the declared physical plane',
        z_F_mm=task.get('gun_z_F_mm',-3112.),p_GeV=task['p_GeV'],full_Calypso_geometry=True,
        tracking='Native Calypso; no independent tracking or simplified detector',
        historical_RDO_local_truth_replaced=False),indent=2)+'\n')


if __name__=='__main__':main()
