#!/usr/bin/env python3
"""Build the real Calypso cluster tool with an opt-in truth-position control."""
import argparse
import json
from pathlib import Path
import shlex
import subprocess

from build_slope_diagnostics import component_metadata


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    source=Path('/eos/home-x/xcheng/FASER/geoRotation_calypso')
    build=source/'build'
    output=args.output.resolve();(output/'lib').mkdir(parents=True,exist_ok=True)
    package=build/'Tracker/TrackerRecTools/FaserSiClusterizationTool'
    records=[r for r in json.loads((build/'compile_commands.json').read_text())
             if Path(r['directory'])==package]
    recipes=[]
    objects={k:[] for k in ['FaserSiClusterizationToolLib','FaserSiClusterizationTool']}
    for record in records:
        target=next(k for k in objects if '/'+k+'.dir/' in record['output'])
        obj=Path.cwd()/Path(record['file']).with_suffix('.o').name
        command=shlex.split(record['command'])
        command[command.index('-o')+1]=str(obj)
        command=[x for x in command if x!='-g']
        command.insert(1,'-I'+str(source/'Tracker/TrackerSimEvent'))
        command.insert(1,'-I/cvmfs/atlas.cern.ch/repo/sw/software/24.0/Athena/24.0.41/InstallArea/x86_64-el9-gcc13-opt/include')
        subprocess.run(command,check=True);recipes.append(command);objects[target].append(str(obj))
    links=[]
    for target in objects:
        command=shlex.split((package/f'CMakeFiles/{target}.dir/link.txt').read_text())
        if command[0].endswith('atlas_build_run.sh'):command=command[1:]
        command=[x for x in command if not x.endswith('.o')]
        command[command.index('-o')+1]=str(output/'lib'/('lib'+target+'.so'))
        command[command.index('-o'):command.index('-o')]=objects[target]
        for i,token in enumerate(command):
            if token.startswith('-') or not token.endswith(('.so','.a')):continue
            if 'libFaserSiClusterizationToolLib.so' in token:
                command[i]=str(output/'lib/libFaserSiClusterizationToolLib.so')
            elif not Path(token).is_absolute():command[i]=str((package/token).resolve())
        if target.endswith('Lib'):
            command.append(str((build/'x86_64-el9-gcc13-opt/lib/libTrackerSimEvent.so').resolve()))
        subprocess.run(command,check=True);links.append(command)
    component_metadata(output,'FaserSiClusterizationTool','libFaserSiClusterizationTool')
    scope='Tracker/TrackerRecTools/FaserSiClusterizationTool'
    patch=subprocess.check_output(['git','-C',str(source),'diff','--',scope],text=True)
    (output/'source.patch').write_text(patch)
    (output/'software.json').write_text(json.dumps(dict(source=str(source),scope=scope,
        source_commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip(),
        compile_commands=recipes,link_commands=links,
        control='UseTruthClusterPosition defaults false; only measured local strip coordinate is replaced',
        tracking='Native Calypso clustering, seeding, CKF and refit'),indent=2)+'\n')


if __name__=='__main__':main()
