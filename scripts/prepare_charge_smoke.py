#!/usr/bin/env python3
"""Compose complete native libraries and freeze observer ON/OFF validation."""
import argparse
import json
from pathlib import Path
import shutil


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args()
    root=args.root.resolve();repo=Path(__file__).resolve().parents[1]
    for variant in ['beforeGeo','afterGeo']:
        for job in [variant,'observer']:
            assert (root/'build_jobs'/job/'exit_code.txt').read_text().strip()=='0'
        sources=[root/('software_'+variant),root/'observer']
        destination=root/('overlay_'+variant);destination.mkdir(exist_ok=False)
        provenance={}
        for folder in ['lib','python']:
            (destination/folder).mkdir()
            for source in sources:
                for item in sorted((source/folder).iterdir()):
                    target=destination/folder/item.name
                    assert not target.exists()
                    target.symlink_to(item.resolve(),target_is_directory=item.is_dir())
                    provenance[str(target.relative_to(destination))]=str(item.resolve())
        (destination/'software.json').write_text(json.dumps(dict(variant=variant,
            sources=[str(s) for s in sources],files=provenance,isolated_smoothing=True,
            observer_changes_reconstruction=False),indent=2)+'\n')
    code=root/'charge_smoke_code'
    for folder in ['configs','scripts']:(code/folder).mkdir(parents=True,exist_ok=False)
    for name in ['run_slope_task.py','ntuple_tracking_check.py','charge_diagnostics.py']:
        shutil.copy2(repo/'scripts'/name,code/'scripts'/name)
    shutil.copy2(repo/'configs/slope_tracking.py',code/'configs/slope_tracking.py')
    samples=json.loads((repo/'configs/samples.json').read_text())['samples']
    for variant in ['beforeGeo','afterGeo']:
        tasks=[]
        for dsid in ['100044','100116','100117']:
            sample=next(s for s in samples if s['dsid']==dsid)
            for mode in ['on','off']:
                options=['--no-seed-writer','--actual-covariance','--isolate-smoothing']
                if mode=='on':options.append('--charge-diagnostics')
                source=dict(input=sample['files'][0]['path'],skip=0,events=20,dsid=dsid,
                    file_index=0,sample_entry=0,variant=variant,
                    material='/eos/home-x/xcheng/FASER/geoRotation_test_results/stage02_closure/material-current-source.json'
                        if variant=='beforeGeo' else '/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json',
                    runner=str(code/'configs/slope_tracking.py'),options=options,
                    output=str(root/'charge_smoke'/variant/dsid/mode))
                path=root/f'charge_smoke_source_{variant}_{dsid}_{mode}.json'
                path.write_text(json.dumps(source,indent=2)+'\n')
                tasks.append(dict(kind='analysis',command=['python3',str(code/'scripts/run_slope_task.py'),
                    '--task',str(path)],output=source['output']))
        (root/('charge_smoke_'+variant+'_tasks.json')).write_text(json.dumps(tasks,indent=2)+'\n')


if __name__=='__main__':
    main()
