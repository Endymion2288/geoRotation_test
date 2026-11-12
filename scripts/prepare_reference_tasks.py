#!/usr/bin/env python3
"""Freeze native geometry comparisons and the GitHub-macro truth control."""
import argparse
import datetime
import json
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inputs', type=Path,
        default=Path(__file__).resolve().parents[1]/'configs/reference_110051.json',
        help='Previously measured full-input manifest (CollectionTree entries per file)')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(args.inputs.read_text())
    sample, = [s for s in manifest['samples'] if s['dsid'] == '110051']
    files = sample['files']
    actual = sorted(Path(files[0]['path']).parent.glob('*.root'))
    assert [str(p) for p in actual] == [f['path'] for f in files]
    assert all(p.stat().st_size == f['size_bytes'] for p, f in zip(actual, files))
    for relative in ['configs/reference_tracking.py', 'scripts/reference_truth_observer.py']:
        dest = root/'code'/relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, dest)
    closure = root.parent/'stage02_closure'
    cases = dict(
        beforeGeo=dict(variant='beforeGeo', material=str(closure/'material-current-source.json'),
            software='software_before', commit='1eb885ae22299fb39bb352bead6005d41b0f953e'),
        afterGeo=dict(variant='afterGeo',
            material='/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json',
            software='software_after', commit='88c457c1e4bf4154e67b96114ee8e7020d9b1a97'),
        afterGeo_strict=dict(variant='afterGeo', material=str(closure/'material-current-rotated.json'),
            software='software_after', commit='88c457c1e4bf4154e67b96114ee8e7020d9b1a97'))
    for case, settings in cases.items():
        tasks = []
        for index, source in enumerate(files):
            tasks.append(dict(kind='reference', dsid='110051', case=case,
                variant=settings['variant'], file_index=index, skip=0, events=source['entries'],
                input=source['path'], material=settings['material'], mask='0,1,2,3,4,5',
                runner=str(root/'code/configs/reference_tracking.py'),
                macro=str(root/'reference/adapted/ckf_residual_perLayer.C'),
                output=str(root/'shards'/case/f'f{index:02d}')))
        (root/(case+'_tasks.json')).write_text(json.dumps(tasks,indent=2)+'\n')
    first = json.loads((root/'afterGeo_tasks.json').read_text())[0]
    (root/'afterGeo_f00_tasks.json').write_text(json.dumps([first],indent=2)+'\n')
    controls=[]
    for case,mask in [('afterGeo_mask012','0,1,2'),('afterGeo_4station','')]:
        task=dict(first,case=case,mask=mask,output=str(root/'shards'/case/'f00'))
        controls.append(task)
        (root/(case+'_tasks.json')).write_text(json.dumps([task],indent=2)+'\n')
    (root/'configuration_controls_tasks.json').write_text(json.dumps(controls,indent=2)+'\n')
    manifest = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        samples=[sample], previous_entry_measurement=str(args.inputs.resolve()),
        directory_and_file_sizes_rechecked=True, source_events=sum(f['entries'] for f in files),
        cases=cases, scope='Full 110051 sample; overall native Calypso before/after geometry rotation',
        source_macro_commit='33ca2281b961935456dfd4dc9fc06235e24b551e',
        cuts='chi2Sum/NDF<20 && rMax_smt<90 && NDF>9',
        residual='(1/abs(eQOP_stage)-1/abs(truth_qop_MeV_inv)/1000)/(1/abs(truth_qop_MeV_inv)/1000)',
        binning=dict(bins=100, low=-1.2, high=1.2, normalization='raw Entries'),
        truth_branches=['native','production','common_native','common_production'],
        audit=False, reset_event_maps=False, isolate_smoothing=False,
        mask=[0,1,2,3,4,5], geometry='FASERNU-04', conditions='OFLCOND-FASER-05')
    (root/'inputs.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(files=len(files), events=manifest['source_events'], cases=list(cases))))


if __name__ == '__main__':
    main()
