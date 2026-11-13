#!/usr/bin/env python3
"""Reuse the real Calypso SiHit-centroid tool for controlled charge replays."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['smoke', 'replay'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    code = root/'centroid_control_v1_code'
    if args.phase == 'smoke':
        overlay = root/'overlay_centroid_v1'
        overlay.mkdir(exist_ok=False)
        sources = [root/'overlay_hitcontrol_v2', root.parent.parent/'stage05_slopes/software_cluster_control_v2']
        links = {}
        for folder in ['lib', 'python']:
            (overlay/folder).mkdir()
            for source in sources:
                for item in (source/folder).iterdir():
                    target = overlay/folder/item.name
                    assert not target.exists()
                    target.symlink_to(item.resolve(), target_is_directory=item.is_dir())
                    links[str(target.relative_to(overlay))] = str(item.resolve())
        dump(overlay/'software.json', dict(sources=list(map(str, sources)), files=links,
            control='Actual Calypso UseTruthClusterPosition; preserve widths/errors/noise/cluster IDs and unmatched-truth fallback'))
        for folder in ['scripts', 'configs']:
            (code/folder).mkdir(parents=True, exist_ok=False)
        shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
        for name in ['run_slope_task.py', 'ntuple_tracking_check.py', 'slope_stage_extract.py',
                     'audit_reduce.py', 'charge_audit.py']:
            shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
        output = root/'centroid_smoke_v1'
        originals = [json.loads((root/'charge_smoke/afterGeo'/dsid/'off/source_task.json').read_text())
                     for dsid in ['100044', '100116', '100117']]
        modes = ['default', 'centroid']
    else:
        valid = json.loads((root/'centroid_smoke_v1/validation.json').read_text())
        assert valid['all_default_58_fields_exact'] and valid['all_actual_cluster_libraries_match']
        output = root/'centroid_replay_v1'
        originals = [json.loads((Path(job['output'])/'source_task.json').read_text()) for job in
                     json.loads((root/'charge_replay_v1/tasks.json').read_text())]
        modes = ['centroid', 'centroid_primary_filter']
    output.mkdir(exist_ok=False)
    jobs, sources = [], []
    for original in originals:
        for mode in modes:
            options = ['--no-seed-writer', '--actual-covariance', '--isolate-smoothing', '--stage-audit', '--charge-audit']
            if mode != 'default':
                options.append('--truth-cluster-position')
            if mode == 'centroid_primary_filter':
                with np.load(Path(original['output'])/'truth.npz') as truth:
                    primary = (truth['status'] == 1) & (truth['barcode'] > 0) & (truth['barcode'] < 200000)
                    barcodes = np.unique(truth['barcode'][primary])
                    assert len(barcodes) == 1
                options += ['--truth-source-link-barcode', str(int(barcodes[0]))]
            name = f"f{original['file_index']:02d}_e{original['skip']:08d}"
            task = dict(original, runner=str(code/'configs/slope_tracking.py'), options=options,
                original_output=original['output'], output=str(output/original['dsid']/mode/name))
            source = output/(original['dsid']+'_'+mode+'_'+name+'.json')
            dump(source, task)
            jobs.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                                                       '--task', str(source)], output=task['output']))
            sources.append(task)
    dump(output/'tasks.json', jobs)
    dump(output/'manifest.json', dict(phase=args.phase, jobs=len(jobs), source_tasks=sources,
        scope='Native Calypso clustering, segment/seed finding, CKF, refit; isolated smoothing ON',
        control='Only measured strip coordinate replaced by SDO-majority energy-weighted same-wafer SiHit centroid; seed and branch responses may change'))
    print(output/'tasks.json', len(jobs), 'jobs')


if __name__ == '__main__':
    main()
