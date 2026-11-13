#!/usr/bin/env python3
"""Freeze native CKF source-link and simultaneous q/B reversal controls."""
import argparse
import json
from pathlib import Path
import shutil


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['smoke', 'replay'], required=True)
    parser.add_argument('--revision', default='v1')
    parser.add_argument('--observer', type=Path,
                        help='New observer build; existing native tracking build is reused')
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    if args.phase == 'smoke':
        for name in ['hitcontrol_build_v1_job', 'hitcontrol_observer_build_v1_job']:
            assert (root/name/'exit_code.txt').read_text().strip() == '0'
        overlay = root/('overlay_hitcontrol' if args.revision == 'v1'
                        else 'overlay_hitcontrol_'+args.revision)
        overlay.mkdir(exist_ok=False)
        sources = [root/'software_hitcontrol',
                   args.observer.resolve() if args.observer else root/'observer_hitcontrol']
        linked = {}
        for folder in ['lib', 'python']:
            (overlay/folder).mkdir()
            for source in sources:
                for item in (source/folder).iterdir():
                    target = overlay/folder/item.name
                    assert not target.exists()
                    target.symlink_to(item.resolve(), target_is_directory=item.is_dir())
                    linked[str(target.relative_to(overlay))] = str(item.resolve())
        dump(overlay/'software.json', dict(sources=list(map(str, sources)), files=linked,
            control='Actual native CKF source-link filtering and consistent q/B reversal, default OFF'))
        code = root/('hitcontrol_'+args.revision+'_code')
        for folder in ['scripts', 'configs']:
            (code/folder).mkdir(parents=True)
        shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
        for name in ['run_slope_task.py', 'ntuple_tracking_check.py', 'slope_stage_extract.py',
                     'audit_reduce.py', 'charge_audit.py']:
            shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
        output = root/('hitcontrol_smoke_'+args.revision)
        output.mkdir()
        original = []
        for dsid in ['100044', '100116', '100117']:
            task = json.loads((root/'charge_smoke/afterGeo'/dsid/'off/source_task.json').read_text())
            original.append(task)
    else:
        validation = json.loads((root/('hitcontrol_smoke_'+args.revision)/'validation.json').read_text())
        assert validation['all_default_fields_exact'] and validation['all_qB_paths_closed']
        code = root/('hitcontrol_'+args.revision+'_code')
        output = root/('hitcontrol_replay_'+args.revision)
        output.mkdir(exist_ok=False)
        original = [json.loads((Path(job['output'])/'source_task.json').read_text()) for job in
                    json.loads((root/'charge_replay_v1/tasks.json').read_text())]
    tasks, frozen = [], []
    for original_task in original:
        dsid = original_task['dsid']
        modes = ['default', 'truth_filter', 'reverse_qB'] if args.phase == 'smoke' else ['truth_filter', 'reverse_qB']
        for mode in modes:
            options = ['--no-seed-writer', '--actual-covariance', '--isolate-smoothing', '--stage-audit', '--charge-audit']
            if mode == 'truth_filter':
                # The sample has one primary; verify its actual barcode before freezing.
                import numpy as np
                with np.load(Path(original_task['output'])/'truth.npz') as truth:
                    primary = (truth['status'] == 1) & (truth['barcode'] > 0) & (truth['barcode'] < 200000)
                    barcodes = np.unique(truth['barcode'][primary])
                    assert len(barcodes) == 1
                    barcode = int(barcodes[0])
                options.extend(['--truth-source-link-barcode', str(barcode)])
            if mode == 'reverse_qB':
                options.append('--reverse-charge-and-field')
            name = f"f{original_task['file_index']:02d}_e{original_task['skip']:08d}"
            task = dict(original_task, runner=str(code/'configs/slope_tracking.py'), options=options,
                output=str(output/dsid/mode/name), original_output=original_task['output'])
            path = output/(dsid+'_'+mode+'_'+name+'.json')
            dump(path, task)
            tasks.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                                                       '--task', str(path)], output=task['output']))
            frozen.append(task)
    dump(output/'tasks.json', tasks)
    dump(output/'manifest.json', dict(phase=args.phase, controls=['native CKF source-link SDO filter',
        'seed q/p + full covariance cross terms and CKF/refit/extrapolation B reversed together'],
        scope='SDO filter leaves seed means/covariances and CKF algorithm unchanged; subsequent branch/candidate selection may change',
        jobs=len(tasks), source_tasks=frozen))
    print(output/'tasks.json', len(tasks), 'jobs')


if __name__ == '__main__':
    main()
