#!/usr/bin/env python3
"""Recover completed counter runs and retry only confirmed failed runs."""
import argparse
from collections import Counter
import json
from pathlib import Path
import shutil


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    original = root/'mode_counter_100k_v1'
    jobs = [job for variant in ['beforeGeo', 'afterGeo'] for job in
            json.loads((original/(variant+'_tasks.json')).read_text())]
    pending = [job['output'] for job in jobs if not (Path(job['output'])/'exit_code.txt').exists()]
    assert not pending, f'{len(pending)} original jobs have not returned: {pending[:3]}'
    repo = Path(__file__).resolve().parents[1]
    code = root/'mode_counter_v2_code'
    for folder in ['scripts', 'configs']:
        (code/folder).mkdir(parents=True, exist_ok=False)
    for name in ['run_charge_mode_counter.py', 'recover_charge_mode_counter.py', 'run_slope_task.py',
                 'ntuple_tracking_check.py', 'charge_diagnostics.py', 'validate_charge_audit.py']:
        shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
    shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
    output = root/'mode_counter_recovery_v2'
    output.mkdir(exist_ok=False)
    effective, plans, counts = [], {variant: [] for variant in ['beforeGeo', 'afterGeo']}, Counter()
    for job in jobs:
        folder = Path(job['output'])
        variant, dsid, mode, shard = folder.relative_to(original).parts
        exit_code = int((folder/'exit_code.txt').read_text())
        task = json.loads(Path(job['command'][-1]).read_text())
        destination = output/variant/dsid/mode/shard
        validation = folder/'counter_validation.json'
        if validation.exists() and json.loads(validation.read_text())['all_58_fields_exact']:
            action = 'reuse' if exit_code == 0 else 'validated_before_cleanup_failure'
            data, checked = folder, folder
        else:
            assert exit_code != 0, folder
            arrays = all((folder/(name+'.npz')).is_file() for name in ['events', 'truth', 'tracks', 'segments'])
            if (folder/'summary.json').exists() and arrays:
                action, data, checked = 'validate_completed_reconstruction', folder, destination
                command = ['python3', str(code/'scripts/recover_charge_mode_counter.py'),
                           '--source', str(folder), '--output', str(destination)]
            else:
                action, data, checked = 'retry_failed_reconstruction', destination, destination
                retry = dict(task, runner=str(code/'configs/slope_tracking.py'), output=str(destination),
                             failed_output=str(folder))
                sourcefile = output/(variant+'_'+dsid+'_'+mode+'_'+shard+'.json')
                dump(sourcefile, retry)
                command = ['python3', str(code/'scripts/run_charge_mode_counter.py'), '--task', str(sourcefile)]
            plans[variant].append(dict(kind='analysis', command=command, output=str(destination)))
        counts[action] += 1
        effective.append(dict(variant=variant, dsid=dsid, mode=mode, shard=shard,
            original_counter_output=str(folder), counter_source=str(data), validation_output=str(checked),
            action=action, original_exit_code=exit_code, task=task))
    assert len(effective) == 400
    dump(output/'effective_sources.json', effective)
    dump(output/'manifest.json', dict(actions=dict(counts), original_jobs=400,
        scope='Keep every completed run; validate on scratch; retry only confirmed nonzero, incomplete outputs',
        first_failure='EOS shared-library reads or immediate EOS NPZ readback; original failures retained'))
    for variant, tasks in plans.items():
        dump(output/(variant+'_tasks.json'), tasks)
    print(json.dumps(dict(actions=dict(counts), planned_jobs={v: len(j) for v, j in plans.items()}), indent=2))


if __name__ == '__main__':
    main()
