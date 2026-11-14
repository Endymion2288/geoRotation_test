#!/usr/bin/env python3
"""Replay the frozen charge cohort with corrected material, gating exact output."""
import argparse
from collections import Counter
import json
from pathlib import Path
import shutil


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def prepare(root):
    control = root/'material_encoding_control_v1'
    assert json.loads((control/'full_validation.json').read_text())['all_exact']
    output = control/'replay_v1'
    output.mkdir(exist_ok=False)
    (output/'code').mkdir()
    frozen = output/'code'/Path(__file__).name
    shutil.copy2(__file__, frozen)
    source = root/'charge_replay_v1'
    jobs, smoke, counts = [], [], Counter()
    source_jobs = json.loads((source/'tasks.json').read_text())
    assert len(source_jobs) == 119
    for dsid in ['100116', '100117']:
        old_targets = [r for r in json.loads((source/'targets.json').read_text())
                       if r['dsid'] == dsid and r['wrong']]
        measured = json.loads((control/'analysis_v2/full'/dsid/'replay_events.json').read_text())
        assert {r['sample_entry'] for r in old_targets} == {
            r['sample_entry'] for r in measured if r['corrected']['status'] == 2}
        assert all(r['legacy']['status'] == r['corrected']['status'] == 2 for r in measured)
    for original in source_jobs:
        task_path = Path(original['command'][-1])
        task = json.loads(task_path.read_text())
        dsid = task['dsid']
        old_output = task['output']
        task.update(material=str(control/'corrected_afs_map.json'),
            output=str(output/dsid/Path(old_output).name),
            original_output=str(control/dsid/'full/shards/afterGeo'/
                                f"s{task['sample_entry']//5000:04d}"),
            legacy_replay_output=old_output,
            target_note='Frozen original wrong/normal membership; legacy target q/p values are provenance only.')
        assert '--isolate-smoothing' in task['options']
        path = output/task_path.name
        write_json(path, task)
        command = list(original['command'])
        command[-1] = str(path)
        job = dict(kind='analysis', command=command, output=task['output'])
        jobs.append(job)
        wrong = bool(task['replay_targets'][0]['wrong'])
        counts[dsid+'|'+('wrong' if wrong else 'normal')] += 1
        if wrong and not any(Path(r['output']).parent.name == dsid for r in smoke):
            smoke.append(job)
    assert len(smoke) == 2
    write_json(output/'tasks.json', jobs)
    write_json(output/'smoke_tasks.json', smoke)
    for phase in ['smoke', 'full']:
        write_json(output/(phase+'_validation_tasks.json'), [dict(kind='analysis',
            command=['python3', str(frozen), '--root', str(root), '--phase', phase],
            output=str(output/'validation_jobs'/phase))])
    write_json(output/'manifest.json', dict(
        source=str(source), source_tasks=len(source_jobs), counts=dict(counts),
        library=str(root/'overlay_full_audit_afterGeo/lib'),
        material=str(control/'corrected_afs_map.json'),
        changes=['material path', 'output path', 'exact-reference shard path'],
        invariants='Original input, entry, runner, options, audit library and 119-event cohort; isolated smoothing ON.',
        gate='Require all 58 arrays to equal the corrected full shard for two smoke events before submitting the remaining 117.',
        interpretation='Detailed mechanism replay, not a new population misidentification measurement.'))
    print(json.dumps(dict(output=str(output), counts=dict(counts), smoke=len(smoke))))


def validate(root, phase):
    output = root/'material_encoding_control_v1/replay_v1'
    jobs = json.loads((output/('smoke_tasks.json' if phase == 'smoke' else 'tasks.json')).read_text())
    checks = []
    for job in jobs:
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0', str(folder)
        check = json.loads((folder/'replay_validation.json').read_text())
        assert check['all_reconstruction_arrays_identical'] and check['fields_exact'] == 58, str(folder)
        run = json.loads((folder/'run.json').read_text())
        assert run['audit'] and run['isolate_smoothing'] and run['charge_matrix_audit']
        checks.append(dict(output=str(folder), fields_exact=58,
                           updates=check['updates'], smoothing_recursions=check['smoothing_recursions'],
                           max_update_identity_error_GeVinv=check['max_update_identity_error_GeVinv'],
                           max_smoothing_identity_error_GeVinv=check['max_periodic_smoother_identity_error_GeVinv']))
    write_json(output/(phase+'_validation.json'), dict(all_exact=True, events=len(checks), checks=checks))
    if phase == 'smoke':
        done = {r['output'] for r in jobs}
        remaining = [r for r in json.loads((output/'tasks.json').read_text()) if r['output'] not in done]
        assert len(remaining) == 117
        write_json(output/'remaining_tasks.json', remaining)
    print(json.dumps(dict(phase=phase, all_exact=True, events=len(checks))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['prepare', 'smoke', 'full'], required=True)
    args = parser.parse_args()
    if args.phase == 'prepare':
        prepare(args.root.resolve())
    else:
        validate(args.root.resolve(), args.phase)


if __name__ == '__main__':
    main()
