#!/usr/bin/env python3
"""Retain failed attempts and freeze identical Calypso process-control retries."""
import argparse
import datetime
import json
from pathlib import Path
import shutil

from prepare_cluster_measurements import dump


def prepare(control, output):
    repo = Path(__file__).resolve().parents[1]
    observations = json.loads((control/'failure_observation_v1.json').read_text())
    failed = {row['output'] for row in observations}
    jobs = json.loads((control/'all_tasks.json').read_text())
    assert len(failed) == 16 and len(jobs) == 792
    output.mkdir(exist_ok=False)
    for folder in ['tasks', 'batches']:
        (output/folder).mkdir()
    code = output/'sequence_code_v1'
    code.mkdir()
    for name in ['momentum_task.py', 'run_analysis_sequence.py']:
        shutil.copy2(repo/'scripts'/name, code/name)
    retries, evidence = [], []
    for job in jobs:
        original_folder = Path(job['output'])
        if str(original_folder) not in failed:
            continue
        assert (original_folder/'exit_code.txt').read_text().strip() != '0'
        original = json.loads((original_folder/'source_task.json').read_text())
        assert job['command'][:2] == ['python3', str(control/'code/scripts/run_slope_task.py')]
        frozen = json.loads(Path(job['command'][-1]).read_text())
        assert original == frozen
        destination = output/'reconstruction'/original['dsid']/original_folder.parent.name/original_folder.name
        task = dict(original, output=str(destination))
        assert {k:v for k,v in task.items() if k != 'output'} == {k:v for k,v in original.items() if k != 'output'}
        path = output/'tasks'/Path(job['command'][-1]).name
        dump(path, task)
        retries.append(dict(kind='analysis', command=[*job['command'][:-1], str(path)], output=str(destination)))
        evidence.append(dict(original_output=str(original_folder), recovery_output=str(destination),
                             source_task=str(path), only_output_path_changed=True))
    assert len(retries) == 16
    dump(output/'recovery_jobs.json', retries)
    batches = []
    for index in range(0, len(retries), 2):
        path = output/'batches'/f'batch{index//2:02d}.json'
        dump(path, retries[index:index+2])
        batches.append(dict(kind='analysis', command=['python3', str(code/'run_analysis_sequence.py'),
            '--task-files', str(path), '--worker', str(code/'momentum_task.py'), '--keep-going'],
            output=str(output/'jobs'/f'batch{index//2:02d}')))
    dump(output/'batch_tasks.json', batches)
    dump(output/'manifest.json', dict(time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        original_jobs=str(control/'all_tasks.json'), attempts=evidence,
        scope='Retry only 16 closed failed child attempts in new directories, with frozen runner, libraries, inputs and physics unchanged. Two original running attempts are retained.'))


def resolve(control, output):
    attempts = [output]
    if (control/'recovery_v2/manifest.json').exists():
        attempts.append(control/'recovery_v2')
    replacements, retries = {}, {}
    for attempt in attempts:
        replacements.update({row['original_output']:row['recovery_output']
            for row in json.loads((attempt/'manifest.json').read_text())['attempts']})
        retries.update({row['output']:row for row in json.loads((attempt/'recovery_jobs.json').read_text())})
    resolved = []
    for job in json.loads((control/'all_tasks.json').read_text()):
        row = retries[replacements[job['output']]] if job['output'] in replacements else job
        assert (Path(row['output'])/'exit_code.txt').read_text().strip() == '0', row['output']
        resolved.append(row)
    assert len(resolved) == len({row['output'] for row in resolved}) == 792
    dump(control/'resolved_tasks_v1.json', resolved)


def prepare_stalled(control, output):
    evidence = control/'metadata_stall_observation_v1'
    closure = json.loads((evidence/'removal_closure.json').read_text())
    assert not closure['queue']
    assert {(row['ProcId'], row['JobStatus']) for row in closure['history']} == {(445, 3), (604, 3)}
    originals = json.loads((control/'remaining_tasks.json').read_text())
    output.mkdir(exist_ok=False)
    (output/'tasks').mkdir()
    jobs, attempts = [], []
    for index in [445, 604]:
        job = originals[index]
        old = Path(job['output'])
        assert not (old/'exit_code.txt').exists()
        assert 'Py:MetaReader' in (evidence/str(index)/'tracking.log').read_text()
        original = json.loads(Path(job['command'][-1]).read_text())
        destination = output/'reconstruction'/original['dsid']/old.parent.name/old.name
        task = dict(original, output=str(destination))
        path = output/'tasks'/Path(job['command'][-1]).name
        dump(path, task)
        jobs.append(dict(kind='analysis', command=[*job['command'][:-1], str(path)], output=str(destination)))
        attempts.append(dict(original_output=str(old), recovery_output=str(destination),
            original_job_id=f'1578842.{index}', source_task=str(path), only_output_path_changed=True))
    dump(output/'recovery_jobs.json', jobs)
    dump(output/'manifest.json', dict(attempts=attempts, closure=str(evidence/'removal_closure.json'),
        scope='Identical physics after two documented pre-CKF metadata stalls were explicitly removed; no timeout is treated as an exit.'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['prepare', 'prepare-stalled', 'resolve'], required=True)
    args = parser.parse_args()
    control = args.root.resolve()/'fixed_hit_process_v1'
    output = control/('recovery_v2' if args.phase == 'prepare-stalled' else 'recovery_v1')
    {'prepare':prepare, 'prepare-stalled':prepare_stalled, 'resolve':resolve}[args.phase](control, output)


if __name__ == '__main__':
    main()
