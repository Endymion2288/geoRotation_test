#!/usr/bin/env python3
"""Recover terminal joint-regression failures without changing frozen jobs."""
import argparse
import json
from pathlib import Path
import subprocess


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def prepare(root, cluster, schedd):
    source = root/'joint_regression_full_v1'
    output = root/'joint_regression_recovery_v1'
    for dsid in ['100116', '100117']:
        check = json.loads((root/'joint_environment_validation_v1'/dsid/'validation.json').read_text())
        assert check['all_58_fields_exact']
    queue = subprocess.check_output(['condor_q', '-name', schedd, '-constraint',
        f'ClusterId == {cluster}', '-json'], text=True)
    live = {row['ProcId'] for row in json.loads(queue or '[]')}
    original_jobs = json.loads((source/'tasks.json').read_text())
    failures = []
    for proc, job in enumerate(original_jobs):
        code = Path(job['output'])/'exit_code.txt'
        if code.exists() and code.read_text().strip() != '0':
            assert proc not in live, f'{cluster}.{proc} is still queued'
            failures.append((proc, job, int(code.read_text())))
    assert failures, 'No terminal failures to recover'
    output.mkdir(exist_ok=False)
    (output/'original_queue.json').write_text(queue)
    jobs, records = [], []
    for proc, job, code in failures:
        command = job['command']
        path = Path(command[command.index('--task')+1])
        task = json.loads(path.read_text())
        assert task['output'] == job['output']
        replacement = dict(task, output=str(output/'jobs'/f'p{proc:04d}'))
        replacement_path = output/f'p{proc:04d}_source.json'
        dump(replacement_path, replacement)
        new_command = command.copy()
        new_command[new_command.index('--task')+1] = str(replacement_path)
        recovery = dict(job, command=new_command, output=replacement['output'])
        jobs.append(recovery)
        records.append(dict(original_proc=proc, original_exit_code=code,
            original_job=job, original_source=str(path), recovery_job=recovery))
    dump(output/'tasks.json', jobs)
    dump(output/'manifest.json', dict(original_cluster=cluster, original_schedd=schedd,
        records=records, original_source=str(source),
        change='Explicit analysis dependency path in worker; reconstruction input/config/libraries unchanged',
        environment_validation='Two actual 20-event Calypso runs; all 58 output fields exact'))
    print(output/'tasks.json', len(jobs), 'terminal failures')


def effective(root):
    recovery = root/'joint_regression_recovery_v1'
    records = json.loads((recovery/'manifest.json').read_text())['records']
    replacements = {r['original_job']['output']: r for r in records}
    source = root/'joint_regression_full_v1'
    output = root/'joint_regression_effective_v1'
    manifest = json.loads((source/'manifest.json').read_text())
    groups, coverage = [], []
    # Check all sources before creating the effective view.
    for group in manifest['groups']:
        jobs = json.loads((Path(group['group'])/'beforeGeo_tasks.json').read_text())
        for job in jobs:
            replacement = replacements.get(job['output'])
            current = replacement['recovery_job'] if replacement else job
            folder = Path(current['output'])
            assert (folder/'exit_code.txt').read_text().strip() == '0', str(folder)
            original_task = json.loads((Path(job['command'][job['command'].index('--task')+1])).read_text())
            actual_task = json.loads((folder/'source_task.json').read_text())
            assert {k: v for k, v in original_task.items() if k != 'output'} == \
                   {k: v for k, v in actual_task.items() if k != 'output'}, str(folder)
        groups.append((group, jobs))
    output.mkdir(exist_ok=False)
    effective_groups = []
    for group, jobs in groups:
        path = output/group['dsid']/group['mode']/group['control']
        path.mkdir(parents=True)
        updated = [replacements[j['output']]['recovery_job'] if j['output'] in replacements else j for j in jobs]
        dump(path/'beforeGeo_tasks.json', updated)
        dump(path/'manifest.json', json.loads((Path(group['group'])/'manifest.json').read_text()))
        effective_groups.append(dict(group, group=str(path)))
        coverage.append(dict(dsid=group['dsid'], mode=group['mode'], control=group['control'],
            jobs=len(jobs), recovered=sum(j['output'] in replacements for j in jobs)))
    manifest.update(groups=effective_groups, original_manifest=str(source/'manifest.json'),
        recovery_manifest=str(recovery/'manifest.json'), effective_coverage=coverage)
    dump(output/'manifest.json', manifest)
    print(output, len(effective_groups), 'complete groups')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['prepare', 'effective'], required=True)
    parser.add_argument('--cluster', type=int, default=1578686)
    parser.add_argument('--schedd', default='bigbird20.cern.ch')
    args = parser.parse_args()
    if args.phase == 'prepare':
        prepare(args.root.resolve(), args.cluster, args.schedd)
    else:
        effective(args.root.resolve())


if __name__ == '__main__':
    main()
