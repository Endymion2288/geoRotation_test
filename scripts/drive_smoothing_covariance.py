#!/usr/bin/env python3
"""Monitor recorded builds, validate smoke runs, submit full Calypso and reduction jobs."""
import argparse
import datetime
import json
from pathlib import Path
import subprocess
import time


def completed(destinations):
    for destination in destinations:
        status = destination/'exit_code.txt'
        if status.exists() and status.read_text().strip() != '0':
            raise RuntimeError(f'Failed task: {destination}')
    return all((p/'exit_code.txt').exists() for p in destinations)


def wait_for(destinations, label):
    while not completed(destinations):
        done = sum((p/'exit_code.txt').exists() for p in destinations)
        print(f'{datetime.datetime.now(datetime.timezone.utc).isoformat()} {label}: {done}/{len(destinations)} complete', flush=True)
        time.sleep(30)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--tag', default='20261112')
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]

    def submit(tasks, library, name, memory=3500):
        worker = root/(name+'_worker')
        if not (worker/'submission.json').exists():
            subprocess.run(['python3', str(repo/'scripts/submit_momentum.py'),
                '--tasks', str(tasks), '--library', str(library), '--name', name,
                '--memory-mb', str(memory), '--disk-gb', '12', '--flavour', 'workday'], check=True)
        return json.loads((worker/'submission.json').read_text())

    wait_for([root/'build_before', root/'build_after'], 'builds')
    cases = ['beforeGeo', 'afterGeo', 'afterGeo_originalVariance']
    submissions = {}
    for case in cases:
        software = root/('software_before' if case == 'beforeGeo' else 'software_after')
        assert (software/'lib/FaserActsKalmanFilter.confdb2').is_file()
        submissions[case+'_smoke'] = submit(root/(case+'_smoke_tasks.json'), software/'lib',
            f's04_smcov_{case}_smoke_{args.tag}')
    wait_for([root/'smoke'/case for case in cases], 'smoke')
    expected = json.loads((root/'inputs.json').read_text())['cases']
    for case in cases:
        run = json.loads((root/'smoke'/case/'run.json').read_text())
        assert run['success'] and run['isolate_smoothing']
        assert not run['audit'] and not run['reset_event_maps']
        assert run['args']['initial_variance_inflation'] == expected[case]['initial_variance_inflation']
        software = root/expected[case]['software']
        assert str(software/'lib/libFaserActsKalmanFilter.so') in run['loaded_acts_libraries']
        for key, factor, sigma in zip(['covLoc0', 'covLoc1', 'covPhi', 'covTheta', 'covQOverP'],
                expected[case]['initial_variance_inflation'], [.9, .019, .033, .00029, .0001]):
            assert abs(run['seed_covariance'][key]/(factor*sigma*sigma)-1) < 1e-12
        submissions[case] = submit(root/(case+'_tasks.json'), software/'lib',
            f's04_smcov_{case}_full_{args.tag}')
    (root/'reconstruction_submissions.json').write_text(json.dumps(submissions, indent=2)+'\n')
    destinations = [Path(task['output']) for case in cases
        for task in json.loads((root/(case+'_tasks.json')).read_text())]
    wait_for(destinations, 'full reconstruction')
    tasks = [dict(kind='analysis', command=['python3',
        str(root/'code/analysis/smoothing_covariance.py'), '--root', str(root)],
        output=str(root/'analysis_job'))]
    from shutil import copytree
    copytree(repo/'analysis', root/'code/analysis', dirs_exist_ok=True,
        ignore=lambda directory, names: [n for n in names if n == '__pycache__'])
    taskfile = root/'analysis_tasks.json'
    taskfile.write_text(json.dumps(tasks, indent=2)+'\n')
    submissions['analysis'] = submit(taskfile, root/'software_after/lib',
        f's04_smcov_analysis_{args.tag}', memory=6000)
    wait_for([root/'analysis_job'], 'analysis')
    (root/'all_submissions.json').write_text(json.dumps(submissions, indent=2)+'\n')
    print('All reconstruction and figure tasks completed.', flush=True)


if __name__ == '__main__':
    main()
