#!/usr/bin/env python3
"""Execute a frozen sequence of dependent analysis tasks on one batch worker."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task-files', type=Path, nargs='+', required=True)
    parser.add_argument('--worker', type=Path, required=True)
    parser.add_argument('--keep-going', action='store_true',
                        help='Run independent attempts after a failure; still return a failure status')
    args = parser.parse_args()
    failures = []
    for path in [path.resolve() for path in args.task_files]:
        jobs = json.loads(path.read_text())
        assert jobs and all(job['kind'] == 'analysis' for job in jobs)
        for index in range(len(jobs)):
            print('Running dependent analysis', path, index, flush=True)
            with tempfile.TemporaryDirectory(prefix='analysis-step-') as directory:
                result = subprocess.run([sys.executable, str(args.worker.resolve()), '--tasks', str(path),
                                         '--index', str(index)], cwd=directory)
            if result.returncode:
                failures.append(dict(tasks=str(path), index=index, exit_code=result.returncode))
                if not args.keep_going:
                    raise SystemExit(result.returncode)
    if failures:
        print(json.dumps(dict(failed_independent_attempts=failures)), flush=True)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
