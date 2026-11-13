#!/usr/bin/env python3
"""Check the actual Calypso direction-chart header on real writer angles."""
import argparse
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--recipes', type=Path, required=True)
    parser.add_argument('--header-directory', type=Path, required=True)
    parser.add_argument('--angles', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    driver = Path(__file__).resolve().parents[1] / 'configs/refit_direction_chart_check.cxx'
    records = json.loads(args.recipes.read_text())
    recipe = next(r for r in records if r['source'].endswith('TrackFitterFunction.cxx'))
    command = [token for token in shlex.split(recipe['command']) if token != '-c']
    command[-1] = str(driver)
    command.append('-I' + str(args.header_directory))
    compiler_lib = Path(command[0]).parent.parent / 'lib64'
    env = dict(os.environ, LD_LIBRARY_PATH=str(compiler_lib) + ':' + os.environ.get('LD_LIBRARY_PATH', ''))
    with tempfile.TemporaryDirectory(prefix='s05-chart-math-') as scratch:
        executable = Path(scratch) / 'check'
        command[command.index('-o') + 1] = str(executable)
        (args.output / 'compile_command.json').write_text(json.dumps(command, indent=2) + '\n')
        with (args.output / 'compile.log').open('w') as log:
            status = subprocess.run(command, cwd=scratch, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
        assert status == 0, str(args.output / 'compile.log')
        rows = []
        with np.load(args.angles) as angles:
            for stage in ['prt', 'flt', 'smt']:
                pair = np.column_stack([angles['ePHI_' + stage], angles['eTHETA_' + stage]])
                rows.append(pair[np.isfinite(pair).all(axis=1)])
        stream = io.StringIO()
        np.savetxt(stream, np.concatenate(rows), fmt='%.17g')
        result = subprocess.run([str(executable)], input=stream.getvalue(), capture_output=True, text=True, env=env)
        measured = json.loads(result.stdout)
        measured.update(actual_header=str(args.header_directory / 'RefitDirectionChart.h'),
                        input=str(args.angles), compiler_runtime=str(compiler_lib),
                        passed=result.returncode == 0,
                        scope='Exact Calypso header: numerical derivative, direction roundtrip and inverse Jacobian; no tracking')
        (args.output / 'validation.json').write_text(json.dumps(measured, indent=2) + '\n')
        print(json.dumps(measured))
        assert result.returncode == 0


if __name__ == '__main__':
    main()
