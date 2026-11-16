#!/usr/bin/env python3
"""Compile a Calypso observer with the tracking package's exact Eigen ABI."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    task = json.loads(args.task.read_text())
    output = Path(task['build_output'])
    output.mkdir(parents=True, exist_ok=False)
    source = output/'LatestSeedFieldClosure.cxx'
    shutil.copy2(task['helper_source'], source)
    command = json.loads(Path(task['compile_recipes']).read_text())[0]
    command.insert(1, '-I'+str(Path(task['helper_source']).parent))
    command[command.index('-o')+1] = str(Path.cwd()/'seed_field_closure.o')
    command[-1] = str(source)
    subprocess.run(command, check=True)
    link = json.loads(Path(task['link_recipe']).read_text())
    link = [arg for arg in link if not arg.endswith('.o')]
    link[link.index('-o')+1] = task['helper_library']
    link[1:1] = [str(Path.cwd()/'seed_field_closure.o')]
    link += ['-L'+str(Path(task['compile_recipes']).parent/'lib'), '-lFaserActsKalmanFilter']
    subprocess.run(link, check=True)
    (output/'build.json').write_text(json.dumps(dict(
        source=str(source), compile=command, link=link,
        scope='Observer using native Calypso ToolHandle/extrapolator and CircleFit; no new tracking implementation'),
        indent=2)+'\n')
    subprocess.run(['python3', task['runner'], '--task', str(args.task)], check=True)


if __name__ == '__main__':
    main()
