#!/usr/bin/env python3
"""Build the real Calypso conditions package with its recorded CMake recipes."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--source-patch', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    base = Path('/eos/home-x/xcheng/FASER/geoRotation_calypso')
    relative = Path('Tracker/TrackerConditions/FaserSCT_ConditionsAlgorithms')
    package = relative.name
    build = base/'build'
    diff = subprocess.check_output(['git', '-C', str(args.source), 'diff', 'HEAD'], text=True)
    assert diff == args.source_patch.read_text()
    commit = subprocess.check_output(['git', '-C', str(args.source), 'rev-parse', 'HEAD'], text=True).strip()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    lib = output/'lib'
    lib.mkdir()
    scratch = Path.cwd()/'alignment_build'
    directory = scratch/relative
    directory.mkdir(parents=True)
    local_lib = scratch/'lib'
    local_lib.mkdir()
    jobs = []
    for row in json.loads((build/'compile_commands.json').read_text()):
        if Path(row['directory']) != build/relative:
            continue
        command = row['command'].replace(str(base/relative), str(args.source/relative)).replace(' -g ', ' ')
        tokens = shlex.split(command)
        (directory/tokens[tokens.index('-o')+1]).parent.mkdir(parents=True, exist_ok=True)
        jobs.append(dict(command=command, directory=str(directory), source=row['file'].replace(str(base/relative), str(args.source/relative))))
    assert len(jobs) == 11
    (output/'compile_recipes.json').write_text(json.dumps(jobs, indent=2)+'\n')
    def compile_one(row):
        print('Compiling', row['source'], flush=True)
        subprocess.run(row['command'], cwd=row['directory'], shell=True, check=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(compile_one, jobs))
    recipe = build/relative/'CMakeFiles'/(package+'.dir')/'link.txt'
    tokens = shlex.split(recipe.read_text())
    if tokens[0].endswith('atlas_build_run.sh'):
        tokens = tokens[1:]
    target = local_lib/('lib'+package+'.so')
    tokens[tokens.index('-o')+1] = str(target)
    for i, token in enumerate(tokens):
        if not token.startswith('-') and token.endswith(('.so', '.a')) and not Path(token).is_absolute():
            tokens[i] = str((build/relative/token).resolve())
    (output/'link_recipe.json').write_text(json.dumps(tokens, indent=2)+'\n')
    subprocess.run(tokens, cwd=directory, check=True)
    shutil.copy2(target, lib/target.name)
    python = output/'python'/package
    python.mkdir(parents=True)
    (python/'__init__.py').touch()
    tools = Path('/cvmfs/atlas.cern.ch/repo/sw/software/24.0/AthenaExternals/24.0.41/InstallArea/x86_64-el9-gcc13-opt/bin')
    generated = scratch/'genconf'
    generated.mkdir()
    env = dict(os.environ, LD_LIBRARY_PATH=str(lib)+':'+os.environ.get('LD_LIBRARY_PATH', ''))
    library = 'lib'+package
    subprocess.run([str(tools/'listcomponents'), '--output', str(lib/(library+'.components')), library+'.so'], env=env, check=True)
    subprocess.run([str(tools/'genconf'), '-o', str(generated), '-p', package, '--no-init', '-i', library+'.so'], env=env, check=True)
    shutil.copy2(generated/(package+'.confdb'), lib/(package+'.confdb'))
    merge = Path('/cvmfs/atlas.cern.ch/repo/sw/software/24.0/Athena/24.0.41/InstallArea/x86_64-el9-gcc13-opt/cmake/modules/scripts/mergeConfdb2.py')
    inputs = generated/(package+'_inputs.txt')
    inputs.write_text(str(generated/(package+'.confdb2_part'))+'\n')
    subprocess.run([str(merge), str(lib/(package+'.confdb2')), str(inputs)], env=env, check=True)
    shutil.copy2(generated/(package+'Conf.py'), python/(package+'Conf.py'))
    manifest = dict(source=str(args.source), commit=commit, source_patch=diff,
        compiled_translation_units=len(jobs), compiler_recipe_base=str(build/'compile_commands.json'),
        actual_package=package, shared_dependencies=str(base/'run'), jobs=args.jobs,
        athena='24.0.41', completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (output/'software.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(manifest), flush=True)


if __name__ == '__main__':
    main()
