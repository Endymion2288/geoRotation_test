#!/usr/bin/env python3
"""Compile the actual latest Calypso tracking package against its installed dependencies."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--latest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    old = Path('/eos/home-x/xcheng/FASER/geoRotation_calypso')
    latest, source, output = args.latest.resolve(), args.source.resolve(), args.output.resolve()
    package = Path('Tracking/Acts/FaserActsKalmanFilter')
    recipe_file = old/'build/compile_commands.json'
    records = json.loads(recipe_file.read_text())
    template = next(r for r in records if r['file'].endswith('/FaserActsKalmanFilter/src/CKF2.cxx'))
    cmake = (source/package/'CMakeLists.txt').read_text()
    block = cmake.split('atlas_add_component(FaserActsKalmanFilter', 1)[1].split('INCLUDE_DIRS', 1)[0]
    patterns = re.findall(r'^\s*(src/[^\s#]+\.cxx)\s*$', block, re.M)
    # atlas_add_component expands source patterns with file(GLOB); legacy
    # entries that no longer exist in official master do not add sources.
    units = [unit for unit in patterns if (source/package/unit).is_file()]
    assert len(units) >= 12 and 'src/CKF2.cxx' in units
    scratch = Path.cwd()/'latest_native_build'
    scratch.mkdir(exist_ok=False)
    lib = output/'lib'
    lib.mkdir(parents=True, exist_ok=False)
    commands = []
    for i, unit in enumerate(units):
        tokens = shlex.split(template['command'])
        tokens = [t.replace(str(old), str(latest)) for t in tokens if t != '-g']
        tokens[tokens.index('-o')+1] = str(scratch/f'unit_{i}.o')
        tokens[-1] = str(source/package/unit)
        tokens.insert(1, '-I'+str(latest/'run/include'))
        tokens.insert(1, '-I'+str(source/package))
        commands.append(tokens)
    output.mkdir(parents=True, exist_ok=True)
    (output/'compile_recipes.json').write_text(json.dumps(commands, indent=2)+'\n')
    def compile_one(command):
        print('Compiling', command[-1], flush=True)
        subprocess.run(command, check=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(compile_one, commands))
    original = shlex.split((old/'build'/package/'CMakeFiles/FaserActsKalmanFilter.dir/link.txt').read_text())
    if original[0].endswith('atlas_build_run.sh'):
        original = original[1:]
    linked = []
    skip = False
    for token in original:
        if skip:
            skip = False
            continue
        if token == '-o':
            skip = True
            continue
        if token.endswith('.o'):
            continue
        if not token.startswith('-') and token.endswith(('.so', '.a')) and not Path(token).is_absolute():
            token = str((old/'build'/package/token).resolve())
        if token.startswith(str(old)) and token.endswith(('.so', '.a')):
            token = str(latest/'run/lib'/Path(token).name)
            assert Path(token).is_file(), token
        token = token.replace(str(old/'build'), str(latest/'run'))
        linked.append(token)
    linked += ['-o', str(lib/'libFaserActsKalmanFilter.so')]
    # --as-needed resolves only symbols already referenced at each library.
    linked[1:1] = [str(scratch/f'unit_{i}.o') for i in range(len(units))]
    linked += ['-L'+str(latest/'run/lib'), '-lRootWriterToolsLib', '-lTrackingTruthUtilsLib']
    (output/'link_recipe.json').write_text(json.dumps(linked, indent=2)+'\n')
    subprocess.run(linked, check=True)
    python = output/'python/FaserActsKalmanFilter'
    shutil.copytree(source/package/'python', python)
    (python/'__init__.py').touch()
    tools = Path('/cvmfs/atlas.cern.ch/repo/sw/software/24.0/AthenaExternals/24.0.41/InstallArea/x86_64-el9-gcc13-opt/bin')
    env = dict(os.environ, LD_LIBRARY_PATH=str(lib)+':'+os.environ.get('LD_LIBRARY_PATH', ''))
    generated = scratch/'genconf'
    generated.mkdir()
    subprocess.run([str(tools/'listcomponents'), '--output', str(lib/'libFaserActsKalmanFilter.components'),
                    'libFaserActsKalmanFilter.so'], env=env, check=True)
    subprocess.run([str(tools/'genconf'), '-o', str(generated), '-p', 'FaserActsKalmanFilter',
                    '--no-init', '-i', 'libFaserActsKalmanFilter.so'], env=env, check=True)
    shutil.copy2(generated/'FaserActsKalmanFilter.confdb', lib/'FaserActsKalmanFilter.confdb')
    inputs = generated/'inputs.txt'
    inputs.write_text(str(generated/'FaserActsKalmanFilter.confdb2_part')+'\n')
    merge = Path('/cvmfs/atlas.cern.ch/repo/sw/software/24.0/Athena/24.0.41/InstallArea/x86_64-el9-gcc13-opt/cmake/modules/scripts/mergeConfdb2.py')
    subprocess.run([str(merge), str(lib/'FaserActsKalmanFilter.confdb2'), str(inputs)], env=env, check=True)
    shutil.copy2(generated/'FaserActsKalmanFilterConf.py', python/'FaserActsKalmanFilterConf.py')
    manifest = dict(source=str(source), latest=str(latest),
        commit=subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip(),
        source_patch=subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True),
        source_units=units, unmatched_cmake_patterns=[p for p in patterns if p not in units],
        compilation_recipe_template=str(recipe_file),
        shared_dependencies=str(latest/'run'), jobs=args.jobs,
        scope='Full actual latest Calypso FaserActsKalmanFilter sources, native ACTS, installed latest dependencies')
    (output/'software.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    main()
