#!/usr/bin/env python3
"""Rebuild native Calypso packages using their existing CMake compile/link recipes."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess


PACKAGES = ('FaserActsGeometry', 'FaserActsKalmanFilter')
TARGETS = ('FaserActsGeometryLib', 'FaserActsGeometry', 'FaserActsKalmanFilter')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--source-patch', type=Path,
        help='Permit exactly this recorded tracked source diff')
    args = parser.parse_args()
    base = Path('/eos/home-x/xcheng/FASER/geoRotation_calypso')
    original_build = base/'build'
    source = args.source.resolve()
    diff = subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True)
    if args.source_patch:
        assert diff == args.source_patch.read_text()
    else:
        assert not diff
    commit = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    lib = output/'lib'
    lib.mkdir(exist_ok=True)
    scratch = Path.cwd()/'native_build'
    scratch.mkdir(exist_ok=False)
    local_lib = scratch/'lib'
    local_lib.mkdir()
    commands = json.loads((original_build/'compile_commands.json').read_text())
    replacements = {str(base/'Tracking/Acts'/p): str(source/'Tracking/Acts'/p)
                    for p in (*PACKAGES, 'FaserActsGeometryInterfaces')}

    def relocate(text):
        for old, new in replacements.items():
            text = text.replace(old, new)
        return text

    jobs = []
    for record in commands:
        if Path(record['directory']).name not in PACKAGES:
            continue
        relative = Path(record['directory']).relative_to(original_build)
        if not any('CMakeFiles/'+target+'.dir/' in record['output'] for target in TARGETS):
            continue
        directory = scratch/relative
        directory.mkdir(parents=True, exist_ok=True)
        command = relocate(record['command'])
        tokens = shlex.split(command)
        target_object = directory/tokens[tokens.index('-o')+1]
        target_object.parent.mkdir(parents=True, exist_ok=True)
        # Debug information has no reconstruction effect and needlessly enlarges EOS products.
        command = command.replace(' -g ', ' ')
        jobs.append(dict(command=command, directory=str(directory), source=relocate(record['file'])))
    assert len(jobs) > 30
    (output/'compile_recipes.json').write_text(json.dumps(jobs, indent=2)+'\n')

    link_recipes = []
    for target in TARGETS:
        package = 'FaserActsGeometry' if target.startswith('FaserActsGeometry') else 'FaserActsKalmanFilter'
        relative = Path('Tracking/Acts')/package
        recipe = original_build/relative/'CMakeFiles'/(target+'.dir')/'link.txt'
        tokens = shlex.split(recipe.read_text())
        if tokens[0].endswith('atlas_build_run.sh'):
            tokens = tokens[1:]
        dest = local_lib/('lib'+target+'.so')
        tokens[tokens.index('-o')+1] = str(dest)
        for i, token in enumerate(tokens):
            if token.startswith('-'):
                continue
            if token.endswith(('.so', '.a')) and not Path(token).is_absolute():
                tokens[i] = str((original_build/relative/token).resolve())
            for dependency in TARGETS:
                name = 'lib'+dependency+'.so'
                if token.endswith('/'+name) and token != str(dest):
                    candidate = local_lib/name
                    tokens[i] = str(candidate)
            if tokens[i].endswith(('.so', '.a')) and Path(tokens[i]).parent != local_lib:
                assert Path(tokens[i]).is_file(), tokens[i]
        tokens += ['-Wl,-rpath,$ORIGIN']
        link_recipes.append(dict(target=target, command=tokens, directory=str(scratch/relative)))
    (output/'link_recipes.json').write_text(json.dumps(link_recipes, indent=2)+'\n')

    def compile_one(item):
        name = Path(item['source']).name
        print('Compiling', name, flush=True)
        status = subprocess.run(item['command'], cwd=item['directory'], shell=True).returncode
        if status:
            raise RuntimeError(f'Native compilation failed: {name}: {status}')
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(compile_one, jobs))
    for recipe in link_recipes:
        subprocess.run(recipe['command'], cwd=recipe['directory'], check=True)
        name = 'lib'+recipe['target']+'.so'
        shutil.copy2(local_lib/name, lib/name)
    for package in PACKAGES:
        dest = output/'python'/package
        shutil.copytree(source/'Tracking/Acts'/package/'python', dest, dirs_exist_ok=True)
        (dest/'__init__.py').touch()
    manifest = dict(source=str(source), commit=commit, native_packages=list(PACKAGES), targets=list(TARGETS),
        source_clean=not diff, source_patch=diff,
        compiler_recipe_base=str(original_build/'compile_commands.json'),
        compiled_translation_units=len(jobs), shared_unchanged_dependencies=str(base/'run'),
        athena='24.0.41', debug_symbols=False, jobs=args.jobs,
        completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (output/'software.json').write_text(json.dumps(manifest, indent=2)+'\n')
    if args.source_patch:
        # Register the additional Gaudi property using Athena's own metadata tools.
        tools = Path('/cvmfs/atlas.cern.ch/repo/sw/software/24.0/AthenaExternals/24.0.41/InstallArea/x86_64-el9-gcc13-opt/bin')
        env = dict(os.environ, LD_LIBRARY_PATH=str(lib)+':'+os.environ.get('LD_LIBRARY_PATH', ''))
        generated = scratch/'genconf'
        generated.mkdir()
        subprocess.run([str(tools/'listcomponents'), '--output',
            str(lib/'libFaserActsKalmanFilter.components'), 'libFaserActsKalmanFilter.so'], env=env, check=True)
        subprocess.run([str(tools/'genconf'), '-o', str(generated), '-p', 'FaserActsKalmanFilter',
            '--no-init', '-i', 'libFaserActsKalmanFilter.so'], env=env, check=True)
        shutil.copy2(generated/'FaserActsKalmanFilter.confdb', lib/'FaserActsKalmanFilter.confdb')
        inputs = generated/'inputs.txt'
        inputs.write_text(str(generated/'FaserActsKalmanFilter.confdb2_part')+'\n')
        merge = Path('/cvmfs/atlas.cern.ch/repo/sw/software/24.0/Athena/24.0.41/InstallArea/x86_64-el9-gcc13-opt/cmake/modules/scripts/mergeConfdb2.py')
        subprocess.run([str(merge), str(lib/'FaserActsKalmanFilter.confdb2'), str(inputs)], env=env, check=True)
        shutil.copy2(generated/'FaserActsKalmanFilterConf.py',
            output/'python/FaserActsKalmanFilter/FaserActsKalmanFilterConf.py')
    print(json.dumps(manifest), flush=True)


if __name__ == '__main__':
    main()
