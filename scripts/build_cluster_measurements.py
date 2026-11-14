#!/usr/bin/env python3
"""Build the actual Calypso cluster package from an isolated audited worktree."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import json
from pathlib import Path
import shlex
import subprocess

from build_slope_diagnostics import component_metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--source-patch', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    base = Path('/eos/home-x/xcheng/FASER/geoRotation_calypso')
    relative = Path('Tracker/TrackerRecTools/FaserSiClusterizationTool')
    package = base/'build'/relative
    diff = subprocess.check_output(['git', '-C', str(args.source), 'diff', 'HEAD'], text=True)
    assert diff == args.source_patch.read_text()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output/'lib').mkdir()
    scratch = Path.cwd()/'cluster_measurements_build'
    scratch.mkdir()
    targets = ['FaserSiClusterizationToolLib', 'FaserSiClusterizationTool']
    objects = {name: [] for name in targets}
    recipes = []
    for row in json.loads((base/'build/compile_commands.json').read_text()):
        if Path(row['directory']) != package:
            continue
        tokens = shlex.split(row['command'].replace(str(base/relative), str(args.source/relative)))
        obj = scratch/row['output']
        obj.parent.mkdir(parents=True, exist_ok=True)
        tokens[tokens.index('-o')+1] = str(obj)
        tokens = [token for token in tokens if token != '-g']
        tokens[1:1] = ['-I'+str(base/'Tracker/TrackerSimEvent'),
            '-I/cvmfs/atlas.cern.ch/repo/sw/software/24.0/Athena/24.0.41/InstallArea/x86_64-el9-gcc13-opt/include']
        target = next(name for name in targets if '/'+name+'.dir/' in row['output'])
        objects[target].append(str(obj))
        recipes.append(dict(command=tokens, source=tokens[tokens.index('-c')+1]))
    assert recipes and all(objects.values())
    (output/'compile_recipes.json').write_text(json.dumps(recipes, indent=2)+'\n')

    def compile_one(row):
        print('Compiling', row['source'], flush=True)
        subprocess.run(row['command'], check=True)

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(compile_one, recipes))
    links = []
    for target in targets:
        tokens = shlex.split((package/f'CMakeFiles/{target}.dir/link.txt').read_text())
        if tokens[0].endswith('atlas_build_run.sh'):
            tokens = tokens[1:]
        tokens = [token for token in tokens if not token.endswith('.o')]
        tokens[tokens.index('-o')+1] = str(output/'lib'/('lib'+target+'.so'))
        tokens[tokens.index('-o'):tokens.index('-o')] = objects[target]
        for i, token in enumerate(tokens):
            if token.startswith('-') or not token.endswith(('.so', '.a')):
                continue
            if token.endswith('libFaserSiClusterizationToolLib.so'):
                tokens[i] = str(output/'lib/libFaserSiClusterizationToolLib.so')
            elif not Path(token).is_absolute():
                tokens[i] = str((package/token).resolve())
        if target.endswith('Lib'):
            tokens.append(str((base/'build/x86_64-el9-gcc13-opt/lib/libTrackerSimEvent.so').resolve()))
        subprocess.run(tokens, check=True)
        links.append(tokens)
    (output/'link_recipes.json').write_text(json.dumps(links, indent=2)+'\n')
    component_metadata(output, 'FaserSiClusterizationTool', 'libFaserSiClusterizationTool')
    (output/'source.patch').write_text(diff)
    manifest = dict(source=str(args.source), source_commit=subprocess.check_output(
        ['git', '-C', str(args.source), 'rev-parse', 'HEAD'], text=True).strip(),
        source_patch=str(args.source_patch), compiled_translation_units=len(recipes),
        jobs=args.jobs, athena='24.0.41', recipe_base=str(package),
        control='Read-only native cluster/RDO/SiHit geometry audit; optional local-X measurement variance scale defaults to one',
        completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (output/'software.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(manifest), flush=True)


if __name__ == '__main__':
    main()
