#!/usr/bin/env python3
"""Build the actual Calypso NtupleDumper using its existing Athena recipes."""
import argparse
import json
from pathlib import Path
import shlex
import shutil
import subprocess

from build_slope_diagnostics import component_metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path, default=Path('/eos/home-x/xcheng/FASER/geoRotation_calypso_main'))
    args = parser.parse_args()
    rotated = Path('/eos/home-x/xcheng/FASER/geoRotation_calypso')
    master = args.source.resolve()
    output = args.output.resolve()
    lib = output/'lib'
    lib.mkdir(parents=True, exist_ok=True)
    build = rotated/'build'
    source_dir = rotated/'PhysicsAnalysis/NtupleDumper'
    target_dir = master/'PhysicsAnalysis/NtupleDumper'
    records = [r for r in json.loads((build/'compile_commands.json').read_text())
               if Path(r['directory']).name == 'NtupleDumper' and 'NtupleDumper.dir/' in r['output']]
    objects = []
    recipes = []
    for record in records:
        tokens = shlex.split(record['command'].replace(str(source_dir), str(target_dir)))
        obj = Path.cwd()/Path(record['file']).with_suffix('.o').name
        objects.append(str(obj))
        tokens[tokens.index('-o')+1] = str(obj)
        tokens = [x for x in tokens if x != '-g']
        recipes.append(tokens)
        subprocess.run(tokens, check=True)
    package_build = build/'PhysicsAnalysis/NtupleDumper'
    link = shlex.split((package_build/'CMakeFiles/NtupleDumper.dir/link.txt').read_text())
    if link[0].endswith('atlas_build_run.sh'):
        link = link[1:]
    link = [x for x in link if not x.endswith('.o')]
    link[link.index('-o')+1] = str(lib/'libNtupleDumper.so')
    link[link.index('-o'):link.index('-o')] = objects
    for i, token in enumerate(link):
        if not token.startswith('-') and token.endswith(('.so', '.a')) and not Path(token).is_absolute():
            link[i] = str((package_build/token).resolve())
    subprocess.run(link, check=True)
    py = output/'python/NtupleDumper'
    shutil.copytree(target_dir/'python', py, dirs_exist_ok=True)
    (py/'__init__.py').touch()
    component_metadata(output, 'NtupleDumper', 'libNtupleDumper')
    (output/'software.json').write_text(json.dumps(dict(source=str(master),
        native_package='PhysicsAnalysis/NtupleDumper', source_changes=True,
        tracking_changes=False,
        angle_control_available='RotatedPhysicalAngles' in (target_dir/'src/NtupleDumperAlg.h').read_text(),
        compatibility_change='Source diff is archived; master waveform accessor compatibility is separate from the opt-in rotated angle control',
        source_commit=subprocess.check_output(['git','-C',str(master),'rev-parse','HEAD'],text=True).strip(),
        compile_commands=recipes, link_command=link), indent=2)+'\n')
    (output/'source.patch').write_text(subprocess.check_output(['git','-C',str(master),'diff','--','PhysicsAnalysis/NtupleDumper'],text=True))


if __name__ == '__main__':
    main()
