#!/usr/bin/env python3
"""Build an opt-in observer from Calypso source with its existing Athena recipes."""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess


def component_metadata(output, package, library):
    """Use the same listcomponents/genconf/confdb2 pipeline as Athena CMake."""
    lib = output/'lib'
    env = dict(os.environ, LD_LIBRARY_PATH=str(lib)+':'+os.environ.get('LD_LIBRARY_PATH', ''))
    tools = Path('/cvmfs/atlas.cern.ch/repo/sw/software/24.0/AthenaExternals/24.0.41/InstallArea/x86_64-el9-gcc13-opt/bin')
    subprocess.run([str(tools/'listcomponents'), '--output', str(lib/(library+'.components')),
                    library+'.so'], env=env, check=True)
    generated = Path.cwd()/('genConf_'+package)
    generated.mkdir()
    subprocess.run([str(tools/'genconf'), '-o', str(generated), '-p', package,
                    '--no-init', '-i', library+'.so'], env=env, check=True)
    shutil.copy2(generated/(package+'.confdb'), lib/(library+'.confdb'))
    part = generated/(package+'.confdb2_part')
    inputs = generated/'inputs.txt'
    inputs.write_text(str(part)+'\n')
    merge = Path('/cvmfs/atlas.cern.ch/repo/sw/software/24.0/Athena/24.0.41/InstallArea/x86_64-el9-gcc13-opt/cmake/modules/scripts/mergeConfdb2.py')
    subprocess.run([str(merge), str(lib/(package+'.confdb2')), str(inputs)], env=env, check=True)
    py = output/'python'/package
    py.mkdir(parents=True, exist_ok=True)
    shutil.copy2(generated/(package+'Conf.py'), py/(package+'Conf.py'))
    (py/'__init__.py').touch()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = Path('/eos/home-x/xcheng/FASER/geoRotation_calypso')
    build = source/'build'
    output = args.output.resolve()
    lib = output/'lib'
    lib.mkdir(parents=True, exist_ok=True)
    unit = source/'Tracking/Acts/FaserActsKalmanFilter/src/CalypsoSlopeDiagnosticsAlg.cxx'
    compile_commands = json.loads((build/'compile_commands.json').read_text())
    template = next(x for x in compile_commands if x['file'].endswith('/NtupleDumperAlg.cxx'))
    tokens = shlex.split(template['command'])
    tokens[tokens.index('-o')+1] = str(Path.cwd()/'slope_observer.o')
    tokens[tokens.index('-c')+1] = str(unit)
    tokens = [x for x in tokens if x != '-g']
    subprocess.run(tokens, check=True)
    link = shlex.split((Path(template['directory'])/'CMakeFiles/NtupleDumper.dir/link.txt').read_text())
    if link[0].endswith('atlas_build_run.sh'):
        link = link[1:]
    link[link.index('-o')+1] = str(lib/'libCalypsoSlopeDiagnostics.so')
    link = [x for x in link if not x.endswith('.o')]
    link.insert(link.index('-o'), str(Path.cwd()/'slope_observer.o'))
    for i, token in enumerate(link):
        if token.startswith('-Wl,-soname,'):
            link[i] = '-Wl,-soname,libCalypsoSlopeDiagnostics.so'
        elif not token.startswith('-') and token.endswith(('.so', '.a')) and not Path(token).is_absolute():
            link[i] = str((Path(template['directory'])/token).resolve())
    subprocess.run(link, check=True)
    component_metadata(output, 'CalypsoSlopeDiagnostics', 'libCalypsoSlopeDiagnostics')
    shutil.copy2(unit, output/unit.name)
    manifest = dict(source=str(source), unit=str(unit), tracking_changes=False,
                    compiler_recipe=template, compile_command=tokens, link_command=link,
                    component='CalypsoSlopeDiagnosticsAlg', athena='24.0.41',
                    source_commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip())
    manifest['registered_components']=(lib/'libCalypsoSlopeDiagnostics.components').read_text().splitlines()
    (output/'software.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print('Built CalypsoSlopeDiagnosticsAlg', flush=True)


if __name__ == '__main__':
    main()
