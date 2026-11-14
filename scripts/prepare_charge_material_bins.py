#!/usr/bin/env python3
"""Freeze native material-bin queries and a row-major rotation control."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys


def dump(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    common = root/'software_material_bins_v4'
    assert (common/'software.json').exists()
    previous = root/'physical_prior_material_trace_v3'
    prior = root/'physical_prior_transport_v2'
    old = json.loads((previous/'manifest.json').read_text())
    output = root/'physical_prior_material_bins_v4'
    output.mkdir(exist_ok=False)
    code = output/'code'
    for folder, names in {'scripts': ['run_slope_transport.py', 'rotate_material.py'],
                          'configs': ['charge_prior_transport.py']}.items():
        (code/folder).mkdir(parents=True)
        for name in names:
            shutil.copy2(repo/folder/name, code/folder/name)
    source_map = root.parents[1]/'stage02_closure/material-current-source.json'
    corrected_map = output/'material-rowmajor-rotated.json'
    subprocess.run([sys.executable, str(code/'scripts/rotate_material.py'),
        '--source', str(source_map), '--output', str(corrected_map),
        '--report', str(output/'rotation_report.json')], check=True)
    libraries = {}
    for variant, library in old['libraries'].items():
        original = Path(library).parent
        overlay = output/('software_'+variant)
        (overlay/'lib').mkdir(parents=True)
        (overlay/'python').mkdir()
        for path in (original/'lib').iterdir():
            target = common/'lib'/path.name if 'CalypsoSlopeDiagnostics' in path.name else path.resolve()
            (overlay/'lib'/path.name).symlink_to(target)
        for path in (original/'python').iterdir():
            target = common/'python'/path.name if path.name == 'CalypsoSlopeDiagnostics' else path.resolve()
            (overlay/'python'/path.name).symlink_to(target)
        libraries[variant] = str(overlay/'lib')
    jobs = {name: [] for name in ['beforeGeo_trace', 'afterGeo_trace',
        'afterGeo_corrected_trace', 'afterGeo_corrected_full']}
    sources = {}
    for dsid in old['samples']:
        for case in jobs:
            variant = 'beforeGeo' if case.startswith('beforeGeo') else 'afterGeo'
            original = (prior/dsid/'material_on'/(variant+'_source.json')
                if case.endswith('full') else previous/dsid/(variant+'_source.json'))
            task = json.loads(original.read_text())
            destination = output/dsid/case
            destination.parent.mkdir(parents=True, exist_ok=True)
            task.update(runner=str(code/'configs/charge_prior_transport.py'),
                options=['--trace-material'] if case.endswith('trace') else [],
                output=str(destination))
            if 'corrected' in case:
                task['material'] = str(corrected_map)
            source = destination.with_name(case+'_source.json')
            dump(source, task)
            jobs[case].append(dict(kind='analysis', command=['python3',
                str(code/'scripts/run_slope_transport.py'), '--task', str(source)],
                output=task['output']))
            sources[dsid+'|'+case] = str(source)
    for case, tasks in jobs.items():
        dump(output/(case+'_tasks.json'), tasks)
    dump(output/'manifest.json', dict(libraries=libraries, source_tasks=sources,
        previous_manifest=str(previous/'manifest.json'),
        prior_manifest=str(prior/'manifest.json'), source_map=str(source_map),
        corrected_map=str(corrected_map), observer=str(common/'software.json'),
        selection=old['samples'],
        scope='Native Calypso material queries on the frozen 44 selected physical points; '
              'corrected-material transport also covers all 7680 common finite physical prior points. '
              'The legacy map and all existing results remain frozen.'))
    print(output)


if __name__ == '__main__':
    main()
