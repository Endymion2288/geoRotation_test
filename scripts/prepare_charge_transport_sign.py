#!/usr/bin/env python3
"""Freeze endpoint-sign observers and material on/off native transport jobs."""
import argparse
import json
from pathlib import Path
import shutil


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    old = root/'physical_prior_transport_v1'
    common = root/'software_transport_sign_v2'
    assert (common/'software.json').exists()
    output = root/'physical_prior_transport_v2'
    output.mkdir(exist_ok=False)
    code = output/'code'
    for folder in ['configs', 'scripts']:
        (code/folder).mkdir(parents=True)
    shutil.copy2(repo/'configs/charge_prior_transport.py', code/'configs/charge_prior_transport.py')
    shutil.copy2(repo/'scripts/run_slope_transport.py', code/'scripts/run_slope_transport.py')
    previous = json.loads((old/'manifest.json').read_text())
    libraries, jobs, native_sources = {}, {}, {}
    for variant in ['beforeGeo', 'afterGeo']:
        overlay = output/('software_'+variant)
        (overlay/'lib').mkdir(parents=True)
        (overlay/'python').mkdir()
        baseline = Path(previous['libraries'][variant]).parent
        for path in (baseline/'lib').iterdir():
            target = common/'lib'/path.name if 'CalypsoSlopeDiagnostics' in path.name else path.resolve()
            (overlay/'lib'/path.name).symlink_to(target)
        # Freeze native geometry objects that previously came from the release
        # path; the observer is the only changed C++ component.
        run = json.loads((old/'100116'/variant/'run.json').read_text())
        native_sources[variant] = run['loaded_acts_libraries']
        for name in run['loaded_acts_libraries']:
            path = Path(name)
            destination = overlay/'lib'/path.name
            if 'FaserActsGeometry' in path.name and not destination.exists():
                shutil.copy2(path, destination)
        for path in (baseline/'python').iterdir():
            if path.name != 'CalypsoSlopeDiagnostics':
                (overlay/'python'/path.name).symlink_to(path.resolve())
        (overlay/'python/CalypsoSlopeDiagnostics').symlink_to(common/'python/CalypsoSlopeDiagnostics')
        libraries[variant] = str(overlay/'lib')
        jobs[variant] = []
    for dsid in previous['samples']:
        for variant in libraries:
            for control in ['material_on', 'material_off']:
                path = output/dsid/control
                path.mkdir(parents=True, exist_ok=True)
                task = json.loads((old/dsid/(variant+'_source.json')).read_text())
                task.update(runner=str(code/'configs/charge_prior_transport.py'),
                    options=['--disable-material'] if control == 'material_off' else [],
                    output=str(path/variant))
                source = path/(variant+'_source.json')
                dump(source, task)
                jobs[variant].append(dict(kind='analysis', command=['python3',
                    str(code/'scripts/run_slope_transport.py'), '--task', str(source)], output=task['output']))
    for variant, tasks in jobs.items():
        dump(output/(variant+'_tasks.json'), tasks)
    manifest = dict(previous, libraries=libraries, previous=str(old),
        native_library_sources=native_sources, endpoint_observer=str(common/'software.json'),
        controls=['material_on', 'material_off'],
        output_charge='Actual result->charge() and result->parameters()[4], read after native propagation',
        comparison='Material-on existing fields must equal v1; no propagation or input changes')
    dump(output/'manifest.json', manifest)
    print(output)


if __name__ == '__main__':
    main()
