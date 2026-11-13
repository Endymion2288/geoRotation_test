#!/usr/bin/env python3
"""Combine the validated S05 refit chart with S06 native beforeGeo observers."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.root.resolve()
    source = repo.parent/'geoRotation_calypso_charge_refit_before'
    reference = repo.parent/'geoRotation_test_results/stage05_slopes/chart_update_refit_v1/build/software.json'
    recorded = json.loads(reference.read_text())
    assert recorded['commit'] == '1eb885ae22299fb39bb352bead6005d41b0f953e'
    subprocess.run(['git', '-C', str(repo.parent/'geoRotation_calypso'), 'worktree', 'add',
                    '--detach', '--no-checkout', str(source), recorded['commit']], check=True)
    subprocess.run(['git', '-C', str(source), 'sparse-checkout', 'init', '--cone'], check=True)
    subprocess.run(['git', '-C', str(source), 'sparse-checkout', 'set',
        'Tracking/Acts/FaserActsKalmanFilter', 'Tracking/Acts/FaserActsGeometry',
        'Tracking/Acts/FaserActsGeometryInterfaces'], check=True)
    subprocess.run(['git', '-C', str(source), 'read-tree', '-mu', 'HEAD'], check=True)
    subprocess.run(['git', '-C', str(source), 'apply', str(repo/'patches/0022_charge_diagnostics_before.patch')], check=True)
    relative = Path('Tracking/Acts/FaserActsKalmanFilter/src')
    parent = Path(recorded['source'])
    files = ['KalmanFitterTool.cxx', 'KalmanFitterTool.h', 'TrackFitterFunction.cxx',
             'TrackFitterFunction.h', 'RefitDirectionChart.h']
    # These are the actual, already measured S05 native refit components.
    for name in files:
        shutil.copy2(parent/relative/name, source/relative/name)
    subprocess.run(['git', '-C', str(source), 'add', '-N', str(relative/'ChargeDiagnostics.h'),
                    str(relative/'RefitDirectionChart.h')], check=True)
    patch = subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True)
    (root/'refit_control_source.patch').write_text(patch)
    # The refit-only increment is the same five files as S05 patch0017.
    increment = subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD', '--',
                                       *[str(relative/name) for name in files]], text=True)
    (repo/'patches/0027_charge_refit_chart_before.patch').write_text(increment)
    code = root/'refit_control_build_v1_code'
    code.mkdir(exist_ok=False)
    shutil.copy2(repo/'scripts/build_native_calypso.py', code/'build_native_calypso.py')
    tasks = [dict(kind='analysis', command=['python3', str(code/'build_native_calypso.py'),
        '--source', str(source), '--output', str(root/'software_refit_before'), '--jobs', '4',
        '--source-patch', str(root/'refit_control_source.patch'),
        '--reuse-geometry', str(root/'software_beforeGeo')], output=str(root/'refit_control_build_v1_job'))]
    (root/'refit_control_build_v1_tasks.json').write_text(json.dumps(tasks, indent=2)+'\n')
    (root/'refit_control_provenance.json').write_text(json.dumps(dict(source=str(source),
        commit=recorded['commit'], refit_source=str(parent), refit_source_reference=str(reference),
        files=files, baseline='patch0022 native beforeGeo isolated smoothing and read-only charge CSV'), indent=2)+'\n')
    print(root/'refit_control_build_v1_tasks.json')


if __name__ == '__main__':
    main()
