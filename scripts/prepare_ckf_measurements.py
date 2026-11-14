#!/usr/bin/env python3
"""Freeze CKF measurement controls with actual native Calypso covariance."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess

from prepare_cluster_measurements import dump, sources, override


def task(root, destination, original, mode):
    options = ['--no-seed-writer', '--actual-covariance', '--isolate-smoothing',
               '--stage-audit', '--charge-audit', '--cluster-coordinate-audit']
    if mode != 'native':
        seedfile = destination/'seeds'/(original['dsid']+f"_{original['sample_entry']:08d}.txt")
        if not seedfile.exists():
            override(Path(original['output']), seedfile)
        options += ['--seed-override-file', str(seedfile), '--seed-override-mode', 'mean']
    if mode == 'centroid_fixed':
        options += ['--truth-cluster-position']
    elif mode == 'ckf_quarter':
        options += ['--ckf-measurement-variance-scale', '.25']
    elif mode == 'ckf_four':
        options += ['--ckf-measurement-variance-scale', '4']
    elif mode == 'ckf_cluster':
        options += ['--ckf-cluster-measurement-variance']
    output = destination/'reconstruction'/original['dsid']/mode/f"e{original['sample_entry']:08d}"
    current = dict(original, runner=str(destination/'code/configs/slope_tracking.py'),
        original_output=original['output'], output=str(output), options=options, coordinate_mode=mode)
    source = destination/'tasks'/(original['dsid']+'_'+mode+f"_{original['sample_entry']:08d}.json")
    dump(source, current)
    runner = 'run_charge_replay.py' if mode == 'native' else 'run_slope_task.py'
    return dict(kind='analysis', command=['python3', str(destination/'code/scripts'/runner),
                '--task', str(source)], output=str(output))


def prepare_build(root, destination, source):
    repo = Path(__file__).resolve().parents[1]
    destination.mkdir(exist_ok=False)
    for folder in ['tasks', 'seeds', 'code/scripts', 'code/configs']:
        (destination/folder).mkdir(parents=True)
    for name in ['build_native_calypso.py', 'run_slope_task.py', 'run_charge_replay.py',
                 'validate_charge_audit.py', 'audit_reduce.py', 'charge_audit.py',
                 'slope_stage_extract.py', 'ntuple_tracking_check.py',
                 'prepare_cluster_measurements.py', 'prepare_ckf_measurements.py', 'cluster_coordinate_io.py']:
        shutil.copy2(repo/'scripts'/name, destination/'code/scripts'/name)
    shutil.copy2(repo/'configs/slope_tracking.py', destination/'code/configs/slope_tracking.py')
    diff = subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True)
    assert 'UseClusterMeasurementVariance' in diff and 'MeasurementVarianceScale' in diff
    (destination/'source.patch').write_text(diff)
    dump(destination/'build_tasks.json', [dict(kind='analysis', command=['python3',
        str(destination/'code/scripts/build_native_calypso.py'), '--source', str(source),
        '--source-patch', str(destination/'source.patch'), '--output', str(destination/'software'),
        '--jobs', '4'], output=str(destination/'jobs/build'))])
    dump(destination/'manifest.json', dict(source=str(source), source_commit='88c457c1e4bf4154e67b96114ee8e7020d9b1a97',
        cluster_control=str(root/'cluster_measurement_v1/software'),
        material=str(root/'material_encoding_control_v1/corrected_afs_map.json'),
        scope='Six frozen original purity-one wrong events plus 60 original normal controls; isolated smoothing ON',
        native_measurement_variance_mm2=.0231**2,
        controls={'native':'Exact corrected-material default; read-only coordinate audit',
            'centroid_fixed':'Same-wafer energy-weighted SiHit midpoint, original six seed means and full covariance',
            'ckf_quarter':'Actual CKF measurement V multiplied by 0.25; original full seed and cluster covariance',
            'ckf_four':'Actual CKF measurement V multiplied by 4; original full seed and cluster covariance',
            'ckf_cluster':'Actual CKF uses the existing cluster local-X variance; original full seed'},
        gates='Default must exactly reproduce all 58 corrected-reference arrays. Audit library and full seed identity checked; selected hits can change and are reported.',
        interpretation='Mechanism controls on a fixed outcome cohort; no population rate or optimized production error claimed.'))


def smoke(root, destination):
    assert (destination/'jobs/build/exit_code.txt').read_text().strip() == '0'
    overlay = destination/'overlay'
    overlay.mkdir(exist_ok=False)
    links = {}
    for folder in ['lib', 'python']:
        (overlay/folder).mkdir()
        for source in [destination/'software', root/'observer_hitcontrol_v2', root/'cluster_measurement_v1/software']:
            for item in (source/folder).iterdir():
                target = overlay/folder/item.name
                assert not target.exists()
                target.symlink_to(item.resolve(), target_is_directory=item.is_dir())
                links[str(target.relative_to(overlay))] = str(item.resolve())
    dump(overlay/'software.json', dict(files=links, native_tracking=True))
    selected = []
    original_sources = sources(root)
    for dsid in ['100116', '100117']:
        for wrong in [True, False]:
            selected.append(next(r for r in original_sources if r['dsid'] == dsid
                and bool(r['replay_targets'][0]['wrong']) == wrong))
    jobs = [task(root, destination, row, mode) for row in selected for mode in
            ['native', 'centroid_fixed', 'ckf_quarter', 'ckf_four', 'ckf_cluster']]
    dump(destination/'smoke_tasks.json', jobs)
    dump(destination/'smoke_validation_tasks.json', [dict(kind='analysis', command=['python3',
        str(destination/'code/scripts/prepare_ckf_measurements.py'), '--root', str(root),
        '--revision', destination.name, '--phase', 'validate-smoke'], output=str(destination/'jobs/validate_smoke'))])


def validate_smoke(destination):
    import numpy as np
    from validate_charge_audit import compare_tables
    from cluster_coordinate_io import read_coordinates
    checks = []
    for job in json.loads((destination/'smoke_tasks.json').read_text()):
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0', str(folder)
        current = json.loads((folder/'source_task.json').read_text())
        mode = current['coordinate_mode']
        run = json.loads((folder/'run.json').read_text())
        assert run['isolate_smoothing'] and run['args']['cluster_coordinate_audit']
        with gzip.open(folder/'charge_audit.jsonl.gz', 'rt') as file:
            audit = [json.loads(line) for line in file]
        with gzip.open(Path(current['original_output'])/'charge_audit.jsonl.gz', 'rt') as file:
            original = [json.loads(line) for line in file]
        seeds = [{r['seed']: r for r in rows if r['stage'] == 'seed'} for rows in [original, audit]]
        assert seeds[0].keys() == seeds[1].keys()
        for index in seeds[0]:
            a, b = seeds[0][index], seeds[1][index]
            for key in ['bound', 'cov', 'surface']:
                assert a[key] == b[key], (str(folder), key)
            assert sorted(a['hits']) == sorted(b['hits'])
        record = dict(output=str(folder), mode=mode, full_seed_mean_covariance_surface_exact=True)
        if mode == 'native':
            fields = compare_tables(folder, Path(current['original_output']))
            assert len(fields) == 58
            record['corrected_reference_fields_exact'] = len(fields)
        coordinates = {r['cluster_id']: r for r in read_coordinates(folder/'cluster_coordinates.jsonl.gz')}
        n = 0
        for row in audit:
            if row['stage'] != 'selected_state_trace':
                continue
            for step in row['steps']:
                if not step['measurement_flag']:
                    continue
                hit = coordinates[step['hit']]
                factor = {'ckf_quarter': .25, 'ckf_four': 4.}.get(mode, 1.)
                expected = hit['used_covariance_mm2'][0] if mode == 'ckf_cluster' else .0231**2*factor
                np.testing.assert_array_equal(step['measurement_cov'], [[expected]])
                np.testing.assert_array_equal(step['measurement'], [hit['used_local_mm'][0]])
                n += 1
        record['actual_ckf_measurement_covariances_checked'] = n
        assert n > 0
        checks.append(record)
    dump(destination/'smoke_validation.json', dict(all_exact=True, jobs=len(checks), checks=checks))


def replay(root, destination):
    assert json.loads((destination/'smoke_validation.json').read_text())['all_exact']
    done = {r['output'] for r in json.loads((destination/'smoke_tasks.json').read_text())}
    jobs = [task(root, destination, original, mode) for original in sources(root) for mode in
            ['native', 'centroid_fixed', 'ckf_quarter', 'ckf_four', 'ckf_cluster']]
    assert len(jobs) == 330
    dump(destination/'all_tasks.json', jobs)
    remaining = [job for job in jobs if job['output'] not in done]
    assert len(remaining) == 310
    dump(destination/'replay_tasks.json', remaining)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--revision', default='ckf_measurement_v1')
    parser.add_argument('--phase', choices=['build', 'smoke', 'validate-smoke', 'replay'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    destination = root/args.revision
    if args.phase == 'build':
        prepare_build(root, destination, args.source.resolve())
    elif args.phase == 'smoke':
        smoke(root, destination)
    elif args.phase == 'validate-smoke':
        validate_smoke(destination)
    else:
        replay(root, destination)
    print(destination, args.phase)


if __name__ == '__main__':
    main()
