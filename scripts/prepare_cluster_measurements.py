#!/usr/bin/env python3
"""Freeze actual Calypso coordinate/error controls for the existing charge cohort."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def sources(root):
    retained = {(r['dsid'], r['sample_entry']) for r in json.loads(
        (root/'fixed_seed_analysis_v3/events.json').read_text())}
    assert len(retained) == 66
    rows = []
    for job in json.loads((root/'material_encoding_control_v1/replay_v1/tasks.json').read_text()):
        task = json.loads((Path(job['output'])/'source_task.json').read_text())
        if (task['dsid'], task['sample_entry']) in retained:
            rows.append(task)
    assert len(rows) == 66
    return rows


def override(folder, output):
    with gzip.open(folder/'charge_audit.jsonl.gz', 'rt') as file:
        seeds = [r for line in file if (r := json.loads(line))['stage'] == 'seed']
    task = json.loads((folder/'source_task.json').read_text())
    output.write_text('\n'.join(' '.join([str(r['source_entry']-task['skip']), str(r['seed']),
        *[format(x, '.17g') for x in r['bound']]]) for r in seeds)+'\n')
    return len(seeds)


def task(root, destination, original, mode):
    code = destination/'code'
    options = ['--no-seed-writer', '--actual-covariance', '--isolate-smoothing',
               '--stage-audit', '--charge-audit']
    if mode != 'audit_off':
        options += ['--cluster-coordinate-audit']
    if mode == 'centroid_fixed' or mode.startswith('variance_'):
        override_file = destination/'seeds'/(original['dsid']+f"_{original['sample_entry']:08d}.txt")
        if not override_file.exists():
            override(Path(original['output']), override_file)
        options += ['--seed-override-file', str(override_file), '--seed-override-mode', 'mean']
    if mode == 'centroid_fixed':
        options += ['--truth-cluster-position']
    if mode.startswith('variance_'):
        options += ['--measurement-variance-scale', '.25' if mode == 'variance_quarter' else '4.']
    output = destination/'reconstruction'/original['dsid']/mode/f"e{original['sample_entry']:08d}"
    current = dict(original, runner=str(code/'configs/slope_tracking.py'),
        original_output=original['output'], output=str(output), options=options, coordinate_mode=mode)
    source = destination/'tasks'/(original['dsid']+'_'+mode+f"_{original['sample_entry']:08d}.json")
    dump(source, current)
    runner = 'run_charge_replay.py' if mode in ['native', 'audit_off'] else 'run_slope_task.py'
    return dict(kind='analysis', command=['python3', str(code/'scripts'/runner), '--task', str(source)],
                output=str(output))


def prepare_build(root, destination, source):
    repo = Path(__file__).resolve().parents[1]
    destination.mkdir(exist_ok=False)
    for folder in ['tasks', 'seeds', 'code/scripts', 'code/configs']:
        (destination/folder).mkdir(parents=True)
    for name in ['build_cluster_measurements.py', 'build_slope_diagnostics.py', 'run_slope_task.py',
                 'run_charge_replay.py', 'validate_charge_audit.py', 'audit_reduce.py',
                 'charge_audit.py', 'slope_stage_extract.py', 'ntuple_tracking_check.py',
                 'prepare_cluster_measurements.py']:
        shutil.copy2(repo/'scripts'/name, destination/'code/scripts'/name)
    shutil.copy2(repo/'configs/slope_tracking.py', destination/'code/configs/slope_tracking.py')
    patch = subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True)
    assert patch and 'CoordinateAuditFile' in patch and 'MeasurementVarianceScale' in patch
    (destination/'source.patch').write_text(patch)
    jobs = [dict(kind='analysis', command=['python3',
        str(destination/'code/scripts/build_cluster_measurements.py'), '--source', str(source),
        '--source-patch', str(destination/'source.patch'), '--output', str(destination/'software'),
        '--jobs', '4'], output=str(destination/'jobs/build'))]
    dump(destination/'build_tasks.json', jobs)
    dump(destination/'manifest.json', dict(source=str(source), material=str(root/'material_encoding_control_v1/corrected_afs_map.json'),
        scope='Frozen six original purity-one wrong events and 60 normal controls; isolated smoothing ON',
        invariants='Actual Calypso cluster, segment/seed, CKF and refit. Audit defaults off; variance factor defaults one.',
        quantities='Native strip geometry, actual cluster channel/SDO lists and same-wafer SiHit endpoints; no invented sensor momentum.',
        gates='Require exact 58-field default/audit and corrected-reference closure before expansion; check complete seed and selected-hit identities separately.'))


def prepare_smoke(root, destination):
    assert (destination/'jobs/build/exit_code.txt').read_text().strip() == '0'
    overlay = destination/'overlay'
    overlay.mkdir(exist_ok=False)
    links = {}
    for folder in ['lib', 'python']:
        (overlay/folder).mkdir()
        for source in [root/'overlay_hitcontrol_v2', destination/'software']:
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
    jobs = [task(root, destination, original, mode) for original in selected for mode in ['audit_off', 'native']]
    dump(destination/'smoke_tasks.json', jobs)
    dump(destination/'smoke_validation_tasks.json', [dict(kind='analysis', command=['python3',
        str(destination/'code/scripts/prepare_cluster_measurements.py'), '--root', str(root),
        '--revision', destination.name, '--phase', 'validate-smoke'], output=str(destination/'jobs/validate_smoke'))])


def validate_smoke(destination):
    from cluster_coordinate_io import read_coordinates
    from validate_charge_audit import compare_tables
    checks = []
    jobs = json.loads((destination/'smoke_tasks.json').read_text())
    expected = {str((destination/'software/lib'/name).resolve()) for name in
                ['libFaserSiClusterizationTool.so', 'libFaserSiClusterizationToolLib.so']}
    for job in jobs:
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        task = json.loads((folder/'source_task.json').read_text())
        run = json.loads((folder/'run.json').read_text())
        assert run['isolate_smoothing'] and set(run['loaded_cluster_libraries']) == expected
        replay = json.loads((folder/'replay_validation.json').read_text())
        assert replay['all_reconstruction_arrays_identical'] and replay['fields_exact'] == 58
        record = dict(output=str(folder), corrected_reference_fields_exact=58)
        if task['coordinate_mode'] == 'native':
            fields = compare_tables(folder, folder.parent.parent/'audit_off'/folder.name)
            assert len(fields) == 58
            rows = list(read_coordinates(folder/'cluster_coordinates.jsonl.gz'))
            assert rows and all(row['variance_scale'] == 1 for row in rows)
            assert all(row['native_local_mm'] == row['used_local_mm'] for row in rows)
            record.update(audit_on_off_fields_exact=58, clusters=len(rows))
        checks.append(record)
    dump(destination/'smoke_validation.json', dict(all_exact=True, checks=checks, jobs=len(jobs)))


def prepare_replay(root, destination):
    assert json.loads((destination/'smoke_validation.json').read_text())['all_exact']
    done = {r['output'] for r in json.loads((destination/'smoke_tasks.json').read_text())}
    jobs = []
    for original in sources(root):
        modes = ['native', 'centroid_fixed']
        if original['replay_targets'][0]['wrong']:
            assert original['replay_targets'][0]['purity'] == 1
            modes += ['variance_quarter', 'variance_four']
        for mode in modes:
            job = task(root, destination, original, mode)
            if job['output'] not in done:
                jobs.append(job)
    assert len(jobs) == 140
    dump(destination/'replay_tasks.json', jobs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--revision', default='cluster_measurement_v1')
    parser.add_argument('--phase', choices=['build', 'smoke', 'validate-smoke', 'replay'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    destination = root/args.revision
    if args.phase == 'build':
        prepare_build(root, destination, args.source.resolve())
    elif args.phase == 'smoke':
        prepare_smoke(root, destination)
    elif args.phase == 'validate-smoke':
        validate_smoke(destination)
    else:
        prepare_replay(root, destination)
    print(destination, args.phase)


if __name__ == '__main__':
    main()
