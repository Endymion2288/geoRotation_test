#!/usr/bin/env python3
"""Freeze native Calypso strip-centre controls that retain strip quantization."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_cluster_measurements import dump, sources, override


def build(root, destination, source):
    repo = Path(__file__).resolve().parents[1]
    destination.mkdir(exist_ok=False)
    code = destination/'build_code'
    code.mkdir()
    for name in ['build_cluster_measurements.py', 'build_slope_diagnostics.py']:
        shutil.copy2(repo/'scripts'/name, code/name)
    diff = subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True)
    assert 'UseSDOSupportedStripPosition' in diff
    (destination/'source.patch').write_text(diff)
    dump(destination/'build_tasks.json', [dict(kind='analysis', command=['python3',
        str(code/'build_cluster_measurements.py'), '--source', str(source),
        '--source-patch', str(destination/'source.patch'), '--output', str(destination/'software'),
        '--jobs', '4'], output=str(destination/'jobs/build'))])
    dump(destination/'manifest.json', dict(source=str(source),
        source_commit='88c457c1e4bf4154e67b96114ee8e7020d9b1a97',
        material=str(root/'material_encoding_control_v1/corrected_afs_map.json'),
        coordinate='Unweighted mean of actual native strip centres with a positive deposit from the cluster majority barcode; original Lorentz shift added',
        invariants='Original RDO identifiers, cluster width/covariance, local Y, all six seed means/full covariance/surface; isolated smoothing ON; native fixed CKF variance',
        fallback='Keep the original cluster coordinate when no supported strip exists',
        scope='Frozen six pure wrong events and 60 original normal controls; truth-assisted mechanism diagnostic, not a production selection or population rate'))


def task(root, destination, original, mode):
    options = ['--no-seed-writer', '--actual-covariance', '--isolate-smoothing',
               '--stage-audit', '--charge-audit', '--cluster-coordinate-audit']
    if mode == 'sdo_strips_fixed':
        seed = destination/'seeds'/(original['dsid']+f"_{original['sample_entry']:08d}.txt")
        if not seed.exists():
            override(Path(original['output']), seed)
        options += ['--seed-override-file', str(seed), '--seed-override-mode', 'mean',
                    '--sdo-supported-strip-position']
    output = destination/'reconstruction'/original['dsid']/mode/f"e{original['sample_entry']:08d}"
    current = dict(original, original_output=original['output'], output=str(output),
        runner=str(destination/'reconstruction_code/configs/slope_tracking.py'),
        coordinate_mode=mode, options=options)
    config = destination/'tasks'/(original['dsid']+'_'+mode+f"_{original['sample_entry']:08d}.json")
    dump(config, current)
    worker = 'run_charge_replay.py' if mode == 'native' else 'run_slope_task.py'
    return dict(kind='analysis', command=['python3', str(destination/'reconstruction_code/scripts'/worker),
                '--task', str(config)], output=str(output))


def smoke(root, destination):
    assert (destination/'jobs/build/exit_code.txt').read_text().strip() == '0'
    repo = Path(__file__).resolve().parents[1]
    code = destination/'reconstruction_code'
    for folder in ['scripts', 'configs']:
        (code/folder).mkdir(parents=True, exist_ok=False)
    for name in ['run_slope_task.py', 'run_charge_replay.py', 'validate_charge_audit.py',
                 'audit_reduce.py', 'charge_audit.py', 'slope_stage_extract.py',
                 'ntuple_tracking_check.py', 'cluster_coordinate_io.py',
                 'prepare_cluster_measurements.py', 'prepare_sdo_strip_coordinates.py']:
        shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
    shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
    assert '--sdo-supported-strip-position' in (code/'configs/slope_tracking.py').read_text()
    for folder in ['tasks', 'seeds']:
        (destination/folder).mkdir()
    overlay = destination/'overlay'
    links = {}
    for folder in ['lib', 'python']:
        (overlay/folder).mkdir(parents=True)
        for source in [root/'ckf_measurement_v1/software', root/'observer_hitcontrol_v2', destination/'software']:
            for item in (source/folder).iterdir():
                target = overlay/folder/item.name
                assert not target.exists()
                target.symlink_to(item.resolve(), target_is_directory=item.is_dir())
                links[str(target.relative_to(overlay))] = str(item.resolve())
    dump(overlay/'software.json', dict(files=links, native_tracking=True))
    selected = []
    original = sources(root)
    for dsid in ['100116', '100117']:
        for wrong in [True, False]:
            selected.append(next(r for r in original if r['dsid'] == dsid and
                            bool(r['replay_targets'][0]['wrong']) == wrong))
    jobs = [task(root, destination, row, mode) for row in selected for mode in ['native', 'sdo_strips_fixed']]
    dump(destination/'smoke_tasks.json', jobs)
    dump(destination/'smoke_validation_tasks.json', [dict(kind='analysis', command=['python3',
        str(code/'scripts/prepare_sdo_strip_coordinates.py'), '--root', str(root),
        '--phase', 'validate-smoke'], output=str(destination/'jobs/validate_smoke'))])


def coordinate_check(native, current, enabled):
    import numpy as np
    assert native.keys() == current.keys()
    changed, supported = 0, 0
    for key, old in native.items():
        new = current[key]
        for field, value in old.items():
            if field not in ['used_local_mm', 'sdo_supported_position_enabled']:
                assert new[field] == value, (key, field)
        assert new['sdo_supported_position_enabled'] == enabled
        channels = [r for r in old['channels'] if old['majority_barcode'] and
                    any(b == old['majority_barcode'] and q > 0 for b, q in r['deposits'])]
        assert len(channels) == new['sdo_supported_channels']
        expected = sum(r['centre_local_mm'][0] for r in channels)/len(channels)+old['lorentz_shift_mm'] if channels else None
        assert new['sdo_supported_local_x_mm'] == expected
        used = expected if enabled and expected is not None else old['native_local_mm'][0]+old['lorentz_shift_mm']
        np.testing.assert_array_equal(new['used_local_mm'], [used, old['native_local_mm'][1]])
        changed += new['used_local_mm'] != old['used_local_mm']
        supported += bool(channels)
    return dict(clusters=len(native), supported_clusters=supported, coordinates_changed=changed,
                native_width_covariance_channels_sihits_exact=True)


def validate_smoke(root, destination):
    import gzip
    import numpy as np
    from cluster_coordinate_io import read_coordinates
    from validate_charge_audit import compare_tables
    checks = []
    expected = {str((destination/'software/lib'/name).resolve()) for name in
                ['libFaserSiClusterizationTool.so', 'libFaserSiClusterizationToolLib.so']}
    for job in json.loads((destination/'smoke_tasks.json').read_text()):
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0', str(folder)
        current = json.loads((folder/'source_task.json').read_text())
        mode = current['coordinate_mode']
        run = json.loads((folder/'run.json').read_text())
        assert run['success'] and run['isolate_smoothing']
        assert set(run['loaded_cluster_libraries']) == expected
        with gzip.open(folder/'charge_audit.jsonl.gz', 'rt') as file:
            audit = [json.loads(line) for line in file]
        with gzip.open(Path(current['original_output'])/'charge_audit.jsonl.gz', 'rt') as file:
            original = [json.loads(line) for line in file]
        seeds = [{r['seed']: r for r in rows if r['stage'] == 'seed'} for rows in [original, audit]]
        assert seeds[0].keys() == seeds[1].keys()
        for index in seeds[0]:
            for key in ['bound', 'cov', 'surface']:
                assert seeds[0][index][key] == seeds[1][index][key]
            assert sorted(seeds[0][index]['hits']) == sorted(seeds[1][index]['hits'])
        rows = {r['cluster_id']: r for r in read_coordinates(folder/'cluster_coordinates.jsonl.gz')}
        native_folder = folder.parent.parent/'native'/folder.name
        native = {r['cluster_id']: r for r in read_coordinates(native_folder/'cluster_coordinates.jsonl.gz')}
        record = dict(output=str(folder), full_seed_exact=True,
            coordinates=coordinate_check(native, rows, mode == 'sdo_strips_fixed'))
        if mode == 'native':
            assert len(compare_tables(folder, Path(current['original_output']))) == 58
            record['corrected_reference_fields_exact'] = 58
        n = 0
        for row in audit:
            if row['stage'] != 'selected_state_trace':
                continue
            for step in row['steps']:
                if step['measurement_flag']:
                    np.testing.assert_array_equal(step['measurement_cov'], [[.0231**2]])
                    np.testing.assert_array_equal(step['measurement'], [rows[step['hit']]['used_local_mm'][0]])
                    n += 1
        assert n > 0
        record['actual_ckf_coordinates_variances_checked'] = n
        checks.append(record)
    dump(destination/'smoke_validation.json', dict(all_exact=True, jobs=len(checks), checks=checks))


def replay(root, destination):
    assert json.loads((destination/'smoke_validation.json').read_text())['all_exact']
    done = {r['output'] for r in json.loads((destination/'smoke_tasks.json').read_text())}
    jobs = [task(root, destination, original, 'sdo_strips_fixed') for original in sources(root)]
    assert len(jobs) == 66
    dump(destination/'all_tasks.json', jobs)
    remaining = [r for r in jobs if r['output'] not in done]
    assert len(remaining) == 62
    dump(destination/'replay_tasks.json', remaining)


def analysis(root, destination):
    assert json.loads((destination/'smoke_validation.json').read_text())['all_exact']
    repo = Path(__file__).resolve().parents[1]
    code = destination/'analysis_code_v1'
    for folder in ['analysis', 'scripts']:
        shutil.copytree(repo/folder, code/folder, ignore=shutil.ignore_patterns('__pycache__'))
    for label in ['smoke', 'full']:
        command = ['python3', str(code/'analysis/sdo_strip_coordinates.py'), '--root', str(root),
                   '--output', str(destination/'analysis_v1'/label)]
        if label == 'smoke':
            command += ['--smoke']
        dump(destination/(label+'_analysis_tasks.json'), [dict(kind='analysis', command=command,
             output=str(destination/'jobs'/('analysis_v1_'+label)))])
    figures = destination/'figures_v1'
    figures.mkdir()
    tools = '/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/software_tools'
    for phase in ['plot', 'validate']:
        dump(figures/(phase+'_tasks.json'), [dict(kind='analysis', command=['python3',
            str(code/'analysis/sdo_strip_figures.py'), '--phase', phase, '--input',
            str(destination/'analysis_v1/full'), '--output', str(figures/'results'),
            '--pdf-tools', tools], output=str(figures/'jobs'/phase))])
    dump(figures/'report_tasks.json', [dict(kind='analysis', command=['python3',
        str(code/'analysis/charge_report.py'), '--figures', str(root/'report_v8'),
        str(figures/'results'), '--output', str(root/'report_v9'), '--pdf-tools', tools],
        output=str(figures/'jobs/report'))])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--phase', choices=['build', 'smoke', 'validate-smoke', 'replay', 'analysis'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    destination = root/'sdo_strip_coordinate_v1'
    if args.phase == 'build':
        build(root, destination, args.source.resolve())
    elif args.phase == 'smoke':
        smoke(root, destination)
    elif args.phase == 'validate-smoke':
        validate_smoke(root, destination)
    elif args.phase == 'analysis':
        analysis(root, destination)
    else:
        replay(root, destination)
    print(destination, args.phase)


if __name__ == '__main__':
    main()
