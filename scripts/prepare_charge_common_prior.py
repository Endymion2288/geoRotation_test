#!/usr/bin/env python3
"""Freeze native Calypso observation and complete-seed gates for Q5."""
import argparse
import gzip
import json
from pathlib import Path
import shutil

from prepare_cluster_measurements import dump


VARIANTS = ['beforeGeo', 'afterGeo']
SAMPLES = ['100044', '100116', '100117']


def seed_rows(folder):
    entry = -1
    rows = []
    with gzip.open(folder/'prior_audit.jsonl.gz', 'rt') as source:
        for line in source:
            row = json.loads(line)
            if row['stage'] == 'event':
                entry += 1
            if row['stage'] == 'seed':
                rows.append(dict(row, local_event_index=entry))
    assert entry == 19
    return rows


def smoke(root, output):
    assert (output/'jobs/build_before/exit_code.txt').read_text().strip() == '0'
    repo = Path(__file__).resolve().parents[1]
    code = output/'smoke_code_v1'
    for folder in ['scripts', 'configs']:
        (code/folder).mkdir(parents=True, exist_ok=False)
    for name in ['run_slope_task.py', 'validate_charge_audit.py', 'prepare_cluster_measurements.py',
                 'prepare_charge_common_prior.py', 'ntuple_tracking_check.py']:
        shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
    shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
    overlays = {}
    for variant in VARIANTS:
        overlay = output/('overlay_'+variant)
        inputs = ([output/'software_before', root/'observer'] if variant == 'beforeGeo'
                  else [root/'fixed_hit_process_v1/overlay'])
        links = {}
        for folder in ['lib', 'python']:
            (overlay/folder).mkdir(parents=True)
            for package in inputs:
                for item in (package/folder).iterdir():
                    target = overlay/folder/item.name
                    assert not target.exists()
                    target.symlink_to(item.resolve(), target_is_directory=item.is_dir())
                    links[str(target.relative_to(overlay))] = str(item.resolve())
        dump(overlay/'software.json', dict(files=links, variant=variant, actual_calypso=True))
        overlays[variant] = str(overlay/'lib')
    (output/'tasks').mkdir()
    for variant in VARIANTS:
        jobs = []
        for dsid in SAMPLES:
            baseline = root/'charge_smoke'/variant/dsid/'on'
            original = json.loads((baseline/'source_task.json').read_text())
            for mode in ['default', 'audit']:
                folder = output/'smoke'/variant/dsid/mode
                options = ['--no-seed-writer', '--actual-covariance', '--isolate-smoothing']
                if mode == 'audit':
                    options += ['--prior-component-audit']
                task = dict(original, runner=str(code/'configs/slope_tracking.py'), options=options,
                    output=str(folder), reference=str(baseline), component_mode=mode)
                path = output/'tasks'/(variant+'_'+dsid+'_'+mode+'.json')
                dump(path, task)
                jobs.append(dict(kind='analysis', command=['python3',
                    str(code/'scripts/run_slope_task.py'), '--task', str(path)], output=str(folder)))
        dump(output/(variant+'_smoke_tasks.json'), jobs)
    dump(output/'smoke_manifest.json', dict(libraries=overlays,
        gates='Exact reconstruction arrays versus recorded native baseline; audit versus default; explicit native-frame before observer',
        scope='Observation gates only. Original baseline material is retained for exact replay; the common physical-prior experiment uses its separately validated source/rotated map pair.'))
    dump(output/'smoke_validation_tasks.json', [dict(kind='analysis', command=['python3',
        str(code/'scripts/prepare_charge_common_prior.py'), '--root', str(root),
        '--phase', 'validate-smoke'], output=str(output/'jobs/smoke_validation'))])


def validate_smoke(root, output):
    from validate_charge_audit import compare_tables
    checks = []
    for variant in VARIANTS:
        for dsid in SAMPLES:
            folder = output/'smoke'/variant/dsid
            for mode in ['default', 'audit']:
                assert (folder/mode/'exit_code.txt').read_text().strip() == '0'
            fields = compare_tables(folder/'default', root/'charge_smoke'/variant/dsid/'on')
            audit_fields = compare_tables(folder/'audit', folder/'default')
            assert len(fields) == len(audit_fields) == 58
            seeds = seed_rows(folder/'audit')
            assert seeds and all(len(row['bound']) == len(row['cov']) == 6 and
                len(row['surface']) == 4 for row in seeds)
            if variant == 'beforeGeo':
                assert all(row['parameter_frame'] == 'FASER' and
                    'position_native_mm' in row and 'momentum_native_GeV' in row for row in seeds)
            checks.append(dict(variant=variant, dsid=dsid, default_fields_exact=len(fields),
                audit_fields_exact=len(audit_fields), complete_seed_rows=len(seeds)))
    dump(output/'smoke_validation.json', dict(all_exact=True, checks=checks))


def identity(root, output):
    assert json.loads((output/'smoke_validation.json').read_text())['all_exact']
    (output/'identity_inputs').mkdir()
    code = output/'smoke_code_v1'
    for variant in VARIANTS:
        jobs = []
        for dsid in SAMPLES:
            folder = output/'smoke'/variant/dsid/'audit'
            original = json.loads((folder/'source_task.json').read_text())
            lines = []
            for row in seed_rows(folder):
                values = [*row['bound'], *[v for r in row['cov'] for v in r],
                          *[v for r in row['surface'] for v in r]]
                assert len(values) == 58
                lines.append(' '.join([str(row['local_event_index']), str(row['seed']),
                                      *[format(v, '.17g') for v in values]]))
            inputs = output/'identity_inputs'/(variant+'_'+dsid+'.txt')
            inputs.write_text('\n'.join(lines)+'\n')
            destination = folder.parent/'full_identity'
            task = dict(original, output=str(destination),
                options=[*original['options'], '--full-seed-input-file', str(inputs)])
            path = output/'tasks'/(variant+'_'+dsid+'_full_identity.json')
            dump(path, task)
            jobs.append(dict(kind='analysis', command=['python3',
                str(code/'scripts/run_slope_task.py'), '--task', str(path)], output=str(destination)))
        dump(output/(variant+'_identity_tasks.json'), jobs)
    dump(output/'identity_validation_tasks.json', [dict(kind='analysis', command=['python3',
        str(code/'scripts/prepare_charge_common_prior.py'), '--root', str(root),
        '--phase', 'validate-identity'], output=str(output/'jobs/identity_validation'))])


def validate_identity(root, output):
    from validate_charge_audit import compare_tables
    checks = []
    for variant in VARIANTS:
        for dsid in SAMPLES:
            folder = output/'smoke'/variant/dsid
            assert (folder/'full_identity/exit_code.txt').read_text().strip() == '0'
            assert len(compare_tables(folder/'full_identity', folder/'audit')) == 58
            a, b = seed_rows(folder/'audit'), seed_rows(folder/'full_identity')
            assert len(a) == len(b)
            for old, new in zip(a, b):
                for key in ['local_event_index', 'seed', 'bound', 'cov', 'surface']:
                    assert old[key] == new[key], (variant, dsid, key)
            checks.append(dict(variant=variant, dsid=dsid, fields_exact=58,
                complete_mean_covariance_surface_exact=True, seed_rows=len(a)))
    dump(output/'identity_validation.json', dict(all_exact=True, checks=checks))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['smoke', 'validate-smoke', 'identity', 'validate-identity'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output = root/'common_prior_ckf_v1'
    {'smoke':smoke, 'validate-smoke':validate_smoke, 'identity':identity,
     'validate-identity':validate_identity}[args.phase](root, output)
    print(output, args.phase)


if __name__ == '__main__':
    main()
