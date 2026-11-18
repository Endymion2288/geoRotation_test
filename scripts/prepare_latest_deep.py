#!/usr/bin/env python3
"""Prepare isolated build and smoke validation for the latest-Calypso deep momentum study."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    latest = repo.parent / 'geoRotation_calypso_latest'
    source = args.source.resolve()
    assert source.exists(), source
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)

    patch = repo / 'patches/0049_latest_seed_and_cov_audit.patch'
    diff = subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True)
    assert diff == patch.read_text(), 'Worktree diff does not match 0049_latest_seed_and_cov_audit.patch'
    subprocess.run(['git', '-C', str(source), 'apply', '--reverse', '--check', str(patch)], check=True)
    (output / 'source.patch').write_text(diff)
    shutil.copy2(repo / 'patches/LatestSeedAndStateAudit.h', output / 'LatestSeedAndStateAudit.h')

    compressed = latest / 'Tracking/Acts/FaserActsGeometry/data/material-maps-alma9-rotated-rowmajor.json.gz'
    material = output / 'material-latest-rowmajor.json'
    material.write_bytes(gzip.decompress(compressed.read_bytes()))

    reference = output / 'reference'
    old_ref = repo.parent / 'geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/reference/adapted'
    reference.mkdir()
    for name in ['ckf_residual_perLayer.C', 'CommonUtils.h', 'SetStyle.hpp']:
        shutil.copy2(old_ref / name, reference / name)

    code = output / 'code_setup'
    for relative in (
        'scripts/build_latest_momentum.py',
        'scripts/run_latest_deep.py',
        'scripts/reference_truth_observer.py',
        'configs/latest_deep_tracking.py',
        'analysis/latest_deep_reduce.py',
    ):
        target = code / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / relative, target)

    build_task = dict(
        kind='analysis',
        command=[
            'python3', str(code / 'scripts/build_latest_momentum.py'),
            '--source', str(source),
            '--latest', str(latest),
            '--output', str(output / 'software_deep_v1'),
            '--jobs', '4',
        ],
        output=str(output / 'jobs/build_deep_v1'),
    )
    (output / 'build_tasks.json').write_text(json.dumps([build_task], indent=2) + '\n')

    rdo0 = '/eos/experiment/faser/data0/sim/mc24/foresee/110051/rdo/master_dev/FaserMC-MC24_FS_Aee_010MeV_1Em4-110051-00000-00009-master_dev-RDO.root'
    for mode, audit_flag, retain in [('native_smoke', False, False), ('audit_smoke', True, True)]:
        src_task = dict(
            input=rdo0, file_index=0, events=20, skip=0,
            output=str(output / 'runs' / mode),
            material=str(material),
            audit=audit_flag, retain_audit=retain, options=[],
            macro=str(reference / 'ckf_residual_perLayer.C') if audit_flag else None,
        )
        src_file = output / f'{mode}_source.json'
        src_file.write_text(json.dumps(src_task, indent=2) + '\n')
        job = dict(
            kind='analysis',
            command=['python3', str(code / 'scripts/run_latest_deep.py'), '--task', str(src_file)],
            output=str(output / 'runs' / mode),
        )
        (output / f'{mode}_tasks.json').write_text(json.dumps([job], indent=2) + '\n')

    manifest = dict(
        latest=str(latest),
        source=str(source),
        commit=subprocess.check_output(['git', '-C', str(latest), 'rev-parse', 'HEAD'], text=True).strip(),
        material=str(material),
        material_source=str(compressed),
        material_bytes=material.stat().st_size,
        source_patch=str(patch),
        isolate_smoothing=True,
        sample='110051 all 20 files / 20000 entries',
    )
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(output)


if __name__ == '__main__':
    main()
