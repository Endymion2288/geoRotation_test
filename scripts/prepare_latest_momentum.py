#!/usr/bin/env python3
"""Freeze latest-Calypso instrumentation for the station-two momentum investigation."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess


def replace(text, before, after):
    assert text.count(before) == 1, before
    return text.replace(before, after)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True, help='Fresh isolated worktree path')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    latest = repo.parent/'geoRotation_calypso_latest'
    source = args.source.resolve()
    if not source.exists():
        subprocess.run(['git','-C',str(latest),'worktree','add','--detach','--no-checkout',str(source),'HEAD'],check=True)
        subprocess.run(['git','-C',str(source),'sparse-checkout','set','Tracking/Acts/FaserActsKalmanFilter'],check=True)
        subprocess.run(['git','-C',str(source),'checkout','HEAD'],check=True)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    package = source/'Tracking/Acts/FaserActsKalmanFilter/src'
    # The patch and its added header are the complete versioned source change.
    patch = repo/'patches/0048_latest_momentum_audit.patch'
    subprocess.run(['git','-C',str(source),'apply','--check',str(patch)],check=True)
    subprocess.run(['git','-C',str(source),'apply',str(patch)],check=True)
    shutil.copy2(repo/'patches/LatestMomentumAudit.h',package/'LatestMomentumAudit.h')
    # Explicit material bytes: use exactly the map packaged by latest, not either historical JSON.
    compressed = latest/'Tracking/Acts/FaserActsGeometry/data/material-maps-alma9-rotated-rowmajor.json.gz'
    material = output/'material-latest-rowmajor.json'
    material.write_bytes(gzip.decompress(compressed.read_bytes()))
    diff = subprocess.check_output(['git','-C',str(source),'diff','HEAD'],text=True)
    patch = repo/'patches/0048_latest_momentum_audit.patch'
    patch.write_text(diff)
    (output/'source.patch').write_text(diff)
    shutil.copy2(package/'LatestMomentumAudit.h',output/'LatestMomentumAudit.h')
    subprocess.run(['git','-C',str(source),'apply','--reverse','--check',str(patch)],check=True)
    code = output/'code'
    for relative in ('scripts/build_latest_momentum.py','configs/reference_tracking.py','scripts/reference_truth_observer.py'):
        target=code/relative
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(repo/relative,target)
    task=dict(kind='analysis',command=['python3',str(code/'scripts/build_latest_momentum.py'),
        '--source',str(source),'--latest',str(latest),'--output',str(output/'software'),'--jobs','4'],
        output=str(output/'jobs/build'))
    (output/'build_tasks.json').write_text(json.dumps([task],indent=2)+'\n')
    (output/'manifest.json').write_text(json.dumps(dict(latest=str(latest),source=str(source),
        commit=subprocess.check_output(['git','-C',str(latest),'rev-parse','HEAD'],text=True).strip(),
        material=str(material), material_source=str(compressed),material_bytes=material.stat().st_size,
        source_patch=str(patch),isolate_smoothing=True, sample='110051 all 20 files / 20000 entries'),indent=2)+'\n')
    print(output)


if __name__ == '__main__':
    main()
