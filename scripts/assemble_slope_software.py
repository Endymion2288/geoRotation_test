#!/usr/bin/env python3
"""Compose recorded native libraries and the separate S05 observer by symlink."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--variant', choices=['beforeGeo','afterGeo'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    short = 'before' if args.variant == 'beforeGeo' else 'after'
    base = root.parent/'stage05_geometry_comparison'/('software_'+short)
    sources = [base, root/'software_common']
    if short == 'before':
        sources.append(root/'software_native_ntuple_before')
    for source in sources:
        assert (source/'software.json').is_file(), source
    destination = root/('software_'+short)
    provenance = {}
    for folder in ['lib','python']:
        (destination/folder).mkdir(parents=True, exist_ok=True)
        for source in sources:
            for item in sorted((source/folder).iterdir()):
                target = destination/folder/item.name
                if target.is_symlink():
                    target.unlink()
                assert not target.exists(), target
                target.symlink_to(item.resolve(), target_is_directory=item.is_dir())
                provenance[str(target.relative_to(destination))] = str(item.resolve())
    (destination/'software.json').write_text(json.dumps(dict(variant=args.variant,
        native_sources=[str(s) for s in sources], files=provenance,
        tracking_implementation='Actual native Calypso CKF2 and KalmanFitterTool',
        observer_changes_reconstruction=False), indent=2)+'\n')
    print(destination/'lib')


if __name__ == '__main__':
    main()
