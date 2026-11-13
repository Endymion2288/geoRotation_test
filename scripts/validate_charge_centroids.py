#!/usr/bin/env python3
"""Require default closure and actual centroid-library use before expansion."""
import argparse
import gzip
import json
from pathlib import Path
import re

from validate_charge_audit import compare_tables


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root/'centroid_smoke_v1'
    samples = {}
    expected = {str((args.root.parent.parent/'stage05_slopes/software_cluster_control_v2/lib'/name).resolve())
                for name in ['libFaserSiClusterizationTool.so', 'libFaserSiClusterizationToolLib.so']}
    for dsid in ['100044', '100116', '100117']:
        outputs = {mode: root/dsid/mode/'f00_e00000000' for mode in ['default', 'centroid']}
        fields = compare_tables(outputs['default'], args.root/'charge_smoke/afterGeo'/dsid/'off')
        assert len(fields) == 58
        sample = {}
        for mode, folder in outputs.items():
            assert (folder/'exit_code.txt').read_text().strip() == '0'
            run = json.loads((folder/'run.json').read_text())
            assert run['success'] and run['isolate_smoothing']
            assert set(run['loaded_cluster_libraries']) == expected
            assert run['args']['truth_cluster_position'] == (mode == 'centroid')
            log = gzip.open(folder/'tracking.log.gz', 'rt').read()
            values = re.findall(r'Truth cluster position: enabled=.*?(True|False)\s+replaced=(\d+) fallback=(\d+)', log)
            assert len(values) == 1
            enabled, replaced, fallback = values[0]
            assert (enabled == 'True') == (mode == 'centroid')
            replaced, fallback = int(replaced), int(fallback)
            assert (replaced > 0) if mode == 'centroid' else (replaced == 0)
            sample[mode] = dict(replaced=replaced, fallback=fallback, actual_libraries=run['loaded_cluster_libraries'])
        samples[dsid] = dict(default_fields_exact=len(fields), modes=sample)
    result = dict(samples=samples, all_default_58_fields_exact=True, all_actual_cluster_libraries_match=True,
        scope='Actual Calypso default vs same-wafer majority-SiHit centroid; 20 source events for each of three DSIDs')
    (root/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
