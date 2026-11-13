#!/usr/bin/env python3
"""Require unchanged native reconstruction and identical CKF before refit."""
import argparse
import json
from pathlib import Path

import numpy as np

from validate_charge_audit import compare_tables


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    output = args.root/'refit_control_smoke_v1'
    samples = {}
    for dsid in ['100116', '100117']:
        native = output/dsid/'native/shards/beforeGeo/s0000'
        original = args.root/'charge_smoke/beforeGeo'/dsid/'on'
        fields = compare_tables(native, original)
        checks = {}
        with np.load(native/'charge_lineage.npz') as reference:
            for mode in ['native', 'smoother', 'updater', 'both']:
                folder = output/dsid/mode/'shards/beforeGeo/s0000'
                assert (folder/'exit_code.txt').read_text().strip() == '0'
                with np.load(folder/'charge_lineage.npz') as current:
                    assert current.files == reference.files
                    a, b = current['stage'] != 7, reference['stage'] != 7
                    for field in current.files:
                        np.testing.assert_array_equal(current[field][a], reference[field][b],
                                                      err_msg=dsid+'/'+mode+'/'+field)
                    checks[mode] = dict(CKF_rows_exact=int(a.sum()), fields=len(current.files))
        samples[dsid] = dict(native_fields_exact=len(fields), modes=checks)
    result = dict(samples=samples, native_58_fields_exact=all(
        r['native_fields_exact'] == 58 for r in samples.values()), all_CKF_rows_exact=True,
        scope='20 actual RDO events per charge; only native refit chart switches differ')
    (output/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
