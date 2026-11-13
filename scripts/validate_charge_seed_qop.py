#!/usr/bin/env python3
"""Validate exact default closure and q/p-only changes in actual Calypso seeds."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from validate_charge_audit import compare_tables


def audit(folder):
    with gzip.open(folder/'charge_audit.jsonl.gz', 'rt') as source:
        return [json.loads(line) for line in source]


def seeds(folder):
    return {(row['source_entry'], row['seed']): row for row in audit(folder) if row['stage'] == 'seed'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--mean-smoke', action='store_true', help='Check identical complete-mean overrides')
    args = parser.parse_args()
    root = args.root/('seed_mean_smoke_v1' if args.mean_smoke else 'seed_qop_smoke_v1')
    records = []
    for dsid in ['100044', '100116', '100117']:
        original = args.root/'hitcontrol_smoke_v2'/dsid/'default/f00_e00000000'
        old = seeds(original)
        for mode in (['identical_mean'] if args.mean_smoke else ['identical', 'production_qop']):
            folder = root/dsid/mode/'f00_e00000000'
            assert (folder/'exit_code.txt').read_text().strip() == '0'
            run = json.loads((folder/'run.json').read_text())
            assert run['success'] and run['isolate_smoothing']
            assert run['seed_override_mode'] == ('mean' if args.mean_smoke else 'qop')
            new = seeds(folder)
            assert old.keys() == new.keys()
            override = np.atleast_2d(np.loadtxt(run['seed_override_file']))
            expected = {(int(row[0]), int(row[1])): row[6] for row in override}
            changed = 0
            for key in old:
                a, b = old[key], new[key]
                for field in ['hits', 'surface', 'cov']:
                    assert a[field] == b[field]
                np.testing.assert_array_equal(np.delete(a['bound'], 4), np.delete(b['bound'], 4))
                assert b['bound'][4] == expected[key]
                changed += a['bound'][4] != b['bound'][4]
            if mode in ['identical', 'identical_mean']:
                fields = compare_tables(folder, original)
                assert len(fields) == 58 and changed == 0
            else:
                assert changed > 0
            records.append(dict(dsid=dsid, mode=mode, seed_rows=len(old), qop_changed=changed,
                fields_exact=58 if mode in ['identical', 'identical_mean'] else None))
    result = dict(samples=records, all_identical_58_fields_exact=True,
        all_changes_qop_only_with_identical_covariance=True,
        scope='Three DSIDs, 20 events each; actual Calypso existing seed override, no tracking replacement')
    if args.mean_smoke:
        result.pop('all_changes_qop_only_with_identical_covariance')
        result['all_full_mean_overrides_exact'] = True
    (root/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
