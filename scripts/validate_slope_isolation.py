#!/usr/bin/env python3
"""Compare actual Calypso observer outputs for default, explicit ON and OFF."""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    result = {}
    for dsid in ['100044', '100116', '100117']:
        base = args.root/dsid
        runs = {mode: json.loads((base/mode/'run.json').read_text())
                for mode in ['default', 'isolated', 'shared']}
        for mode, run in runs.items():
            assert (base/mode/'exit_code.txt').read_text().strip() == '0'
            assert run['success'] and run['isolate_smoothing'] == (mode != 'shared')
            assert len(run['loaded_acts_libraries']) == 4
        equal, changes, counts = [], {}, {}
        for name in ['events', 'truth', 'tracks', 'segments']:
            with np.load(base/'default'/(name+'.npz')) as a, \
                 np.load(base/'isolated'/(name+'.npz')) as b, \
                 np.load(base/'shared'/(name+'.npz')) as c:
                assert a.files == b.files == c.files
                counts[name] = len(a['source_entry'])
                for key in a.files:
                    np.testing.assert_array_equal(a[key], b[key], err_msg=dsid+'/'+name+'/'+key)
                    equal.append(name+'/'+key)
                    same = (a[key] == c[key]) | (np.isnan(a[key]) & np.isnan(c[key]))
                    if not same.all():
                        changes[name+'/'+key] = int((~same).sum())
        result[dsid] = dict(counts=counts, default_explicit_on_equal_arrays=equal,
                            isolated_vs_shared_changed_elements=changes,
                            actual_loaded_libraries=runs['default']['loaded_acts_libraries'])
    (args.root/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:dict(counts=v['counts'], equal_arrays=len(v['default_explicit_on_equal_arrays']),
          changed=v['isolated_vs_shared_changed_elements']) for k,v in result.items()},indent=2))


if __name__ == '__main__':
    main()
