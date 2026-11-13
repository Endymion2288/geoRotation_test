#!/usr/bin/env python3
"""Verify detailed audit output against the frozen actual Calypso baseline."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np


def compare_tables(current, original, entry=None):
    checked = []
    for name in ['events', 'truth', 'tracks', 'segments']:
        with np.load(current/(name+'.npz')) as a, np.load(original/(name+'.npz')) as b:
            assert a.files == b.files
            mask = np.ones(len(b['source_entry']), dtype=bool) if entry is None else b['source_entry'] == entry
            if name == 'tracks':
                rows = np.flatnonzero(mask)
                offsets = b['hit_offsets']
                hits = np.concatenate([b['hit_identifiers'][offsets[i]:offsets[i+1]] for i in rows]) if len(rows) else np.array([], dtype=b['hit_identifiers'].dtype)
                expected_offsets = np.r_[0, np.cumsum(np.diff(offsets)[mask])]
            for field in b.files:
                expected = (hits if field == 'hit_identifiers' else expected_offsets) if name == 'tracks' and field in ['hit_identifiers', 'hit_offsets'] else b[field][mask]
                np.testing.assert_array_equal(a[field], expected, err_msg=str(current)+'/'+name+'/'+field)
                checked.append(name+'/'+field)
    return checked


def inspect_audit(folder):
    update_errors, smoothing_errors = [], []
    traces = 0
    with gzip.open(folder/'charge_audit.jsonl.gz', 'rt') as source:
        for line in source:
            row = json.loads(line)
            if row['stage'] not in ['selected_state_trace', 'refit_state_trace']:
                continue
            traces += 1
            for state in row['steps']:
                if 'update_identity_error' in state:
                    update_errors.append(abs(state['update_identity_error']))
                if 'smoothing_identity_error' in state:
                    smoothing_errors.append(abs(state['smoothing_identity_error']))
    return dict(traces=traces, updates=len(update_errors), smoothing_recursions=len(smoothing_errors),
        max_update_identity_error_GeVinv=max(update_errors, default=0.),
        max_periodic_smoother_identity_error_GeVinv=max(smoothing_errors, default=0.))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    result = {}
    for dsid in ['100116', '100117']:
        current = args.root/'full_audit_smoke_v1'/dsid
        assert (current/'exit_code.txt').read_text().strip() == '0'
        fields = compare_tables(current, args.root/'charge_smoke/afterGeo'/dsid/'off')
        run = json.loads((current/'run.json').read_text())
        assert run['audit'] and run['isolate_smoothing'] and run['charge_matrix_audit']
        result[dsid] = dict(fields_exact=len(fields), all_reconstruction_arrays_identical=True,
                            **inspect_audit(current))
    (args.root/'full_audit_smoke_v1/validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
