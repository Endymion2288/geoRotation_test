#!/usr/bin/env python3
"""Validate native replay and native decoded material objects before physics jobs."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--native', default='native_smoke_v3')
    parser.add_argument('--audit', default='audit_smoke_v4')
    args = parser.parse_args()
    root = args.root
    native = root/'runs'/args.native
    audit = root/'runs'/args.audit
    a, b = np.load(native/'writer.npz'), np.load(audit/'writer.npz')
    for field in a.files:
        np.testing.assert_equal(a[field], b[field], err_msg=field)
    material = json.loads((root/'material-latest-rowmajor.json').read_text())
    entries = material['Surfaces']['entries']
    queries = []
    with gzip.open(audit/'audit.jsonl.gz', 'rt') as stream:
        for line in stream:
            row = json.loads(line)
            if 'material_bin_transform' not in row:
                continue
            gid = row['geometry_id']
            ids = dict(volume=gid>>56, layer=(gid>>36)&0xfff,
                       approach=(gid>>28)&0xff, sensitive=(gid>>8)&0xfffff)
            candidates = [e for e in entries if all(ids[k]==v for k,v in e.items() if k!='value')]
            assert len(candidates) == 1, ids
            slab = candidates[0]['value']['material']
            transform = slab['binUtility']['transform']
            matrix = np.asarray(row['material_bin_transform'])
            np.testing.assert_allclose(matrix[:3,:3], np.asarray(transform['rotation']).reshape(3,3), atol=1e-15)
            np.testing.assert_allclose(matrix[:3,3], transform['translation'], atol=1e-12)
            bins = row['material_bin_indices']
            raw = slab['data'][bins[1]][bins[0]]
            np.testing.assert_allclose(row['material_full_thickness_mm'], raw['thickness'], rtol=1e-6)
            if raw['material'] is not None:
                np.testing.assert_allclose(row['material_parameters'], raw['material'], rtol=1e-6)
            queries.append(dict(geometry_id=gid, native_matrix=matrix.tolist(),
                native_local=row['material_bin_local'], native_indices=bins,
                native_raw_thickness_mm=row['material_full_thickness_mm'],
                json_raw_thickness_mm=raw['thickness']))
    assert queries, 'No native binned-material query was observed'
    states = np.load(audit/'states.npz')
    values = states['values']
    fields = list(states['fields'])
    errors = values[:,fields.index('update_identity_error')]
    finite = errors[np.isfinite(errors)]
    assert len(finite)>0 and np.max(np.abs(finite))<1e-10
    result = dict(native=str(native), audit=str(audit),
        replay_fields=a.files, exact_equal_nan=True,
        native_material_queries=len(queries), material_queries=queries,
        update_identity_count=len(finite), update_identity_max_abs=float(np.max(np.abs(finite))))
    (root/'smoke_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k!='material_queries'})


if __name__ == '__main__':
    main()
