#!/usr/bin/env python3
"""Validate real Calypso default closure, seed isolation and q/B conjugation."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path

import numpy as np

from validate_charge_audit import compare_tables


def audit(folder):
    with gzip.open(folder/'charge_audit.jsonl.gz', 'rt') as source:
        return [json.loads(line) for line in source]


def compare_qb_tables(current, original):
    fields, maximum = 0, 0.
    for name in ['events', 'truth', 'tracks', 'segments']:
        with np.load(current/(name+'.npz')) as a, np.load(original/(name+'.npz')) as b:
            assert a.files == b.files
            for field in a.files:
                expected = b[field]
                if name == 'tracks' and field in ['charge', 'qop_GeVinv']:
                    expected = -expected
                if expected.dtype.kind == 'f':
                    np.testing.assert_allclose(a[field], expected, rtol=1e-10, atol=1e-10,
                                               err_msg=name+'/'+field)
                    finite = np.isfinite(expected) & np.isfinite(a[field])
                    if finite.any():
                        maximum = max(maximum, float(np.max(np.abs(a[field][finite]-expected[finite]))))
                else:
                    np.testing.assert_array_equal(a[field], expected, err_msg=name+'/'+field)
                fields += 1
    return dict(fields=fields, maximum_absolute_float_difference=maximum,
                exact_event_hit_candidate_counts=True)


def compare_qb_states(current, original):
    def selected(rows):
        return {(r['source_entry'], r['stage'], r['seed'], r['candidate']): r for r in rows
                if r['stage'] in ['seed', 'selected_state_trace', 'refit_state_trace']}
    a, b = selected(audit(current)), selected(audit(original))
    assert a.keys() == b.keys()
    t = np.diag([1., 1., 1., 1., -1., 1.])
    counts = Counter()
    for key, expected in b.items():
        observed = a[key]
        if key[1] == 'seed':
            np.testing.assert_allclose(observed['bound'], t@expected['bound'], rtol=1e-11, atol=1e-11)
            np.testing.assert_allclose(observed['cov'], t@expected['cov']@t, rtol=1e-11, atol=1e-11)
            counts['seed_means_and_full_covariance'] += 1
            continue
        assert len(observed['steps']) == len(expected['steps'])
        for left, right in zip(observed['steps'], expected['steps']):
            for field in ['hit', 'geometry_id', 'measurement_flag', 'outlier_flag', 'hole_flag']:
                assert left.get(field) == right.get(field), (key, field)
            for stage in ['predicted', 'filtered', 'smoothed']:
                if stage not in right:
                    continue
                np.testing.assert_allclose(left[stage], t@right[stage], rtol=1e-10, atol=1e-10)
                np.testing.assert_allclose(left[stage+'_cov'], t@right[stage+'_cov']@t,
                                           rtol=1e-9, atol=1e-9)
                counts[stage+'_mean_and_full_covariance'] += 1
            counts['same_path_states'] += 1
    return dict(counts, all_paths_and_covariance_signs_closed=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--revision', default='v1')
    args = parser.parse_args()
    root = args.root/('hitcontrol_smoke_'+args.revision)
    result = {}
    for dsid in ['100044', '100116', '100117']:
        folders = {mode: root/dsid/mode/'f00_e00000000' for mode in ['default', 'truth_filter', 'reverse_qB']}
        for folder in folders.values():
            assert (folder/'exit_code.txt').read_text().strip() == '0'
        baseline = args.root/'charge_smoke/afterGeo'/dsid/'off'
        default = compare_tables(folders['default'], baseline)
        qb = compare_qb_tables(folders['reverse_qB'], folders['default'])
        states = compare_qb_states(folders['reverse_qB'], folders['default'])
        def seeds(folder):
            return {(r['source_entry'], r['seed']): (r['bound'], r['cov'])
                    for r in audit(folder) if r['stage'] == 'seed'}
        assert seeds(folders['truth_filter']) == seeds(folders['default'])
        filtered = audit(folders['truth_filter'])
        summaries = [r for r in filtered if r['stage'] == 'truth_source_link_filter']
        assert len(summaries) == 20
        for row in filtered:
            if row['stage'] != 'selected':
                continue
            assert dict(row['matches']).get(summaries[0]['barcode']) == row['measurements']
        result[dsid] = dict(default_fields_exact=len(default), qB_tables=qb, qB_states=states,
            filter_seed_means_and_covariances_exact=True,
            source_links_original=sum(r['original'] for r in summaries),
            source_links_retained=sum(r['retained'] for r in summaries))
    summary = dict(samples=result, all_default_fields_exact=True, all_qB_paths_closed=True,
        scope='20 events per DSID; actual native CKF/refit/propagator and full SDO-supported source-link filter')
    (root/'validation.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
