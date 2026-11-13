#!/usr/bin/env python3
"""Measure full-sample response to the actual Jacobian-matched seed prior."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from charge_analysis import paired_difference, summarize
from charge_tables import STATUS, derive
from slope_tables import merge_tables, write_json


def seed_closure(root, source, dsid, tracks, selected_rows):
    counts = Counter()
    selected_keys = set(zip(map(int, tracks['sample_entry'][selected_rows]),
                            map(int, tracks['ordinal'][selected_rows])))
    for job in json.loads((source/'afterGeo_tasks.json').read_text()):
        output = Path(job['output'])
        task = json.loads((output/'source_task.json').read_text())
        baseline = root/f'main_{dsid}_v1/shards/afterGeo'/f"s{task['sample_entry']//5000:04d}"
        with np.load(baseline/'charge_lineage.npz') as a, np.load(output/'charge_lineage.npz') as b:
            sa, sb = a['stage'] == 0, b['stage'] == 0
            for field in ['source_entry', 'seed', 'qop_GeVinv', 'phi', 'theta', 'loc0', 'loc1']:
                np.testing.assert_array_equal(a[field][sa], b[field][sb], err_msg='Matched prior changed seed mean: '+field)
            np.testing.assert_allclose(b['c44_GeVinv2'][sb], 1e-5, rtol=1e-12, atol=0)
            counts['seeds'] += int(sb.sum())
            final = b['stage'] == 7
            status = b['status'][final]
            counts['native_final'] += len(status)
            counts['refit_success'] += int((status == 1).sum())
            counts['refit_fallback'] += int((status == 0).sum())
            counts['invalid_final'] += int((status < 0).sum())
            entries = b['source_entry'][final]-task['skip']+task['sample_entry']
            chosen = np.array([(int(e), int(o)) in selected_keys
                for e, o in zip(entries, b['ordinal'][final])], dtype=bool)
            counts['selected_refit_success'] += int((chosen & (status == 1)).sum())
            counts['selected_refit_fallback'] += int((chosen & (status == 0)).sum())
            counts['selected_invalid_final'] += int((chosen & (status < 0)).sum())
    assert counts['selected_refit_success']+counts['selected_refit_fallback'] == len(selected_rows)
    return dict(counts, all_seed_means_exact=True, covariance_qop_GeVinv2=1e-5)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    groups, coverage, migrations, seeds = {}, {}, {}, {}
    for dsid in ['100116', '100117']:
        source = args.root/'matched_prior_500k_v1'/dsid
        tables, coverage[dsid] = merge_tables(source, 'afterGeo')
        d = derive(tables, dsid)
        groups[dsid] = summarize(tables, d)
        primary = d['primary']
        current = {key: d[key][primary] for key in ['status', 'selected', 'wrong', 'fiducial']}
        current.update(sample_entry=tables['truth']['sample_entry'][primary],
            barcode=tables['truth']['barcode'][primary], truth_p_GeV=d['truth_p'][primary])
        rows = d['chosen'][primary]
        values = np.full(len(rows), np.nan)
        values[rows >= 0] = tables['tracks']['qop_GeVinv'][rows[rows >= 0]]
        current['qop_GeVinv'] = values
        np.savez_compressed(args.output/(dsid+'_outcomes.npz'), **current)
        migrations[dsid] = {}
        for variant in ['afterGeo', 'beforeGeo']:
            with np.load(args.root/f'main_analysis_v1/{dsid}_{variant}_outcomes.npz') as archive:
                baseline = {key: archive[key] for key in ['status', 'selected', 'wrong']}
                for field in ['sample_entry', 'barcode', 'truth_p_GeV', 'fiducial']:
                    np.testing.assert_array_equal(archive[field], current[field], err_msg=field)
            matrix = np.bincount(7*baseline['status']+current['status'], minlength=49).reshape(7, 7)
            migrations[dsid][variant] = dict(states=STATUS, counts=matrix.tolist(),
                rate_difference=paired_difference(baseline, current, np.ones(len(rows), dtype=bool)))
        seeds[dsid] = seed_closure(args.root, source, dsid, tables['tracks'], d['chosen'][d['selected']])
        print(dsid, groups[dsid]['conditional_misid'], flush=True)
        del tables, d
    write_json(args.output/'summary.json', dict(groups=groups, coverage=coverage,
        migrations=migrations, seed_validation=seeds,
        prior='beforeGeo full seed covariance mapped by actual local/angular Jacobian; J C_F J^T',
        nonlinear_scope='Matches the first-order covariance, not the finite-width physical spherical prior distribution',
        reconstruction='Actual geoRotation Calypso; native means, field, geometry, material, finding and refit; isolated smoothing ON'))


if __name__ == '__main__':
    main()
