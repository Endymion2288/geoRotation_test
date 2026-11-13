#!/usr/bin/env python3
"""Measure refit-only charge response with identical native CKF rows."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from charge_analysis import paired_difference, summarize
from charge_modes import compact
from charge_tables import STATUS, derive
from slope_tables import merge_tables, write_json


def check_lineage(source, tracks, selected_rows, compare=False):
    selected = set(zip(map(int, tracks['sample_entry'][selected_rows]),
                       map(int, tracks['ordinal'][selected_rows])))
    result = Counter()
    for job in json.loads((source/'beforeGeo_tasks.json').read_text()):
        folder = Path(job['output'])
        task = json.loads((folder/'source_task.json').read_text())
        with np.load(folder/'charge_lineage.npz') as current:
            if compare:
                with np.load(Path(task['original_output'])/'charge_lineage.npz') as original:
                    a, b = current['stage'] != 7, original['stage'] != 7
                    for field in current.files:
                        np.testing.assert_array_equal(current[field][a], original[field][b], err_msg=field)
                    result['all_CKF_rows_exact'] += int(a.sum())
            final = current['stage'] == 7
            entries = current['source_entry'][final]-task['skip']+task['sample_entry']
            chosen = np.array([(int(e), int(o)) in selected
                for e, o in zip(entries, current['ordinal'][final])])
            status = current['status'][final]
            qop = current['qop_GeVinv'][final]
            converted = current['qop_GeVinv'][current['stage'] == 6]
            assert len(converted) == len(qop)
            charge = 1 if task['dsid'] == '100116' else -1
            for label, mask in [('all_final', np.ones(len(qop), dtype=bool)), ('selected', chosen)]:
                result[label+'_total'] += int(mask.sum())
                result[label+'_refit_success'] += int((mask & (status == 1)).sum())
                result[label+'_fallback'] += int((mask & (status == 0)).sum())
                result[label+'_correct_to_wrong_in_refit'] += int((mask & (charge*converted > 0) & (charge*qop < 0)).sum())
                result[label+'_wrong_to_correct_in_refit'] += int((mask & (charge*converted < 0) & (charge*qop > 0)).sum())
    assert result['selected_total'] == len(selected_rows)
    return dict(result)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    groups, coverage, counts, migrations = {}, {}, {}, {}
    for dsid in ['100116', '100117']:
        baseline = None
        for mode in ['native', 'smoother', 'updater', 'both']:
            source = (args.root/'mode_analysis_v1/forward_references'/dsid if mode == 'native'
                      else args.root/'refit_control_100k_v1'/dsid/mode)
            tables, cov = merge_tables(source, 'beforeGeo')
            d = derive(tables, dsid)
            key = dsid+'|'+mode
            groups[key] = summarize(tables, d)
            coverage[key] = cov
            counts[key] = check_lineage(source, tables['tracks'], d['chosen'][d['selected']], mode != 'native')
            current = compact(tables, d)
            np.savez_compressed(args.output/(key.replace('|', '_')+'_outcomes.npz'), **current)
            if baseline is None:
                baseline = current
            for field in ['sample_entry', 'barcode', 'truth_p_GeV', 'fiducial']:
                np.testing.assert_array_equal(current[field], baseline[field])
            migrations[key] = dict(states=STATUS, counts=np.bincount(
                7*baseline['status']+current['status'], minlength=49).reshape(7, 7).tolist(),
                conditional_difference=paired_difference(baseline, current, np.ones(len(current['status']), dtype=bool)))
            print(key, groups[key]['conditional_misid'], counts[key], flush=True)
            del tables, d
    write_json(args.output/'summary.json', dict(groups=groups, coverage=coverage, refit_counts=counts,
        migrations=migrations, source_events_per_charge_mode=100000, all_nonfinal_CKF_rows_exact=True,
        scope='Same full input coverage as mode control: first 20,000 of each of five files',
        control='Actual beforeGeo Calypso refit only: regular smoother, regular updater, both; candidate-isolated CKF smoothing stays ON'))


if __name__ == '__main__':
    main()
