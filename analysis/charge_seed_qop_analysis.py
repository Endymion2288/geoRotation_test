#!/usr/bin/env python3
"""Separate native seed curvature means from cluster-coordinate controls."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from charge_centroid_analysis import read_tables
from charge_replay_analysis import read_audit
from charge_tables import STATUS, derive
from slope_tables import write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    jobs = json.loads((args.root/'seed_qop_replay_v1/tasks.json').read_text())
    targets = {(row['dsid'], row['sample_entry']): row for row in
               json.loads((args.root/'charge_replay_v1/targets.json').read_text())}
    events, groups = [], defaultdict(Counter)
    for job in jobs:
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        task = json.loads((folder/'source_task.json').read_text())
        run = json.loads((folder/'run.json').read_text())
        assert run['success'] and run['isolate_smoothing'] and run['seed_override_mode'] == 'qop'
        target = targets[task['dsid'], task['sample_entry']]
        original = read_audit(Path(task['original_output'])/'charge_audit.jsonl.gz')
        current = read_audit(folder/'charge_audit.jsonl.gz')
        a = {row['seed']: row for row in original if row['stage'] == 'seed'}
        b = {row['seed']: row for row in current if row['stage'] == 'seed'}
        assert a.keys() == b.keys()
        tables = read_tables(folder, task)
        decision = derive(tables, task['dsid'])
        primary = np.flatnonzero(decision['primary'])
        assert len(primary) == 1
        i = primary[0]
        expected_qop = decision['charge']/decision['truth_p'][i]
        for seed in a:
            for field in ['hits', 'surface', 'cov']:
                assert a[seed][field] == b[seed][field]
            np.testing.assert_array_equal(np.delete(a[seed]['bound'], 4), np.delete(b[seed]['bound'], 4))
            np.testing.assert_allclose(b[seed]['bound'][4], expected_qop, rtol=1e-14, atol=0.)
        chosen = int(decision['chosen'][i])
        outcome = STATUS[int(decision['status'][i])]
        cohort = 'wrong' if target['wrong'] else 'normal'
        group = groups[task['dsid']+'|'+cohort]
        group['events'] += 1
        group[outcome] += 1
        group['seeds'] += len(a)
        group['qop_changed'] += sum(a[key]['bound'][4] != b[key]['bound'][4] for key in a)
        tracks = tables['tracks']
        events.append(dict(dsid=task['dsid'], sample_entry=task['sample_entry'],
            original_wrong=target['wrong'], original_purity=target['purity'],
            original_final_qop_GeVinv=target['final_qop_GeVinv'],
            outcome=outcome, final_qop_GeVinv=float(tracks['qop_GeVinv'][chosen]) if chosen >= 0 else None,
            final_variance_GeVinv2=float(tracks['qop_variance_GeVinv2'][chosen]) if chosen >= 0 else None,
            final_purity=float(tracks['purity'][chosen]) if chosen >= 0 else None,
            native_seed_qop_GeVinv={str(key): a[key]['bound'][4] for key in a},
            replacement_qop_GeVinv=float(expected_qop),
            all_other_seed_means_covariances_hits_exact=True,
            audit=str(folder/'charge_audit.jsonl.gz')))
        print(task['dsid'], task['sample_entry'], cohort, outcome, flush=True)
    assert len(events) == len(jobs) == 119
    write_json(args.output/'events.json', events)
    write_json(args.output/'summary.json', dict(groups={key: dict(value) for key, value in groups.items()},
        jobs=119, all_jobs_exit_zero=True, all_changes_qop_mean_only=True,
        scope='All 59 afterGeo wrong signs and 60 frozen normal controls; actual CKF SeedOverrideMode=qop; isolated smoothing ON',
        reference='Same primary stored production charge / momentum, not sensor local G4 momentum',
        interpretation='Native cluster coordinates and seed locations/directions/full covariance/hit sets retained. CKF branches and winners may change in response to the mean curvature. Quality failure is separate from restored charge.'))


if __name__ == '__main__':
    main()
