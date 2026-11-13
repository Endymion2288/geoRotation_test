#!/usr/bin/env python3
"""Measure centroid response with the original complete native Calypso seeds."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from charge_centroid_analysis import measurements, read_tables
from charge_replay_analysis import read_audit, operation
from charge_tables import STATUS, derive
from slope_tables import write_json


def quality_record(tracks, chosen):
    if chosen < 0:
        return None
    return {key: tracks[key][chosen].item() for key in
            ['chi2', 'ndof', 'rmax_mm', 'n_hits', 'n_layers', 'station_bits', 'purity']}


def failed_quality(record):
    if record is None:
        return ['no_matched_track']
    return [name for name, passed in [
        ('three_downstream_stations', record['station_bits'] & 14 == 14),
        ('positive_ndof', record['ndof'] > 0),
        ('chi2_ndof_below_25', record['ndof'] > 0 and record['chi2']/record['ndof'] < 25),
        ('at_least_seven_layers', record['n_layers'] >= 7),
        ('ndof_plus_five_at_least_twelve', record['ndof']+5 >= 12),
        ('rmax_below_95_mm', record['rmax_mm'] < 95)] if not passed]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    valid = json.loads((args.root/'seed_mean_smoke_v1/validation.json').read_text())
    assert valid['all_identical_58_fields_exact'] and valid['all_full_mean_overrides_exact']
    jobs = json.loads((args.root/'centroid_fixed_seed_v1/tasks.json').read_text())
    targets = {(row['dsid'], row['sample_entry']): row for row in
               json.loads((args.root/'charge_replay_v1/targets.json').read_text())}
    free = {(row['dsid'], row['sample_entry']): row for row in
            json.loads((args.root/'centroid_analysis_v1/events.json').read_text()) if row['mode'] == 'centroid'}
    expected = {str((args.root.parent.parent/'stage05_slopes/software_cluster_control_v2/lib'/name).resolve())
                for name in ['libFaserSiClusterizationTool.so', 'libFaserSiClusterizationToolLib.so']}
    events, groups = [], defaultdict(Counter)
    for job in jobs:
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        task = json.loads((folder/'source_task.json').read_text())
        run = json.loads((folder/'run.json').read_text())
        assert run['success'] and run['isolate_smoothing'] and run['args']['truth_cluster_position']
        assert run['seed_override_mode'] == 'mean'
        assert set(run['loaded_cluster_libraries']) == expected
        key = task['dsid'], task['sample_entry']
        target, centroid = targets[key], free[key]
        original = read_audit(Path(task['original_output'])/'charge_audit.jsonl.gz')
        current = read_audit(folder/'charge_audit.jsonl.gz')
        old_seeds = {row['seed']: row for row in original if row['stage'] == 'seed'}
        new_seeds = {row['seed']: row for row in current if row['stage'] == 'seed'}
        assert old_seeds.keys() == new_seeds.keys()
        for seed in old_seeds:
            assert Counter(old_seeds[seed]['hits']) == Counter(new_seeds[seed]['hits'])
            for field in ['surface', 'cov', 'bound']:
                assert old_seeds[seed][field] == new_seeds[seed][field]
        seed_hit_order_changed = sum(old_seeds[seed]['hits'] != new_seeds[seed]['hits'] for seed in old_seeds)
        a, b = measurements(original), measurements(current)
        common = a.keys() & b.keys()
        delta = [float(np.max(np.abs(np.asarray(a[hit])-b[hit]))) for hit in common]
        tables = read_tables(folder, task)
        d = derive(tables, task['dsid'])
        primary = np.flatnonzero(d['primary'])
        assert len(primary) == 1
        chosen = int(d['chosen'][primary[0]])
        status = STATUS[int(d['status'][primary[0]])]
        quality = quality_record(tables['tracks'], chosen)
        cohort = 'pure_wrong' if target['wrong'] else 'normal'
        group = groups[task['dsid']+'|'+cohort]
        group['events'] += 1
        group[status] += 1
        group['native_seeds_exact'] += len(old_seeds)
        group['seed_hit_order_changed'] += seed_hit_order_changed
        row = dict(dsid=task['dsid'], sample_entry=task['sample_entry'], cohort=cohort,
            outcome=status, free_seed_centroid_outcome=centroid['outcome'],
            native_final_qop_GeVinv=target['final_qop_GeVinv'],
            free_seed_centroid_qop_GeVinv=centroid['final_qop_GeVinv'],
            final_qop_GeVinv=float(tables['tracks']['qop_GeVinv'][chosen]) if chosen >= 0 else None,
            quality=quality, failed_quality_conditions=failed_quality(quality),
            native_seed_means_full_covariance_surfaces_exact=True,
            native_seed_hit_multisets_exact=True, seed_hit_order_changed=seed_hit_order_changed,
            common_audited_measurements=len(common), changed_common_measurements=sum(value > 0 for value in delta),
            maximum_common_measurement_change_mm=max(delta, default=None),
            audit=str(folder/'charge_audit.jsonl.gz'))
        if target['wrong']:
            assert target['purity'] == 1 and len(old_seeds) == 1
            traces = [record for record in current if record['stage'] == 'selected_state_trace']
            assert len(traces) == 1
            charge = 1 if task['dsid'] == '100116' else -1
            # Matrix identity is checked for every actual measurement update.
            # SDO classification belongs to the original frozen replay table.
            updates = [operation(step, [], 0, charge) for step in traces[0]['steps'] if step['measurement_flag']]
            for update in updates:
                update.pop('truth_barcodes')
                update.pop('hit_truth')
            old_trace = [record for record in original if record['stage'] == 'selected_state_trace']
            assert len(old_trace) == 1
            old_hits = [step['hit'] for step in old_trace[0]['steps'] if step['measurement_flag']]
            new_hits = [step['hit'] for step in traces[0]['steps'] if step['measurement_flag']]
            row.update(same_selected_ordered_measurement_hits=old_hits == new_hits,
                       native_selected_measurement_hits=old_hits, selected_measurement_hits=new_hits,
                       per_update=updates)
        events.append(row)
        print(key, cohort, status, flush=True)
    assert len(events) == len(jobs) == 66
    write_json(args.output/'events.json', events)
    write_json(args.output/'summary.json', dict(groups={key: dict(value) for key, value in groups.items()},
        jobs=66, all_jobs_exit_zero=True, all_full_native_seeds_exact=True,
        all_actual_cluster_libraries_match=True,
        scope='Six native afterGeo purity-one wrong signs plus 60 frozen normal controls; isolated smoothing ON',
        control='Actual Calypso SiHit-centroid measurements with full original native seed means, covariance, surfaces and hits',
        interpretation='The full numerical seed means/covariances/surfaces and constituent hit multisets are exact; seed hit-vector order can differ. Native CKF still chooses branches and hits; selected path equality and quality failures are reported separately.'))


if __name__ == '__main__':
    main()
