#!/usr/bin/env python3
"""Measure native Calypso charge response to actual SiHit-centroid clusters."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import re

import numpy as np

from charge_replay_analysis import read_audit
from charge_tables import STATUS, derive
from slope_tables import write_json


def read_tables(folder, task):
    tables = {}
    for name in ['events', 'truth', 'tracks']:
        with np.load(folder/(name+'.npz')) as archive:
            tables[name] = {key: archive[key] for key in archive.files}
        tables[name]['sample_entry'] = tables[name]['source_entry']-task['skip']+task['sample_entry']
        tables[name]['file_index'] = np.full(len(tables[name]['source_entry']), task['file_index'])
    return tables


def compare_seeds(original, current):
    old = {r['seed']: r for r in original if r['stage'] == 'seed'}
    new = {r['seed']: r for r in current if r['stage'] == 'seed'}
    old_hits = Counter(tuple(sorted(row['hits'])) for row in old.values())
    new_hits = Counter(tuple(sorted(row['hits'])) for row in new.values())
    result = dict(original_seeds=len(old), centroid_seeds=len(new),
        same_seed_hit_sets=sum((old_hits & new_hits).values()),
        removed_seed_hit_sets=sum((old_hits-new_hits).values()),
        added_seed_hit_sets=sum((new_hits-old_hits).values()),
        same_id_hits_surface=0, mean_changed=0, covariance_changed=0)
    for seed in old.keys() & new.keys():
        a, b = old[seed], new[seed]
        if sorted(a['hits']) != sorted(b['hits']) or a['surface'] != b['surface']:
            continue
        result['same_id_hits_surface'] += 1
        result['mean_changed'] += a['bound'] != b['bound']
        result['covariance_changed'] += a['cov'] != b['cov']
    return result


def measurements(rows):
    values = {}
    for row in rows:
        if row['stage'] != 'selected_state_trace':
            continue
        for state in row['steps']:
            if not state['measurement_flag']:
                continue
            hit = int(state['hit'])
            value = state['measurement']
            if hit in values:
                assert values[hit] == value
            values[hit] = value
    return values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    targets = {(row['dsid'], row['sample_entry']): row for row in
               json.loads((args.root/'charge_replay_v1/targets.json').read_text())}
    jobs = json.loads((args.root/'centroid_replay_v1/tasks.json').read_text())
    expected = {str((args.root.parent.parent/'stage05_slopes/software_cluster_control_v2/lib'/name).resolve())
                for name in ['libFaserSiClusterizationTool.so', 'libFaserSiClusterizationToolLib.so']}
    events, groups = [], defaultdict(Counter)
    for job in jobs:
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        task = json.loads((folder/'source_task.json').read_text())
        run = json.loads((folder/'run.json').read_text())
        assert run['success'] and run['isolate_smoothing'] and run['args']['truth_cluster_position']
        assert set(run['loaded_cluster_libraries']) == expected
        target = targets[task['dsid'], task['sample_entry']]
        mode = folder.parent.name
        assert mode in ['centroid', 'centroid_primary_filter']
        cohort = 'wrong' if target['wrong'] else 'normal'
        original = read_audit(Path(task['original_output'])/'charge_audit.jsonl.gz')
        current = read_audit(folder/'charge_audit.jsonl.gz')
        old_m, new_m = measurements(original), measurements(current)
        common = old_m.keys() & new_m.keys()
        delta = [float(np.max(np.abs(np.asarray(old_m[hit])-new_m[hit]))) for hit in common]
        tables = read_tables(folder, task)
        decision = derive(tables, task['dsid'])
        primary = np.flatnonzero(decision['primary'])
        assert len(primary) == 1
        i = primary[0]
        chosen = int(decision['chosen'][i])
        status = STATUS[int(decision['status'][i])]
        track = tables['tracks']
        if mode == 'centroid_primary_filter':
            barcode = int(tables['truth']['barcode'][i])
            assert run['truth_source_link_barcode'] == barcode
            assert all(dict(row['matches']).get(barcode) == row['measurements']
                       for row in current if row['stage'] == 'selected')
        log = gzip.open(folder/'tracking.log.gz', 'rt').read()
        replacement = re.findall(r'Truth cluster position: enabled=.*?True\s+replaced=(\d+) fallback=(\d+)', log)
        assert len(replacement) == 1
        replaced, fallback = map(int, replacement[0])
        seeds = compare_seeds(original, current)
        group = groups['|'.join([task['dsid'], mode, cohort])]
        group['events'] += 1
        group[status] += 1
        group['replaced_clusters'] += replaced
        group['fallback_clusters'] += fallback
        for key, value in seeds.items():
            group[key] += value
        events.append(dict(dsid=task['dsid'], sample_entry=task['sample_entry'],
            file_index=task['file_index'], mode=mode, original_wrong=target['wrong'],
            original_final_purity=target['purity'], original_final_qop_GeVinv=target['final_qop_GeVinv'],
            outcome=status, final_qop_GeVinv=float(track['qop_GeVinv'][chosen]) if chosen >= 0 else None,
            final_purity=float(track['purity'][chosen]) if chosen >= 0 else None,
            chosen_ordinal=int(track['ordinal'][chosen]) if chosen >= 0 else None,
            seed_response=seeds, replaced_clusters=replaced, fallback_clusters=fallback,
            common_audited_measurement_hits=len(common),
            changed_common_measurements=sum(value > 0 for value in delta),
            maximum_common_measurement_change_mm=max(delta, default=None),
            audit=str(folder/'charge_audit.jsonl.gz')))
        print(task['dsid'], task['sample_entry'], mode, cohort, status, flush=True)
    assert len(events) == len(jobs) == 238
    write_json(args.output/'events.json', events)
    write_json(args.output/'summary.json', dict(groups={key: dict(value) for key, value in groups.items()},
        jobs=238, all_jobs_exit_zero=True, all_actual_cluster_libraries_match=True,
        scope='All 59 afterGeo wrong signs and 60 frozen normal controls; two centroid modes; isolated smoothing ON',
        control='Actual Calypso cluster maker replaces measured locX with same-wafer majority-SiHit energy-weighted midpoint; original width/error/noise/IDs and missing-truth fallback retained',
        interpretation='Measurement changes can alter segments, seeds, CKF branches and final selection. Quality failures remain separate from recovered correct signs; this is not a fixed-seed or digitization-off experiment.'))


if __name__ == '__main__':
    main()
