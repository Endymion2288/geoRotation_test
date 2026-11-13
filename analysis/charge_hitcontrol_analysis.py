#!/usr/bin/env python3
"""Summarize controlled native CKF hit filtering and full q/B conjugation."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from validate_charge_hitcontrol import audit, compare_qb_states, compare_qb_tables
from charge_tables import STATUS, derive
from slope_tables import write_json


def read_tables(folder, task):
    result = {}
    for name in ['events', 'truth', 'tracks', 'segments']:
        with np.load(folder/(name+'.npz')) as archive:
            result[name] = {key: archive[key] for key in archive.files}
        result[name]['sample_entry'] = result[name]['source_entry']-task['skip']+task['sample_entry']
        result[name]['file_index'] = np.full(len(result[name]['source_entry']), task['file_index'])
    return result


def seed_map(rows):
    return {(r['source_entry'], r['seed']): (r['bound'], r['cov'])
            for r in rows if r['stage'] == 'seed'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--revision', default='v2')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    replay = args.root/('hitcontrol_replay_'+args.revision)
    jobs = json.loads((replay/'tasks.json').read_text())
    targets = {(r['dsid'], r['sample_entry']): r for r in
               json.loads((args.root/'charge_replay_v1/targets.json').read_text())}
    records, groups, closure = [], defaultdict(Counter), Counter()
    for job in jobs:
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        task = json.loads((folder/'source_task.json').read_text())
        original = Path(task['original_output'])
        target = targets[task['dsid'], task['sample_entry']]
        mode = folder.parent.name
        cohort = 'wrong' if target['wrong'] else 'normal'
        rows = audit(folder)
        baseline_seeds = seed_map(audit(original))
        current_seeds = seed_map(rows)
        if mode == 'reverse_qB':
            table_check = compare_qb_tables(folder, original)
            state_check = compare_qb_states(folder, original)
            closure['qB_events'] += 1
            closure['qB_fields'] += table_check['fields']
            for key, value in state_check.items():
                if isinstance(value, int):
                    closure['qB_'+key] += value
            records.append(dict(dsid=task['dsid'], sample_entry=task['sample_entry'],
                file_index=task['file_index'], mode=mode, original_wrong=target['wrong'],
                all_58_fields_closed=True, **state_check))
            continue
        assert mode == 'truth_filter'
        assert current_seeds == baseline_seeds
        summary = [r for r in rows if r['stage'] == 'truth_source_link_filter']
        assert len(summary) == 1
        barcode = summary[0]['barcode']
        selected = [r for r in rows if r['stage'] == 'selected']
        assert all(dict(r['matches']).get(barcode) == r['measurements'] for r in selected)
        tables = read_tables(folder, task)
        decision = derive(tables, task['dsid'])
        primary = np.flatnonzero(decision['primary'])
        assert len(primary) == 1
        i = primary[0]
        chosen = int(decision['chosen'][i])
        valid = chosen >= 0
        track = tables['tracks']
        status = STATUS[int(decision['status'][i])]
        counts = groups[task['dsid']+'|'+cohort]
        counts[status] += 1
        counts['events'] += 1
        counts['source_links_original'] += summary[0]['original']
        counts['source_links_retained'] += summary[0]['retained']
        counts['seed_means_covariances_exact'] += 1
        records.append(dict(dsid=task['dsid'], sample_entry=task['sample_entry'],
            file_index=task['file_index'], mode=mode, original_wrong=target['wrong'],
            primary_barcode=barcode, original_final_qop_GeVinv=target['final_qop_GeVinv'],
            outcome=status, final_qop_GeVinv=float(track['qop_GeVinv'][chosen]) if valid else None,
            final_purity=float(track['purity'][chosen]) if valid else None,
            chosen_ordinal=int(track['ordinal'][chosen]) if valid else None,
            chosen_hits=list(map(int, track['hit_identifiers'][
                track['hit_offsets'][chosen]:track['hit_offsets'][chosen+1]])) if valid else [],
            seed_means_covariances_exact=True, all_selected_measurements_primary_supported=True,
            original_source_links=summary[0]['original'], retained_source_links=summary[0]['retained']))
        print(task['dsid'], task['sample_entry'], cohort, status, flush=True)
    assert len(records) == len(jobs) == 238
    assert closure['qB_events'] == 119
    write_json(args.output/'events.json', records)
    write_json(args.output/'summary.json', dict(filter_groups={k: dict(v) for k, v in groups.items()},
        qB_closure=dict(closure), all_238_jobs_succeeded=True,
        scope='All 59 final afterGeo wrong signs and 60 frozen normal controls; same input events and seeds',
        control='Actual native combinatorial CKF; only primary-SDO source-link availability changes',
        interpretation='Filter may change branches and winners. Lost and quality-failed tracks are separate from recovered correct signs. Q/B preserves the recorded trajectory and tests sign equivariance, not charge-reflected stochastic G4 physics.'))


if __name__ == '__main__':
    main()
