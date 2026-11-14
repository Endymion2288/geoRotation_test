#!/usr/bin/env python3
"""Compare charge-crossing operations in the fixed legacy/corrected cohort."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from charge_centroid_analysis import measurements, read_tables
from charge_lineage_analysis import classify_path
from charge_replay_analysis import operation, read_audit
from charge_tables import STATUS, derive
from slope_tables import write_json


def target_path(folder, task, hit_truth):
    records = read_audit(folder/'charge_audit.jsonl.gz')
    tables = read_tables(folder, task)
    d = derive(tables, task['dsid'])
    indices = np.flatnonzero(d['primary'])
    assert len(indices) == 1
    i = int(indices[0])
    chosen = int(d['chosen'][i])
    assert chosen >= 0
    tracks = tables['tracks']
    offsets = tracks['hit_offsets']
    hits = tracks['hit_identifiers'][offsets[chosen]:offsets[chosen+1]].tolist()
    final_qop = float(tracks['qop_GeVinv'][chosen])
    finals = [r for r in records if r['stage'] == 'final' and r['hits'] == hits
              and np.isclose(1000*r['bound_trk'][4], final_qop, rtol=1e-12, atol=0.)]
    assert len(finals) == 1, (str(folder), len(finals))
    final = finals[0]
    candidate, seed = final['candidate'], final['seed']

    def one(stage):
        rows = [r for r in records if r['stage'] == stage and r['seed'] == seed
                and (stage == 'seed' or r['candidate'] == candidate)]
        assert len(rows) == 1, (str(folder), stage, len(rows))
        return rows[0]

    start, trace = one('seed'), one('selected_state_trace')
    selected, converted = one('selected'), one('converted')
    np.testing.assert_allclose(selected['bound'][4], 1000*converted['bound_trk'][4], rtol=1e-12, atol=0.)
    charge = d['charge']
    primary = int(tables['truth']['barcode'][i])
    qops, stages, updates = [start['bound'][4]], [1], []
    for step in trace['steps']:
        if not step['measurement_flag']:
            continue
        hit = int(step['hit'])
        update = operation(step, hit_truth.get(hit, []), primary, charge)
        if hit not in hit_truth:
            update['hit_truth'] = 'not_in_legacy_measurement_audit'
        updates.append(update)
        stages.extend([2, 3])
        qops.extend([step['predicted'][4], step['filtered'][4]])
    stages.extend([4, 5, 6])
    qops.extend([selected['bound'][4], 1000*converted['bound_trk'][4], final_qop])
    first, last = classify_path(stages, qops, charge)
    bad = [r for r in updates if r['correct_to_wrong']]
    if last == 3:
        assert bad
    return dict(outcome=STATUS[int(d['status'][i])], candidate=candidate, seed=seed,
        ordinal=int(tracks['ordinal'][chosen]), primary_barcode=primary,
        seed_qop_GeVinv=float(start['bound'][4]), native_reference_qop_GeVinv=float(selected['bound'][4]),
        final_qop_GeVinv=final_qop, final_purity=float(tracks['purity'][chosen]),
        first_wrong=first, final_wrong_entry=last,
        last_update_into_wrong=bad[-1] if bad else None,
        seed_bound=start['bound'], seed_covariance=start['cov'], seed_surface=start['surface'],
        seed_hits=start['hits'], ordered_measurement_hits=[u['hit_id'] for u in updates],
        per_update=updates, audit=str(folder/'charge_audit.jsonl.gz')), records


def reduce(root, replay, output, jobs_path=None):
    jobs = json.loads((jobs_path or replay/'tasks.json').read_text())
    if jobs_path is None:
        assert json.loads((replay/'full_validation.json').read_text())['all_exact']
    grouped = defaultdict(list)
    for job in jobs:
        folder = Path(job['output'])
        validation = json.loads((folder/'replay_validation.json').read_text())
        assert validation['all_reconstruction_arrays_identical'] and validation['fields_exact'] == 58
        task = json.loads((folder/'source_task.json').read_text())
        original = json.loads((Path(task['legacy_replay_output'])/'source_task.json').read_text())
        grouped[original['original_output']].append((folder, task))
    events, groups = [], defaultdict(Counter)
    for original, tasks in grouped.items():
        entries = {task['skip'] for _, task in tasks}
        support = {entry: {} for entry in entries}
        with np.load(Path(original)/'charge_lineage.npz') as archive:
            mask = (archive['stage'] == 8) & np.isin(archive['source_entry'], list(entries))
            for entry, hit, barcodes in zip(archive['source_entry'][mask], archive['hit_id'][mask],
                                           archive['truth_barcodes'][mask]):
                truth = [int(b) for b in barcodes.split('|') if b]
                previous = support[int(entry)].setdefault(int(hit), truth)
                assert previous == truth
        for folder, task in tasks:
            old, a = target_path(Path(task['legacy_replay_output']), task, support[task['skip']])
            new, b = target_path(folder, task, support[task['skip']])
            cohort = 'wrong' if task['replay_targets'][0]['wrong'] else 'normal'
            assert (old['outcome'] == 'wrong') == (cohort == 'wrong')
            same_path = old['ordered_measurement_hits'] == new['ordered_measurement_hits']
            old_m, new_m = measurements(a), measurements(b)
            common = old_m.keys() & new_m.keys()
            assert all(old_m[hit] == new_m[hit] for hit in common)
            seed_equal = all(old[k] == new[k] for k in
                             ['seed', 'seed_bound', 'seed_covariance', 'seed_surface', 'seed_hits'])
            old_last, new_last = old['last_update_into_wrong'], new['last_update_into_wrong']
            same_last_hit = old_last is not None and new_last is not None and old_last['hit_id'] == new_last['hit_id']
            row = dict(dsid=task['dsid'], sample_entry=task['sample_entry'], file_index=task['file_index'],
                source_entry=task['skip'], cohort=cohort, legacy=old, corrected=new,
                same_selected_ordered_hits=same_path, same_complete_seed=seed_equal,
                common_audited_measurements=len(common), common_measurements_exact=True,
                same_last_crossing_hit=same_last_hit)
            events.append(row)
            group = groups[task['dsid']+'|'+cohort]
            group['events'] += 1
            group['same_ordered_hits'] += same_path
            group['same_complete_seed'] += seed_equal
            group['same_first_wrong_stage'] += old['first_wrong'] == new['first_wrong']
            group['same_final_wrong_entry'] += old['final_wrong_entry'] == new['final_wrong_entry']
            group['same_last_crossing_hit'] += same_last_hit
            group['final_qop_changed'] += old['final_qop_GeVinv'] != new['final_qop_GeVinv']
            for label, value in [('legacy', old), ('corrected', new)]:
                group[label+'_outcome_'+value['outcome']] += 1
                if cohort == 'wrong':
                    group[label+'_last_entry_'+str(value['final_wrong_entry'])] += 1
                    if value['final_wrong_entry'] == 3:
                        group[label+'_crossing_'+value['last_update_into_wrong']['hit_truth']] += 1
            print(task['dsid'], task['sample_entry'], cohort, old['outcome'], new['outcome'], flush=True)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/'events.json', events)
    write_json(output/'summary.json', dict(events=len(events), groups={k: dict(v) for k, v in groups.items()},
        all_current_reconstruction_fields_exact=True, all_common_measurements_exact=True,
        support_source='Original unmodified RDO SDO support retained in legacy compact measurement audit; missing coverage is explicitly separate.',
        scope='Frozen original 59 wrong signs and 60 normal controls; material encoding only, isolated smoothing ON.',
        interpretation='Recompute actual matrix updates and candidate mapping in each encoding; no legacy numerical q/p is imposed on corrected traces.'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs', type=Path, help='Two-event smoke analysis before the full validator')
    args = parser.parse_args()
    reduce(args.root, args.root/'material_encoding_control_v1/replay_v1', args.output, args.jobs)


if __name__ == '__main__':
    main()
