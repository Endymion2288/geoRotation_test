#!/usr/bin/env python3
"""Locate charge crossings using actual replayed Calypso matrices and SDOs."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path

import numpy as np

from charge_lineage_analysis import classify_path
from slope_tables import write_json


def read_audit(path):
    with gzip.open(path, 'rt') as source:
        return [json.loads(line) for line in source]


def hit_category(barcodes, primary):
    if not barcodes:
        return 'no_SDO'
    if primary not in barcodes:
        return 'no_primary_support'
    return 'primary_only' if len(barcodes) == 1 else 'primary_and_other'


def operation(step, barcodes, primary, charge):
    p, f, c = [np.asarray(step[k], dtype=float) for k in
               ['predicted', 'filtered', 'predicted_cov']]
    h, v = [np.asarray(step[k], dtype=float) for k in ['projector', 'measurement_cov']]
    residual = np.asarray(step['measurement'])-h@p
    innovation = h@c@h.T+v
    gain = np.linalg.solve(innovation, h@c).T
    delta = f-p
    sigma = np.sqrt(np.maximum(np.diag(c), 0.))
    correlations = np.divide(c[4], sigma[4]*sigma,
                             out=np.full(6, np.nan), where=sigma[4]*sigma > 0)
    error = float(delta[4]-(gain@residual)[4])
    assert abs(error) < 1e-8
    return dict(hit_id=int(step['hit']), geometry_id=int(step['geometry_id']),
        station=(int(step['hit']) >> 59) & 3, layer=(int(step['hit']) >> 57) & 3,
        z_F_mm=step['z_F_mm'], state=step['state'],
        truth_barcodes=barcodes, hit_truth=hit_category(barcodes, primary),
        predicted_qop_GeVinv=float(p[4]), filtered_qop_GeVinv=float(f[4]),
        sigma_predicted_qop_GeVinv=float(sigma[4]),
        predicted_signed_significance=float(charge*p[4]/sigma[4]) if sigma[4] else None,
        filtered_sigma_qop_GeVinv=float(np.sqrt(max(0., step['filtered_cov'][4][4]))),
        residual_mm=residual.tolist(), innovation_mm2=innovation.tolist(),
        innovation_chi2=float(residual@np.linalg.solve(innovation, residual)),
        gain_qop_GeVinv_per_mm=gain[4].tolist(), delta_qop_GeVinv=float(delta[4]),
        identity_error_GeVinv=error, correlation_qop_loc_phi_theta_time=correlations.tolist(),
        correct_to_wrong=bool(charge*p[4] > 0 and charge*f[4] < 0),
        wrong_to_correct=bool(charge*p[4] < 0 and charge*f[4] > 0))


def inspect_target(task, target, original):
    output = Path(task['output'])
    assert (output/'exit_code.txt').read_text().strip() == '0'
    validation = json.loads((output/'replay_validation.json').read_text())
    assert validation['all_reconstruction_arrays_identical'] and validation['fields_exact'] == 58
    audit = read_audit(output/'charge_audit.jsonl.gz')
    candidate, seed = target['candidate'], target['seed']
    selected = [r for r in audit if r['stage'] == 'selected_state_trace' and
                r['candidate'] == candidate and r['seed'] == seed]
    assert len(selected) == 1
    seed_rows = [r for r in audit if r['stage'] == 'seed' and r['seed'] == seed]
    assert len(seed_rows) == 1
    seed_row = seed_rows[0]
    assert seed_row['bound'][4] == target['seed_qop_GeVinv']
    charge = 1 if task['dsid'] == '100116' else -1
    truths = [r for r in audit if r['stage'] == 'truth' and r['pdg'] == -13*charge and
              r['status'] == 1 and 0 < r['barcode'] < 200000]
    assert len(truths) == 1
    primary = truths[0]['barcode']
    mask = (original['stage'] == 8) & (original['source_entry'] == task['skip'])
    hit_truth, own = {}, {}
    for i in np.flatnonzero(mask):
        hit = int(original['hit_id'][i])
        barcodes = [int(b) for b in original['truth_barcodes'][i].split('|') if b]
        if hit in hit_truth:
            assert hit_truth[hit] == barcodes
        hit_truth[hit] = barcodes
        if int(original['candidate'][i]) == candidate:
            assert hit not in own
            own[hit] = float(original['qop_GeVinv'][i])
    updates, stages, qops = [], [1], [seed_row['bound'][4]]
    for step in selected[0]['steps']:
        if not step['measurement_flag']:
            continue
        hit = int(step['hit'])
        assert hit in own and own[hit] == step['filtered'][4]
        updates.append(operation(step, hit_truth[hit], primary, charge))
        stages.extend([2, 3])
        qops.extend([step['predicted'][4], step['filtered'][4]])
    assert len(updates) == len(own)
    stages.extend([4, 5, 6])
    qops.extend([target[k] for k in ['reference_qop_GeVinv', 'converted_qop_GeVinv', 'final_qop_GeVinv']])
    first, last = classify_path(stages, qops, charge)
    assert (first, last) == (target['first_wrong'], target['final_wrong_entry'])
    bad = [r for r in updates if r['correct_to_wrong']]
    first_update = bad[0] if bad else None
    last_update = bad[-1] if bad else None
    seed_categories = Counter(hit_category(hit_truth.get(int(hit), []), primary)
                              for hit in seed_row['hits'])
    circles = [r for r in audit if r['stage'] == 'seed_circle' and r['seed'] == seed]
    return dict(dsid=task['dsid'], sample_entry=task['sample_entry'], source_entry=task['skip'],
        file_index=task['file_index'], candidate=candidate, seed=seed, primary_barcode=primary,
        wrong=target['wrong'], final_purity=target['purity'], final_truth_p_GeV=target['truth_p_GeV'],
        seed_qop_GeVinv=target['seed_qop_GeVinv'], first_wrong=first, last_wrong_entry=last,
        final_qop_GeVinv=target['final_qop_GeVinv'],
        final_variance_GeVinv2=target['final_variance_GeVinv2'],
        first_update_into_wrong=first_update, last_update_into_wrong=last_update,
        seed_hit_categories=dict(seed_categories), seed_circle=circles,
        per_update=updates, audit=str(output/'charge_audit.jsonl.gz'),
        exact_original_reconstruction_fields=validation['fields_exact'])


def pruning(root, targets):
    result = []
    for dsid in ['100116', '100117']:
        path = root/'lineage_summary_v1'/(dsid+'_afterGeo_candidates.npz')
        with np.load(path) as archive:
            rows = {k: archive[k] for k in archive.files}
        for target in targets:
            if target['dsid'] != dsid or not target['wrong']:
                continue
            event = rows['sample_entry'] == target['sample_entry']
            correct = rows['qop_GeVinv']*(1 if dsid == '100116' else -1) > 0
            mask = event & correct & rows['primary_supported'] & (rows['rejection_reason'] > 0)
            candidates = [{k: rows[k][i].item() for k in rows} for i in np.flatnonzero(mask)]
            if candidates:
                result.append(dict(dsid=dsid, sample_entry=target['sample_entry'],
                    final_candidate=target['candidate'], primary_supported_correct_pruned=candidates))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    replay = args.root/'charge_replay_v1'
    tasks = json.loads((replay/'tasks.json').read_text())
    targets = json.loads((replay/'targets.json').read_text())
    target_by_key = {(r['dsid'], r['sample_entry']): r for r in targets}
    grouped = defaultdict(list)
    for job in tasks:
        task = json.loads((Path(job['output'])/'source_task.json').read_text())
        grouped[task['original_output']].append(task)
    records = []
    for original, tasks in grouped.items():
        with np.load(Path(original)/'charge_lineage.npz') as archive:
            fields = {key: archive[key] for key in
                      ['stage', 'source_entry', 'candidate', 'hit_id', 'qop_GeVinv', 'truth_barcodes']}
        for task in tasks:
            target = target_by_key[task['dsid'], task['sample_entry']]
            records.append(inspect_target(task, target, fields))
        print('Closed original shard', original, flush=True)
    assert len(records) == len(targets) == 119
    groups = {}
    for dsid in ['100116', '100117']:
        for cohort in ['wrong', 'normal']:
            rows = [r for r in records if r['dsid'] == dsid and r['wrong'] == (cohort == 'wrong')]
            updates = [u for row in rows for u in row['per_update']]
            last = [r['last_update_into_wrong'] for r in rows if r['last_wrong_entry'] == 3]
            groups[dsid+'|'+cohort] = dict(events=len(rows),
                final_purity_one=sum(r['final_purity'] == 1 for r in rows),
                seed_wrong=sum(r['seed_qop_GeVinv']*(1 if dsid == '100116' else -1) < 0 for r in rows),
                last_entry=dict(Counter(r['last_wrong_entry'] for r in rows)),
                last_crossing_hit_truth=dict(Counter(r['hit_truth'] for r in last)),
                last_crossing_station=dict(Counter(r['station'] for r in last)),
                last_crossing_innovation_chi2=[r['innovation_chi2'] for r in last],
                update_identity_max_error_GeVinv=max(abs(u['identity_error_GeVinv']) for u in updates),
                all_update_hit_truth=dict(Counter(u['hit_truth'] for u in updates)))
    pruned = pruning(args.root, targets)
    write_json(args.output/'events.json', records)
    write_json(args.output/'pruning.json', pruned)
    write_json(args.output/'summary.json', dict(groups=groups, events=len(records),
        all_58_arrays_exact=True, original_per_hit_filtered_qop_exact=True,
        all_stage_crossings_match_original_lineage=True, pruning_events=len(pruned),
        pruned_candidate_ties=sum(c['primary_majority_tied'] for r in pruned
            for c in r['primary_supported_correct_pruned']),
        scope='All 59 afterGeo wrong signs and 60 file-stratified modulus-997 normal controls',
        inference='SDO support describes the crossing hit; a primary-supported hit alone does not prove a correct trajectory hypothesis. Causality requires controlled hit/seed replay.'))


if __name__ == '__main__':
    main()
