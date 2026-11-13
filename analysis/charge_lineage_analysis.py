#!/usr/bin/env python3
"""Reduce actual Calypso seed/candidate/update/conversion/refit charge traces."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from charge_tables import derive
from slope_tables import write_json


ORIGINS = ['no_wrong_sign', 'seed', 'CKF_prediction', 'CKF_measurement_update',
           'native_CKF_reference_fit_or_target_propagation', 'conversion', 'refit',
           'incomplete_path']


def classify_path(stages, qops, charge):
    """Native reference drives conversion; isolated state smoothing is separate."""
    signs = np.sign(qops)*charge
    valid = np.isfinite(signs) & (signs != 0)
    wrong = np.flatnonzero(valid & (signs < 0))
    first = int(stages[wrong[0]]) if len(wrong) else 0
    last_entry = 0
    previous = 1
    for stage, sign, ok in zip(stages, signs, valid):
        if not ok:
            continue
        if sign < 0 and previous > 0:
            last_entry = int(stage)
        previous = sign
    return first, last_entry if previous < 0 else 0


def read_local(output, source):
    tables = {}
    for name in ['events', 'truth', 'tracks']:
        with np.load(output/(name+'.npz')) as archive:
            tables[name] = {k: archive[k] for k in archive.files}
        n = len(tables[name]['source_entry'])
        tables[name]['file_index'] = np.full(n, source['file_index'], dtype=np.int32)
        tables[name]['sample_entry'] = tables[name]['source_entry']-source['skip']+source['sample_entry']
    return tables


def reduce(output, destination):
    assert (output/'exit_code.txt').read_text().strip() == '0'
    source = json.loads((output/'source_task.json').read_text())
    tables = read_local(output, source)
    d = derive(tables, source['dsid'])
    h, t = tables['truth'], tables['tracks']
    primary = {int(e): int(b) for e, b in zip(h['source_entry'][d['primary']], h['barcode'][d['primary']])}
    with np.load(output/'charge_lineage.npz') as archive:
        fields = {k: archive[k] for k in ['stage', 'source_entry', 'seed', 'candidate',
            'selected', 'ordinal', 'state_order', 'geometry_id', 'hit_id', 'qop_GeVinv',
            'c44_GeVinv2', 'chi2', 'n_measurements', 'status', 'truth_barcodes']}
    max_barcode_length = int(np.char.str_len(fields['truth_barcodes']).max(initial=0))
    assert max_barcode_length < 256, 'Fixed-width barcode strings may have been truncated'
    stage = fields['stage']
    seed_rows, reference_rows, measurement_rows, selected_rows = {}, {}, defaultdict(list), defaultdict(list)
    selected_refs, converted, rejected = {}, {}, {}
    final_rows = np.flatnonzero(stage == 7)
    for row in np.flatnonzero(np.isin(stage, [0, 1, 2, 3, 4, 5, 6, 8, 10])):
        event = int(fields['source_entry'][row])
        candidate = int(fields['candidate'][row])
        key = (event, candidate)
        s = int(stage[row])
        if s == 0:
            seed_rows[event, int(fields['seed'][row])] = row
        elif s == 1:
            assert key not in reference_rows
            reference_rows[key] = row
        elif s == 8:
            measurement_rows[key].append(row)
        elif s in [2, 3, 4]:
            selected_rows[key].append(row)
        elif s == 5:
            selected_refs[key] = row
        elif s == 6:
            converted[key] = row
        elif s == 10:
            rejected[key] = row
    selected_keys = {(int(fields['source_entry'][r]), int(fields['candidate'][r])) for r in final_rows}
    barcode_cache = {s: tuple(int(b) for b in s.split('|')) if s else ()
                     for s in np.unique(fields['truth_barcodes'][stage == 8])}
    candidates = []
    by_event = defaultdict(list)
    for key, row in reference_rows.items():
        counts = Counter()
        states = measurement_rows[key]
        for state in states:
            counts.update(barcode_cache[fields['truth_barcodes'][state]])
        n = len(states)
        barcode = primary.get(key[0], 0)
        matched = counts[barcode] if barcode else 0
        maximum = max(counts.values(), default=0)
        tied = sum(v == maximum for v in counts.values()) if maximum else 0
        pure = n > 0 and matched/n > .7 and matched == maximum
        rejection = rejected.get(key)
        record = dict(source_entry=key[0], sample_entry=key[0]-source['skip']+source['sample_entry'],
            candidate=key[1], seed=int(fields['seed'][row]), barcode=barcode,
            n_measurements=n, primary_fraction=matched/n if n else 0.,
            primary_majority_tied=tied > 1 and matched == maximum,
            primary_supported=pure, qop_GeVinv=float(fields['qop_GeVinv'][row]),
            selected=key in selected_keys,
            rejection_reason=int(fields['status'][rejection]) if rejection is not None else 0,
            pruning_winner=int(fields['selected'][rejection]) if rejection is not None else -1)
        candidates.append(record)
        by_event[key[0]].append(record)
    del fields['truth_barcodes']
    track_lookup = {(int(e), int(o)): i for i, (e, o) in enumerate(zip(t['source_entry'], t['ordinal']))}
    best = {int(row): i for i, row in enumerate(d['chosen']) if row >= 0}
    finals, examples = [], []
    conversion_max = 0.
    counts = Counter()
    for row in final_rows:
        event, candidate, ordinal = (int(fields[name][row]) for name in ['source_entry', 'candidate', 'ordinal'])
        key = (event, candidate)
        seed = int(fields['seed'][row])
        assert (event, seed) in seed_rows
        sr = seed_rows[event, seed]
        rr = reference_rows[key]
        cr = converted[key]
        own = selected_rows[key]
        forward = sorted([i for i in own if stage[i] in [2, 3]],
                         key=lambda i: (-fields['state_order'][i], stage[i]))
        sequence = [sr, *forward, rr, cr, row]
        labels = [1, *[2 if stage[i] == 2 else 3 for i in forward], 4, 5, 6]
        first, final_entry = classify_path(labels, fields['qop_GeVinv'][sequence], d['charge'])
        native_ref, copy_ref = fields['qop_GeVinv'][rr], fields['qop_GeVinv'][selected_refs[key]]
        np.testing.assert_equal(native_ref, copy_ref, err_msg='Copy/isolation changed reference q/p')
        conv_qop = fields['qop_GeVinv'][cr]
        if np.isfinite(native_ref) and np.isfinite(conv_qop):
            error = abs(conv_qop-native_ref)/max(abs(native_ref), 1e-30)
            conversion_max = max(conversion_max, float(error))
            assert np.sign(conv_qop) == np.sign(native_ref)
            assert error < 1e-12
        observer = track_lookup.get((event, ordinal))
        if int(fields['status'][row]) >= 0:
            assert observer is not None
            np.testing.assert_equal(fields['qop_GeVinv'][row], t['qop_GeVinv'][observer])
        truth_row = best.get(observer)
        selected = truth_row is not None and bool(d['selected'][truth_row])
        wrong = selected and bool(d['wrong'][truth_row])
        supported = [c for c in by_event[event] if c['primary_supported']]
        correct = [c for c in supported if c['qop_GeVinv']*d['charge'] > 0]
        pure_correct_pruned = [c for c in correct if c['rejection_reason'] > 0]
        filtered = [i for i in forward if stage[i] == 3]
        smoothed = [i for i in own if stage[i] == 4]
        record = dict(source_entry=event, sample_entry=event-source['skip']+source['sample_entry'],
            ordinal=ordinal, candidate=candidate, seed=seed, selected=selected, wrong=wrong,
            matched_primary_best=truth_row is not None,
            truth_status=int(d['status'][truth_row]) if truth_row is not None else -1,
            purity=float(t['purity'][observer]) if observer is not None else 0.,
            first_wrong=first, final_wrong_entry=final_entry,
            seed_qop_GeVinv=float(fields['qop_GeVinv'][sr]), reference_qop_GeVinv=float(native_ref),
            converted_qop_GeVinv=float(conv_qop), final_qop_GeVinv=float(fields['qop_GeVinv'][row]),
            final_variance_GeVinv2=float(fields['c44_GeVinv2'][row]),
            refit_status=int(fields['status'][row]),
            n_filtered_wrong=int(np.count_nonzero(fields['qop_GeVinv'][filtered]*d['charge'] < 0)),
            n_isolated_smoothed_wrong=int(np.count_nonzero(fields['qop_GeVinv'][smoothed]*d['charge'] < 0)),
            pure_correct_candidate_exists=bool(correct),
            pure_correct_candidate_pruned=bool(pure_correct_pruned),
            primary_supported_candidates=len(supported))
        finals.append(record)
        counts['native_final'] += 1
        counts['invalid_final_parameters'] += fields['status'][row] < 0
        if selected:
            counts['selected'] += 1
            counts['selected_refit_success' if fields['status'][row] == 1 else 'selected_fallback'] += 1
            counts['selected_seed_wrong'] += fields['qop_GeVinv'][sr]*d['charge'] < 0
            counts['selected_reference_wrong'] += native_ref*d['charge'] < 0
            counts['selected_first_'+ORIGINS[first]] += 1
            if wrong:
                counts['wrong'] += 1
                counts['wrong_final_entry_'+ORIGINS[final_entry]] += 1
                counts['wrong_first_'+ORIGINS[first]] += 1
                counts['wrong_with_pure_correct_candidate'] += bool(correct)
                counts['wrong_with_pure_correct_pruned_candidate'] += bool(pure_correct_pruned)
        if wrong or (selected and record['sample_entry'] % 997 == 0):
            first_bad = next((i for i in forward if fields['qop_GeVinv'][i]*d['charge'] < 0), None)
            details = dict(record, input=source['input'], file_index=source['file_index'],
                run=int(t['run'][observer]), event=int(t['event'][observer]),
                truth_p_GeV=float(d['truth_p'][truth_row]),
                native_reference_variance_GeVinv2=float(fields['c44_GeVinv2'][rr]),
                seed_variance_GeVinv2=float(fields['c44_GeVinv2'][sr]),
                first_wrong_forward_state=None if first_bad is None else dict(
                    stage=int(stage[first_bad]), state_order=int(fields['state_order'][first_bad]),
                    geometry_id=int(fields['geometry_id'][first_bad]), hit_id=int(fields['hit_id'][first_bad]),
                    qop_GeVinv=float(fields['qop_GeVinv'][first_bad]),
                    variance_GeVinv2=float(fields['c44_GeVinv2'][first_bad])))
            examples.append(details)
    assert counts['selected'] == int(d['selected'].sum())
    assert counts['wrong'] == int(d['wrong'].sum())
    destination.mkdir(parents=True, exist_ok=False)
    for name, records in [('candidates', candidates), ('finals', finals)]:
        if records:
            np.savez_compressed(destination/(name+'.npz'), **{k: np.array([r[k] for r in records]) for k in records[0]})
    write_json(destination/'summary.json', dict(counts={k: int(v) for k, v in counts.items()},
        source_events=len(tables['events']['source_entry']), input=source,
        native_candidates=len(candidates), barcode_string_max_length=max_barcode_length,
        candidate_primary_majority_ties=sum(c['primary_majority_tied'] for c in candidates),
        conversion_qop_max_relative_error=conversion_max, origin_codes=ORIGINS,
        interpretation='P/F in forward measurement order; native target fit then conversion and refit. Isolated S is a separate candidate-state diagnostic, not the native reference mean.',
        measurement_order='state_order descending for forward filtering; ascending for smoothing',
        all_final_qop_and_identity_closed=True))
    # These rows are chosen by a frozen rule, including all selected wrong signs.
    write_json(destination/'replay_events.json', examples)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    reduce(args.source, args.output)


if __name__ == '__main__':
    main()
