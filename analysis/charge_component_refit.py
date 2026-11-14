#!/usr/bin/env python3
"""Compare finite physical-prior components along the actual Calypso chain."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from validate_charge_components import input_rows, validate


SAMPLES = ('100044', '100116', '100117')
SERIES = ('matched_n32_f256_v1', 'matched_n64_f1024_v1')


def audit(path):
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            yield json.loads(line)


def components(folder):
    records, truths = {}, {}
    event = -1
    for row in audit(folder/'prior_audit.jsonl.gz'):
        stage = row['stage']
        if stage == 'event':
            event += 1
            assert row.get('local_event_index', event) == event
        elif stage == 'truth' and row.get('signal') and row['status'] == 1 and abs(row['pdg']) == 13:
            assert event not in truths
            truths[event] = int(row['pdg'])
        if row.get('input_component', -1) < 0:
            continue
        key = event, int(row['input_component'])
        item = records.setdefault(key, dict(seed=None, selected={}, states=defaultdict(list),
                                             converted={}, final={}))
        if stage == 'seed':
            assert item['seed'] is None
            item['seed'] = row
        elif stage == 'selected_state':
            item['states'][row['candidate']].append(row)
        elif stage in ('selected', 'converted', 'final'):
            assert row['candidate'] not in item[stage]
            item[stage][row['candidate']] = row
    assert event == 19
    for (entry, _), item in records.items():
        assert entry in truths and item['seed'] is not None
        item['truth_charge'] = 1 if truths[entry] < 0 else -1
    return records


def first(item):
    ordered = sorted(item['selected'].values(),
                     key=lambda row: (-row['measurements'], row['chi2'], row['candidate']))
    return ordered[0]['candidate'] if ordered else None


def qop(row, stage):
    return float(row['bound'][4] if stage == 'selected' else row['bound_trk'][4]*1000.)


def sign(value):
    return (1 if value > 0 else -1) if value is not None and math.isfinite(value) and value else 0


def metrics(values):
    if not values:
        return dict(n=0, max=None, median=None, p95=None)
    x = np.abs(values)
    return dict(n=len(x), max=float(np.max(x)), median=float(np.median(x)),
                p95=float(np.quantile(x, .95)))


def gates(reference, new):
    """The extra conversion/refit must not change any previously recorded CKF value."""
    assert reference.keys() == new.keys()
    counts = Counter()
    for key, a in reference.items():
        b = new[key]
        for stage in ('seed', 'selected', 'states'):
            assert a[stage] == b[stage], (key, stage)
        counts['components'] += 1
        counts['selected_candidates'] += len(a['selected'])
        counts['state_records'] += sum(map(len, a['states'].values()))
        assert b['converted'].keys() == b['final'].keys() == b['selected'].keys()
        for candidate, row in b['selected'].items():
            delta = qop(b['converted'][candidate], 'converted') - qop(row, 'selected')
            assert abs(delta) < 5e-13 * max(1., abs(qop(row, 'selected'))), (key, delta)
            counts['conversion_sign_changes'] += sign(qop(b['converted'][candidate], 'converted')) != sign(qop(row, 'selected'))
    return dict(counts, CKF_seed_selected_and_state_records_exact=True)


def compare(a, b):
    assert a.keys() == b.keys()
    counts = Counter(components=len(a))
    stage_differences = defaultdict(list)
    normalized_differences = defaultdict(list)
    examples, rows = [], []
    for key in sorted(a):
        before, after = a[key], b[key]
        ia, ib = first(before), first(after)
        counts['selected_in_both'] += ia is not None and ib is not None
        counts['before_only_selected'] += ia is not None and ib is None
        counts['after_only_selected'] += ia is None and ib is not None
        record = dict(entry=key[0], component=key[1], truth_charge=before['truth_charge'],
                      native_seed=before['seed']['seed'],
                      before_selected=ia is not None, after_selected=ib is not None)
        for label, item, candidate in [('before', before, ia), ('after', after, ib)]:
            record[label+'_refit_success'] = (item['final'][candidate]['refit'] == 1
                if candidate in item['final'] else None)
            record[label+'_final_hits'] = (item['final'][candidate]['hits']
                if candidate in item['final'] else None)
            for stage in ('selected', 'converted', 'final'):
                row = item[stage].get(candidate)
                value = qop(row, stage) if row else None
                record[label+'_'+stage+'_qop_GeVinv'] = value
                record[label+'_'+stage+'_qop_variance_GeVinv2'] = row.get('variance_qop') if row else None
                record[label+'_'+stage+'_correct'] = (sign(value) == item['truth_charge']) if sign(value) else None
            if candidate is not None:
                counts[label+'_selected'] += 1
                for stage in ('selected', 'converted', 'final'):
                    if record[label+'_'+stage+'_correct'] is not None:
                        counts[label+'_'+stage+('_correct' if record[label+'_'+stage+'_correct'] else '_wrong')] += 1
                if candidate in item['final']:
                    counts[label+'_refit_success'] += item['final'][candidate]['refit'] == 1
                    counts[label+'_refit_fallback'] += item['final'][candidate]['refit'] == 0
                values = [record[label+'_'+stage+'_qop_GeVinv'] for stage in ('selected', 'converted', 'final')]
                if all(sign(v) for v in values):
                    counts[label+'_CKF_to_final_sign_changes'] += sign(values[0]) != sign(values[2])
        if ia is not None and ib is not None:
            same_path = before['selected'][ia]['hits'] == after['selected'][ib]['hits']
            counts['same_selected_hit_sequence'] += same_path
            record['same_selected_hit_sequence'] = same_path
            for stage in ('selected', 'converted', 'final'):
                va, vb = (record[k+'_'+stage+'_qop_GeVinv'] for k in ('before', 'after'))
                if va is not None and vb is not None:
                    stage_differences[stage].append(vb-va)
                    ca, cb = (record[k+'_'+stage+'_qop_variance_GeVinv2'] for k in ('before', 'after'))
                    if ca is not None and cb is not None and math.isfinite(ca+cb) and ca+cb > 0:
                        normalized_differences[stage].append((vb-va)/math.sqrt((ca+cb)/2))
                    counts[stage+'_sign_migrations'] += sign(va) != sign(vb)
                    if stage == 'final' and sign(va) != sign(vb):
                        examples.append(record.copy())
            if record['before_refit_success'] is not None and record['after_refit_success'] is not None:
                for left, right, category in [(True, True, 'refit_success_both'),
                                              (True, False, 'refit_before_only'),
                                              (False, True, 'refit_after_only'),
                                              (False, False, 'refit_fallback_both')]:
                    counts[category] += record['before_refit_success'] is left and record['after_refit_success'] is right
                counts['same_final_hit_sequence'] += record['before_final_hits'] == record['after_final_hits']
                if record['before_refit_success'] and record['after_refit_success']:
                    stage_differences['final_refit_success_both'].append(
                        record['after_final_qop_GeVinv']-record['before_final_qop_GeVinv'])
                    counts['final_refit_success_both_sign_migrations'] += (
                        sign(record['before_final_qop_GeVinv']) != sign(record['after_final_qop_GeVinv']))
            sa, sb = before['states'][ia], after['states'][ib]
            state_keys = lambda seq: [(r.get('hit'), r['measurement_flag'], r['outlier_flag'], r['hole_flag']) for r in seq]
            if state_keys(sa) == state_keys(sb):
                counts['same_full_state_path'] += 1
                for ra, rb in zip(sa, sb):
                    if not ra['measurement_flag']:
                        continue
                    counts['paired_measurement_states'] += 1
                    for state_stage in ('predicted', 'filtered', 'smoothed'):
                        va, vb = ra[state_stage][4], rb[state_stage][4]
                        stage_differences[state_stage].append(vb-va)
                        counts[state_stage+'_state_sign_migrations'] += sign(va) != sign(vb)
        rows.append(record)
    return dict(counts=dict(counts), absolute_qop_differences_GeVinv={k: metrics(v) for k, v in stage_differences.items()},
                absolute_qop_difference_over_reported_sigma={k: metrics(v) for k, v in normalized_differences.items()},
                final_sign_migration_examples=examples), rows


def width_input_checks(root, series_list, samples):
    """Check the same standardized cloud while preserving the total prior C."""
    results = []
    manifests = {name: json.loads((root/'common_prior_components_v4'/name/'manifest.json').read_text())
                 for name in series_list}
    for i, first_series in enumerate(series_list):
        first_manifest = manifests[first_series]
        for second_series in series_list[i+1:]:
            second_manifest = manifests[second_series]
            if first_manifest['n_components'] != second_manifest['n_components']:
                continue
            f1 = first_manifest['within_component_covariance_fraction']
            f2 = second_manifest['within_component_covariance_fraction']
            scale = np.sqrt((1-f2)/(1-f1))
            counts = 0
            mean_error, covariance_error = 0., 0.
            for dsid in samples:
                a = input_rows(root/'common_prior_components_v4'/first_series/dsid/'beforeGeo_components.txt')
                b = input_rows(root/'common_prior_components_v4'/second_series/dsid/'beforeGeo_components.txt')
                assert a.keys() == b.keys()
                seeds = defaultdict(list)
                for key, row in a.items():
                    seeds[(key[0], row['seed'])].append(key)
                for keys in seeds.values():
                    centers = np.asarray([a[key]['bound'] for key in keys])
                    mean = centers.mean(axis=0)
                    for key in keys:
                        expected = mean+scale*(a[key]['bound']-mean)
                        np.testing.assert_allclose(b[key]['bound'], expected, rtol=2e-13, atol=2e-12)
                        np.testing.assert_allclose(b[key]['cov'], a[key]['cov']*f2/f1, rtol=2e-13, atol=1e-20)
                        mean_error = max(mean_error, float(np.max(np.abs(b[key]['bound']-expected))))
                        covariance_error = max(covariance_error, float(np.max(np.abs(b[key]['cov']-a[key]['cov']*f2/f1))))
                        counts += 1
            results.append(dict(first=first_series, second=second_series, paired_components=counts,
                same_standardized_cloud=True, center_rescaling_for_fixed_total_C=float(scale),
                mean_relation_max_absolute_error=mean_error, covariance_relation_max_absolute_error=covariance_error,
                scope='The same standardized cloud and total mean/C; finite mixture distributions are approximations, not identical posteriors.'))
    return results


def sampling_checks(results, rows, samples):
    comparisons = []
    series_list = list(results)
    for i, first_series in enumerate(series_list):
        a = results[first_series]
        for second_series in series_list[i+1:]:
            b = results[second_series]
            same_count = a['n_components'] == b['n_components']
            same_fraction = a['within_fraction'] == b['within_fraction']
            if not (same_count or same_fraction):
                continue
            for dsid in samples:
                first_rows = [r for r in rows if r['series'] == first_series and r['dsid'] == dsid and r['comparison'] == 'CKF']
                second_rows = [r for r in rows if r['series'] == second_series and r['dsid'] == dsid and r['comparison'] == 'CKF']
                item = dict(first=first_series, second=second_series, dsid=dsid,
                            same_component_count=same_count, same_within_fraction=same_fraction,
                            kind='within_covariance_split' if same_count else 'component_count')
                for name, current in [('first', first_rows), ('second', second_rows)]:
                    selected = sum(row['after_selected'] for row in current)
                    wrong = sum(row['after_selected_correct'] is False for row in current)
                    item[name] = dict(components=len(current), selected=selected, wrong=wrong,
                        selected_fraction=selected/len(current), wrong_among_selected=wrong/selected if selected else None)
                if same_count:
                    first = {(r['entry'], r['component']): r for r in first_rows}
                    second = {(r['entry'], r['component']): r for r in second_rows}
                    assert first.keys() == second.keys()
                    counts = Counter()
                    differences = []
                    for key, left in first.items():
                        right = second[key]
                        la, lb = left['after_selected'], right['after_selected']
                        counts['selected_both'] += la and lb
                        counts['first_only_selected'] += la and not lb
                        counts['second_only_selected'] += lb and not la
                        if la and lb:
                            va, vb = left['after_selected_qop_GeVinv'], right['after_selected_qop_GeVinv']
                            differences.append(vb-va)
                            counts['selected_qop_sign_changes'] += sign(va) != sign(vb)
                    item['same_cloud_paired_outcomes'] = dict(counts)
                    item['same_cloud_absolute_qop_change_GeVinv'] = metrics(differences)
                comparisons.append(item)
    return comparisons


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--full-chain', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--series', nargs='+', default=list(SERIES))
    parser.add_argument('--samples', nargs='+', default=list(SAMPLES))
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    results, validations, rows, ckf_gates = {}, [], [], []
    for series in args.series:
        prior = args.root/'common_prior_components_v4'/series
        manifest = json.loads((prior/'manifest.json').read_text())
        info = dict(n_components=manifest['n_components'], within_fraction=manifest['within_component_covariance_fraction'],
                    groups={}, input_covariance_scaled_error_max=0., input_direction_error_max=0.)
        for dsid in args.samples:
            checks = manifest['samples'][dsid]['input_checks']
            info['input_covariance_scaled_error_max'] = max(info['input_covariance_scaled_error_max'],
                max(c['total_F_covariance_scaled_error'] for c in checks))
            info['input_direction_error_max'] = max(info['input_direction_error_max'],
                max(c['maximum_inner_sigma_direction_error'] for c in checks))
            for variant in ('beforeGeo', 'afterGeo'):
                validations.append(validate(prior/dsid/variant))
            a = components(prior/dsid/'beforeGeo')
            b = components(prior/dsid/'afterGeo')
            comparison, paired = compare(a, b)
            info['groups'][dsid] = dict(CKF=comparison)
            for row in paired:
                rows.append(dict(series=series, dsid=dsid, comparison='CKF', **row))
            chains = {}
            for variant, modes in [('beforeGeo', ('native', 'regular_both')), ('afterGeo', ('native',))]:
                for mode in modes:
                    folder = args.full_chain/series/dsid/(variant+'_'+mode)
                    validations.append(validate(folder))
                    chains[variant+'_'+mode] = components(folder)
                    gate = gates(a if variant == 'beforeGeo' else b, chains[variant+'_'+mode])
                    ckf_gates.append(dict(series=series, dsid=dsid, variant=variant, mode=mode, **gate))
            for name, first_key, second_key in [
                    ('full_native', 'beforeGeo_native', 'afterGeo_native'),
                    ('before_refit_chart', 'beforeGeo_native', 'beforeGeo_regular_both'),
                    ('full_regular_before', 'beforeGeo_regular_both', 'afterGeo_native')]:
                comparison, paired = compare(chains[first_key], chains[second_key])
                info['groups'][dsid][name] = comparison
                for row in paired:
                    rows.append(dict(series=series, dsid=dsid, comparison=name, **row))
        results[series] = info
    width_checks = width_input_checks(args.root, args.series, args.samples)
    (output/'summary.json').write_text(json.dumps(dict(series=results, CKF_gates=ckf_gates,
        same_cloud_width_input_checks=width_checks,
        sampling_comparisons=sampling_checks(results, rows, args.samples),
        scope='Common moment-matched finite prior components in actual Calypso CKF/conversion/refit; '
              '20 source events per sample. Component counts are diagnostic ensembles, not production efficiencies. '
              'First-order finite component C mapping and no posterior likelihood aggregation.'), indent=2)+'\n')
    (output/'validation.json').write_text(json.dumps(dict(all_complete=True, checks=validations,
        CKF_and_conversion_gates_passed=True, reference_CKF_gates=ckf_gates), indent=2)+'\n')
    (output/'paired_components.json').write_text(json.dumps(rows, indent=2)+'\n')
    print(json.dumps(dict(output=str(output), validated_jobs=len(validations),
                         full_chain_jobs=len(ckf_gates), paired_rows=len(rows)), indent=2))


if __name__ == '__main__':
    main()
