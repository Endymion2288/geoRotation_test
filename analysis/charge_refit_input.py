#!/usr/bin/env python3
"""Audit physical refit priors and native state differences in actual Calypso."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from charge_component_refit import audit, components, first, metrics, qop, sign
from slope_tables import write_json
from validate_charge_components import validate


ROTATE_F_TO_A = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
SAMPLES = ('100044', '100116', '100117')


def physical(row, variant, state_stage=None):
    """Express a recorded bound mean/C as position, unit direction, q/p, time.

    This is an analytic coordinate Jacobian, not track transport or fitting.
    Its 8D covariance has the expected two surface/direction constraints.
    """
    if state_stage is None:
        bound, covariance = np.asarray(row['bound']), np.asarray(row['cov'])
    else:
        bound = np.asarray(row[state_stage])
        covariance = np.asarray(row[state_stage+'_cov'])
    surface = np.asarray(row['surface'])
    phi, theta = bound[2:4]
    cp, sp, ct, st = np.cos(phi), np.sin(phi), np.cos(theta), np.sin(theta)
    direction = np.array([cp*st, sp*st, ct])
    rotation = np.eye(3) if variant == 'beforeGeo' else ROTATE_F_TO_A.T
    position = surface[:3, :2] @ bound[:2]+surface[:3, 3]
    jacobian = np.zeros((8, 6))
    jacobian[:3, :2] = rotation @ surface[:3, :2]
    jacobian[3:6, 2] = rotation @ np.array([-sp*st, cp*st, 0.])
    jacobian[3:6, 3] = rotation @ np.array([cp*ct, sp*ct, -st])
    jacobian[6, 4] = jacobian[7, 5] = 1.
    mean = np.r_[rotation @ position, rotation @ direction, bound[4:]]
    physical_covariance = jacobian @ covariance @ jacobian.T
    if state_stage is None:
        suffix = 'native' if variant == 'beforeGeo' else 'A'
        np.testing.assert_allclose(position, row['position_'+suffix+'_mm'], rtol=0, atol=2e-8)
        momentum = np.asarray(row['momentum_'+suffix+'_GeV'])
        np.testing.assert_allclose(direction, momentum/np.linalg.norm(momentum), rtol=0, atol=2e-12)
    return mean, physical_covariance


def refits(folder):
    result = {}
    event = -1
    for row in audit(folder/'prior_audit.jsonl.gz'):
        if row['stage'] == 'event':
            event += 1
        if row.get('input_component', -1) < 0 or row['stage'] not in ('refit_input', 'refit', 'refit_state'):
            continue
        key = event, row['input_component'], row['candidate']
        item = result.setdefault(key, dict(states=[]))
        if row['stage'] == 'refit_state':
            item['states'].append(row)
        else:
            assert row['stage'] not in item
            item[row['stage']] = row
    assert event == 19
    return result


def observer_gate(old, new, compare_final=True):
    assert old.keys() == new.keys()
    counts = Counter()
    for key, previous in old.items():
        current = new[key]
        stages = ['seed', 'selected', 'states', 'converted', 'truth_charge']
        if compare_final:
            stages.append('final')
        for stage in stages:
            assert previous[stage] == current[stage], (key, stage)
        counts['components'] += 1
        counts['selected_candidates'] += len(current['selected'])
        counts['selected_state_records'] += sum(map(len, current['states'].values()))
    return dict(counts, existing_CKF_and_conversion_records_exact=True,
                final_records_compared_and_exact=compare_final)


def input_contract(folder, tracks, records):
    """Check the actual native copied-C/diagonal inflation contract and coverage."""
    run = json.loads((folder/'run.json').read_text())
    scale = run['refit_covariance_scale']
    count, maximum = 0, 0.
    for key, item in tracks.items():
        for candidate, final in item['final'].items():
            refit = records.get((*key, candidate), {})
            if final['refit'] == 1:
                assert 'refit_input' in refit and 'refit' in refit
            if 'refit_input' not in refit:
                continue
            row = refit['refit_input']
            converted = item['converted'][candidate]
            expected = np.eye(6)
            expected[:5, :5] = np.asarray(converted['cov_trk'])
            units = np.array([1., 1., 1., 1., 1000., 1.])
            expected *= np.outer(units, units)
            expected[np.diag_indices(6)] *= scale
            expected *= 10.
            np.testing.assert_allclose(row['cov'], expected, rtol=2e-13, atol=1e-18)
            np.testing.assert_allclose(row['bound'][4], qop(converted, 'converted'), rtol=2e-13, atol=1e-15)
            maximum = max(maximum, float(np.max(np.abs(np.asarray(row['cov'])-expected))))
            count += 1
    return dict(refit_inputs=count, diagonal_scale=scale, uniform_multiplier=10.,
        copied_Trk_covariance_unit_conversion_and_inflation_checked=True,
        maximum_covariance_absolute_error=maximum,
        every_successful_final_has_native_input_and_refit_record=True)


def scale_gate(default, control, default_refit, control_refit):
    check = observer_gate(default, control, compare_final=False)
    assert default_refit.keys() == control_refit.keys()
    checked = 0
    for key, old in default_refit.items():
        current = control_refit[key]
        assert ('refit_input' in old) == ('refit_input' in current)
        if 'refit_input' not in old:
            continue
        a, b = old['refit_input'], current['refit_input']
        for field in ('bound', 'surface'):
            assert a[field] == b[field], (key, field)
        expected = np.asarray(a['cov']).copy()
        expected[np.diag_indices(6)] /= 10.
        np.testing.assert_allclose(b['cov'], expected, rtol=2e-13, atol=1e-18)
        checked += 1
    return dict(check, refit_input_pairs=checked, means_and_surfaces_exact=True,
                covariance_only_diagonal_divided_by_10=True)


def input_difference(before, after):
    ma, ca = physical(before, 'beforeGeo')
    mb, cb = physical(after, 'afterGeo')
    diagonal = np.maximum((np.abs(np.diag(ca))+np.abs(np.diag(cb)))/2., 1e-30)
    normalized = np.abs(cb-ca)/np.sqrt(np.outer(diagonal, diagonal))
    result = dict(position_difference_mm=float(np.linalg.norm(mb[:3]-ma[:3])),
        direction_difference_norm=float(np.linalg.norm(mb[3:6]-ma[3:6])),
        qop_difference_GeVinv=float(abs(mb[6]-ma[6])),
        physical_covariance_scaled_max=float(np.max(normalized)),
        position_covariance_scaled_max=float(np.max(normalized[:3, :3])),
        direction_covariance_scaled_max=float(np.max(normalized[3:6, 3:6])),
        direction_qop_covariance_scaled_max=float(np.max(normalized[3:6, 6])),
        before_direction_variance_trace=float(np.trace(ca[3:6, 3:6])),
        after_direction_variance_trace=float(np.trace(cb[3:6, 3:6])))
    return result


def state_differences(before, after):
    def keyed(states):
        rows = [r for r in states if r.get('hit') is not None]
        result = {(r['hit'], r['measurement_flag'], r['outlier_flag'], r['hole_flag']): r for r in rows}
        assert len(result) == len(rows)
        return result
    a, b = keyed(before), keyed(after)
    rows = []
    for key in a.keys() & b.keys():
        ra, rb = a[key], b[key]
        item = dict(hit=key[0], measurement=key[1], outlier=key[2], hole=key[3],
                    path_length_mm=ra['path_length_mm'])
        for stage in ('predicted', 'filtered', 'smoothed'):
            if ra.get(stage) is None or rb.get(stage) is None:
                continue
            qa, qb = ra[stage][4], rb[stage][4]
            item[stage+'_qop_before_GeVinv'] = qa
            item[stage+'_qop_after_GeVinv'] = qb
            item[stage+'_sign_difference'] = sign(qa) != sign(qb)
        rows.append(item)
    return sorted(rows, key=lambda r: r['path_length_mm'])


def compare(before, after, before_refit, after_refit):
    assert before.keys() == after.keys()
    counts, rows, examples, seed_checks = Counter(), [], [], []
    for key, a in before.items():
        b = after[key]
        seed_checks.append(input_difference(a['seed'], b['seed']))
        ia, ib = first(a), first(b)
        if ia is None or ib is None:
            continue
        counts['selected_both'] += 1
        ra = before_refit.get((*key, ia), {})
        rb = after_refit.get((*key, ib), {})
        same_hits = a['selected'][ia]['hits'] == b['selected'][ib]['hits']
        qa, qb = qop(a['final'][ia], 'final'), qop(b['final'][ib], 'final')
        item = dict(entry=key[0], component=key[1], native_seed=a['seed']['seed'],
            same_CKF_hits=same_hits, before_refit_success=a['final'][ia]['refit'] == 1,
            after_refit_success=b['final'][ib]['refit'] == 1,
            final_sign_difference=sign(qa) != sign(qb),
            before_final_qop_GeVinv=qa, after_final_qop_GeVinv=qb,
            before_has_refit_input='refit_input' in ra, after_has_refit_input='refit_input' in rb)
        counts['refit_input_before'] += item['before_has_refit_input']
        counts['refit_input_after'] += item['after_has_refit_input']
        if item['before_has_refit_input'] and item['after_has_refit_input']:
            item.update(input_difference(ra['refit_input'], rb['refit_input']))
            counts['refit_input_both'] += 1
        if item['final_sign_difference']:
            examples.append(dict(item, paired_refit_states=state_differences(ra.get('states', []), rb.get('states', []))))
        rows.append(item)
    cohorts = {}
    for name, predicate in [
            ('all_input_pairs', lambda r: True),
            ('same_CKF_hits', lambda r: r['same_CKF_hits']),
            ('both_refit_success', lambda r: r['before_refit_success'] and r['after_refit_success']),
            ('same_hits_both_refit_success', lambda r: r['same_CKF_hits'] and r['before_refit_success'] and r['after_refit_success'])]:
        selected = [r for r in rows if 'physical_covariance_scaled_max' in r and predicate(r)]
        cohorts[name] = dict(n=len(selected), metrics={field: metrics([r[field] for r in selected])
            for field in ('position_difference_mm', 'direction_difference_norm', 'qop_difference_GeVinv',
                          'physical_covariance_scaled_max', 'position_covariance_scaled_max',
                          'direction_covariance_scaled_max', 'direction_qop_covariance_scaled_max')})
    seed_metrics = {field: metrics([r[field] for r in seed_checks]) for field in seed_checks[0]}
    return dict(counts=dict(counts), cohorts=cohorts, seed_physical_prior_metrics=seed_metrics,
                final_sign_difference_examples=examples), rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--replay', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scale-control', type=Path)
    parser.add_argument('--series', nargs='+', default=['matched_n32_f256_v1'])
    args = parser.parse_args()
    original = args.root/'common_prior_component_refit_v1'
    args.output.mkdir(parents=True, exist_ok=False)
    summary, checks, rows, controls, control_checks = {}, [], [], {}, []
    for series in args.series:
        for dsid in SAMPLES:
            after_folder = original/series/dsid/'afterGeo_native'
            after, after_refit = components(after_folder), refits(after_folder)
            if args.scale_control:
                after_control_folder = args.scale_control/series/dsid/'afterGeo_native'
                after_control = components(after_control_folder)
                after_control_refit = refits(after_control_folder)
                control_checks.append(dict(series=series, dsid=dsid, mode='afterGeo',
                    runtime=validate(after_control_folder),
                    contract=input_contract(after_control_folder, after_control, after_control_refit),
                    scale=scale_gate(after, after_control, after_refit, after_control_refit)))
            for mode in ('native', 'regular_both'):
                name = 'beforeGeo_'+mode
                folder = args.replay/series/dsid/name
                before = components(folder)
                before_refit = refits(folder)
                gate = observer_gate(components(original/series/dsid/name), before)
                checks.append(dict(series=series, dsid=dsid, mode=mode, **gate,
                    runtime=validate(folder), contract=input_contract(folder, before, before_refit)))
                result, pairs = compare(before, after, before_refit, after_refit)
                summary[series+'|'+dsid+'|'+mode] = result
                rows.extend(dict(series=series, dsid=dsid, mode=mode, **r) for r in pairs)
                if args.scale_control:
                    control_folder = args.scale_control/series/dsid/name
                    control = components(control_folder)
                    control_refit = refits(control_folder)
                    control_checks.append(dict(series=series, dsid=dsid, mode=name,
                        runtime=validate(control_folder),
                        contract=input_contract(control_folder, control, control_refit),
                        scale=scale_gate(before, control, before_refit, control_refit)))
                    result, pairs = compare(control, after_control, control_refit, after_control_refit)
                    controls[series+'|'+dsid+'|'+mode] = result
                    rows.extend(dict(series=series, dsid=dsid, mode=mode+'_scale1', **r) for r in pairs)
    write_json(args.output/'summary.json', dict(groups=summary, scale1_groups=controls,
        scope='Recorded native refit inputs in a common physical frame; conditional diagnostic components, not production rates.'))
    write_json(args.output/'paired_refit_inputs.json', rows)
    write_json(args.output/'validation.json', dict(observer_gates=checks, scale_control_gates=control_checks,
        complete_existing_records_exact=True, analytic_physical_coordinate_Jacobian=True,
        original_full_chain=str(original), replay=str(args.replay)))
    print(json.dumps(dict(groups=len(summary), pairs=len(rows), observer_gates=len(checks)), indent=2))


if __name__ == '__main__':
    main()
