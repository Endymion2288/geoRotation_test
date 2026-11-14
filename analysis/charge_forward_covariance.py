#!/usr/bin/env python3
"""Locate forward differences in recorded, complete native Calypso paths.

Coordinate changes and covariance balances below are diagnostics on native
outputs. They perform no track propagation, measurement update or smoothing.
The balance C_pred - J C_prev J^T is not labelled native process noise: material
updates and any other operations absent from the saved J can contribute.
"""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path

import numpy as np

from charge_component_refit import components, first, metrics
from bound_surface_coordinates import direction_f, surface_position, ROTATE_F_TO_A
from charge_smoothing_operations import scaled
from slope_tables import write_json


STAGES = ('predicted', 'filtered', 'smoothed')
SAMPLES = ('100044', '100116', '100117')
THRESHOLDS = (1e-8, 1e-6, 1e-4)


def chart(bound, variant):
    """Native bound -> common (local0, local1, phi_A, latitude_A, q/p, t).

    The before Jacobian is the analytic map used by RefitDirectionChart.h.
    Local coordinates are retained only after the paired surfaces are checked.
    """
    result = np.asarray(bound, dtype=float).copy()
    jacobian = np.eye(6)
    if variant == 'afterGeo':
        result[3] = np.pi / 2. - result[3]
        jacobian[3, 3] = -1.
        return result, jacobian
    phi, theta = result[2:4]
    sp, cp, st, ct = np.sin(phi), np.cos(phi), np.sin(theta), np.cos(theta)
    direction = ROTATE_F_TO_A @ np.array([cp*st, sp*st, ct])
    transverse2 = direction[0]**2 + direction[1]**2
    transverse = np.sqrt(transverse2)
    result[2:4] = np.arctan2(direction[1], direction[0]), np.arctan2(direction[2], transverse)
    derivatives = ROTATE_F_TO_A @ np.array([[-sp*st, cp*ct], [cp*st, sp*ct], [0., -st]])
    gradient = np.array([[-direction[1]/transverse2, direction[0]/transverse2, 0.],
                         [-direction[0]*direction[2]/transverse,
                          -direction[1]*direction[2]/transverse, transverse]])
    jacobian[2:4, 2:4] = gradient @ derivatives
    return result, jacobian


def common(row, variant, stage=None):
    bound = row['bound'] if stage is None else row[stage]
    covariance = row['cov'] if stage is None else row[stage+'_cov']
    mean, jacobian = chart(bound, variant)
    return mean, jacobian @ np.asarray(covariance) @ jacobian.T, jacobian


def matrix_metrics(a, b, reference_a=None, reference_b=None):
    delta = scaled(a, b,
        np.diag(reference_a) if reference_a is not None else None,
        np.diag(reference_b) if reference_b is not None else None)
    return dict(all_scaled_max=float(delta.max()),
        direction_scaled_max=float(delta[2:4, 2:4].max()),
        direction_absolute_max_rad2=float(np.abs(a[2:4, 2:4]-b[2:4, 2:4]).max()),
        direction_frobenius_rad2=float(np.linalg.norm(a[2:4, 2:4]-b[2:4, 2:4])),
        qop_variance_absolute_GeVinv2=float(abs(a[4, 4]-b[4, 4])))


def state_metrics(a, b, stage):
    ma, ca, _ = common(a, 'beforeGeo', stage)
    mb, cb, _ = common(b, 'afterGeo', stage)
    da, ga = direction_f(np.asarray(a[stage]), 'beforeGeo')
    db, gb = direction_f(np.asarray(b[stage]), 'afterGeo')
    pca = ga @ np.asarray(a[stage+'_cov']) @ ga.T
    pcb = gb @ np.asarray(b[stage+'_cov']) @ gb.T
    position_difference = None
    if a.get('surface_type') is not None and b.get('surface_type') is not None:
        assert a['surface_type'] == b['surface_type']
        pa, _ = surface_position(a[stage], a['surface'], a['surface_type'])
        pb, _ = surface_position(b[stage], b['surface'], b['surface_type'])
        position_difference = float(np.linalg.norm(pa-ROTATE_F_TO_A.T @ pb))
    angle_delta = mb[2:4]-ma[2:4]
    angle_delta[0] = np.arctan2(np.sin(angle_delta[0]), np.cos(angle_delta[0]))
    return dict(**matrix_metrics(ca, cb),
        position_difference_mm=position_difference,
        direction_difference_norm=float(np.linalg.norm(da-db)),
        angular_difference_rad=float(np.abs(angle_delta).max()),
        qop_difference_GeVinv=float(abs(ma[4]-mb[4])),
        physical_direction_C_scaled_max=float(scaled(pca, pcb).max()),
        before_direction_variance_trace_rad2=float(np.trace(ca[2:4, 2:4])),
        after_direction_variance_trace_rad2=float(np.trace(cb[2:4, 2:4])))


def covariance_balance(ca, cb, ja, jb, pa, pb):
    """Exact accounting for after-minus-before C_pred in one common chart.

    This separates the saved previous-C and J contractions, plus the remaining
    covariance balance. It is not a counterfactual native transport experiment.
    All term comparisons use the same target predicted-C normalization.
    """
    propagated_a, propagated_b = ja @ ca @ ja.T, jb @ cb @ jb.T
    balance_a, balance_b = pa-propagated_a, pb-propagated_b
    inherited = ja @ (cb-ca) @ ja.T
    jacobian = propagated_b-ja @ cb @ ja.T
    balance = balance_b-balance_a
    delta = pb-pa
    terms = dict(inherited_C=inherited, saved_J=jacobian, remaining_balance=balance)
    residual = delta-sum(terms.values())
    denominator = np.sqrt(np.outer(
        np.maximum((np.abs(np.diag(pa))+np.abs(np.diag(pb)))/2., 1e-30),
        np.maximum((np.abs(np.diag(pa))+np.abs(np.diag(pb)))/2., 1e-30)))
    error = float(np.abs(residual/denominator).max())
    scale = max(1., *(float(np.abs(v/denominator).max()) for v in terms.values()))
    assert error < 1e-8*scale, (error, scale)
    result = dict(decomposition_scaled_residual=error,
        before_balance_direction_rad2=balance_a[2:4, 2:4].tolist(),
        after_balance_direction_rad2=balance_b[2:4, 2:4].tolist(),
        target_delta_direction_rad2=delta[2:4, 2:4].tolist())
    for name, value in terms.items():
        result[name] = dict(**matrix_metrics(value, np.zeros((6, 6)), pa, pb),
                           direction_matrix_rad2=value[2:4, 2:4].tolist())
    source_std = np.sqrt(np.maximum((np.abs(np.diag(ca))+np.abs(np.diag(cb)))/2., 1e-30))
    target_std = np.sqrt(np.maximum((np.abs(np.diag(pa))+np.abs(np.diag(pb)))/2., 1e-30))
    result['saved_J_difference_scaled_max'] = float(np.abs((jb-ja)*source_std[None, :]/target_std[:, None]).max())
    result['saved_J_difference_absolute_max'] = float(np.abs(jb-ja).max())
    return result


def forward_states(item, candidate):
    states = list(reversed(item['states'][candidate]))
    assert states
    assert states[0]['acts_previous_index'] == 4294967295
    for previous, current in zip(states, states[1:]):
        assert current['acts_previous_index'] == previous['acts_state_index']
    return states


def signature(row):
    return tuple(row.get(field) for field in
        ('geometry_id', 'hit', 'measurement_flag', 'outlier_flag', 'hole_flag'))


def paired_path(a, b, ia, ib, tag, surface_catalog=None):
    before, after = forward_states(a, ia), forward_states(b, ib)
    assert list(map(signature, before)) == list(map(signature, after)), tag
    rows = []
    rotation = np.eye(4)
    rotation[:3, :3] = ROTATE_F_TO_A
    previous = (a['seed'], b['seed'])
    previous_stage = None
    for ordinal, (ra, rb) in enumerate(zip(before, after)):
        for variant, state in [('beforeGeo', ra), ('afterGeo', rb)]:
            metadata = (surface_catalog or {}).get(variant+'|'+str(state['geometry_id']))
            if metadata is not None:
                np.testing.assert_allclose(state['surface'], metadata['surface'], rtol=0, atol=2e-9)
                state['surface_type'] = metadata['surface_type']
        surface_error = float(np.abs(rotation @ np.asarray(ra['surface'])-np.asarray(rb['surface'])).max())
        assert surface_error < 2e-9, (tag, ordinal, surface_error)
        assert all(stage in ra and stage in rb for stage in STAGES)
        item = dict(**tag, forward_index=ordinal, path_length_mm=ra['path_length_mm'],
            geometry_id=ra['geometry_id'], hit=ra.get('hit'),
            surface_type=ra.get('surface_type'),
            measurement_flag=ra['measurement_flag'], outlier_flag=ra['outlier_flag'],
            hole_flag=ra['hole_flag'], surface_rotation_error_mm=surface_error,
            path_length_difference_mm=rb['path_length_mm']-ra['path_length_mm'])
        item['path_step_mm'] = (ra['path_length_mm']-before[ordinal-1]['path_length_mm']
                                if ordinal else None)
        for stage in STAGES:
            item[stage] = state_metrics(ra, rb, stage)
        item['material'] = {name: {field: r.get(field) for field in
            ('material_full_X0', 'material_full_thickness_mm', 'material_path_correction')}
            for name, r in [('before', ra), ('after', rb)]}
        item['material']['raw_slab_equal'] = all(ra.get(field) == rb.get(field)
            for field in ('material_full_X0', 'material_full_thickness_mm'))
        _, ca, ba = common(previous[0], 'beforeGeo', previous_stage)
        _, cb, bb = common(previous[1], 'afterGeo', previous_stage)
        _, pa, ta = common(ra, 'beforeGeo', 'predicted')
        _, pb, tb = common(rb, 'afterGeo', 'predicted')
        # Stored J links the preceding filtered state to this prediction;
        # the first stored state instead uses the explicit input seed.
        ja = ta @ np.asarray(ra['jacobian']) @ np.linalg.inv(ba)
        jb = tb @ np.asarray(rb['jacobian']) @ np.linalg.inv(bb)
        item['link'] = dict(from_seed=ordinal == 0,
            previous_filtered_difference=matrix_metrics(ca, cb),
            **covariance_balance(ca, cb, ja, jb, pa, pb))
        # This same-surface balance is descriptive: for passive states it
        # can include material operations rather than a measurement update.
        _, fa, _ = common(ra, 'beforeGeo', 'filtered')
        _, fb, _ = common(rb, 'afterGeo', 'filtered')
        item['predicted_to_filtered_change'] = matrix_metrics(fa-pa, fb-pb, pa, pb)
        rows.append(item)
        previous, previous_stage = (ra, rb), 'filtered'
    return rows, before, after


def analyze(root, output, surface_catalog_path=None):
    manifest = json.loads((root/'manifest.json').read_text())
    surface_catalog = json.loads(surface_catalog_path.read_text()) if surface_catalog_path else {}
    for mode in ('off', 'on'):
        assert json.loads((root/('gate_'+mode+'_v1')/'validation.json').read_text())['all_exact']
    points = {(p['dsid'], p['variant']): p for p in manifest['points']}
    summary, all_rows, first_crossings, extremes, coverage = {}, [], [], [], []
    for dsid in SAMPLES:
        folders = {v: Path(points[dsid, v]['on']) for v in ('beforeGeo', 'afterGeo')}
        before, after = (components(folders[v]) for v in ('beforeGeo', 'afterGeo'))
        assert before.keys() == after.keys()
        counts = Counter(input_components=len(before))
        sample, paths = [], {}
        for key, a in before.items():
            b = after[key]
            ia, ib = first(a), first(b)
            counts['before_selected'] += ia is not None
            counts['after_selected'] += ib is not None
            counts['neither_selected'] += ia is None and ib is None
            counts['one_sided_selected'] += (ia is None) != (ib is None)
            if ia is None or ib is None:
                continue
            counts['common_selected'] += 1
            if a['selected'][ia]['hits'] != b['selected'][ib]['hits']:
                counts['different_ordered_hits'] += 1
                continue
            counts['same_ordered_hits'] += 1
            rows, raw_a, raw_b = paired_path(a, b, ia, ib,
                dict(dsid=dsid, entry=key[0], component=key[1], candidate_before=ia, candidate_after=ib), surface_catalog)
            counts['terminal_only'] += len(rows) == 1
            counts['measurement_states'] += sum(r['measurement_flag'] for r in rows)
            counts['passive_states'] += sum(not r['measurement_flag'] for r in rows)
            counts['material_states'] += sum(r['material']['before']['material_full_X0'] is not None for r in rows)
            counts['different_raw_material_slabs'] += sum(not r['material']['raw_slab_equal'] for r in rows)
            counts['paired_states'] += len(rows)
            counts['position_measured_states'] += sum(r['surface_type'] is not None for r in rows)
            counts['position_unmeasured_states'] += sum(r['surface_type'] is None for r in rows)
            counts['negative_path_steps'] += sum(r['path_step_mm'] is not None and r['path_step_mm'] < -1e-8 for r in rows)
            paths[key] = dict(rows=rows, before_seed=a['seed'], after_seed=b['seed'],
                              before_states=raw_a, after_states=raw_b)
            for threshold in THRESHOLDS:
                found = next(((r, stage) for r in rows for stage in ('predicted', 'filtered')
                    if r[stage]['direction_scaled_max'] > threshold), None)
                if found:
                    r, stage = found
                    first_crossings.append(dict(threshold=threshold, first_stage=stage, **r))
            sample.extend(rows)
        fields = ('direction_scaled_max', 'direction_absolute_max_rad2', 'direction_frobenius_rad2',
                  'physical_direction_C_scaled_max', 'position_difference_mm',
                  'direction_difference_norm', 'qop_difference_GeVinv')
        cohorts = {}
        for name, predicate in [('all', lambda r: True), ('measurement', lambda r: r['measurement_flag']),
                                ('passive', lambda r: not r['measurement_flag'])]:
            selected = [r for r in sample if predicate(r)]
            cohorts[name] = dict(states=len(selected), stages={stage: {
                field: metrics([r[stage][field] for r in selected if r[stage][field] is not None])
                for field in fields} for stage in STAGES})
        summary[dsid] = dict(counts=counts, cohorts=cohorts,
            link_terms={term: metrics([r['link'][term]['direction_scaled_max'] for r in sample])
                        for term in ('inherited_C', 'saved_J', 'remaining_balance')})
        selected_keys = {(1, 24)} if dsid == '100117' else set()
        for stage in STAGES:
            worst = max(sample, key=lambda r: r[stage]['direction_scaled_max'])
            selected_keys.add((worst['entry'], worst['component']))
        for key in sorted(selected_keys):
            extremes.append(dict(dsid=dsid, entry=key[0], component=key[1], **paths[key]))
        coverage.append(dict(dsid=dsid, source_before=str(folders['beforeGeo']),
            source_after=str(folders['afterGeo']), **dict(counts)))
        all_rows.extend(sample)
    output.mkdir(parents=True, exist_ok=False)
    for name, rows in [('paired_forward_states', all_rows), ('first_forward_crossings', first_crossings)]:
        with gzip.open(output/(name+'.jsonl.gz'), 'wt', compresslevel=4) as stream:
            for row in rows:
                stream.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
    write_json(output/'extreme_forward_paths.json', extremes)
    write_json(output/'summary.json', dict(groups=summary, source_events=60,
        surface_catalog=str(surface_catalog_path) if surface_catalog_path else None,
        scope='All stored states of first-ranked common ordered-hit native Calypso paths; '
              'isolated smoothing ON, n32 C/256 analytic input. Common bound angular chart and '
              'physical direction metrics are separate. Position requires native surface type; '
              'unknown positions are unmeasured, never interpreted as planar. '
              'Covariance balance is not measured native Q. '
              'No new tracking or production-rate estimate; full posterior remains open.'))
    write_json(output/'validation.json', dict(coverage=coverage,
        paired_states=len(all_rows), full_stream_replay_gates_passed=True,
        position_measured_states=sum(r['surface_type'] is not None for r in all_rows),
        position_unmeasured_states=sum(r['surface_type'] is None for r in all_rows),
        linked_forward_state_indices_and_surface_rotation_checked=True,
        different_raw_material_slabs=sum(not r['material']['raw_slab_equal'] for r in all_rows),
        covariance_balance_algebra_checked=True,
        maximum_decomposition_scaled_residual=max(r['link']['decomposition_scaled_residual'] for r in all_rows)))
    print(json.dumps(dict(paired_forward_states=len(all_rows), first_crossings=len(first_crossings),
                          extreme_paths=len(extremes))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--surface-catalog', type=Path)
    args = parser.parse_args()
    analyze(args.input, args.output, args.surface_catalog)


if __name__ == '__main__':
    main()
