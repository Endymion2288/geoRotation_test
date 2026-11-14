#!/usr/bin/env python3
"""Read fixed-observation innovation likelihoods from actual native CKF states.

No propagation, update, smoothing or refit is implemented here. A posterior is
only normalized when every input component has one complete observation path.
"""
import argparse
from collections import Counter, defaultdict
from itertools import zip_longest
import gzip
import json
import math
from pathlib import Path

import numpy as np

from charge_component_refit import audit, components, metrics
from charge_forward_covariance import common, forward_states
from bound_surface_coordinates import ROTATE_F_TO_A
from slope_tables import write_json
from validate_charge_components import input_rows, validate


SAMPLES = ('100044', '100116', '100117')
VARIANTS = ('beforeGeo', 'afterGeo')


def innovation(row):
    """Normalized Gaussian innovation density in actual measurement units."""
    y = np.asarray(row['measurement'], dtype=float)
    h = np.asarray(row['projector'], dtype=float)
    mean = np.asarray(row['predicted'], dtype=float)
    cov = np.asarray(row['predicted_cov'], dtype=float)
    variance = np.asarray(row['measurement_cov'], dtype=float)
    residual = y-h @ mean
    s = h @ cov @ h.T+variance
    assert np.isfinite(residual).all() and np.isfinite(s).all()
    np.testing.assert_allclose(s, s.T, rtol=2e-12, atol=1e-18)
    chol = np.linalg.cholesky(s)
    pull = np.linalg.solve(chol, residual)
    chi2 = float(pull @ pull)
    logdet = float(2*np.log(np.diag(chol)).sum())
    return dict(chi2=chi2, logdet=logdet, dimension=len(y),
        log_likelihood=-.5*(chi2+logdet+len(y)*math.log(2*math.pi)),
        residual=residual.tolist(), innovation_covariance=s.tolist())


def normalized_weights(log_likelihoods):
    """Equal finite-prior component masses; stable even for very small L."""
    values = np.asarray(log_likelihoods, dtype=float)
    assert len(values) and np.isfinite(values).all()
    maximum = float(values.max())
    shifted = np.exp(values-maximum)
    weights = shifted/shifted.sum()
    evidence = maximum+math.log(float(shifted.sum()))-math.log(len(values))
    return weights, evidence


def mixture(means, covariances, weights):
    """Total C includes native within-component C and between-component spread."""
    means, covariances = np.asarray(means), np.asarray(covariances)
    weights = np.asarray(weights)
    np.testing.assert_allclose(weights.sum(), 1., rtol=0, atol=2e-15)
    mean = weights @ means
    within = np.einsum('i,ijk->jk', weights, covariances)
    delta = means-mean
    between = np.einsum('i,ij,ik->jk', weights, delta, delta)
    return dict(mean=mean.tolist(), within_cov=within.tolist(),
        between_cov=between.tolist(), cov=(within+between).tolist())


def component_path(item, expected_hits, variant):
    """Keep all failures/multiplicity; never choose a surviving first track."""
    result = dict(selected_candidates=len(item['selected']), complete=False)
    if len(item['selected']) != 1:
        return dict(result, reason='no_selected' if not item['selected'] else 'multiple_selected'), None
    candidate = next(iter(item['selected']))
    path = forward_states(item, candidate)
    measurements = [r for r in path if r['measurement_flag']]
    hits = [r['hit'] for r in measurements]
    result.update(candidate=candidate, states=len(path), measurements=len(measurements),
        terminal_only=len(path) == 1, actual_hits=hits,
        outliers=sum(r['outlier_flag'] for r in path), holes=sum(r['hole_flag'] for r in path),
        refit_success=item['final'].get(candidate, {}).get('refit') == 1)
    if hits != expected_hits or result['outliers']:
        return dict(result, reason='different_observation_path'), None
    if not measurements:
        return dict(result, reason='empty_observation_path'), None
    try:
        innovations = [innovation(row) for row in measurements]
    except (AssertionError, np.linalg.LinAlgError, ValueError):
        return dict(result, reason='invalid_innovation'), None
    mean, cov, _ = common(measurements[0], variant, 'smoothed')
    if not (np.isfinite(mean).all() and np.isfinite(cov).all()):
        return dict(result, reason='invalid_smoothed_gaussian'), None
    result.update(complete=True, reason='complete',
        log_likelihood=sum(r['log_likelihood'] for r in innovations),
        innovation_chi2=sum(r['chi2'] for r in innovations),
        innovation_dimensions=sum(r['dimension'] for r in innovations),
        native_chi2=item['selected'][candidate]['chi2'],
        smoothed_mean=mean.tolist(), smoothed_cov=cov.tolist(),
        innovations=innovations)
    final = item['final'].get(candidate)
    if final is not None:
        result.update(final_qop_GeVinv=final['bound_trk'][4]*1000.,
            final_qop_variance_GeVinv2=final['variance_qop'])
    return result, measurements


def observation_signature(row):
    return {k:row[k] for k in ('geometry_id', 'hit', 'measurement', 'measurement_cov', 'projector')}


def event_posterior(rows, truth_charge):
    complete = [row for row in rows if row['complete']]
    result = dict(input_components=len(rows), complete_components=len(complete),
        missing_prior_mass=(len(rows)-len(complete))/len(rows),
        complete=len(rows) == len(complete), posterior=None)
    if not result['complete']:
        return result
    weights, log_evidence = normalized_weights([row['log_likelihood'] for row in rows])
    gaussian = mixture([r['smoothed_mean'] for r in rows],
        [r['smoothed_cov'] for r in rows], weights)
    probability = 0.
    for w, row in zip(weights, rows):
        variance = row['smoothed_cov'][4][4]
        assert variance > 0
        probability += w*.5*math.erfc(truth_charge*row['smoothed_mean'][4]/math.sqrt(2*variance))
        row['posterior_weight'] = float(w)
    result['posterior'] = dict(gaussian, log_evidence=log_evidence,
        effective_components=float(1/(weights @ weights)), max_weight=float(weights.max()),
        wrong_charge_probability=float(probability),
        qop_mean_GeVinv=gaussian['mean'][4], qop_variance_GeVinv2=gaussian['cov'][4][4])
    # The native refit reuses the observations and inflates C. Its output is
    # recorded separately, not treated as a second Bayesian likelihood factor.
    if all('final_qop_GeVinv' in r for r in rows):
        qop = np.array([r['final_qop_GeVinv'] for r in rows])
        variance = np.array([r['final_qop_variance_GeVinv2'] for r in rows])
        mean = float(weights @ qop)
        result['native_final_output_mixture'] = dict(qop_mean_GeVinv=mean,
            qop_variance_GeVinv2=float(weights @ (variance+(qop-mean)**2)),
            refit_success_prior_fraction=sum(r['refit_success'] for r in rows)/len(rows),
            refit_success_CKF_posterior_weight=float(sum(w for w, r in zip(weights, rows) if r['refit_success'])))
    return result


def gate(root, output):
    checks = []
    for point in json.loads((root/'manifest.json').read_text())['points']:
        if point['mode'] != 'default':
            continue
        old, new = Path(point['reference']), Path(point['output'])
        counts = {}
        for name in ('prior_audit', 'smoothing_operations'):
            count = 0
            for a, b in zip_longest(audit(old/(name+'.jsonl.gz')), audit(new/(name+'.jsonl.gz'))):
                assert a == b, (point['dsid'], point['variant'], name, count)
                count += 1
            counts[name+'_rows_exact'] = count
        arrays = 0
        for name in ('events', 'truth', 'tracks', 'segments', 'refit_outcomes'):
            if not (old/(name+'.npz')).exists() and not (new/(name+'.npz')).exists():
                continue
            with np.load(old/(name+'.npz')) as a, np.load(new/(name+'.npz')) as b:
                assert a.files == b.files
                for field in a.files:
                    np.testing.assert_array_equal(a[field], b[field])
                    arrays += 1
        run = json.loads((new/'run.json').read_text())
        assert run['ckf_chi2_max'] is None and run['ckf_nmax'] is None
        # CKF2Config overrides the C++ property's 15 with 25 in these runs.
        assert run['actual_ckf_chi2_max'] == 25. and run['actual_ckf_nmax'] == 10
        counts.update(array_fields_exact=arrays)
        checks.append(dict(dsid=point['dsid'], variant=point['variant'], **counts, runtime=validate(new)))
    output.mkdir(exist_ok=False)
    write_json(output/'validation.json', dict(default_replay_exact=True, checks=checks))
    print(json.dumps(dict(default_replay_exact=True, jobs=len(checks))))


def analyze(root, output, gate_validation):
    assert json.loads(gate_validation.read_text())['default_replay_exact']
    manifest = json.loads((root/'manifest.json').read_text())
    observations = json.loads((root/'observations.json').read_text())
    points = {(p['dsid'], p['variant'], p['mode']): p for p in manifest['points']}
    groups, rows_out, event_rows, runtime, comparison_rows = {}, [], [], [], []
    for dsid in SAMPLES:
        common_observations, common_surfaces = {}, {}
        event_results, component_results = {}, {}
        for variant in VARIANTS:
            point = points[dsid, variant, 'fixed_observation']
            folder = Path(point['output'])
            runtime.append(validate(folder))
            run = json.loads((folder/'run.json').read_text())
            assert run['actual_ckf_chi2_max'] == 1e30 and run['actual_ckf_nmax'] == 1
            assert run['actual_ckf_max_step_mm'] == 10.
            assert run['source_link_selection_file'] == str(root/'observations'/(dsid+'_hits.txt'))
            records = components(folder)
            task = json.loads((folder/'source_task.json').read_text())
            assert records.keys() == input_rows(Path(task['component_file'])).keys()
            group, by_event = Counter(input_components=len(records)), defaultdict(list)
            for key, item in records.items():
                expected = observations[dsid]['events'][key[0]]['forward_ordered_hits']
                result, measurements = component_path(item, expected, variant)
                row = dict(dsid=dsid, variant=variant, entry=key[0], component=key[1],
                    truth_charge=item['truth_charge'], **result)
                rows_out.append(row)
                by_event[key[0]].append(row)
                component_results[variant, key] = row
                group[result['reason']] += 1
                group['selected_candidates'] += len(item['selected'])
                group['terminal_only'] += result.get('terminal_only', False)
                group['refit_success'] += result.get('refit_success', False)
                if measurements is not None:
                    signature = list(map(observation_signature, measurements))
                    if key[0] not in common_observations:
                        common_observations[key[0]] = signature
                    assert signature == common_observations[key[0]], (dsid, variant, key)
                    surface = np.asarray(measurements[0]['surface'])
                    if variant == 'afterGeo':
                        rotate = np.eye(4)
                        rotate[:3, :3] = ROTATE_F_TO_A.T
                        surface = rotate @ surface
                    if key[0] not in common_surfaces:
                        common_surfaces[key[0]] = surface
                    np.testing.assert_allclose(surface, common_surfaces[key[0]], rtol=0, atol=2e-9)
            for entry in range(20):
                rows = by_event[entry]
                if rows:
                    assert len({r['truth_charge'] for r in rows}) == 1
                    result = event_posterior(rows, rows[0]['truth_charge'])
                    event_results[variant, entry] = result
                    event_rows.append(dict(dsid=dsid, variant=variant, entry=entry, **result))
                    group['complete_event_posteriors'] += result['complete']
                    group['incomplete_event_posteriors'] += not result['complete']
                else:
                    assert not observations[dsid]['events'][entry]['forward_ordered_hits']
                    group['events_without_prior_components'] += 1
                    event_rows.append(dict(dsid=dsid, variant=variant, entry=entry, input_components=0,
                        complete=False, posterior=None, reason='no_native_seed'))
            groups[dsid+'|'+variant] = dict(group)
        assert {key for v, key in component_results if v == 'beforeGeo'} == {
            key for v, key in component_results if v == 'afterGeo'}
        paired = []
        for entry in range(20):
            if ('beforeGeo', entry) not in event_results:
                continue
            a, b = (event_results[v, entry] for v in VARIANTS)
            row = dict(dsid=dsid, entry=entry, complete_both=a['complete'] and b['complete'],
                before_complete_components=a['complete_components'], after_complete_components=b['complete_components'],
                input_components=a['input_components'])
            if row['complete_both']:
                pa, pb = a['posterior'], b['posterior']
                for name in ('qop_mean_GeVinv', 'qop_variance_GeVinv2', 'log_evidence',
                             'wrong_charge_probability', 'effective_components'):
                    row[name+'_difference'] = pb[name]-pa[name]
                weights = {v:np.array([component_results[v, (entry, i)]['posterior_weight']
                    for i in sorted(k[1] for vv, k in component_results if vv == v and k[0] == entry)]) for v in VARIANTS}
                row['weight_total_variation'] = float(.5*np.abs(weights['afterGeo']-weights['beforeGeo']).sum())
                row['qop_mean_difference_over_posterior_sigma'] = abs(row['qop_mean_GeVinv_difference'])/math.sqrt(
                    (pa['qop_variance_GeVinv2']+pb['qop_variance_GeVinv2'])/2.)
                paired.append(row)
            comparison_rows.append(row)
        groups[dsid+'|comparison'] = dict(input_events=len(event_results)//2, complete_both=len(paired),
            differences={name:metrics([r[name] for r in paired]) for name in
                ('qop_mean_GeVinv_difference', 'log_evidence_difference', 'wrong_charge_probability_difference',
                 'weight_total_variation', 'qop_mean_difference_over_posterior_sigma')})
    output.mkdir(exist_ok=False)
    for filename, rows in [('components.jsonl.gz', rows_out), ('events.jsonl.gz', event_rows),
                           ('comparisons.jsonl.gz', comparison_rows)]:
        with gzip.open(output/filename, 'wt', compresslevel=4) as stream:
            for row in rows:
                stream.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
    write_json(output/'summary.json', dict(groups=groups, runtime=runtime,
        source_events=60, finite_series=manifest['component_series'],
        scope='Actual native CKF Gaussian innovation likelihood for a common frozen observation sequence. '
              'All finite prior components retained, equal prior component masses. Full posterior only when '
              'every component completes the observations. Gaussian moment C includes within/between terms '
              'at the same first sensor in a common regular chart. Missing mass is not normalized away. '
              'This is a finite-mixture EKF diagnostic, not a production charge-misID rate; native refit '
              'reuses data and inflates C, so its weighted output is separately labelled.'))
    write_json(output/'validation.json', dict(default_gate_exact=True,
        common_measurement_values_covariance_projector_exact=True,
        common_first_sensor_surface_checked=True, all_input_components_retained=True,
        incomplete_events_have_no_normalized_posterior=all(r['posterior'] is None for r in event_rows if not r['complete']),
        component_records=len(rows_out), source_event_records=len(event_rows)))
    print(json.dumps(groups))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--phase', choices=['gate', 'analyze'], required=True)
    parser.add_argument('--gate-validation', type=Path)
    args = parser.parse_args()
    if args.phase == 'gate':
        gate(args.input, args.output)
    else:
        analyze(args.input, args.output, args.gate_validation)


if __name__ == '__main__':
    main()
