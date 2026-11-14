#!/usr/bin/env python3
"""Validate complete component inputs and actual native Calypso CKF audits."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path

import numpy as np


def direction(phi, theta):
    return np.array([np.cos(phi) * np.sin(theta),
                     np.sin(phi) * np.sin(theta), np.cos(theta)])


def input_rows(path):
    rows = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split()
        assert len(fields) == 61
        entry, seed, component = map(int, fields[:3])
        values = np.asarray(list(map(float, fields[3:])))
        key = entry, component
        assert key not in rows
        rows[key] = dict(seed=seed, bound=values[:6], cov=values[6:42].reshape(6, 6),
                         surface=values[42:].reshape(4, 4))
    return rows


def validate(folder):
    assert (folder/'exit_code.txt').read_text().strip() == '0'
    task = json.loads((folder/'source_task.json').read_text())
    run = json.loads((folder/'run.json').read_text())
    assert run['success'] and run['isolate_smoothing'] and run['prior_component_audit']
    assert run['seed_component_input_file'] == task['component_file']
    assert run.get('component_refit', False) == task.get('component_refit', False)
    expected_libraries = {'libCalypsoSlopeDiagnostics.so', 'libFaserActsGeometry.so',
                          'libFaserActsGeometryLib.so', 'libFaserActsKalmanFilter.so'}
    loaded = [Path(name) for name in run['loaded_acts_libraries']]
    assert {name.name for name in loaded} == expected_libraries
    assert all(name.parent.resolve() == Path(run['library']).resolve() for name in loaded)
    options = task.get('options', [])
    fixed_fitter = '--component-fixed-hit-fitter' in options
    assert run.get('component_fixed_hit_fitter', False) == fixed_fitter
    assert run.get('actual_component_fixed_hit_fitter', False) == fixed_fitter
    observations = {}
    if fixed_fitter:
        for line in Path(run['source_link_selection_file']).read_text().splitlines():
            entry, count, *hits = map(int, line.split())
            assert count == len(hits) == len(set(hits))
            observations[entry] = set(hits)
    assert run['refit_update_chart_regular'] == ('--regular-refit-update-chart' in options)
    assert run['refit_smoothing_chart_regular'] == ('--regular-refit-smoothing-chart' in options)
    assert run.get('ckf_update_chart_regular', False) == ('--regular-ckf-update-chart' in options)
    assert run.get('ckf_smoothing_chart_regular', False) == ('--regular-ckf-smoothing-chart' in options)
    expected = input_rows(Path(task['component_file']))
    seeds, terminals, stages = set(), set(), Counter()
    found, selected, converted, finals = Counter(), Counter(), Counter(), Counter()
    fitter_inputs = set()
    events, empty, zero_native = [], [], []
    current = None
    updates, maximum_error = 0, 0.
    with gzip.open(folder/'prior_audit.jsonl.gz', 'rt') as stream:
        for line in stream:
            row = json.loads(line)
            stage = row['stage']
            stages[stage] += 1
            if stage == 'event':
                current = row.get('local_event_index', row['event'])
                events.append(current)
                if row['seeds'] == 0:
                    zero_native.append(current)
            if stage == 'component_event_end':
                assert row['components'] == 0 and not any(e == current for e, _ in expected)
                empty.append(current)
            if 'input_component' not in row:
                continue
            key = current, row['input_component']
            assert key in expected and row['seed'] == expected[key]['seed']
            if stage == 'seed':
                assert key not in seeds
                seeds.add(key)
                # Acts canonicalizes bound (phi,theta) at construction.  A
                # representation such as (phi+pi,-theta) is identical as a
                # physical direction, so compare the invariant direction and
                # the remaining injected quantities rather than raw angles.
                np.testing.assert_allclose(row['bound'][:2], expected[key]['bound'][:2],
                                           rtol=0, atol=2e-12)
                np.testing.assert_allclose(direction(row['bound'][2], row['bound'][3]),
                                           direction(expected[key]['bound'][2], expected[key]['bound'][3]),
                                           rtol=0, atol=2e-12)
                np.testing.assert_allclose(row['bound'][4:], expected[key]['bound'][4:],
                                           rtol=0, atol=2e-12)
                np.testing.assert_allclose(row['surface'], expected[key]['surface'],
                                           rtol=0, atol=2e-12)
                np.testing.assert_allclose(row['cov'], expected[key]['cov'],
                                           rtol=2e-14, atol=1e-18)
            elif stage == 'found':
                found[key] += 1
            elif stage == 'selected':
                selected[key] += 1
            elif stage == 'converted':
                converted[key] += 1
            elif stage == 'final':
                finals[key] += 1
            elif stage == 'component_fixed_observations':
                assert fixed_fitter and key not in fitter_inputs
                assert row['count'] == len(row['hit_ids']) == len(set(row['hit_ids']))
                assert set(row['hit_ids']) == observations[current], key
                fitter_inputs.add(key)
            elif stage in ['component_end', 'component_finder_error']:
                assert key not in terminals
                terminals.add(key)
                if stage == 'component_end':
                    assert row['found'] == found[key] and row['selected'] == selected[key]
                elif fixed_fitter:
                    assert row['error_category'] and row['error_message']
            elif stage == 'selected_state' and row['measurement_flag']:
                predicted = np.asarray(row['predicted'])
                covariance = np.asarray(row['predicted_cov'])
                projector = np.asarray(row['projector'])
                measurement = np.asarray(row['measurement'])
                variance = np.asarray(row['measurement_cov'])
                residual = measurement - projector @ predicted
                innovation = projector @ covariance @ projector.T + variance
                gain = np.linalg.solve(innovation, projector @ covariance).T
                expected_qop = predicted[4] + (gain @ residual)[4]
                error = abs(row['filtered'][4] - expected_qop)
                maximum_error = max(maximum_error, error)
                assert error < 1e-11 * max(1., abs(expected_qop)), (key, row['state'], error)
                updates += 1
    assert events == list(range(task['events']))
    assert seeds == terminals == set(expected)
    assert fitter_inputs == (set(expected) if fixed_fitter else set())
    assert empty == [e for e in events if not any(k[0] == e for k in expected)]
    if run.get('component_refit', False):
        assert converted == finals == selected
    else:
        assert not converted and not finals
    return dict(dsid=task['dsid'], variant=task['variant'], source_events=len(events),
                input_components=len(expected), complete_mean_covariance_surface_exact=True,
                all_components_have_native_terminal_record=True,
                zero_native_seed_events=zero_native, empty_component_events=empty,
                native_selected_candidates=sum(selected.values()),
                native_conversion_refit_candidates=sum(finals.values()),
                native_finder_errors=stages['component_finder_error'],
                actual_loaded_acts_and_observer_libraries_checked=True,
                configured_component_input_refit_and_chart_checked=True,
                actual_native_solver='fixed_source_link_KalmanFitter' if fixed_fitter else 'CKF_finder',
                actual_fitter_observation_sets_checked=len(fitter_inputs),
                measurement_updates_checked=updates, maximum_gain_identity_error_GeVinv=maximum_error,
                stages=dict(stages), run=str(folder/'run.json'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checks = []
    for variant in ['beforeGeo', 'afterGeo']:
        for job in json.loads((args.input/(variant+'_tasks.json')).read_text()):
            checks.append(validate(Path(job['output'])))
    args.output.write_text(json.dumps(dict(all_complete=True, checks=checks,
        scope='Input, event coverage, native component completion and measurement-update gate; '
              'does not establish finite-prior equivalence or an aggregated posterior'), indent=2)+'\n')
    print(json.dumps(dict(jobs=len(checks), components=sum(r['input_components'] for r in checks),
                         updates=sum(r['measurement_updates_checked'] for r in checks))))


if __name__ == '__main__':
    main()
