#!/usr/bin/env python3
"""Measure paired native Calypso reconstructions of symmetric full-G4 inputs."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta

from charge_analysis import paired_difference, summarize
from charge_modes import compact
from charge_tables import STATUS, derive, measurement
from slope_tables import histogram, merge_tables, slope, write_json


MODES = ['beforeGeo', 'beforeGeo_regular_both', 'afterGeo']
DSIDS = ['100116', '100117']


def charge_difference(plus, minus):
    """Bonferroni binomial bounds retain uncertainty for two zero counts."""
    def interval(row):
        k, n = row['wrong'], row['total']
        if not n:
            return None
        return [float(beta.ppf(.0125, k, n-k+1)) if k else 0.,
                float(beta.ppf(.9875, k+1, n-k)) if k < n else 1.]
    a, b = interval(plus), interval(minus)
    return dict(value=minus['value']-plus['value'] if a and b else None,
        interval95=[b[0]-a[1], b[1]-a[0]] if a and b else [None, None],
        method='Conservative difference of simultaneous 97.5% Clopper-Pearson intervals',
        mu_plus=plus, mu_minus=minus)


def refit_counts(jobs, tables, decision):
    selected = set(zip(tables['tracks']['sample_entry'][decision['chosen'][decision['selected']]],
                       tables['tracks']['ordinal'][decision['chosen'][decision['selected']]]))
    counts = Counter()
    for job in jobs:
        folder = Path(job['output'])
        task = json.loads((folder/'source_task.json').read_text())
        with np.load(folder/'charge_lineage.npz') as rows:
            final = rows['stage'] == 7
            entries = rows['source_entry'][final]-task['skip']+task['sample_entry']
            chosen = np.array([(entry, ordinal) in selected
                               for entry, ordinal in zip(entries, rows['ordinal'][final])], dtype=bool)
            qop = rows['qop_GeVinv'][final]
            converted = rows['qop_GeVinv'][rows['stage'] == 6]
            assert len(qop) == len(converted)
            status = rows['status'][final]
            charge = decision['charge']
            for label, mask in [('all_final', np.ones(len(qop), dtype=bool)), ('selected', chosen)]:
                counts[label+'_total'] += int(mask.sum())
                for outcome, value in [('success', 1), ('fallback', 0), ('invalid', -1)]:
                    counts[label+'_'+outcome] += int(np.count_nonzero(mask & (status == value)))
                counts[label+'_correct_to_wrong'] += int(np.count_nonzero(mask & (charge*converted > 0) & (charge*qop < 0)))
                counts[label+'_wrong_to_correct'] += int(np.count_nonzero(mask & (charge*converted < 0) & (charge*qop > 0)))
    assert counts['selected_total'] == int(decision['selected'].sum())
    return dict(counts)


def angular_histograms(tables, decision):
    tracks, truth = tables['tracks'], tables['truth']
    rows = decision['chosen'][decision['selected']]
    p = tracks['projection'][rows, 0]
    valid = tracks['projection_valid'][rows, 0] == 1
    launch = slope(truth['production_momentum_F_GeV'][decision['selected']])
    sy = slope(p[:, 3:6])
    residual = sy-launch
    finite = valid & np.isfinite(residual)
    pull_valid = valid & np.isfinite(p[:, 6]) & (p[:, 6] > 0)
    moments = dict(finite=int(finite.sum()),
        quantile_probabilities=[.025, .16, .5, .84, .975],
        residual_quantiles=np.quantile(residual[finite], [.025, .16, .5, .84, .975]).tolist() if finite.any() else None,
        absolute_residual_68=float(np.quantile(np.abs(residual[finite]), .68)) if finite.any() else None,
        absolute_residual_above_1mrad=int(np.count_nonzero(finite & (np.abs(residual) > .001))))
    return dict(st1_sy=histogram(sy[valid], np.linspace(-.005, .005, 201)),
        st1_sy_minus_launch=histogram(residual[valid], np.linspace(-.005, .005, 201)),
        st1_sy_minus_launch_fine=histogram(residual[valid], np.linspace(-.0005, .0005, 201)),
        residual_summary=moments,
        st1_sy_pull=histogram(residual[pull_valid]/np.sqrt(p[pull_valid, 6]), np.linspace(-20, 20, 201)),
        invalid_pull_variance=int(np.count_nonzero(valid & ~pull_valid)),
        failed_st1_projections=int((~valid).sum()),
        reference='Physical sy projected by actual Calypso to z_F=0 mm; exact G4 launch at the same plane')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', default='symmetric_gun_20k_v1')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source = args.root/args.source
    manifest = json.loads((source/'manifest.json').read_text())
    tasks = manifest['source_tasks']
    generated, results, closures, outcomes, launch_precision = {}, {}, {}, {}, {}
    for dsid in DSIDS:
        momenta = sorted({task['p_GeV'] for task in tasks if task['dsid'] == dsid})
        for momentum in momenta:
            sources = [task for task in tasks if task['dsid'] == dsid and task['p_GeV'] == momentum]
            sources.sort(key=lambda task: task['sample_entry'])
            count = sum(task['events'] for task in sources)
            assert count == manifest['events_per_charge_and_momentum']
            for task in sources:
                folder = Path(task['output'])
                for label in ['', 'simulate', 'digitize', *[mode+'_wrapper' for mode in MODES]]:
                    assert (folder/label/'exit_code.txt').read_text().strip() == '0'
            reference_truth = None
            for mode in MODES:
                variant = 'afterGeo' if mode == 'afterGeo' else 'beforeGeo'
                jobs = [dict(output=str(Path(task['output'])/mode),
                             exit_code_file=str(Path(task['output'])/(mode+'_wrapper')/'exit_code.txt'))
                        for task in sources]
                reference = args.output/'input_references'/f'{dsid}_p{momentum:g}'/mode
                reference.mkdir(parents=True)
                write_json(reference/(variant+'_tasks.json'), jobs)
                write_json(reference/'manifest.json', dict(source_events_per_variant=count,
                    runtime_invariants=dict(audit=False, reset_event_maps=False, isolate_smoothing=True)))
                tables, coverage = merge_tables(reference, variant)
                decision = derive(tables, dsid)
                primary = decision['primary']
                assert primary.sum() == count
                np.testing.assert_array_equal(tables['truth']['sample_entry'][primary], np.arange(count))
                truth = {field: tables['truth'][field][primary] for field in [
                    'sample_entry', 'production_position_F_mm', 'production_momentum_F_GeV',
                    'truth_hit_counts', 'truth_sensor_direction_F', 'truth_station_position_F_mm']}
                if reference_truth is None:
                    reference_truth = truth
                    magnitude = np.linalg.norm(truth['production_momentum_F_GeV'].astype(float), axis=1)
                    relative_error = np.max(np.abs(magnitude/momentum-1.))
                    assert relative_error < 1e-6
                    launch_precision[f'{dsid}|{momentum:g}'] = dict(
                        requested_momentum_GeV=momentum,
                        observed_momentum_range_GeV=[float(magnitude.min()), float(magnitude.max())],
                        maximum_relative_difference=float(relative_error),
                        nominal_momentum_relative_budget=1e-6,
                        reference='Actual stored generated momentum is used in residuals; charge-to-charge equality is checked exactly')
                    np.testing.assert_allclose(truth['production_position_F_mm'][:, 2], 0., atol=1e-12)
                    generated[dsid, momentum] = truth
                else:
                    for field in truth:
                        np.testing.assert_array_equal(truth[field], reference_truth[field], err_msg=field)
                expected_software = 'software_afterGeo' if variant == 'afterGeo' else 'software_refit_before'
                expected = {str((args.root/expected_software/'lib'/name).resolve()) for name in [
                    'libFaserActsGeometry.so', 'libFaserActsGeometryLib.so', 'libFaserActsKalmanFilter.so']}
                expected.add(str((args.root/'observer/lib/libCalypsoSlopeDiagnostics.so').resolve()))
                for job in jobs:
                    run = json.loads((Path(job['output'])/'run.json').read_text())
                    assert set(run['loaded_acts_libraries']) == expected
                key = f'{dsid}|{momentum:g}|{mode}'
                result = summarize(tables, decision)
                result['angular_histograms'] = angular_histograms(tables, decision)
                result['refit_counts'] = refit_counts(jobs, tables, decision)
                result['coverage'] = coverage
                result['reconstruction_wall_seconds_per_source_event'] = sum(r['wall_s'] for r in coverage['resources'])/count
                position = tables['truth']['production_position_F_mm']
                result['production_regions'] = {}
                for xpos in [False, True]:
                    for ypos in [False, True]:
                        mask = ((position[:, 0] >= 0) == xpos) & ((position[:, 1] >= 0) == ypos)
                        result['production_regions'][f'x{"plus" if xpos else "minus"}_y{"plus" if ypos else "minus"}'] = measurement(tables, decision, mask)
                results[key] = result
                current = compact(tables, decision)
                outcomes[key] = current
                np.savez_compressed(args.output/(key.replace('|', '_')+'_outcomes.npz'), **current)
                if mode != 'beforeGeo':
                    baseline = outcomes[f'{dsid}|{momentum:g}|beforeGeo']
                    closures[key] = dict(states=STATUS, counts=np.bincount(
                        7*baseline['status']+current['status'], minlength=49).reshape(7, 7).tolist(),
                        conditional_difference=paired_difference(baseline, current, np.ones(count, dtype=bool)))
                print(key, result['conditional_misid'], flush=True)
                del tables, decision
    differences = {}
    for momentum in sorted({task['p_GeV'] for task in tasks}):
        a, b = [generated[dsid, momentum] for dsid in DSIDS]
        for field in ['sample_entry', 'production_position_F_mm', 'production_momentum_F_GeV']:
            np.testing.assert_array_equal(a[field], b[field], err_msg='Actual mu+/- launch: '+field)
        for mode in MODES:
            plus, minus = [results[f'{dsid}|{momentum:g}|{mode}']['conditional_misid'] for dsid in DSIDS]
            differences[f'{momentum:g}|{mode}'] = charge_difference(plus, minus)
    write_json(args.output/'summary.json', dict(groups=results, migrations=closures,
        charge_differences=differences, all_actual_generated_conditions_exact=True,
        nominal_launch_precision=launch_precision,
        all_versions_same_truth=True, all_actual_libraries_match=True,
        events_per_charge_momentum=manifest['events_per_charge_and_momentum'],
        scope='Independent full Calypso G4 and digitization; isolated smoothing ON in all three real reconstructions',
        interpretation='Generated launch spectra are identical; actual G4 interactions and detector acceptance can differ by charge. Spatial regions test the real apparatus without assuming exact reflection symmetry.'))


if __name__ == '__main__':
    main()
