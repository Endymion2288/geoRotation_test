#!/usr/bin/env python3
"""Check native closure and observer invariance for actual Calypso refit controls."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from slope_controls import archive
from slope_tables import write_json


def equal_archives(left, right, names):
    checked = {}
    for name in names:
        a, b = (archive(folder / (name + '.npz')) for folder in [left, right])
        assert a.keys() == b.keys(), name
        for key in a:
            equal = np.array_equal(a[key], b[key]) if a[key].dtype.kind in 'US' else np.array_equal(a[key], b[key], equal_nan=True)
            assert equal, (name, key, left, right)
        checked[name] = len(a)
    return checked


def check_basis(args):
    modes = ['native', 'observed', 'basis', 'plane', 'basis_plane']
    folders = {mode: args.root / 'smoke' / mode for mode in modes}
    result = dict(runs={}, comparisons={})
    identity = None
    for mode, folder in folders.items():
        assert (folder / 'exit_code.txt').read_text().strip() == '0', folder
        run = json.loads((folder / 'run.json').read_text())
        task = json.loads((folder / 'source_task.json').read_text())
        source = tuple(task[k] for k in ['input', 'skip', 'events', 'variant', 'material'])
        if identity is None:
            identity = source
        assert source == identity
        expected = dict(conversion_covariance_basis_corrected=mode in ['basis', 'basis_plane'],
                        refit_covariance_basis_corrected=mode in ['basis', 'basis_plane'],
                        refit_uses_input_reference_plane=mode in ['plane', 'basis_plane'],
                        refit_outcome_observer=mode != 'native')
        for key, value in expected.items():
            assert run[key] == value, (mode, key)
        assert run['success'] and run['refit_covariance_scale'] == 10
        assert not run['refit_smoothing_phi_periodic']
        assert not run['refit_smoothing_angles_canonical']
        assert Path(run['library']).resolve() == (args.root / 'software/lib').resolve()
        result['runs'][mode] = dict(source_events=task['events'], library=run['library'], **expected)
    tables = ['events', 'truth', 'tracks', 'segments', 'native_states']
    result['comparisons']['native_vs_original'] = equal_archives(folders['native'], args.original, tables)
    result['comparisons']['observed_vs_native'] = equal_archives(
        folders['observed'], folders['native'], tables + ['native_angles'])
    base = archive(folders['native'] / 'native_states.npz')
    for mode in modes[1:]:
        checks = dict(source_invariance=equal_archives(
            folders['native'], folders[mode], ['events', 'truth', 'segments']))
        data = archive(folders[mode] / 'native_states.npz')
        masks = [np.isin(table['stage'], [1, 3, 4, 5]) for table in [base, data]]
        for key in base:
            np.testing.assert_array_equal(base[key][masks[0]], data[key][masks[1]],
                                          err_msg=mode + ' seed/CKF ' + key)
        checks['seed_CKF_exact_rows'] = int(masks[0].sum())
        outcomes = archive(folders[mode] / 'refit_outcomes.npz')
        tracks = archive(folders[mode] / 'tracks.npz')
        for key in ['source_entry', 'run', 'event', 'ordinal']:
            np.testing.assert_array_equal(outcomes[key], tracks[key], err_msg=mode + ' outcome ' + key)
        assert np.isin(outcomes['refit_success'], [0, 1]).all()
        checks['refit_outcomes'] = dict(selected=len(outcomes['ordinal']),
            succeeded=int(outcomes['refit_success'].sum()),
            fallback=int((outcomes['refit_success'] == 0).sum()))
        result['comparisons'][mode] = checks
    result['definition'] = ('Actual Calypso saved-field checks; observer leaves all native fields unchanged; '
                            'controls keep source events, truth, segments, seed and CKF states fixed. '
                            'Refit forward and final states may change.')
    return result


def check_angles(args):
    folders = {mode: args.root / 'smoke' / mode for mode in ['native', 'periodic', 'canonical']}
    runs = {}
    identity = None
    for mode, folder in folders.items():
        assert (folder / 'exit_code.txt').read_text().strip() == '0', folder
        run = json.loads((folder / 'run.json').read_text())
        task = json.loads((folder / 'source_task.json').read_text())
        source = tuple(task[k] for k in ['input', 'skip', 'events', 'variant', 'material'])
        if identity is None:
            identity = source
        assert source == identity
        assert run['success'] and run['refit_covariance_scale'] == 10
        assert run['refit_smoothing_phi_periodic'] == (mode != 'native')
        assert run['refit_smoothing_angles_canonical'] == (mode == 'canonical')
        assert Path(run['library']).resolve() == (args.root / 'software/lib').resolve()
        runs[mode] = dict(source_events=task['events'], library=run['library'],
                          periodic=run['refit_smoothing_phi_periodic'],
                          canonical=run['refit_smoothing_angles_canonical'])
    result = dict(runs=runs, comparisons={})
    tables = ['events', 'truth', 'tracks', 'segments', 'native_states']
    result['comparisons']['native_vs_original'] = equal_archives(folders['native'], args.original, tables)
    for mode in ['native', 'periodic']:
        result['comparisons'][mode + '_vs_previous'] = equal_archives(
            folders[mode], args.previous / 'smoke' / mode, tables + ['native_angles'])
    base = archive(folders['native'] / 'native_states.npz')
    for mode in ['periodic', 'canonical']:
        result['comparisons'][mode + '_source_invariance'] = equal_archives(
            folders['native'], folders[mode], ['events', 'truth', 'segments'])
        data = archive(folders[mode] / 'native_states.npz')
        for key in ['entry', 'track_index', 'stage', 'station', 'layer']:
            np.testing.assert_array_equal(base[key], data[key])
        checked = {}
        for label, stages in [('seed_CKF', [1, 3, 4, 5]), ('refit_forward', [8, 9])]:
            mask = np.isin(base['stage'], stages)
            for key in base:
                np.testing.assert_array_equal(base[key][mask], data[key][mask], err_msg=mode + ' ' + key)
            checked[label] = int(mask.sum())
        smoothed = base['stage'] == 10
        changed = ~np.all(np.isclose(base['values'][smoothed], data['values'][smoothed], rtol=0, atol=0, equal_nan=True), axis=1)
        checked.update(refit_smoothed=int(smoothed.sum()), refit_smoothed_changed=int(changed.sum()))
        result['comparisons'][mode + '_native_states'] = checked
    angles = archive(folders['canonical'] / 'native_angles.npz')
    phi, theta = angles['ePHI_smt'], angles['eTHETA_smt']
    selected = (angles['writer'] == 1) & np.isfinite(phi) & np.isfinite(theta)
    tolerance = 8 * np.finfo(np.float32).eps
    outside = selected & ((np.abs(phi) > np.pi + tolerance) | (theta < -tolerance) | (theta > np.pi + tolerance))
    assert not outside.any()
    result['canonical_writer'] = dict(finite_refit_measurement_states=int(selected.sum()),
                                      angular_chart_violations=int(outside.sum()),
                                      float_storage_tolerance_rad=tolerance)
    result['definition'] = 'Actual Calypso outputs; exact saved-field and state identity checks, no tracking recalculation'
    return result


def check_update(args):
    modes = ['native', 'canonical', 'updater', 'update_canonical']
    folders = {mode: args.root/'smoke'/mode for mode in modes}
    result = dict(runs={}, comparisons={})
    identity = None
    for mode, folder in folders.items():
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        run = json.loads((folder/'run.json').read_text())
        task = json.loads((folder/'source_task.json').read_text())
        source = tuple(task[k] for k in ['input', 'skip', 'events', 'variant', 'material'])
        if identity is None:
            identity = source
        assert source == identity
        assert run['success'] and run['refit_covariance_scale'] == 10
        assert run['refit_outcome_observer']
        assert run['refit_update_angles_canonical'] == (mode in ['updater', 'update_canonical'])
        assert run['refit_smoothing_angles_canonical'] == (mode in ['canonical', 'update_canonical'])
        assert run['refit_smoothing_phi_periodic'] == (mode in ['canonical', 'update_canonical'])
        assert not any(run[k] for k in ['conversion_covariance_basis_corrected',
                                       'refit_covariance_basis_corrected', 'refit_uses_input_reference_plane'])
        assert Path(run['library']).resolve() == (args.root/'software/lib').resolve()
        result['runs'][mode] = dict(source_events=task['events'], library=run['library'])
    tables = ['events', 'truth', 'tracks', 'segments', 'native_states', 'native_angles']
    result['comparisons']['native_vs_original'] = equal_archives(
        folders['native'], args.original, tables[:-1])
    for mode in ['native', 'canonical']:
        result['comparisons'][mode+'_vs_previous'] = equal_archives(
            folders[mode], args.previous/'smoke'/mode, tables)
    base = archive(folders['native']/'native_states.npz')
    for mode in modes[1:]:
        checks = dict(source_invariance=equal_archives(
            folders['native'], folders[mode], ['events', 'truth', 'segments']))
        data = archive(folders[mode]/'native_states.npz')
        masks = [np.isin(a['stage'], [1, 3, 4, 5]) for a in [base, data]]
        for key in base:
            np.testing.assert_array_equal(base[key][masks[0]], data[key][masks[1]],
                                          err_msg=mode+' seed/CKF '+key)
        checks['seed_CKF_exact_rows'] = int(masks[0].sum())
        if mode in ['updater', 'update_canonical']:
            angles = archive(folders[mode]/'native_angles.npz')
            phi, theta = angles['ePHI_flt'], angles['eTHETA_flt']
            mask = (angles['writer'] == 1) & np.isfinite(phi) & np.isfinite(theta)
            tolerance = 8*np.finfo(np.float32).eps
            bad = mask & ((np.abs(phi) > np.pi+tolerance) | (theta < -tolerance) | (theta > np.pi+tolerance))
            assert not bad.any()
            checks['finite_filtered_refit_measurements'] = int(mask.sum())
            checks['filtered_angular_chart_violations'] = int(bad.sum())
        result['comparisons'][mode] = checks
    result['definition'] = 'Actual Calypso canonical update control; exact off/smoother-only closure and seed/CKF invariance; refit forward may change'
    return result


def check_chart(args):
    folders = {mode: args.root/'smoke'/mode for mode in ['native', 'regular_chart']}
    result = dict(runs={}, comparisons={})
    source = None
    for mode, folder in folders.items():
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        run = json.loads((folder/'run.json').read_text())
        task = json.loads((folder/'source_task.json').read_text())
        identity = tuple(task[k] for k in ['input', 'skip', 'events', 'variant', 'material'])
        if source is None:
            source = identity
        assert source == identity
        assert run['success'] and run['refit_covariance_scale'] == 10
        assert run['refit_smoothing_chart_regular'] == (mode == 'regular_chart')
        assert run['refit_outcome_observer']
        assert not any(run[k] for k in ['refit_smoothing_phi_periodic',
                                        'refit_smoothing_angles_canonical',
                                        'refit_update_angles_canonical',
                                        'conversion_covariance_basis_corrected',
                                        'refit_covariance_basis_corrected',
                                        'refit_uses_input_reference_plane'])
        assert Path(run['library']).resolve() == (args.root/'software/lib').resolve()
        result['runs'][mode] = dict(source_events=task['events'], library=run['library'])
        outcomes, tracks = (archive(folder/(name+'.npz')) for name in ['refit_outcomes', 'tracks'])
        for key in ['source_entry', 'run', 'event', 'ordinal']:
            np.testing.assert_array_equal(outcomes[key], tracks[key])
    result['comparisons']['native_vs_original'] = equal_archives(folders['native'], args.original,
        ['events', 'truth', 'tracks', 'segments', 'native_states', 'native_angles', 'refit_outcomes'])
    result['comparisons']['source_invariance'] = equal_archives(folders['native'], folders['regular_chart'],
        ['events', 'truth', 'segments'])
    before, after = (archive(folders[mode]/'native_states.npz') for mode in ['native', 'regular_chart'])
    for key in ['entry', 'track_index', 'stage', 'station', 'layer']:
        np.testing.assert_array_equal(before[key], after[key])
    for label, stages in [('seed_CKF', [1, 3, 4, 5]), ('refit_forward', [8, 9])]:
        mask = np.isin(before['stage'], stages)
        for key in before:
            np.testing.assert_array_equal(before[key][mask], after[key][mask], err_msg=label+' '+key)
        result['comparisons'][label+'_exact_rows'] = int(mask.sum())
    result['definition'] = 'Actual Calypso beam-regular smoother; native exact closure and unchanged seed/CKF/refit forward states; source and outcome identities checked'
    return result


def check_chart_update(args):
    modes = ['native', 'regular_chart', 'regular_update', 'regular_both']
    folders = {mode: args.root/'smoke'/mode for mode in modes}
    result = dict(runs={}, comparisons={})
    for mode, folder in folders.items():
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        run = json.loads((folder/'run.json').read_text())
        assert run['success'] and run['refit_covariance_scale'] == 10
        assert run['refit_smoothing_chart_regular'] == (mode in ['regular_chart', 'regular_both'])
        assert run['refit_update_chart_regular'] == (mode in ['regular_update', 'regular_both'])
        assert run['refit_outcome_observer']
        assert Path(run['library']).resolve() == (args.root/'software/lib').resolve()
        actual = {str(Path(name).resolve()) for name in run['loaded_acts_libraries']}
        for name in ['libFaserActsGeometry.so', 'libFaserActsGeometryLib.so', 'libFaserActsKalmanFilter.so']:
            assert str((args.root/'software/lib'/name).resolve()) in actual
        tracks, outcomes = (archive(folder/(name+'.npz')) for name in ['tracks', 'refit_outcomes'])
        for key in ['source_entry', 'run', 'event', 'ordinal']:
            np.testing.assert_array_equal(tracks[key], outcomes[key])
        result['runs'][mode] = dict(selected=len(tracks['ordinal']),
            succeeded=int(outcomes['refit_success'].sum()), loaded_acts_libraries=sorted(actual))
    tables = ['events', 'truth', 'tracks', 'segments', 'native_states', 'native_angles', 'refit_outcomes']
    result['comparisons']['native_vs_original'] = equal_archives(folders['native'], args.original, tables)
    result['comparisons']['smoother_vs_previous'] = equal_archives(
        folders['regular_chart'], args.previous/'smoke/regular_chart', tables)
    base = archive(folders['native']/'native_states.npz')
    for mode in modes[1:]:
        checks = dict(source_invariance=equal_archives(
            folders['native'], folders[mode], ['events', 'truth', 'segments']))
        actual = archive(folders[mode]/'native_states.npz')
        for stages, label in [([1, 3, 4, 5], 'seed_CKF')]:
            a, b = [np.isin(table['stage'], stages) for table in [base, actual]]
            for key in base:
                np.testing.assert_array_equal(base[key][a], actual[key][b], err_msg=mode+' '+key)
            checks[label+'_exact_rows'] = int(a.sum())
        result['comparisons'][mode] = checks
    result['definition'] = 'Actual Calypso regular update control; exact native and smoother-only closure, source/seed/CKF unchanged, actual loaded libraries and fit return identities checked'
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--suite', choices=['angles', 'basis', 'update', 'chart', 'chart_update'], default='angles')
    parser.add_argument('--previous', type=Path)
    parser.add_argument('--original', type=Path, required=True)
    args = parser.parse_args()
    if args.suite in ['angles', 'update', 'chart_update'] and args.previous is None:
        parser.error('--previous is required for the angular controls')
    result = {'basis': check_basis, 'angles': check_angles, 'update': check_update,
              'chart': check_chart, 'chart_update': check_chart_update}[args.suite](args)
    write_json(args.root / 'smoke_validation.json', result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
