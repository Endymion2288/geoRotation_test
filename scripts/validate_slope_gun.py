#!/usr/bin/env python3
"""Validate completed full-Calypso gun pipelines and fixed-source controls."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from slope_controls import archive
from slope_tables import write_json


def validate(folder, require_library_map=False):
    assert (folder / 'exit_code.txt').read_text().strip() == '0', folder
    task = json.loads((folder / 'source_task.json').read_text())
    modes = task['reconstruction_modes']
    result = dict(source_events=task['events'], pdg=task['pdg'], seed=task['seed'],
                  source_field_comparisons={}, runs={})
    for step in ['simulate', 'digitize']:
        assert (folder / step / 'exit_code.txt').read_text().strip() == '0', step
    native = folder / modes[0]['label']
    source = {name: archive(native / (name + '.npz')) for name in ['events', 'truth', 'segments']}
    for mode in modes:
        label, variant = mode['label'], mode['variant']
        path = folder / label
        assert (folder / (label + '_wrapper') / 'exit_code.txt').read_text().strip() == '0', label
        run = json.loads((path / 'run.json').read_text())
        summary = json.loads((path / 'summary.json').read_text())
        assert run['success'] and summary['source_events'] == task['events'], label
        assert Path(run['library']).resolve() == Path(task['libraries'][variant]).resolve(), label
        assert run['args']['material'] == task['materials'][variant], label
        assert run['args']['variant'] == variant
        assert run['geometry'] == 'FASERNU-04' and run['conditions'] == 'OFLCOND-FASER-05'
        for option in ['refit-outcomes', 'correct-conversion-basis', 'correct-refit-basis',
                       'refit-input-reference-plane', 'wrap-refit-smoothing-phi',
                       'canonical-refit-smoothing-angles', 'canonical-refit-update-angles',
                       'regular-refit-smoothing-chart', 'regular-refit-update-chart']:
            assert bool(run['args'].get(option.replace('-', '_'), False)) == ('--' + option in mode.get('options', [])), (label, option)
        result['runs'][label] = dict(library=run['library'], source_events=summary['source_events'])
        if require_library_map or 'loaded_acts_libraries' in run:
            actual = {str(Path(name).resolve()) for name in run['loaded_acts_libraries']}
            expected = [Path(run['library']) / name for name in [
                'libFaserActsGeometry.so', 'libFaserActsGeometryLib.so', 'libFaserActsKalmanFilter.so']]
            assert all(path.is_file() and str(path.resolve()) in actual for path in expected), (label, actual)
            result['runs'][label]['loaded_acts_libraries'] = run['loaded_acts_libraries']
            result['runs'][label]['actual_geometry_and_tracking_libraries_verified'] = True
        # Before/after geometry can change reconstructed segments. Only the
        # controls of the same native configuration must preserve all three.
        checks = {}
        for name in ['events', 'truth', 'segments'] if variant == modes[0]['variant'] else ['events', 'truth']:
            actual = archive(path / (name + '.npz'))
            assert source[name].keys() == actual.keys(), (label, name)
            keys = list(source[name])
            if variant != modes[0]['variant']:
                # Candidate counts and propagated truth models belong to the
                # reconstruction; source identities and G4 hit truth do not.
                excluded = ['n_tracks'] if name == 'events' else ['truth_model_momentum_F_GeV', 'truth_model_valid']
                keys = [key for key in keys if key not in excluded]
            for key in keys:
                values = source[name][key]
                equal = np.array_equal(values, actual[key]) if values.dtype.kind in 'US' else np.array_equal(values, actual[key], equal_nan=True)
                assert equal, (label, name, key)
            checks[name] = keys
        result['source_field_comparisons'][label] = checks
        if '--refit-outcomes' in mode.get('options', []):
            outcomes, tracks = (archive(path / (name + '.npz')) for name in ['refit_outcomes', 'tracks'])
            for key in ['source_entry', 'run', 'event', 'ordinal']:
                np.testing.assert_array_equal(outcomes[key], tracks[key], err_msg=label + ' ' + key)
            assert np.isin(outcomes['refit_success'], [0, 1]).all()
            result['runs'][label]['refit_outcomes'] = summary['refit_outcomes']
    result['all_steps_exit_zero'] = True
    result['definition'] = 'Same full Calypso G4/digitization RDO; actual tracking controls and observer identities; no tracking recomputation'
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--require-library-map', action='store_true')
    args = parser.parse_args()
    result = validate(args.folder, args.require_library_map)
    write_json(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
