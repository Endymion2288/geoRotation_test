#!/usr/bin/env python3
"""Check actual default and nominal Calypso surfaces and reconstruction tables."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def sensor_transforms(path):
    rows = list(csv.DictReader(path.open()))
    ids = [int(row['wafer_id']) for row in rows]
    assert len(set(ids)) == len(ids) and len(ids) > 100
    matrices = {label: np.asarray([[[float(row[f'{label}{i}{j}']) for j in range(4)]
                                  for i in range(4)] for row in rows])
                for label in ['actual', 'nominal']}
    delta = matrices['actual']-matrices['nominal']
    summary = dict(sensors=len(ids), frame='Physical FASER frame; local-to-global sensor transforms',
        maximum_translation_difference_mm=float(np.linalg.norm(delta[:, :3, 3], axis=1).max()),
        maximum_rotation_matrix_difference=float(np.abs(delta[:, :3, :3]).max()),
        maximum_absolute_matrix_difference=float(np.abs(delta).max()))
    return ids, matrices, summary


def equal_tables(current, original):
    checked = []
    for name in ['events', 'truth', 'tracks', 'segments']:
        with np.load(current/(name+'.npz')) as a, np.load(original/(name+'.npz')) as b:
            assert set(a.files) == set(b.files)
            for field in a.files:
                np.testing.assert_array_equal(a[field], b[field], err_msg=str(current)+'/'+name+'/'+field)
                checked.append(name+'/'+field)
    assert len(checked) == 58
    return checked


def check_run(folder, root, nominal):
    run = json.loads((folder/'run.json').read_text())
    assert run['success'] and run['isolate_smoothing']
    assert run['alignment']['nominal_alignment'] == nominal
    assert run['alignment']['conditions_displacements_disabled'] == nominal
    expected = str((root/'software_nominal_alignment_v1/lib/libFaserSCT_ConditionsAlgorithms.so').resolve())
    assert run['loaded_conditions_libraries'] == [expected]
    assert bool(run['alignment']['tracker_alignment_folders']) != nominal
    ids, matrices, summary = sensor_transforms(folder/'alignment.csv')
    if nominal:
        assert summary['maximum_absolute_matrix_difference'] < 1e-9
    return ids, matrices, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['compatibility', 'gun-smoke'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output = root/('alignment_compatibility_v1' if args.phase == 'compatibility' else 'alignment_gun_smoke_v1')
    tasks = json.loads((output/'manifest.json').read_text())['source_tasks']
    recovery = root/'alignment_gun_recovery_v2/manifest.json'
    if args.phase == 'gun-smoke' and recovery.exists():
        replacements = {(row['dsid'], row['nominal_alignment']): row for row in
                        json.loads(recovery.read_text())['source_tasks']}
        tasks = [replacements.get((row['dsid'], row['nominal_alignment']), row) for row in tasks]
    samples = {}
    if args.phase == 'compatibility':
        reference = None
        for task in tasks:
            folder = Path(task['output'])
            assert (folder/'exit_code.txt').read_text().strip() == '0'
            checked = equal_tables(folder, Path(task['original_output']))
            ids, matrices, summary = check_run(folder, root, False)
            if reference:
                assert ids == reference[0]
                for key in matrices:
                    np.testing.assert_array_equal(matrices[key], reference[1][key])
            else:
                reference = ids, matrices
            samples[task['variant']+'|'+task['dsid']] = dict(fields_exact=len(checked), surfaces=summary)
        result = dict(samples=samples, all_58_fields_exact=True, actual_libraries_match=True,
            all_default_sensor_transforms_equal=True, controls=len(tasks))
    else:
        generated = {}
        reconstructed = {}
        for task in tasks:
            folder = Path(task['output'])
            nominal = task['nominal_alignment']
            audit_reference = None
            for name in ['simulate', 'digitize']:
                assert (folder/name/'exit_code.txt').read_text().strip() == '0'
                run = json.loads((folder/(name+'_run.json')).read_text())
                assert run['success'] and run['alignment']['nominal_alignment'] == nominal
                assert run['alignment']['conditions_displacements_disabled'] == nominal
                assert bool(run['alignment']['tracker_alignment_folders']) != nominal
                if nominal:
                    ids, matrices, summary = sensor_transforms(folder/name/'alignment.csv')
                    assert summary['maximum_absolute_matrix_difference'] < 1e-9
                    if audit_reference:
                        assert ids == audit_reference[0]
                        np.testing.assert_allclose(matrices['actual'], audit_reference[1]['actual'], rtol=0, atol=1e-9)
                    else:
                        audit_reference = ids, matrices
            sample = {}
            truth_reference = None
            for mode in task['reconstruction_modes']:
                label = mode['label']
                assert (folder/(label+'_wrapper')/'exit_code.txt').read_text().strip() == '0'
                ids, matrices, summary = check_run(folder/label, root, nominal)
                if nominal:
                    assert ids == audit_reference[0]
                    np.testing.assert_allclose(matrices['actual'], audit_reference[1]['actual'], rtol=0, atol=1e-9)
                else:
                    equal_tables(folder/label, Path(task['original_output'])/label)
                with np.load(folder/label/'truth.npz') as t:
                    primary = (t['pdg'] == task['pdg']) & (t['status'] == 1) & (t['barcode'] > 0) & (t['barcode'] < 200000)
                    assert int(primary.sum()) == task['events']
                    truth = {key: t[key][primary] for key in ['production_position_F_mm', 'production_momentum_F_GeV', 'truth_hit_counts']}
                if truth_reference:
                    for key in truth:
                        np.testing.assert_array_equal(truth[key], truth_reference[key])
                else:
                    truth_reference = truth
                sample[label] = summary
            generated[task['dsid'], nominal] = truth_reference
            reconstructed[task['dsid'], nominal] = folder
            samples[task['dsid']+('|nominal' if nominal else '|default')] = sample
        default_nominal_tables = {}
        for dsid in ['100116', '100117']:
            for key in ['production_position_F_mm', 'production_momentum_F_GeV']:
                np.testing.assert_array_equal(generated[dsid, False][key], generated[dsid, True][key])
            for label in samples[dsid+'|default']:
                fields = equal_tables(reconstructed[dsid, True]/label,
                                      reconstructed[dsid, False]/label)
                default_nominal_tables[dsid+'|'+label] = dict(fields_exact=len(fields))
        for nominal in [False, True]:
            for key in ['production_position_F_mm', 'production_momentum_F_GeV']:
                np.testing.assert_array_equal(generated['100116', nominal][key], generated['100117', nominal][key])
        result = dict(samples=samples, all_steps_succeeded=True, actual_libraries_match=True,
            all_default_58_fields_exact=True, nominal_G4_digitization_reconstruction_surfaces_match=True,
            identical_generated_conditions=True, all_reconstruction_modes_same_truth=True,
            default_nominal_tables=default_nominal_tables,
            all_default_nominal_58_fields_exact=True,
            source_events_per_charge_alignment=tasks[0]['events'])
    (output/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
