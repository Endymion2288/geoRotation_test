#!/usr/bin/env python3
"""Check the complete generated charge-control chain and loaded Calypso libs."""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    output = args.root/'symmetric_gun_smoke_v1'
    tasks = json.loads((output/'manifest.json').read_text())['source_tasks']
    generated = []
    samples = {}
    for task in tasks:
        folder = Path(task['output'])
        for name in ['simulate', 'digitize', 'beforeGeo_wrapper', 'beforeGeo_regular_both_wrapper', 'afterGeo_wrapper']:
            assert (folder/name/'exit_code.txt').read_text().strip() == '0'
        values = []
        for mode in ['beforeGeo', 'beforeGeo_regular_both', 'afterGeo']:
            current = folder/mode
            run = json.loads((current/'run.json').read_text())
            assert run['isolate_smoothing']
            software = 'software_afterGeo' if mode == 'afterGeo' else 'software_refit_before'
            expected = {str((args.root/software/'lib'/name).resolve()) for name in [
                'libFaserActsGeometry.so', 'libFaserActsGeometryLib.so', 'libFaserActsKalmanFilter.so']}
            expected.add(str((args.root/'observer/lib/libCalypsoSlopeDiagnostics.so').resolve()))
            assert set(run['loaded_acts_libraries']) == expected
            with np.load(current/'truth.npz') as truth:
                mask = (truth['pdg'] == task['pdg']) & (truth['status'] == 1) & (truth['barcode'] > 0) & (truth['barcode'] < 200000)
                assert mask.sum() == task['events']
                values.append({field: truth[field][mask] for field in [
                    'production_position_F_mm', 'production_momentum_F_GeV', 'truth_hit_counts']})
        for row in values[1:]:
            for key in row:
                np.testing.assert_array_equal(row[key], values[0][key])
        generated.append(values[0])
        samples[task['dsid']] = dict(events=task['events'], all_steps_succeeded=True,
            actual_libraries_match=True, all_versions_same_truth=True)
    for key in ['production_position_F_mm', 'production_momentum_F_GeV']:
        np.testing.assert_array_equal(generated[0][key], generated[1][key], err_msg=key)
    result = dict(samples=samples, all_steps_succeeded=True, actual_libraries_match=True,
        identical_generated_conditions=True, all_versions_same_truth=True,
        scope='Full Calypso G4 -> digitization -> before/native, before/chart, after/native; 20 events per charge; isolation ON')
    (output/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
