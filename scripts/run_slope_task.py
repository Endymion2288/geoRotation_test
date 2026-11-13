#!/usr/bin/env python3
"""Run native Calypso on scratch and retain compact angle analysis tables."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import awkward as ak
import numpy as np
import uproot


def read_ntuple(output):
    with uproot.open('native_ntuple.root') as file:
        key = next(k for k, v in file.classnames().items() if v == 'TTree')
        tree = file[key]
        branches = ['Track_px0', 'Track_py0', 'Track_pz0', 'Track_ThetaY_atVetoNu',
                    'Track_ThetaX_atVetoNu', 'Track_X_atVetoNu', 'Track_Y_atVetoNu']
        native = tree.arrays(branches, library='ak')
    reference = {k: ak.to_numpy(ak.flatten(native[k])) for k in branches}
    reference['track_offsets'] = np.r_[0,np.cumsum(ak.to_numpy(ak.num(native[branches[0]])))]
    np.savez_compressed('native_ntuple_branches.npz', **reference)
    shutil.copy2('native_ntuple_branches.npz',output/'native_ntuple_branches.npz')
    return reference


def check_ntuple(arrays, output, reference, physical_angles=False, rotated=False):
    tracks = arrays['tracks']
    long = tracks['station_bits'] & 14 == 14
    raw = tracks['upstream_raw'][long]
    assert len(raw) == len(reference['Track_py0'])
    max_momentum_error = 0.
    for i, name in enumerate(['Track_px0', 'Track_py0', 'Track_pz0']):
        delta = np.abs(reference[name]*.001-raw[:, i+3])
        max_momentum_error = max(max_momentum_error, float(delta.max(initial=0)))
        assert np.allclose(reference[name]*.001, raw[:, i+3], rtol=2e-14, atol=1e-13)
    projection = tracks['projection'][long]
    valid = tracks['projection_valid'][long, 1] == 1
    expected = projection[:, 1, 7]
    if physical_angles:
        expected = np.arctan2(projection[:,1,4],projection[:,1,5])
    observed = reference['Track_ThetaY_atVetoNu']
    assert np.all(observed[~valid] == 999999.)
    delta = np.abs(observed[valid]-expected[valid])
    assert np.allclose(observed[valid], expected[valid], rtol=2e-10, atol=2e-12)
    momentum=projection[:,1,3:6]
    with np.errstate(divide='ignore',invalid='ignore'):
        expected_x = np.arctan2(momentum[:,0],momentum[:,2]) if physical_angles else np.arctan(-momentum[:,2]/momentum[:,0]) if rotated else np.arctan(momentum[:,0]/momentum[:,2])
    assert np.allclose(reference['Track_ThetaX_atVetoNu'][valid],expected_x[valid],rtol=2e-10,atol=2e-12)
    result = dict(native_long_tracks=int(long.sum()),
                  max_native_upstream_momentum_error_GeV=max_momentum_error,
                  common_vetonu_success=int(valid.sum()),
                  native_and_observer_failures_equal=True,
                  physical_angle_control=physical_angles,
                  max_vetonu_legacy_angle_error_rad=float(delta.max(initial=0)))
    (output/'ntuple_validation.json').write_text(json.dumps(result, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    task = json.loads(args.task.read_text())
    output = Path(task['output'])
    output.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, task['runner'], '--input', task['input'], '--events', str(task['events']),
               '--skip', str(task['skip']), '--material', task['material'], '--variant', task['variant'],
               *task.get('options', [])]
    begin = time.monotonic()
    with Path('tracking.log').open('w') as log:
        status = subprocess.run(['/usr/bin/time', '-v', '-o', 'tracking_resource.txt', *command],
                                stdout=log, stderr=subprocess.STDOUT).returncode
    for name in ['run.json', 'tracking_resource.txt']:
        if Path(name).exists():
            shutil.copy2(name, output/name)
    with Path('tracking.log').open('rb') as source, gzip.open(output/'tracking.log.gz', 'wb') as target:
        shutil.copyfileobj(source, target)
    (output/'source_task.json').write_text(json.dumps(task, indent=2)+'\n')
    if status:
        raise RuntimeError(f'Calypso reconstruction failed: exit {status}; see {output}/tracking.log.gz')
    reference = read_ntuple(output) if '--ntuple-check' in task.get('options', []) else None
    if '--no-observer' in task.get('options', []):
        assert reference is not None
        (output/'summary.json').write_text(json.dumps(dict(source_events=len(reference['track_offsets'])-1,
            native_long_tracks=len(reference['Track_py0']), native_only=True), indent=2)+'\n')
        return
    arrays = {}
    with uproot.open('slopes.root') as file:
        for name in ['events', 'truth', 'tracks', 'segments']:
            values = file[name].arrays(library='ak')
            arrays[name] = {key: ak.to_numpy(values[key]) for key in values.fields if key != 'hit_identifiers'}
            if name == 'tracks':
                lengths = ak.to_numpy(ak.num(values['hit_identifiers']))
                arrays[name]['hit_offsets'] = np.r_[0, np.cumsum(lengths)]
                arrays[name]['hit_identifiers'] = ak.to_numpy(ak.flatten(values['hit_identifiers']))
    assert np.array_equal(arrays['events']['source_entry'], np.arange(task['skip'], task['skip']+task['events']))
    summary = dict(task=task, source_events=len(arrays['events']['source_entry']),
                   primary_truth=len(arrays['truth']['source_entry']),
                   native_tracks=len(arrays['tracks']['source_entry']),
                   segments=len(arrays['segments']['source_entry']),
                   wall_s=time.monotonic()-begin, units='mm, GeV, radians; segments.parameters_F momentum in MeV',
                   planes=['fixed st1 z_F=0 mm', 'fixed VetoNu z_F=-3112 mm', 'same primary SiHit midpoint z_F'],
                   native_tracking=True, standalone_tracking=False)
    if '--ntuple-check' in task.get('options', []):
        check_ntuple(arrays, output, reference, '--physical-ntuple-angles' in task.get('options', []),task['variant']!='beforeGeo')
    for name, data in arrays.items():
        # Truth is retained once per primary; secondary-majority tracks retain their own reference.
        if name == 'tracks':
            for branch in ['production_position_F_mm', 'truth_hit_counts', 'truth_station_position_F_mm',
                           'truth_station_direction_F', 'truth_sensor_position_F_mm', 'truth_sensor_direction_F']:
                data.pop(branch)
        np.savez_compressed(name+'.npz', **data)
        shutil.copy2(name+'.npz',output/(name+'.npz'))
    if '--refit-outcomes' in task.get('options', []):
        outcomes = np.atleast_1d(np.genfromtxt('refit_outcomes.csv', delimiter=',', names=True, dtype=np.int64))
        fields = {key: outcomes[key] for key in outcomes.dtype.names}
        fields['source_entry'] = fields['local_event_index'] + task['skip']
        for key in ['source_entry', 'run', 'event', 'ordinal']:
            np.testing.assert_array_equal(fields[key], arrays['tracks'][key], err_msg='refit outcome identity: '+key)
        assert np.isin(fields['refit_success'], [0, 1]).all()
        np.savez_compressed('refit_outcomes.npz', **fields)
        shutil.copy2('refit_outcomes.npz', output/'refit_outcomes.npz')
        summary['refit_outcomes'] = dict(selected=len(outcomes),
            succeeded=int(fields['refit_success'].sum()),
            fallback=int((fields['refit_success'] == 0).sum()),
            identity='Checked against observer source entry, run, event and final ordinal')
    if '--stage-audit' in task.get('options', []):
        from slope_stage_extract import extract
        summary['audit_state_rows'] = extract('stages.jsonl', output, task)
        Path('stages.jsonl').unlink()
    if '--native-states' in task.get('options', []):
        from slope_native_states import extract
        summary['native_state_rows'] = extract(output, task, arrays)
    summary['compact_bytes'] = sum(p.stat().st_size for p in output.glob('*.npz'))
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
