#!/usr/bin/env python3
"""Read native Calypso fits with frozen seed, admitted cluster IDs and process terms."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

import numpy as np

from charge_measurement_coordinates import coordinates, outcome
from charge_replay_analysis import operation, read_audit
from material_encoding_replay import target_path
from slope_tables import write_json

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_charge_process_noise import COORDINATES, PROCESSES


def native_inputs(root):
    rows = json.loads((root/'ckf_measurement_v1/analysis_v3/full/events.json').read_text())
    old = {(r['dsid'], r['sample_entry']): r for r in rows}
    sdo = json.loads((root/'sdo_strip_coordinate_v1/analysis_v1/full/events.json').read_text())
    for row in sdo:
        old[row['dsid'], row['sample_entry']]['modes']['sdo_strips_fixed'] = row['modes']['sdo_strips_fixed']
    return old


def coordinate_checks(native, current, mode):
    assert native.keys() == current.keys()
    for hit, old in native.items():
        row = current[hit]
        for key in ['native_local_mm', 'lorentz_shift_mm', 'phi_pitch_mm', 'width_mm',
                    'width_cells', 'error_strategy', 'majority_barcode', 'channels',
                    'same_wafer_hits', 'used_covariance_mm2']:
            assert row[key] == old[key], (hit, key)
        expected = old['used_local_mm'][0]
        if mode == 'sdo_strips_fixed' and row['sdo_supported_channels']:
            expected = row['sdo_supported_local_x_mm']
        if mode == 'centroid_fixed' and row['majority_centroid_local_mm'] is not None:
            expected = row['majority_centroid_local_mm'][0]+row['lorentz_shift_mm']
        assert row['used_local_mm'] == [expected, old['used_local_mm'][1]]


def resource(folder):
    row = {}
    for line in (folder/'tracking_resource.txt').read_text().splitlines():
        for name, key in [('User time (seconds)', 'user_s'), ('System time (seconds)', 'system_s'),
                          ('Maximum resident set size (kbytes)', 'rss_KiB')]:
            if line.strip().startswith(name+':'):
                row[key] = float(line.split(':', 1)[1])
    assert len(row) == 3
    return row


def reduce(root, output, smoke=False, jobs_path=None):
    control = root/'fixed_hit_process_v1'
    assert json.loads((control/'smoke_validation.json').read_text())['all_exact']
    previous = native_inputs(root)
    jobs = json.loads((jobs_path or control/('fixed_smoke_tasks.json' if smoke else 'all_tasks.json')).read_text())
    assert len(jobs) == (48 if smoke else 792)
    expected_libs = {str((control/'software/lib'/name).resolve()) for name in
        ['libFaserActsGeometry.so', 'libFaserActsGeometryLib.so', 'libFaserActsKalmanFilter.so']}
    expected_libs.add(str((root/'observer_hitcontrol_v2/lib/libCalypsoSlopeDiagnostics.so').resolve()))
    expected_clusters = {str((root/'sdo_strip_coordinate_v1/software/lib'/name).resolve()) for name in
        ['libFaserSiClusterizationTool.so', 'libFaserSiClusterizationToolLib.so']}
    groups, events, updates, resources = defaultdict(Counter), {}, 0, []
    for job in jobs:
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0', str(folder)
        task = json.loads((folder/'source_task.json').read_text())
        run = json.loads((folder/'run.json').read_text())
        coordinate, process = task['coordinate_mode'], task['process_mode']
        assert run['success'] and run['isolate_smoothing'] and run['run_only_input_seeds']
        assert run['ckf_material_scattering'] == (process not in ['no_scattering', 'no_material'])
        assert run['ckf_material_energy_loss'] == (process not in ['no_energy_loss', 'no_material'])
        assert set(run['loaded_acts_libraries']) == expected_libs
        assert set(run['loaded_cluster_libraries']) == expected_clusters
        key = task['dsid'], task['sample_entry']
        old = previous[key]
        original_folder = root/'ckf_measurement_v1/reconstruction'/task['dsid']/'native'/folder.name
        original = read_audit(original_folder/'charge_audit.jsonl.gz')
        original_hits = old['native_ordered_measurement_hits']
        original_traces = [r for r in original if r['stage'] == 'selected_state_trace' and
            [s['hit'] for s in r['steps'] if s['measurement_flag']] == original_hits]
        assert len(original_traces) == 1
        seed_index = original_traces[0]['seed']
        seed = next(r for r in original if r['stage'] == 'seed' and r['seed'] == seed_index)
        audit = read_audit(folder/'charge_audit.jsonl.gz')
        inputs = [r for r in audit if r['stage'] == 'seed']
        assert len(inputs) == 1 and inputs[0]['seed'] == seed_index
        for field in ['bound', 'cov', 'surface']:
            assert inputs[0][field] == seed[field], (key, field)
        assert sorted(inputs[0]['hits']) == sorted(seed['hits'])
        source = [r for r in audit if r['stage'] == 'source_link_selection']
        assert len(source) == 1 and source[0]['hit_ids'] == sorted(original_hits)
        assert source[0]['retained'] == len(original_hits)
        native_coords, actual_coords = coordinates(original_folder), coordinates(folder)
        coordinate_checks(native_coords, actual_coords, coordinate)
        primary = old['modes']['native']['primary_barcode']
        hit_truth = {hit: sorted({barcode for channel in row['channels'] for barcode, _ in channel['deposits']})
                     for hit, row in actual_coords.items()}
        charge = 1 if task['dsid'] == '100116' else -1
        for row in audit:
            if row['stage'] != 'selected_state_trace':
                continue
            assert row['seed'] == seed_index
            for step in row['steps']:
                if not step['measurement_flag']:
                    continue
                assert step['hit'] in original_hits
                np.testing.assert_array_equal(step['measurement'], [actual_coords[step['hit']]['used_local_mm'][0]])
                np.testing.assert_array_equal(step['measurement_cov'], [[.0231**2]])
                operation(step, hit_truth[step['hit']], primary, charge)
                updates += 1
        result = outcome(folder, task)
        path = target_path(folder, task, hit_truth)[0] if result['final_qop_GeVinv'] is not None else None
        ordered = path['ordered_measurement_hits'] if path else []
        assert set(ordered) <= set(original_hits)
        result.update(seed=seed_index, same_ordered_hits=ordered == original_hits,
            selected_measurement_hits=ordered, missing_native_hit_count=len(set(original_hits)-set(ordered)),
            native_reference_qop_GeVinv=path['native_reference_qop_GeVinv'] if path else None,
            last_update_into_wrong=path['last_update_into_wrong'] if path else None,
            final_wrong_entry=path['final_wrong_entry'] if path else None,
            audit=str(folder/'charge_audit.jsonl.gz'))
        if path and old['cohort'] == 'pure_wrong':
            trace = next(r for r in audit if r['stage'] == 'selected_state_trace' and
                         r['seed'] == path['seed'] and r['candidate'] == path['candidate'])
            result['per_update'] = path['per_update']
            result['qop_trace'] = [dict(hit=step.get('hit'), z_F_mm=step['z_F_mm'],
                measurement=step['measurement_flag'], **{name:step[name][4] for name in ['predicted', 'filtered', 'smoothed']},
                **{name+'_variance':step[name+'_cov'][4][4] for name in ['predicted', 'filtered', 'smoothed']})
                for step in trace['steps'] if step['measurement_flag']]
        row = events.setdefault(key, dict(dsid=task['dsid'], sample_entry=task['sample_entry'],
            input=task['input'], source_entry=task['skip'], cohort=old['cohort'],
            original_selected_seed=seed_index, original_ordered_hits=original_hits,
            unrestricted_modes={mode:old['modes'][mode] for mode in COORDINATES}, modes={}))
        mode = coordinate+'|'+process
        assert mode not in row['modes']
        row['modes'][mode] = result
        counts = groups['|'.join([task['dsid'], old['cohort'], coordinate, process])]
        counts['events'] += 1
        counts[result['status']] += 1
        counts['same_ordered_hits'] += result['same_ordered_hits']
        counts['final_correct_sign'] += result['final_qop_GeVinv'] is not None and charge*result['final_qop_GeVinv'] > 0
        counts['ckf_correct_sign'] += result['native_reference_qop_GeVinv'] is not None and charge*result['native_reference_qop_GeVinv'] > 0
        resources.append(dict(dsid=task['dsid'], sample_entry=task['sample_entry'], mode=mode, **resource(folder)))
        print(*key, mode, result['status'], result['same_ordered_hits'], flush=True)
    for row in events.values():
        assert set(row['modes']) == {c+'|'+p for c in COORDINATES for p in PROCESSES}
    assert len(events) == (4 if smoke else 66)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/'events.json', list(events.values()))
    write_json(output/'resources.json', dict(rows=resources, jobs=len(resources),
        total_cpu_s=sum(r['user_s']+r['system_s'] for r in resources),
        maximum_rss_KiB=max(r['rss_KiB'] for r in resources)))
    write_json(output/'summary.json', dict(events=len(events), jobs=len(jobs),
        groups={key:dict(value) for key,value in groups.items()},
        actual_measurement_updates_checked=updates, full_seed_and_surface_exact=True,
        original_ancestor_seed_fixed=True, all_admitted_hit_ids_exact=True,
        all_selected_hits_are_subset=True, all_actual_measurement_values_and_variances_verified=True,
        scope='Actual geoRotation Calypso CKF/refit with isolated smoothing ON and corrected material encoding; six original pure wrong + 60 normal controls',
        interpretation='One original ancestor seed and exact admitted cluster IDs. A native finder may still reject hits or fail; ordered path equality and missing hits are counted. Process switches affect only CKF material terms; refit and observer stay native. Outcome-selected mechanism cohort, not a production-rate estimate.'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--jobs', type=Path,
                        help='Frozen resolved attempts; retain the same event and mode matrix')
    args = parser.parse_args()
    reduce(args.root, args.output, args.smoke, args.jobs)


if __name__ == '__main__':
    main()
