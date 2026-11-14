#!/usr/bin/env python3
"""Compare actual Calypso SDO-supported strip positions with centroid controls."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

import numpy as np

from charge_measurement_coordinates import coordinates, coordinate_features, full_seed_checks, outcome
from charge_measurement_details import quantiles
from charge_replay_analysis import operation, read_audit
from material_encoding_replay import target_path
from slope_tables import write_json

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_sdo_strip_coordinates import coordinate_check


def unlinked_features(row):
    hits = [hit for hit in row['same_wafer_hits'] if hit['barcode'] == 0]
    channels = [channel for channel in row['channels'] if not channel['deposits']]
    half_pitch = .5*row['phi_pitch_mm']
    overlaps = []
    for channel in channels:
        centre = channel['centre_local_mm'][0]
        if any(max(hit['start_local_mm'][0], hit['end_local_mm'][0]) >= centre-half_pitch and
               min(hit['start_local_mm'][0], hit['end_local_mm'][0]) <= centre+half_pitch for hit in hits):
            overlaps.append(channel['strip'])
    result = dict(barcode_zero_sihits=len(hits),
        barcode_zero_energy_deposit_MeV=sum(hit['energy_loss_MeV'] for hit in hits),
        no_sdo_strips_with_barcode_zero_endpoint_x_overlap=overlaps)
    if hits:
        points = np.array([point for hit in hits for point in [hit['start_local_mm'], hit['end_local_mm']]])
        result.update(barcode_zero_local_min_mm=points.min(axis=0).tolist(),
                      barcode_zero_local_max_mm=points.max(axis=0).tolist())
    return result


def reduce(root, output, smoke=False):
    control = root/'sdo_strip_coordinate_v1'
    old_control = root/'ckf_measurement_v1'
    assert json.loads((control/'smoke_validation.json').read_text())['all_exact']
    source = json.loads((old_control/'analysis_v3/full/events.json').read_text())
    previous = {(r['dsid'], r['sample_entry']): r for r in source}
    jobs = json.loads((control/('smoke_tasks.json' if smoke else 'all_tasks.json')).read_text())
    jobs = [job for job in jobs if '/sdo_strips_fixed/' in job['output']]
    assert len(jobs) == (4 if smoke else 66)
    events, spatial, evidence, counts = [], [], [], defaultdict(Counter)
    source_inputs = set()
    total_updates = 0
    expected_acts = {str((old_control/'software/lib'/name).resolve()) for name in
        ['libFaserActsGeometry.so', 'libFaserActsGeometryLib.so', 'libFaserActsKalmanFilter.so']}
    expected_acts.add(str((root/'observer_hitcontrol_v2/lib/libCalypsoSlopeDiagnostics.so').resolve()))
    expected_cluster = {str((control/'software/lib'/name).resolve()) for name in
        ['libFaserSiClusterizationTool.so', 'libFaserSiClusterizationToolLib.so']}
    for job in jobs:
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0', str(folder)
        task = json.loads((folder/'source_task.json').read_text())
        source_inputs.add(task['input'])
        run = json.loads((folder/'run.json').read_text())
        assert run['success'] and run['isolate_smoothing'] and run['args']['sdo_supported_strip_position']
        assert set(run['loaded_acts_libraries']) == expected_acts
        assert set(run['loaded_cluster_libraries']) == expected_cluster
        dsid, entry = task['dsid'], task['sample_entry']
        old = previous[dsid, entry]
        native_folder = old_control/'reconstruction'/dsid/'native'/f'e{entry:08d}'
        original = read_audit(native_folder/'charge_audit.jsonl.gz')
        audit = read_audit(folder/'charge_audit.jsonl.gz')
        seed = full_seed_checks(original, audit)
        native_coords, new_coords = coordinates(native_folder), coordinates(folder)
        coordinate_validation = coordinate_check(native_coords, new_coords, True)
        hit_truth = {hit: sorted({b for channel in row['channels'] for b, _ in channel['deposits']})
                     for hit, row in new_coords.items()}
        stats = outcome(folder, task)
        target = target_path(folder, task, hit_truth)[0] if stats['final_qop_GeVinv'] is not None else None
        primary = old['modes']['native']['primary_barcode']
        hits = target['ordered_measurement_hits'] if target else []
        modes = {name: old['modes'][name] for name in ['native', 'centroid_fixed']}
        modes['sdo_strips_fixed'] = dict(stats, same_ordered_hits=hits == old['native_ordered_measurement_hits'],
            ordered_hits=hits, seed_checks=seed,
            native_reference_qop_GeVinv=target['native_reference_qop_GeVinv'] if target else None,
            last_update_into_wrong=target['last_update_into_wrong'] if target else None,
            final_wrong_entry=target['final_wrong_entry'] if target else None,
            per_update=target['per_update'] if target else [],
            audit=str(folder/'charge_audit.jsonl.gz'))
        for row in audit:
            if row['stage'] != 'selected_state_trace':
                continue
            for step in row['steps']:
                if not step['measurement_flag']:
                    continue
                hit = new_coords[step['hit']]
                np.testing.assert_array_equal(step['measurement'], [hit['used_local_mm'][0]])
                np.testing.assert_array_equal(step['measurement_cov'], [[.0231**2]])
                operation(step, hit_truth[step['hit']], primary, 1 if dsid == '100116' else -1)
                total_updates += 1
        for hit in old['native_ordered_measurement_hits']:
            original_row, new_row = native_coords[hit], new_coords[hit]
            row = coordinate_features(original_row, primary)
            row.update(unlinked_features(original_row))
            row.update(dsid=dsid, sample_entry=entry, cohort=old['cohort'],
                sdo_supported_channels=new_row['sdo_supported_channels'],
                sdo_used_x_mm=new_row['used_local_mm'][0],
                native_minus_sdo_mm=original_row['used_local_mm'][0]-new_row['used_local_mm'][0],
                selected_in_sdo=hit in hits)
            if 'primary_centroid_local_mm' in row:
                row['sdo_minus_primary_centroid_mm'] = new_row['used_local_mm'][0]-row['primary_centroid_local_mm'][0]
            spatial.append(row)
        if old['cohort'] == 'pure_wrong':
            rows = [row for row in spatial if row['dsid'] == dsid and row['sample_entry'] == entry
                    and 'measured_minus_primary_centroid_mm' in row]
            largest = max(rows, key=lambda row: abs(row['measured_minus_primary_centroid_mm']))
            hit = largest['hit_id']
            evidence.append(dict(dsid=dsid, sample_entry=entry, features=largest,
                native_cluster=native_coords[hit], supported_used_x_mm=new_coords[hit]['used_local_mm'][0],
                native_last_crossing=old['modes']['native']['last_update_into_wrong']))
        for mode, value in modes.items():
            group = counts['|'.join([dsid, old['cohort'], mode])]
            group['events'] += 1
            group[value['status']] += 1
            group['same_ordered_hits'] += value['same_ordered_hits']
            group['final_correct_sign'] += value['final_qop_GeVinv'] is not None and value['final_qop_GeVinv']*(1 if dsid == '100116' else -1) > 0
        events.append(dict(dsid=dsid, sample_entry=entry, input=task['input'], source_entry=task['skip'],
            cohort=old['cohort'], modes=modes,
            coordinates=coordinate_validation, native_ordered_hits=old['native_ordered_measurement_hits']))
        print(dsid, entry, {name: value['status'] for name, value in modes.items()}, flush=True)
    spatial_summary = {}
    for cohort in ['pure_wrong', 'normal']:
        rows = [row for row in spatial if row['cohort'] == cohort]
        valid = [row for row in rows if 'sdo_minus_primary_centroid_mm' in row]
        spatial_summary[cohort] = dict(native_selected_states=len(rows),
            coordinates_changed=sum(row['native_minus_sdo_mm'] != 0 for row in rows),
            coordinates_changed_above_1e_minus12_mm=sum(abs(row['native_minus_sdo_mm']) > 1e-12 for row in rows),
            fallback_no_majority_support=sum(row['sdo_supported_channels'] == 0 for row in rows),
            native_absolute_error_mm=quantiles([abs(row['measured_minus_primary_centroid_mm']) for row in valid]),
            supported_absolute_error_mm=quantiles([abs(row['sdo_minus_primary_centroid_mm']) for row in valid]),
            states_with_barcode_zero_sihits=sum(row['barcode_zero_sihits'] > 0 for row in rows),
            no_sdo_strips_with_barcode_zero_endpoint_x_overlap=sum(len(row['no_sdo_strips_with_barcode_zero_endpoint_x_overlap']) for row in rows))
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/'events.json', events)
    write_json(output/'coordinates.json', spatial)
    write_json(output/'cluster_evidence.json', evidence)
    write_json(output/'summary.json', dict(events=len(events), source_files=len(source_inputs),
        source_inputs=sorted(source_inputs), actual_measurement_updates_checked=total_updates,
        groups={key: dict(value) for key, value in counts.items()}, spatial=spatial_summary,
        all_full_seeds_exact=True, all_native_widths_errors_channel_lists_sihits_exact=True,
        all_actual_ckf_coordinates_and_fixed_variances_verified=True,
        scope='Same six pure wrong and 60 normal controls; native Calypso strip centres, corrected material map, isolated smoothing ON',
        interpretation='SDO controls remove unsupported-strip influence on measured X while retaining strip quantization and the original cluster identifiers/width/error. Comparisons use each original native selected measurement as the spatial denominator; selected-hit migration is reported separately.'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()
    reduce(args.root, args.output, args.smoke)


if __name__ == '__main__':
    main()
