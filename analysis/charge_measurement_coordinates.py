#!/usr/bin/env python3
"""Inspect real strip geometry and actual CKF measurement-error controls."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import sys

import numpy as np

from charge_centroid_analysis import read_tables
from charge_fixed_seed_analysis import quality_record, failed_quality
from charge_replay_analysis import operation, read_audit
from charge_tables import STATUS, derive
from material_encoding_replay import target_path
from slope_tables import write_json
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from cluster_coordinate_io import read_coordinates


def coordinates(folder):
    result = {}
    for row in read_coordinates(folder/'cluster_coordinates.jsonl.gz'):
        key = int(row['cluster_id'])
        assert key not in result
        result[key] = row
    return result


def coordinate_features(row, primary):
    """Use native sensor-local endpoints and native cell lookup, without transport."""
    hits = [hit for hit in row['same_wafer_hits'] if hit['barcode'] == primary]
    result = dict(hit_id=row['cluster_id'], station=row['station'], layer=row['layer'],
        strip_side=row['side'], n_strips=len(row['channels']), phi_pitch_mm=row['phi_pitch_mm'],
        cluster_sigma_mm=float(np.sqrt(row['used_covariance_mm2'][0])),
        native_measured_x_mm=row['native_local_mm'][0],
        majority_barcode=row['majority_barcode'], primary_sihit_steps=len(hits),
        native_channel_strip_indices=[channel['strip'] for channel in row['channels']])
    support = [set(int(b) for b, charge in channel['deposits']) for channel in row['channels']]
    result.update(channels_with_primary=sum(primary in s for s in support),
        channels_without_sdo=sum(not s for s in support),
        channels_with_other=sum(bool(s-{primary}) for s in support))
    if hits:
        starts = np.array([hit['start_local_mm'] for hit in hits])
        ends = np.array([hit['end_local_mm'] for hit in hits])
        midpoints = .5*(starts+ends)
        weights = np.array([hit['energy_loss_MeV'] if hit['energy_loss_MeV'] > 0 else 1. for hit in hits])
        centroid = np.average(midpoints, axis=0, weights=weights)
        measured = row['native_local_mm'][0]
        points = []
        for start, end in zip(starts, ends):
            dz = end[2]-start[2]
            if dz != 0:
                t = -start[2]/dz
                if 0 <= t <= 1:
                    points.append((start+t*(end-start)).tolist())
        result.update(primary_centroid_local_mm=centroid.tolist(),
            measured_minus_primary_centroid_mm=float(measured-centroid[0]),
            primary_midpoint_x_range_mm=float(np.ptp(midpoints[:, 0])),
            primary_endpoints_x_range_mm=float(np.ptp(np.r_[starts[:, 0], ends[:, 0]])),
            primary_time_range_ns=float(np.ptp([hit['time_ns'] for hit in hits])),
            primary_sensor_plane_crossings=points,
            primary_energy_deposit_MeV=float(sum(hit['energy_loss_MeV'] for hit in hits)))
        # Endpoint connectivity describes saved G4 step pieces. Persisted SiHits
        # lack G4 track IDs, so disconnected pieces are not labelled secondaries.
        neighbors = [set() for _ in hits]
        for i in range(len(hits)):
            for j in range(i):
                if min(np.linalg.norm(ends[i]-starts[j]), np.linalg.norm(ends[j]-starts[i])) < 1e-5:
                    neighbors[i].add(j)
                    neighbors[j].add(i)
        remaining = set(range(len(hits)))
        components = 0
        while remaining:
            components += 1
            pending = [remaining.pop()]
            while pending:
                reachable = neighbors[pending.pop()] & remaining
                remaining -= reachable
                pending.extend(reachable)
        result['primary_endpoint_connected_components_10nm'] = components
        if row['majority_barcode'] == primary:
            np.testing.assert_allclose(centroid[:2], row['majority_centroid_local_mm'], rtol=0, atol=1e-12)
            cell = row['centroid_cell']
            if cell is not None:
                centre = cell['centre_local_mm'][0]
                result.update(native_centroid_cell_strip=cell['strip'],
                    centroid_cell_is_cluster_channel=cell['strip'] in result['native_channel_strip_indices'],
                    measured_minus_centroid_cell_mm=float(measured-centre),
                    centroid_cell_minus_truth_centroid_mm=float(centre-centroid[0]))
    return result


def full_seed_checks(original, current):
    old = {r['seed']: r for r in original if r['stage'] == 'seed'}
    new = {r['seed']: r for r in current if r['stage'] == 'seed'}
    assert old.keys() == new.keys()
    for index in old:
        a, b = old[index], new[index]
        for key in ['bound', 'cov', 'surface']:
            assert a[key] == b[key], key
        assert sorted(a['hits']) == sorted(b['hits'])
    return dict(seed_count=len(old), full_mean_covariance_surface_exact=True,
                hit_multisets_exact=True,
                hit_vector_order_changed=sum(old[i]['hits'] != new[i]['hits'] for i in old))


def outcome(folder, task):
    tables = read_tables(folder, task)
    decision = derive(tables, task['dsid'])
    primary = np.flatnonzero(decision['primary'])
    assert len(primary) == 1
    i = int(primary[0])
    chosen = int(decision['chosen'][i])
    quality = quality_record(tables['tracks'], chosen)
    return dict(status=STATUS[int(decision['status'][i])],
        primary_barcode=int(tables['truth']['barcode'][i]),
        final_qop_GeVinv=float(tables['tracks']['qop_GeVinv'][chosen]) if chosen >= 0 else None,
        final_qop_variance_GeVinv2=float(tables['tracks']['qop_variance_GeVinv2'][chosen]) if chosen >= 0 else None,
        quality=quality, failed_quality_conditions=failed_quality(quality))


def reduce(root, output, smoke=False):
    control = root/'ckf_measurement_v1'
    jobs = json.loads((control/('smoke_tasks.json' if smoke else 'all_tasks.json')).read_text())
    grouped = defaultdict(dict)
    for job in jobs:
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0', str(folder)
        task = json.loads((folder/'source_task.json').read_text())
        run = json.loads((folder/'run.json').read_text())
        assert run['success'] and run['isolate_smoothing'] and run['charge_matrix_audit']
        expected = {str((control/'software/lib'/name).resolve()) for name in
                    ['libFaserActsGeometry.so', 'libFaserActsGeometryLib.so', 'libFaserActsKalmanFilter.so']}
        expected.add(str((root/'observer_hitcontrol_v2/lib/libCalypsoSlopeDiagnostics.so').resolve()))
        assert set(run['loaded_acts_libraries']) == expected
        cluster_libraries = {str((root/'cluster_measurement_v1/software/lib'/name).resolve()) for name in
                             ['libFaserSiClusterizationTool.so', 'libFaserSiClusterizationToolLib.so']}
        assert set(run['loaded_cluster_libraries']) == cluster_libraries
        if task['coordinate_mode'] == 'native':
            validation = json.loads((folder/'replay_validation.json').read_text())
            assert validation['all_reconstruction_arrays_identical'] and validation['fields_exact'] == 58
        key = task['dsid'], task['sample_entry']
        grouped[key][task['coordinate_mode']] = (folder, task)
    events, counts, spatial = [], defaultdict(Counter), []
    for (dsid, entry), modes in sorted(grouped.items()):
        assert set(modes) == {'native', 'centroid_fixed', 'ckf_quarter', 'ckf_four', 'ckf_cluster'}
        native_folder, task = modes['native']
        native, original = target_path(native_folder, task, {})
        native_coordinates = coordinates(native_folder)
        primary = native['primary_barcode']
        wrong = bool(task['replay_targets'][0]['wrong'])
        cohort = 'pure_wrong' if wrong else 'normal'
        assert (native['outcome'] == 'wrong') == wrong
        result = dict(dsid=dsid, sample_entry=entry, cohort=cohort,
            native_ordered_measurement_hits=native['ordered_measurement_hits'], modes={})
        for step in native['per_update']:
            hit = native_coordinates[step['hit_id']]
            row = coordinate_features(hit, primary)
            row.update(dsid=dsid, sample_entry=entry, cohort=cohort,
                actual_ckf_sigma_mm=.0231, z_F_mm=step['z_F_mm'],
                update_residual_mm=step['residual_mm'][0],
                update_gain_qop_GeVinv_per_mm=step['gain_qop_GeVinv_per_mm'][0],
                predicted_qop_GeVinv=step['predicted_qop_GeVinv'],
                filtered_qop_GeVinv=step['filtered_qop_GeVinv'])
            spatial.append(row)
        for mode, (folder, current_task) in modes.items():
            current = read_audit(folder/'charge_audit.jsonl.gz')
            seed = full_seed_checks(original, current)
            stats = outcome(folder, current_task)
            trace = target_path(folder, current_task, {})[0] if stats['final_qop_GeVinv'] is not None else None
            coords = coordinates(folder)
            updates_checked = 0
            for record in current:
                if record['stage'] != 'selected_state_trace':
                    continue
                for step in record['steps']:
                    if not step['measurement_flag']:
                        continue
                    hit = coords[step['hit']]
                    scale = {'ckf_quarter': .25, 'ckf_four': 4.}.get(mode, 1.)
                    expected = hit['used_covariance_mm2'][0] if mode == 'ckf_cluster' else .0231**2*scale
                    np.testing.assert_array_equal(step['measurement_cov'], [[expected]])
                    np.testing.assert_array_equal(step['measurement'], [hit['used_local_mm'][0]])
                    operation(step, [], primary, 1 if dsid == '100116' else -1)
                    updates_checked += 1
            if mode in ['ckf_quarter', 'ckf_four', 'ckf_cluster']:
                assert native_coordinates == coords
            same_hits = trace is not None and trace['ordered_measurement_hits'] == native['ordered_measurement_hits']
            result['modes'][mode] = dict(stats, seed_checks=seed, same_ordered_hits=same_hits,
                ordered_hits=trace['ordered_measurement_hits'] if trace else [],
                native_reference_qop_GeVinv=trace['native_reference_qop_GeVinv'] if trace else None,
                final_wrong_entry=trace['final_wrong_entry'] if trace else None,
                last_update_into_wrong=trace['last_update_into_wrong'] if trace else None,
                per_update=trace['per_update'] if trace else [], updates_checked=updates_checked,
                audit=str(folder/'charge_audit.jsonl.gz'))
            group = counts['|'.join([dsid, cohort, mode])]
            group['events'] += 1
            group[stats['status']] += 1
            group['same_ordered_hits'] += same_hits
            if stats['final_qop_GeVinv'] is not None:
                group['final_correct_sign'] += stats['final_qop_GeVinv']*(1 if dsid == '100116' else -1) > 0
        events.append(result)
        print(dsid, entry, cohort, {m: v['status'] for m, v in result['modes'].items()}, flush=True)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/'events.json', events)
    write_json(output/'coordinates.json', spatial)
    write_json(output/'summary.json', dict(events=len(events), source_jobs=len(jobs),
        groups={key: dict(value) for key, value in counts.items()},
        all_native_corrected_reference_fields_exact=True, all_full_seeds_exact=True,
        all_actual_ckf_measurement_variances_verified=True,
        scope='Six pure wrong and 60 original normal controls; corrected material encoding, isolated smoothing ON',
        interpretation='Outcome-selected mechanism cohort. Sensor endpoints are G4 SiHits from the actual input; no local G4 momentum inferred. CKF variance controls retain identical cluster coordinates and errors.'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()
    reduce(args.root, args.output, args.smoke)


if __name__ == '__main__':
    main()
