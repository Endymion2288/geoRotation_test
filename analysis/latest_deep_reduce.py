"""Reduce latest-Calypso states, full 18-sensor covariances, and 52k-seed forensics."""
from collections import defaultdict
import json
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np


SEED_FIELDS = [
    'seedIndex', 'majorityParticleId', 'track_loc0', 'track_loc1',
    'track_phi', 'track_theta', 'track_QoverP', 'track_time',
    't_loc0', 't_loc1', 't_phi', 't_theta', 't_QoverP', 't_time',
]
WRITER_SCALARS = ['event_nr', 'track_nr', 'majorityParticleId', 'chi2Sum', 'NDF', 'rMax_smt']
WRITER_VECTORS = [
    'volume_id', 'layer_id', 't_eQOP', 'eQOP_prt', 'eQOP_flt', 'eQOP_smt',
    'err_eQOP_prt', 'err_eQOP_flt', 'err_eQOP_smt',
]
STATE_FIELDS = [
    'entry', 'candidate', 'seed', 'station', 'layer', 'measurement', 'z_F_mm',
    'p_production_GeV', 'native_truth_qop_GeV_inv', 'q_predicted', 'q_filtered', 'q_smoothed',
    'var_predicted', 'var_filtered', 'var_smoothed', 'previous_q_filtered',
    'dq_transport', 'dq_update', 'dq_smoothing', 'measurement_residual_mm', 'innovation_mm2',
    'gain_qop_per_mm', 'update_identity_error', 'seed_qop', 'seed_truth_qop',
    'github', 'purity', 'good_seed', 'same_seed_truth', 'material_X0', 'sigma_over_truth',
    'truth_charge', 'seed_barcode', 'track_barcode', 'predicted_corr_min_eigenvalue',
    'filtered_corr_min_eigenvalue', 'corr_qop_phi', 'variance_update_identity_error',
    'wafer_side', 'sensor_ordinal_in_layer', 'stereo_sin_alpha',
    'var_pred_loc0', 'var_pred_loc1', 'var_pred_phi', 'var_pred_theta',
    'corr_qop_loc0', 'corr_qop_loc1', 'corr_qop_theta',
]
TRACK_FIELDS = [
    'entry', 'candidate', 'seed', 'seed_barcode', 'track_barcode', 'final_barcode',
    'github', 'purity', 'good_seed', 'p_production_GeV', 'seed_qop', 'seed_truth_qop',
    'refit_ok', 'final_qop_GeV_inv', 'final_p_GeV', 'final_charge',
]
SEED_FORENSIC_FIELDS = [
    'entry', 'seed', 'stations', 'n_clusters', 'chi2', 'chi2_is_nan', 'min_radicand',
    'circle_r_mm', 'circle_cx_mm', 'circle_cy_mm', 'p_fake9_GeV', 'p_centers3_GeV',
    'charge_seed', 'max_shared_prev', 'n_segs', 'all_segs_same_barcode',
    'n_distinct_barcodes', 'has_opp_primary_mix', 'has_secondary_seg', 'min_seg_purity',
    'seg0_barcode', 'seg1_barcode', 'seg2_barcode',
    'phi_chord', 'theta_chord', 'phi_seg0', 'theta_seg0', 'phi_circle_tan', 'theta_line',
    'writer_majority_barcode', 'writer_reco_qop', 'writer_truth_qop',
    'writer_reco_phi', 'writer_truth_phi', 'writer_reco_theta', 'writer_truth_theta',
    'truth_prod_p_GeV', 'truth_pdg',
]
SEED_EVENT_FIELDS = [
    'entry', 'raw_s0', 'raw_s1', 'raw_s2', 'raw_s3',
    'gb_s0', 'gb_s1', 'gb_s2', 'gb_s3',
    'raw_seed_combinations', 'raw_nan_chi2_seeds', 'selected_seeds',
]


def grouped_audit(path):
    if not path.exists():
        return
    rows = []
    with path.open() as stream:
        for line in stream:
            row = json.loads(line)
            if row['stage'] == 'event' and rows:
                yield rows
                rows = []
            rows.append(row)
    if rows:
        yield rows


def seed_good(values):
    q, tq, theta, truth_theta = values[6], values[12], values[5], values[11]
    return (all(map(math.isfinite, [q, tq, theta, truth_theta])) and
            q * tq > 0 and abs(abs(tq / q) - 1) < 0.1 and abs(theta - truth_theta) < 1e-4)


def covariance_health(matrix):
    matrix = np.asarray(matrix, dtype=float)
    diag = np.diag(matrix)
    if not np.all(np.isfinite(matrix)) or np.any(diag <= 0):
        return math.nan
    corr = matrix / np.sqrt(np.outer(diag, diag))
    return float(np.linalg.eigvalsh((corr + corr.T) / 2)[0])


def safe_corr(cov, i, j):
    vi, vj = cov[i, i], cov[j, j]
    if vi > 0 and vj > 0 and math.isfinite(cov[i, j]):
        return float(cov[i, j] / math.sqrt(vi * vj))
    return math.nan


def reduce_run(directory, destination, task):
    import uproot
    directory, destination = Path(directory), Path(destination)
    with uproot.open(directory / 'reference_3station_forward_seed_summary_circleFitTrackSeedTool.root') as source:
        seed_columns = source['tree'].arrays(['event_number', *SEED_FIELDS], library='np')
    seed_events, seed_rows = [], []
    for entry in range(len(seed_columns['event_number'])):
        columns = [seed_columns[name][entry] for name in SEED_FIELDS]
        assert len(set(map(len, columns))) == 1
        seeds = {}
        for values in zip(*columns):
            values = list(map(float, values))
            seeds[int(values[0])] = values
            seed_rows.append([task['skip'] + entry, *values])
        seed_events.append((int(seed_columns['event_number'][entry]), seeds))
    assert len(seed_events) == task['events']

    with uproot.open(directory / 'track_states.root') as source:
        native = source['tree'].arrays(WRITER_SCALARS + WRITER_VECTORS, library='np')
    native_entries = len(native['event_nr'])

    audits = iter(grouped_audit(directory / 'momentum_audit.jsonl'))
    observed = (directory / 'momentum_audit.jsonl').exists()
    writer_scalars, writer_vectors, offsets = [], [], [0]
    states_out, geometry_ids, hit_ids, tracks_out = [], [], [], []
    seed_forensics_out, seed_events_out = [], []
    purity_values, seed_values_out = [], []
    ordinal = 0

    with (directory / 'production_truth.jsonl').open() as stream:
        for local_entry, line in enumerate(stream):
            final_event = json.loads(line)
            event_number, seeds = seed_events[local_entry]
            assert final_event['event'] == event_number
            audit = next(audits) if observed else []
            truth = {r['barcode']: r for r in audit if r['stage'] == 'truth'}
            initial = {r['seed']: r for r in audit if r['stage'] == 'seed'}
            selected = {r['candidate']: r for r in audit if r['stage'] == 'selected'}
            final_trk = {r['candidate']: r for r in audit if r['stage'] == 'final_trk'}
            state_map = defaultdict(list)
            seed_for_map = {}
            for row in audit:
                assert row['source_entry'] == local_entry
                if row['stage'] == 'selected_state':
                    state_map[row['candidate']].append(row)
                elif row['stage'] == 'seed_forensics':
                    seed_for_map[row['seed']] = row
                elif row['stage'] == 'seed_event_summary':
                    r0123 = row['raw_segments_s0123']
                    g0123 = row['gb_segments_s0123']
                    seed_events_out.append([
                        task['skip'] + local_entry,
                        *r0123, *g0123,
                        row['raw_seed_combinations'],
                        row['raw_nan_chi2_seeds'],
                        row['selected_seeds'],
                    ])

            if observed:
                for sidx, svals in sorted(seeds.items()):
                    sf = seed_for_map.get(sidx)
                    if sf is None:
                        continue
                    seg_bcs = [int(x) for x in sf['seg_barcodes']]
                    seg_mh = [int(x) for x in sf['seg_maj_hits']]
                    seg_nc = [int(x) for x in sf['seg_nclusters']]
                    n_segs = len(seg_bcs)
                    distinct_bcs = len({b for b in seg_bcs if b > 0})
                    all_same_bc = int(n_segs >= 2 and all(b > 0 and b == seg_bcs[0] for b in seg_bcs))
                    pdgs = {truth[b]['pdg'] for b in seg_bcs if b in truth and 0 < b < 200000}
                    has_opp_prim = int(11 in pdgs and -11 in pdgs)
                    has_sec = int(any(b == 0 or b % 1000000 >= 200000 for b in seg_bcs))
                    min_pur = min((mh / nc if nc > 0 else 0.0) for mh, nc in zip(seg_mh, seg_nc)) if n_segs else 0.0
                    padded_bcs = (seg_bcs + [-1, -1, -1])[:3]
                    maj_bc = int(svals[1]) if math.isfinite(svals[1]) else -1
                    tr_s = truth.get(maj_bc)
                    prod_p = float(np.linalg.norm(tr_s['momentum_F_MeV']) * 0.001) if tr_s else math.nan
                    prod_pdg = float(tr_s['pdg']) if tr_s else math.nan
                    chi2_val = float(sf['chi2']) if sf['chi2'] is not None else math.nan
                    min_rad = float(sf['min_radicand']) if sf['min_radicand'] is not None else math.nan
                    p_c3 = float(sf['p_centers3_GeV']) if sf['p_centers3_GeV'] is not None else math.nan
                    phi_s0 = float(sf['phi_seg0']) if sf['phi_seg0'] is not None else math.nan
                    th_s0 = float(sf['theta_seg0']) if sf['theta_seg0'] is not None else math.nan
                    phi_ct = float(sf['phi_circle_tan']) if sf['phi_circle_tan'] is not None else math.nan
                    th_ln = float(sf['theta_line']) if sf['theta_line'] is not None else math.nan
                    seed_forensics_out.append([
                        task['skip'] + local_entry, sidx,
                        sf['stations'], sf['n_clusters'], chi2_val, sf['chi2_is_nan'], min_rad,
                        float(sf['circle_r_mm']) if sf['circle_r_mm'] is not None else math.nan,
                        float(sf['circle_cx_mm']) if sf['circle_cx_mm'] is not None else math.nan,
                        float(sf['circle_cy_mm']) if sf['circle_cy_mm'] is not None else math.nan,
                        float(sf['p_fake9_GeV']) if sf['p_fake9_GeV'] is not None else math.nan,
                        p_c3,
                        float(sf['charge_seed']) if sf['charge_seed'] is not None else math.nan,
                        sf['max_shared_prev'], n_segs, all_same_bc,
                        distinct_bcs, has_opp_prim, has_sec, min_pur,
                        *padded_bcs,
                        float(sf['phi_chord']) if sf['phi_chord'] is not None else math.nan,
                        float(sf['theta_chord']) if sf['theta_chord'] is not None else math.nan,
                        phi_s0, th_s0, phi_ct, th_ln,
                        maj_bc, svals[6], svals[12], svals[4], svals[10], svals[5], svals[11],
                        prod_p, prod_pdg,
                    ])

            assert not observed or len(selected) == len(final_event['final'])
            for candidate, final in enumerate(final_event['final']):
                assert ordinal < native_entries
                scalar = [float(native[name][ordinal]) for name in WRITER_SCALARS]
                assert int(scalar[0]) == event_number and int(scalar[1]) == candidate
                arrays = [np.asarray(native[name][ordinal], dtype=np.float64) for name in WRITER_VECTORS]
                assert len(set(map(len, arrays))) == 1
                writer_scalars.append(scalar)
                writer_vectors.extend(np.array(arrays).T.tolist())
                offsets.append(len(writer_vectors))
                github = scalar[4] > 9 and scalar[3] / scalar[4] < 20 and scalar[5] < 90
                purity = final['majority_hits'] / final['clusters'] if final['clusters'] else 0.0
                seed = selected[candidate]['seed'] if observed else -1
                seed_values = seeds[seed] if observed else None
                good = seed_good(seed_values) if observed else False
                purity_values.append(purity)
                seed_values_out.append(int(good))
                if observed:
                    barcode, seed_barcode = int(scalar[2]), int(seed_values[1])
                    tr = truth.get(barcode)
                    p = np.linalg.norm(tr['momentum_F_MeV']) * 0.001 if tr else math.nan
                    charge = (-1.0 if tr['pdg'] > 0 else 1.0) if tr else math.nan
                    ft = final_trk.get(candidate)
                    refit_ok = float(ft['refit']) if ft else math.nan
                    final_mom_vec = ft.get('momentum_A_GeV') if ft else None
                    final_p = float(np.linalg.norm(final_mom_vec)) if final_mom_vec else math.nan
                    final_ch = float(ft['charge']) if ft and ft.get('charge') is not None else math.nan
                    final_qop = final_ch / final_p if final_p > 0 and math.isfinite(final_ch) else math.nan
                    tracks_out.append([
                        task['skip'] + local_entry, candidate, seed, seed_barcode, barcode,
                        final['majority_barcode'], github, purity, good, p,
                        initial[seed]['bound'][4], seed_values[12],
                        refit_ok, final_qop, final_p, final_ch,
                    ])
                    reverse_states = state_map[candidate]
                    measurements = [s for s in reverse_states if s['measurement_flag']]
                    assert len(measurements) == len(arrays[0])
                    native_truth = {}
                    for j, s in enumerate(measurements):
                        assert (s['geometry_id'] >> 56) == int(arrays[0][j])
                        for k, name in enumerate(('predicted', 'filtered', 'smoothed')):
                            np.testing.assert_allclose(s[name][4], arrays[3 + k][j], rtol=1e-6, atol=1e-12)
                        native_truth[s['state']] = arrays[2][j] * 1000.0
                    previous = initial[seed]['bound'][4]
                    layer_ordinals = defaultdict(int)
                    for s in reversed(reverse_states):
                        effective = dict(s)
                        if 'filtered' not in effective:
                            effective['filtered'] = s['predicted']
                            effective['filtered_cov'] = s['predicted_cov']
                        values = [effective.get(n, [math.nan] * 6)[4] for n in ('predicted', 'filtered', 'smoothed')]
                        variances = [effective.get(n, [[math.nan] * 6] * 6)[4][4] for n in
                                     ('predicted_cov', 'filtered_cov', 'smoothed_cov')]
                        qp, qf, qs = values
                        residual, innovation, gain = math.nan, math.nan, math.nan
                        if s.get('measurement_flag') and 'measurement' in s:
                            H = np.asarray(s['projector'], dtype=float)
                            C = np.asarray(s['predicted_cov'], dtype=float)
                            innovation = float((H @ C @ H.T + np.asarray(s['measurement_cov']))[0, 0])
                            residual = float((np.asarray(s['measurement']) - H @ np.asarray(s['predicted']))[0])
                            gain = float((C @ H.T)[4, 0] / innovation)
                        error = (qf - qp) - gain * residual
                        covariance_error = math.nan
                        if math.isfinite(gain):
                            covariance_error = variances[1] - (variances[0] - gain * gain * innovation)
                        cp = np.asarray(effective.get('predicted_cov', [[math.nan] * 6] * 6), dtype=float)
                        cqphi = safe_corr(cp, 4, 2)
                        cqloc0 = safe_corr(cp, 4, 0)
                        cqloc1 = safe_corr(cp, 4, 1)
                        cqtheta = safe_corr(cp, 4, 3)
                        geo = s['geometry_id']
                        station = s.get('wafer_station', {25: 1, 27: 2, 29: 3}.get(geo >> 56, 0))
                        layer = s.get('wafer_layer', ((geo >> 36) & 0xfff) // 2 - 1)
                        wside = float(s['wafer_side']) if 'wafer_side' in s else math.nan
                        if s['measurement_flag']:
                            ord_in_layer = float(layer_ordinals[(station, layer)])
                            layer_ordinals[(station, layer)] += 1
                        else:
                            ord_in_layer = math.nan
                        # Surface local X axis in ACTS A is col(0) = (s00, s10, s20).
                        # Since (x_A, y_A, z_A) = (z_F, y_F, -x_F), component along +x_F is -s20.
                        surf = s['surface']
                        stereo_sin = -float(surf[2][0]) if surf else math.nan
                        states_out.append([
                            task['skip'] + local_entry, candidate, seed, station, layer,
                            s['measurement_flag'], surf[0][3], p, native_truth.get(s['state'], math.nan),
                            *values, *variances, previous, qp - previous, qf - qp, qs - qf,
                            residual, innovation, gain, error,
                            initial[seed]['bound'][4], seed_values[12], github, purity, good, seed_barcode == barcode,
                            s.get('material_full_X0', 0.0),
                            math.sqrt(variances[0]) * p if variances[0] >= 0 else math.nan,
                            charge, seed_barcode, barcode,
                            covariance_health(effective['predicted_cov']),
                            covariance_health(effective['filtered_cov']),
                            cqphi, covariance_error,
                            wside, ord_in_layer, stereo_sin,
                            float(cp[0, 0]), float(cp[1, 1]), float(cp[2, 2]), float(cp[3, 3]),
                            cqloc0, cqloc1, cqtheta,
                        ])
                        geometry_ids.append(geo)
                        hit_ids.append(s.get('hit', 0))
                        previous = qf
                ordinal += 1
    assert ordinal == native_entries
    assert not observed or next(audits, None) is None

    with uproot.recreate(directory / 'production.root') as friend:
        friend.mktree('production', {'final_purity': 'float64', 'seed_good': 'int32'})
        if ordinal:
            friend['production'].extend(dict(
                final_purity=np.asarray(purity_values),
                seed_good=np.asarray(seed_values_out, dtype=np.int32),
            ))

    np.savez_compressed(
        destination / 'writer.npz',
        scalars=np.asarray(writer_scalars),
        values=np.asarray(writer_vectors),
        offsets=offsets,
        scalar_names=WRITER_SCALARS,
        vector_names=WRITER_VECTORS,
        seeds=np.asarray(seed_rows),
        seed_fields=['entry', *SEED_FIELDS],
    )
    if observed:
        np.savez_compressed(
            destination / 'states.npz',
            values=np.asarray(states_out, dtype=np.float64),
            fields=STATE_FIELDS,
            geometry_ids=np.asarray(geometry_ids, dtype=np.uint64),
            hit_ids=np.asarray(hit_ids, dtype=np.uint64),
            tracks=np.asarray(tracks_out, dtype=np.float64),
            track_fields=TRACK_FIELDS,
        )
        np.savez_compressed(
            destination / 'seed_forensics.npz',
            seeds=np.asarray(seed_forensics_out, dtype=np.float64),
            seed_fields=SEED_FORENSIC_FIELDS,
            events=np.asarray(seed_events_out, dtype=np.float64),
            event_fields=SEED_EVENT_FIELDS,
        )

    summary = dict(
        source_events=task['events'],
        writer_tracks=ordinal,
        seed_count=len(seed_rows),
        seed_forensics_count=len(seed_forensics_out),
        audit_states=len(states_out),
        audit_observed=observed,
        first_source_entry=task['skip'],
        last_source_entry=task['skip'] + len(seed_events) - 1,
        seed_writer_event_entries=len(seed_events),
        reduction_version='deep-uproot-v1',
        truth='Native writer SiHit-linked HepMC production p; explicit production and local continuation bounds for diagnostics',
    )
    if task.get('macro') and observed:
        link = directory / 'reference_3station_forward.root'
        if not link.exists():
            link.symlink_to('track_states.root')
        label = task.get('macro_label', 'afterGeo_latest')
        invocation = f'{task["macro"]}("reference_3station_forward.root","native_counts","{label}","t_eQOP",false)'
        subprocess.run(['root', '-l', '-b', '-q', invocation], check=True, cwd=directory)
        shutil.copy2(directory / 'native_counts.root', destination / 'native_counts.root')
    (destination / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary
