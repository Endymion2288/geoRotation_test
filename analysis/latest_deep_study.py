#!/usr/bin/env python3
"""Deep source-code attribution of latest-Calypso seed, 18-sensor S1/S2, and smoothed momentum residuals."""
import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import subprocess

import numpy as np


SENSOR_LABELS = [f'S{st}L{ly}_s{sub}' for st in (1, 2, 3) for ly in (0, 1, 2) for sub in (0, 1)]


def describe(values):
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if not len(finite):
        return dict(raw=int(len(values)), finite=0)
    return dict(
        raw=int(len(values)),
        finite=int(len(finite)),
        mean=float(np.mean(finite)),
        std=float(np.std(finite)),
        quantiles=dict(zip(['p05', 'p16', 'p50', 'p84', 'p95'],
                           np.quantile(finite, [.05, .16, .5, .84, .95]).tolist())),
        below_minus01=float(np.mean(finite < -0.1)),
        above_plus01=float(np.mean(finite > 0.1)),
        below_minus05=float(np.mean(finite < -0.5)),
    )


def describe_p(values):
    result = describe(values)
    edges = np.linspace(-1.2, 1.2, 101)
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite):
        counts, _ = np.histogram(finite, edges)
        if np.any(counts):
            idx = int(np.argmax(counts))
            result['modal_bin'] = [float(edges[idx]), float(edges[idx + 1])]
    return result


def p_residual(q, p):
    with np.errstate(divide='ignore', invalid='ignore'):
        return 1.0 / (np.abs(q) * p) - 1.0


def native_formula_residual(q, truth_q):
    def divide(a, b):
        a, b = np.broadcast_arrays(a, b)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.divide(a, b, out=np.zeros_like(a, dtype=float), where=b != 0)
    p = divide(1.0, np.abs(q))
    truth = divide(1.0, np.abs(truth_q))
    return divide(p - truth, truth)


def event_interval(indicator, keys_2col, source_events, bounded=True):
    keys = np.asarray(keys_2col, dtype=int)[:, :2]
    _, inverse = np.unique(keys, axis=0, return_inverse=True)
    count = np.bincount(inverse)
    passed = np.bincount(inverse, weights=np.asarray(indicator, dtype=float))
    ratio = float(np.sum(passed) / np.sum(count))
    variance = source_events / max(1, source_events - 1) * np.sum((passed - ratio * count) ** 2) / np.sum(count) ** 2
    half = 1.96 * float(np.sqrt(variance))
    lo, hi = ratio - half, ratio + half
    if bounded:
        lo, hi = max(0.0, lo), min(1.0, hi)
    return dict(fraction=ratio, ci95=[lo, hi], unit='source-event cluster delta interval')


def load_chains(chains_dir, allowed_files=None):
    chains = defaultdict(list)
    paths = sorted(Path(chains_dir).glob('file*'))
    assert paths, f'No truth_vertices files in {chains_dir}'
    for p in paths:
        fi = int(p.name[4:])
        if allowed_files is not None and fi not in allowed_files:
            continue
        assert (p / 'exit_code.txt').read_text().strip() == '0', p
        for entry, barcode, pdg, z, mom, pz in np.load(p / 'vertices.npz')['values']:
            chains[(fi, int(entry), int(barcode) % 1000000, int(pdg))].append((float(z), float(mom), float(pz)))
    return chains


def load_deep_group(root, name, files=None):
    paths = sorted((root / 'runs' / name).glob('file*'))
    if files is not None:
        fset = set(files)
        paths = [p for p in paths if int(p.name[4:]) in fset]
    assert paths, f'No run directories found for {name}'
    values, tracks, seeds, identities, hits = [], [], [], [], []
    forensic_seeds, forensic_events = [], []
    coverage = []
    fields = track_fields = fseed_fields = fevent_fields = None
    for path in paths:
        task = json.loads((path / 'source_task.json').read_text())
        assert (path / 'exit_code.txt').read_text().strip() == '0', path
        summary = json.loads((path / 'summary.json').read_text())
        assert summary['source_events'] == task['events']
        run = json.loads((path / 'run.json').read_text())
        assert run['success'] and run['isolate_smoothing']
        assert run['args']['material'] == str(root / 'material-latest-rowmajor.json')
        assert run['calypso_dir'].endswith('/geoRotation_calypso_latest/run')
        assert any('/software_deep_v1/lib/libFaserActsKalmanFilter.so' in p for p in run['loaded_acts_libraries'])
        st = np.load(path / 'states.npz')
        wr = np.load(path / 'writer.npz')
        sf = np.load(path / 'seed_forensics.npz')
        assert fields is None or list(fields) == list(st['fields'])
        fields = list(st['fields'])
        track_fields = list(st['track_fields'])
        fseed_fields = list(sf['seed_fields'])
        fevent_fields = list(sf['event_fields'])
        fi = task['file_index']
        values.append(np.column_stack([np.full(len(st['values']), fi), st['values']]))
        tracks.append(np.column_stack([np.full(len(st['tracks']), fi), st['tracks']]))
        seeds.append(np.column_stack([np.full(len(wr['seeds']), fi), wr['seeds']]))
        forensic_seeds.append(np.column_stack([np.full(len(sf['seeds']), fi), sf['seeds']]))
        forensic_events.append(np.column_stack([np.full(len(sf['events']), fi), sf['events']]))
        identities.append(st['geometry_ids'])
        hits.append(st['hit_ids'])
        coverage.append(dict(file=task['input'], file_index=fi, skip=task['skip'], events=task['events']))
    assert len({c['file'] for c in coverage}) == len(coverage)
    return dict(
        values=np.concatenate(values),
        fields=['file', *fields],
        tracks=np.concatenate(tracks),
        track_fields=['file', *track_fields],
        seeds=np.concatenate(seeds),
        forensic_seeds=np.concatenate(forensic_seeds),
        forensic_seed_fields=['file', *fseed_fields],
        forensic_events=np.concatenate(forensic_events),
        forensic_event_fields=['file', *fevent_fields],
        geometry=np.concatenate(identities),
        hits=np.concatenate(hits),
        paths=paths,
        coverage=coverage,
    )


def subset_files(data, files):
    fset = set(files)
    file_list = list(fset)
    smask = np.isin(data['values'][:, 0], file_list)
    return {
        **data,
        'values': data['values'][smask],
        'geometry': data['geometry'][smask],
        'hits': data['hits'][smask],
        'tracks': data['tracks'][np.isin(data['tracks'][:, 0], file_list)],
        'seeds': data['seeds'][np.isin(data['seeds'][:, 0], file_list)],
        'forensic_seeds': data['forensic_seeds'][np.isin(data['forensic_seeds'][:, 0], file_list)],
        'forensic_events': data['forensic_events'][np.isin(data['forensic_events'][:, 0], file_list)],
        'paths': [p for p in data['paths'] if int(p.name[4:]) in fset],
        'coverage': [c for c in data['coverage'] if c['file_index'] in fset],
    }


def analyze_seed_population(data, chains):
    fs = data['forensic_seeds']
    sf = {n: i for i, n in enumerate(data['forensic_seed_fields'])}
    fe = data['forensic_events']
    ef = {n: i for i, n in enumerate(data['forensic_event_fields'])}

    with np.errstate(divide='ignore', invalid='ignore'):
        q_reco = fs[:, sf['writer_reco_qop']]
        q_truth = fs[:, sf['writer_truth_qop']]
        res_p = np.where((q_reco != 0) & (q_truth != 0) & np.isfinite(q_reco) & np.isfinite(q_truth),
                         np.abs(q_truth / q_reco) - 1.0, np.nan)

    # Match each seed's majority truth particle to continuation chain
    p_min_ratio = np.full(len(fs), np.nan)
    for idx, row in enumerate(fs):
        bc = int(row[sf['writer_majority_barcode']]) if np.isfinite(row[sf['writer_majority_barcode']]) else -1
        pdg = int(row[sf['truth_pdg']]) if np.isfinite(row[sf['truth_pdg']]) else 0
        if bc <= 0 or abs(pdg) != 11:
            continue
        ch = sorted(chains.get((int(row[sf['file']]), int(row[sf['entry']]), bc % 1000000, pdg), []))
        if ch and ch[0][1] > 0:
            in_trk = [pt[1] for pt in ch if pt[0] <= 2500.0]
            p_min_ratio[idx] = min(in_trk) / ch[0][1] if in_trk else 1.0

    valid_truth = np.isfinite(res_p)
    pure_3st = (
        valid_truth
        & (fs[:, sf['stations']] == 3)
        & (fs[:, sf['all_segs_same_barcode']] == 1)
        & (fs[:, sf['has_secondary_seg']] == 0)
        & (fs[:, sf['min_seg_purity']] >= 0.8)
    )
    cat_pure_low_brem = pure_3st & (np.isnan(p_min_ratio) | (p_min_ratio >= 0.95))
    cat_pure_hard_brem = pure_3st & (p_min_ratio < 0.95)
    cat_opp_primary = valid_truth & (~pure_3st) & (fs[:, sf['has_opp_primary_mix']] == 1)
    cat_other_mixed = valid_truth & (~pure_3st) & (fs[:, sf['has_opp_primary_mix']] == 0)
    assert int(cat_pure_low_brem.sum() + cat_pure_hard_brem.sum() +
               cat_opp_primary.sum() + cat_other_mixed.sum()) == int(valid_truth.sum())

    left_spike = valid_truth & (res_p < -0.5)
    n_left = max(1, int(left_spike.sum()))
    categories = {}
    for name, mask in [
        ('all_valid_truth', valid_truth),
        ('pure_3st_low_brem_ge095', cat_pure_low_brem),
        ('pure_3st_hard_brem_lt095', cat_pure_hard_brem),
        ('cross_track_opposite_primary_e_plus_e_minus', cat_opp_primary),
        ('secondary_or_impure_or_2station', cat_other_mixed),
    ]:
        categories[name] = dict(
            count=int(mask.sum()),
            fraction_of_valid=float(mask.sum() / max(1, valid_truth.sum())),
            share_of_below_minus05_spike=float((mask & left_spike).sum() / n_left),
            residual=describe_p(res_p[mask]),
        )

    # Compare 9 fakePositions vs 3 station centers on pure_3st_low_brem
    p9 = fs[cat_pure_low_brem, sf['p_fake9_GeV']]
    p3 = fs[cat_pure_low_brem, sf['p_centers3_GeV']]
    fake_vs_centers = describe(p9 / p3 - 1.0)

    # Shared clusters with earlier selected seeds (CircleFitTrackSeedTool.cxx:162 > 6 check)
    shared = fs[:, sf['max_shared_prev']]
    sharing_summary = dict(
        zero_shared=int(np.sum(shared == 0)),
        shared_1_to_5=int(np.sum((shared >= 1) & (shared <= 5))),
        shared_exact_6_full_station=int(np.sum(shared == 6)),
        shared_gt_6=int(np.sum(shared > 6)),
        exact_6_below_minus05_fraction=float(np.mean(res_p[valid_truth & (shared == 6)] < -0.5))
        if np.any(valid_truth & (shared == 6)) else 0.0,
    )

    # NaN chi2 / negative radicand in Seed::getY (CircleFitTrackSeedTool.cxx:356)
    nan_chi2 = fs[:, sf['chi2_is_nan']] > 0
    neg_rad = fs[:, sf['min_radicand']] < 0
    nan_summary = dict(
        raw_combinations_total=int(np.sum(fe[:, ef['raw_seed_combinations']])),
        raw_nan_chi2_total=int(np.sum(fe[:, ef['raw_nan_chi2_seeds']])),
        selected_seeds_total=int(len(fs)),
        selected_missing_truth=int(np.sum(~valid_truth)),
        selected_nan_chi2=int(np.sum(nan_chi2)),
        selected_negative_getY_radicand=int(np.sum(neg_rad)),
        events_ge3_gb_segments_any_station=int(np.sum(
            (fe[:, ef['gb_s1']] >= 3) | (fe[:, ef['gb_s2']] >= 3) | (fe[:, ef['gb_s3']] >= 3)
        )),
        mean_raw_segments_s123=[float(np.mean(fe[:, ef[k]])) for k in ('raw_s1', 'raw_s2', 'raw_s3')],
        mean_gb_segments_s123=[float(np.mean(fe[:, ef[k]])) for k in ('gb_s1', 'gb_s2', 'gb_s3')],
    )

    # Seed direction bias: S1->S2 chord vs S1 segment vs circle tangent at origin
    t_phi = fs[:, sf['writer_truth_phi']]
    t_theta = fs[:, sf['writer_truth_theta']]
    ele_mask = cat_pure_low_brem & (fs[:, sf['truth_pdg']] == 11)
    pos_mask = cat_pure_low_brem & (fs[:, sf['truth_pdg']] == -11)
    direction_summary = dict(
        pure_3st_low_brem=dict(
            chord_dphi_mrad=describe((fs[cat_pure_low_brem, sf['phi_chord']] - t_phi[cat_pure_low_brem]) * 1e3),
            seg0_dphi_mrad=describe((fs[cat_pure_low_brem, sf['phi_seg0']] - t_phi[cat_pure_low_brem]) * 1e3),
            circle_tan_dphi_mrad=describe((fs[cat_pure_low_brem, sf['phi_circle_tan']] - t_phi[cat_pure_low_brem]) * 1e3),
            chord_dtheta_mrad=describe((fs[cat_pure_low_brem, sf['theta_chord']] - t_theta[cat_pure_low_brem]) * 1e3),
            seg0_dtheta_mrad=describe((fs[cat_pure_low_brem, sf['theta_seg0']] - t_theta[cat_pure_low_brem]) * 1e3),
            linear_fit_dtheta_mrad=describe((fs[cat_pure_low_brem, sf['theta_line']] - t_theta[cat_pure_low_brem]) * 1e3),
        ),
        electron_chord_dphi_mrad=describe((fs[ele_mask, sf['phi_chord']] - t_phi[ele_mask]) * 1e3),
        positron_chord_dphi_mrad=describe((fs[pos_mask, sf['phi_chord']] - t_phi[pos_mask]) * 1e3),
        electron_circle_tan_dphi_mrad=describe((fs[ele_mask, sf['phi_circle_tan']] - t_phi[ele_mask]) * 1e3),
        positron_circle_tan_dphi_mrad=describe((fs[pos_mask, sf['phi_circle_tan']] - t_phi[pos_mask]) * 1e3),
    )

    return dict(
        categories=categories,
        fake9_vs_centers3_relative_diff=fake_vs_centers,
        cluster_sharing=sharing_summary,
        combinatorics_and_nan_chi2=nan_summary,
        direction_estimators=direction_summary,
        masks=dict(
            valid_truth=valid_truth,
            pure_low_brem=cat_pure_low_brem,
            pure_hard_brem=cat_pure_hard_brem,
            opp_primary=cat_opp_primary,
            other_mixed=cat_other_mixed,
            ele=ele_mask,
            pos=pos_mask,
        ),
        res_p=res_p,
    )


def analyze_18_sensors(data, require_good=True, same_truth=True):
    v = data['values']
    f = {n: i for i, n in enumerate(data['fields'])}
    mask = (v[:, f['github']] > 0) & (v[:, f['purity']] > 0.7) & (v[:, f['measurement']] > 0)
    if require_good:
        mask &= v[:, f['good_seed']] > 0
    if same_truth:
        mask &= v[:, f['same_seed_truth']] > 0
    rows = v[mask]
    out = {}
    for st in (1, 2, 3):
        for ly in (0, 1, 2):
            for sub in (0, 1):
                key = f'S{st}L{ly}_s{sub}'
                sub_rows = rows[
                    (rows[:, f['station']] == st)
                    & (rows[:, f['layer']] == ly)
                    & (rows[:, f['sensor_ordinal_in_layer']] == sub)
                ]
                if not len(sub_rows):
                    out[key] = dict(count=0)
                    continue
                p = sub_rows[:, f['p_production_GeV']]
                tq = sub_rows[:, f['truth_charge']] / p
                rp = p_residual(sub_rows[:, f['q_predicted']], p)
                rf = p_residual(sub_rows[:, f['q_filtered']], p)
                rs = p_residual(sub_rows[:, f['q_smoothed']], p)
                res_mm = sub_rows[:, f['measurement_residual_mm']]
                inn_mm = np.sqrt(np.maximum(0.0, sub_rows[:, f['innovation_mm2']]))
                with np.errstate(divide='ignore', invalid='ignore'):
                    pull = np.where(inn_mm > 0, res_mm / inn_mm, np.nan)
                out[key] = dict(
                    count=int(len(sub_rows)),
                    predicted_p_res=describe_p(rp),
                    filtered_p_res=describe_p(rf),
                    smoothed_p_res=describe_p(rs),
                    predicted_curv_err=describe(sub_rows[:, f['q_predicted']] / tq - 1.0),
                    filtered_curv_err=describe(sub_rows[:, f['q_filtered']] / tq - 1.0),
                    smoothed_curv_err=describe(sub_rows[:, f['q_smoothed']] / tq - 1.0),
                    wrong_charge_prt=float(np.mean(sub_rows[:, f['q_predicted']] * tq < 0)),
                    wrong_charge_flt=float(np.mean(sub_rows[:, f['q_filtered']] * tq < 0)),
                    wrong_charge_smt=float(np.mean(sub_rows[:, f['q_smoothed']] * tq < 0)),
                    residual_mm=describe(res_mm),
                    innovation_sigma_mm=describe(inn_mm),
                    pull=describe(pull),
                    gain_qop_per_mm=describe(sub_rows[:, f['gain_qop_per_mm']]),
                    sigma_qop_over_truth=describe(np.sqrt(np.maximum(0.0, sub_rows[:, f['var_predicted']])) / np.abs(tq)),
                    sigma_loc0_mm=describe(np.sqrt(np.maximum(0.0, sub_rows[:, f['var_pred_loc0']]))),
                    sigma_loc1_mm=describe(np.sqrt(np.maximum(0.0, sub_rows[:, f['var_pred_loc1']]))),
                    corr_qop_phi=describe(sub_rows[:, f['corr_qop_phi']]),
                    corr_qop_loc1=describe(sub_rows[:, f['corr_qop_loc1']]),
                )
    return out


def transitions_table(data, require_good=True, same_truth=False):
    v = data['values']
    f = {n: i for i, n in enumerate(data['fields'])}
    keep = (v[:, f['github']] > 0) & (v[:, f['purity']] > 0.7)
    if require_good:
        keep &= v[:, f['good_seed']] > 0
    if same_truth:
        keep &= v[:, f['same_seed_truth']] > 0
    indices = np.flatnonzero(keep)
    if not len(indices):
        return np.zeros((0, 17)), [], []
    boundaries = np.r_[0, 1 + np.flatnonzero(np.any(np.diff(v[indices, :3], axis=0) != 0, axis=1)), len(indices)]
    output, paths = [], []
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        ii = indices[start:end]
        a = v[ii]
        meas = a[:, f['measurement']] > 0
        s1 = np.flatnonzero(meas & (a[:, f['station']] == 1))
        s2 = np.flatnonzero(meas & (a[:, f['station']] == 2))
        s3 = np.flatnonzero(meas & (a[:, f['station']] == 3))
        if not len(s1) or not len(s2):
            continue
        first2, last1 = int(s2[0]), int(s1[-1])
        row1, row2 = a[last1], a[first2]
        p = row2[f['p_production_GeV']]
        tq = row2[f['truth_charge']] / p
        s3_last_flt = a[int(s3[-1]), f['q_filtered']] if len(s3) else np.nan
        s3_last_smt = a[int(s3[-1]), f['q_smoothed']] if len(s3) else np.nan
        output.append([
            row2[0], row2[f['entry']], row2[f['candidate']], row2[f['seed']],
            p, tq, row2[f['seed_qop']], row1[f['q_filtered']],
            row2[f['q_predicted']], row2[f['q_filtered']],
            a[s2[1], f['q_predicted']] if len(s2) > 1 else np.nan,
            row2[f['q_smoothed']], s3_last_flt, s3_last_smt,
            row2[f['same_seed_truth']], row2[f['track_barcode']],
            a[0, f['z_F_mm']],
        ])
        paths.append((
            tuple(v[ii[meas], f['station']].astype(int)),
            tuple(int(x) for x in data['hits'][ii[meas]]),
        ))
    cols = [
        'file', 'entry', 'candidate', 'seed', 'p', 'truth_qop', 'seed_qop', 's1_last_flt',
        's2_first_prt', 's2_first_flt', 's2_second_prt', 's2_first_smt',
        's3_last_flt', 's3_last_smt', 'same_truth', 'track_barcode', 'z_first_mm',
    ]
    return np.asarray(output, dtype=float).reshape(-1, len(cols)), cols, paths


def analyze_ungated_and_refit(data, chains):
    rows, cols, _ = transitions_table(data, require_good=False, same_truth=True)
    c = {n: i for i, n in enumerate(cols)}
    p_min_all, p_pre_s2, n_vtx = [], [], []
    for r in rows:
        tq = r[c['truth_qop']]
        if not np.isfinite(tq) or tq == 0:
            p_min_all.append(np.nan)
            p_pre_s2.append(np.nan)
            n_vtx.append(0)
            continue
        pdg = int(-11 * np.sign(tq))
        ch = sorted(chains.get((int(r[0]), int(r[1]), int(r[c['track_barcode']]) % 1000000, pdg), []))
        n_vtx.append(len(ch))
        if ch and ch[0][1] > 0:
            in_trk = [pt[1] for pt in ch if pt[0] <= 2500.0]
            pre2 = [pt[1] for pt in ch if pt[0] <= 1200.0]
            p_min_all.append(min(in_trk) / ch[0][1] if in_trk else 1.0)
            p_pre_s2.append(min(pre2) / ch[0][1] if pre2 else 1.0)
        else:
            p_min_all.append(np.nan)
            p_pre_s2.append(np.nan)
    p_min_all = np.asarray(p_min_all)
    p_pre_s2 = np.asarray(p_pre_s2)
    n_vtx = np.asarray(n_vtx)

    s2_smt_prod = p_residual(rows[:, c['s2_first_smt']], rows[:, c['p']])
    s2_smt_local = p_residual(rows[:, c['s2_first_smt']], rows[:, c['p']] * p_pre_s2)
    curv_err = rows[:, c['s2_first_smt']] / rows[:, c['truth_qop']] - 1.0

    strata = {}
    for name, m in [
        ('all_ungated_same_truth', np.ones(len(rows), bool)),
        ('no_brem_pmin_ge_0995', (n_vtx > 1) & (p_min_all >= 0.995)),
        ('mild_brem_095_to_0995', (n_vtx > 1) & (p_min_all >= 0.95) & (p_min_all < 0.995)),
        ('hard_brem_pmin_lt_095', (n_vtx > 1) & (p_min_all < 0.95)),
        ('pre_s2_hard_brem_lt_090', (n_vtx > 1) & (p_pre_s2 < 0.90)),
    ]:
        strata[name] = dict(
            tracks=int(m.sum()),
            s2_smt_vs_production=describe_p(s2_smt_prod[m]),
            s2_smt_vs_pre_s2_upper=describe_p(s2_smt_local[m]),
            signed_qop_err_vs_production=describe(curv_err[m]),
        )

    # Compare CKF selected smoothed vs post-refit final Trk::Track
    tr = data['tracks']
    tf = {n: i for i, n in enumerate(data['track_fields'])}
    lookup_smt = {(int(r[0]), int(r[1]), int(r[2])): r[c['s2_first_smt']] for r in rows}
    keep_tr = (tr[:, tf['github']] > 0) & (tr[:, tf['purity']] > 0.7) & (tr[:, tf['seed_barcode']] == tr[:, tf['track_barcode']])
    ckf_q, fin_q, prod_p, ref_ok, good_flag = [], [], [], [], []
    for row in tr[keep_tr]:
        key = (int(row[0]), int(row[tf['entry']]), int(row[tf['candidate']]))
        if key not in lookup_smt:
            continue
        ckf_q.append(lookup_smt[key])
        fin_q.append(row[tf['final_qop_GeV_inv']])
        prod_p.append(row[tf['p_production_GeV']])
        ref_ok.append(row[tf['refit_ok']])
        good_flag.append(row[tf['good_seed']])
    ckf_q = np.asarray(ckf_q)
    fin_q = np.asarray(fin_q)
    prod_p = np.asarray(prod_p)
    ref_ok = np.asarray(ref_ok)
    good_flag = np.asarray(good_flag) > 0

    refit_summary = dict(
        ungated_same_truth=dict(
            tracks=int(len(ckf_q)),
            refit_success_fraction=float(np.mean(ref_ok == 1)),
            ckf_s2_smt_p_res=describe_p(p_residual(ckf_q, prod_p)),
            final_refit_p_res=describe_p(p_residual(fin_q, prod_p)),
            relative_p_diff_final_vs_ckf=describe(np.abs(ckf_q / fin_q) - 1.0),
            same_charge_fraction=float(np.mean(ckf_q * fin_q > 0)),
        ),
        good_seed_same_truth=dict(
            tracks=int(good_flag.sum()),
            refit_success_fraction=float(np.mean(ref_ok[good_flag] == 1)) if good_flag.any() else 0.0,
            ckf_s2_smt_p_res=describe_p(p_residual(ckf_q[good_flag], prod_p[good_flag])),
            final_refit_p_res=describe_p(p_residual(fin_q[good_flag], prod_p[good_flag])),
            relative_p_diff_final_vs_ckf=describe(np.abs(ckf_q[good_flag] / fin_q[good_flag]) - 1.0),
        ),
    )
    return dict(
        radiation_strata=strata,
        refit_comparison=refit_summary,
        arrays=dict(
            s2_smt_prod=s2_smt_prod,
            s2_smt_local=s2_smt_local,
            n_vtx=n_vtx,
            p_min_all=p_min_all,
            ckf_p_res_good=p_residual(ckf_q[good_flag], prod_p[good_flag]),
            fin_p_res_good=p_residual(fin_q[good_flag], prod_p[good_flag]),
        ),
    )


def summarize_group(data, chains, files_label):
    tr = data['tracks']
    tf = {n: i for i, n in enumerate(data['track_fields'])}
    source_events = sum(c['events'] for c in data['coverage'])
    gh = tr[:, tf['github']] > 0
    pur = gh & (tr[:, tf['purity']] > 0.7)
    same = pur & (tr[:, tf['seed_barcode']] == tr[:, tf['track_barcode']])
    good_all = pur & (tr[:, tf['good_seed']] > 0)
    good_same = same & (tr[:, tf['good_seed']] > 0)

    sensors_good_same = analyze_18_sensors(data, require_good=True, same_truth=True)
    t_good_same, tcols, _ = transitions_table(data, require_good=True, same_truth=True)
    tc = {n: i for i, n in enumerate(tcols)}
    t_ungated_same, _, _ = transitions_table(data, require_good=False, same_truth=True)

    # Also compute no-brem ungated median
    ungated_refit = analyze_ungated_and_refit(data, chains)

    return dict(
        files_label=files_label,
        file_count=len(data['coverage']),
        source_events=source_events,
        track_counts=dict(
            writer_tracks=int(len(tr)),
            github_cuts=int(gh.sum()),
            purity_gt_07=int(pur.sum()),
            same_seed_track_truth=int(same.sum()),
            good_seed_all=int(good_all.sum()),
            good_seed_same_truth=int(good_same.sum()),
            seed_final_majority_changed_in_good_seed=int((good_all & (tr[:, tf['seed_barcode']] != tr[:, tf['track_barcode']])).sum()),
        ),
        good_seed_same_truth_transitions={
            name: describe_p(p_residual(t_good_same[:, tc[name]], t_good_same[:, tc['p']]))
            for name in ('seed_qop', 's1_last_flt', 's2_first_prt', 's2_first_flt', 's2_second_prt', 's2_first_smt', 's3_last_flt', 's3_last_smt')
        } if len(t_good_same) else {},
        ungated_same_truth_transitions={
            name: describe_p(p_residual(t_ungated_same[:, tc[name]], t_ungated_same[:, tc['p']]))
            for name in ('seed_qop', 's1_last_flt', 's2_first_prt', 's2_first_flt', 's2_second_prt', 's2_first_smt', 's3_last_flt', 's3_last_smt')
        } if len(t_ungated_same) else {},
        eighteen_sensors_good_same=sensors_good_same,
        ungated_and_refit=ungated_refit['radiation_strata'],
        refit_comparison=ungated_refit['refit_comparison'],
    )


def paired_path_comparison(base_data, ctrl_data, require_good_in_base=True):
    # Match on (file, entry, seed, track_barcode, hit_sequence) using UNGATED transitions
    # so we can also compare identical trajectories even when seed B changes the seed-p gate!
    a, cols, paths_a = transitions_table(base_data, require_good=False, same_truth=True)
    b, _, paths_b = transitions_table(ctrl_data, require_good=False, same_truth=True)
    c = {n: i for i, n in enumerate(cols)}
    base_good_keys = set()
    if require_good_in_base:
        ag, _, pg = transitions_table(base_data, require_good=True, same_truth=True)
        for r, pth in zip(ag, pg):
            base_good_keys.add((int(r[0]), int(r[1]), int(r[c['seed']]), int(r[c['track_barcode']]), pth))

    def index_by_key(rows, paths):
        d = {}
        for r, pth in zip(rows, paths):
            k = (int(r[0]), int(r[1]), int(r[c['seed']]), int(r[c['track_barcode']]), pth)
            if not require_good_in_base or k in base_good_keys:
                d[k] = r
        return d

    da, db = index_by_key(a, paths_a), index_by_key(b, paths_b)
    common = sorted(da.keys() & db.keys())
    source_events = sum(x['events'] for x in base_data['coverage'])
    out = dict(
        base_eligible=len(da),
        ctrl_matched_same_hits=len(common),
        retention_fraction=float(len(common) / max(1, len(da))),
    )
    if common:
        xa = np.asarray([da[k] for k in common])
        xb = np.asarray([db[k] for k in common])
        for comp in ('s2_first_prt', 's2_first_flt', 's2_second_prt', 's2_first_smt'):
            ra = p_residual(xa[:, c[comp]], xa[:, c['p']])
            rb = p_residual(xb[:, c[comp]], xb[:, c['p']])
            out[comp] = dict(
                baseline=describe_p(ra),
                control=describe_p(rb),
                paired_diff=describe(rb - ra),
                left_minus01_diff=event_interval((rb < -0.1).astype(float) - (ra < -0.1).astype(float),
                                                 xa, source_events, bounded=False),
            )
    return out


def make_figures(root, output, full_groups, heldout_groups, seed_pop, ungated_arrays, summaries):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    output.mkdir(exist_ok=False)
    bins = {}
    histogram_metadata = {}
    edges = np.linspace(-1.2, 1.2, 101)

    def record_hist(key, values, bin_edges):
        arr = np.asarray(values, dtype=float)
        fin = arr[np.isfinite(arr)]
        counts, _ = np.histogram(fin, bin_edges)
        bins[key + '_counts'] = counts
        bins[key + '_edges'] = np.asarray(bin_edges)
        bins[key + '_sumw2'] = counts.copy()
        histogram_metadata[key] = dict(
            raw=int(len(arr)),
            finite=int(len(fin)),
            underflow=int(np.sum(fin < bin_edges[0])),
            overflow=int(np.sum(fin > bin_edges[-1])),
        )
        return counts, len(fin)

    def save(fig, name, pdf):
        fig.savefig(output / f'{name}.png', dpi=170)
        fig.savefig(output / f'{name}.pdf')
        pdf.savefig(fig)
        plt.close(fig)

    with PdfPages(output / 'diagnostics.pdf') as pdf:
        # Figure 1: Complete 52,338-seed forensic decomposition + B=0.4737 calibration + phi chord vs tangent
        fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), layout='constrained')
        res_p = seed_pop['res_p']
        m = seed_pop['masks']
        c_low, _ = record_hist('seed_cat_pure_low_brem', res_p[m['pure_low_brem']], edges)
        c_brem, _ = record_hist('seed_cat_pure_hard_brem', res_p[m['pure_hard_brem']], edges)
        c_opp, _ = record_hist('seed_cat_opp_primary', res_p[m['opp_primary']], edges)
        c_mix, _ = record_hist('seed_cat_other_mixed', res_p[m['other_mixed']], edges)

        centers = 0.5 * (edges[:-1] + edges[1:])
        width = edges[1] - edges[0]
        b0 = np.zeros_like(c_low, dtype=float)
        for counts_arr, label, color in [
            (c_opp, f"Cross-track e+/e- mix (N={m['opp_primary'].sum()})", '#d62728'),
            (c_mix, f"Secondary / impure / 2-st (N={m['other_mixed'].sum()})", '#ff7f0e'),
            (c_brem, f"Pure 3-st, hard brem >5% (N={m['pure_hard_brem'].sum()})", '#9467bd'),
            (c_low, f"Pure 3-st, low brem <=5% (N={m['pure_low_brem'].sum()})", '#1f77b4'),
        ]:
            axes[0].bar(centers, counts_arr, width=width, bottom=b0, label=label, color=color, alpha=0.88)
            b0 += counts_arr
        for cut in (-0.1, 0.1):
            axes[0].axvline(cut, color='black', ls='--', lw=0.9)
        axes[0].set(
            xlim=(-1.1, 1.1),
            xlabel='(p_seed - p_truth) / p_truth',
            ylabel='Selected seeds / 0.024 bin',
            title='All 50,364 valid-truth seeds stacked by origin',
        )
        axes[0].legend(fontsize=7.2, loc='upper right')

        # Panel B: Held-out files 1-19 pure 3-st seeds: native B=0.55 T vs calibrated B=0.473662 T
        fs_base19 = heldout_groups['baseline']['forensic_seeds']
        fs_cal19 = heldout_groups['seed_b_04737']['forensic_seeds']
        sf = {n: i for i, n in enumerate(heldout_groups['baseline']['forensic_seed_fields'])}
        for arr_fs, key_name, label, color in [
            (fs_base19, 'heldout19_seed_B055', 'Native B = 0.55 T (9 fake pts)', '#1f77b4'),
            (fs_cal19, 'heldout19_seed_B04737', 'Calibrated B = 0.4737 T (files 1-19)', '#2ca02c'),
        ]:
            qr, qt = arr_fs[:, sf['writer_reco_qop']], arr_fs[:, sf['writer_truth_qop']]
            pure = (
                (arr_fs[:, sf['stations']] == 3)
                & (arr_fs[:, sf['all_segs_same_barcode']] == 1)
                & (arr_fs[:, sf['has_secondary_seg']] == 0)
                & (arr_fs[:, sf['min_seg_purity']] >= 0.8)
                & (qr != 0) & (qt != 0) & np.isfinite(qr) & np.isfinite(qt)
            )
            rp = np.abs(qt[pure] / qr[pure]) - 1.0
            cnt, nfin = record_hist(key_name, rp, np.linspace(-0.5, 0.5, 101))
            axes[1].stairs(cnt / max(1, nfin), np.linspace(-0.5, 0.5, 101), label=f'{label} (N={nfin})', color=color, lw=1.4)
        for cut in (-0.1, 0.1):
            axes[1].axvline(cut, color='black', ls='--', lw=0.9)
        axes[1].axvline(0.0, color='gray', ls=':', lw=0.8)
        axes[1].set(
            xlim=(-0.45, 0.45),
            xlabel='(p_seed - p_truth) / p_truth',
            ylabel='Fraction of pure 3-st seeds / 0.01 bin',
            title='Held-out 19 files: CircleFit B-scale calibration',
        )
        axes[1].legend(fontsize=7.5)

        # Panel C: Seed bending angle phi_A bias (S1->S2 chord vs circle tangent at origin)
        fs_all = full_groups['baseline']['forensic_seeds']
        phi_edges = np.linspace(-2.5, 2.5, 101)
        for mask_k, col_name, key_name, label, color, ls in [
            (m['ele'], 'phi_chord', 'phi_chord_ele', 'e- native S1-S2 chord', '#d62728', '-'),
            (m['pos'], 'phi_chord', 'phi_chord_pos', 'e+ native S1-S2 chord', '#1f77b4', '-'),
            (m['pure_low_brem'], 'phi_circle_tan', 'phi_tan_all', 'e+/e- circle tangent at seed plane', '#2ca02c', '--'),
        ]:
            dphi = (fs_all[mask_k, sf[col_name]] - fs_all[mask_k, sf['writer_truth_phi']]) * 1e3
            cnt, nfin = record_hist(key_name, dphi, phi_edges)
            axes[2].stairs(cnt / max(1, nfin), phi_edges, label=f'{label} (med={np.nanmedian(dphi):+.2f} mrad)', color=color, ls=ls, lw=1.3)
        axes[2].axvline(0.0, color='gray', ls=':', lw=0.8)
        axes[2].set(
            xlim=(-2.5, 2.5),
            xlabel='Seed phi_A - propagated truth phi_A [mrad]',
            ylabel='Fraction of seeds / 0.05 mrad bin',
            title='Bending angle at seed surface: chord vs tangent',
        )
        axes[2].legend(fontsize=7.2)
        fig.suptitle('DSID 110051 (A\' -> e+e-) | latest Calypso + row-major material | ACTS frame A, seed plane x_A = z_min - 10 mm\nSource-code decomposition of CircleFitTrackSeedTool seed momentum and angle residuals')
        save(fig, 'seed_population_forensics', pdf)

        # Figure 2: 18-sensor forward evolution across S1, S2, S3 (baseline vs covLoc fix vs qop 1000)
        fig, axes = plt.subplots(1, 3, figsize=(17, 5.4), layout='constrained')
        x_idx = np.arange(len(SENSOR_LABELS))
        style_map = [
            ('baseline', 'Baseline [covLoc swapped, q/p x1e5]', '#d62728', '-o'),
            ('fix_cov_loc_swap', 'Fix covLoc0/covLoc1 swap only [q/p x1e5]', '#ff7f0e', '-s'),
            ('qop_1000_only', 'q/p inflation x1e3 only [covLoc swapped]', '#9467bd', '--^'),
            ('qop_1000_loc_fixed', 'Fix covLoc swap + q/p inflation x1e3', '#2ca02c', '-D'),
        ]
        for gname, glabel, color, fmt in style_map:
            sdict = summaries['full_20files'][gname]['eighteen_sensors_good_same']
            prt_left = [sdict[k]['predicted_p_res']['below_minus01'] * 100 if sdict[k]['count'] else np.nan for k in SENSOR_LABELS]
            pull_std = [sdict[k]['pull']['std'] if sdict[k]['count'] else np.nan for k in SENSOR_LABELS]
            wrg_q = [sdict[k]['wrong_charge_flt'] * 100 if sdict[k]['count'] else np.nan for k in SENSOR_LABELS]
            axes[0].plot(x_idx, prt_left, fmt, ms=4, lw=1.3, color=color, label=glabel)
            axes[1].plot(x_idx, pull_std, fmt, ms=4, lw=1.3, color=color, label=glabel)
            axes[2].plot(x_idx, wrg_q, fmt, ms=4, lw=1.3, color=color, label=glabel)

        for ax in axes:
            ax.axvline(5.5, color='gray', ls='--', lw=0.8)
            ax.axvline(11.5, color='gray', ls='--', lw=0.8)
            ax.set_xticks(x_idx)
            ax.set_xticklabels(SENSOR_LABELS, rotation=65, fontsize=7.5)
            ax.grid(alpha=0.25)
        axes[1].axhline(1.0, color='black', ls=':', lw=0.9, label='Unit pull RMS = 1')
        axes[1].set_ylim(0.0, 1.58)
        axes[0].set(ylabel='Predicted p residual < -0.1 [%]', title='Predicted left-tail fraction by wafer')
        axes[1].set(ylabel='Measurement pull std(r / sigma_inn)', title='Innovation pull RMS by wafer')
        axes[2].set(ylabel='Filtered wrong-charge fraction [%]', title='Filtered wrong-charge rate by wafer')
        axes[0].legend(fontsize=6.8, loc='upper left')
        axes[1].legend(fontsize=6.8, loc='upper right')
        fig.suptitle('DSID 110051 | 20 files, 20,000 source events | latest Calypso | GitHub + purity > 0.7 + good seed + same truth\nAll 18 individual SCT wafer sides in forward CKF order (vertical dashed lines mark S1|S2 and S2|S3)')
        save(fig, 'eighteen_sensors_evolution', pdf)

        # Figure 3: Station 2 first vs second sensor & curvature-to-momentum folding
        fig, axes = plt.subplots(1, 3, figsize=(16, 5.1), layout='constrained')
        for gname, title_tag, ax in [
            ('baseline', 'Baseline (covLoc swapped, q/p x1e5)', axes[0]),
            ('fix_cov_loc_swap', 'Fix covLoc0/covLoc1 swap (q/p x1e5)', axes[1]),
        ]:
            t_arr, tcols, _ = transitions_table(full_groups[gname], require_good=True, same_truth=True)
            tc = {n: i for i, n in enumerate(tcols)}
            for cname, color, lbl in [
                ('s2_first_prt', '#2ca02c', '1st S2 wafer predicted'),
                ('s2_first_flt', '#d62728', '1st S2 wafer filtered'),
                ('s2_second_prt', '#ff7f0e', '2nd S2 wafer predicted'),
                ('s2_first_smt', '#1f77b4', '1st S2 wafer smoothed'),
            ]:
                rp = p_residual(t_arr[:, tc[cname]], t_arr[:, tc['p']])
                cnt, nfin = record_hist(f'{gname}_{cname}', rp, edges)
                ax.stairs(cnt, edges, color=color, label=f'{lbl} (med={np.nanmedian(rp)*100:+.1f}%)', lw=1.3)
            ax.set(xlim=(-1.1, 1.1), xlabel='(p - p_prod) / p_prod', ylabel='Tracks / 0.024 bin', title=f'{title_tag} (N={len(t_arr)})')
            ax.legend(fontsize=7)

        # Panel C: Signed curvature error vs folded p residual on baseline 1st S2 filtered
        t_base, tcols, _ = transitions_table(full_groups['baseline'], require_good=True, same_truth=True)
        tc = {n: i for i, n in enumerate(tcols)}
        curv_prt = t_base[:, tc['s2_first_prt']] / t_base[:, tc['truth_qop']] - 1.0
        curv_flt = t_base[:, tc['s2_first_flt']] / t_base[:, tc['truth_qop']] - 1.0
        curv_smt = t_base[:, tc['s2_first_smt']] / t_base[:, tc['truth_qop']] - 1.0
        c_edges = np.linspace(-18, 18, 91)
        for arr_c, key_c, lbl, col in [
            (curv_prt, 's2_curv_prt', f'1st S2 prt (std={np.nanstd(curv_prt):.2f})', '#2ca02c'),
            (curv_flt, 's2_curv_flt', f'1st S2 flt (std={np.nanstd(curv_flt):.2f})', '#d62728'),
            (curv_smt, 's2_curv_smt', f'1st S2 smt (std={np.nanstd(curv_smt):.2f})', '#1f77b4'),
        ]:
            cnt, _ = record_hist(key_c, arr_c, c_edges)
            axes[2].stairs(cnt, c_edges, color=col, label=lbl, lw=1.3)
        axes[2].axvline(-1.0, color='black', ls='--', lw=0.9, label='q/p = 0 (charge flip boundary)')
        axes[2].set(xlim=(-18, 18), xlabel='Signed curvature error (q/p)/(q/p)_truth - 1', ylabel='Tracks / 0.4 bin',
                    title='Symmetric q/p broadening before 1/|q/p| fold')
        axes[2].legend(fontsize=7)
        fig.suptitle('DSID 110051 | 20 files, 20,000 source events | latest Calypso | same seed/track truth\nWhy Station 2 Layer 0 green curve has a -0.92 spike: two-wafer merging + wide curvature folded by p = 1/|q/p|')
        save(fig, 'station2_wafer_and_folding', pdf)

        # Figure 4: Smoothed bias causal decoupling (held-out 19 files) + ungated bremsstrahlung attribution + CKF vs refit
        fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), layout='constrained')
        narrow_edges = np.linspace(-0.6, 0.6, 101)
        for gname, lbl, col, ls in [
            ('baseline', 'Baseline [B=0.55, q/p x1e5]', '#d62728', '-'),
            ('qop_1000_loc_fixed', 'Fix covLoc + q/p x1e3 [B=0.55]', '#ff7f0e', '--'),
            ('seed_b_04737', 'C++ Seed B=0.4737 T only [q/p x1e5]', '#1f77b4', '-.'),
            ('combined_phys_fix', 'Combined: Seed B=0.4737 + covLoc + q/p x1e3', '#2ca02c', '-'),
        ]:
            t19, tcols, _ = transitions_table(heldout_groups[gname], require_good=True, same_truth=True)
            tc = {n: i for i, n in enumerate(tcols)}
            rp = p_residual(t19[:, tc['s2_first_smt']], t19[:, tc['p']])
            cnt, nfin = record_hist(f'heldout19_smt_{gname}', rp, narrow_edges)
            axes[0].stairs(cnt / max(1, nfin), narrow_edges, color=col, ls=ls, lw=1.4,
                           label=f'{lbl} (N={nfin}, med={np.nanmedian(rp)*100:+.1f}%)')
        axes[0].axvline(0.0, color='gray', ls=':', lw=0.8)
        axes[0].set(xlim=(-0.5, 0.5), xlabel='First S2 smoothed (p - p_prod) / p_prod',
                    ylabel='Fraction of tracks / 0.012 bin',
                    title='Held-out 19 files: original good-seed gate')
        axes[0].legend(fontsize=6.8)

        # Panel B: Ungated 20-file smoothed residual decomposed by electron bremsstrahlung
        ua = ungated_arrays
        nv = ua['n_vtx']
        pm = ua['p_min_all']
        for mask_u, arr_u, key_u, lbl, col, ls in [
            (np.ones(len(nv), bool), ua['s2_smt_prod'], 'ungated_all_prod', 'All ungated vs p_prod', '#1f77b4', '-'),
            ((nv > 1) & (pm >= 0.995), ua['s2_smt_prod'], 'ungated_nobrem_prod', 'No-brem (p_min/p_prod>=0.995) vs p_prod', '#2ca02c', '-'),
            ((nv > 1) & (pm < 0.95), ua['s2_smt_prod'], 'ungated_hardbrem_prod', 'Hard-brem (p_min/p_prod<0.95) vs p_prod', '#d62728', '--'),
            ((nv > 1) & (pm < 0.95), ua['s2_smt_local'], 'ungated_hardbrem_local', 'Hard-brem vs pre-S2 local truth upper', '#9467bd', ':'),
        ]:
            vals = arr_u[mask_u]
            cnt, nfin = record_hist(key_u, vals, narrow_edges)
            axes[1].stairs(cnt / max(1, nfin), narrow_edges, color=col, ls=ls, lw=1.4,
                           label=f'{lbl} (N={nfin}, med={np.nanmedian(vals)*100:+.1f}%)')
        axes[1].axvline(0.0, color='gray', ls=':', lw=0.8)
        axes[1].set(xlim=(-0.55, 0.55), xlabel='First S2 smoothed momentum residual',
                    ylabel='Fraction of tracks / 0.012 bin',
                    title='Ungated 20 files: electron bremsstrahlung effect')
        axes[1].legend(fontsize=6.7)

        # Panel C: CKF selected smoothed vs post-refit final Trk::Track
        for arr_r, key_r, lbl, col, ls in [
            (ua['ckf_p_res_good'], 'ckf_vs_refit_ckf', 'CKF selected S2 smoothed', '#1f77b4', '-'),
            (ua['fin_p_res_good'], 'ckf_vs_refit_final', 'Post-refit final Trk::Track (KalmanFitterTool1)', '#ff7f0e', '--'),
        ]:
            cnt, nfin = record_hist(key_r, arr_r, narrow_edges)
            axes[2].stairs(cnt / max(1, nfin), narrow_edges, color=col, ls=ls, lw=1.4,
                           label=f'{lbl} (N={nfin}, med={np.nanmedian(arr_r)*100:+.2f}%)')
        axes[2].axvline(0.0, color='gray', ls=':', lw=0.8)
        axes[2].set(xlim=(-0.45, 0.45), xlabel='(p - p_prod) / p_prod',
                    ylabel='Fraction of tracks / 0.012 bin',
                    title='CKF smoothed state vs post-refit final track')
        axes[2].legend(fontsize=7.2)
        fig.suptitle('DSID 110051 | latest Calypso | Decoupling smoothed momentum shift: seed-B selection bias vs electron bremsstrahlung vs refit')
        save(fig, 'smoothed_bias_decomposition', pdf)

    np.savez_compressed(output / 'bin_counts.npz', default_edges=edges, **bins)
    (output / 'histogram_metadata.json').write_text(json.dumps(histogram_metadata, indent=2) + '\n')


def export_reference_nine_layers(root, output, full_groups, heldout_groups, summaries):
    import ROOT
    ROOT.gROOT.SetBatch(True)

    specs = [
        ('baseline', full_groups['baseline'], 'latest baseline', '20 files / 20000 source events'),
        ('fix_cov_loc_swap', full_groups['fix_cov_loc_swap'], 'covLoc swap fix', '20 files / 20000 source events'),
        ('qop_1000_only', full_groups['qop_1000_only'], 'q/p prior x1000', '20 files / 20000 source events'),
        ('qop_1000_loc_fixed', full_groups['qop_1000_loc_fixed'], 'covLoc+q/p x1000', '20 files / 20000 source events'),
        ('seed_b_04737_heldout19', heldout_groups['seed_b_04737'], 'Seed B=0.4737T', '19 held-out files / 19000 source events (file 0 excluded)'),
        ('combined_phys_fix_heldout19', heldout_groups['combined_phys_fix'], 'Combined (19f)', '19 held-out files / 19000 source events (file 0 excluded)'),
    ]
    names = [f'h_{component}_{station}_{layer}' for component in range(3)
             for station in (1, 2, 3) for layer in (0, 1, 2)]
    reference_figures = []
    raw_entries = {}
    reference_bins = {}

    for key, data, legend_label, scope_str in specs:
        counts_path = output / f'{key}_counts_3station_forward.root'
        prefix = output / f'110051_{key}_nine_layers'
        merged_histograms, expected_entries, expected_counts = {}, {}, {}
        for path in data['paths']:
            source_file = ROOT.TFile.Open(str(path / 'native_counts.root'))
            for name in names:
                sh = source_file.Get(name)
                assert sh and sh.GetNbinsX() == 100
                vals = np.asarray([sh.GetBinContent(i) for i in range(102)])
                if name not in merged_histograms:
                    merged_histograms[name] = sh.Clone(name)
                    merged_histograms[name].SetDirectory(0)
                    expected_entries[name] = sh.GetEntries()
                    expected_counts[name] = vals
                else:
                    assert merged_histograms[name].Add(sh)
                    expected_entries[name] += sh.GetEntries()
                    expected_counts[name] += vals
            source_file.Close()

        merged_file = ROOT.TFile(str(counts_path), 'RECREATE')
        for name, hist in merged_histograms.items():
            assert hist.GetEntries() == expected_entries[name]
            np.testing.assert_equal([hist.GetBinContent(i) for i in range(102)], expected_counts[name])
            hist.Write(name)
        merged_file.Close()

        env = {
            **os.environ,
            'REFERENCE_SCOPE_LABEL': f'110051 | {scope_str} | ACTS A | HepMC p | GitHub + good seed + purity>0.7',
        }
        invocation = (
            str(root / 'reference/ckf_residual_perLayer.C')
            + '(' + ','.join(json.dumps(x) for x in [str(counts_path), str(prefix), legend_label, 't_eQOP']) + ',true)'
        )
        subprocess.run(['root', '-l', '-b', '-q', invocation], env=env, check=True)
        reference_figures.append(prefix.with_suffix('.pdf'))

        actual = ROOT.TFile.Open(str(prefix.with_suffix('.root')))
        source = ROOT.TFile.Open(str(counts_path))
        raw_entries[key] = {}
        for name in names:
            a, b = actual.Get(name), source.Get(name)
            assert a.GetEntries() == b.GetEntries()
            assert all(a.GetBinContent(i) == b.GetBinContent(i) for i in range(102))
            raw_entries[key][name] = a.GetEntries()
            reference_bins[f'{key}_{name}_counts_with_flow'] = np.asarray([a.GetBinContent(i) for i in range(102)])
            reference_bins[f'{key}_{name}_sumw2_with_flow'] = np.asarray([a.GetBinError(i) ** 2 for i in range(102)])
        actual.Close()
        source.Close()

    np.savez_compressed(output / 'reference_bin_counts.npz', edges=np.linspace(-1.2, 1.2, 101), **reference_bins)
    return reference_figures, raw_entries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--chains-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--groups', nargs='+', default=[
        'baseline', 'fix_cov_loc_swap', 'qop_1000_only',
        'qop_1000_loc_fixed', 'seed_b_04737', 'combined_phys_fix',
    ])
    args = parser.parse_args()

    frozen = json.loads((args.root.parent.parent / 'stage04_seed_purity/20261113_seed_good_v2/inputs.json').read_text())
    expected = {x['path']: x['entries'] for x in frozen['samples'][0]['files']}
    assert len(expected) == 20

    full_groups = {g: load_deep_group(args.root, g) for g in args.groups}
    for g, data in full_groups.items():
        assert len(data['coverage']) == 20, (g, len(data['coverage']))
        assert all(c['skip'] == 0 and c['events'] == expected[c['file']] for c in data['coverage'])

    chains_all = load_chains(args.chains_dir)
    heldout_files = set(range(1, 20))
    chains_heldout = load_chains(args.chains_dir, allowed_files=heldout_files)
    heldout_groups = {g: subset_files(d, heldout_files) for g, d in full_groups.items()}

    seed_pop = analyze_seed_population(full_groups['baseline'], chains_all)
    ungated_full = analyze_ungated_and_refit(full_groups['baseline'], chains_all)

    full_non_calibrated = ['baseline', 'fix_cov_loc_swap', 'qop_1000_only', 'qop_1000_loc_fixed']
    summaries = dict(
        seed_population_forensics_20files={k: v for k, v in seed_pop.items() if k not in ('masks', 'res_p')},
        full_20files={
            g: summarize_group(full_groups[g], chains_all, 'all 20 files [0..19], 20000 events')
            for g in full_non_calibrated
        },
        heldout_19files_excluding_file0={
            g: summarize_group(heldout_groups[g], chains_heldout, 'held-out 19 files [1..19], 19000 events (file 0 excluded)')
            for g in args.groups
        },
        paired_same_path_20files={
            g: paired_path_comparison(full_groups['baseline'], full_groups[g], require_good_in_base=True)
            for g in full_non_calibrated if g != 'baseline'
        },
        paired_same_path_heldout19_base_good_gate={
            g: paired_path_comparison(heldout_groups['baseline'], heldout_groups[g], require_good_in_base=True)
            for g in args.groups if g != 'baseline'
        },
        paired_same_path_heldout19_all_ungated={
            g: paired_path_comparison(heldout_groups['baseline'], heldout_groups[g], require_good_in_base=False)
            for g in args.groups if g != 'baseline'
        },
    )

    make_figures(args.root, args.output, full_groups, heldout_groups, seed_pop, ungated_full['arrays'], summaries)
    (args.output / 'summary.json').write_text(json.dumps(summaries, indent=2) + '\n')

    ref_pdfs, raw_entries = export_reference_nine_layers(args.root, args.output, full_groups, heldout_groups, summaries)
    import sys
    sys.path.insert(0, str(args.root.parent.parent / 'stage04_seed_purity/20261113_seed_good_v2/software_tools'))
    from pypdf import PdfReader, PdfWriter
    from PIL import Image

    merged = PdfWriter()
    for p in ref_pdfs:
        merged.append(str(p))
    merged.append(str(args.output / 'diagnostics.pdf'))
    with (args.output / 'report.pdf').open('wb') as stream:
        merged.write(stream)

    index = []
    for png in sorted(args.output.glob('*.png')):
        with Image.open(png) as img:
            img.load()
            sz = list(img.size)
        pdf_p = png.with_suffix('.pdf')
        assert len(PdfReader(pdf_p).pages) == 1
        index.append(dict(png=png.name, pdf=pdf_p.name, pixels=sz))
    (args.output / 'figure_validation.json').write_text(json.dumps(dict(
        figure_index=index,
        pdf_pages=len(PdfReader(args.output / 'report.pdf').pages),
        reference_histograms=27 * len(ref_pdfs),
        reference_bins_with_flow=2754 * len(ref_pdfs),
        reference_raw_entries=raw_entries,
        reference_histograms_exact=True,
    ), indent=2) + '\n')


if __name__ == '__main__':
    main()
