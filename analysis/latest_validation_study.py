#!/usr/bin/env python3
"""Multi-core aggregation and English figure suite for latest Calypso seed fixes, particle hypothesis, and seeding checks."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


ALL_GROUPS = [
    'C0_baseline',
    'C1_cov_loc',
    'C2_orig_var',
    'C3_seed_b',
    'C4_no_nan',
    'C5_local_tan',
    'H1_electron_hypo',
    'S3_strip_only',
    'S1_dedup_prefilter',
    'S2_strip_field_fit',
]

SAMPLE_ORDER = ['110051', '110051_heldout19', '110033', '100116', '100117', '100130']
RAW_SAMPLES = ['110051', '110033', '100116', '100117', '100130']


def wrap_pi(dphi):
    return (dphi + np.pi) % (2.0 * np.pi) - np.pi


def robust_stats(values):
    arr = np.asarray(values, dtype=np.float64)
    fin = arr[np.isfinite(arr)]
    if len(fin) == 0:
        return {
            'raw': int(len(arr)), 'finite': 0,
            'mean': None, 'std': None, 'median': None,
            'p16': None, 'p84': None, 'half68': None,
            'p05': None, 'p95': None,
            'below_m01': None, 'above_p01': None, 'below_m05': None,
        }
    q05, q16, q50, q84, q95 = np.quantile(fin, [0.05, 0.16, 0.50, 0.84, 0.95])
    return {
        'raw': int(len(arr)),
        'finite': int(len(fin)),
        'mean': float(np.mean(fin)),
        'std': float(np.std(fin)),
        'median': float(q50),
        'p16': float(q16),
        'p84': float(q84),
        'half68': float(0.5 * (q84 - q16)),
        'p05': float(q05),
        'p95': float(q95),
        'below_m01': float(np.mean(fin < -0.1)),
        'above_p01': float(np.mean(fin > +0.1)),
        'below_m05': float(np.mean(fin < -0.5)),
    }


def pull_stats(values):
    arr = np.asarray(values, dtype=np.float64)
    fin = arr[np.isfinite(arr)]
    if len(fin) == 0:
        return {'finite': 0, 'median': None, 'mean': None, 'rms': None, 'sigma68': None}
    q16, q50, q84 = np.quantile(fin, [0.16, 0.50, 0.84])
    clip = fin[(fin > -8.0) & (fin < 8.0)]
    return {
        'finite': int(len(fin)),
        'median': float(q50),
        'mean_clipped8': float(np.mean(clip)) if len(clip) else float(np.mean(fin)),
        'rms_clipped8': float(np.std(clip)) if len(clip) else float(np.std(fin)),
        'rms_raw': float(np.std(fin)),
        'sigma68': float(0.5 * (q84 - q16)),
    }


def _load_one_npz(path):
    with np.load(path) as d:
        return (
            str(path.parent.name),
            d['events'].copy(),
            d['signals'].copy(),
            d['seeds'].copy(),
            d['tracks'].copy(),
            [str(x) for x in d['event_cols']],
            [str(x) for x in d['signal_cols']],
            [str(x) for x in d['seed_cols']],
            [str(x) for x in d['track_cols']],
        )


def load_all_data(root, groups, workers=8):
    tasks = []
    for g in groups:
        for s in RAW_SAMPLES:
            sdir = root / 'runs' / g / s
            if not sdir.exists():
                continue
            shards = sorted(sdir.glob('*/shard_data.npz'))
            for sp in shards:
                tasks.append((g, s, sp))

    def loader(item):
        g, s, sp = item
        return g, s, sp.parent, _load_one_npz(sp)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        loaded = list(pool.map(loader, tasks))

    by_gs = {}
    cols_meta = None
    for g, s, pdir, tup in loaded:
        sid, ev, sig, sd, trk, ec, sc, sdc, tc = tup
        if cols_meta is None:
            cols_meta = {
                'event': {n: i for i, n in enumerate(ec)},
                'signal': {n: i for i, n in enumerate(sc)},
                'seed': {n: i for i, n in enumerate(sdc)},
                'track': {n: i for i, n in enumerate(tc)},
            }
        by_gs.setdefault((g, s), []).append((sid, pdir, ev, sig, sd, trk))

    merged = {}
    for (g, s), items in by_gs.items():
        items.sort(key=lambda x: x[0])
        merged[(g, s)] = {
            'paths': [x[1] for x in items],
            'events': np.vstack([x[2] for x in items]),
            'signals': np.vstack([x[3] for x in items]) if any(len(x[3]) for x in items) else np.zeros((0, len(cols_meta['signal']))),
            'seeds': np.vstack([x[4] for x in items]) if any(len(x[4]) for x in items) else np.zeros((0, len(cols_meta['seed']))),
            'tracks': np.vstack([x[5] for x in items]) if any(len(x[5]) for x in items) else np.zeros((0, len(cols_meta['track']))),
        }
        if s == '110051':
            # Strictly exclude calibration file00 (entries < 100000) for held-out 19-file evaluation
            items19 = [x for x in items if x[0] != 'file00']
            merged[(g, '110051_heldout19')] = {
                'paths': [x[1] for x in items19],
                'events': np.vstack([x[2] for x in items19]),
                'signals': np.vstack([x[3] for x in items19]),
                'seeds': np.vstack([x[4] for x in items19]),
                'tracks': np.vstack([x[5] for x in items19]),
            }
    return merged, cols_meta


def summarize_group_sample(data, cols):
    ec, sc, sdc, tc = cols['event'], cols['signal'], cols['seed'], cols['track']
    ev, sig, sd, trk = data['events'], data['signals'], data['seeds'], data['tracks']

    # Event & combinatorics summary
    n_events = len(ev)
    raw_comb = int(np.sum(ev[:, ec['raw_comb']]))
    raw_nan = int(np.sum(ev[:, ec['raw_nan_chi2']]))
    n_seeds = len(sd)
    sel_nan = int(np.sum(ev[:, ec['n_seeds_nan_chi2']]))

    # Signal efficiency & 2-track pair efficiency
    fid_mask = sig[:, sc['fiducial']] > 0.5
    fid_sig = sig[fid_mask]
    n_fid = len(fid_sig)
    fid_seed_eff = float(np.mean(fid_sig[:, sc['n_seed_match']] >= 1)) if n_fid else None
    fid_ckf_eff = float(np.mean(fid_sig[:, sc['n_ckf_match']] >= 1)) if n_fid else None
    fid_dup_rate = float(np.mean(fid_sig[:, sc['n_ckf_match']] >= 2)) if n_fid else None

    # Pair efficiency on events that have >=2 fiducial truth particles
    pair_eff = None
    n_pair_events = 0
    if n_fid > 0:
        f_entries = fid_sig[:, sc['entry']].astype(np.int64)
        f_matched = (fid_sig[:, sc['n_ckf_match']] >= 1).astype(np.int32)
        uniq_ent, inv, cnts = np.unique(f_entries, return_inverse=True, return_counts=True)
        matched_per_ent = np.bincount(inv, weights=f_matched)
        two_mask = cnts >= 2
        n_pair_events = int(np.sum(two_mask))
        if n_pair_events > 0:
            pair_eff = float(np.mean(matched_per_ent[two_mask] >= 2))

    # Seed population composition & residuals
    valid_seed_p = np.isfinite(sd[:, sdc['p_prod_GeV']]) & (sd[:, sdc['p_prod_GeV']] > 0)
    sd_v = sd[valid_seed_p]
    p_seed_res_prod = (sd_v[:, sdc['p_seed_GeV']] - sd_v[:, sdc['p_prod_GeV']]) / sd_v[:, sdc['p_prod_GeV']]
    valid_st1 = np.isfinite(sd_v[:, sdc['p_st1_GeV']]) & (sd_v[:, sdc['p_st1_GeV']] > 0)
    p_seed_res_st1 = (sd_v[valid_st1, sdc['p_seed_GeV']] - sd_v[valid_st1, sdc['p_st1_GeV']]) / sd_v[valid_st1, sdc['p_st1_GeV']]

    clean_3st = (
        (sd_v[:, sdc['all_same_bc']] > 0.5)
        & (sd_v[:, sdc['stations']] == 3)
        & (sd_v[:, sdc['purity']] >= 0.7)
        & ((sd_v[:, sdc['p_pre_s2_GeV']] / np.maximum(sd_v[:, sdc['p_st1_GeV']], 1e-6)) >= 0.95)
    )
    opp_mix = sd_v[:, sdc['opp_primary_mix']] > 0.5
    clean_seed_res_st1 = (sd_v[clean_3st, sdc['p_seed_GeV']] - sd_v[clean_3st, sdc['p_st1_GeV']]) / sd_v[clean_3st, sdc['p_st1_GeV']]

    # Seed angle residuals (clean 3-station seeds with valid S1 hit truth)
    ang_mask = clean_3st & (sd_v[:, sdc['seed_truth_valid']] > 0.5)
    dphi_mrad = wrap_pi(sd_v[ang_mask, sdc['phi_A']] - sd_v[ang_mask, sdc['t_phi_A']]) * 1e3
    dtheta_mrad = (sd_v[ang_mask, sdc['theta_A']] - sd_v[ang_mask, sdc['t_theta_A']]) * 1e3
    q_ang = sd_v[ang_mask, sdc['truth_q']]
    dphi_neg = dphi_mrad[q_ang < 0]
    dphi_pos = dphi_mrad[q_ang > 0]

    # Strip-fit chi2/ndf separation on real vs cross-track seeds
    has_strip = np.isfinite(sd_v[:, sdc['chi2_strip']]) & (sd_v[:, sdc['ndf_strip']] > 0)
    strip_red_clean = sd_v[clean_3st & has_strip, sdc['chi2_strip']] / sd_v[clean_3st & has_strip, sdc['ndf_strip']]
    strip_red_opp = sd_v[opp_mix & has_strip, sdc['chi2_strip']] / sd_v[opp_mix & has_strip, sdc['ndf_strip']]

    # Build good-seed lookup for tracks
    good_seed_keys = set()
    if len(sd):
        good_m = (
            (sd[:, sdc['seed_truth_valid']] > 0.5)
            & np.isfinite(sd[:, sdc['p_prod_GeV']]) & (sd[:, sdc['p_prod_GeV']] > 0)
            & (np.abs(sd[:, sdc['p_seed_GeV']] / sd[:, sdc['p_prod_GeV']] - 1.0) < 0.1)
            & (sd[:, sdc['q_seed']] * sd[:, sdc['truth_q']] > 0.0)
            & (np.abs(sd[:, sdc['theta_A']] - sd[:, sdc['t_theta_A']]) < 1e-4)
        )
        for r in sd[good_m]:
            good_seed_keys.add((int(round(r[sdc['entry']])), int(round(r[sdc['iseed']]))))

    # Track selection: GitHub quality cuts + purity > 0.7
    gh_mask = (
        (trk[:, tc['ndf_ckf']] > 9.0)
        & ((trk[:, tc['chi2_ckf']] / np.maximum(trk[:, tc['ndf_ckf']], 1.0)) < 20.0)
        & (trk[:, tc['rmax_ckf']] < 90.0)
    ) if len(trk) else np.zeros(0, dtype=bool)
    pur_mask = gh_mask & (trk[:, tc['trk_purity']] >= 0.7)
    fake_gh_fraction = float(np.mean(trk[gh_mask, tc['trk_purity']] < 0.7)) if np.any(gh_mask) else 0.0

    same_mask = pur_mask & (trk[:, tc['seed_maj_bc']] == trk[:, tc['trk_maj_bc']]) & (trk[:, tc['trk_maj_bc']] > 0)
    hijack_rate_pur = float(np.mean(trk[pur_mask, tc['seed_maj_bc']] != trk[pur_mask, tc['trk_maj_bc']])) if np.any(pur_mask) else 0.0

    trk_good_flag = np.asarray(
        [(int(round(r[tc['entry']])), int(round(r[tc['iseed']]))) in good_seed_keys for r in trk],
        dtype=bool,
    ) if len(trk) else np.zeros(0, dtype=bool)
    good_all_mask = pur_mask & trk_good_flag
    good_same_mask = same_mask & trk_good_flag
    hijack_rate_good = float(np.mean(trk[good_all_mask, tc['seed_maj_bc']] != trk[good_all_mask, tc['trk_maj_bc']])) if np.any(good_all_mask) else 0.0

    def p_res_from_qop(qop_arr, p_ref):
        with np.errstate(divide='ignore', invalid='ignore'):
            p_rec = np.where(np.isfinite(qop_arr) & (np.abs(qop_arr) > 0.0), 1.0 / np.abs(qop_arr), np.nan)
            return np.where(np.isfinite(p_ref) & (p_ref > 0.0), (p_rec - p_ref) / p_ref, np.nan)

    def cohort_metrics(mask):
        sub = trk[mask]
        if len(sub) == 0:
            return {'count': 0}
        p_prod = sub[:, tc['p_prod_GeV']]
        p_st1 = sub[:, tc['p_st1_GeV']]
        p_pre_s2 = sub[:, tc['p_pre_s2_GeV']]
        p_min = sub[:, tc['p_min_trk_GeV']]
        tq = sub[:, tc['truth_q']]

        ref_ok = sub[:, tc['refit_ok']] > 0.5
        final_qop = np.where(ref_ok & np.isfinite(sub[:, tc['ref_w0_qop']]), sub[:, tc['ref_w0_qop']], sub[:, tc['ckf_w0_qop']])
        wrong_q_ckf = float(np.mean(sub[:, tc['ckf_w0_qop']] * tq < 0.0))
        wrong_q_final = float(np.mean(final_qop * tq < 0.0))

        # Wafer-0 5-parameter covariance pulls (exact SiHit local truth on w0!)
        w0_ok = (sub[:, tc['hit_w0_valid']] > 0.5) & (sub[:, tc['ckf_w0_geo']] > 0)
        w0_ref_ok = w0_ok & ref_ok & (sub[:, tc['ref_w0_geo']] == sub[:, tc['ckf_w0_geo']])

        def calc_pulls(sub_p, prefix):
            if len(sub_p) == 0:
                return {}
            l0 = (sub_p[:, tc[f'{prefix}_loc0']] - sub_p[:, tc['hit_w0_loc0']]) / sub_p[:, tc[f'{prefix}_s0']]
            l1 = (sub_p[:, tc[f'{prefix}_loc1']] - sub_p[:, tc['hit_w0_loc1']]) / sub_p[:, tc[f'{prefix}_s1']]
            ph = wrap_pi(sub_p[:, tc[f'{prefix}_phi']] - sub_p[:, tc['hit_w0_phi']]) / sub_p[:, tc[f'{prefix}_s2']]
            th = (sub_p[:, tc[f'{prefix}_theta']] - sub_p[:, tc['hit_w0_theta']]) / sub_p[:, tc[f'{prefix}_s3']]
            qp_prod = (sub_p[:, tc[f'{prefix}_qop']] - sub_p[:, tc['truth_q']] / sub_p[:, tc['p_prod_GeV']]) / sub_p[:, tc[f'{prefix}_s4']]
            qp_st1 = (sub_p[:, tc[f'{prefix}_qop']] - sub_p[:, tc['truth_q']] / sub_p[:, tc['p_st1_GeV']]) / sub_p[:, tc[f'{prefix}_s4']]
            # Non-brem subset (p_min / p_st1 >= 0.995) for clean Gaussian pull calibration check
            nb = (sub_p[:, tc['p_min_trk_GeV']] / np.maximum(sub_p[:, tc['p_st1_GeV']], 1e-6)) >= 0.995
            return {
                'loc0': pull_stats(l0),
                'loc1': pull_stats(l1),
                'phi': pull_stats(ph),
                'theta': pull_stats(th),
                'qop_vs_prod': pull_stats(qp_prod),
                'qop_vs_st1': pull_stats(qp_st1),
                'qop_nobrem_vs_st1': pull_stats(qp_st1[nb]),
                'phi_nobrem': pull_stats(ph[nb]),
            }

        nobrem = (p_min / np.maximum(p_st1, 1e-6)) >= 0.995
        hardbrem = (p_min / np.maximum(p_st1, 1e-6)) < 0.95

        return {
            'count': int(len(sub)),
            'refit_success_rate': float(np.mean(ref_ok)),
            'wrong_charge_ckf': wrong_q_ckf,
            'wrong_charge_final': wrong_q_final,
            'chi2_ndf_ckf': robust_stats(sub[:, tc['chi2_ckf']] / np.maximum(sub[:, tc['ndf_ckf']], 1.0)),
            'chi2_ndf_refit': robust_stats(sub[ref_ok, tc['chi2_refit']] / np.maximum(sub[ref_ok, tc['ndf_refit']], 1.0)),
            's2_first_pull': pull_stats(sub[:, tc['s2_first_pull']]),
            'vs_prod': {
                's2_first_prt': robust_stats(p_res_from_qop(sub[:, tc['s2_first_prt_qop']], p_prod)),
                's2_first_flt': robust_stats(p_res_from_qop(sub[:, tc['s2_first_flt_qop']], p_prod)),
                's2_first_smt': robust_stats(p_res_from_qop(sub[:, tc['s2_first_smt_qop']], p_prod)),
                'ckf_w0_smt': robust_stats(p_res_from_qop(sub[:, tc['ckf_w0_qop']], p_prod)),
                'final_refit': robust_stats(p_res_from_qop(final_qop, p_prod)),
            },
            'vs_st1_local': {
                's2_first_prt': robust_stats(p_res_from_qop(sub[:, tc['s2_first_prt_qop']], p_st1)),
                's2_first_flt': robust_stats(p_res_from_qop(sub[:, tc['s2_first_flt_qop']], p_st1)),
                's2_first_smt': robust_stats(p_res_from_qop(sub[:, tc['s2_first_smt_qop']], p_st1)),
                'ckf_w0_smt': robust_stats(p_res_from_qop(sub[:, tc['ckf_w0_qop']], p_st1)),
                'final_refit': robust_stats(p_res_from_qop(final_qop, p_st1)),
                'final_nobrem': robust_stats(p_res_from_qop(final_qop[nobrem], p_st1[nobrem])),
                'final_hardbrem_vs_st1': robust_stats(p_res_from_qop(final_qop[hardbrem], p_st1[hardbrem])),
                'final_hardbrem_vs_pres2': robust_stats(p_res_from_qop(final_qop[hardbrem], p_pre_s2[hardbrem])),
            },
            'pulls_ckf_w0': calc_pulls(sub[w0_ok], 'ckf_w0'),
            'pulls_refit_w0': calc_pulls(sub[w0_ref_ok], 'ref_w0'),
        }

    return {
        'events': int(n_events),
        'combinatorics': {
            'raw_combinations': raw_comb,
            'raw_nan_chi2': raw_nan,
            'selected_seeds': int(n_seeds),
            'selected_nan_chi2': sel_nan,
            'seeds_per_event': float(n_seeds / max(n_events, 1)),
        },
        'efficiency': {
            'fiducial_signals': int(n_fid),
            'seed_efficiency': fid_seed_eff,
            'ckf_efficiency': fid_ckf_eff,
            'duplicate_rate': fid_dup_rate,
            'two_track_pair_events': n_pair_events,
            'two_track_pair_efficiency': pair_eff,
            'fake_track_fraction_in_github_cuts': fake_gh_fraction,
        },
        'seed_distributions': {
            'valid_truth_seeds': int(len(sd_v)),
            'clean_3st_fraction': float(np.mean(clean_3st)) if len(sd_v) else 0.0,
            'cross_opp_primary_fraction': float(np.mean(opp_mix)) if len(sd_v) else 0.0,
            'shared_exact6_count': int(np.sum(sd[:, sdc['max_shared_prev']] == 6)) if len(sd) else 0,
            'shared_union_gt6_count': int(np.sum(sd[:, sdc['shared_union_prev']] > 6)) if len(sd) else 0,
            'all_seed_p_res_vs_prod': robust_stats(p_seed_res_prod),
            'all_seed_p_res_vs_st1': robust_stats(p_seed_res_st1),
            'clean_3st_p_res_vs_st1': robust_stats(clean_seed_res_st1),
            'clean_dphi_all_mrad': robust_stats(dphi_mrad),
            'clean_dphi_qneg_mrad': robust_stats(dphi_neg),
            'clean_dphi_qpos_mrad': robust_stats(dphi_pos),
            'clean_dtheta_mrad': robust_stats(dtheta_mrad),
            'strip_red_chi2_clean': robust_stats(strip_red_clean),
            'strip_red_chi2_cross_opp': robust_stats(strip_red_opp),
        },
        'track_counts': {
            'writer_tracks': int(len(trk)),
            'github_cuts': int(np.sum(gh_mask)),
            'purity_ge_07': int(np.sum(pur_mask)),
            'same_seed_track_truth': int(np.sum(same_mask)),
            'hijack_rate_in_purity07': hijack_rate_pur,
            'good_seed_all': int(np.sum(good_all_mask)),
            'good_seed_same_truth': int(np.sum(good_same_mask)),
            'hijack_rate_in_good_seed': hijack_rate_good,
        },
        'ungated_purity07': cohort_metrics(pur_mask),
        'ungated_same_truth': cohort_metrics(same_mask),
        'good_seed_same_truth': cohort_metrics(good_same_mask),
    }


def export_figures(root, output, merged, cols, summaries, field_closure):
    sdc, tc = cols['seed'], cols['track']
    bins_store = {}
    hist_meta = {}
    edges_wide = np.linspace(-1.2, 1.2, 101)
    edges_core = np.linspace(-0.5, 0.5, 81)

    def rec_hist(key, values, edges):
        arr = np.asarray(values, dtype=np.float64)
        fin = arr[np.isfinite(arr)]
        cnt, _ = np.histogram(fin, bins=edges)
        bins_store[f'{key}_counts'] = cnt
        bins_store[f'{key}_edges'] = edges
        hist_meta[key] = {
            'raw': int(len(arr)),
            'finite': int(len(fin)),
            'underflow': int(np.sum(fin < edges[0])),
            'overflow': int(np.sum(fin >= edges[-1])),
        }
        return cnt, len(fin)

    def save_fig(fig, stem, pdf_handle):
        fig.savefig(output / f'{stem}.png', dpi=170)
        fig.savefig(output / f'{stem}.pdf')
        pdf_handle.savefig(fig)
        plt.close(fig)

    with PdfPages(output / 'diagnostics.pdf') as pdf:
        # FIGURE 1: Multi-sample B_eff deterministic closure + RDO clean seed momentum calibration across 5 species
        fig, axes = plt.subplots(1, 3, figsize=(17.4, 5.4), layout='constrained')
        smp_labels = [
            ('110051', "110051 A'->e+e- (1f, 250 ev)", '#1f77b4'),
            ('110033', "110033 A'->e+e- (1f, 250 ev)", '#17becf'),
            ('100116', '100116 PG mu+ 2D (1f, 250 ev)', '#2ca02c'),
            ('100117', '100117 PG mu- 2D (1f, 250 ev)', '#ff7f0e'),
            ('100130', '100130 PG K_S hadrons (1f, 250 ev)', '#9467bd'),
        ]
        b_edges = np.linspace(0.465, 0.482, 69)
        for smp, lbl, col in smp_labels:
            fc_file = root / 'field_closure' / smp / 'field_closure.json'
            if fc_file.exists():
                pts = json.loads(fc_file.read_text())['particles']
                beff = np.array([p['effective_B_T'] for p in pts])
                rmax = np.array([p['r_max_mm'] for p in pts])
                cnt, nf = rec_hist(f'fc_beff_{smp}', beff, b_edges)
                axes[0].stairs(cnt / max(nf, 1), b_edges, color=col, lw=1.4,
                               label=f'{lbl}: N={nf}, med={np.median(beff):.5f} T')
                rbins = np.linspace(0, 100, 6)
                xs_b, ys_b = [], []
                for bi in range(len(rbins) - 1):
                    m = (rmax >= rbins[bi]) & (rmax < rbins[bi + 1])
                    if np.sum(m) >= 5:
                        xs_b.append(0.5 * (rbins[bi] + rbins[bi + 1]))
                        ys_b.append(np.median(beff[m]))
                if xs_b:
                    axes[1].plot(xs_b, ys_b, 'o-', color=col, lw=1.3, ms=4.5, label=lbl)

        axes[0].axvline(0.473662, color='black', ls='--', lw=1.1, label='Calibrated B = 0.473662 T')
        axes[0].set(xlabel='Effective 3-station circle field B_eff [T]', ylabel='Normalized fraction / bin',
                    title='Full-field CircleFit B_eff (entries [0,250), FASER z-y circle)')
        axes[0].legend(fontsize=6.4, loc='upper right')
        axes[1].axhline(0.473662, color='black', ls='--', lw=1.0, label='0.473662 T constant')
        axes[1].set(xlim=(0, 100), ylim=(0.468, 0.480), xlabel='Max transverse radius r_max in S1-S3 [mm]',
                    ylabel='Median B_eff [T]', title='Transverse bore uniformity (r < 100 mm, S1-S3)')
        axes[1].grid(alpha=0.25)
        axes[1].legend(fontsize=6.6)

        # Panel C: Actual RDO clean 3-station seed residual (vs pre-S1 local truth) for B=0.55 T (C2) vs B=0.4737 T (C3)
        for smp, lbl, col in [
            ('110051_heldout19', "110051 e+e- (19f, 19k)", '#1f77b4'),
            ('100116', '100116 mu+ (1f [0,10k))', '#2ca02c'),
            ('100130', '100130 K_S hadrons (1f [0,10k))', '#9467bd'),
        ]:
            for gname, ls, tag in [('C2_orig_var', '--', '0.55T'), ('C3_seed_b', '-', '0.4737T')]:
                sd = merged[(gname, smp)]['seeds']
                cl = (
                    (sd[:, sdc['all_same_bc']] > 0.5) & (sd[:, sdc['stations']] == 3)
                    & (sd[:, sdc['purity']] >= 0.7)
                    & np.isfinite(sd[:, sdc['p_st1_GeV']]) & (sd[:, sdc['p_st1_GeV']] > 0)
                    & ((sd[:, sdc['p_pre_s2_GeV']] / np.maximum(sd[:, sdc['p_st1_GeV']], 1e-6)) >= 0.95)
                )
                vals = (sd[cl, sdc['p_seed_GeV']] - sd[cl, sdc['p_st1_GeV']]) / sd[cl, sdc['p_st1_GeV']]
                cnt, nf = rec_hist(f'rdo_seed_{smp}_{gname}', vals, edges_core)
                axes[2].stairs(cnt / max(nf, 1), edges_core, color=col, ls=ls, lw=1.35,
                               label=f'{lbl} [{tag}] (med={np.nanmedian(vals)*100:+.2f}%)')
        axes[2].axvline(0.0, color='black', ls=':', lw=0.9)
        axes[2].set(xlim=(-0.45, 0.45), xlabel='Clean 3-st seed (p_seed - p_st1) / p_st1 [dimensionless]',
                    ylabel='Fraction of seeds / 0.0125 bin',
                    title='RDO seed targetSurface (ACTS A): 0.55 T vs 0.4737 T')
        axes[2].legend(fontsize=6.2, loc='upper right')
        fig.suptitle('SeedBFieldTesla = 0.473662 T universality: deterministic closure (5 samples x 250 ev) & full RDO seeding (ACTS frame A)')
        save_fig(fig, 'seed_b_multisample_closure', pdf)

        # FIGURE 2: Step-by-step momentum distribution evolution C0 -> C1 -> C2 -> C3 -> C4 -> C5 (all 6 steps)
        fig, axes = plt.subplots(1, 3, figsize=(17.6, 5.4), layout='constrained')
        step_styles = [
            ('C0_baseline', 'C0: latest baseline', '#d62728', '-'),
            ('C1_cov_loc', 'C1: + fix covLoc swap', '#ff7f0e', '--'),
            ('C2_orig_var', 'C2: + orig var [1e2..1e3]', '#9467bd', '-.'),
            ('C3_seed_b', 'C3: + Seed B=0.4737 T', '#1f77b4', ':'),
            ('C4_no_nan', 'C4: + getY NaN guard', '#17becf', '--'),
            ('C5_local_tan', 'C5: + local tangent (all 5)', '#2ca02c', '-'),
        ]
        for gname, glbl, col, ls in step_styles:
            # Panel A: 110051 held-out 19f S2 first filtered p residual (ungated same_truth)
            t19 = merged[(gname, '110051_heldout19')]['tracks']
            m19 = (
                (t19[:, tc['ndf_ckf']] > 9) & ((t19[:, tc['chi2_ckf']] / np.maximum(t19[:, tc['ndf_ckf']], 1)) < 20)
                & (t19[:, tc['rmax_ckf']] < 90) & (t19[:, tc['trk_purity']] >= 0.7)
                & (t19[:, tc['seed_maj_bc']] == t19[:, tc['trk_maj_bc']])
            )
            qflt = t19[m19, tc['s2_first_flt_qop']]
            p19 = t19[m19, tc['p_prod_GeV']]
            r_flt = (1.0 / np.abs(qflt) - p19) / p19
            cnt, nf = rec_hist(f'step_s2flt_{gname}', r_flt, edges_wide)
            axes[0].stairs(cnt / max(nf, 1), edges_wide, color=col, ls=ls, lw=1.35,
                           label=f'{glbl} (med={np.nanmedian(r_flt)*100:+.1f}%)')

            # Panel B: 110051 held-out 19f good-seed CKF smoothed p residual
            s19 = merged[(gname, '110051_heldout19')]['seeds']
            gm = (
                (s19[:, sdc['seed_truth_valid']] > 0.5) & (s19[:, sdc['p_prod_GeV']] > 0)
                & (np.abs(s19[:, sdc['p_seed_GeV']] / s19[:, sdc['p_prod_GeV']] - 1.0) < 0.1)
                & (s19[:, sdc['q_seed']] * s19[:, sdc['truth_q']] > 0)
                & (np.abs(s19[:, sdc['theta_A']] - s19[:, sdc['t_theta_A']]) < 1e-4)
            )
            gkeys = {(int(round(r[0])), int(round(r[1]))) for r in s19[gm]}
            tg = np.array([(int(round(r[0])), int(round(r[2]))) in gkeys for r in t19])
            mg = m19 & tg
            qsmt = t19[mg, tc['s2_first_smt_qop']]
            r_smt = (1.0 / np.abs(qsmt) - t19[mg, tc['p_prod_GeV']]) / t19[mg, tc['p_prod_GeV']]
            cnt_g, nf_g = rec_hist(f'step_goodsmt_{gname}', r_smt, edges_core)
            axes[1].stairs(cnt_g / max(nf_g, 1), edges_core, color=col, ls=ls, lw=1.35,
                           label=f'{glbl} (N={nf_g}, med={np.nanmedian(r_smt)*100:+.2f}%)')

            # Panel C: 100116 + 100117 muons ungated final refit p residual vs pre-S1 local truth
            tmu = np.vstack([merged[(gname, '100116')]['tracks'], merged[(gname, '100117')]['tracks']])
            mmu = (
                (tmu[:, tc['ndf_ckf']] > 9) & ((tmu[:, tc['chi2_ckf']] / np.maximum(tmu[:, tc['ndf_ckf']], 1)) < 20)
                & (tmu[:, tc['rmax_ckf']] < 90) & (tmu[:, tc['trk_purity']] >= 0.7)
            )
            rok = tmu[mmu, tc['refit_ok']] > 0.5
            qfin = np.where(rok & np.isfinite(tmu[mmu, tc['ref_w0_qop']]), tmu[mmu, tc['ref_w0_qop']], tmu[mmu, tc['ckf_w0_qop']])
            pst1 = tmu[mmu, tc['p_st1_GeV']]
            r_mu = (1.0 / np.abs(qfin) - pst1) / pst1
            cnt_m, nf_m = rec_hist(f'step_mu_final_{gname}', r_mu, edges_core)
            axes[2].stairs(cnt_m / max(nf_m, 1), edges_core, color=col, ls=ls, lw=1.35,
                           label=f'{glbl} (N={nf_m}, med={np.nanmedian(r_mu)*100:+.2f}%)')

        axes[0].set(xlim=(-1.1, 1.1), xlabel='1st S2 wafer filtered (p - p_prod) / p_prod',
                    ylabel='Fraction of tracks / 0.024 bin',
                    title="110051 (19f [01..19], 19k ev): 1st S2 wafer filtered")
        axes[0].legend(fontsize=6.2)
        axes[1].axvline(0.0, color='gray', ls=':', lw=0.8)
        axes[1].set(xlim=(-0.45, 0.45), xlabel='1st S2 wafer smoothed (p - p_prod) / p_prod',
                    ylabel='Fraction of tracks / 0.0125 bin',
                    title='110051 (19f, 19k ev): good-seed S2 smoothed')
        axes[1].legend(fontsize=6.2)
        axes[2].axvline(0.0, color='gray', ls=':', lw=0.8)
        axes[2].set(xlim=(-0.40, 0.40), xlabel='1st S1 wafer w0 refit (p_final - p_st1) / p_st1',
                    ylabel='Fraction of tracks / 0.0125 bin',
                    title='100116+100117 mu+/- (2f [0,10k)x2, 20k ev): final refit')
        axes[2].legend(fontsize=6.2)
        fig.suptitle('Step-by-step momentum impact C0 -> C1 -> C2 -> C3 -> C4 -> C5 (Rotated ACTS frame A; GeV p from |1/(q/p)|)')
        save_fig(fig, 'ablation_step_by_step_momentum', pdf)

        # FIGURE 3: Kinematic necessity of EVERY individual step C0 -> C1 -> C2 -> C3 -> C4 -> C5
        fig, axes = plt.subplots(1, 3, figsize=(17.4, 5.4), layout='constrained')
        c_steps = ['C0_baseline', 'C1_cov_loc', 'C2_orig_var', 'C3_seed_b', 'C4_no_nan', 'C5_local_tan']
        short_lbls = ['C0\nbaseline', 'C1\n+covLoc', 'C2\n+origVar', 'C3\n+seedB', 'C4\n+noNaN', 'C5\n+locTan']
        x_idx = np.arange(len(c_steps))

        hijack_110051 = [summaries[(g, '110051_heldout19')]['track_counts']['hijack_rate_in_purity07'] * 100.0 for g in c_steps]
        s2_pull_110051 = [summaries[(g, '110051_heldout19')]['ungated_same_truth']['s2_first_pull']['rms_clipped8'] for g in c_steps]
        nan_seeds_110051 = [summaries[(g, '110051')]['combinatorics']['selected_nan_chi2'] for g in c_steps]

        ax0_t = axes[0].twinx()
        l1 = axes[0].plot(x_idx, hijack_110051, 'o-', color='#d62728', lw=1.6, label='Cross-track hijack rate [%] (left)')
        l2 = ax0_t.plot(x_idx, s2_pull_110051, 's--', color='#1f77b4', lw=1.6, label='1st S2 wafer pull RMS (right)')
        ax0_t.axhline(1.0, color='#1f77b4', ls=':', lw=0.8)
        axes[0].set_xticks(x_idx)
        axes[0].set_xticklabels(short_lbls, fontsize=8)
        axes[0].set_ylabel('CKF seed->track cross-hijack rate [%]', color='#d62728')
        ax0_t.set_ylabel('1st S2 wafer innovation pull RMS', color='#1f77b4')
        axes[0].set_ylim(0, max(max(hijack_110051) * 1.25, 1.0))
        ax0_t.set_ylim(0, 1.35)
        axes[0].set_title('C1 (covLoc) & C2 (origVar): 110051 19f (19k ev)')
        axes[0].legend(l1 + l2, [x.get_label() for x in (l1 + l2)], fontsize=6.8, loc='center right')
        axes[0].grid(alpha=0.25)

        # Panel B: Selected NaN-chi2 seeds & raw NaN-chi2 seeds across C0..C5 on 110051 + 110033 + 100130
        for smp, lbl, col in [
            ('110051', "110051 A' (20f, 20k ev)", '#1f77b4'),
            ('110033', "110033 A' (5f, 10k ev)", '#ff7f0e'),
            ('100130', '100130 K_S (1f [0,10k))', '#9467bd'),
        ]:
            vals_nan = [summaries[(g, smp)]['combinatorics']['raw_nan_chi2'] for g in c_steps]
            axes[1].plot(x_idx, vals_nan, 'o-', color=col, lw=1.5, label=f'{lbl} raw NaN chi2')
        axes[1].plot(x_idx, nan_seeds_110051, 's--', color='#d62728', lw=1.5, label='110051 (20f) selected NaN chi2')
        axes[1].set_xticks(x_idx)
        axes[1].set_xticklabels(short_lbls, fontsize=8)
        axes[1].set(ylabel='Seeds with NaN getY()/chi2 [count]', title='C4 (getY NaN guard + strict-weak sort)')
        axes[1].grid(alpha=0.25)
        axes[1].legend(fontsize=6.8)

        # Panel C: Seed bending angle residual delta_phi_A [mrad] for q<0 vs q>0 in C4 (chord) vs C5 (local tangent)
        phi_edges = np.linspace(-0.6, 0.6, 81)
        for gname, ls, gtag in [('C4_no_nan', '--', 'C4 S1->S2 chord'), ('C5_local_tan', '-', 'C5 circle tangent')]:
            sd = merged[(gname, '110051_heldout19')]['seeds']
            cl = (sd[:, sdc['all_same_bc']] > 0.5) & (sd[:, sdc['stations']] == 3) & (sd[:, sdc['purity']] >= 0.7) & (sd[:, sdc['seed_truth_valid']] > 0.5)
            dphi = wrap_pi(sd[cl, sdc['phi_A']] - sd[cl, sdc['t_phi_A']]) * 1e3
            q = sd[cl, sdc['truth_q']]
            for qsign, qlbl, col in [(-1, 'e-', '#1f77b4'), (+1, 'e+', '#d62728')]:
                v = dphi[q * qsign > 0]
                cnt, nf = rec_hist(f'dphi_{gname}_{qlbl}', v, phi_edges)
                axes[2].stairs(cnt / max(nf, 1), phi_edges, color=col, ls=ls, lw=1.35,
                               label=f'{gtag} {qlbl} (med={np.median(v):+.3f}, s68={0.5*(np.quantile(v,0.84)-np.quantile(v,0.16)):.3f} mrad)')
        axes[2].axvline(0.0, color='black', ls=':', lw=0.8)
        axes[2].set(xlim=(-0.55, 0.55), xlabel='Seed targetSurface phi_A - phi_truth,A [mrad] (ACTS frame A)',
                    ylabel='Fraction of clean seeds / 0.015 mrad',
                    title='C5 (local tangent): 110051 19f (19k ev)')
        axes[2].legend(fontsize=6.3)
        fig.suptitle('Observable-specific necessity of C1 (covLoc), C2 (origVar), C4 (no-NaN), and C5 (local tangent) in Rotated ACTS frame A')
        save_fig(fig, 'ablation_kinematics_necessity', pdf)

        # FIGURE 4: Check 1 - ParticleHypothesis::muon() [C5] vs ParticleHypothesis::electron() [H1] under identical material
        fig, axes = plt.subplots(1, 3, figsize=(17.6, 5.4), layout='constrained')
        for smp, slbl, col_m, col_e in [
            ('110051_heldout19', "110051 e+/- (19f, 19k)", '#1f77b4', '#d62728'),
            ('100116', '100116 mu+ (1f [0,10k))', '#2ca02c', '#ff7f0e'),
            ('100130', '100130 K_S hadrons (1f [0,10k))', '#9467bd', '#8c564b'),
        ]:
            for gname, htag, col, ls in [('C5_local_tan', 'muon()', col_m, '-'), ('H1_electron_hypo', 'electron()', col_e, '--')]:
                trk = merged[(gname, smp)]['tracks']
                m = (
                    (trk[:, tc['ndf_ckf']] > 9) & ((trk[:, tc['chi2_ckf']] / np.maximum(trk[:, tc['ndf_ckf']], 1)) < 20)
                    & (trk[:, tc['rmax_ckf']] < 90) & (trk[:, tc['trk_purity']] >= 0.7)
                )
                rok = trk[m, tc['refit_ok']] > 0.5
                qfin = np.where(rok & np.isfinite(trk[m, tc['ref_w0_qop']]), trk[m, tc['ref_w0_qop']], trk[m, tc['ckf_w0_qop']])
                pst1 = trk[m, tc['p_st1_GeV']]
                r_p = (1.0 / np.abs(qfin) - pst1) / pst1
                cnt, nf = rec_hist(f'hypo_pres_{smp}_{gname}', r_p, edges_core)
                axes[0].stairs(cnt / max(nf, 1), edges_core, color=col, ls=ls, lw=1.35,
                               label=f'{slbl} {htag} (med={np.nanmedian(r_p)*100:+.2f}%)')

                c2_edges = np.linspace(0, 8, 65)
                c2_vals = trk[m][rok, tc['chi2_refit']] / np.maximum(trk[m][rok, tc['ndf_refit']], 1.0)
                cnt2, nf2 = rec_hist(f'hypo_chi2_{smp}_{gname}', c2_vals, c2_edges)
                axes[1].stairs(cnt2 / max(nf2, 1), c2_edges, color=col, ls=ls, lw=1.35,
                               label=f'{slbl} {htag} (med={np.nanmedian(c2_vals):.2f})')

                p_edges = np.linspace(-5, 5, 81)
                sub_r = trk[m][rok & (trk[m, tc['hit_w0_valid']] > 0.5) & (trk[m, tc['ref_w0_geo']] == trk[m, tc['ckf_w0_geo']])]
                pull_qop = (sub_r[:, tc['ref_w0_qop']] - sub_r[:, tc['truth_q']] / sub_r[:, tc['p_st1_GeV']]) / sub_r[:, tc['ref_w0_s4']]
                cnt3, nf3 = rec_hist(f'hypo_pull_{smp}_{gname}', pull_qop, p_edges)
                s68 = 0.5 * (np.nanquantile(pull_qop, 0.84) - np.nanquantile(pull_qop, 0.16))
                axes[2].stairs(cnt3 / max(nf3, 1), p_edges, color=col, ls=ls, lw=1.35,
                               label=f'{slbl} {htag} (med={np.nanmedian(pull_qop):+.2f}, s68={s68:.2f})')

        axes[0].axvline(0.0, color='gray', ls=':', lw=0.8)
        axes[0].set(xlim=(-0.45, 0.45), xlabel='1st wafer w0 refit (p_final - p_st1) / p_st1',
                    ylabel='Fraction of tracks / 0.0125 bin', title='Momentum residual: muon() vs electron()')
        axes[0].legend(fontsize=5.8, loc='upper right')
        axes[1].set(xlim=(0, 8), xlabel='Post-refit chi2 / ndf [dimensionless]', ylabel='Fraction of tracks / 0.125 bin',
                    title='Post-refit chi2/ndf (same rotated material map)')
        axes[1].legend(fontsize=6.0)
        axes[2].axvline(0.0, color='gray', ls=':', lw=0.8)
        axes[2].set(xlim=(-5, 5), xlabel='1st wafer w0 pull (q/p - q/p_st1) / sigma(q/p) (ACTS A)',
                    ylabel='Fraction of tracks / 0.125 bin', title='Wafer-0 q/p covariance pull')
        axes[2].legend(fontsize=5.8)
        fig.suptitle('Exploratory Check 1: ParticleHypothesis::muon() vs electron() on identical events & identical material map (ACTS frame A)')
        save_fig(fig, 'particle_hypothesis_comparison', pdf)

        # FIGURE 5: Check 2 - Seed deduplication, Cartesian prefilter, and 18-cluster 2-dipole strip field fit
        fig, axes = plt.subplots(1, 3, figsize=(17.4, 5.4), layout='constrained')
        sd_c5 = merged[('C5_local_tan', '110051_heldout19')]['seeds']
        hs = np.isfinite(sd_c5[:, sdc['chi2_strip']]) & (sd_c5[:, sdc['ndf_strip']] > 0)
        cl_s = hs & (sd_c5[:, sdc['all_same_bc']] > 0.5) & (sd_c5[:, sdc['stations']] == 3) & (sd_c5[:, sdc['purity']] >= 0.7)
        op_s = hs & (sd_c5[:, sdc['opp_primary_mix']] > 0.5)
        log_edges = np.linspace(-1.5, 5.5, 71)
        for m_s, k_s, lbl_s, col_s in [
            (cl_s, 'strip_chi2_clean', 'Clean single-particle 3-st seeds', '#2ca02c'),
            (op_s, 'strip_chi2_cross', 'Cross-track opposite-sign e+/e- seeds', '#d62728'),
        ]:
            lvals = np.log10(np.maximum(sd_c5[m_s, sdc['chi2_strip']] / sd_c5[m_s, sdc['ndf_strip']], 1e-3))
            cnt, nf = rec_hist(k_s, lvals, log_edges)
            axes[0].stairs(cnt / max(nf, 1), log_edges, color=col_s, lw=1.4, label=f'{lbl_s} (N={nf})')
        axes[0].axvline(math.log10(35.0), color='black', ls='--', lw=1.0, label='Cut chi2_strip/ndf = 35')
        axes[0].set(xlabel='log10(18-cluster 2-dipole strip chi2 / ndf)', ylabel='Fraction of seeds / 0.1 dex',
                    title='110051 (19f, 19k ev): 18-cluster strip fit vs fakes')
        axes[0].legend(fontsize=6.6)

        for gname, glbl, col, ls in [
            ('C0_baseline', 'C0: latest baseline', '#d62728', ':'),
            ('C5_local_tan', 'C5: 5 core fixes (circle)', '#ff7f0e', '--'),
            ('S3_strip_only', 'S3: + 18-cl strip fit only', '#9467bd', '-'),
            ('S1_dedup_prefilter', 'S1: + prefilter + S23 dedup', '#1f77b4', '-.'),
            ('S2_strip_field_fit', 'S2: strip + prefilter + dedup', '#2ca02c', '-'),
        ]:
            sd = merged[(gname, '110051_heldout19')]['seeds']
            vs = np.isfinite(sd[:, sdc['p_prod_GeV']]) & (sd[:, sdc['p_prod_GeV']] > 0)
            r_sd = (sd[vs, sdc['p_seed_GeV']] - sd[vs, sdc['p_prod_GeV']]) / sd[vs, sdc['p_prod_GeV']]
            cnt, nf = rec_hist(f'seed_all_{gname}', r_sd, edges_wide)
            spike = np.mean(r_sd < -0.5) * 100.0
            axes[1].stairs(cnt, edges_wide, color=col, ls=ls, lw=1.35,
                           label=f'{glbl} (N={nf}, <-0.5: {spike:.1f}%)')
        axes[1].set(xlim=(-1.1, 1.1), xlabel='Seed targetSurface (p_seed - p_prod) / p_prod',
                    ylabel='Seeds / 0.024 bin', title='110051 (19f, 19k ev): fake-seed spike comparison')
        axes[1].legend(fontsize=6.0, loc='upper right')

        smps_bar = ['110051_heldout19', '110033', '100130', '100116', '100117']
        s_short = ["110051\n19f (19k)", "110033\n5f (10k)", "100130\n1f (10k)", "100116\n1f (10k)", "100117\n1f (10k)"]
        xb = np.arange(len(smps_bar))
        w = 0.20
        for gi, (gname, glbl, col) in enumerate([
            ('C5_local_tan', 'C5 (std circle)', '#ff7f0e'),
            ('S3_strip_only', 'S3 (18-cl strip only)', '#9467bd'),
            ('S1_dedup_prefilter', 'S1 (prefilter+dedup)', '#1f77b4'),
            ('S2_strip_field_fit', 'S2 (strip+pref+dedup)', '#2ca02c'),
        ]):
            effs = [summaries[(gname, sm)]['efficiency']['ckf_efficiency'] * 100.0 for sm in smps_bar]
            axes[2].bar(xb + (gi - 1.5) * w, effs, width=w, color=col, alpha=0.88, label=f'{glbl} eff [%]')
        axes[2].set_xticks(xb)
        axes[2].set_xticklabels(s_short, fontsize=7.6)
        axes[2].set(ylim=(15, 105), ylabel='Fiducial signal CKF track efficiency [%]',
                    title='5-sample fiducial CKF efficiency (ACTS frame A)')
        axes[2].grid(axis='y', alpha=0.25)
        axes[2].legend(fontsize=6.2, loc='lower right')
        fig.suptitle('Exploratory Check 2: Pure 18-cluster strip fit (S3) vs aggressive prefilter + union/S2-S3 deduplication (S1, S2) (ACTS frame A)')
        save_fig(fig, 'seeding_combinatorics_checks', pdf)

        # FIGURE 6: Multi-sample synthesis across all 5 samples and key configurations
        fig, axes = plt.subplots(2, 2, figsize=(16.8, 9.6), layout='constrained')
        syn_groups = [
            ('C0_baseline', 'C0 baseline', '#d62728', 'o-'),
            ('C2_orig_var', 'C2 +covLoc+var', '#9467bd', 's--'),
            ('C3_seed_b', 'C3 +seedB', '#1f77b4', '^:'),
            ('C5_local_tan', 'C5 all 5 fixes', '#2ca02c', 'D-'),
            ('H1_electron_hypo', 'H1 electron()', '#ff7f0e', 'v--'),
            ('S3_strip_only', 'S3 strip fit only', '#8c564b', 'X-'),
            ('S2_strip_field_fit', 'S2 strip+pref+dedup', '#17becf', 'P-.'),
        ]
        syn_samples = [
            ('110051_heldout19', "110051 A' e+e-\n19f / 19k ev"),
            ('110033', "110033 A' e+e-\n5f / 10k ev"),
            ('100116', "100116 PG mu+\n1f [0,10k)"),
            ('100117', "100117 PG mu-\n1f [0,10k)"),
            ('100130', "100130 PG K_S\n1f [0,10k)"),
        ]
        xs_syn = np.arange(len(syn_samples))
        for gname, glbl, col, fmt in syn_groups:
            med_b = [summaries[(gname, sm)]['ungated_purity07']['vs_st1_local']['final_refit']['median'] * 100.0 for sm, _ in syn_samples]
            res68 = [summaries[(gname, sm)]['ungated_purity07']['vs_st1_local']['final_refit']['half68'] * 100.0 for sm, _ in syn_samples]
            pull_s68 = [summaries[(gname, sm)]['ungated_purity07']['pulls_refit_w0']['qop_vs_st1']['sigma68'] for sm, _ in syn_samples]
            fid_eff = [summaries[(gname, sm)]['efficiency']['ckf_efficiency'] * 100.0 for sm, _ in syn_samples]
            axes[0, 0].plot(xs_syn, med_b, fmt, color=col, lw=1.4, ms=5.5, label=glbl)
            axes[0, 1].plot(xs_syn, res68, fmt, color=col, lw=1.4, ms=5.5, label=glbl)
            axes[1, 0].plot(xs_syn, pull_s68, fmt, color=col, lw=1.4, ms=5.5, label=glbl)
            axes[1, 1].plot(xs_syn, fid_eff, fmt, color=col, lw=1.4, ms=5.5, label=glbl)

        axes[0, 0].axhline(0.0, color='black', ls=':', lw=0.9)
        axes[0, 0].set(ylabel='Median (p_final - p_st1) / p_st1 [%]',
                       title='A. Final refit momentum bias vs pre-S1 truth (wafer w0, ACTS A)')
        axes[0, 1].set(ylabel='68% half-width sigma_68(dp/p_st1) [%]',
                       title='B. Final refit momentum resolution sigma_68 (wafer w0, ACTS A)')
        axes[1, 0].axhline(1.0, color='black', ls='--', lw=0.9, label='Unit Gaussian sigma=1')
        axes[1, 0].set(ylabel='Wafer-0 refit q/p pull sigma_68',
                       title='C. Wafer-0 q/p covariance pull width sigma_68 (ACTS A)')
        axes[1, 1].set(ylabel='Fiducial signal CKF track efficiency [%]',
                       title='D. Uniform single-track fiducial CKF efficiency across all 5 samples')
        for ax in axes.ravel():
            ax.set_xticks(xs_syn)
            ax.set_xticklabels([lbl for _, lbl in syn_samples], fontsize=8)
            ax.grid(alpha=0.25)
            ax.legend(fontsize=7.0, ncol=2)
        fig.suptitle('Cross-sample generalization summary (59,000 independent events per group; Rotated ACTS frame A, wafer w0)')
        save_fig(fig, 'multisample_summary_comparison', pdf)

    np.savez_compressed(output / 'bin_counts.npz', **bins_store)
    (output / 'histogram_metadata.json').write_text(json.dumps(hist_meta, indent=2) + '\n')


def export_reference_nine_layers(root, output, merged):
    import ROOT
    ROOT.gROOT.SetBatch(True)

    specs = [
        ('C0_baseline_all20_ungated', merged[('C0_baseline', '110051')]['paths'], 'native_counts_all.root',
         'C0 baseline', '20 files / 20000 events | GitHub + purity>0.7 (ungated)'),
        ('C0_baseline_19f_ungated', merged[('C0_baseline', '110051_heldout19')]['paths'], 'native_counts_all.root',
         'C0 base (19f)', '19 held-out files / 19000 events | GitHub + purity>0.7'),
        ('C1_cov_loc_19f_ungated', merged[('C1_cov_loc', '110051_heldout19')]['paths'], 'native_counts_all.root',
         'C1 +covLoc', '19 held-out files / 19000 events | GitHub + purity>0.7'),
        ('C2_orig_var_19f_ungated', merged[('C2_orig_var', '110051_heldout19')]['paths'], 'native_counts_all.root',
         'C2 +origVar', '19 held-out files / 19000 events | GitHub + purity>0.7'),
        ('C3_seed_b_19f_ungated', merged[('C3_seed_b', '110051_heldout19')]['paths'], 'native_counts_all.root',
         'C3 +seedB', '19 held-out files / 19000 events | GitHub + purity>0.7'),
        ('C4_no_nan_19f_ungated', merged[('C4_no_nan', '110051_heldout19')]['paths'], 'native_counts_all.root',
         'C4 +noNaN', '19 held-out files / 19000 events | GitHub + purity>0.7'),
        ('C5_local_tan_19f_ungated', merged[('C5_local_tan', '110051_heldout19')]['paths'], 'native_counts_all.root',
         'C5 all 5 fixes', '19 held-out files / 19000 events | GitHub + purity>0.7'),
        ('C0_baseline_19f_goodseed', merged[('C0_baseline', '110051_heldout19')]['paths'], 'native_counts_good.root',
         'C0 base good', '19 held-out files / 19000 events | GitHub + good seed + purity>0.7'),
        ('C5_local_tan_19f_goodseed', merged[('C5_local_tan', '110051_heldout19')]['paths'], 'native_counts_good.root',
         'C5 5-fix good', '19 held-out files / 19000 events | GitHub + good seed + purity>0.7'),
        ('H1_electron_19f_goodseed', merged[('H1_electron_hypo', '110051_heldout19')]['paths'], 'native_counts_good.root',
         'H1 e-hypo good', '19 held-out files / 19000 events | electron() + good seed + purity>0.7'),
        ('S3_strip_only_19f_ungated', merged[('S3_strip_only', '110051_heldout19')]['paths'], 'native_counts_all.root',
         'S3 strip only', '19 held-out files / 19000 events | S3 strip fit + GitHub + purity>0.7'),
        ('S2_strip_fit_19f_ungated', merged[('S2_strip_field_fit', '110051_heldout19')]['paths'], 'native_counts_all.root',
         'S2 strip+dedup', '19 held-out files / 19000 events | S2 seeding + GitHub + purity>0.7'),
    ]

    names = [f'h_{c}_{s}_{l}' for c in range(3) for s in (1, 2, 3) for l in (0, 1, 2)]
    ref_pdfs = []
    ref_bins = {}
    raw_entries = {}

    for key, paths, root_name, leg_lbl, scope_txt in specs:
        merged_h, exp_entries, exp_counts = {}, {}, {}
        for p in paths:
            sf = ROOT.TFile.Open(str(p / root_name))
            for nm in names:
                sh = sf.Get(nm)
                assert sh and sh.GetNbinsX() == 100
                vals = np.asarray([sh.GetBinContent(i) for i in range(102)])
                if nm not in merged_h:
                    cl = sh.Clone(nm)
                    cl.SetDirectory(0)
                    merged_h[nm] = cl
                    exp_entries[nm] = float(sh.GetEntries())
                    exp_counts[nm] = vals
                else:
                    merged_h[nm].Add(sh)
                    exp_entries[nm] += float(sh.GetEntries())
                    exp_counts[nm] += vals
            sf.Close()

        cpath = output / f'{key}_counts_3station_forward.root'
        tf = ROOT.TFile.Open(str(cpath), 'RECREATE')
        for nm in names:
            merged_h[nm].Write()
        tf.Close()

        prefix = output / f'110051_{key}_nine_layers'
        header = f'110051 | {scope_txt} | ACTS A | HepMC p'
        macro = root / 'reference/ckf_residual_perLayer.C'
        env = {**os.environ, 'REFERENCE_SCOPE_LABEL': header}
        inv = f'{macro}({",".join(json.dumps(x) for x in [str(cpath), str(prefix), leg_lbl, "t_eQOP"])},true)'
        subprocess.run(['root', '-l', '-b', '-q', inv], env=env, check=True)

        drawn = ROOT.TFile.Open(f'{prefix}.root')
        raw_entries[key] = {}
        for nm in names:
            h = drawn.Get(nm)
            assert h
            vals = np.asarray([h.GetBinContent(i) for i in range(102)])
            assert float(h.GetEntries()) == exp_entries[nm]
            assert np.array_equal(vals, exp_counts[nm])
            ref_bins[f'{key}_{nm}'] = vals
            raw_entries[key][nm] = float(h.GetEntries())
        drawn.Close()
        ref_pdfs.append(f'{prefix}.pdf')

    np.savez_compressed(output / 'reference_bin_counts.npz', **ref_bins)
    return ref_pdfs, raw_entries


def validate_and_compile_pdf(root, output, ref_pdfs, raw_entries):
    import sys
    sys.path.insert(0, str(root.parent.parent / 'stage04_seed_purity/20261113_seed_good_v2/software_tools'))
    from PIL import Image
    from pypdf import PdfReader, PdfWriter

    writer = PdfWriter()
    for p in [output / 'diagnostics.pdf', *map(Path, ref_pdfs)]:
        for page in PdfReader(str(p)).pages:
            writer.add_page(page)
    with (output / 'report.pdf').open('wb') as f:
        writer.write(f)

    fig_idx = []
    for png in sorted(output.glob('*.png')):
        with Image.open(png) as img:
            w, h = img.size
            arr = np.asarray(img.convert('L'))
        assert w >= 900 and h >= 500
        assert float(arr.std()) > 5.0
        pdf_p = png.with_suffix('.pdf')
        assert pdf_p.exists() and len(PdfReader(str(pdf_p)).pages) == 1
        fig_idx.append({'png': png.name, 'pdf': pdf_p.name, 'pixels': [w, h]})

    n_pages = len(PdfReader(str(output / 'report.pdf')).pages)
    val = {
        'figure_index': fig_idx,
        'pdf_pages': n_pages,
        'reference_figures': len(ref_pdfs),
        'reference_histograms': len(ref_pdfs) * 27,
        'reference_bins_with_flow': len(ref_pdfs) * 27 * 102,
        'reference_histograms_exact': True,
        'reference_raw_entries': raw_entries,
    }
    (output / 'figure_validation.json').write_text(json.dumps(val, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--groups', nargs='+', default=ALL_GROUPS)
    parser.add_argument('--workers', type=int, default=max(1, min(16, os.cpu_count() or 4)))
    args = parser.parse_args()

    root = args.root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    merged, cols = load_all_data(root, args.groups, workers=args.workers)

    pairs = sorted(merged.keys())
    def sum_one(gs):
        return gs, summarize_group_sample(merged[gs], cols)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        sum_items = list(pool.map(sum_one, pairs))
    summaries = dict(sum_items)

    fc_path = root / 'field_closure_summary.json'
    field_closure = json.loads(fc_path.read_text()) if fc_path.exists() else {}

    nested_summary = {'field_closure_multisample': field_closure, 'groups': {}}
    for (g, s), val in summaries.items():
        nested_summary['groups'].setdefault(g, {})[s] = val
    (output / 'summary.json').write_text(json.dumps(nested_summary, indent=2) + '\n')

    export_figures(root, output, merged, cols, summaries, field_closure)
    ref_pdfs, raw_entries = export_reference_nine_layers(root, output, merged)
    validate_and_compile_pdf(root, output, ref_pdfs, raw_entries)
    print(f'Completed study in {output}: {len(summaries)} (group, sample) cells')


if __name__ == '__main__':
    main()
