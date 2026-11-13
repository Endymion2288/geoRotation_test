#!/usr/bin/env python3
"""Measure charge response on complete native Calypso before/after samples."""
import argparse
import json
from pathlib import Path

import numpy as np

from charge_history import PBINS, rate
from charge_tables import STATUS, derive, measurement, truth_balance, weighted_rate
from slope_tables import histogram, merge_tables, write_json


DSIDS = ['100116', '100117']
VARIANTS = ['beforeGeo', 'afterGeo']


def paired_difference(a, b, mask):
    """Ratio difference with the covariance of the same source truth events."""
    selected = [d['selected'][mask].astype(float) for d in [a, b]]
    wrong = [d['wrong'][mask].astype(float) for d in [a, b]]
    totals = [s.sum() for s in selected]
    if not all(totals):
        return dict(value=None, standard_error=None, interval95=[None, None])
    values = [y.sum()/n for y, n in zip(wrong, totals)]
    influence = [(y-p*s)/n for y, p, s, n in zip(wrong, values, selected, totals)]
    delta = influence[1]-influence[0]
    error = float(np.sqrt(delta@delta*len(delta)/max(1, len(delta)-1)))
    value = float(values[1]-values[0])
    return dict(value=value, standard_error=error,
                interval68=[value-error, value+error],
                interval95=[value-1.96*error, value+1.96*error],
                source_truth=int(mask.sum()), paired_source_event_covariance=True)


def histograms(tables, d):
    t = tables['tracks']
    rows = d['chosen'][d['selected']]
    p = np.linalg.norm(t['upstream_F'][rows, 3:], axis=1)
    truth_p = d['truth_p'][d['selected']]
    qop = t['qop_GeVinv'][rows]
    wrong = d['wrong'][d['selected']]
    sigma = np.sqrt(np.maximum(t['qop_variance_GeVinv2'][rows], 0))
    truth_qop = d['charge']/truth_p
    with np.errstate(divide='ignore', invalid='ignore'):
        pull = (qop-truth_qop)/sigma
        relative = (qop-truth_qop)/np.abs(truth_qop)
    result = {}
    for name, values, edges in [
        ('signed_qop', qop, np.linspace(-.03, .03, 301)),
        ('wrong_signed_qop', qop[wrong], np.linspace(-.03, .03, 301)),
        ('qop_residual', relative, np.linspace(-6, 6, 241)),
        ('qop_pull', pull, np.linspace(-20, 20, 201)),
        ('reco_p', p, np.linspace(0, 5000, 201)),
        ('truth_p', truth_p, np.linspace(0, 5000, 201)),
        ('relative_p', p/truth_p-1, np.linspace(-1.2, 1.2, 101)),
    ]:
        result[name] = histogram(values, edges)
    edges = np.linspace(0, 2000, 101)
    counts, xe, ye = np.histogram2d(truth_p, p, bins=[edges, edges])
    result['truth_reco_p'] = dict(counts=counts.astype(int).tolist(),
        truth_edges=xe.tolist(), reco_edges=ye.tolist(), selected=len(p),
        outside=int(len(p)-counts.sum()), sumw2=counts.astype(int).tolist())
    return result


def summarize(tables, d):
    t, h, e = tables['tracks'], tables['truth'], tables['events']
    result = measurement(tables, d)
    result['histograms'] = histograms(tables, d)
    result['truth_p_bins'] = []
    for lo, hi in zip(PBINS[:-1], PBINS[1:]):
        mask = (d['truth_p'] >= lo) & (d['truth_p'] < hi)
        result['truth_p_bins'].append(dict(low_GeV=float(lo), high_GeV=float(hi),
            **measurement(tables, d, mask)))
    result['truth_p_outside_bins'] = int(np.count_nonzero(
        d['fiducial'] & ((d['truth_p'] < PBINS[0]) | (d['truth_p'] >= PBINS[-1]))))
    good = d['track_quality']
    unmatched = ~d['track_matched']
    result['quality_tracks'] = dict(total=int(good.sum()),
        unmatched_primary=int(np.count_nonzero(good & unmatched)),
        low_purity=int(np.count_nonzero(good & (t['purity'] <= .7))),
        invalid_curvature=int(np.count_nonzero(good & ~d['track_valid'])))
    result['source_events'] = len(e['sample_entry'])
    result['observed_final_tracks'] = len(t['source_entry'])
    present = d['chosen'] >= 0
    pureless = np.zeros(len(h['sample_entry']), dtype=bool)
    pureless[present] = d['fiducial'][present] & good[d['chosen'][present]] & \
        d['track_valid'][d['chosen'][present]]
    wrong = np.zeros_like(pureless)
    wrong[present] = t['charge'][d['chosen'][present]]*d['charge'] < 0
    result['no_purity_cut'] = dict(
        inclusive=rate(np.count_nonzero(pureless & wrong), pureless.sum()),
        exactly_one_good=rate(np.count_nonzero(pureless & wrong & d['one_good']),
                             np.count_nonzero(pureless & d['one_good'])))
    for name, mask in [('discovery', h['file_index'] < 3), ('held_out', h['file_index'] >= 3)]:
        result[name] = measurement(tables, d, mask)
    return result


def outcomes(output, dsid, variant, tables, d):
    """Compact truth keyed decisions for independent lineage and report joins."""
    h, t = tables['truth'], tables['tracks']
    present = d['chosen'] >= 0
    payload = {k: h[k] for k in ['sample_entry', 'file_index', 'source_entry', 'run', 'event', 'barcode']}
    payload.update(status=d['status'], selected=d['selected'], wrong=d['wrong'],
        fiducial=d['fiducial'], truth_p_GeV=d['truth_p'], one_good=d['one_good'],
        duplicates=d['duplicates'], truth_features=d['features'])
    for name, field, fill in [('ordinal', 'ordinal', -1), ('qop_GeVinv', 'qop_GeVinv', np.nan),
                             ('purity', 'purity', np.nan)]:
        values = np.full(len(present), fill, dtype=t[field].dtype)
        values[present] = t[field][d['chosen'][present]]
        payload[name] = values
    np.savez_compressed(output/f'{dsid}_{variant}_outcomes.npz', **payload)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    data, decisions, coverage, groups, migrations = {}, {}, {}, {}, {}
    for dsid in DSIDS:
        data[dsid], decisions[dsid] = {}, {}
        for variant in VARIANTS:
            tables, cov = merge_tables(args.root/f'main_{dsid}_v1', variant)
            # Segments and hit arrays have no role in this event-level estimator.
            tables.pop('segments')
            tables['tracks'].pop('hit_identifiers')
            tables['tracks'].pop('hit_offsets')
            data[dsid][variant] = tables
            key = dsid+'|'+variant
            coverage[key] = cov
            decisions[dsid][variant] = {r: derive(tables, dsid, r) for r in [95, 100]}
            groups[key] = {f'r{r}': summarize(tables, decisions[dsid][variant][r]) for r in [95, 100]}
            outcomes(args.output, dsid, variant, tables, decisions[dsid][variant][95])
            print(key, groups[key]['r95']['conditional_misid'], flush=True)
        a, b = [data[dsid][v]['truth'] for v in VARIANTS]
        for field in ['sample_entry', 'barcode', 'pdg', 'production_momentum_F_GeV',
                      'production_position_F_mm', 'truth_sensor_direction_F', 'truth_station_position_F_mm']:
            np.testing.assert_array_equal(a[field], b[field], err_msg='before/after truth: '+field)
        migrations[dsid] = {}
        for radius in [95, 100]:
            before, after = [decisions[dsid][v][radius] for v in VARIANTS]
            mask = before['primary'] | after['primary']
            matrix = np.bincount((7*before['status'][mask]+after['status'][mask]), minlength=49).reshape(7, 7)
            migrations[dsid][f'r{radius}'] = dict(states=STATUS, counts=matrix.tolist(),
                truth_identity_and_physical_input_exact=True,
                conditional_rate_difference=paired_difference(before, after, mask))
    balances, weighted = {}, {}
    for variant in VARIANTS:
        tables = {k: data[k][variant] for k in DSIDS}
        d = {k: decisions[k][variant][95] for k in DSIDS}
        weights, balances[variant] = truth_balance(tables, d)
        weighted[variant] = {}
        for dsid in DSIDS:
            weighted[variant][dsid] = {}
            for split in ['all', 'discovery', 'held_out']:
                files = tables[dsid]['truth']['file_index']
                mask = np.ones(len(files), dtype=bool) if split == 'all' else \
                    (files < 3 if split == 'discovery' else files >= 3)
                w = weights[dsid]*mask
                selected = d[dsid]['selected'] & mask
                supported = w > 0
                row = dict(raw=rate(np.count_nonzero(d[dsid]['wrong'] & selected), selected.sum()),
                    common_support=rate(np.count_nonzero(d[dsid]['wrong'] & selected & supported),
                                        np.count_nonzero(selected & supported)),
                    reweighted=weighted_rate(d[dsid]['wrong'], selected, w),
                    truth_correct_efficiency=weighted_rate(d[dsid]['status'] == 1, d[dsid]['fiducial'] & mask, w),
                    truth_wrong_efficiency=weighted_rate(d[dsid]['wrong'], d[dsid]['fiducial'] & mask, w))
                weighted[variant][dsid][split] = row
    differences = {}
    for variant in VARIANTS:
        differences[variant] = {}
        for split in ['all', 'held_out']:
            minus, plus = [weighted[variant][k][split] for k in ['100117', '100116']]
            row = {}
            for scope in ['raw', 'common_support', 'reweighted']:
                a, b = minus[scope], plus[scope]
                value = a['value']-b['value']
                variances = [r['standard_error']**2 if scope == 'reweighted' else
                    r['value']*(1-r['value'])/r['total'] for r in [a, b]]
                error = float(np.sqrt(sum(variances)))
                row[scope] = dict(mu_minus_minus_mu_plus=value, standard_error=error,
                    interval95=[value-1.96*error, value+1.96*error],
                    ratio_mu_minus_over_mu_plus=a['value']/b['value'] if b['value'] else None)
            differences[variant][split] = row
    result = dict(groups=groups, coverage=coverage, migrations=migrations,
        spectrum_balance=balances, weighted_measurements=weighted, charge_differences=differences,
        configuration=json.loads((Path(__file__).resolve().parents[1]/'configs/charge_study.json').read_text())
            if (Path(__file__).resolve().parents[1]/'configs/charge_study.json').exists() else None,
        uncertainty='One primary muon per source entry; exact binomial intervals; paired ratio delta method for rotation',
        comparison='Native before/after geometry versions, both isolated smoothing ON; native seed priors differ')
    write_json(args.output/'summary.json', result)


if __name__ == '__main__':
    main()
