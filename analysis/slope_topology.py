#!/usr/bin/env python3
"""Condition actual Calypso A-prime tracks on their generator-pair topology."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_analysis import derive
from slope_controls import save
from slope_tables import histogram, merge_tables, write_json


GROUPS = {
    'decay_z_mm': ([-np.inf, -1000, -500, 0, np.inf], 'Decay vertex z_F [mm]'),
    'opening_mrad': ([0, .05, .1, .25, np.inf], 'Production opening angle [mrad]'),
    'st1_centroid_separation_mm': ([0, .1, .5, 2, np.inf], 'ST1 truth SiHit centroid separation [mm]'),
}


def pair_topology(truth):
    """Use the two status-1 e+/e- at a common production vertex per source event."""
    n = len(truth['barcode'])
    fields = {name: np.full(n, np.nan) for name in GROUPS}
    selected = (np.abs(truth['pdg']) == 11) & (truth['status'] == 1)
    indices = np.flatnonzero(selected)
    order = indices[np.argsort(truth['sample_entry'][indices], kind='stable')]
    entries, starts, counts = np.unique(truth['sample_entry'][order], return_index=True, return_counts=True)
    counters = dict(events_with_electrons=len(entries), valid_pairs=0, wrong_multiplicity=0,
                    wrong_charge_pair=0, different_vertices=0, missing_st1_centroids=0)
    for start, count in zip(starts, counts):
        if count != 2:
            counters['wrong_multiplicity'] += 1
            continue
        pair = order[start:start+2]
        if truth['pdg'][pair].sum() != 0:
            counters['wrong_charge_pair'] += 1
            continue
        vertices = truth['production_position_F_mm'][pair]
        if not np.all(np.isfinite(vertices)) or not np.allclose(vertices[0], vertices[1], rtol=0, atol=1e-6):
            counters['different_vertices'] += 1
            continue
        momenta = truth['production_momentum_F_GeV'][pair]
        directions = momenta / np.linalg.norm(momenta, axis=1)[:, None]
        opening = np.arctan2(np.linalg.norm(np.cross(*directions)), np.dot(*directions)) * 1000
        fields['decay_z_mm'][pair] = vertices[:, 2]
        fields['opening_mrad'][pair] = opening
        positions = truth['truth_station_position_F_mm'][pair, 1, :2]
        valid = (truth['truth_hit_counts'][pair, 1] > 0).all() and np.isfinite(positions).all()
        if valid:
            fields['st1_centroid_separation_mm'][pair] = np.linalg.norm(positions[0] - positions[1])
        else:
            counters['missing_st1_centroids'] += 1
        counters['valid_pairs'] += 1
    return fields, counters


def fraction(events, numerator, denominator, count):
    """Ratio with source-event cluster delta uncertainty, including paired tracks."""
    den = np.bincount(events[denominator], minlength=count)
    num = np.bincount(events[numerator & denominator], minlength=count)
    total = den.sum()
    if not total:
        return dict(numerator=0, denominator=0, estimate=None, interval_95=None)
    estimate = num.sum()/total
    error = np.sqrt(count/(count-1) * np.sum((num-estimate*den)**2))/total if count > 1 else 0
    return dict(numerator=int(num.sum()), denominator=int(total), estimate=float(estimate),
                interval_95=[float(max(0, estimate-1.95996398454*error)), float(min(1, estimate+1.95996398454*error))],
                method='Source-event cluster delta method; descriptive interval, not a zero-count upper limit')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--dsid', choices=['110033', '110051'], required=True)
    parser.add_argument('--variant', choices=['beforeGeo', 'afterGeo'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    tables, coverage = merge_tables(args.root, args.variant)
    d = derive(tables, args.variant, args.dsid)
    truth, tracks = tables['truth'], tables['tracks']
    topology, pairs = pair_topology(truth)
    n = coverage['source_events']
    base = d['long'] & d['primary'] & d['pure']
    result = dict(dsid=args.dsid, variant=args.variant, coverage=coverage, pair_counts=pairs,
        selection='Matched status-1 primary e+/-; long ST1-3 and SDO purity >= 0.7; duplicates retained',
        definition='Actual Calypso saved reconstruction and truth; topology grouping only, no reconstruction or transport recalculated',
        truth_reference='Fixed-plane propagated production-p/SiHit-chord model; not local G4 momentum',
        separation_reference='Transverse separation of actual station-1 truth SiHit centroids, not a common-z extrapolation',
        groups={})
    groups = dict(GROUPS)
    if args.dsid == '110051':
        # This production has a much narrower opening spectrum. The first
        # pass left only six reco-p selected tracks above 0.05 mrad.
        # A second count check found four entries below 5 microradians;
        # merge that interval with 5-10 rather than interpret its shape.
        groups['opening_mrad'] = ([0, .01, .02, np.inf], GROUPS['opening_mrad'][1])
    with PdfPages(args.output/'topology_report.pdf') as report:
        for name, (edges, label) in groups.items():
            values = topology[name][d['index']]
            finite = np.isfinite(values)
            group = dict(edges=[float(x) if np.isfinite(x) else ('-inf' if x < 0 else 'inf') for x in edges],
                label=label, missing_topology=int((base & ~finite).sum()), bins=[])
            fig, axes = plt.subplots(2, 3, figsize=(14, 8))
            for lo, hi in zip(edges[:-1], edges[1:]):
                selected = base & finite & (values >= lo) & (values < hi)
                item = dict(lower=str(lo), upper=str(hi), selected=int(selected.sum()), selections={})
                for cut, mask in [('all_p', selected), ('reco_p_50_100', selected & d['selections']['reco_p_window']),
                                  ('truth_p_50_100', selected & d['selections']['truth_p_window'])]:
                    planes = {}
                    for pi, plane in enumerate(['st1', 'VetoNu', 'sensor']):
                        reco, reference = d['reconstructed'][:, pi], d['truth_slopes'][d['index'], pi]
                        common = mask & np.isfinite(reco) & np.isfinite(reference)
                        residual = reco-reference
                        planes[plane] = dict(selected=int(mask.sum()), common_finite=int(common.sum()),
                            failed=int((mask & (tracks['projection_valid'][:, pi] == 0)).sum()),
                            missing_reference=int((mask & (tracks['projection_valid'][:, pi] < 0)).sum()),
                            near_zero=fraction(tracks['sample_entry'], np.abs(reco) < .0001, common, n),
                            momentum_tail=fraction(tracks['sample_entry'], d['momentum_tail'], common, n),
                            wrong_charge=fraction(tracks['sample_entry'], d['wrong_charge'], common, n),
                            residual_mean=float(np.mean(residual[common])) if common.any() else None,
                            residual_std=float(np.std(residual[common])) if common.any() else None,
                            reco=histogram(reco[mask], np.linspace(-.01, .01, 201)),
                            truth=histogram(reference[mask], np.linspace(-.01, .01, 201)),
                            residual=histogram(residual[mask], np.linspace(-.005, .005, 201)))
                    item['selections'][cut] = planes
                group['bins'].append(item)
                for pi, plane in enumerate(['st1', 'VetoNu']):
                    for ci, field in enumerate(['truth', 'reco', 'residual']):
                        h = item['selections']['reco_p_50_100'][plane][field]
                        axes[pi, ci].stairs(np.asarray(h['counts'])/max(1, h['finite']), np.asarray(h['edges'])*1000,
                            label=f'{lo:g} to {hi:g}; N={h["finite"]}')
                        axes[pi, ci].set(title=plane+' '+field, xlabel='Physical sy [mrad]' if field != 'residual' else 'Reco - truth model sy [mrad]',
                                         ylabel='Entries / finite entries')
                        axes[pi, ci].grid(alpha=.2)
                        axes[pi, ci].legend(fontsize=7)
            assert sum(row['selected'] for row in group['bins']) + group['missing_topology'] == int(base.sum())
            result['groups'][name] = group
            fig.suptitle(f'{args.dsid} | {args.variant} | {n:,} source events | {label}\n'
                         'Actual native Calypso; matched pure long tracks; reco p=50-100 GeV; conditional shapes')
            save(fig, args.output, name, report)
    write_json(args.output/'summary.json', result)
    print(json.dumps(dict(dsid=args.dsid, variant=args.variant, source_events=n, pairs=pairs)))


if __name__ == '__main__':
    main()
