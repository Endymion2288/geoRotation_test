#!/usr/bin/env python3
"""Validate frozen response curves against inputs and retain all flow sumw2."""
import argparse
import json
from pathlib import Path

import numpy as np

from slope_analysis import derive, EDGES
from slope_response_systematics import common_coverage, MODELS, PERIOD
from slope_statistics import response_cells
from slope_tables import merge_tables, write_json


def observed(values):
    values = np.asarray(values)
    finite = values[np.isfinite(values)]
    counts = np.histogram(finite, EDGES)[0]
    under, over = int((finite < EDGES[0]).sum()), int((finite > EDGES[-1]).sum())
    assert int(counts.sum())+under+over == len(finite)
    return dict(counts=counts.tolist(), sumw2=counts.tolist(), edges=EDGES.tolist(),
                raw_entries=len(values), finite=len(finite), invalid=len(values)-len(finite),
                underflow=under, overflow=over, exact_right_edge=int((finite == EDGES[-1]).sum()),
                right_edge_inclusive=True)


def prediction_flows(ts, reco, cells, index, train_truth, val_truth, train_reco, reference):
    rng = np.random.default_rng(reference['seed'])
    counts, sumw2 = np.zeros(len(EDGES)-1), np.zeros(len(EDGES)-1)
    flows = np.zeros(4)
    expected = 0.
    for cell in np.unique(cells[val_truth]):
        t = np.flatnonzero(train_truth & (cells == cell))
        v = np.flatnonzero(val_truth & (cells == cell))
        r = np.flatnonzero(train_reco & (cells[index] == cell))
        if len(t) < 30 or len(r) < 20:
            continue
        residual = reco[r]-ts[index[r]]
        samples = ts[v, None]+rng.choice(residual, size=(len(v), 8))
        weight = len(r)/len(t)/8
        h = np.histogram(samples, EDGES)[0]
        counts += h*weight
        sumw2 += h*weight**2
        under, over = (samples < EDGES[0]).sum(), (samples >= EDGES[-1]).sum()
        flows += [under*weight, over*weight, under*weight**2, over*weight**2]
        expected += len(v)*len(r)/len(t)
    np.testing.assert_allclose(counts, reference['predicted_counts'], rtol=1e-12, atol=1e-10)
    np.testing.assert_allclose(sumw2, reference['mc_sumw2'], rtol=1e-12, atol=1e-10)
    np.testing.assert_allclose(flows[:2], [reference['prediction_underflow'], reference['prediction_overflow']], rtol=1e-12, atol=1e-10)
    np.testing.assert_allclose(expected, reference['prediction_expected_total'], rtol=1e-12)
    np.testing.assert_allclose(counts.sum()+flows[:2].sum(), expected, rtol=1e-12)
    return dict(counts_and_sumw2_match_frozen=True, expected_total=expected,
                underflow=float(flows[0]), overflow=float(flows[1]),
                underflow_sumw2=float(flows[2]), overflow_sumw2=float(flows[3]),
                in_range_sumw2=float(sumw2.sum()), total_sumw2=float(sumw2.sum()+flows[2:].sum()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--variant', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source = args.root/'main_100044_v1'
    tables, coverage = merge_tables(source, args.variant)
    d = derive(tables, args.variant, '100044')
    t, r = tables['truth'], tables['tracks']
    split = json.loads((source/'manifest.json').read_text())['held_out_split_sample_entry']
    phase = json.loads((args.root/'phase_response_v2'/args.variant/'summary.json').read_text())['responses']
    systems = json.loads((args.root/'response_systematics_v1'/args.variant/'summary.json').read_text())['planes']
    p = np.linalg.norm(t['production_momentum_F_GeV'], axis=1)
    x, y = t['truth_sensor_position_F_mm'][:, :2].T
    selected = d['primary'] & d['selections']['reco_p_window'] & d['pure']
    result = dict(variant=args.variant, source_events=coverage['source_events'], planes={},
                  validation='Input-level reproduction of all frozen response bins, Monte Carlo sumw2 and flows; no reconstruction rerun')
    for i, plane in enumerate(['st1', 'VetoNu']):
        ts, reco = d['truth_slopes'][:, i], d['reconstructed'][:, i]
        valid = d['signal_truth'] & np.isfinite(ts) & np.isfinite(p) & np.isfinite(x) & np.isfinite(y)
        tt = valid & (t['sample_entry'] < split)
        vt = valid & (t['sample_entry'] >= split)
        tr = selected & np.isfinite(reco) & (r['sample_entry'] < split) & valid[d['index']]
        vr = selected & np.isfinite(reco) & (r['sample_entry'] >= split) & valid[d['index']]
        common, support = common_coverage(ts, p, x, y, t['sample_entry'], d['signal_truth'],
            reco, d['index'], selected, r['sample_entry'], split)
        assert support == systems[plane]['support']
        references = {'nominal_phase_v2': (phase[plane+'|reco_p_window'], {}, vt)}
        references.update({name: (reference, MODELS[name], vt & common)
                           for name, reference in systems[plane]['responses'].items()})
        rows = {}
        for name, (reference, config, target) in references.items():
            cells = response_cells(ts, p, x, y, PERIOD, config)
            covered_cells = [v['cell'] for v in reference['cells'] if v['covered']]
            target_reco = vr & target[d['index']]
            covered = target_reco & np.isin(cells[d['index']], covered_cells)
            obs = observed(reco[target_reco])
            good = observed(reco[covered])
            missing = observed(reco[target_reco & ~covered])
            assert obs['counts'] == reference['observed_counts']
            assert good['counts'] == reference['covered_observed_counts']
            assert missing['counts'] == reference['uncovered_validation_reco_counts']
            assert obs['finite'] == reference['validation_reco']
            for key in ['underflow', 'overflow', 'finite']:
                assert good[key]+missing[key] == obs[key]
            rows[name] = dict(observed=obs, covered=good, uncovered=missing,
                prediction=prediction_flows(ts, reco, cells, d['index'], tt, target, tr, reference))
        result['planes'][plane] = rows
        print(args.variant, plane, 'all frozen curve counts and sumw2 reproduced', flush=True)
    write_json(args.output/'validation.json', result)
    print(args.output/'validation.json')


if __name__ == '__main__':
    main()
