#!/usr/bin/env python3
"""Test empirical-response grouping on a common held-out truth population."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_analysis import derive, EDGES
from slope_controls import save
from slope_statistics import conditional_response, response_cells
from slope_tables import merge_tables, write_json


MODELS = {
    'phase8': {},
    'phase4': {'phase_bins': 4},
    'phase16': {'phase_bins': 16},
    'phase8_fine_p': {'momentum_edges_GeV': [50., 75., 100., 150., 250., 500., 1000.]},
    'phase8_fine_position': {'transverse_position_edges_mm': [-30., 0., 30.]},
}
PERIOD = .08/63.


def common_coverage(truth_slope, truth_p, x, y, entries, selected_truth,
                    reco, index, selected_reco, reco_entries, split):
    valid_truth = selected_truth & np.isfinite(truth_slope) & np.isfinite(truth_p) & np.isfinite(x) & np.isfinite(y)
    train_truth = valid_truth & (entries < split)
    train_reco = selected_reco & np.isfinite(reco) & (reco_entries < split) & valid_truth[index]
    common = valid_truth.copy()
    support = {}
    for name, config in MODELS.items():
        cells = response_cells(truth_slope, truth_p, x, y, PERIOD, config)
        size = int(cells.max(initial=0))+1
        nt = np.bincount(cells[train_truth], minlength=size)
        nr = np.bincount(cells[index[train_reco]], minlength=size)
        covered = (nt >= 30) & ((nr == 0) | (nr >= 20))
        mask = covered[cells]
        common &= mask
        validation = selected_reco & np.isfinite(reco) & (reco_entries >= split) & valid_truth[index]
        support[name] = dict(covered_validation_reco=int((validation & mask[index]).sum()),
                             uncovered_validation_reco=int((validation & ~mask[index]).sum()))
    return common, support


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--variant', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    tables, coverage = merge_tables(args.root, args.variant)
    d = derive(tables, args.variant, '100044')
    truth, tracks = tables['truth'], tables['tracks']
    split = json.loads((args.root/'manifest.json').read_text())['held_out_split_sample_entry']
    windows = json.loads((args.reference/'peak_windows.json').read_text())['windows']
    p = np.linalg.norm(truth['production_momentum_F_GeV'], axis=1)
    x, y = truth['truth_sensor_position_F_mm'][:, :2].T
    selected = d['primary'] & d['selections']['reco_p_window'] & d['pure']
    result = dict(variant=args.variant, source_events=coverage['source_events'], models=MODELS,
                  discovery_validation_split=split, phase_period=PERIOD,
                  definition='Five predeclared grouping choices; identical common covered validation truths; training always uses the full discovery population.',
                  planes={})
    with PdfPages(args.output/'response_systematics.pdf') as report:
        for i, plane in enumerate(['st1', 'VetoNu']):
            key = plane+'|reco_p_window'
            reco, ts = d['reconstructed'][:, i], d['truth_slopes'][:, i]
            common, support = common_coverage(ts, p, x, y, truth['sample_entry'], d['signal_truth'],
                reco, d['index'], selected, tracks['sample_entry'], split)
            responses = {}
            for name, config in MODELS.items():
                response = conditional_response(ts, p, x, y, truth['sample_entry'], reco, d['index'], selected,
                    tracks['sample_entry'], split, EDGES, windows=windows[key], truth_selected=d['signal_truth'],
                    phase_period=PERIOD, conditioning=config, validation_truth_mask=common)
                assert response['uncovered_validation_reco'] == 0
                if responses:
                    prior = next(iter(responses.values()))
                    assert response['covered_observed_counts'] == prior['covered_observed_counts']
                    assert response['validation_reco'] == prior['validation_reco']
                responses[name] = response
                print(args.variant, plane, name, response['validation_reco'], flush=True)
            valid_truth = d['signal_truth'] & np.isfinite(ts) & np.isfinite(p) & np.isfinite(x) & np.isfinite(y)
            validation = selected & np.isfinite(reco) & (tracks['sample_entry'] >= split) & valid_truth[d['index']]
            missing = validation & ~common[d['index']]
            nominal = responses['phase8']
            assert nominal['validation_reco']+int(missing.sum()) == int(validation.sum())
            windows_result = []
            for j, w in enumerate(windows[key]):
                rows = [r['covered_peak_comparison'][j] for r in responses.values()]
                total_window = int((validation & (np.abs(reco-w['center']) < w['half_width'])).sum())
                uncovered_window = int((missing & (np.abs(reco-w['center']) < w['half_width'])).sum())
                upper = max(0., max(r['difference_entries_interval_95'][1] for r in rows))
                windows_result.append(dict(window=w, full_validation_window_entries=total_window,
                    common_covered_window_entries=rows[0]['observed_entries'], uncovered_window_entries=uncovered_window,
                    predicted_entries_range=[min(r['predicted_entries'] for r in rows), max(r['predicted_entries'] for r in rows)],
                    observed_minus_predicted_entries_range=[min(r['observed_minus_predicted_entries'] for r in rows), max(r['observed_minus_predicted_entries'] for r in rows)],
                    per_model_interval_envelope_95=[min(r['difference_entries_interval_95'][0] for r in rows), max(r['difference_entries_interval_95'][1] for r in rows)],
                    covered_positive_excess_upper_entries=upper,
                    full_positive_excess_upper_entries_including_uncovered=upper+uncovered_window,
                    full_positive_excess_upper_fraction_including_uncovered=(upper+uncovered_window)/validation.sum(),
                    interpretation='Per frozen window, envelope across five measured grouping choices; all uncovered window entries assigned to unexplained excess. Not a simultaneous or global peak limit.'))
            result['planes'][plane] = dict(support=support, full_validation_reco=int(validation.sum()),
                common_covered_reco=nominal['validation_reco'], uncovered_reco=int(missing.sum()),
                windows=windows_result, responses=responses)
            fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
            axes[0].stairs(nominal['covered_observed_counts'], EDGES, label='Held-out reconstruction: common covered cells', color='black')
            for name, response in responses.items():
                axes[0].stairs(response['predicted_counts'], EDGES, label=name, alpha=.8)
                axes[1].stairs(np.array(nominal['covered_observed_counts'])-response['predicted_counts'], EDGES, label=name)
            axes[0].stairs(np.histogram(reco[missing], EDGES)[0], EDGES, label='Uncovered held-out reconstruction', color='gray', linestyle=':')
            for ax in axes:
                ax.set(xlim=(-.006, .004), ylabel='Entries / 0.2 mrad')
                ax.grid(alpha=.2)
                ax.legend(fontsize=7)
            axes[1].axhline(0, color='k', lw=.6)
            axes[1].set_xlabel('Physical sy')
            fig.suptitle(f'{args.variant} | {plane} | Response grouping systematics | {coverage["source_events"]:,} source events\nNative p = 50-250 GeV; discovery < {split:,}; same common covered held-out truths')
            save(fig, args.output, plane+'_response_systematics', report)
    write_json(args.output/'summary.json', result)
    print(args.output/'summary.json')


if __name__ == '__main__':
    main()
