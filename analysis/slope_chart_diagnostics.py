#!/usr/bin/env python3
"""Check direction-chart effects against fixed-input truth and afterGeo."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_centroids import collect
from slope_controls import save
from slope_periodic_refit import CONTROLS
from slope_statistics import event_ratio_interval
from slope_tables import histogram, join_truth, slope, write_json


def moments(values):
    finite = values[np.isfinite(values)]
    return dict(finite=len(finite), mean=float(np.mean(finite)),
                rms=float(np.sqrt(np.mean(finite**2))),
                quantile_probabilities=[.025, .16, .5, .84, .975],
                quantiles=np.quantile(finite, [.025, .16, .5, .84, .975]).tolist(),
                absolute_error_68=float(np.quantile(np.abs(finite), .68)))


def analyse(native_paths, control_paths, output, report, label, generated=False,
            validation_entry=None, after_paths=None, control_name='regular_chart'):
    native, count = collect(native_paths)
    control, other_count = collect(control_paths)
    assert count == other_count
    r, q, truth = native['tracks'], control['tracks'], native['truth']
    for key in ['source_entry', 'ordinal']:
        assert np.array_equal(r[key], q[key]), key
    same_hits = r['hit_signatures'] == q['hit_signatures']
    same_truth = r['barcode'] == q['barcode']
    if control_name == 'regular_chart':
        assert same_hits.all() and same_truth.all()
    index, matched = join_truth(native)
    selected = matched & (np.abs(truth['pdg'][index]) == 13) & (r['purity'] >= .7) & (r['station_bits'] & 14 == 14)
    if validation_entry is not None:
        selected &= r['source_entry'] >= validation_entry
    momentum = np.linalg.norm(r['upstream_F'][:, 3:], axis=1)
    truth_p = np.linalg.norm(truth['production_momentum_F_GeV'][index], axis=1)
    result = dict(source_events=count, selected_candidates=int(selected.sum()),
                  control=control_name, control_definition=CONTROLS[control_name][1],
                  native_paths=[str(path) for path in native_paths],
                  control_paths=[str(path) for path in control_paths],
                  after_paths=[str(path) for path in after_paths] if after_paths else [],
                  validation_source_entry_min=validation_entry,
                  generated_st1_truth=generated,
                  candidate_identity_unchanged=True,
                  same_complete_hit_sets=int(same_hits.sum()),
                  majority_changes=int((~same_truth).sum()), planes={},
                  momentum_relative_change=moments((np.linalg.norm(q['upstream_F'][:, 3:], axis=1)-momentum)/momentum))
    after = collect(after_paths)[0] if after_paths else None
    for plane, i in [('st1', 0), ('sensor', 2)]:
        a, b = [np.where(t['projection_valid'][:, i] == 1,
                         slope(t['projection'][:, i, 3:6]), np.nan) for t in [r, q]]
        reference = (slope(truth['production_momentum_F_GeV'][index]) if generated and i == 0
                     else slope(truth['truth_model_momentum_F_GeV'][index, i]) if i == 0
                     else slope(truth['truth_sensor_direction_F'][index]))
        if generated and i == 0:
            assert np.max(np.abs(truth['production_position_F_mm'][:, 2])) < 1.e-9
        finite = selected & np.isfinite(a) & np.isfinite(b) & np.isfinite(reference)
        modes = dict(native=a, truth=reference)
        modes[control_name] = b
        item = dict(reference=('Exact generated st1 direction' if generated and i == 0
                               else 'SiHit-chord/production-p same-plane model' if i == 0
                               else 'Same-sensor SiHit chord direction'), selections={})
        windows = [('all_p', np.ones(len(a), bool)),
                   ('native_p', (momentum > 50) & (momentum < 250)),
                   ('truth_p', (truth_p > 50) & (truth_p < 250))]
        windows += [(name+'|same_hits_truth', window & same_hits & same_truth)
                    for name, window in windows.copy()]
        for name, window in windows:
            mask = finite & window
            c = dict(common_finite=int(mask.sum()), windows={}, residuals={})
            for width in [25, 50, 100, 200]:
                inside = {mode: np.abs(values) < width*1.e-6 for mode, values in modes.items()}
                c['windows'][str(width)] = dict(
                    half_width_microradians=width,
                    counts={mode: int((mask & valid).sum()) for mode, valid in inside.items()},
                    native_minus_truth=event_ratio_interval(r['sample_entry'],
                        mask*(inside['native'].astype(int)-inside['truth'].astype(int)), mask, (0, count)),
                    control_minus_truth=event_ratio_interval(r['sample_entry'],
                        mask*(inside[control_name].astype(int)-inside['truth'].astype(int)), mask, (0, count)),
                    control_minus_native=event_ratio_interval(r['sample_entry'],
                        mask*(inside[control_name].astype(int)-inside['native'].astype(int)), mask, (0, count)))
            for mode, values in [('native', a), (control_name, b)]:
                residual = values-reference
                variance = (r if mode == 'native' else q)['projection'][:, i, 6]
                pull_valid = mask & np.isfinite(variance) & (variance > 0)
                c['residuals'][mode] = dict(
                    summary=moments(residual[mask]),
                    histogram=histogram(residual[mask], np.linspace(-.003, .003, 301)),
                    pull_summary=moments(residual[pull_valid]/np.sqrt(variance[pull_valid])),
                    pull_histogram=histogram(residual[pull_valid]/np.sqrt(variance[pull_valid]), np.linspace(-10, 10, 201)),
                    invalid_pull_variance=int((mask & ~pull_valid).sum()),
                    large_error={str(width): event_ratio_interval(r['sample_entry'],
                        mask & (np.abs(residual) > width*.001), mask, (0, count))
                        for width in [1, 3]})
            item['selections'][name] = c
        result['planes'][plane] = item
        fixed = finite & (momentum > 50) & (momentum < 250)
        item['figure_histograms'] = {}
        for values, name in [(a, 'Native'), (b, CONTROLS[control_name][0]), (reference, 'Truth')]:
            h = histogram(values[fixed], np.linspace(-.0015, .0015, 151))
            item['figure_histograms'][name] = h
        if after is not None:
            ar = after['tracks']
            ti, linked = join_truth(after)
            chosen = linked & (np.abs(after['truth']['pdg'][ti]) == 13) & (ar['purity'] >= .7) & (ar['station_bits'] & 14 == 14)
            ap = np.linalg.norm(ar['upstream_F'][:, 3:], axis=1)
            chosen &= (ap > 50) & (ap < 250) & (ar['projection_valid'][:, i] == 1)
            sy = slope(ar['projection'][:, i, 3:6])
            h = histogram(sy[chosen], np.linspace(-.0015, .0015, 151))
            item['figure_histograms']['afterGeo_own_cohort'] = h
            item['afterGeo'] = dict(selected=int(chosen.sum()), near_zero=int((chosen & (np.abs(sy) < .0001)).sum()),
                                   definition='Same generated source; own afterGeo candidates, hits and reconstructed-p cut; descriptive version reference')
        render_plane(item, label, plane, count, control_name, output, report)
    return result


def render_plane(item, label, plane, count, control, output, report):
    """Render saved counts so display fixes do not repeat the event bootstrap."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    for name, h in item['figure_histograms'].items():
        after = name == 'afterGeo_own_cohort'
        axes[0].stairs(h['counts'], np.array(h['edges'])*1000,
            label='AfterGeo (own cohort)' if after else name, linestyle='--' if after else '-')
    for mode, name in [('native', 'Native'), (control, CONTROLS[control][0])]:
        for ax, key, scale in [(axes[1], 'histogram', 1000), (axes[2], 'pull_histogram', 1)]:
            h = item['selections']['native_p']['residuals'][mode][key]
            ax.stairs(h['counts'], np.array(h['edges'])*scale, label=name)
    for axis, key in [(axes[1], 'histogram'), (axes[2], 'pull_histogram')]:
        flows = []
        for mode, name in [('native', 'Native'), (control, 'Control')]:
            h = item['selections']['native_p']['residuals'][mode][key]
            flows.append(f"{name} outside range: {h['underflow']+h['overflow']}; invalid: {h['invalid']}")
        axis.text(.03, 1.02, '\n'.join(flows), transform=axis.transAxes,
                  va='bottom', fontsize=7)
    for ax, xlabel in zip(axes, ['Physical sy [mrad]', 'Reco sy - reference sy [mrad]', 'Residual / propagated sigma(sy)']):
        ax.set(xlabel=xlabel, ylabel='Candidates / bin')
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    axes[1].set(yscale='log', ylim=(.8, None))
    fig.suptitle(f'{label} | {plane} | {count:,} source events | Actual Calypso\nFixed native candidates and p = 50-250 GeV; purity >= 0.7; complete-hit changes reported')
    save(fig, output, label+'_'+plane+'_'+control+'_diagnostics', report)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--regular-root', type=Path)
    parser.add_argument('--native-root', type=Path)
    parser.add_argument('--gun-root', type=Path)
    parser.add_argument('--display-summary', type=Path,
                        help='Render frozen counts without reading events or recalculating metrics')
    parser.add_argument('--control', choices=[k for k in CONTROLS if k.startswith('regular_')],
                        default='regular_chart')
    parser.add_argument('--mode', choices=['all', 'rdo', 'gun'], default='all')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.display_summary:
        data = json.loads(args.display_summary.read_text())
        with PdfPages(args.output/'chart_diagnostics_report.pdf') as report:
            for label, case in data['cases'].items():
                for plane, item in case['planes'].items():
                    render_plane(item, label, plane, case['source_events'], case['control'], args.output, report)
        data['display_source_summary'] = str(args.display_summary)
        write_json(args.output/'summary.json', data)
        return
    cases = {}
    with PdfPages(args.output/'chart_diagnostics_report.pdf') as report:
        if args.mode in ['all', 'rdo']:
            cases['production_muplus'] = analyse(
                sorted((args.native_root/'production/observed/shards/beforeGeo').glob('s*')),
                sorted((args.regular_root/('production/'+args.control+'/shards/beforeGeo')).glob('s*')),
                args.output, report, 'production_muplus', validation_entry=250000,
                control_name=args.control)
        if args.mode in ['all', 'gun']:
            for charge in ['muplus', 'muminus']:
                p = args.gun_root/'shards'/charge
                cases[charge] = analyse([p/'beforeGeo_native'], [p/('beforeGeo_'+args.control)],
                    args.output, report, charge, generated=True, after_paths=[p/'afterGeo'],
                    control_name=args.control)
    write_json(args.output/'summary.json', dict(cases=cases,
        definition='Actual Calypso outputs; no transport or tracking recomputation. Windows fixed at 25/50/100/200 microradians; 100 is the primary window. Intervals use source-event empirical bootstrap; narrower/wider windows are diagnostics.'))
    print(args.output/'summary.json')


if __name__ == '__main__':
    main()
