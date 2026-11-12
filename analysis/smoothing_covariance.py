#!/usr/bin/env python3
"""Event-paired seed-variance tail statistics and reference-macro figures."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

CASES = ['beforeGeo', 'afterGeo', 'afterGeo_originalVariance']
STAGES = ['prt', 'flt', 'smt']
LAYERS = [(s, l) for s in [1, 2, 3] for l in [0, 1, 2]]


def paired_rate(first, second, field):
    """Ratio difference with source-event clustering and paired event covariance."""
    n_first, n_second = first[:, 0].sum(), second[:, 0].sum()
    p_first = first[:, field].sum()/n_first
    p_second = second[:, field].sum()/n_second
    influence = ((second[:, field]-p_second*second[:, 0])/n_second
                 -(first[:, field]-p_first*first[:, 0])/n_first)
    error = np.sqrt(len(first)/(len(first)-1)*np.square(influence).sum())
    delta = p_second-p_first
    return dict(first_rate=float(p_first), second_rate=float(p_second),
        delta_pp=float(100*delta), lower95_pp=float(100*(delta-1.96*error)),
        upper95_pp=float(100*(delta+1.96*error)), method='paired source-event cluster delta method')


def load_events(root, case):
    tasks = json.loads((root/(case+'_tasks.json')).read_text())
    settings = json.loads((root/'inputs.json').read_text())['cases'][case]
    counts, tracks, keys = [], [], []
    for task in tasks:
        run = json.loads((Path(task['output'])/'run.json').read_text())
        assert str(root/settings['software']/'lib/libFaserActsKalmanFilter.so') in run['loaded_acts_libraries']
        for prop, factor, sigma in zip(['covLoc0', 'covLoc1', 'covPhi', 'covTheta', 'covQOverP'],
                settings['initial_variance_inflation'], [.9, .019, .033, .00029, .0001]):
            assert abs(run['seed_covariance'][prop]/(factor*sigma*sigma)-1) < 1e-12
        with np.load(Path(task['output'])/'event_counts.npz') as data:
            entries = data['source_entry']
            assert entries.tolist() == list(range(task['skip'], task['skip']+task['events']))
            keys.extend((task['input'], int(entry)) for entry in entries)
            counts.append(data['counts'])
            tracks.append(data['tracks'])
    assert len(keys) == len(set(keys)) == 20000
    return np.concatenate(counts), np.concatenate(tracks), keys


def statistics(root):
    data = {case: load_events(root, case) for case in CASES}
    assert all(data[case][2] == data[CASES[0]][2] for case in CASES)
    summary = dict(source_events=20000, files=20, isolate_smoothing=True,
        reference='SiHit-linked HepMC production momentum, same majority barcode',
        denominator='finite reconstructed and truth q/p, per layer and stage',
        thresholds={'tail_half': 'residual < -0.5', 'tail_point8': 'residual < -0.8'}, cases={}, comparisons={})
    for case, (counts, tracks, _) in data.items():
        total = counts.sum(axis=0)
        frozen = json.loads((root/'reductions'/(case+'_truth_summary.json')).read_text())
        assert int(tracks[:, 0].sum()) == frozen['counts']['writer_tracks']
        assert int(tracks[:, 1].sum()) == frozen['counts']['github_selected_tracks']
        layers = {}
        for stage_index, stage in enumerate(STAGES):
            for layer_index, (station, layer) in enumerate(LAYERS):
                values = total[stage_index, layer_index]
                source = frozen['per_layer'][f'{station}_{layer}']
                for field_index, field in enumerate(['finite', 'below_minus_half', 'below_minus_point8']):
                    assert int(values[field_index]) == source[stage+'_native_'+field]
                layers[f'{stage}_S{station}L{layer}'] = dict(finite=int(values[0]),
                    raw_states=int(source['states']),
                    invalid_reco=int(source.get(stage+'_native_invalid_reco', 0)),
                    invalid_truth=int(source.get(stage+'_native_invalid_truth', 0)),
                    tail_half=int(values[1]), tail_point8=int(values[2]),
                    tail_half_pct=float(100*values[1]/values[0]),
                    tail_point8_pct=float(100*values[2]/values[0]))
        summary['cases'][case] = dict(writer_tracks=int(tracks[:, 0].sum()),
            selected_tracks=int(tracks[:, 1].sum()), layers=layers)
    for first, second in [(CASES[1], CASES[2]), (CASES[0], CASES[1]), (CASES[0], CASES[2])]:
        comparison = {}
        for stage_index, stage in enumerate(STAGES):
            for layer_index, (station, layer) in enumerate(LAYERS):
                comparison[f'{stage}_S{station}L{layer}'] = {
                    label: paired_rate(data[first][0][:, stage_index, layer_index],
                        data[second][0][:, stage_index, layer_index], field)
                    for label, field in [('tail_half', 1), ('tail_point8', 2)]}
        summary['comparisons'][first+'_vs_'+second] = comparison
    (root/'reductions/tail_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    with (root/'reductions/tail_comparisons.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=['comparison', 'layer_stage', 'threshold',
            'first_rate', 'second_rate', 'delta_pp', 'lower95_pp', 'upper95_pp', 'method'])
        writer.writeheader()
        for pair, rows in summary['comparisons'].items():
            for layer, thresholds in rows.items():
                for threshold, value in thresholds.items():
                    writer.writerow(dict(comparison=pair, layer_stage=layer, threshold=threshold, **value))
    return summary


def tail_figures(root, summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from reference_figures import read_histograms
    histograms = {case: read_histograms(root/'reductions'/(case+'_native_3station.root'),
        'tail_'+case) for case in CASES}
    labels = ['beforeGeo: original variance', 'afterGeo: enlarged variance', 'afterGeo: original variance']
    colors = ['#222222', '#d62728', '#1f77b4']
    fig, axes = plt.subplots(3, 3, figsize=(15, 10), sharex=True)
    bins = {}
    for axis, (station, layer) in zip(axes.flat, LAYERS):
        for case, label, color in zip(CASES, labels, colors):
            hist = histograms[case][(0, station, layer)]
            edges = np.array([hist.GetBinLowEdge(i) for i in range(1, 102)])
            counts = np.array([hist.GetBinContent(i) for i in range(1, 101)])
            axis.stairs(counts/hist.GetEntries(), edges, label=label, color=color, linewidth=1.5)
            bins[f'{case}_S{station}L{layer}'] = dict(edges=edges.tolist(), counts=counts.tolist(),
                sumw2=[hist.GetBinError(i)**2 for i in range(1, 101)],
                underflow=hist.GetBinContent(0), overflow=hist.GetBinContent(101),
                raw_entries=hist.GetEntries())
        axis.set_title(f'Station {station}, layer {layer}')
        axis.set_yscale('log')
        axis.set_xlim(-1.2, .3)
        axis.set_ylim(1e-5, 1)
        axis.axvline(-.5, color='.6', linewidth=.7, linestyle=':')
        axis.axvline(-.8, color='.6', linewidth=.7, linestyle=':')
        axis.set_xlabel('(reco p - truth p) / truth p')
        axis.set_ylabel('Entries / raw Entries')
        entries = [int(histograms[case][(0, station, layer)].GetEntries()) for case in CASES]
        axis.text(.03, .96, 'Raw N: '+' / '.join(f'{n:,}' for n in entries),
            transform=axis.transAxes, va='top', fontsize=8)
    handles, legend_labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc='upper center', bbox_to_anchor=(.5, .90),
        ncol=3, frameon=False, fontsize=10)
    fig.suptitle(f"110051 | {summary['files']} files | {summary['source_events']:,} source events | candidate-isolated smoothing\n"
        'Selected CKF smoothed sensor states; native FASER / rotated ACTS frames; SiHit-linked HepMC production truth\n'
        'Cuts: chi2/NDF < 20, NDF > 9, rmax < 90 mm; normalization includes NaN / flow entries')
    fig.tight_layout(rect=(0, 0, 1, .91))
    prefix = root/'figures/110051_smoothed_left_tail'
    fig.savefig(str(prefix)+'.png', dpi=150)
    fig.savefig(str(prefix)+'.pdf')
    plt.close(fig)
    (root/'reductions/smoothed_bin_counts.json').write_text(json.dumps(bins, indent=2)+'\n')
    import sys
    sys.path.insert(0, str(root/'software_tools'))
    from pypdf import PdfReader, PdfWriter
    writer = PdfWriter()
    for path in [root/'figures/reference_comparison.pdf', prefix.with_suffix('.pdf')]:
        for page in PdfReader(path).pages:
            writer.add_page(page)
    with (root/'figures/smoothing_covariance_comparison.pdf').open('wb') as stream:
        writer.write(stream)
    index_path = root/'figures/figure_index.json'
    index = json.loads(index_path.read_text())
    index.append(dict(name=prefix.name, description=f"110051 | {summary['files']} files | {summary['source_events']:,} events | smoothed left tail; raw-Entries normalization; isolated smoothing",
        reference='geometryRotation.pdf p17; supplemental logarithmic tail view'))
    index_path.write_text(json.dumps(index, indent=2)+'\n')
    with (root/'figures/figure_index.md').open('a') as stream:
        stream.write(f"\n- [{prefix.name}]({prefix.name}.png): isolated smoothing, {summary['source_events']:,} source events, smoothed left tail\n")
    combined = PdfReader(root/'figures/smoothing_covariance_comparison.pdf')
    assert len(combined.pages) == len(index)
    for page, item in zip(combined.pages, index):
        original = PdfReader(root/'figures'/(item['name']+'.pdf')).pages[0]
        assert page.get_contents().get_data() == original.get_contents().get_data()
        assert list(page.mediabox) == list(original.mediabox)
        assert list(page.cropbox) == list(original.cropbox)
    (root/'figures/combined_figure_validation.json').write_text(json.dumps(dict(
        page_count=len(index), contents_and_page_boxes_preserved=True,
        file='smoothing_covariance_comparison.pdf'), indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    from reference_figures import merge, plot
    merge(args.root, CASES)
    summary = statistics(args.root)
    plot(args.root, CASES)
    tail_figures(args.root, summary)
    print(json.dumps(dict(complete=True, cases={k: v['selected_tracks']
        for k, v in summary['cases'].items()})), flush=True)


if __name__ == '__main__':
    main()
