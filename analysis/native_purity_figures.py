#!/usr/bin/env python3
"""Validate and plot native-only geometry comparison with final-track purity."""
import argparse
from collections import Counter
import csv
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile

import numpy as np

from reference_figures import COMPONENTS, LAYERS, histogram_name, overlay, read_histograms
from smoothing_covariance import paired_rate

CASES = ['beforeGeo', 'afterGeo']
STAGES = ['prt', 'flt', 'smt']


def validate_smoke(root):
    import ROOT
    reference = Path(json.loads((root/'inputs.json').read_text())['reference'])
    checks = {}
    for case in CASES:
        source = root/'smoke'/case
        assert (source/'exit_code.txt').read_text().strip() == '0'
        old = read_histograms(reference/'smoke'/case/'native_3station_counts.root', 'old_'+case)
        baseline = read_histograms(source/'baseline_3station_counts.root', 'base_'+case)
        selected = read_histograms(source/'native_3station_counts.root', 'new_'+case)
        for key, hist in baseline.items():
            assert hist.GetEntries() == old[key].GetEntries()
            for i in range(102):
                assert hist.GetBinContent(i) == old[key].GetBinContent(i), (case, key, i)
                assert selected[key].GetBinContent(i) <= hist.GetBinContent(i)
        summary = json.loads((source/'native_purity.json').read_text())
        assert summary['counts']['source_events'] == 20
        assert summary['counts']['selected_tracks'] < summary['counts']['github_selected_tracks']
        checks[case] = dict(unchanged_baseline_histograms=27, bins_and_flow_checked=2754,
            raw_entries_equal=True, purity_counts=summary['purity_counts'], counts=summary['counts'])
    (root/'smoke_validation.json').write_text(json.dumps(checks, indent=2)+'\n')
    print(json.dumps(checks, indent=2))


def merge(root):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    destination = root/'reductions'
    destination.mkdir(exist_ok=True)
    controls = json.loads((root/'inputs.json').read_text())
    coverage, summaries, event_data, bins, purity_effects = {}, {}, {}, {}, {}
    for case in CASES:
        tasks = json.loads((root/(case+'_tasks.json')).read_text())
        histograms, counts, purity_counts = {}, Counter(), Counter()
        layers = {f'{s}_{l}': Counter() for s, l in LAYERS}
        events, tracks, keys, resources = [], [], [], []
        for task in tasks:
            source = Path(task['output'])
            assert (source/'exit_code.txt').read_text().strip() == '0', source
            assert json.loads((source/'task.json').read_text()) == task
            metadata = json.loads((source/'input_metadata.json').read_text())
            assert metadata == dict(input=task['input'], collection_tree_entries=1000, skip=0, events=1000)
            run = json.loads((source/'run.json').read_text())
            assert run['success'] and run['isolate_smoothing'] and run['args']['native_final_purity']
            assert not run['audit'] and not run['reset_event_maps']
            software = Path(controls['cases'][case]['software'])
            assert run['loaded_acts_libraries'] == sorted(str(software/'lib'/name) for name in
                ['libFaserActsKalmanFilter.so', 'libFaserActsGeometry.so', 'libFaserActsGeometryLib.so'])
            for prop, factor, sigma in zip(['covLoc0', 'covLoc1', 'covPhi', 'covTheta', 'covQOverP'],
                    task['initial_variance_inflation'], [.9, .019, .033, .00029, .0001]):
                assert abs(run['seed_covariance'][prop]/(factor*sigma*sigma)-1) < 1e-12
            observed = json.loads((source/'native_purity.json').read_text())
            assert observed['counts']['source_events'] == 1000
            counts.update(observed['counts'])
            purity_counts.update(observed['purity_counts'])
            for key, value in observed['per_layer'].items():
                layers[key].update(value)
            with np.load(source/'event_counts.npz') as data:
                assert data['source_entry'].tolist() == list(range(1000))
                keys.extend((task['input'], int(e)) for e in data['source_entry'])
                events.append(data['counts'])
                tracks.append(data['tracks'])
            current = read_histograms(source/'native_3station_counts.root', case+str(task['file_index']))
            for key, hist in current.items():
                if key in histograms:
                    histograms[key].Add(hist)
                else:
                    histograms[key] = hist
            resources.append(dict(file_index=task['file_index'], hostname=run['hostname'],
                wall_s=run['wall_s'], **json.loads((source/'resources.json').read_text())))
        assert len(tasks) == 20 and len(keys) == len(set(keys)) == 20000
        event_data[case] = np.concatenate(events), np.concatenate(tracks), keys
        totals = event_data[case][0].sum(axis=0)
        assert event_data[case][1].sum(axis=0).tolist() == [counts[k] for k in
            ['writer_tracks', 'github_selected_tracks', 'selected_tracks']]
        reference = Path(controls['reference'])
        old_summary = json.loads((reference/'reductions'/(case+'_truth_summary.json')).read_text())
        for field in ['writer_tracks', 'github_selected_tracks']:
            assert counts[field] == old_summary['counts'][field]
        old_events, old_keys = [], []
        for task in json.loads((reference/(case+'_tasks.json')).read_text()):
            with np.load(Path(task['output'])/'event_counts.npz') as old:
                old_keys.extend((task['input'], int(e)) for e in old['source_entry'])
                old_events.append(old['counts'])
        assert old_keys == keys
        baseline_events = np.concatenate(old_events)
        purity_effects[case] = {f'smt_S{s}L{l}': {label: paired_rate(
            baseline_events[:, 2, li], event_data[case][0][:, 2, li], field)
            for label, field in [('tail_half', 1), ('tail_point8', 2)]}
            for li, (s, l) in enumerate(LAYERS)}
        flow_checks = []
        for (component, station, layer), hist in histograms.items():
            values = layers[f'{station}_{layer}']
            prefix = COMPONENTS[component]+'_native'
            assert hist.GetEntries() == values['states']
            visible_finite = values[prefix+'_finite']-values[prefix+'_underflow']-values[prefix+'_overflow']
            invalid_visible = hist.Integral()-visible_finite
            assert 0 <= invalid_visible <= values[prefix+'_invalid_reco']
            assert hist.GetBinContent(0) == values[prefix+'_underflow']
            assert hist.GetBinContent(101) == values['states']-values[prefix+'_finite']+values[prefix+'_overflow']-invalid_visible
            flow_checks.append(dict(histogram=histogram_name(component, station, layer),
                raw_entries=hist.GetEntries(), visible_invalid_reco=invalid_visible, closure=True))
            bins[f'{case}_{histogram_name(component, station, layer)}'] = dict(
                edges=[hist.GetBinLowEdge(i) for i in range(1, 102)],
                flow_counts=[int(hist.GetBinContent(i)) for i in range(102)],
                flow_sumw=[hist.GetBinContent(i) for i in range(102)],
                flow_sumw2=[hist.GetBinError(i)**2 for i in range(102)], raw_entries=hist.GetEntries())
        metrics = {}
        for si, stage in enumerate(STAGES):
            for li, (station, layer) in enumerate(LAYERS):
                value = totals[si, li]
                source = layers[f'{station}_{layer}']
                for fi, field in enumerate(['finite', 'below_minus_half', 'below_minus_point8']):
                    assert int(value[fi]) == source[stage+'_native_'+field]
                metrics[f'{stage}_S{station}L{layer}'] = dict(raw_states=source['states'],
                    finite=int(value[0]), tail_half=int(value[1]), tail_point8=int(value[2]),
                    tail_half_pct=float(100*value[1]/value[0]),
                    tail_point8_pct=float(100*value[2]/value[0]))
        with tempfile.TemporaryDirectory() as scratch:
            name = case+'_native_3station.root'
            output = ROOT.TFile(str(Path(scratch)/name), 'RECREATE')
            for key, hist in histograms.items():
                assert hist.Write(histogram_name(*key)) > 0
            assert not output.TestBit(ROOT.TFile.kWriteError)
            output.Close()
            shutil.copy2(Path(scratch)/name, destination/name)
        coverage[case] = dict(files=20, source_events=20000, entries_per_file=1000,
            source_ranges='each frozen file: [0,1000)', unique_file_entry_keys=20000)
        summaries[case] = dict(counts=dict(counts), layers=metrics, per_layer=layers,
            purity_counts=dict(purity_counts), resources=resources, histogram_flow_checks=flow_checks)
    assert event_data['beforeGeo'][2] == event_data['afterGeo'][2]
    comparisons = {}
    for si, stage in enumerate(STAGES):
        for li, (station, layer) in enumerate(LAYERS):
            comparisons[f'{stage}_S{station}L{layer}'] = {label: paired_rate(
                event_data['beforeGeo'][0][:, si, li], event_data['afterGeo'][0][:, si, li], field)
                for label, field in [('tail_half', 1), ('tail_point8', 2)]}
    summary = dict(cases=summaries, comparisons=comparisons, selection=controls['final_purity_cut'],
        truth=controls['truth'], isolate_smoothing=True, coverage=coverage,
        purity_effects_vs_isolation_only=purity_effects)
    (destination/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (destination/'bin_counts.json').write_text(json.dumps(bins, indent=2)+'\n')
    (destination/'histogram_schema.json').write_text(json.dumps(dict(
        histogram_count=54, components={'0': 'smoothed', '1': 'filtered', '2': 'predicted'},
        bins=100, range=[-1.2, 1.2], residual='p_reco/p_native-1 (dimensionless)',
        flow_index={'0': 'underflow', '1..100': 'visible bins', '101': 'overflow including native ROOT NaN behavior'},
        weights='unit weights; raw bin counts = sumw = sumw2',
        normalization='Main GitHub figures: raw Entries; supplemental tail figures: Entries / raw Entries'),
        indent=2)+'\n')
    with (destination/'layer_metrics.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=['case', 'stage_layer', *next(iter(summaries['beforeGeo']['layers'].values())).keys()])
        writer.writeheader()
        for case in CASES:
            for key, value in summaries[case]['layers'].items():
                writer.writerow(dict(case=case, stage_layer=key, **value))
    print(json.dumps({case: v['counts'] for case, v in summaries.items()}, indent=2))


def plot(root):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    summary = json.loads((root/'reductions/summary.json').read_text())
    figures = root/'figures'
    figures.mkdir(exist_ok=True)
    data = {case: read_histograms(root/'reductions'/(case+'_native_3station.root'), 'plot_'+case) for case in CASES}
    for panel, (station, layer) in enumerate(LAYERS):
        os.environ[f'REFERENCE_YMAX_{panel}'] = str(max(h[(c, station, layer)].GetMaximum()
            for h in data.values() for c in range(3)))
    ROOT.gROOT.ProcessLine('.L '+str(root/'reference/plot/ckf_residual_perLayer.C'))
    first_scope = summary['coverage']['beforeGeo']
    assert first_scope == summary['coverage']['afterGeo']
    entries_per_file = first_scope['source_events']//first_scope['files']
    scope = (f"110051 | {first_scope['files']} files, entries [0,{entries_per_file}) each"
             f" | {first_scope['source_events']:,} source events"
             ' | native SiHit-linked HepMC truth')
    cuts = 'Isolated smoothing | final purity>0.7 | #chi^{2}/NDF<20, NDF>9, r_{max}<90 mm'
    frames = 'CKF sensor states | FASER / rotated ACTS frames'
    index = []
    for case in CASES:
        prefix = figures/('110051_'+case+'_native_nine_layers')
        os.environ['REFERENCE_SCOPE_LABEL'] = '#splitline{'+scope+' | '+case+'}{'+cuts+' | '+frames+'}'
        ROOT.ckf_residual_perLayer(str(root/'reductions'/(case+'_native_3station.root')),
            str(prefix), case, 't_eQOP', True)
        index.append(dict(name=prefix.name, description=scope+' | '+case+' | '+cuts))
    prefix = figures/'110051_geometry_overlay'
    canvas, retained = overlay(data['beforeGeo'], data['afterGeo'], prefix,
        '#splitline{'+scope+' | beforeGeo vs afterGeo}{'+cuts+' | '+frames+'}', CASES)
    canvas.Close()
    index.append(dict(name=prefix.name, description='Original P/F/S colors; solid beforeGeo, dashed afterGeo; raw Entries'))
    for log in [False, True]:
        fig, axes = plt.subplots(3, 3, figsize=(15, 10), sharex=True)
        for ax, (station, layer) in zip(axes.flat, LAYERS):
            for case, color in zip(CASES, ['#222222', '#d62728']):
                h = data[case][(0, station, layer)]
                edges = [h.GetBinLowEdge(i) for i in range(1, 102)]
                values = np.array([h.GetBinContent(i) for i in range(1, 101)])/h.GetEntries()
                ax.stairs(values, edges, color=color, label=f'{case} (N={int(h.GetEntries()):,})')
            ax.set_title(f'Station {station}, layer {layer}')
            ax.set_xlim(-1.2, .3)
            ax.set_xlabel('(reco p - truth p) / truth p')
            ax.set_ylabel('Entries / raw Entries')
            ax.legend(fontsize=8, frameon=False)
            if log:
                ax.set_yscale('log')
                ax.set_ylim(1e-5, 1)
            else:
                ax.set_ylim(bottom=0)
        fig.suptitle(scope+'\nIsolated smoothing; final purity > 0.7; chi2/NDF < 20; NDF > 9; rmax < 90 mm\n'
            'Selected CKF smoothed sensor states; native FASER / rotated ACTS frames; normalization includes flow / NaN entries')
        fig.tight_layout()
        prefix = figures/('110051_smoothed_left_tail'+('_log' if log else ''))
        fig.savefig(str(prefix)+'.png', dpi=150)
        fig.savefig(str(prefix)+'.pdf')
        plt.close(fig)
        index.append(dict(name=prefix.name, description='Smoothed native residual, common bins, raw-Entries normalization'))
    sys.path.insert(0, str(root/'software_tools'))
    from pypdf import PdfReader, PdfWriter
    entries_checks = {}
    for case in CASES:
        text = PdfReader(figures/('110051_'+case+'_native_nine_layers.pdf')).pages[0].extract_text()
        displayed = list(map(int, re.findall(r'N\s*=\s*(\d+)', text)))
        expected = [int(data[case][(component, station, layer)].GetEntries())
            for station, layer in LAYERS for component in range(3)]
        assert displayed == expected, (case, displayed, expected)
        entries_checks[case] = dict(labels_checked=len(expected), all_equal_to_frozen_histograms=True)
    writer = PdfWriter()
    for item in index:
        reader = PdfReader(figures/(item['name']+'.pdf'))
        assert len(reader.pages) == 1
        writer.add_page(reader.pages[0])
    with (figures/'native_purity_comparison.pdf').open('wb') as stream:
        writer.write(stream)
    combined = PdfReader(figures/'native_purity_comparison.pdf')
    assert len(combined.pages) == len(index)
    for page, item in zip(combined.pages, index):
        original = PdfReader(figures/(item['name']+'.pdf')).pages[0]
        assert page.get_contents().get_data() == original.get_contents().get_data()
        assert list(page.mediabox) == list(original.mediabox)
        assert list(page.cropbox) == list(original.cropbox)
        assert page.rotation == original.rotation
    (figures/'figure_validation.json').write_text(json.dumps(dict(pages=len(index),
        contents_and_page_boxes_preserved=True, raw_entries_labels=entries_checks), indent=2)+'\n')
    (figures/'figure_index.json').write_text(json.dumps(index, indent=2)+'\n')
    (figures/'figure_index.md').write_text('# Native truth, isolated smoothing, final purity > 0.7\n\n'+
        '\n'.join(f"- [{x['name']}]({x['name']}.png): {x['description']}" for x in index)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['smoke', 'merge', 'plot', 'all'], default='all')
    args = parser.parse_args()
    if args.phase == 'smoke':
        validate_smoke(args.root)
    if args.phase in ['all', 'merge']:
        merge(args.root)
    if args.phase in ['all', 'plot']:
        plot(args.root)


if __name__ == '__main__':
    main()
