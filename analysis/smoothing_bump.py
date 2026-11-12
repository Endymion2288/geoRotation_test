#!/usr/bin/env python3
"""Compare saved Calypso candidate states before/after isolated smoothing."""
import argparse
from collections import Counter, defaultdict
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile

import numpy as np

from pilot import event_groups, process_event, read_records, truth_fiducial


LAYERS = [(s, l) for s in (1, 2, 3) for l in (0, 1, 2)]
STAGES = ('smoothed', 'filtered', 'predicted')
EDGES = np.linspace(-1.2, 1.2, 101)


def residual(step, stage, momentum):
    bound = step.get(stage)
    qop = bound[4] if bound else None
    return 1. / abs(qop) / momentum - 1. if qop and math.isfinite(qop) else math.nan


def histogram(values):
    values = np.asarray(values, dtype=float)
    valid = values[np.isfinite(values)]
    bins = np.bincount(np.searchsorted(EDGES, valid, side='right'), minlength=102)
    return dict(flow_bins=bins.tolist(), finite=len(valid), invalid=len(values)-len(valid),
                tail_half=int(np.count_nonzero(valid < -.5)),
                bump_fifth=int(np.count_nonzero(valid < -.8)))


def reduce(task, output):
    counts = Counter()
    values = defaultdict(list)
    pairs = {f'{s}_{l}': Counter() for s, l in LAYERS}
    for entry, records in enumerate(event_groups(read_records(Path(task['trace'])))):
        counts['source_events'] += 1
        truth = {r['barcode']: r for r in records if r['stage'] == 'truth'}
        lookup = {(r['stage'], r['seed'], r['candidate']): r for r in records}
        particles, _ = process_event((records[0]['run'], records[0]['lumi'], records[0]['event']),
                                     records, '110051')
        final = {(p['final']['seed'], p['final']['candidate']): p for p in particles
                 if p['fiducial'] and p['truth_H18'] and p['final']}
        for meta in (r for r in records if r['stage'] == 'selected_shared'):
            identity = (meta['seed'], meta['candidate'])
            before = lookup[('selected_shared_state_trace', *identity)]['steps']
            after = lookup[('selected_state_trace', *identity)]['steps']
            assert len(before) == len(after)
            counts['selected_candidates'] += 1
            barcode = meta['matches'][0][0] if meta['matches'] else None
            particle = truth.get(barcode)
            cohorts = {}
            if particle:
                cohorts['all_truth'] = particle
                if particle['signal'] and particle['status'] == 1 and abs(particle['pdg']) == 11:
                    cohorts['signal'] = particle
                    if truth_fiducial(particle) and sum(particle['hit_counts'][1:]) >= 18:
                        cohorts['fiducial_H18'] = particle
                        if meta['matches'][0][1] / len(meta['hits']) >= .7:
                            cohorts['purity70'] = particle
            else:
                counts['missing_selected_majority_truth'] += 1
            chosen = final.get(identity)
            if chosen:
                cohorts['final_unique'] = truth[chosen['barcode']]
                if chosen['final']['good95']:
                    cohorts['final_good95'] = truth[chosen['barcode']]
                    counts['old_good95_candidates'] += 1
                    if chosen['barcode'] == barcode:
                        cohorts['final_good95_same_majority'] = truth[barcode]
                    else:
                        counts['old_good95_majority_mismatch'] += 1
            for name, particle in cohorts.items():
                counts[name+'_candidates'] += 1
                momentum = math.sqrt(sum(x*x for x in particle['momentum_F_MeV'])) * .001
                for left, right in zip(before, after):
                    assert (left['hit'], left['measurement_flag']) == (right['hit'], right['measurement_flag'])
                    if not left['measurement_flag']:
                        continue
                    station, layer = (left['hit'] >> 59) & 3, (left['hit'] >> 57) & 3
                    assert (station, layer) in LAYERS
                    assert left.get('predicted') == right.get('predicted')
                    assert left.get('filtered') == right.get('filtered')
                    key = f'{station}_{layer}'
                    for mode, step in [('shared', left), ('isolated', right)]:
                        for stage in STAGES:
                            values[(name, mode, key, stage)].append(residual(step, stage, momentum))
                    if name == 'final_good95':
                        a, b = (residual(step, 'smoothed', momentum) for step in (left, right))
                        p = pairs[key]
                        p['states'] += 1
                        p['smoothed_qop_changed'] += left.get('smoothed', [None]*5)[4] != right.get('smoothed', [None]*5)[4]
                        p['shared_smoother_nonclosure'] += abs(left.get('smoothing_identity_error', 0.)) * momentum > 1e-6
                        p['isolated_smoother_nonclosure'] += abs(right.get('smoothing_identity_error', 0.)) * momentum > 1e-6
                        if math.isfinite(a) and math.isfinite(b):
                            p['common_finite'] += 1
                            for label, threshold in [('tail', -.5), ('bump', -.8)]:
                                p[label+'_recovered'] += a < threshold <= b
                                p[label+'_introduced'] += b < threshold <= a
        if (entry+1) % 250 == 0:
            print('Reduced', task['file_index'], entry+1, flush=True)
    assert counts['source_events'] == task['events']
    result = dict(task=task, counts=dict(counts), pairs={k: dict(v) for k, v in pairs.items()},
                  histograms={'|'.join(k): histogram(v) for k, v in values.items()},
                  bin_edges=EDGES.tolist(),
                  definition='Same saved Calypso candidate, hits and truth; shared vs candidate-isolated smoothing')
    (output/'counts.json').write_text(json.dumps(result, indent=2)+'\n')


def merge(root):
    tasks = json.loads((root/'tasks.json').read_text())
    counts = Counter()
    pairs = {f'{s}_{l}': Counter() for s, l in LAYERS}
    histograms = {}
    for task in tasks:
        folder = Path(task['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        value = json.loads((folder/'counts.json').read_text())
        assert value['task'] == task['source']
        counts.update(value['counts'])
        for key, row in value['pairs'].items():
            pairs[key].update(row)
        for key, row in value['histograms'].items():
            if key not in histograms:
                histograms[key] = dict(flow_bins=[0]*102, finite=0, invalid=0, tail_half=0, bump_fifth=0)
            dest = histograms[key]
            dest['flow_bins'] = [a+b for a, b in zip(dest['flow_bins'], row['flow_bins'])]
            for name in ['finite', 'invalid', 'tail_half', 'bump_fifth']:
                dest[name] += row[name]
    old = json.loads((root/'manifest.json').read_text())
    reference = json.loads(Path(old['old_histograms']).read_text())['variants']['strict']
    assert counts['old_good95_candidates'] == reference['counts']['good95_particles']
    for station, layer in LAYERS:
        for stage in STAGES:
            original = reference['layers']['good95'][station-1][layer][stage]
            actual = histograms[f'final_good95|isolated|{station}_{layer}|{stage}']
            assert actual['flow_bins'][1:101] == original['counts'], (station, layer, stage)
            assert (actual['finite'], actual['tail_half'], actual['invalid']) == (original['n'], original['tail_half'], original['invalid'])
    result = dict(counts=dict(counts), pairs={k: dict(v) for k, v in pairs.items()},
                  histograms=histograms, old_figure_27_histograms_exactly_reproduced=True,
                  bin_edges=EDGES.tolist(), scope=old)
    native_counts = Counter()
    native_layers = {f'{s}_{l}': Counter() for s, l in LAYERS}
    for file_index in old['file_indices']:
        folder = Path(old['native_root'])/'shards/afterGeo_strict'/f'f{file_index:02d}'
        source = json.loads((folder/'truth_comparison.json').read_text())
        native_counts.update(source['counts'])
        for key, row in source['per_layer'].items():
            native_layers[key].update(row)
    result['native_same_files_counts'] = dict(native_counts)
    result['native_same_files_per_layer'] = {k: dict(v) for k, v in native_layers.items()}
    (root/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def render(root, summary, destination):
    import ROOT
    from reference_figures import histogram_name, overlay, worker_histogram
    ROOT.gROOT.SetBatch(True)
    settings = summary['scope']
    data = {}
    for mode in ['shared', 'isolated']:
        histograms = {}
        for component, stage in enumerate(STAGES):
            for station, layer in LAYERS:
                name = histogram_name(component, station, layer)
                h = ROOT.TH1F(mode+'_'+name, '', 100, -1.2, 1.2)
                h.SetDirectory(0)
                row = summary['histograms'][f'final_good95|{mode}|{station}_{layer}|{stage}']
                for index, value in enumerate(row['flow_bins']):
                    h.SetBinContent(index, value)
                h.SetEntries(row['finite'])
                histograms[(component, station, layer)] = h
        data[mode] = histograms
    native = {}
    for index in settings['file_indices']:
        path = Path(settings['native_root'])/'shards/afterGeo_strict'/f'f{index:02d}'/'production_3station_counts.root'
        file = ROOT.TFile.Open(str(path))
        for component in range(3):
            for station, layer in LAYERS:
                h, _ = worker_histogram(file, component, station, layer)
                key = component, station, layer
                if key not in native:
                    native[key] = h.Clone('native_'+histogram_name(*key))
                    native[key].SetDirectory(0)
                else:
                    native[key].Add(h)
        file.Close()
    data['native'] = native
    macro = Path(settings['native_root'])/'reference/adapted/ckf_residual_perLayer.C'
    ROOT.gROOT.ProcessLine('.L '+str(macro))
    for panel, key in enumerate(LAYERS):
        os.environ[f'REFERENCE_YMAX_{panel}'] = str(max(h[(c, *key)].GetMaximum() for h in data.values() for c in range(3)))
    index = []
    for mode, title in [('native', 'Native CKF, GitHub cuts'), ('shared', 'S04 selection, original smoothed slots'),
                        ('isolated', 'S04 selection, candidate-isolated smoothing')]:
        prefix = destination/('110051_'+mode+'_nine_layers')
        counts_path = destination/(prefix.name+'_3station_counts.root')
        file = ROOT.TFile(str(counts_path), 'RECREATE')
        for key, hist in data[mode].items():
            hist.Write(histogram_name(*key))
        file.Close()
        os.environ['REFERENCE_SCOPE_LABEL'] = '110051 | Same 10 files, 10,000 events | Production truth | '+title
        ROOT.ckf_residual_perLayer(str(counts_path), str(prefix), mode, 't_eQOP', True)
        index.append(prefix.name)
    for first, second, name, text in [
            ('native', 'shared', 'selection_comparison', 'Native GitHub selection vs S04 selection; original smoothed slots'),
            ('shared', 'isolated', 'smoothing_only', 'Same tracks, hits, truth and cuts; only isolated smoothing changes')]:
        prefix = destination/('110051_'+name)
        canvas, retained = overlay(data[first], data[second], prefix,
            '110051 | Same 10 files, 10,000 events | '+text, [first, second])
        canvas.Close()
        index.append(prefix.name)
    sys.path.insert(0, str(Path(settings['native_root'])/'software_tools'))
    from pypdf import PdfReader, PdfWriter
    pdf = PdfWriter()
    for name in index:
        pdf.add_page(PdfReader(destination/(name+'.pdf')).pages[0])
    with (destination/'smoothing_bump_comparison.pdf').open('wb') as stream:
        pdf.write(stream)
    (destination/'figure_index.json').write_text(json.dumps(index, indent=2)+'\n')


def plot(root, summary):
    destination = root/'figures'
    destination.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='calypso_bump_figures_') as scratch:
        render(root, summary, Path(scratch))
        for path in Path(scratch).iterdir():
            shutil.copy2(path, destination/path.name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', choices=['reduce', 'merge', 'plot'], required=True)
    parser.add_argument('--task', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.phase == 'reduce':
        reduce(json.loads(args.task.read_text()), args.output)
    elif args.phase == 'merge':
        plot(args.output, merge(args.output))
    else:
        plot(args.output, json.loads((args.output/'summary.json').read_text()))


if __name__ == '__main__':
    main()
