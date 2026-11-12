#!/usr/bin/env python3
"""Measure legacy selection effects on actual archived Calypso isolated states."""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path

import numpy as np

from pilot import candidate_rank, event_groups, pnorm, read_records, truth_fiducial
from smoothing_bump import EDGES, LAYERS, STAGES
from smoothing_covariance import paired_rate

CUTS = ('signal', 'fiducial', 'H18', 'purity', 'unique', 'chi2', 'layers', 'ndof', 'radius')
QUALITY = ('chi2', 'layers', 'ndof', 'radius')
FULL = set(CUTS)
CUMULATIVE = [
    ('base', set()), ('signal', {'signal'}),
    ('fiducial', {'signal', 'fiducial'}),
    ('H18', {'signal', 'fiducial', 'H18'}),
    ('purity', {'signal', 'fiducial', 'H18', 'purity'}),
    ('unique', {'signal', 'fiducial', 'H18', 'purity', 'unique'}),
    ('good95', FULL),
]
SPECS = dict(CUMULATIVE)
SPECS.update({f'minus_{cut}': FULL-{cut} for cut in CUTS})
SPECS['minus_good95'] = FULL-set(QUALITY)
SPECS.update({f'only_{cut}': {cut} for cut in CUTS})
SPECS['only_good95'] = set(QUALITY)
COHORTS = ['selected_all', 'github', 'github_selected_purity', 'original_figure']
COHORTS += [prefix+name for prefix in ('legacy_', 'github_') for name in SPECS]
MODES = ('shared', 'isolated')


def identity(row):
    return row['seed'], row['candidate']


def final_selection(finals, truth, enabled):
    """Recompute matching/ranking for each ablation; rank before quality as in pilot."""
    eligible = []
    for row in finals:
        matches = row.get('matches', [])
        if not matches or pnorm(row) is None:
            continue
        particle = truth.get(matches[0][0])
        if particle is None:
            continue
        tests = dict(
            signal=particle['signal'] and particle['status'] == 1 and abs(particle['pdg']) == 11,
            fiducial=truth_fiducial(particle),
            H18=sum(particle['hit_counts'][1:]) >= 18,
            purity=bool(row['hits']) and matches[0][1]/len(row['hits']) >= .7)
        if not all(tests[cut] for cut in enabled & tests.keys()):
            continue
        eligible.append(row)
    if 'unique' in enabled:
        groups = defaultdict(list)
        for row in eligible:
            groups[row['matches'][0][0]].append(row)
        eligible = [max(group, key=candidate_rank) for group in groups.values()]
    result = set()
    for row in eligible:
        ndof = row.get('ndof', 0)
        chi2 = row.get('chi2')
        tests = dict(chi2=ndof > 0 and chi2 is not None and math.isfinite(chi2) and chi2/ndof < 25,
            layers=row.get('layers', 0) >= 7, ndof=ndof+5 >= 12,
            radius=row.get('rmax_mm', math.inf) < 95)
        if all(tests[cut] for cut in enabled & tests.keys()):
            result.add(identity(row))
    return result


def candidate_arrays(before, after, momentum):
    """One candidate's P/F/S histograms and event-cluster sufficient counts."""
    bins = np.zeros((2, 3, 9, 102), dtype=np.int64)
    metrics = np.zeros((2, 3, 9, 4), dtype=np.int64)
    event = np.zeros((2, 9, 5), dtype=np.int32)
    for mode_index, steps in enumerate((before, after)):
        layers = np.array([LAYERS.index(((r['hit'] >> 59) & 3, (r['hit'] >> 57) & 3)) for r in steps])
        for stage_index, stage in enumerate(STAGES):
            qop = np.array([r.get(stage, [math.nan]*5)[4] for r in steps], dtype=float)
            with np.errstate(divide='ignore', invalid='ignore'):
                residual = 1/np.abs(qop)/momentum-1
            valid = np.isfinite(residual)
            bin_index = np.searchsorted(EDGES, residual, side='right')
            np.add.at(bins[mode_index, stage_index], (layers[valid], bin_index[valid]), 1)
            for layer in range(9):
                take = layers == layer
                finite = take & valid
                n = int(finite.sum())
                metrics[mode_index, stage_index, layer] = [n, int((take & ~valid).sum()),
                    int((finite & (residual < -.5)).sum()), int((finite & (residual < -.8)).sum())]
                if stage_index == 0:
                    event[mode_index, layer] = [n, metrics[mode_index, 0, layer, 2],
                        metrics[mode_index, 0, layer, 3],
                        int((finite & (bin_index >= 10) & (bin_index < 18)).sum()),
                        int((finite & (bin_index >= 26) & (bin_index < 34)).sum())]
    return bins, metrics, event


def reduce(task, output):
    geometry = json.loads(Path(task['geometry']).read_text())
    assert geometry['geometry'] == 'FASERNU-04' and geometry['conditions'] == 'OFLCOND-FASER-05'
    sensors = {row['id']: np.array(row['local_to_F']) for row in geometry['tracker']}
    hist = np.zeros((len(COHORTS), 2, 3, 9, 102), dtype=np.int64)
    metrics = np.zeros((len(COHORTS), 2, 3, 9, 4), dtype=np.int64)
    events = np.zeros((task['events'], len(COHORTS), 2, 9, 5), dtype=np.int32)
    tracks = np.zeros((task['events'], len(COHORTS)), dtype=np.int32)
    counters = Counter()
    geometry_ids = {}
    max_z_error = 0.
    for entry, records in enumerate(event_groups(read_records(Path(task['trace'])))):
        counters['source_events'] += 1
        truth = {r['barcode']: r for r in records if r['stage'] == 'truth'}
        finals = [r for r in records if r['stage'] == 'final']
        final_by_id = {identity(row): row for row in finals}
        selected = [r for r in records if r['stage'] == 'selected']
        selected_by_id = {identity(row): row for row in selected}
        lookup = {(r['stage'], *identity(r)): r for r in records if r['stage'].endswith('_state_trace')}
        selections = {name: final_selection(finals, truth, enabled) for name, enabled in SPECS.items()}
        github = set()
        pairs = {}
        for metadata in selected:
            key = identity(metadata)
            before = [s for s in lookup[('selected_shared_state_trace', *key)]['steps'] if s['measurement_flag']]
            after = [s for s in lookup[('selected_state_trace', *key)]['steps'] if s['measurement_flag']]
            assert len(before) == len(after) == len(metadata['hits'])
            positions = []
            for left, right in zip(before, after):
                assert (left['hit'], left['geometry_id']) == (right['hit'], right['geometry_id'])
                assert left['predicted'] == right['predicted'] and left['filtered'] == right['filtered']
                wafer = right['hit'] >> 53 << 53
                transform = sensors[wafer]
                max_z_error = max(max_z_error, abs(transform[2, 3]-right['z_F_mm']))
                assert abs(transform[2, 3]-right['z_F_mm']) < 1e-8
                old_id = geometry_ids.setdefault(wafer, right['geometry_id'])
                assert old_id == right['geometry_id']
                loc = right['smoothed'][:2]
                position = transform[:2, :2] @ loc + transform[:2, 3]
                positions.append(position)
            rmax = float(np.max(np.linalg.norm(positions, axis=1)))
            ndof = metadata['ndof']
            if ndof > 9 and metadata['chi2']/ndof < 20 and rmax < 90:
                github.add(key)
            pairs[key] = before, after
        sets = dict(selected_all=set(selected_by_id), github=github,
            github_selected_purity={key for key in github if selected_by_id[key]['matches'] and
                selected_by_id[key]['matches'][0][1]/len(selected_by_id[key]['hits']) >= .7},
            original_figure=selections['good95'])
        for name, chosen in selections.items():
            sets['legacy_'+name] = chosen
            sets['github_'+name] = github & chosen
        for metadata in selected:
            key = identity(metadata)
            counters['selected_candidates'] += 1
            majority = metadata['matches'][0][0] if metadata['matches'] else None
            particle = truth.get(majority)
            final = final_by_id.get(key)
            if final and final['matches']:
                counters['selected_final_majority_mismatches'] += majority != final['matches'][0][0]
            if particle:
                momentum = math.sqrt(sum(x*x for x in particle['momentum_F_MeV']))*.001
                data = candidate_arrays(*pairs[key], momentum)
                indices = [i for i, name in enumerate(COHORTS) if name != 'original_figure' and key in sets[name]]
                hist[indices] += data[0]
                metrics[indices] += data[1]
                events[entry, indices] += data[2]
                tracks[entry, indices] += 1
            else:
                counters['missing_selected_majority_truth'] += 1
            if key in sets['original_figure']:
                legacy_truth = truth[final['matches'][0][0]]
                momentum = math.sqrt(sum(x*x for x in legacy_truth['momentum_F_MeV']))*.001
                data = candidate_arrays(*pairs[key], momentum)
                i = COHORTS.index('original_figure')
                hist[i] += data[0]
                metrics[i] += data[1]
                events[entry, i] += data[2]
                tracks[entry, i] += 1
                counters['old_figure_candidates'] += 1
                counters['old_figure_majority_mismatches'] += majority != final['matches'][0][0]
        if (entry+1) % 250 == 0:
            print('Reduced', task['file_index'], entry+1, flush=True)
    assert counters['source_events'] == task['events']
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'counts.npz', hist=hist, metrics=metrics, events=events, tracks=tracks)
    (output/'counts.json').write_text(json.dumps(dict(task=task, counts=dict(counters),
        cohorts=COHORTS, modes=MODES, stages=STAGES, max_sensor_z_error_mm=max_z_error,
        geometry_ids={str(k): v for k, v in geometry_ids.items()}), indent=2)+'\n')


def ratio_interval(values, numerator, denominator):
    n = values[:, denominator].sum()
    p = values[:, numerator].sum()/n
    influence = (values[:, numerator]-p*values[:, denominator])/n
    error = np.sqrt(len(values)/(len(values)-1)*np.square(influence).sum())
    return dict(value=float(p), lower95=float(p-1.96*error), upper95=float(p+1.96*error))


def merge(root):
    tasks = json.loads((root/'tasks.json').read_text())
    total_hist = total_metrics = None
    event_parts, track_parts, keys = [], [], []
    counters = Counter()
    geometry_ids = {}
    for task in tasks:
        source = Path(task['output'])
        assert (source/'exit_code.txt').read_text().strip() == '0'
        info = json.loads((source/'counts.json').read_text())
        assert info['cohorts'] == COHORTS and info['task'] == task['source']
        counters.update(info['counts'])
        for wafer, gid in info['geometry_ids'].items():
            assert geometry_ids.setdefault(wafer, gid) == gid
        keys.extend((task['source']['input'], entry) for entry in range(task['source']['events']))
        with np.load(source/'counts.npz') as data:
            if total_hist is None:
                total_hist, total_metrics = data['hist'].copy(), data['metrics'].copy()
            else:
                total_hist += data['hist']
                total_metrics += data['metrics']
            event_parts.append(data['events']); track_parts.append(data['tracks'])
    assert len(keys) == len(set(keys)) == 10000
    events, tracks = np.concatenate(event_parts), np.concatenate(track_parts)
    original = json.loads((Path(json.loads((root/'manifest.json').read_text())['source_bump'])/'summary.json').read_text())
    for ci, old_name in [(COHORTS.index('original_figure'), 'final_good95'), (COHORTS.index('selected_all'), 'all_truth')]:
        for mi, mode in enumerate(MODES):
            for si, stage in enumerate(STAGES):
                for li, (station, layer) in enumerate(LAYERS):
                    old = original['histograms'][f'{old_name}|{mode}|{station}_{layer}|{stage}']
                    assert total_hist[ci, mi, si, li].tolist() == old['flow_bins']
                    assert total_metrics[ci, mi, si, li].tolist() == [old['finite'], old['invalid'], old['tail_half'], old['bump_fifth']]
    validation = json.loads((root/'geometry_check/writer_smoothing_check.json').read_text())
    assert validation['sensor_geometry_matches_actual_acts_surfaces']
    assert validation['maximum_sensor_transform_error'] < 1e-10
    summaries = {}
    for ci, cohort in enumerate(COHORTS):
        summaries[cohort] = dict(candidates=int(tracks[:, ci].sum()), modes={})
        for mi, mode in enumerate(MODES):
            layers = {}
            for li, (station, layer) in enumerate(LAYERS):
                data = events[:, ci, mi, li]
                assert np.array_equal(data[:, :3].sum(axis=0), total_metrics[ci, mi, 0, li, [0, 2, 3]])
                assert data[:, 3].sum() == total_hist[ci, mi, 0, li, 10:18].sum()
                assert data[:, 4].sum() == total_hist[ci, mi, 0, li, 26:34].sum()
                layers[f'S{station}L{layer}'] = dict(finite=int(data[:, 0].sum()),
                    tail_half=ratio_interval(data, 1, 0), extreme=ratio_interval(data, 2, 0),
                    left_rise=ratio_interval(data, 3, 4),
                    window_counts=data[:, 3:5].sum(axis=0).tolist())
            summaries[cohort]['modes'][mode] = layers
    comparisons = {}
    pairs = [('github', 'github_'+name) for name, _ in CUMULATIVE]
    pairs += [('github_'+first, 'github_'+second) for (first, _), (second, _) in zip(CUMULATIVE, CUMULATIVE[1:])]
    pairs += [(prefix+'good95', prefix+name) for prefix in ('legacy_', 'github_') for name in SPECS if name.startswith('minus_')]
    pairs += [('github', 'github_'+name) for name in SPECS if name.startswith('only_')]
    for first, second in pairs:
        a, b = COHORTS.index(first), COHORTS.index(second)
        comparisons[first+'_vs_'+second] = {f'S{s}L{l}': {
            label: paired_rate(events[:, a, 1, li], events[:, b, 1, li], field)
            for label, field in [('tail_half', 1), ('extreme', 2)]}
            for li, (s, l) in enumerate(LAYERS)}
    result = dict(counts=dict(counters), source_files=10, source_events=10000,
        input_entries=keys, cohorts=summaries, comparisons=comparisons,
        cuts={name: sorted(enabled) for name, enabled in SPECS.items()},
        exact_original_figure_54_histograms=True, selected_all_54_histograms_reproduced=True,
        sensor_geometry_validation=validation, sensor_wafer_count=len(geometry_ids),
        left_rise_definition='equal 8-bin windows: residual [-0.984,-0.792) / [-0.600,-0.408); ratio > 1 means an excess in the left window',
        truth='Selected-majority production truth fixed per candidate for all ablations; original_figure alone uses exact legacy final-majority truth',
        selection='GitHub cuts evaluated on ON smoothed positions from actual sensor transforms; freeze passing candidates for shared/isolated comparison')
    (root/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    np.savez_compressed(root/'histograms.npz', hist=total_hist, metrics=total_metrics)
    # Per-source-event sufficient counts support paired uncertainties without state storage.
    np.savez_compressed(root/'event_counts.npz', events=events, tracks=tracks)
    print(json.dumps(dict(complete=True, candidates={name: summaries[name]['candidates'] for name in
        ['github', 'github_purity', 'github_unique', 'github_good95', 'legacy_good95', 'original_figure']})), flush=True)
    return result


def plot(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import ROOT
    ROOT.gROOT.SetBatch(True)
    from reference_figures import overlay
    import sys
    sys.path.insert(0, str(Path('/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_smoothing_covariance/20261112_full110051/software_tools')))
    from pypdf import PdfReader, PdfWriter
    summary = json.loads((root/'summary.json').read_text())
    with np.load(root/'histograms.npz') as source:
        hist = source['hist']; metrics = source['metrics']
    figures = root/'figures'; figures.mkdir(exist_ok=True)
    index = []
    prefix_title = '110051 | exact legacy sample: 10 files, 10,000 source entries; each file entries 0-999\n'
    suffix_title = '\nActual Calypso sensor states; rotated ACTS frame; production truth; finite-residual normalization; 100 source bins'
    def shape(names, labels, name, mode='isolated'):
        fig, axes = plt.subplots(3, 3, figsize=(15, 10), sharex=True)
        colors = ['#222222', '#e69f00', '#009e73', '#0072b2', '#cc79a7']
        for li, (axis, (station, layer)) in enumerate(zip(axes.flat, LAYERS)):
            ns = []
            for j, cohort in enumerate(names):
                ci = COHORTS.index(cohort); mi = MODES.index(mode)
                values = hist[ci, mi, 0, li]
                n = metrics[ci, mi, 0, li, 0]; ns.append(int(n))
                axis.stairs(values[1:101]/n, EDGES, label=labels[j], color=colors[j], linewidth=1.5)
            axis.set_title(f'Station {station}, layer {layer}')
            axis.set_xlim(-1.2, 1.2); axis.set_ylim(bottom=0)
            axis.set_xlabel('(reco p - truth p) / truth p'); axis.set_ylabel('Entries / finite residuals')
            axis.text(.02, .96, 'N: '+' / '.join(f'{n:,}' for n in ns), transform=axis.transAxes, va='top', fontsize=8)
            axis.axvspan(-.984, -.792, color='.8', alpha=.2)
            axis.axvspan(-.600, -.408, color='#0072b2', alpha=.07)
        handles, labels_ = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels_, loc='upper center', bbox_to_anchor=(.5, .915), ncol=2, frameon=False)
        fig.suptitle(prefix_title+'Candidate-isolated smoothing ON; '+('GitHub cuts + legacy selections' if names[0].startswith('github') else 'legacy selection ablations')+suffix_title, fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, .87))
        fig.savefig(figures/(name+'.png'), dpi=150); fig.savefig(figures/(name+'.pdf')); plt.close(fig)
        index.append(name)
    shape(['github', 'github_H18', 'github_purity', 'github_good95'],
        ['GitHub only', '+ signal, fiducial, H18', '+ final SDO purity >= 0.7', '+ one per truth and final Good95'], '110051_cumulative_isolated')
    shape(['github', 'github_only_purity', 'github_only_good95', 'github_good95'],
        ['GitHub only', 'GitHub + final purity only', 'GitHub + final Good95 only', 'GitHub + all legacy cuts'], '110051_individual_isolated')
    shape(['github_good95', 'github_minus_purity', 'github_minus_good95'],
        ['All legacy cuts', 'Remove final purity', 'Remove final Good95'], '110051_remove_one_isolated')
    # Exact reproduction of the image supplied by the user, including legacy final truth.
    ci = COHORTS.index('original_figure')
    data = []
    for mi, mode in enumerate(MODES):
        group = {}
        for si, stage in enumerate(STAGES):
            for li, (station, layer) in enumerate(LAYERS):
                h = ROOT.TH1D(f'{mode}_{si}_{station}_{layer}', '', 100, -1.2, 1.2); h.SetDirectory(0)
                for bi, value in enumerate(hist[ci, mi, si, li]):
                    h.SetBinContent(bi, float(value)); h.SetBinError(bi, math.sqrt(value))
                h.SetEntries(float(metrics[ci, mi, si, li, 0])); group[(si, station, layer)] = h
        data.append(group)
    name = '110051_uploaded_figure_reproduced'
    canvas, retained = overlay(*data, figures/name,
        '#splitline{110051 | Exact Image 1 reproduction | 10 files | 10,000 source entries | S04 full legacy selection}{Same candidates, hits, truth and cuts | shared / isolated smoothing | sensor states | rotated ACTS frame | production truth | finite Entries}',
        ['shared', 'isolated'])
    canvas.Close(); index.append(name)
    # Remove one condition from the complete selection. Effects are not additive.
    names = ['minus_'+cut for cut in CUTS]+['minus_good95']
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for axis, field, title in zip(axes, ['tail_half', 'extreme'], ['Residual < -0.5', 'Residual < -0.8']):
        values = np.array([[100*(summary['cohorts']['github_'+name]['modes']['isolated'][f'S{s}L{l}'][field]['value']-
            summary['cohorts']['github_good95']['modes']['isolated'][f'S{s}L{l}'][field]['value']) for name in names] for s,l in LAYERS])
        image = axis.imshow(values, cmap='RdBu_r', vmin=-max(1., abs(values).max()), vmax=max(1., abs(values).max()), aspect='auto')
        axis.set_xticks(range(len(names)), [n.removeprefix('minus_') for n in names], rotation=60, ha='right')
        axis.set_yticks(range(9), [f'S{s} L{l}' for s,l in LAYERS]); axis.set_title(title)
        for i in range(9):
            for j in range(len(names)):
                axis.text(j, i, f'{values[i,j]:+.2f}', ha='center', va='center', fontsize=8)
        fig.colorbar(image, ax=axis, label='Tail-rate increase when condition is removed [percentage points]')
    fig.suptitle(prefix_title+'GitHub cuts fixed; isolated smoothing; remove one legacy condition from the complete selection'+suffix_title, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, .88))
    name='110051_leave_one_out'; fig.savefig(figures/(name+'.png'),dpi=150); fig.savefig(figures/(name+'.pdf')); plt.close(fig); index.append(name)
    # Paired source-event confidence intervals for tail rates and left-window excess.
    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    names = ['github']+['github_'+name for name,_ in CUMULATIVE]
    labels = ['GitHub', '+final valid', '+signal', '+fiducial', '+H18', '+purity', '+one/truth', '+Good95']
    for col, layer in enumerate(['S1L0', 'S3L2']):
        for row, field in enumerate(['tail_half', 'left_rise']):
            axis = axes[row,col]
            values = [summary['cohorts'][name]['modes']['isolated'][layer][field] for name in names]
            scale = 100 if row == 0 else 1
            y=np.array([v['value'] for v in values])*scale
            errors=np.array([[v['value']-v['lower95'] for v in values], [v['upper95']-v['value'] for v in values]])*scale
            axis.errorbar(range(len(names)), y, yerr=errors, fmt='o-', capsize=3)
            axis.set_xticks(range(len(names)), labels, rotation=35, ha='right')
            axis.set_title(layer); axis.set_ylabel('Tail rate [%]' if row == 0 else 'Left-window / shoulder-window entries')
            if row == 1: axis.axhline(1, color='.6', linestyle='--')
            axis.grid(axis='y', alpha=.2)
    fig.suptitle(prefix_title+'Isolated smoothing ON; cumulative legacy cuts; 95% source-event cluster delta intervals\n'
        'Left rise: [-0.984,-0.792) / [-0.600,-0.408); equal-width windows; same candidate-majority truth throughout',fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .88))
    name='110051_cut_rates_and_left_rise';fig.savefig(figures/(name+'.png'),dpi=150);fig.savefig(figures/(name+'.pdf'));plt.close(fig);index.append(name)
    writer=PdfWriter()
    for name in index: writer.add_page(PdfReader(figures/(name+'.pdf')).pages[0])
    with (figures/'selection_ablation.pdf').open('wb') as stream:writer.write(stream)
    combined=PdfReader(figures/'selection_ablation.pdf')
    for page,name in zip(combined.pages,index):
        original=PdfReader(figures/(name+'.pdf')).pages[0]
        assert page.get_contents().get_data()==original.get_contents().get_data()
        for box in ['mediabox','cropbox','artbox','bleedbox','trimbox']:
            assert list(getattr(page,box))==list(getattr(original,box))
        assert page.rotation==original.rotation
    from PIL import Image
    for name in index:
        with Image.open(figures/(name+'.png')) as im:im.verify()
    (root/'figure_index.json').write_text(json.dumps(index,indent=2)+'\n')
    (root/'validation.json').write_text(json.dumps(dict(pdf_pages=len(index),png_integrity=True,
        pdf_contents_boxes_rotation_equal=True,exact_uploaded_figure_histograms=True,
        source_events=summary['source_events'],sensor_geometry_validated=True),indent=2)+'\n')


def export_tables(root):
    """Export compact frozen results; unit-weight counts are also sumw and sumw2."""
    summary = json.loads((root/'summary.json').read_text())
    with np.load(root/'histograms.npz') as source:
        metrics = source['metrics']
    (root/'histogram_schema.json').write_text(json.dumps(dict(
        file='histograms.npz', hist_axes=['cohort', 'mode', 'stage', 'layer', 'flow_bin'],
        cohorts=COHORTS, modes=MODES, stages=STAGES, layers=LAYERS, edges=EDGES.tolist(),
        underflow_index=0, overflow_index=101, weights='unit weights; hist stores identical raw counts, sumw and sumw2',
        metrics_axes=['cohort', 'mode', 'stage', 'layer', 'field'],
        metrics_fields=['finite', 'invalid', 'tail_half', 'extreme'],
        event_counts_file='event_counts.npz', event_fields=['finite', 'tail_half', 'extreme', 'left_window', 'shoulder_window']), indent=2)+'\n')
    fields = ['cohort', 'mode', 'layer', 'candidates', 'finite', 'invalid', 'tail_half_count',
        'extreme_count', 'tail_half_pct', 'extreme_pct', 'left_rise', 'left_window', 'shoulder_window']
    with (root/'selection_metrics.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for ci, name in enumerate(COHORTS):
            for mi, mode in enumerate(MODES):
                for li, (station, layer) in enumerate(LAYERS):
                    key = f'S{station}L{layer}'
                    item = summary['cohorts'][name]['modes'][mode][key]
                    writer.writerow(dict(cohort=name, mode=mode, layer=key,
                        candidates=summary['cohorts'][name]['candidates'], finite=item['finite'],
                        invalid=int(metrics[ci, mi, 0, li, 1]),
                        tail_half_count=int(metrics[ci, mi, 0, li, 2]), extreme_count=int(metrics[ci, mi, 0, li, 3]),
                        tail_half_pct=100*item['tail_half']['value'], extreme_pct=100*item['extreme']['value'],
                        left_rise=item['left_rise']['value'], left_window=item['window_counts'][0],
                        shoulder_window=item['window_counts'][1]))
    with (root/'paired_cut_effects.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=['comparison', 'layer', 'threshold',
            'first_rate', 'second_rate', 'delta_pp', 'lower95_pp', 'upper95_pp', 'method'])
        writer.writeheader()
        for comparison, layers in summary['comparisons'].items():
            for layer, thresholds in layers.items():
                for threshold, result in thresholds.items():
                    writer.writerow(dict(comparison=comparison, layer=layer, threshold=threshold, **result))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['reduce','merge','plot','tables'],required=True)
    parser.add_argument('--task',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.phase=='reduce': reduce(json.loads(args.task.read_text()),args.output)
    elif args.phase=='merge': merge(args.output);plot(args.output);export_tables(args.output)
    elif args.phase=='tables':export_tables(args.output)
    else:plot(args.output);export_tables(args.output)


if __name__=='__main__':main()
