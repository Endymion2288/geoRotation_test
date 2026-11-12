#!/usr/bin/env python3
"""Verify audit-to-writer state ownership and plot GitHub-cut OFF/ON counts."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

LAYERS = [(s, l) for s in (1, 2, 3) for l in (0, 1, 2)]


def audit_events(path):
    rows = []
    with path.open() as stream:
        for line in stream:
            row = json.loads(line)
            if row['stage'] == 'event' and rows:
                yield rows
                rows = []
            rows.append(row)
    if rows:
        yield rows


def check_writer(scratch, destination):
    """Compare actual writer parameters to old audit snapshots, before scratch cleanup."""
    import ROOT
    source = ROOT.TFile.Open(str(scratch/'track_states.root'))
    tree = source.Get('tree')
    names = ['eLOC0', 'eLOC1', 'ePHI', 'eTHETA', 'eQOP', 'eT']
    branches = ['event_nr', 'track_nr', 'volume_id', 'layer_id', 'module_id',
                't_eQOP', 'NDF', 'chi2Sum', 'rMax_smt']
    branches += [name+'_'+stage for name in names for stage in ('prt', 'flt', 'smt')]
    tree.SetBranchStatus('*', 0)
    for branch in branches:
        tree.SetBranchStatus(branch, 1)
    counts = Counter()
    tails = {f'{s}_{l}': {mode: Counter() for mode in ('shared', 'isolated')} for s, l in LAYERS}
    histograms = {mode: np.zeros((9, 102), dtype=np.int64) for mode in ('shared', 'isolated')}
    ordinal = 0
    max_rmax_error = 0.
    max_relative_writer_error = 0.
    for event_index, rows in enumerate(audit_events(scratch/'stages.jsonl')):
        counts['source_events'] += 1
        selected = [r for r in rows if r['stage'] == 'selected']
        states = {}
        for row in rows:
            if row['stage'] in ('selected_shared_state', 'selected_state'):
                key = row['stage'], row['seed'], row['candidate']
                states.setdefault(key, []).append(row)
        for track_index, metadata in enumerate(selected):
            tree.GetEntry(ordinal)
            assert int(tree.track_nr) == track_index
            assert int(tree.event_nr) == metadata['event']
            ordinal += 1
            counts['writer_tracks'] += 1
            key = metadata['seed'], metadata['candidate']
            shared_all = states[('selected_shared_state', *key)]
            isolated_all = states[('selected_state', *key)]
            assert len(shared_all) == len(isolated_all)
            shared = [r for r in shared_all if r['measurement_flag']]
            isolated = [r for r in isolated_all if r['measurement_flag']]
            assert len(isolated) == len(tree.t_eQOP)
            passing = tree.NDF > 9 and tree.chi2Sum/tree.NDF < 20 and tree.rMax_smt < 90
            counts['github_selected_tracks'] += passing
            rmax = 0.
            for i, (old, new) in enumerate(zip(shared, isolated)):
                assert (old['hit'], old['geometry_id']) == (new['hit'], new['geometry_id'])
                gid = new['geometry_id']
                assert gid >> 56 == int(tree.volume_id[i])
                assert (gid >> 36) & 0xfff == int(tree.layer_id[i])
                assert (gid >> 8) & 0xfffffff == int(tree.module_id[i])
                for stage, field in [('prt', 'predicted'), ('flt', 'filtered'), ('smt', 'smoothed')]:
                    if stage != 'smt':
                        assert old[field] == new[field]
                    for j, name in enumerate(names):
                        value = float(getattr(tree, name+'_'+stage)[i])
                        expected = new[field][j]
                        error = abs(value-expected)/max(abs(expected), 1e-25)
                        max_relative_writer_error = max(max_relative_writer_error, error)
                        # ROOT writer persists bound parameters as float32.
                        assert value == float(np.float32(expected)), (event_index, track_index, i, stage, name, value, expected)
                        counts['writer_parameters_equal_after_float32_cast'] += 1
                counts['measurement_states'] += 1
                counts['smoothed_qop_changed'] += old['smoothed'][4] != new['smoothed'][4]
                transform = np.array(new['surface'])
                global_position = transform @ np.array([*new['smoothed'][:2], 0., 1.])
                rmax = max(rmax, float(np.hypot(global_position[1], global_position[2])))
                if not passing:
                    continue
                station = {25: 1, 27: 2, 29: 3}[int(tree.volume_id[i])]
                layer = int(tree.layer_id[i])//2-1
                layer_index = LAYERS.index((station, layer))
                truth = float(tree.t_eQOP[i])
                for mode, state in [('shared', old), ('isolated', new)]:
                    qop = state['smoothed'][4]
                    if not truth or not qop or not np.isfinite(truth/qop):
                        continue
                    residual = abs(1000*truth/qop)-1
                    counter = tails[f'{station}_{layer}'][mode]
                    counter['finite'] += 1
                    counter['tail_half'] += residual < -.5
                    counter['tail_point8'] += residual < -.8
                    # Same 100 bins and explicit flow counts as the source macro.
                    bin_index = int(np.searchsorted(np.linspace(-1.2, 1.2, 101), residual, side='right'))
                    histograms[mode][layer_index, bin_index] += 1
            max_rmax_error = max(max_rmax_error, abs(rmax-float(tree.rMax_smt)))
            assert abs(rmax-float(tree.rMax_smt)) < max(2e-5, 1e-7*rmax)
    assert ordinal == tree.GetEntries()
    source.Close()
    result = dict(counts=dict(counts), per_layer={k: {m: dict(v) for m, v in values.items()}
        for k, values in tails.items()}, maximum_relative_writer_rounding_error=max_relative_writer_error,
        maximum_rmax_rounding_error_mm=max_rmax_error,
        all_writer_parameters_equal_post_isolation_float32=True,
        predicted_filtered_hits_geometry_identical=True,
        selection='Freeze writer ON GitHub passing candidates; use identical states/truth in shared and isolated snapshots',
        bins={k: v.tolist() for k, v in histograms.items()})
    (destination/'writer_smoothing_check.json').write_text(json.dumps(result, indent=2)+'\n')


def check_smoke_histograms(destination, reference):
    from reference_figures import read_histograms
    comparisons = {}
    for truth in ('native', 'production', 'common_native', 'common_production'):
        name = truth+'_3station_counts.root'
        first = read_histograms(destination/name, 'audit_'+truth)
        second = read_histograms(reference/name, 'minimal_'+truth)
        for key, histogram in first.items():
            other = second[key]
            assert histogram.GetEntries() == other.GetEntries(), key
            assert all(histogram.GetBinContent(i) == other.GetBinContent(i) for i in range(102)), key
        comparisons[truth] = dict(histograms=27, all_bins_flow_and_raw_entries_equal=True)
    (destination/'audit_vs_minimal.json').write_text(json.dumps(dict(
        reference=str(reference), comparisons=comparisons), indent=2)+'\n')


def draw_smoothed(histograms, prefix, title, labels, log=False, denominator='raw Entries'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 3, figsize=(15, 10), sharex=True)
    bins = {}
    for layer_index, (axis, (station, layer)) in enumerate(zip(axes.flat, LAYERS)):
        for group_index, (group, label) in enumerate(zip(histograms, labels)):
            if isinstance(group, dict):
                histogram = group[(0, station, layer)]
                edges = [histogram.GetBinLowEdge(i) for i in range(1, 102)]
                values = [histogram.GetBinContent(i) for i in range(102)]
                n = histogram.GetEntries()
                sumw2 = [histogram.GetBinError(i)**2 for i in range(102)]
            else:
                edges = np.linspace(-1.2, 1.2, 101)
                values = group[layer_index]
                n = sum(values)
                sumw2 = values
            values = np.asarray(values)
            axis.stairs(values[1:-1]/n, edges, label=f'{label}: N={int(n):,}',
                color=['#555555', '#1f77b4'][group_index], linewidth=1.6)
            bins[f'{group_index}_S{station}L{layer}'] = dict(edges=list(edges),
                counts=values.tolist(), sumw2=list(map(float, sumw2)),
                **{'raw_entries' if denominator == 'raw Entries' else 'finite_entries': float(n)})
        axis.set_title(f'Station {station}, layer {layer}')
        axis.set_xlim(-1.2, 1.2)
        axis.set_xlabel('(reco p - truth p) / truth p')
        axis.set_ylabel('Entries / '+denominator)
        if log:
            axis.set_yscale('log')
            axis.set_ylim(1e-5, 1)
        else:
            axis.set_ylim(bottom=0)
        axis.axvline(-.5, color='.65', linestyle=':', linewidth=.7)
        axis.legend(fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .91))
    fig.savefig(str(prefix)+'.png', dpi=150)
    fig.savefig(str(prefix)+'.pdf')
    plt.close(fig)
    return bins


def full_comparison(root, off, on):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    from reference_figures import read_histograms, overlay
    figures = root/'figures'
    figures.mkdir(exist_ok=True)
    summary = dict(cases={}, coverage={}, cuts='chi2Sum/NDF < 20 && rMax_smt < 90 && NDF > 9',
        truth='SiHit-linked HepMC particle momentum; native writer truth; same formula as GitHub')
    all_bins = {}
    index = []
    for case in ('beforeGeo', 'afterGeo'):
        task_lists = [json.loads((run/(case+'_tasks.json')).read_text()) for run in (off, on)]
        identities = [[(t['input'], t['skip'], t['events']) for t in tasks] for tasks in task_lists]
        assert identities[0] == identities[1]
        assert len(identities[0]) == 20 and all(skip == 0 and n == 1000 for _, skip, n in identities[0])
        assert len({path for path, _, _ in identities[0]}) == 20
        summary['coverage'][case] = identities[0]
        data = []
        counts = []
        for run_index, (run, tasks) in enumerate(zip((off, on), task_lists)):
            counts.append(json.loads((run/'reductions'/(case+'_truth_summary.json')).read_text()))
            for task in tasks:
                manifest = json.loads((Path(task['output'])/'run.json').read_text())
                assert manifest.get('isolate_smoothing', False) == bool(run_index)
                assert manifest['audit'] is False and manifest['reset_event_maps'] is False
            data.append(read_histograms(run/'reductions'/(case+'_native_3station.root'), f'{case}_{run_index}'))
        assert counts[0]['counts']['writer_tracks'] == counts[1]['counts']['writer_tracks']
        rows = {}
        for station, layer in LAYERS:
            pair = []
            for statistics in counts:
                counter = statistics['per_layer'][f'{station}_{layer}']
                n = counter['smt_native_finite']
                tail = counter['smt_native_below_minus_half']
                pair.append(dict(finite=n, tail_half=tail, tail_half_pct=100*tail/n))
            rows[f'S{station}L{layer}'] = dict(off=pair[0], on=pair[1])
        summary['cases'][case] = dict(writer_tracks=counts[0]['counts']['writer_tracks'],
            selected_tracks_off=counts[0]['counts']['github_selected_tracks'],
            selected_tracks_on=counts[1]['counts']['github_selected_tracks'], layers=rows)
        title = '#splitline{110051 | '+case+' | 20 files | 20,000 source entries | shared OFF / isolated ON}{GitHub cuts: chi2/NDF<20, NDF>9, rMax_smt<90 mm | sensor states | '+('native FASER' if case == 'beforeGeo' else 'rotated ACTS')+' frame | SiHit-linked HepMC truth}'
        name = case+'_github_smoothing_off_on'
        canvas, retained = overlay(*data, figures/name, title, ['OFF', 'ON'])
        canvas.Close()
        index.append(name)
        for log in (False, True):
            name = case+'_smoothed_off_on_'+('log' if log else 'linear')
            title = f'110051 | {case} | 20 files | 20,000 source entries; each file entries 0-999\n'
            title += 'GitHub cuts: chi2/NDF < 20, NDF > 9, rMax_smt < 90 mm; cut membership may change with rMax_smt\n'
            title += 'Selected CKF sensor states; '+('native FASER' if case == 'beforeGeo' else 'rotated ACTS')+' frame; SiHit-linked HepMC truth; normalization includes flow/invalid Entries'
            all_bins[name] = draw_smoothed(data, figures/name, title, ['Smoothing OFF', 'Isolated ON'], log)
            index.append(name)
    summary['selection_note'] = 'Identical GitHub cut expression; ON/OFF rMax_smt changes candidate membership. Full plots are not a frozen-candidate comparison.'
    (root/'full_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (root/'bin_counts.json').write_text(json.dumps(all_bins, indent=2)+'\n')
    (root/'figure_index.json').write_text(json.dumps(index, indent=2)+'\n')
    print(json.dumps(summary['cases']), flush=True)


def finish(root, tools):
    import sys
    sys.path.insert(0, str(tools))
    from pypdf import PdfReader, PdfWriter
    index = json.loads((root/'figure_index.json').read_text())
    check = json.loads((root/'audit_check/writer_smoothing_check.json').read_text())
    paired = [np.array(check['bins'][mode]) for mode in ('shared', 'isolated')]
    name = 'afterGeo_same_candidates_audit_check'
    bins = draw_smoothed(paired, root/'figures'/name,
        '110051 | afterGeo | diagnostic subset: first file, entries 0-19 (20 source events)\n'
        'Freeze candidates passing ON writer GitHub cuts; identical hits, P/F parameters, truth and binning\n'
        'Finite smoothed residuals; rotated ACTS frame; SiHit-linked HepMC truth; shared vs old isolated implementation',
        ['Shared snapshot', 'Isolated snapshot'], denominator='finite residuals')
    (root/'paired_bin_counts.json').write_text(json.dumps(bins, indent=2)+'\n')
    index = [item for item in index if item != name] + [name]
    writer = PdfWriter()
    for item in index:
        writer.add_page(PdfReader(root/'figures'/(item+'.pdf')).pages[0])
    output = root/'figures/smoothing_recheck.pdf'
    with output.open('wb') as stream:
        writer.write(stream)
    combined = PdfReader(output)
    assert len(combined.pages) == len(index) == 7
    for page, item in zip(combined.pages, index):
        original = PdfReader(root/'figures'/(item+'.pdf')).pages[0]
        assert page.get_contents().get_data() == original.get_contents().get_data()
        for box in ('mediabox', 'cropbox', 'trimbox', 'bleedbox', 'artbox'):
            assert list(getattr(page, box)) == list(getattr(original, box))
        assert page.rotation == original.rotation
    from PIL import Image
    for item in index:
        with Image.open(root/'figures'/(item+'.png')) as image:
            image.verify()
    (root/'figure_index.json').write_text(json.dumps(index, indent=2)+'\n')
    (root/'validation.json').write_text(json.dumps(dict(pdf_pages=len(combined.pages),
        individual_and_combined_pdf_contents_equal=True, writer_check=check['counts'],
        audit_minimal=json.loads((root/'audit_check/audit_vs_minimal.json').read_text())), indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--off', type=Path)
    parser.add_argument('--on', type=Path)
    parser.add_argument('--finish', type=Path, help='Directory containing pypdf')
    args = parser.parse_args()
    if args.finish:
        finish(args.root, args.finish)
    else:
        full_comparison(args.root, args.off, args.on)


if __name__ == '__main__':
    main()
