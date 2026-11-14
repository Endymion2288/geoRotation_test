#!/usr/bin/env python3
"""Plot native CKF chart interventions and common-state numerical evidence."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from charge_component_refit import metrics
from charge_refit_input_figures import COLORS, LABELS, SAMPLES
from slope_tables import write_json


MODES = ('native', 'update', 'smoothing', 'both')
NAMES = ('Native CKF', 'Regular update', 'Regular smoothing', 'Both regular')


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    pairs = json.loads((source/'paired_components.json').read_text())
    states = []
    with gzip.open(source/'paired_state_metrics.jsonl.gz', 'rt') as stream:
        states.extend(json.loads(line) for line in stream)
    for dsid in SAMPLES:
        keys = []
        for mode in MODES:
            rows = [r for r in pairs if r['dsid'] == dsid and r['mode'] == mode]
            counts = summary['groups'][dsid+'|'+mode]['counts']
            assert len(rows) == counts['components']
            for side in ('before', 'after'):
                assert sum(r[side+'_selected'] for r in rows) == counts[side+'_selected']
                for correct in (True, False):
                    assert sum(r[side+'_selected_correct'] is correct for r in rows) == counts.get(side+'_selected'+('_correct' if correct else '_wrong'), 0)
            keys.append({(r['entry'], r['component'], r['hit']) for r in states
                         if r['dsid'] == dsid and r['mode'] == mode and r['same_ordered_CKF_hits']})
        common = set.intersection(*keys)
        frozen = summary['common_cohorts'][dsid]
        assert frozen['n_states'] == len(common)
        assert frozen['components'] == len({(e, c) for e, c, _ in common})
        for mode in MODES:
            selected = [r for r in states if r['dsid'] == dsid and r['mode'] == mode
                        and (r['entry'], r['component'], r['hit']) in common]
            for stage in ('predicted', 'filtered', 'smoothed'):
                for field, frozen_metric in frozen['modes'][mode][stage].items():
                    assert metrics([r[stage][field] for r in selected]) == frozen_metric
    cases = []
    for dsid, entry, component in [('100044', 1, 23), ('100044', 3, 19), ('100044', 12, 2), ('100117', 0, 25)]:
        rows = [r for r in pairs if r['dsid'] == dsid and r['entry'] == entry and r['component'] == component]
        assert len(rows) == 4
        cases.append(dict(dsid=dsid, entry=entry, component=component, modes={r['mode']:r for r in rows}))
    return dict(summary=summary, cases=cases)


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.ticker import MaxNLocator

    data = values(source)
    output.mkdir(parents=True, exist_ok=False)
    index = []

    def save(fig, name, title, caption, pdf):
        fig.suptitle(title, fontsize=14)
        fig.text(.5, .038, 'Actual beforeGeo CKF vs afterGeo native | Isolated smoothing ON | Same n=32, within C/256 seed components', ha='center', fontsize=8)
        fig.text(.5, .018, caption, ha='center', fontsize=8)
        fig.tight_layout(rect=[0, .085, 1, .94])
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
                          reference='S06 actual Calypso CKF chart control with common finite-prior inputs'))

    with PdfPages(output/'ckf_chart_report.pdf') as pdf:
        fig, axes = plt.subplots(2, 3, figsize=(15, 9))
        for j, dsid in enumerate(SAMPLES):
            counts = [data['summary']['groups'][dsid+'|'+m]['counts'] for m in MODES]
            ax = axes[0, j]
            wrong = [c.get('before_selected_wrong', 0) for c in counts]
            selected = [c['before_selected'] for c in counts]
            bars = ax.bar(np.arange(4), wrong, color=COLORS[j])
            ax.bar_label(bars, labels=[f'{w}/{n}' for w, n in zip(wrong, selected)], padding=3, fontsize=8)
            after_wrong = counts[0].get('after_selected_wrong', 0)
            assert all(c.get('after_selected_wrong', 0) == after_wrong for c in counts)
            ax.axhline(after_wrong, color='black', linestyle='--', linewidth=1,
                       label=f'After native: {after_wrong}/{counts[0]["after_selected"]}')
            ax.set_xticks(np.arange(4), ['Native', 'Update', 'Smooth', 'Both'], fontsize=8)
            ax.set_ylim(0, 1.3*max(max(wrong), after_wrong, 1))
            ax.set_ylabel('Wrong-sign CKF components; label = wrong / selected')
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.set_title(LABELS[j])
            ax.legend(fontsize=8)
            ax.grid(axis='y', alpha=.2)
            ax = axes[1, j]
            for offset, key, label, color in [(-.24, 'before_selected', 'Before selected', '#555555'),
                    (0, 'same_selected_hit_sequence', 'Same before/after ordered hits', '#326ba4'),
                    (.24, 'refit_success_both', 'Both refits successful', '#448557')]:
                bars = ax.bar(np.arange(4)+offset, [c.get(key, 0) for c in counts], .23, label=label, color=color)
                ax.bar_label(bars, fontsize=7, padding=2)
            ax.set_xticks(np.arange(4), ['Native', 'Update', 'Smooth', 'Both'], fontsize=8)
            ax.set_ylim(0, 1.4*max(c['before_selected'] for c in counts))
            ax.set_ylabel('Diagnostic components')
            ax.set_title('Selection, ordered paths and actual refit coverage', fontsize=9)
            if j == 0:
                ax.legend(fontsize=7)
        save(fig, 'ckf_chart_outcomes', 'Native CKF update and smoothing chart: charge and path outcomes',
             '20 source events per sample; component counts are correlated diagnostics | Before refit uses regular update+smoothing; diagonal scale 1 in both versions', pdf)

        fig, axes = plt.subplots(2, 3, figsize=(15, 9))
        mode_colors = ['#555555', '#b7544d', '#d49a37', '#326ba4']
        for j, dsid in enumerate(SAMPLES):
            common = data['summary']['common_cohorts'][dsid]
            for i, (field, label) in enumerate([
                    ('direction_covariance_scaled_max', 'Physical direction-C discrepancy'),
                    ('qop_difference_GeVinv', r'Absolute $q/p$ discrepancy [GeV$^{-1}$]')]):
                ax = axes[i, j]
                for offset, mode, color, name in zip([-.15, -.05, .05, .15], MODES, mode_colors, NAMES):
                    measured = [common['modes'][mode][s][field] for s in ('predicted', 'filtered', 'smoothed')]
                    medians = np.array([m['median'] for m in measured], dtype=float)
                    p95 = np.array([m['p95'] for m in measured], dtype=float)
                    # An exact zero has no logarithmic coordinate. Retain it
                    # in plot_values and state the number in the legend.
                    zero = int(np.count_nonzero(medians == 0))
                    ax.plot(np.arange(3)+offset, np.where(medians > 0, medians, np.nan), 'o-',
                            color=color, label=name+(f'; {zero} zero medians' if zero else ''), markersize=4)
                    for x, lo, hi in zip(np.arange(3)+offset, medians, p95):
                        if lo > 0 and hi > 0:
                            ax.vlines(x, lo, hi, colors=color, alpha=.45)
                        elif hi > 0:
                            ax.plot([x], [hi], '^', color=color, markersize=4)
                ax.set_yscale('log')
                ax.set_xticks(np.arange(3), ['Predicted', 'Filtered', 'Smoothed'])
                ax.set_ylabel(label)
                ax.set_title(f'{LABELS[j]} | {common["components"]} components, {common["n_states"]} states', fontsize=10)
                ax.grid(alpha=.2)
                if j == 0:
                    ax.legend(fontsize=7)
        save(fig, 'ckf_chart_common_states', 'Physical state discrepancies on one common cohort across all four CKF controls',
             'Same ordered hits in all controls | Circles: positive medians; lines / triangles: p95; zero medians omitted | Direction-C metric: dimensionless, FASER frame', pdf)

        fig, axes = plt.subplots(2, 2, figsize=(13, 9))
        for ax, case in zip(axes.flat, data['cases']):
            rows = [case['modes'][m] for m in MODES]
            charge = rows[0]['truth_charge']
            y = [1000*charge*r['before_final_qop_GeVinv'] if r['before_final_qop_GeVinv'] is not None else np.nan for r in rows]
            ax.plot(np.arange(4), y, 'o-', color=COLORS[0], label='Before regular refit')
            after = rows[0]['after_final_qop_GeVinv']
            ax.axhline(1000*charge*after, color=COLORS[2], linestyle='--', label='After native')
            ax.axhline(0, color='gray', linewidth=.8)
            for k, r in enumerate(rows):
                if r['before_selected']:
                    if not r['before_refit_success']:
                        ax.scatter([k], [y[k]], marker='x', color='black', s=75, zorder=5)
                    ax.annotate('same hits' if r['same_selected_hit_sequence'] else 'different hits',
                                (k, y[k]), xytext=(0, 10), textcoords='offset points', ha='center', fontsize=7)
            ax.set_xticks(np.arange(4), ['Native', 'Update', 'Smooth', 'Both'])
            ax.set_ylabel(r'$q_{truth}(q/p)_{final}$ [$10^{-3}$ GeV$^{-1}$]')
            ax.set_title(f'{case["dsid"]} / entry {case["entry"]} / component {case["component"]}', fontsize=10)
            ax.margins(y=.2)
            ax.grid(alpha=.2)
        axes[0, 0].legend(fontsize=8)
        save(fig, 'ckf_chart_original_cases', 'Four original refit sign cases after the CKF chart intervention',
             'Positive = correct sign; negative = wrong sign; crosses = before refit fallback | After reference and all before refit settings are fixed', pdf)
    write_json(output/'plot_values.json', data)
    write_json(output/'figure_index.json', index)


def validate(source, output, pdf_tools):
    import re
    import sys
    from PIL import Image
    sys.path.insert(0, str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics
    assert json.loads((output/'plot_values.json').read_text()) == values(source)
    index = json.loads((output/'figure_index.json').read_text())
    pdf = PdfReader(output/'ckf_chart_report.pdf')
    assert len(index) == len(pdf.pages) == 3
    def normalize(page):
        fonts = page['/Resources']['/Font']
        return re.sub(r'/F\d+(?=\s|-)', lambda m: re.sub(r'^/[A-Z]{6}\+', '/',
            str(fonts[m[0]].get_object()['/BaseFont'])), normalized_graphics(page))
    for item, page in zip(index, pdf.pages):
        with Image.open(output/item['png']) as png:
            png.load()
        original = PdfReader(output/item['pdf']).pages[0]
        assert normalize(original) == normalize(page)
        assert list(original.mediabox) == list(page.mediabox)
        assert list(original.cropbox) == list(page.cropbox) and original.rotation == page.rotation
    write_json(output/'validation.json', dict(pages=3, PNG_decoded=True,
        native_counts_and_common_cohort_metrics_recomputed=True,
        plotted_values_exact=True, normalized_PDF_content_and_boxes_equal=True,
        visual_inspection='Pending original-resolution inspection'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', choices=['plot', 'validate'], required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--pdf-tools', type=Path)
    args = parser.parse_args()
    if args.phase == 'plot':
        plot(args.input, args.output)
    else:
        validate(args.input, args.output, args.pdf_tools)


if __name__ == '__main__':
    main()
