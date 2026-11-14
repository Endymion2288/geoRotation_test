#!/usr/bin/env python3
"""Visualize and independently check the corrected-material charge mechanism."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import numpy as np

from slope_tables import write_json


LABELS = ['Legacy encoding', 'Corrected encoding']
COLORS = ['#c43c39', '#205fb2']
MODES = ['legacy', 'corrected']


def values(source):
    events = json.loads((source/'events.json').read_text())
    assert len(events) == 119
    result = {}
    for dsid in ['100116', '100117']:
        wrong = sorted([r for r in events if r['dsid'] == dsid and r['cohort'] == 'wrong'],
                       key=lambda r: r['sample_entry'])
        assert len(wrong) == (27 if dsid == '100116' else 32)
        result[dsid] = {}
        charge = 1 if dsid == '100116' else -1
        for mode in MODES:
            origins = Counter(r[mode]['final_wrong_entry'] for r in wrong)
            crossing = Counter(r[mode]['last_update_into_wrong']['hit_truth'] for r in wrong
                               if r[mode]['final_wrong_entry'] == 3)
            result[dsid][mode] = dict(
                origins=[origins[i] for i in range(1, 7)],
                crossing_categories={k: crossing[k] for k in
                    ['no_SDO', 'primary_only', 'primary_and_other', 'no_primary_support',
                     'not_in_legacy_measurement_audit']},
                sample_entries=[r['sample_entry'] for r in wrong],
                final_correct_signed_qop=[charge*r[mode]['final_qop_GeVinv'] for r in wrong],
                last_update_sample_entries=[r['sample_entry'] for r in wrong if r[mode]['final_wrong_entry'] == 3],
                last_update_predicted=[charge*r[mode]['last_update_into_wrong']['predicted_qop_GeVinv']
                                       for r in wrong if r[mode]['final_wrong_entry'] == 3],
                last_update_filtered=[charge*r[mode]['last_update_into_wrong']['filtered_qop_GeVinv']
                                      for r in wrong if r[mode]['final_wrong_entry'] == 3])
    return result


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    data = values(source)
    output.mkdir(exist_ok=False)
    index = []

    def save(fig, name, title, pdf):
        fig.suptitle(title, fontsize=13)
        fig.text(.5, .025, 'Actual Calypso | 59 original wrong signs + 60 frozen normal controls | Isolated smoothing ON',
                 ha='center', fontsize=8)
        fig.text(.5, .008, 'Only material rotation encoding changes; all 58 output arrays reproduce the corrected full-sample reference',
                 ha='center', fontsize=8)
        fig.tight_layout(rect=[0, .065, 1, .95])
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
            reference='Frozen 119-event mechanism replay; outcomes are not population rates'))

    with PdfPages(output/'material_replay_report.pdf') as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(13, 8))
        for row, dsid in enumerate(['100116', '100117']):
            for column, family in enumerate(['origins', 'crossing_categories']):
                ax = axes[row, column]
                for offset, mode, color, label in zip([-.18, .18], MODES, COLORS, LABELS):
                    v = data[dsid][mode][family]
                    y = list(v.values()) if isinstance(v, dict) else v
                    x = np.arange(len(y))+offset
                    bars = ax.bar(x, y, width=.35, color=color, alpha=.85, label=label)
                    ax.bar_label(bars, labels=[str(n) if n else '' for n in y], fontsize=8, padding=2)
                labels = (['Seed', 'CKF prediction', 'CKF update', 'Native target', 'Conversion', 'Refit']
                          if column == 0 else ['No SDO', 'Primary only', 'Primary + other',
                                               'No primary', 'Audit missing'])
                ax.set_xticks(np.arange(len(labels)), labels, rotation=20, ha='right', fontsize=8)
                ax.set_ylim(0, max(max(data[dsid][m]['origins']) for m in MODES)+5)
                ax.set_ylabel('Wrong tracks')
                ax.set_title(dsid+' | '+('Last entry into final wrong sign' if column == 0 else 'Crossing-hit SDO support'))
                ax.grid(axis='y', alpha=.2)
        axes[0, 0].legend(fontsize=8)
        save(fig, 'material_replay_origins', 'Charge mechanism after correcting material encoding', pdf)

        fig, axes = plt.subplots(2, 2, figsize=(13, 8))
        for row, dsid in enumerate(['100116', '100117']):
            a, b = [data[dsid][m] for m in MODES]
            x, y = [np.array(d['final_correct_signed_qop']) for d in [a, b]]
            ax = axes[row, 0]
            ax.scatter(x, y, s=23, color=COLORS[1])
            low, high = min(x.min(), y.min()), max(x.max(), y.max())
            ax.plot([low, high], [low, high], '--', color='gray', lw=1)
            ax.set_xlabel(r'Legacy $q_{truth}(q/p)_{final}$ [GeV$^{-1}$]')
            ax.set_ylabel(r'Corrected $q_{truth}(q/p)_{final}$ [GeV$^{-1}$]')
            ax.set_title(dsid+' | Same frozen wrong-event cohort')
            ax.grid(alpha=.2)
            ax = axes[row, 1]
            for mode, color, label, marker in zip(MODES, COLORS, LABELS, ['x', 'o']):
                d = data[dsid][mode]
                ax.scatter(d['last_update_predicted'], d['last_update_filtered'],
                           color=color, marker=marker, s=26, alpha=.7, label=label)
            ax.axhline(0, color='gray', lw=1)
            ax.axvline(0, color='gray', lw=1)
            ax.set_xlabel(r'$q_{truth}(q/p)_{predicted}$ [GeV$^{-1}$]')
            ax.set_ylabel(r'$q_{truth}(q/p)_{filtered}$ [GeV$^{-1}$]')
            ax.set_title(dsid+' | Last correct-to-wrong CKF update')
            ax.grid(alpha=.2)
        axes[0, 1].legend(fontsize=8)
        save(fig, 'material_replay_curvature', 'Final curvature and actual sign-crossing updates', pdf)
    write_json(output/'plot_values.json', data)
    write_json(output/'figure_index.json', index)


def validate(source, output, pdf_tools):
    from PIL import Image
    sys.path.insert(0, str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics
    actual = json.loads((output/'plot_values.json').read_text())
    assert actual == values(source)
    index = json.loads((output/'figure_index.json').read_text())
    report = PdfReader(output/'material_replay_report.pdf')
    assert len(index) == len(report.pages) == 2
    for row, merged in zip(index, report.pages):
        with Image.open(output/row['png']) as png:
            png.load()
        single = PdfReader(output/row['pdf']).pages[0]
        assert normalized_graphics(single) == normalized_graphics(merged)
        assert list(single.mediabox) == list(merged.mediabox)
        assert list(single.cropbox) == list(merged.cropbox) and single.rotation == merged.rotation
    write_json(output/'validation.json', dict(pages=2, source_values_exact=True,
        PNG_decoded=True, normalized_vector_graphics_and_boxes_equal=True,
        visual_inspection='Pending separate image review'))
    print('Verified two corrected-material mechanism pages')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', choices=['plot', 'validate'], required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--figures', type=Path, required=True)
    parser.add_argument('--pdf-tools', type=Path)
    args = parser.parse_args()
    if args.phase == 'plot':
        plot(args.input, args.figures)
    else:
        validate(args.input, args.figures, args.pdf_tools)


if __name__ == '__main__':
    main()
