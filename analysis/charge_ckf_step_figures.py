#!/usr/bin/env python3
"""Plot native CKF step-limit comparisons with complete selection counts."""
import argparse
import json
from pathlib import Path

import numpy as np

from charge_component_refit import audit, metrics
from charge_ckf_step_control import MODES, SAMPLES
from slope_tables import write_json


LABELS = dict(default='Native default', step100mm='CKF limit: 100 mm', step10mm='CKF limit: 10 mm')
COLORS = dict(default='#326ba4', step100mm='#438651', step10mm='#9a6c27')


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    rows = [r for r in audit(source/'paired_states.jsonl.gz') if r['comparison'] == 'between_geometry']
    groups = {}
    for dsid in SAMPLES:
        for mode in MODES:
            key = dsid+'|between_geometry|'+mode
            group = summary['groups'][key]
            selected = [r for r in rows if r['dsid'] == dsid and r['mode'] == mode]
            assert len(selected) == group['counts']['paired_states']
            for stage, fields in group['stages'].items():
                for field, expected in fields.items():
                    assert metrics([r[stage][field] for r in selected]) == expected
            groups[key] = group
    return dict(groups=groups, rows=rows)


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    data = values(source)
    output.mkdir(parents=True, exist_ok=False)
    index = []

    def save(fig, name, title, pdf):
        fig.suptitle(title, fontsize=14, y=.985)
        fig.text(.5, .042, 'Actual Calypso CKF / conversion / refit | Isolated smoothing ON | n32, within-C / 256 | Analytic after seed-C',
            ha='center', fontsize=9)
        fig.text(.5, .022, '60 distinct source events | Regular before CKF/refit angular charts | Finite components are correlated; not production rates',
            ha='center', fontsize=9)
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf', reference='S06 internal native CKF step-limit control'))

    with PdfPages(output/'ckf_step_report.pdf') as pdf:
        fig, ax = plt.subplots(figsize=(16, 8))
        ax.axis('off')
        cells = []
        for dsid in SAMPLES:
            for mode in MODES:
                c = data['groups'][dsid+'|between_geometry|'+mode]['counts']
                cells.append([dsid, LABELS[mode], str(c['inputs']),
                    f"{c['before_selected']} / {c['after_selected']}",
                    str(c.get('different_ordered_hits', 0)), str(c.get('different_full_state_signatures', 0)),
                    str(c.get('selected_sign_difference', 0)), str(c.get('final_sign_difference', 0))])
        table = ax.table(cellText=cells, colLabels=['DSID', 'CKF step setting', 'Input\ncomponents',
            'Selected\nbefore / after', 'Different\nordered hits', 'Same hits, different\nfull state chain',
            'Selected\nsign differences', 'Final\nsign differences'], cellLoc='center',
            bbox=[0., .27, 1., .56])
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        for (row, _), cell in table.get_celld().items():
            if row == 0:
                cell.set_facecolor('#dce6ef')
        ax.text(.5, .89, 'Selection migrations, topology changes, refit failures and terminal-only paths remain in the analysis.',
            ha='center', transform=ax.transAxes, fontsize=10)
        ax.text(.5, .12, 'Default before/after outputs must reproduce the archived full audit streams, smoother operations and observer arrays exactly.\n'
            'CDFs below use common ordered-hit, full-state-signature paths; the coverage table gives the full component denominator.\n'
            'A smaller between-version difference does not by itself establish charge correctness or a converged mixture posterior.',
            ha='center', transform=ax.transAxes, fontsize=10, linespacing=1.7)
        fig.tight_layout(rect=[0, .085, 1, .94])
        save(fig, 'ckf_step_coverage', 'Internal CKF step limits: complete component and selection coverage', pdf)

        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        for column, dsid in enumerate(SAMPLES):
            for row_index, (field, label) in enumerate([
                ('direction_difference_norm', 'Smoothed unit-direction difference'),
                ('direction_scaled_max', 'Smoothed scaled direction-C difference')]):
                ax = axes[row_index, column]
                for mode in MODES:
                    selected = [r for r in data['rows'] if r['dsid'] == dsid and r['mode'] == mode]
                    x = np.sort(np.maximum([r['smoothed'][field] for r in selected], 1e-20))
                    ax.step(x, (np.arange(len(x))+1)/len(x), where='post', color=COLORS[mode], label=LABELS[mode])
                ax.set_xscale('log')
                ax.set_ylim(0, 1.025)
                ax.set_title(dsid+' | Before/after geometry')
                ax.set_xlabel(label+' (zero shown at 1e-20)')
                ax.set_ylabel('Fraction of common full-chain states')
                ax.grid(alpha=.2)
                ax.legend(fontsize=8, loc='lower right')
        fig.tight_layout(rect=[0, .09, 1, .94])
        save(fig, 'ckf_step_distributions', 'Before/after geometry: full native CKF with explicit step limits', pdf)

        fig, axes = plt.subplots(2, 2, figsize=(16, 10))
        for column, stage in enumerate(('predicted', 'smoothed')):
            for row_index, (field, label) in enumerate([
                ('direction_difference_norm', 'Physical unit-direction difference'),
                ('direction_scaled_max', 'Scaled direction-C difference')]):
                ax = axes[row_index, column]
                for mode in MODES:
                    selected = sorted([r for r in data['rows'] if r['dsid'] == '100117' and r['entry'] == 1
                        and r['component'] == 24 and r['mode'] == mode], key=lambda r:r['forward_index'])
                    ax.plot([r['forward_index'] for r in selected],
                        np.maximum([r[stage][field] for r in selected], 1e-20),
                        marker='.', color=COLORS[mode], label=LABELS[mode])
                ax.set_title(stage.capitalize())
                ax.set_yscale('log')
                ax.set_xlabel('Native forward linked-state index')
                ax.set_ylabel(label+' (zero shown at 1e-20)')
                ax.grid(alpha=.2)
                ax.legend(fontsize=8, loc='best')
        fig.tight_layout(rect=[0, .09, 1, .94])
        save(fig, 'ckf_step_extreme', 'DSID 100117, entry 1, component 24: complete CKF step-limit response', pdf)
    write_json(output/'figure_index.json', index)
    write_json(output/'plot_values.json', data)


def validate(source, output, pdf_tools):
    import re
    import sys
    from PIL import Image
    sys.path.insert(0, str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics
    assert json.loads((output/'plot_values.json').read_text()) == values(source)

    def graphics(page):
        fonts = page['/Resources']['/Font']
        return re.sub(r'/F\d+(?=\s|-)', lambda match: re.sub(r'^/[A-Z]{6}\+', '/',
            str(fonts[match[0]].get_object()['/BaseFont'])), normalized_graphics(page))
    index = json.loads((output/'figure_index.json').read_text())
    pages = PdfReader(output/'ckf_step_report.pdf').pages
    assert len(index) == len(pages) == 3
    for row, page in zip(index, pages):
        with Image.open(output/row['png']) as png:
            png.load()
        original = PdfReader(output/row['pdf']).pages[0]
        assert graphics(original) == graphics(page)
        assert (list(original.mediabox), list(original.cropbox), original.rotation) == (
            list(page.mediabox), list(page.cropbox), page.rotation)
    write_json(output/'validation.json', dict(pages=3, plotted_values_exact=True,
        PNG_decoded=True, normalized_PDF_content_and_boxes_equal=True,
        visual_inspection='Pending original-resolution inspection'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--phase', choices=['plot', 'validate'], required=True)
    parser.add_argument('--pdf-tools', type=Path)
    args = parser.parse_args()
    if args.phase == 'plot':
        plot(args.input, args.output)
    else:
        validate(args.input, args.output, args.pdf_tools)


if __name__ == '__main__':
    main()
