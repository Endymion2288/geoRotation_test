#!/usr/bin/env python3
"""Plot complete-C transport controls measured with native Calypso."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from charge_component_refit import metrics
from slope_tables import write_json


SAMPLES = ('100044', '100116', '100117')
CONTROLS = ('native', 'material_off', 'step100mm', 'step10mm')
LABELS = dict(native='Default step limit: 10 m', material_off='Material off: 10 m',
              step100mm='Step limit: 100 mm', step10mm='Step limit: 10 mm')
COLORS = dict(native='#326ba4', material_off='#b9473d', step100mm='#438651', step10mm='#9a6c27')


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    rows = json.loads((source/'paired_cases.json').read_text())
    groups = {}
    for key, group in summary['groups'].items():
        dsid, control = key.split('|')
        selected = [r for r in rows if r['dsid'] == dsid and r['control'] == control]
        for mode, expected in group['modes'].items():
            subset = [r['output'] for r in selected if r['mode'] == mode and 'output' in r]
            for field, result in expected.items():
                assert metrics([r[field] for r in subset]) == result
        groups[key] = dict(group, source_surface_types=dict(Counter(str(r['source_surface_type']) for r in selected)))
    return dict(groups=groups, rows=rows,
        validation=json.loads((source/'validation.json').read_text()))


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
        fig.text(.5, .042, 'Actual Calypso extrapolation | Full bound covariance | Native Plane/Disc surfaces | Isolated CKF input states',
            ha='center', fontsize=9)
        fig.text(.5, .022, 'Four selected extreme paths from 60 source events | Local transport controls; not internal CKF-link replay or posterior closure',
            ha='center', fontsize=9)
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
            reference='S06 native full-covariance bound-to-bound transport controls'))

    with PdfPages(output/'bound_transport_report.pdf') as pdf:
        fig, ax = plt.subplots(figsize=(16, 9))
        ax.axis('off')
        cells = []
        for dsid in SAMPLES:
            for control in CONTROLS:
                c = data['groups'][dsid+'|'+control]['counts']
                cells.append([dsid, LABELS[control], str(c['native_inputs_inputs']),
                    str(c['fixed_before_input_inputs']),
                    str(c['native_inputs_before_failed']+c['fixed_before_input_before_failed']),
                    str(c['native_inputs_after_failed']+c['fixed_before_input_after_failed'])])
        table = ax.table(cellText=cells, colLabels=['DSID', 'Control', 'Native-input\npairs', 'Fixed-input\npairs',
            'Before\nfailed', 'After\nfailed'], cellLoc='center', loc='center', colWidths=[.11, .31, .15, .15, .12, .12])
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2.2)
        for (row, _), cell in table.get_celld().items():
            if row == 0:
                cell.set_facecolor('#dce6ef')
        ax.text(.5, .93, 'Both input modes use the same native source and target geometry IDs. All failures remain in the counts.',
            ha='center', transform=ax.transAxes, fontsize=10)
        ax.text(.5, .065, 'Fixed input: identical native local coordinates; analytic rotation of global angles and the complete tangent covariance.\n'
            'Position is checked against native global coordinates using the recorded surface type.\n'
            'The extension reproduces all 18 legacy array fields and material records exactly in both versions.',
            ha='center', transform=ax.transAxes, fontsize=10, linespacing=1.7)
        fig.tight_layout(rect=[0, .085, 1, .94])
        save(fig, 'bound_transport_coverage', 'Native full-C transport: paired coverage and controls', pdf)

        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        for column, dsid in enumerate(SAMPLES):
            for row_index, (field, label) in enumerate([
                ('direction_difference_norm', 'Physical unit-direction difference'),
                ('direction_scaled_max', 'Scaled direction-C difference')]):
                ax = axes[row_index, column]
                for control in CONTROLS:
                    selected = [r for r in data['rows'] if r['dsid'] == dsid and r['control'] == control
                        and r['mode'] == 'fixed_before_input' and 'output' in r]
                    x = np.sort(np.maximum([r['output'][field] for r in selected], 1e-20))
                    ax.step(x, (np.arange(len(x))+1)/len(x), where='post', label=LABELS[control], color=COLORS[control])
                ax.set_xscale('log')
                ax.set_ylim(0, 1.025)
                ax.set_title(dsid+' | Identical physical tangent input')
                ax.set_xlabel(label+' (zero shown at 1e-20)')
                ax.set_ylabel('Fraction of pairs with both transports valid')
                ax.grid(alpha=.2)
                ax.legend(fontsize=8, loc='lower right')
        fig.tight_layout(rect=[0, .09, 1, .94])
        save(fig, 'bound_transport_fixed_inputs', 'Before/after geometry: fixed complete input, material and step controls', pdf)

        fig, axes = plt.subplots(2, 2, figsize=(16, 10))
        main = [r for r in data['rows'] if (r['dsid'], r['entry'], r['component']) == ('100117', 1, 24)]
        for column, mode in enumerate(('native_inputs', 'fixed_before_input')):
            for row_index, (field, label) in enumerate([
                ('direction_difference_norm', 'Physical unit-direction difference'),
                ('direction_scaled_max', 'Scaled direction-C difference')]):
                ax = axes[row_index, column]
                for control in CONTROLS:
                    selected = sorted([r for r in main if r['mode'] == mode and r['control'] == control
                        and 'output' in r], key=lambda r:r['to_forward_index'])
                    ax.plot([r['to_forward_index'] for r in selected],
                        np.maximum([r['output'][field] for r in selected], 1e-20),
                        marker='.', label=LABELS[control], color=COLORS[control])
                ax.axvline(8, color='#555555', alpha=.4, ls=':')
                ax.set_yscale('log')
                ax.set_title('Own native input' if mode == 'native_inputs' else 'Identical physical tangent input')
                ax.set_xlabel('Target native forward-state index (dotted: long link to index 8)')
                ax.set_ylabel(label+' (zero shown at 1e-20)')
                ax.grid(alpha=.2)
                ax.legend(fontsize=8, loc='best')
        fig.tight_layout(rect=[0, .09, 1, .94])
        save(fig, 'bound_transport_extreme', 'DSID 100117, entry 1, component 24: inherited input and integration controls', pdf)
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
    pages = PdfReader(output/'bound_transport_report.pdf').pages
    assert len(index) == len(pages) == 3
    for row, page in zip(index, pages):
        with Image.open(output/row['png']) as png:
            png.load()
        original = PdfReader(output/row['pdf']).pages[0]
        assert graphics(original) == graphics(page)
        assert (list(original.mediabox), list(original.cropbox), original.rotation) == (
            list(page.mediabox), list(page.cropbox), page.rotation)
    write_json(output/'validation.json', dict(pages=3, native_control_metrics_recomputed=True,
        plotted_values_exact=True, PNG_decoded=True, normalized_PDF_content_and_boxes_equal=True,
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
