#!/usr/bin/env python3
"""Plot actual Calypso material variances and their coordinate arithmetic."""
import argparse
import json
from pathlib import Path

import numpy as np

from charge_component_refit import audit, metrics
from slope_tables import write_json


SAMPLES = ('100044', '100116', '100117')
CONTROLS = ('native', 'step100mm', 'step10mm')
MODES = ('fixed_before_input', 'native_inputs')
LABELS = dict(native='10 m (reference)', step100mm='100 mm', step10mm='10 mm')


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    rows = list(audit(source/'paired_interactions.jsonl.gz'))
    for key, group in summary['groups'].items():
        dsid, control, mode = key.split('|')
        selected = [r for r in rows if (r['dsid'], r['control'], r['mode']) == (dsid, control, mode)]
        assert len(selected) == group['counts'].get('paired_interactions', 0)
        for field, expected in group['differences'].items():
            assert metrics([r['difference'][field] for r in selected]) == expected
        for variant, fields in group['native_relations'].items():
            for field, expected in fields.items():
                assert metrics([r[variant][field] for r in selected if r[variant]['scattering']]) == expected
    return dict(summary=summary, rows=rows,
        coverage=json.loads((source/'validation.json').read_text()))


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
        fig.text(.5, .043, 'Actual Calypso propagationSteps material records | Identical-input and inherited-input controls | FASER direction tensor',
            ha='center', fontsize=9)
        fig.text(.5, .022, 'Same-double-square substitution is an arithmetic diagnostic of recorded terms; fitted track covariance is not rerun.',
            ha='center', fontsize=9)
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf', reference='S06 native material-noise audit'))

    with PdfPages(output/'material_noise_report.pdf') as pdf:
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))
        for ax, mode in zip(axes, MODES):
            ax.axis('off')
            cells = []
            for dsid in SAMPLES:
                for control in CONTROLS:
                    c = data['summary']['groups'][dsid+'|'+control+'|'+mode]['counts']
                    cells.append([dsid, LABELS[control], c['input_links'],
                        c['same_trace_sequence'], c['different_trace_sequence'], c['paired_interactions']])
            table = ax.table(cellText=cells,
                colLabels=['DSID', 'Trace step limit', 'Input\nlinks', 'Same\nsequence', 'Different\nsequence', 'Paired\ninteractions'],
                cellLoc='center', bbox=[0., .23, 1., .59], colWidths=[.12, .29, .14, .15, .15, .15])
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            for (r, _), cell in table.get_celld().items():
                if r == 0:
                    cell.set_facecolor('#dce6ef')
            ax.text(.5, .89, 'Identical physical input' if mode == 'fixed_before_input' else 'Each version\'s frozen native input',
                ha='center', transform=ax.transAxes, fontsize=12)
        coverage = data['coverage']
        fig.text(.5, .15, f"24 transport runs reproduce all original arrays exactly. {coverage['complete_trace_sequence_counts']} paired input links retained.\n"
            f"Material-OFF did not request trace output: {coverage['unrecorded_links']} links are unmeasured for noise; transport replay remains verified.\n"
            f"{coverage['paired_interactions']} paired material interactions; links with no interaction remain in the sequence denominator.",
            ha='center', fontsize=10, linespacing=1.7)
        fig.subplots_adjust(left=.025, right=.985, top=.91, bottom=.17, wspace=.08)
        save(fig, 'material_noise_coverage', 'Native material-noise records: complete transport and trace coverage', pdf)

        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        for i, mode in enumerate(MODES):
            for j, dsid in enumerate(SAMPLES):
                ax = axes[i, j]
                for control, color in [('native', '#326ba4'), ('step10mm', '#9a6c27')]:
                    selected = [r for r in data['rows'] if (r['dsid'], r['mode'], r['control']) == (dsid, mode, control)]
                    for field, style, name in [('physical_direction_Q_scaled_max', '-', 'Recorded Q'),
                            ('same_double_square_Q_scaled_max', '--', 'Same-double-square diagnostic')]:
                        x = np.sort(np.maximum([r['difference'][field] for r in selected], 1e-20))
                        ax.step(x, (np.arange(len(x))+1)/len(x), where='post', color=color, linestyle=style,
                            label=name+' | '+('10 m' if control == 'native' else '10 mm'))
                ax.set_xscale('log')
                ax.set_ylim(0, 1.05)
                ax.set_title(dsid+' | '+('Identical input' if i == 0 else 'Native inputs'))
                ax.set_xlabel('Scaled physical direction-Q difference (zero at 1e-20)')
                ax.set_ylabel('Fraction of paired interactions')
                ax.grid(alpha=.2)
                ax.legend(fontsize=7, loc='lower right')
        fig.tight_layout(rect=[0, .09, 1, .94])
        save(fig, 'material_noise_coordinates', 'Recorded direction noise: rotation difference and arithmetic diagnostic', pdf)

        fig, axes = plt.subplots(1, 2, figsize=(16, 8))
        selected = [r for r in data['rows'] if r['control'] == 'native' and r['mode'] == 'fixed_before_input']
        for variant, color in [('before', '#326ba4'), ('after', '#b54a42')]:
            native = [r[variant] for r in selected if r[variant]['scattering']]
            axes[0].scatter([r['expected_isotropy_from_float_square']*1e8 for r in native],
                [r['isotropy_error']*1e8 for r in native], label=variant+'Geo', color=color, s=24, alpha=.65)
            for field, style, label in [('float_square_relation_error', '-', 'Theta variance vs float square'),
                    ('double_square_relation_error', '--', 'Phi variance vs double square')]:
                axes[1].plot(np.maximum([abs(r[field]) for r in native], 1e-20), linestyle=style,
                    marker='.', color=color, label=variant+'Geo | '+label)
        axes[0].plot([-6, 6], [-6, 6], color='black', linewidth=1, label='Equal')
        axes[0].set_xlabel('Float-square predicted anisotropy / 1e-8')
        axes[0].set_ylabel('Recorded phi/theta anisotropy / 1e-8')
        axes[0].set_title('Isotropy error from recorded native variances')
        axes[1].set_yscale('log')
        axes[1].set_xlabel('Paired interaction index; identical inputs, native limit')
        axes[1].set_ylabel('Absolute relative arithmetic relation error (zero at 1e-20)')
        axes[1].set_title('Float theta-square and double phi-square relations')
        for ax in axes:
            ax.grid(alpha=.2)
            ax.legend(fontsize=8, loc='best')
        fig.tight_layout(rect=[0, .1, 1, .93])
        save(fig, 'material_noise_arithmetic', 'Native scattering variances: measured precision relations', pdf)
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
    pages = PdfReader(output/'material_noise_report.pdf').pages
    index = json.loads((output/'figure_index.json').read_text())
    assert len(pages) == len(index) == 3
    for page, row in zip(pages, index):
        with Image.open(output/row['png']) as png:
            png.load()
        single = PdfReader(output/row['pdf']).pages[0]
        assert graphics(page) == graphics(single)
        assert (list(page.mediabox), list(page.cropbox), page.rotation) == (
            list(single.mediabox), list(single.cropbox), single.rotation)
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
