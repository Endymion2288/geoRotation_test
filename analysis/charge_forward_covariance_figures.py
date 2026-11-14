#!/usr/bin/env python3
"""Plot full-state native forward differences and covariance accounting."""
import argparse
import json
from pathlib import Path

import numpy as np

from charge_component_refit import audit, metrics
from charge_forward_covariance import SAMPLES, STAGES
from slope_tables import write_json


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    validation = json.loads((source/'validation.json').read_text())
    rows = list(audit(source/'paired_forward_states.jsonl.gz'))
    groups = {}
    for dsid in SAMPLES:
        selected = [r for r in rows if r['dsid'] == dsid]
        group = summary['groups'][dsid]
        assert len(selected) == group['counts']['paired_states']
        for stage in STAGES:
            for field, expected in group['cohorts']['all']['stages'][stage].items():
                assert metrics([r[stage][field] for r in selected]) == expected
        groups[dsid] = dict(counts=group['counts'],
            stages={stage: [r[stage]['direction_scaled_max'] for r in selected] for stage in STAGES})
    path = [r for r in rows if (r['dsid'], r['entry'], r['component']) == ('100117', 1, 24)]
    assert path and [r['forward_index'] for r in path] == list(range(len(path)))
    return dict(groups=groups, path=path, validation=validation)


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    data = values(source)
    output.mkdir(parents=True, exist_ok=False)
    index = []
    colors = dict(predicted='#438651', filtered='#b9473d', smoothed='#326ba4')

    def save(fig, name, title, pdf):
        fig.suptitle(title, fontsize=14, y=.985)
        fig.text(.5, .044, 'Actual Calypso | Isolated smoothing ON | n32, within-C / 256 | Analytic seed coordinate Jacobian', ha='center', fontsize=9)
        fig.text(.5, .022, '60 distinct source events | Correlated components, first-ranked common ordered-hit paths | Not production rates', ha='center', fontsize=9)
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
            reference='S06 full native forward-state audit; existing production reference plots retained'))

    with PdfPages(output/'forward_covariance_report.pdf') as pdf:
        fig, ax = plt.subplots(figsize=(15, 7))
        ax.axis('off')
        cells = []
        for dsid in SAMPLES:
            c = data['groups'][dsid]['counts']
            cells.append([dsid, str(c['input_components']), str(c['same_ordered_hits']),
                str(c['terminal_only']), str(c['measurement_states']), str(c['passive_states']),
                str(c['different_raw_material_slabs'])])
        table = ax.table(cellText=cells, colLabels=['DSID', 'Input\ncomponents', 'Common selected\nsame-hit paths',
            'Terminal-only\npaths', 'Measurement\nstates', 'Passive\nstates', 'Different raw\nmaterial slabs'],
            cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 3)
        for (row, _), cell in table.get_celld().items():
            if row == 0:
                cell.set_facecolor('#dce6ef')
        ax.text(.5, .82, 'Complete saved P/F/S states and Jacobians; native linked-state order includes passive surfaces', ha='center', transform=ax.transAxes, fontsize=11)
        ax.text(.5, .16, 'Paired surfaces are checked under the full 3D geometry rotation.\n'
            'Covariances use common (local0, local1, phi_A, latitude_A, q/p, time) coordinates.\n'
            'C_pred - J C_prev J^T is a covariance balance, not an observed native process-noise matrix.',
            ha='center', transform=ax.transAxes, fontsize=10, linespacing=1.7)
        fig.tight_layout(rect=[0, .09, 1, .94])
        save(fig, 'forward_covariance_counts', 'Native forward-state coverage and paired material queries', pdf)

        fig, axes = plt.subplots(1, 3, figsize=(16, 6))
        for ax, dsid in zip(axes, SAMPLES):
            group = data['groups'][dsid]
            for stage in STAGES:
                all_values = np.asarray(group['stages'][stage])
                positive = np.sort(all_values[all_values > 0])
                ax.step(positive, (np.arange(len(positive))+1+(all_values == 0).sum())/len(all_values),
                    where='post', label=stage.capitalize(), color=colors[stage])
            ax.set_xscale('log')
            ax.set_ylim(0, 1.02)
            ax.set_title(f'{dsid} | N = {group["counts"]["paired_states"]:,} states')
            ax.set_xlabel('Scaled direction-C difference')
            ax.set_ylabel('Fraction of all paired native states')
            ax.grid(alpha=.2)
            ax.legend(loc='lower right')
        fig.tight_layout(rect=[0, .10, 1, .93])
        save(fig, 'forward_covariance_distributions', 'Before/after geometry: complete native-state covariance differences', pdf)

        path = data['path']
        x = [r['forward_index'] for r in path]
        fig, axes = plt.subplots(2, 2, figsize=(16, 10))
        for stage in STAGES:
            for ax, field in [(axes[0, 0], 'direction_scaled_max'),
                              (axes[0, 1], 'direction_absolute_max_rad2')]:
                y = np.asarray([r[stage][field] for r in path])
                ax.plot(x, np.where(y > 0, y, np.nan), marker='.', color=colors[stage], label=stage.capitalize())
        axes[0, 0].set_ylabel('Scaled direction-C difference')
        axes[0, 1].set_ylabel(r'Absolute direction-C difference [rad$^2$]')
        for term, label, color in [('inherited_C', 'Previous-C contraction', '#648a45'),
                                  ('saved_J', 'Saved-J contraction', '#b78734'),
                                  ('remaining_balance', 'Remaining covariance balance', '#8165a7')]:
            y = np.asarray([r['link'][term]['direction_scaled_max'] for r in path])
            axes[1, 0].plot(x, np.where(y > 0, y, np.nan), marker='.', color=color, label=label)
        axes[1, 0].set_ylabel('Term size / same target predicted-C scale')
        for stage in ('predicted', 'filtered'):
            y = np.asarray([r[stage]['direction_difference_norm'] for r in path])
            axes[1, 1].plot(x, np.where(y > 0, y, np.nan), marker='.', color=colors[stage], label=stage.capitalize())
        axes[1, 1].set_ylabel('Physical unit-direction difference norm')
        for ax in axes.flat:
            for row in path:
                if not row['measurement_flag']:
                    ax.axvspan(row['forward_index']-.18, row['forward_index']+.18, color='#dddddd', alpha=.4)
            ax.set_yscale('log')
            ax.set_xlabel('Native forward linked-state index (gray: passive surface)')
            ax.grid(alpha=.2)
            ax.legend(fontsize=8, loc='best')
        fig.tight_layout(rect=[0, .09, 1, .94])
        save(fig, 'forward_covariance_extreme', 'DSID 100117, entry 1, component 24: forward covariance and mean differences', pdf)
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
    pages = PdfReader(output/'forward_covariance_report.pdf').pages
    assert len(index) == len(pages) == 3
    for row, page in zip(index, pages):
        with Image.open(output/row['png']) as png:
            png.load()
        original = PdfReader(output/row['pdf']).pages[0]
        assert graphics(original) == graphics(page)
        assert (list(original.mediabox), list(original.cropbox), original.rotation) == (
            list(page.mediabox), list(page.cropbox), page.rotation)
    write_json(output/'validation.json', dict(pages=3, native_state_metrics_recomputed=True,
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
