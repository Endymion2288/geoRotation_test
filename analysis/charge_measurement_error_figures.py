#!/usr/bin/env python3
"""Show actual Calypso measurement-coordinate and covariance controls."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

from slope_tables import write_json

MODES = ['native', 'centroid_fixed', 'ckf_quarter', 'ckf_four', 'ckf_cluster']
LABELS = ['Native', 'SiHit centroid\nfixed seed', 'CKF V x 0.25', 'CKF V x 4', 'CKF cluster V']


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    events = json.loads((source/'events.json').read_text())
    coordinates = json.loads((source/'coordinates.json').read_text())
    assert summary['events'] == len(events) == 66
    return dict(summary=summary, pure=[r for r in events if r['cohort'] == 'pure_wrong'],
                coordinates=coordinates)


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
        fig.text(.5, .028, 'Actual geoRotation Calypso | Corrected material encoding | Isolated smoothing ON', ha='center', fontsize=8)
        fig.text(.5, .009, 'Frozen cohort: 6 originally pure wrong signs + 60 normal controls; outcomes are not population rates', ha='center', fontsize=8)
        fig.tight_layout(rect=[0, .065, 1, .95])
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
            reference='S06 real Calypso measurement mechanism; frozen outcome-selected cohort'))

    with PdfPages(output/'measurement_controls_report.pdf') as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(13, 8))
        for column, dsid in enumerate(['100116', '100117']):
            for row, cohort in enumerate(['pure_wrong', 'normal']):
                ax = axes[row, column]
                bottom = np.zeros(len(MODES))
                statuses = [('correct', '#227f4b'), ('wrong', '#bd3e3e'), ('quality_failed', '#dc9e32'),
                            ('other', '#777777')]
                for status, color in statuses:
                    y = []
                    for mode in MODES:
                        group = data['summary']['groups']['|'.join([dsid, cohort, mode])]
                        y.append(group.get(status, 0) if status != 'other' else
                            group['events']-sum(group.get(s, 0) for s in ['correct', 'wrong', 'quality_failed']))
                    bars = ax.bar(np.arange(len(MODES)), y, bottom=bottom, color=color, label=status.replace('_', ' '))
                    ax.bar_label(bars, labels=[str(n) if n else '' for n in y], label_type='center', fontsize=8)
                    bottom += y
                ax.set_xticks(np.arange(len(MODES)), LABELS, fontsize=8)
                ax.set_ylabel('Events')
                ax.set_title(dsid+' | '+('Original pure wrong' if row == 0 else 'Frozen normal'))
                ax.set_ylim(0, max(bottom)*1.18)
                ax.grid(axis='y', alpha=.2)
        axes[0, 0].legend(fontsize=8, ncol=2)
        save(fig, 'measurement_outcomes', 'Charge response to measured coordinates and actual CKF variance', pdf)

        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
        for cohort, color, label in [('normal', '#777777', 'Normal controls'), ('pure_wrong', '#b83c35', 'Original pure wrong')]:
            rows = [r for r in data['coordinates'] if r['cohort'] == cohort and 'measured_minus_primary_centroid_mm' in r]
            x = [r['n_strips'] for r in rows]
            y = [1000*abs(r['measured_minus_primary_centroid_mm']) for r in rows]
            axes[0].scatter(x, y, color=color, s=18, alpha=.55, label=label)
            axes[1].scatter(x, [1000*r['cluster_sigma_mm'] for r in rows], color=color, s=18, alpha=.55, label=label)
        axes[0].set_ylabel('Absolute cluster - primary SiHit centroid [um]')
        axes[1].set_ylabel('Measurement standard deviation [um]')
        axes[1].axhline(23.1, color='#1f5aa6', linestyle='--', label='Native CKF: 23.1 um')
        for ax in axes:
            ax.set_xlabel('Number of cluster strips')
            ax.grid(alpha=.2)
            ax.legend(fontsize=8)
        axes[0].set_title('Native selected measurement states')
        axes[1].set_title('Cluster covariance versus fixed CKF error')
        save(fig, 'measurement_coordinates', 'Actual strip geometry, SiHit positions and assigned errors', pdf)

        fig, axes = plt.subplots(2, 1, figsize=(13, 7))
        for ax, dsid in zip(axes, ['100116', '100117']):
            charge = 1 if dsid == '100116' else -1
            rows = [r for r in data['pure'] if r['dsid'] == dsid]
            for row in rows:
                y = [charge*row['modes'][mode]['final_qop_GeVinv'] if row['modes'][mode]['final_qop_GeVinv'] is not None else np.nan for mode in MODES]
                ax.plot(np.arange(len(MODES)), y, 'o-', label='Entry '+str(row['sample_entry']), alpha=.85)
                for i, mode in enumerate(MODES):
                    if row['modes'][mode]['status'] == 'quality_failed' and np.isfinite(y[i]):
                        ax.scatter(i, y[i], marker='x', s=100, color='black', zorder=5)
            ax.axhline(0, color='gray', linewidth=1)
            ax.set_xticks(np.arange(len(MODES)), LABELS, fontsize=8)
            ax.set_ylabel(r'$q_{truth}(q/p)_{final}$ [GeV$^{-1}$]')
            ax.set_title(dsid+' | Positive values have the correct sign; black x marks quality failure')
            ax.legend(fontsize=8, ncol=2)
            ax.grid(alpha=.2)
        save(fig, 'measurement_pure_curvature', 'Six pure wrong events with their original complete seeds', pdf)
    write_json(output/'plot_values.json', data)
    write_json(output/'figure_index.json', index)


def validate(source, output, pdf_tools):
    from PIL import Image
    sys.path.insert(0, str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics
    assert json.loads((output/'plot_values.json').read_text()) == values(source)
    index = json.loads((output/'figure_index.json').read_text())
    report = PdfReader(output/'measurement_controls_report.pdf')
    assert len(index) == len(report.pages) == 3
    for row, page in zip(index, report.pages):
        with Image.open(output/row['png']) as png:
            png.load()
        single = PdfReader(output/row['pdf']).pages[0]
        assert normalized_graphics(single) == normalized_graphics(page)
        assert list(single.mediabox) == list(page.mediabox) and list(single.cropbox) == list(page.cropbox)
        assert single.rotation == page.rotation
    write_json(output/'validation.json', dict(pages=3, plotted_source_values_exact=True,
        PNG_decoded=True, normalized_vector_graphics_and_boxes_equal=True,
        visual_inspection='Pending separate full-size PNG review'))


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
