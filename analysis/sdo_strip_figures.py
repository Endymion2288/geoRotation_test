#!/usr/bin/env python3
"""Plot charge outcomes and actual local coordinates for native strip controls."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

from slope_tables import write_json

MODES = ['native', 'sdo_strips_fixed', 'centroid_fixed']
LABELS = ['Native', 'SDO-supported strips\nfixed seed', 'SiHit centroid\nfixed seed']


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    events = json.loads((source/'events.json').read_text())
    coordinates = json.loads((source/'coordinates.json').read_text())
    evidence = json.loads((source/'cluster_evidence.json').read_text())
    assert summary['events'] == len(events) == 66
    assert len(evidence) == 6
    pure = []
    for row in events:
        if row['cohort'] == 'pure_wrong':
            pure.append(dict(dsid=row['dsid'], sample_entry=row['sample_entry'], modes={
                mode: {key: row['modes'][mode][key] for key in ['status', 'final_qop_GeVinv', 'same_ordered_hits']}
                for mode in MODES}))
    histograms = {}
    edges = np.linspace(-.32, .32, 81)
    for cohort in ['pure_wrong', 'normal']:
        rows = [row for row in coordinates if row['cohort'] == cohort and 'sdo_minus_primary_centroid_mm' in row]
        for mode, key in [('native', 'measured_minus_primary_centroid_mm'), ('sdo_strips_fixed', 'sdo_minus_primary_centroid_mm')]:
            data = np.array([row[key] for row in rows])
            assert np.isfinite(data).all()
            bins = np.histogram(data[(data >= edges[0]) & (data < edges[-1])], edges)[0]
            underflow, overflow = int((data < edges[0]).sum()), int((data >= edges[-1]).sum())
            assert bins.sum()+underflow+overflow == len(data)
            histograms[cohort+'|'+mode] = dict(counts=bins.tolist(), entries=len(data),
                underflow=underflow, overflow=overflow)
    return dict(summary=summary, pure=pure, bin_edges_mm=edges.tolist(), histograms=histograms,
                cluster_evidence=evidence)


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    data = values(source)
    output.mkdir(exist_ok=False)
    index = []

    def save(fig, name, title, pdf, top=.95):
        fig.suptitle(title, fontsize=13)
        fig.text(.5, .032, 'Actual geoRotation Calypso | Corrected material map | Isolated smoothing ON | Original complete seeds', ha='center', fontsize=8)
        fig.text(.5, .012, f"66 frozen source events from {data['summary']['source_files']} RDO files: 6 original pure wrong + 60 normal; outcome-selected mechanism cohort", ha='center', fontsize=8)
        fig.tight_layout(rect=[0, .075, 1, top])
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
            reference='S06 native Calypso SDO-supported strip-centre control; strip quantization retained'))

    with PdfPages(output/'sdo_strip_controls_report.pdf') as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(13, 8))
        for column, dsid in enumerate(['100116', '100117']):
            for row, cohort in enumerate(['pure_wrong', 'normal']):
                ax = axes[row, column]
                bottom = np.zeros(len(MODES))
                for status, color in [('correct', '#227f4b'), ('wrong', '#bd3e3e'), ('quality_failed', '#dc9e32'), ('other', '#777777')]:
                    heights = []
                    for mode in MODES:
                        group = data['summary']['groups']['|'.join([dsid, cohort, mode])]
                        heights.append(group.get(status, 0) if status != 'other' else
                            group['events']-sum(group.get(s, 0) for s in ['correct', 'wrong', 'quality_failed']))
                    bars = ax.bar(np.arange(len(MODES)), heights, bottom=bottom, color=color, label=status.replace('_', ' '))
                    ax.bar_label(bars, labels=[str(n) if n else '' for n in heights], label_type='center', fontsize=8)
                    bottom += heights
                ax.set_xticks(np.arange(len(MODES)), LABELS, fontsize=8)
                ax.set_ylabel('Events')
                ax.set_title(dsid+' | '+('Original pure wrong' if cohort == 'pure_wrong' else 'Frozen normal'))
                ax.set_ylim(0, max(bottom)*1.2)
                ax.grid(axis='y', alpha=.2)
        axes[0, 0].legend(fontsize=8, ncol=2)
        save(fig, 'sdo_strip_outcomes', 'Charge response to SDO-supported strip centres and SiHit centroids', pdf)

        fig, axes = plt.subplots(2, 1, figsize=(13, 7))
        for ax, dsid in zip(axes, ['100116', '100117']):
            for row in [r for r in data['pure'] if r['dsid'] == dsid]:
                charge = 1 if dsid == '100116' else -1
                y = [charge*row['modes'][mode]['final_qop_GeVinv'] if row['modes'][mode]['final_qop_GeVinv'] is not None else np.nan for mode in MODES]
                ax.plot(np.arange(len(MODES)), y, 'o-', label='Entry '+str(row['sample_entry']))
                for i, mode in enumerate(MODES):
                    if row['modes'][mode]['status'] == 'quality_failed' and np.isfinite(y[i]):
                        ax.scatter(i, y[i], color='black', marker='x', s=100, zorder=5)
            ax.axhline(0, color='gray', linewidth=1)
            ax.set_xticks(np.arange(len(MODES)), LABELS, fontsize=8)
            ax.set_ylabel(r'$q_{truth}(q/p)_{final}$ [GeV$^{-1}$]')
            ax.set_title(dsid+' | Positive values have the correct sign; black x marks quality failure')
            ax.legend(fontsize=8, ncol=2)
            ax.grid(alpha=.2)
        save(fig, 'sdo_strip_curvature', 'Six original pure wrong events with fixed complete seeds', pdf)

        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
        for ax, cohort in zip(axes, ['pure_wrong', 'normal']):
            for mode, color, label in [('native', '#777777', 'Native cluster'), ('sdo_strips_fixed', '#2459a5', 'SDO-supported strip centres')]:
                row = data['histograms'][cohort+'|'+mode]
                ax.stairs(row['counts'], 1000*np.asarray(data['bin_edges_mm']), color=color,
                    linewidth=1.5, label=f"{label} (N={row['entries']})")
            ax.set_xlabel('Measured X - primary SiHit centroid [um]')
            ax.set_ylabel('Native selected measurement states / 8 um')
            ax.set_yscale('log')
            ax.set_ylim(bottom=.7)
            ax.set_title('Original pure wrong' if cohort == 'pure_wrong' else 'Frozen normal')
            ax.legend(fontsize=8)
            ax.grid(alpha=.2)
        save(fig, 'sdo_strip_coordinates', 'Actual coordinate changes on the fixed native measurement paths', pdf)

        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        fig, axes = plt.subplots(3, 2, figsize=(13, 10))
        for ax, row in zip(axes.flat, data['cluster_evidence']):
            cluster = row['native_cluster']
            origin = cluster['majority_centroid_local_mm']
            centres = np.array([channel['centre_local_mm'][0] for channel in cluster['channels']])
            pitch = cluster['phi_pitch_mm']
            limits = centres.min()-pitch, centres.max()+pitch
            for channel in cluster['channels']:
                x = channel['centre_local_mm'][0]-origin[0]
                color = '#5b9fc1' if channel['deposits'] else '#b3b3b3'
                ax.axvspan(1000*(x-.5*pitch), 1000*(x+.5*pitch), color=color, alpha=.2)
            for hit in cluster['same_wafer_hits']:
                points = np.array([hit['start_local_mm'], hit['end_local_mm']])
                if points[:, 0].max() < limits[0] or points[:, 0].min() > limits[1]:
                    continue
                color = '#d67b18' if hit['barcode'] == 0 else '#222222'
                ax.plot(1000*(points[:, 0]-origin[0]), points[:, 1]-origin[1],
                        color=color, linewidth=1.4, marker='.', markersize=2)
            ax.axvline(1000*(cluster['native_local_mm'][0]-origin[0]), color='#222222', linestyle='--')
            ax.axvline(1000*(row['supported_used_x_mm']-origin[0]), color='#25834a', linestyle=':')
            ax.set_xlim(1000*(limits[0]-origin[0]), 1000*(limits[1]-origin[0]))
            ax.set_xlabel('Sensor local X relative to primary centroid [um]', fontsize=8)
            ax.set_ylabel('Relative sensor local Y [mm]', fontsize=8)
            ax.set_title(f"{row['dsid']} / {row['sample_entry']} | {len(cluster['channels'])} strips | {row['features']['barcode_zero_sihits']} barcode-0 SiHits", fontsize=9)
            ax.tick_params(labelsize=8)
            ax.grid(alpha=.2)
        fig.legend(handles=[Patch(facecolor='#5b9fc1', alpha=.3, label='Strip with SDO contribution'),
            Patch(facecolor='#b3b3b3', alpha=.3, label='Strip without SDO contribution'),
            Line2D([], [], color='#222222', label='G4 SiHit with truth link'),
            Line2D([], [], color='#d67b18', label='G4 SiHit, barcode 0'),
            Line2D([], [], color='#222222', linestyle='--', label='Native measured X'),
            Line2D([], [], color='#25834a', linestyle=':', label='SDO-supported X')],
            ncol=3, loc='upper center', bbox_to_anchor=(.5, .943), fontsize=8)
        save(fig, 'sdo_strip_g4_evidence', 'Same-wafer G4 SiHit segments in each cluster X window', pdf, top=.865)
    write_json(output/'plot_values.json', data)
    write_json(output/'figure_index.json', index)


def validate(source, output, pdf_tools):
    from PIL import Image
    sys.path.insert(0, str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics
    assert json.loads((output/'plot_values.json').read_text()) == values(source)
    index = json.loads((output/'figure_index.json').read_text())
    report = PdfReader(output/'sdo_strip_controls_report.pdf')
    assert len(index) == len(report.pages) == 4
    for row, page in zip(index, report.pages):
        with Image.open(output/row['png']) as image:
            image.load()
        single = PdfReader(output/row['pdf']).pages[0]
        assert normalized_graphics(single) == normalized_graphics(page)
        assert list(single.mediabox) == list(page.mediabox) and list(single.cropbox) == list(page.cropbox)
        assert single.rotation == page.rotation
    write_json(output/'validation.json', dict(pages=4, plotted_values_exact=True,
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
