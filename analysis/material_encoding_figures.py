#!/usr/bin/env python3
"""Plot and validate the complete native Calypso material-encoding control."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

from slope_tables import write_json


DSIDS = ['100044', '100116', '100117']
CONTROLS = ['legacy', 'corrected']
COLORS = ['#c43c39', '#205fb2']
LABELS = ['Legacy material encoding', 'Corrected material encoding']


def read_sources(folder):
    return {dsid: json.loads((folder/'full'/dsid/'summary.json').read_text()) for dsid in DSIDS}


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    data = read_sources(source)
    output.mkdir(parents=True, exist_ok=False)
    index, values = [], {}

    def draw(ax, dsid, family, key, xlabel, log=False):
        notes = []
        for control, color, label in zip(CONTROLS, COLORS, LABELS):
            h = data[dsid][family][control]['histograms'][key]
            values['|'.join([dsid, control, family, key])] = h
            ax.stairs(h['counts'], h['edges'], label=label, color=color)
            notes.append(f"{control}: N={h['raw_entries']}, invalid={h['invalid']}, U/O={h['underflow']}/{h['overflow']}")
        ax.text(.98, .98, '\n'.join(notes), transform=ax.transAxes, ha='right', va='top', fontsize=6)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Entries / bin')
        if log:
            ax.set_yscale('log')
            ax.set_ylim(bottom=.5)
        else:
            ax.set_ylim(bottom=0)
        ax.grid(alpha=.2)

    def save(fig, name, title, pdf, dsid=None, reference=''):
        if dsid:
            n = data[dsid]['coverage']['legacy']['source_events']
            tasks = data[dsid]['coverage']['legacy']['tasks']
            fig.suptitle(f'{dsid}: {title} | {n:,} source events, {tasks} frozen ranges', fontsize=12)
        else:
            fig.suptitle(title, fontsize=12)
        fig.text(.5, .025, 'AfterGeo, isolated smoothing ON | Three-station forward | Only 12 serialized material rotations change',
                 ha='center', fontsize=8)
        fig.text(.5, .009, reference, ha='center', fontsize=7)
        fig.tight_layout(rect=[0, .065, 1, .95])
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf', reference=reference))

    with PdfPages(output/'material_encoding_report.pdf') as pdf:
        for dsid in DSIDS:
            fig, axes = plt.subplots(2, 2, figsize=(12, 8))
            for ax, key, xlabel in zip(axes.flat,
                    ['relative_p', 'signed_qop', 'wrong_signed_qop', 'reco_p'],
                    [r'$(p_{reco}-p_{production})/p_{production}$', r'Signed q/p [GeV$^{-1}$]',
                     r'Wrong-sign q/p [GeV$^{-1}$]', 'Reconstructed p [GeV]']):
                draw(ax, dsid, 'groups', key, xlabel, True)
                ax.set_title(key.replace('_', ' '), fontsize=10)
            axes[0, 0].legend(loc='lower center', fontsize=7)
            save(fig, dsid+'_material_momentum', 'Momentum and signed curvature', pdf, dsid,
                 'Charge-independent truth matching, purity > 0.7, Good95; momentum reference: HepMC production p')

            for residual in [False, True]:
                fig, axes = plt.subplots(2, 3, figsize=(14, 8))
                selections = ['all', 'reco_p_window', 'truth_p_window']
                for row, plane in enumerate(['st1', 'VetoNu']):
                    for col, selection in enumerate(selections):
                        key = '|'.join([plane, selection, 'residual' if residual else 'physical'])
                        xlabel = r'$\Delta\tan(\theta_y)$' if residual else r'$\tan(\theta_y)=p_{y,F}/p_{z,F}$'
                        draw(axes[row, col], dsid, 'slopes', key, xlabel)
                        axes[row, col].set_title(plane+' | '+selection.replace('_', ' '), fontsize=10)
                axes[0, 0].legend(loc='upper left', fontsize=6)
                suffix = 'residual' if residual else 'physical'
                reference = ('Fixed legacy production-p/first-ST1-chord truth model on each physical plane; local G4 truth is not used'
                    if residual else 'S05 / geometryRotation p19 definitions; st1 z_F=0 mm, VetoNu z_F=-3112 mm; p window: 50 < p < 250 GeV')
                save(fig, dsid+'_material_angle_'+suffix, 'Physical angle '+suffix, pdf, dsid, reference)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        metrics = [('conditional_misid', 'Charge misidentification [%]', 100.),
                   ('correct_efficiency', 'Correct-charge truth efficiency [%]', 100.),
                   ('momentum_below_half_production', 'p_reco < 0.5 p_production [%]', 100.),
                   ('angle', r'68% of $|\Delta\tan(\theta_y)|$ [$\mu$rad]', 1e6)]
        for ax, (metric, ylabel, scale) in zip(axes.flat, metrics):
            for offset, control, color, label in zip([-.08, .08], CONTROLS, COLORS, LABELS):
                rows = [data[d]['groups'][control][metric] for d in DSIDS]
                if metric == 'angle':
                    y = np.array([r['absolute_residual_68'] for r in rows])*scale
                    ax.plot(np.arange(3)+offset, y, 'o', color=color, label=label)
                else:
                    y = np.array([r['value'] for r in rows])*scale
                    bounds = np.array([r['interval95'] for r in rows]).T*scale
                    ax.errorbar(np.arange(3)+offset, y, yerr=np.abs(bounds-y), fmt='o', color=color, label=label, capsize=3)
                for dsid, value in zip(DSIDS, rows):
                    values['|'.join([dsid, control, 'metric', metric])] = value
            ax.set_xticks(np.arange(3), DSIDS)
            ax.set_ylabel(ylabel)
            ax.grid(alpha=.2)
        axes[0, 0].legend(fontsize=7)
        save(fig, 'material_joint_metrics', 'Material encoding: joint charge, momentum, angle and efficiency', pdf,
             reference='500,000 source events per sample | Source-event 95% intervals | Native seed priors and all reconstruction cuts fixed')

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ax, dsid in zip(axes, ['100116', '100117']):
            for offset, control, color, label in zip([-.08, .08], CONTROLS, COLORS, LABELS):
                bins = data[dsid]['groups'][control]['truth_p_bins']
                for i, row in enumerate(bins):
                    rate = row['conditional_misid']
                    if rate['value'] is None:
                        continue
                    y = 100*rate['value']
                    errors = np.abs(100*np.array(rate['interval95'])-y).reshape(2, 1)
                    ax.errorbar([i+offset], [y], yerr=errors, fmt='o', color=color, capsize=3,
                                label=label if i == 0 else None)
                values['|'.join([dsid, control, 'metric', 'truth_p_bins'])] = bins
            ax.set_title(dsid)
            ax.set_xticks(range(len(bins)), [f"{r['low_GeV']:g}-{r['high_GeV']:g}" for r in bins], rotation=25)
            ax.set_xlabel('Truth production p [GeV]')
            ax.set_ylabel('Conditional charge misidentification [%]')
            ax.grid(alpha=.2)
        axes[0].legend(fontsize=7)
        save(fig, 'material_truth_p_misid', 'Material encoding: muon charge error versus truth momentum', pdf,
             reference='500,000 source events per charge | Same charge-independent selection and truth bins | Source-event 95% intervals')

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        status = ['Missing', 'Correct', 'Wrong', 'Invalid', 'Impure', 'Quality', 'Nonfid.']
        for ax, dsid in zip(axes, DSIDS):
            matrix = np.array(data[dsid]['paired']['migration'])
            ax.imshow(np.log1p(matrix), cmap='Blues')
            for i, j in zip(*np.nonzero(matrix)):
                ax.text(j, i, str(matrix[i, j]), ha='center', va='center', fontsize=6,
                        color='white' if matrix[i, j] > matrix.max()*.05 else 'black')
            ax.set_xticks(range(7), status, rotation=55, fontsize=7)
            ax.set_yticks(range(7), status, fontsize=7)
            ax.set_xlabel('Corrected encoding')
            ax.set_ylabel('Legacy encoding')
            ax.set_title(dsid)
            values[dsid+'|paired|migration|counts'] = matrix.tolist()
        save(fig, 'material_selection_migration', 'Material encoding: paired truth-status migration', pdf,
             reference='All saved primary/secondary truth rows retained; source file and entry define event identity; shading is log(1+count)')

    write_json(output/'figure_index.json', index)
    write_json(output/'plot_values.json', values)
    write_json(output/'plot_source.json', dict(source=str(source), pages=len(index)))
    differences = {}
    for control in CONTROLS:
        a, b = [data[d]['groups'][control]['conditional_misid'] for d in ['100117', '100116']]
        difference = a['value']-b['value']
        error = float(np.hypot(a['standard_error'], b['standard_error']))
        differences[control] = dict(raw=dict(mu_minus_minus_mu_plus=difference,
            interval95=[difference-1.95996398454*error, difference+1.95996398454*error]))
        for split in ['all', 'held_out']:
            a, b = [data[d]['groups'][control]['production_overlap'][split]['conditional_misid']
                    for d in ['100117', '100116']]
            difference = a['value']-b['value']
            error = float(np.hypot(a['standard_error'], b['standard_error']))
            differences[control]['production_overlap_'+split] = dict(mu_minus_minus_mu_plus=difference,
                interval95=[difference-1.95996398454*error, difference+1.95996398454*error])
    write_json(output/'charge_differences.json', dict(differences=differences,
        interpretation='Independent source samples for the two charges. Weighted intervals are conditional on the previously frozen production-only propensity model.'))


def validate(source, figures, pdf_tools):
    sys.path.insert(0, str(pdf_tools))
    from pypdf import PdfReader
    from PIL import Image
    from validate_charge_gun_figures import normalized_graphics

    data = read_sources(source)
    values = json.loads((figures/'plot_values.json').read_text())
    for key, value in values.items():
        dsid, control, kind, field = key.split('|', 3)
        if control == 'paired':
            expected = data[dsid]['paired']['migration']
        elif kind == 'metric':
            expected = data[dsid]['groups'][control][field]
        else:
            expected = data[dsid][kind][control]['histograms'][field]
            assert sum(value['counts'])+value['underflow']+value['overflow'] == value['finite']
            assert value['sumw'] == value['sumw2'] == value['counts']
        assert value == expected, key
    index = json.loads((figures/'figure_index.json').read_text())
    combined = PdfReader(figures/'material_encoding_report.pdf')
    assert len(combined.pages) == len(index) == 12
    pages = []
    for row, page in zip(index, combined.pages):
        with Image.open(figures/row['png']) as png:
            png.load()
        single = PdfReader(figures/row['pdf'])
        assert len(single.pages) == 1
        original = single.pages[0]
        assert normalized_graphics(page) == normalized_graphics(original)
        assert list(page.mediabox) == list(original.mediabox)
        assert list(page.cropbox) == list(original.cropbox)
        assert page.rotation == original.rotation
        pages.append(dict(figure=row['name'], PNG_decoded=True, vector_graphics_and_boxes_equal=True))
    write_json(figures/'validation.json', dict(pages=12, plotted_values_exact=len(values),
        page_checks=pages, visual_inspection='Pending separate image review'))
    print('Verified material figures:', len(values), 'plotted source records, 12 pages')


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
