#!/usr/bin/env python3
"""Render joint validation metrics and the frozen reference-style distributions."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_tables import write_json


MODES = ['3forward', '3backward', '3merge', '4forward', '4backward', '4merge']
LABELS = ['3F', '3B', '3M', '4F', '4B', '4M']
SPECIES = {'100116': r'$\mu^+$', '100117': r'$\mu^-$',
           '110033': r"$A'\to e^+e^-,\ \epsilon=10^{-5}$",
           '110051': r"$A'\to e^+e^-,\ \epsilon=10^{-4}$"}
COLORS = {'native': '#c43c39', 'both': '#205fb2'}
LEGEND = {'native': 'Native refit', 'both': 'Regular update + smoother'}


def save(fig, output, name, pdf, index, title):
    fig.savefig(output/(name+'.png'), dpi=150)
    fig.savefig(output/(name+'.pdf'))
    pdf.savefig(fig)
    index.append(dict(name=name, title=title, png=name+'.png', pdf=name+'.pdf',
        reference='S06 Q5 native/regular-refit joint regression; frozen source-event selections'))
    plt.close(fig)


def header(fig, title, count, coverage):
    fig.suptitle(title+f' | {count:,} source events', fontsize=13)
    fig.text(.5, .94, coverage, ha='center', fontsize=9)
    fig.text(.5, .012, 'Actual beforeGeo Calypso | isolated CKF smoothing ON | F/B/M: forward/backward/merge | source-event uncertainties',
             ha='center', fontsize=8)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    summary = json.loads(args.input.read_text())
    manifest = json.loads(Path(summary['source_manifest']).read_text())
    coverage = {}
    for dsid, ranges in manifest['ranges'].items():
        files = sorted({r['file_index'] for r in ranges})
        entries = sum(r['events'] for r in ranges)//len(files)
        coverage[dsid] = (f'{len(files)} input files | file indices {files[0]}-{files[-1]} | '
            f'entries [0,{entries:,}) per file')
    args.output.mkdir(exist_ok=False)
    index, counts = [], {}
    with PdfPages(args.output/'joint_regression_report.pdf') as pdf:
        for dsid, label in SPECIES.items():
            source_count = summary['groups'][dsid+'|3forward|native']['source_events']
            fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
            metrics = [('conditional_misid', 'Charge misidentification [%]'),
                       ('correct_efficiency', 'Correct charge / fiducial truth [%]'),
                       ('momentum_below_half_production', r'$p_{reco}/p_{production}<0.5$ [%]'),
                       ('angle_above_1mrad', r'$|\Delta\tan\theta_y|>1$ mrad [%]')]
            for ax, (metric, title) in zip(axes.flat, metrics):
                for control in COLORS:
                    rows = [summary['groups'][f'{dsid}|{mode}|{control}'][metric] for mode in MODES]
                    values = np.array([r['value'] for r in rows], dtype=float)*100
                    interval = np.array([r['interval95'] for r in rows], dtype=float).T*100
                    offset = -.08 if control == 'native' else .08
                    ax.errorbar(np.arange(6)+offset, values, yerr=np.maximum(0, [values-interval[0], interval[1]-values]),
                                marker='o', ms=4, capsize=2, color=COLORS[control], label=LEGEND[control])
                    counts[f'{dsid}|{metric}|{control}'] = rows
                ax.set_title(title, fontsize=10)
                ax.set_xticks(np.arange(6), LABELS)
                ax.set_ylim(bottom=0)
                ax.grid(alpha=.2)
            for control in COLORS:
                widths = [summary['groups'][f'{dsid}|{mode}|{control}']['angle']['absolute_residual_68']*1e6 for mode in MODES]
                axes[1, 1].plot(LABELS, widths, 'o-', color=COLORS[control], label=LEGEND[control])
                fallback = []
                for mode in MODES:
                    row = summary['refit_counters'][f'{dsid}|{mode}|{control}']
                    if mode.endswith('merge'):
                        row = row['merged']['selected_merged_output']
                    else:
                        row = row['selected']
                    fallback.append(100*row.get('fallback', 0)/row['total'] if row['total'] else np.nan)
                axes[1, 2].plot(LABELS, fallback, 'o-', color=COLORS[control], label=LEGEND[control])
                counts[f'{dsid}|angle68|{control}'] = widths
                counts[f'{dsid}|selected_fallback|{control}'] = fallback
            axes[1, 1].set_title(r'68% quantile of $|\Delta\tan\theta_y|$ [$\mu$rad]', fontsize=10)
            axes[1, 2].set_title('Selected refit fallback [%]', fontsize=10)
            for ax in axes[1, 1:]:
                ax.set_ylim(bottom=0)
                ax.grid(alpha=.2)
            axes[0, 0].legend(fontsize=8)
            header(fig, f'{dsid}: {label} | joint regression', source_count, coverage[dsid])
            fig.tight_layout(rect=[0, .04, 1, .91])
            save(fig, args.output, dsid+'_joint_metrics', pdf, index, 'Joint metrics')
            for field, title, xlabel, scale in [
                    ('relative_p', 'Momentum response', r'$(p_{reco}-p_{production})/p_{production}$', 1.),
                    ('angle', r'Angle residual at $z_F=0$ mm', r'$\Delta\tan\theta_y$ [mrad]', 1000.),
                    ('signed_qop', 'Signed curvature', r'$q/p$ [GeV$^{-1}$]', 1.)]:
                fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
                for ax, mode in zip(axes.flat, MODES):
                    for control in COLORS:
                        row = summary['groups'][f'{dsid}|{mode}|{control}']
                        h = row['angle']['residual_fine'] if field == 'angle' else row['histograms'][field]
                        edges = np.asarray(h['edges'])*scale
                        values = np.asarray(h['counts'])
                        assert int(values.sum())+h['underflow']+h['overflow'] == h['finite']
                        assert sum(h['sumw2']) == int(values.sum())
                        ax.stairs(values, edges, color=COLORS[control],
                            label=LEGEND[control]+f"; N={h['finite']:,}, UF/OF={h['underflow']}/{h['overflow']}")
                        counts[f'{dsid}|{mode}|{control}|{field}'] = h
                    ax.set_title(mode)
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel('Entries / bin')
                    ax.grid(alpha=.2)
                    ax.legend(fontsize=6, loc='upper right')
                header(fig, f'{dsid}: {label} | {title}', source_count, coverage[dsid])
                fig.tight_layout(rect=[0, .04, 1, .91])
                save(fig, args.output, dsid+'_'+field, pdf, index, title)
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        for ax, (dsid, label) in zip(axes.flat, SPECIES.items()):
            for mode, text in zip(MODES, LABELS):
                row = summary['comparisons'][dsid+'|'+mode]
                x, y = row['correct_efficiency_difference'], row['conditional_misid_difference']
                counts[f'{dsid}|{mode}|paired_changes'] = dict(correct_efficiency=x, conditional_misid=y)
                ax.errorbar(100*x['value'], 100*y['value'],
                    xerr=100*np.array([[x['value']-x['interval95'][0]], [x['interval95'][1]-x['value']]]),
                    yerr=100*np.array([[y['value']-y['interval95'][0]], [y['interval95'][1]-y['value']]]),
                    marker='o', ms=4, capsize=2)
                ax.annotate(text, (100*x['value'], 100*y['value']), xytext=(4, 4), textcoords='offset points', fontsize=8)
            ax.axhline(0, color='gray', lw=.7)
            ax.axvline(0, color='gray', lw=.7)
            ax.set_title(dsid+': '+label+'\n'+coverage[dsid], fontsize=9)
            ax.set_xlabel('Correct-efficiency change [percentage points]')
            ax.set_ylabel('MisID change [percentage points]')
            ax.grid(alpha=.2)
        fig.suptitle('Paired changes: regular joint refit minus native refit | 95% source-event intervals')
        fig.tight_layout(rect=[0, 0, 1, .95])
        save(fig, args.output, 'paired_charge_efficiency', pdf, index, 'Paired charge and efficiency')
        fig, axes = plt.subplots(2, 2, figsize=(12, 7.5))
        for ax, (dsid, label) in zip(axes.flat, SPECIES.items()):
            for control in COLORS:
                values = [summary['resources'][f'{dsid}|{mode}|{control}']['process_cpu_ms_per_source_event'] for mode in MODES]
                ax.plot(LABELS, values, 'o-', color=COLORS[control], label=LEGEND[control])
                counts[f'{dsid}|CPU_ms|{control}'] = values
            ax.set_title(dsid+': '+label+'\n'+coverage[dsid], fontsize=9)
            ax.set_ylabel('Observed process CPU [ms / source event]')
            ax.grid(alpha=.2)
            ax.set_ylim(bottom=0)
        axes[0, 0].legend(fontsize=8)
        fig.suptitle('Observed CPU including initialization | batches ran on different workers')
        fig.tight_layout(rect=[0, 0, 1, .95])
        save(fig, args.output, 'observed_cpu', pdf, index, 'Observed process CPU')
    write_json(args.output/'bin_counts.json', counts)
    write_json(args.output/'figure_index.json', index)
    (args.output/'figure_index.md').write_text('\n'.join(
        f'{i+1}. [{row["name"]}]({row["pdf"]})' for i, row in enumerate(index))+'\n')
    write_json(args.output/'plot_source.json', dict(summary=str(args.input), pages=len(index),
        axes='All figures are English; raw histogram entries; relative-p uses the reference 100 bins [-1.2,1.2]',
        angle_reference=summary['reference']))


if __name__ == '__main__':
    main()
