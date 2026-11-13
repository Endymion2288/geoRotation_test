#!/usr/bin/env python3
"""Show joint charge-control angle and process-resource measurements."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from charge_history import rate
from slope_tables import write_json


MODES = ['native', 'smoother', 'updater', 'both']
DSIDS = ['100116', '100117']
LABELS = [r'$\mu^+$', r'$\mu^-$']
COLORS = ['#424242', '#ef6c00', '#1565c0', '#388e3c']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    angles = json.loads((args.root/'refit_regression_v1/summary.json').read_text())
    resources = json.loads((args.root/'resource_analysis_v1/summary.json').read_text())
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2})
    index, counts = [], dict(angles=angles, resources=resources)
    with PdfPages(args.output/'refit_joint_regression_report.pdf') as pdf:
        def save(fig, name, title, footer):
            fig.suptitle(title, y=.99, fontsize=14)
            fig.text(.5, .92, 'beforeGeo | 100,000 source events per charge | Native CKF unchanged | Isolated smoothing ON', ha='center', fontsize=9)
            fig.text(.5, .018, footer, ha='center', fontsize=8)
            fig.tight_layout(rect=(.02, .055, .98, .88))
            fig.savefig(args.output/(name+'.png'), dpi=150)
            fig.savefig(args.output/(name+'.pdf'))
            pdf.savefig(fig)
            plt.close(fig)
            index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf', reference='S06 Q5 joint regression of refit controls',
                definition=footer, bin_counts='bin_counts.json'))

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        for column, (dsid, label) in enumerate(zip(DSIDS, LABELS)):
            groups = [angles['groups'][dsid+'|'+mode]['angles'] for mode in MODES]
            axes[0, column].plot(range(4), [g['absolute_residual_68']*1e6 for g in groups], 'o-', color='#1565c0')
            for i, (group, color) in enumerate(zip(groups, COLORS)):
                value = rate(group['absolute_residual_above_1mrad'], group['finite'])
                low, high = np.asarray(value['intervals']['95'])*100
                y = 100*value['value']
                axes[1, column].errorbar(i, y, yerr=[[y-low], [high-y]], marker='o', capsize=4, color=color)
            axes[0, column].set(title=label, ylabel=r'68% quantile of $|\Delta s_y|$ [$\mu$rad]')
            axes[1, column].set(ylabel=r'$|\Delta s_y|>1$ mrad / finite selected [%]', ylim=(0, None))
            for ax in axes[:, column]:
                ax.set_xticks(range(4), MODES)
        save(fig, 'refit_angle_core_and_tail', 'The updater-only charge improvement has an angular-tail cost',
             'Same-plane st1 SiHit-chord truth model with production p; unbinned core and 95% exact tail intervals; mode-specific selected samples')

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        for column, (dsid, label) in enumerate(zip(DSIDS, LABELS)):
            group = angles['groups'][dsid+'|both']
            for row, observable in enumerate(['residual_fine', 'residual_wide']):
                ax = axes[row, column]
                for key, name, color in [('native', 'Native', '#424242'), ('control', 'Regular updater + smoother', '#388e3c')]:
                    hist = group['common_same_candidate_angles'][key][observable]
                    ax.stairs(np.asarray(hist['counts'])/hist['raw_entries'], np.asarray(hist['edges'])*1e6,
                              color=color, label=name)
                ax.set(title=label+f" | {group['same_candidate_and_ordered_hits']:,} same candidates and hits",
                       xlabel=r'$s_y^{reco}-s_y^{truth\ model}$ [$\mu$rad]', ylabel='Fraction / bin')
                if row == 1:
                    ax.set_yscale('log')
                ax.legend()
        save(fig, 'refit_common_candidate_angles', 'Joint chart control improves the same native candidates',
             'Common selected truth, identical CKF candidate ordinal and ordered hit IDs; actual Calypso projections to z_F=0 mm; flow bins retained')

        fig, axes = plt.subplots(1, 2, figsize=(12, 5.8))
        for dsid, label, marker in zip(DSIDS, LABELS, ['o', 's']):
            rows = [resources['refit'][dsid+'|'+mode] for mode in MODES]
            axes[0].plot(range(4), [row['process_cpu_ms_per_source_event'] for row in rows], marker=marker, label=label)
            axes[1].plot(range(4), [row['maximum_process_RSS_KiB']/1024 for row in rows], marker=marker, label=label)
        for ax in axes:
            ax.set_xticks(range(4), MODES)
            ax.legend()
        axes[0].set(ylabel='Complete tracking process CPU / event [ms]')
        axes[1].set(ylabel='Maximum worker process RSS [MiB]')
        save(fig, 'refit_cpu_and_memory', 'Measured process CPU and memory in the refit controls',
             'GNU time user+system includes process initialization; workers and batches differ; these descriptive measurements do not establish a controlled speed ratio')
    write_json(args.output/'bin_counts.json', counts)
    write_json(args.output/'figure_index.json', index)
    (args.output/'figure_index.md').write_text('\n'.join(
        '- ['+row['name']+']('+row['png']+') · [PDF]('+row['pdf']+')' for row in index)+'\n')
    print(json.dumps(dict(figures=len(index), output=str(args.output)), indent=2))


if __name__ == '__main__':
    main()
