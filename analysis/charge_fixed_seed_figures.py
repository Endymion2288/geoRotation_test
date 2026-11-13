#!/usr/bin/env python3
"""Show measured Calypso centroid controls with the complete native seeds."""
import argparse
from collections import Counter
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_tables import write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    rows = json.loads((root/'fixed_seed_analysis_v3/events.json').read_text())
    assert len(rows) == 66
    prior = {(row['dsid'], row['sample_entry']): row for row in
             json.loads((root/'seed_qop_analysis_v1/events.json').read_text())}
    index, counts = [], {}
    args.output.mkdir(parents=True, exist_ok=False)
    labels = ['Native', 'Centroid; free seed', 'Centroid; native seed']
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2})
    with PdfPages(args.output/'fixed_seed_report.pdf') as pdf:
        def save(fig, name, title, footer):
            fig.suptitle(title, y=.99, fontsize=14)
            fig.text(.5, .93, 'Actual afterGeo Calypso | 6 native purity-one wrong signs + 60 frozen normal controls | Isolated smoothing ON',
                     ha='center', fontsize=9)
            fig.text(.5, .018, footer, ha='center', fontsize=8)
            fig.tight_layout(rect=(.02, .06, .98, .89))
            fig.savefig(args.output/(name+'.png'), dpi=150)
            fig.savefig(args.output/(name+'.pdf'))
            pdf.savefig(fig)
            plt.close(fig)
            index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
                reference='S06 Q4: actual measurement changes with native complete seed control',
                selection='Frozen targeted cohorts; quality failures retained separately',
                bin_counts='bin_counts.json', definition=footer))

        fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
        for column, (dsid, charge) in enumerate([('100116', r'$\mu^+$'), ('100117', r'$\mu^-$')]):
            for row_index, cohort in enumerate(['pure_wrong', 'normal']):
                ax = axes[row_index, column]
                chosen = [row for row in rows if row['dsid'] == dsid and row['cohort'] == cohort]
                modes = [Counter({'wrong' if cohort == 'pure_wrong' else 'correct': len(chosen)}),
                         Counter(row['free_seed_centroid_outcome'] for row in chosen),
                         Counter(row['outcome'] for row in chosen)]
                bottom = np.zeros(3)
                for status, color in [('correct', '#43a047'), ('wrong', '#e53935'), ('quality_failed', '#bdbdbd')]:
                    height = np.array([mode[status] for mode in modes])
                    ax.bar(range(3), height, bottom=bottom, color=color, label=status.replace('_', ' '))
                    for i, n in enumerate(height):
                        if n:
                            ax.text(i, bottom[i]+n/2, str(n), ha='center', va='center')
                    bottom += height
                np.testing.assert_array_equal(bottom, len(chosen))
                ax.set(title=charge+' | '+('Native pure wrong' if cohort == 'pure_wrong' else 'Normal controls'),
                       ylabel='Frozen source events', ylim=(0, max(len(chosen)*1.25, 3.1)))
                ax.set_xticks(range(3), labels, rotation=15, ha='right', fontsize=8)
                if row_index == column == 0:
                    ax.legend(fontsize=8, ncol=3, loc='upper center')
                counts[dsid+'|'+cohort] = dict(native=dict(modes[0]), free_seed=dict(modes[1]), native_seed=dict(modes[2]))
        save(fig, 'fixed_seed_outcomes', 'Measurement changes still recover signs when all six seed means are restored',
             'Seed means, full covariance and surfaces are exact; constituent hit multisets match; native CKF may choose different paths; targeted counts are not a global misID rate')

        pure = [row for row in rows if row['cohort'] == 'pure_wrong']
        assert len(pure) == 6
        fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.5))
        for ax, row in zip(axes.flat, pure):
            charge = 1 if row['dsid'] == '100116' else -1
            qops = [charge*row[key] for key in ['native_final_qop_GeVinv', 'free_seed_centroid_qop_GeVinv', 'final_qop_GeVinv']]
            production = charge*prior[row['dsid'], row['sample_entry']]['replacement_qop_GeVinv']
            statuses = ['wrong', row['free_seed_centroid_outcome'], row['outcome']]
            for i, (value, status, color) in enumerate(zip(qops, statuses, ['#424242', '#1565c0', '#ef6c00'])):
                ax.scatter(i, value, color=color, marker='x' if status == 'quality_failed' else 'o', s=55)
            ax.axhline(0, color='black', linewidth=.8)
            ax.axhline(production, color='#00838f', linestyle=':', label='Production |q/p|')
            limit = 1.4*max([abs(value) for value in qops]+[production])
            ax.set(ylim=(-limit, limit), title=f"{row['dsid']} | source entry {row['sample_entry']}",
                   ylabel=r'$q_{truth}(q/p)_{final}$ [GeV$^{-1}$]')
            ax.set_xticks(range(3), ['Native', 'Centroid free', 'Centroid fixed'], rotation=15, ha='right')
            ax.legend(fontsize=8, loc='center', bbox_to_anchor=(.5, .63))
            if row['failed_quality_conditions']:
                ax.text(.04, .06, f"rmax = {row['quality']['rmax_mm']:.4f} mm > 95 mm", transform=ax.transAxes, fontsize=8)
            counts[row['dsid']+'|'+str(row['sample_entry'])] = dict(qop_truth_signed_GeVinv=qops,
                production_abs_qop_GeVinv=production, outcomes=statuses,
                same_native_and_fixed_selected_ordered_measurement_hits=row['same_selected_ordered_measurement_hits'],
                quality=row['quality'], failed_quality_conditions=row['failed_quality_conditions'])
        save(fig, 'fixed_seed_pure_curvature', 'Five selected recoveries and one radius-cut failure with the native seeds',
             'Positive means correct sign; x marks a quality failure; the production line is a reference, not local sensor truth; same selected ordered hits in 2 of 6 events')
    write_json(args.output/'bin_counts.json', counts)
    write_json(args.output/'figure_index.json', index)
    (args.output/'figure_index.md').write_text('\n'.join(
        '- ['+row['name']+']('+row['png']+') · [PDF]('+row['pdf']+')' for row in index)+'\n')
    print(json.dumps(dict(figures=len(index), output=str(args.output)), indent=2))


if __name__ == '__main__':
    main()
