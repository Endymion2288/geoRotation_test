#!/usr/bin/env python3
"""Plot measured response to a first-order, physically matched seed covariance."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from charge_control_figures import intervals
from slope_tables import write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    main_result = json.loads((args.root/'main_analysis_v1/summary.json').read_text())
    prior = json.loads((args.root/'prior_analysis_v1/summary.json').read_text())
    labels = ['Before: native prior', 'After: native prior', 'After: matched prior']
    dsids, colors = ['100116', '100117'], ['#c62828', '#1565c0']
    bins, index = {}, []
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2})
    with PdfPages(args.output/'charge_prior_report.pdf') as report:
        def save(fig, name, title):
            fig.suptitle(title+'\n500,000 RDO events per charge | Actual Calypso | Isolated smoothing ON', fontsize=12)
            fig.text(.5, .01, 'Matched prior: full first-order J C_F J^T at native means; finite-width nonlinear prior equality is not assumed',
                     ha='center', fontsize=8)
            fig.tight_layout(rect=(0, .05, 1, .91))
            for suffix in ['png', 'pdf']:
                fig.savefig(args.output/(name+'.'+suffix), dpi=180)
            report.savefig(fig)
            plt.close(fig)
            index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
                reference='S06 Q5: actual matched seed covariance control',
                definition=title+'; 500k per charge; first-order covariance match only',
                bin_counts='bin_counts.json'))

        groups = {dsid: [main_result['groups'][dsid+'|'+v]['r95'] for v in ['beforeGeo', 'afterGeo']]
                  + [prior['groups'][dsid]] for dsid in dsids}
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
        for ax, key, ylabel in zip(axes, ['conditional_misid', 'truth_correct_efficiency'],
                ['Wrong / selected valid matched [%]', 'Correct / fiducial truth [%]']):
            for dsid, offset, color in zip(dsids, [-.06, .06], colors):
                rows = [g[key] for g in groups[dsid]]
                intervals(ax, np.arange(3)+offset, rows, fmt='o', color=color, capsize=3, label=dsid)
                bins[dsid+'|'+key] = rows
            ax.set(xticks=range(3), xticklabels=labels, ylabel=ylabel)
            ax.tick_params(axis='x', rotation=16)
            ax.legend()
        axes[0].set_ylim(bottom=0)
        save(fig, 'matched_prior_rates_efficiency', 'Native version comparison and matched-prior response | 95% intervals')

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        for i, dsid in enumerate(dsids):
            for ax, observable in zip(axes[i], ['conditional_misid', 'truth_correct_efficiency']):
                for g, label, color, offset in zip(groups[dsid], labels, ['#757575', '#1565c0', '#c62828'], [-.13, 0, .13]):
                    rows = [b[observable] for b in g['truth_p_bins']]
                    intervals(ax, np.arange(len(rows))+offset, rows, fmt='o', color=color, capsize=2, label=label)
                    bins[dsid+'|'+observable+'|'+label+'|truth_p_bins'] = g['truth_p_bins']
                xlabels = [f"{b['low_GeV']:g}-{b['high_GeV']:g}" for b in groups[dsid][0]['truth_p_bins']]
                ax.set(title=dsid, xticks=range(len(xlabels)), xticklabels=xlabels,
                       xlabel='Production truth momentum [GeV]',
                       ylabel='Misidentification [%]' if observable == 'conditional_misid' else 'Correct / fiducial [%]')
                ax.tick_params(axis='x', rotation=25)
                ax.legend(fontsize=7)
                if observable == 'conditional_misid':
                    ax.set_ylim(bottom=0)
        save(fig, 'matched_prior_truth_p', 'Momentum dependence and efficiency | Exact binomial 95% intervals')

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        for i, dsid in enumerate(dsids):
            for j, observable in enumerate(['signed_qop', 'relative_p']):
                ax = axes[i, j]
                for g, label, color in zip(groups[dsid], labels, ['#757575', '#1565c0', '#c62828']):
                    h = g['histograms'][observable]
                    ax.stairs(np.array(h['counts'])/h['raw_entries'], h['edges'], label=label, color=color)
                    bins[dsid+'|'+observable+'|'+label] = h
                ax.set(title=dsid, ylabel='Selected matched fraction / bin',
                       xlabel=r'Signed q/p [GeV$^{-1}$]' if observable == 'signed_qop'
                       else r'$(p_{reco}-p_{production})/p_{production}$')
                ax.set_yscale('log')
                ax.legend(fontsize=7)
        save(fig, 'matched_prior_momentum', 'Signed curvature and production-referenced momentum residual')
    write_json(args.output/'bin_counts.json', bins)
    write_json(args.output/'figure_index.json', index)
    (args.output/'figure_index.md').write_text('\n'.join(
        f"- [{r['name']}]({r['png']}) · [PDF]({r['pdf']})" for r in index)+'\n')
    print(json.dumps(dict(figures=len(index), output=str(args.output)), indent=2))


if __name__ == '__main__':
    main()
