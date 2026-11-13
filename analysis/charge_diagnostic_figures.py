#!/usr/bin/env python3
"""English charge-lineage and accepted overlap-balance report figures."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_tables import write_json


DSIDS = ['100116', '100117']
VARIANTS = ['beforeGeo', 'afterGeo']
COLORS = ['#c62828', '#1565c0']


def weighted_histogram(values, weights, edges):
    valid = np.isfinite(values)
    inside = valid & (values >= edges[0]) & (values < edges[-1])
    counts = np.histogram(values[inside], edges, weights=weights[inside])[0]
    sumw2 = np.histogram(values[inside], edges, weights=weights[inside]**2)[0]
    return dict(edges=edges.tolist(), counts=counts.tolist(), sumw=counts.tolist(), sumw2=sumw2.tolist(),
        raw_entries=len(values), total_sumw=float(weights.sum()), total_sumw2=float(weights@weights),
        underflow_sumw=float(weights[valid & (values < edges[0])].sum()),
        overflow_sumw=float(weights[valid & (values >= edges[-1])].sum()), invalid=int((~valid).sum()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    lineage = json.loads((args.root/'lineage_summary_v1/summary.json').read_text())
    overlap = json.loads((args.root/'overlap_v1/summary.json').read_text())
    assert overlap['balance']['held_out']['passes_frozen_balance_threshold']
    groups = lineage['groups']
    index, bins = [], {}
    scope = 'Actual Calypso | 500,000 source events per charge/configuration | Isolated smoothing ON'
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2, 'legend.fontsize': 8})
    with PdfPages(args.output/'charge_diagnostics_report.pdf') as pdf:
        def save(fig, name, title, reference, footer):
            fig.suptitle(title+'\n'+scope, fontsize=12)
            fig.text(.5, .012, footer, ha='center', fontsize=8)
            fig.tight_layout(rect=(0, .04, 1, .91))
            for suffix in ['png', 'pdf']:
                fig.savefig(args.output/(name+'.'+suffix), dpi=180)
            pdf.savefig(fig)
            plt.close(fig)
            index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf', reference=reference,
                              definition=footer, bin_counts='bin_counts.json'))

        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
        origins = ['seed', 'CKF_measurement_update', 'refit']
        labels = ['Seed', 'CKF measurement update', 'Refit']
        for dsid, ax in zip(DSIDS, axes):
            x = np.arange(2)
            bottom = np.zeros(2)
            for origin, label, color in zip(origins, labels, ['#757575', '#ef6c00', '#6a1b9a']):
                n = np.array([groups[dsid+'|'+v]['final_wrong_entry'][origin] for v in VARIANTS])
                bars = ax.bar(x, n, bottom=bottom, label=label, color=color)
                for bar, count, base in zip(bars, n, bottom):
                    if count and origin != 'seed':
                        ax.text(bar.get_x()+bar.get_width()/2, base+count/2, str(count),
                                ha='center', va='center', color='white', fontsize=10)
                bottom += n
                bins[dsid+'|final_wrong_entry|'+origin] = n.tolist()
            for j, variant in enumerate(VARIANTS):
                assert bottom[j] == groups[dsid+'|'+variant]['counts']['wrong']
                seeds = groups[dsid+'|'+variant]['final_wrong_entry']['seed']
                ax.annotate(f'Total {int(bottom[j])}; seed {seeds}', (j, bottom[j]), xytext=(0, 9),
                            textcoords='offset points', ha='center')
            ax.set(title=dsid, ylabel='Final wrong-sign tracks', xticks=x, xticklabels=VARIANTS,
                   ylim=(0, max(bottom)*1.2))
            ax.legend(loc='upper right')
        save(fig, 'final_wrong_sign_origin', 'Last recorded transition into the final wrong sign',
            'S06 Q4; lineage classification',
            'Matched fiducial Good Tracks, purity > 0.7; native cached reference drives conversion; no conversion-induced sign flips')

        for observable, xlabel in [('curvature_significance', r'$q_{truth}\,\kappa/\sigma_\kappa$'),
                                   ('wrong_curvature_significance', r'Wrong-sign $q_{truth}\,\kappa/\sigma_\kappa$'),
                                   ('production_qop_pull', r'$(\kappa-\kappa_{prod})/\sigma_\kappa$')]:
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
            for dsid, ax in zip(DSIDS, axes):
                for variant, color in zip(VARIANTS, COLORS):
                    h = groups[dsid+'|'+variant]['histograms'][observable]
                    ax.stairs(h['counts'], h['edges'], color=color,
                        label=f"{variant}: N={h['raw_entries']:,}, flows={h['underflow']+h['overflow']:,}")
                    bins[dsid+'|'+variant+'|'+observable] = h
                ax.set(title=dsid, xlabel=xlabel, ylabel='Entries / bin', yscale='log')
                ax.legend()
            save(fig, observable, observable.replace('_', ' '), 'S06 curvature and uncertainty diagnostic',
                 'Final upstream fitted covariance; production truth; finite, invalid and flow counts retained')

        fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
        balance = overlap['balance']['held_out']
        x = np.arange(5)
        axes[0].plot(x, balance['original_SMD'], 'o-', label='Original')
        axes[0].plot(x, balance['overlap_SMD'], 's-', label='Overlap weighted')
        for y in [-.05, .05]:
            axes[0].axhline(y, color='black', ls='--', lw=.8)
        axes[0].set(ylabel='Standardized mean difference', title='Frozen criterion: all |SMD| < 0.05')
        axes[1].plot(x, balance['marginal_weighted_KS_distance'], 'o-', color='#1565c0')
        axes[1].axhline(.05, color='black', ls='--', lw=.8)
        axes[1].set(ylabel='Weighted marginal KS distance', title='Frozen criterion: all KS < 0.05', ylim=(0, .06))
        for ax in axes:
            ax.set(xticks=x, xticklabels=['log p', 'st1 x', 'st1 y', 'sensor sx', 'sensor sy'])
            ax.tick_params(axis='x', rotation=20)
        axes[0].legend()
        bins['accepted_overlap_balance'] = balance
        save(fig, 'accepted_overlap_balance', 'Entrance-condition overlap balance on held-out files 3/4',
             'S06 Q3; discovery-frozen truth-only logistic overlap model',
             'Entrance conditions include true magnetic bending and acceptance; this is distinct from generator-input balancing')

        with np.load(args.root/'overlap_v1/truth_weights.npz') as archive:
            weights = {dsid: archive[dsid] for dsid in DSIDS}
        values, masks = {}, {}
        for dsid in DSIDS:
            with np.load(args.root/f'main_analysis_v1/{dsid}_afterGeo_outcomes.npz') as archive:
                values[dsid] = archive['truth_features']
                masks[dsid] = archive['fiducial'] & (archive['file_index'] >= 3)
        edges = [np.linspace(np.log(5), np.log(5000), 61), np.linspace(-100, 100, 61),
                 np.linspace(-100, 100, 61), np.linspace(-.03, .03, 61), np.linspace(-.03, .03, 61)]
        fig, axes = plt.subplots(5, 2, figsize=(11, 15))
        for j, label in enumerate(['log production p [GeV]', 'st1 truth x [mm]', 'st1 truth y [mm]',
                                   'Sensor truth sx', 'Sensor truth sy']):
            for column, mode in enumerate(['Original', 'Overlap weighted']):
                ax = axes[j, column]
                for dsid, color in zip(DSIDS, COLORS):
                    mask = masks[dsid]
                    w = weights[dsid][mask] if column else np.ones(mask.sum())
                    h = weighted_histogram(values[dsid][mask, j], w, edges[j])
                    ax.stairs(np.array(h['counts'])/h['total_sumw'], h['edges'], color=color, label=dsid)
                    bins[f'heldout|{mode}|{dsid}|{j}'] = h
                ax.set(title=mode, xlabel=label, ylabel='Weighted fraction / bin')
                ax.legend()
        save(fig, 'accepted_overlap_spectra', 'Held-out entrance-condition spectra before and after overlap weighting',
             'S06 Q3; truth-only common support',
             'Weights learned on files 0/1/2; denominator includes flows; sumw, sumw2 and effective sample size retained')

        fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
        for variant, ax in zip(VARIANTS, axes):
            for dsid, color, shift in zip(DSIDS, COLORS, [-.07, .07]):
                rows = [overlap['rates'][split][variant][dsid]['conditional'] for split in ['all', 'held_out']]
                y = np.array([r['value'] for r in rows])*100
                ci = np.array([r['interval95'] for r in rows])*100
                ax.errorbar(np.arange(2)+shift, y, yerr=[y-ci[:, 0], ci[:, 1]-y], fmt='o', color=color, label=dsid)
                bins[variant+'|'+dsid+'|accepted_overlap_rates'] = rows
            ax.set(title=variant, ylabel='Conditional charge misidentification [%]', xticks=range(2),
                   xticklabels=['All files', 'Held-out files 3/4'])
            ax.set_ylim(bottom=0)
            ax.legend()
        save(fig, 'accepted_overlap_rates', 'Charge misidentification after entrance-condition overlap weighting | 95% intervals',
             'S06 Q3; accepted held-out balance',
             'Intervals condition on the discovery-frozen weight model; native before/after seed priors differ')
    write_json(args.output/'bin_counts.json', bins)
    write_json(args.output/'figure_index.json', index)
    (args.output/'figure_index.md').write_text('\n'.join(
        f"- [{r['name']}]({r['png']}) · [PDF]({r['pdf']}) · {r['reference']}" for r in index)+'\n')
    print(json.dumps(dict(figures=len(index), pdf=str(args.output/'charge_diagnostics_report.pdf')), indent=2))


if __name__ == '__main__':
    main()
