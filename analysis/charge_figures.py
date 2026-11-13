#!/usr/bin/env python3
"""English Alma-style charge, momentum and spectrum-balance figures."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
import numpy as np

from charge_tables import truth_balance
from slope_tables import histogram, write_json


DSIDS = ['100116', '100117']
VARIANTS = ['beforeGeo', 'afterGeo']
COLORS = ['#c62828', '#1565c0']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--include-exploratory-balance', action='store_true', help='Append the joint-bin control that did not meet held-out balance criteria')
    args = parser.parse_args()
    result = json.loads((args.input/'summary.json').read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    groups = result['groups']
    index, bins = [], {}
    scope = r'Native Calypso | 500,000 source events per charge/configuration | Isolated smoothing ON'
    selection = 'Fiducial primary; charge-independent SDO match; purity > 0.7; Good Track r < 95 mm'
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2, 'legend.fontsize': 8})
    with PdfPages(args.output/'charge_reference_report.pdf') as pdf:
        def save(fig, name, title, reference, footer=selection):
            fig.suptitle(title+'\n'+scope, fontsize=12)
            fig.text(.5, .012, footer, ha='center', fontsize=8)
            fig.tight_layout(rect=(0, .04, 1, .91))
            for suffix in ['png', 'pdf']:
                fig.savefig(args.output/(name+'.'+suffix), dpi=180)
            pdf.savefig(fig)
            plt.close(fig)
            index.append(dict(name=name, reference=reference, png=name+'.png', pdf=name+'.pdf',
                              selection=footer, bin_counts='bin_counts.json'))

        for observable, xlabel, reference, log in [
            ('signed_qop', r'$q/p$ [GeV$^{-1}$]', 'Alma pp16,23-25', True),
            ('wrong_signed_qop', r'Wrong-sign $q/p$ [GeV$^{-1}$]', 'Alma pp16-18,23-25', False),
            ('qop_residual', r'$(q/p-q/p_{prod})/|q/p_{prod}|$', 'Alma pp23-25; truth residual', True),
            ('qop_pull', r'$(q/p-q/p_{prod})/\sigma_{q/p}$', 'Alma pp23-25; covariance diagnostic', True),
            ('truth_p', 'Production truth momentum [GeV]', 'Alma p21', True),
            ('reco_p', 'Reconstructed momentum [GeV]', 'Alma pp22,27', True),
            ('relative_p', r'$(p_{reco}-p_{prod})/p_{prod}$', 'Alma pp19-22; momentum tail', True),
        ]:
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
            for dsid, ax in zip(DSIDS, axes):
                for variant, color in zip(VARIANTS, COLORS):
                    g = groups[dsid+'|'+variant]['r95']
                    h = g['histograms'][observable]
                    bins[dsid+'|'+variant+'|'+observable] = h
                    ax.stairs(h['counts'], h['edges'], color=color,
                        label=f"{variant}: N={h['raw_entries']:,}, flow={h['underflow']+h['overflow']:,}")
                ax.set(title=r'$\mu^+$' if dsid == '100116' else r'$\mu^-$', xlabel=xlabel,
                       ylabel='Entries / bin', yscale='log' if log else 'linear')
                ax.legend()
            save(fig, observable, observable.replace('_', ' '), reference,
                 selection+('; truth reference is production, not local sensor p' if 'qop_' in observable else ''))

        fig, axes = plt.subplots(2, 2, figsize=(10, 9), sharex=True, sharey=True)
        for i, dsid in enumerate(DSIDS):
            for j, variant in enumerate(VARIANTS):
                h = groups[dsid+'|'+variant]['r95']['histograms']['truth_reco_p']
                ax = axes[i, j]
                ax.pcolormesh(h['truth_edges'], h['reco_edges'], np.asarray(h['counts']).T,
                              norm=LogNorm(vmin=1), rasterized=True)
                ax.plot([0, 2000], [0, 2000], color='red', lw=.7)
                ax.set(title=dsid+' | '+variant+f" | outside={h['outside']:,}",
                       xlabel='Production truth momentum [GeV]', ylabel='Reconstructed momentum [GeV]')
                bins[dsid+'|'+variant+'|truth_reco_p'] = h
        save(fig, 'truth_reco_momentum', 'Truth versus reconstructed momentum', 'Alma pp19-20')

        for radius in [95, 100]:
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
            for dsid, ax in zip(DSIDS, axes):
                for variant, color in zip(VARIANTS, COLORS):
                    rows = groups[dsid+'|'+variant][f'r{radius}']['truth_p_bins']
                    values = [r['conditional_misid'] for r in rows]
                    x = np.arange(len(values))
                    y = np.array([v['value']*100 if v['total'] else np.nan for v in values])
                    ci = np.array([v['intervals']['68'] if v['total'] else [np.nan]*2 for v in values])*100
                    ax.errorbar(x, y, yerr=np.stack([y-ci[:, 0], ci[:, 1]-y]), fmt='o-',
                                label=variant, color=color)
                    bins[dsid+'|'+variant+f'|misid_p_r{radius}'] = rows
                ax.set(title=dsid, xlabel='Truth momentum bin [GeV]', ylabel='Charge misidentification [%]',
                    xticks=range(6), xticklabels=['0-50', '50-100', '100-250', '250-500', '500-1000', '1000+'])
                ax.tick_params(axis='x', rotation=25)
                ax.legend()
            save(fig, f'misid_truth_p_r{radius}', f'Charge misidentification versus truth p | r < {radius} mm | 68% intervals',
                 'Alma pp17-18; truth-p conditioning', selection.replace('95 mm', f'{radius} mm'))

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
        for dsid, ax in zip(DSIDS, axes):
            for variant, color, shift in zip(VARIANTS, COLORS, [-.08, .08]):
                g = groups[dsid+'|'+variant]['r95']
                values = [g['conditional_misid'], g['exactly_one_good'], g['held_out']['conditional_misid']]
                y = np.array([v['value'] for v in values])*100
                ci = np.array([v['intervals']['95'] for v in values])*100
                ax.errorbar(np.arange(3)+shift, y, yerr=np.stack([y-ci[:, 0], ci[:, 1]-y]),
                            fmt='o', color=color, label=variant)
                for x, v in enumerate(values):
                    ax.annotate(f"{v['wrong']}/{v['total']:,}", (x+shift, v['value']*100),
                                xytext=(0, 12), textcoords='offset points', ha='center', fontsize=7)
            ax.set(title=dsid, ylabel='Charge misidentification [%]',
                   xticks=range(3), xticklabels=['Inclusive', 'Exactly one Good Track', 'Held-out files 3/4'])
            ax.set_ylim(0, max(groups[dsid+'|'+v]['r95']['held_out']['conditional_misid']['intervals']['95'][1]
                               for v in VARIANTS)*130)
            ax.margins(x=.12)
            ax.legend(loc='center right')
        save(fig, 'global_charge_rates', 'Charge misidentification: denominators and 95% intervals', 'Alma pp17-18')

        fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
        for dsid, ax in zip(DSIDS, axes):
            migration = result['migrations'][dsid]['r95']
            matrix = np.asarray(migration['counts'])
            labels = ['Unreconstructed', 'Correct', 'Wrong', 'Invalid', 'Low purity', 'Quality failed', 'Not fiducial']
            ax.imshow(matrix, norm=LogNorm(vmin=1, vmax=max(1, matrix.max())))
            for i in range(7):
                for j in range(7):
                    ax.text(j, i, f'{matrix[i,j]:,}', ha='center', va='center', fontsize=7,
                            color='white' if matrix[i,j] > 500 else 'black')
            ax.set(title=dsid, xlabel='afterGeo', ylabel='beforeGeo', xticks=range(7), yticks=range(7),
                   xticklabels=labels, yticklabels=labels)
            ax.tick_params(axis='x', rotation=60)
            bins[dsid+'|migration'] = migration
        save(fig, 'geometry_state_migration', 'Same-truth before/after state migration', 'S06 paired geometry comparison',
             'All source primary muons; includes reconstruction and quality failures; native priors differ')

        if args.include_exploratory_balance:
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
            for variant, ax in zip(VARIANTS, axes):
                for dsid, color, shift in zip(DSIDS, COLORS, [-.08, .08]):
                    row = result['weighted_measurements'][variant][dsid]['held_out']
                    y = np.array([row[k]['value'] for k in ['raw', 'common_support', 'reweighted']])*100
                    ci = np.array([row[k]['intervals']['95'] if k != 'reweighted' else row[k]['interval95']
                                   for k in ['raw', 'common_support', 'reweighted']])*100
                    ax.errorbar(np.arange(3)+shift, y, yerr=np.stack([y-ci[:, 0], ci[:, 1]-y]),
                                fmt='o-', color=color, label=dsid)
                ax.set(title=variant, ylabel='Charge misidentification [%]', xticks=range(3),
                       xticklabels=['Original truth spectrum', 'Common support', 'Reweighted'])
                ax.tick_params(axis='x', rotation=15)
                ax.legend()
            save(fig, 'held_out_spectrum_control', 'Held-out charge comparison after truth-spectrum balance | 95% intervals',
                 'S06 input-spectrum control', 'Weights learned only from discovery files 0/1/2; tested on held-out files 3/4')

            fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
            for variant, ax in zip(VARIANTS, axes):
                balance = result['spectrum_balance'][variant]['balance']['held_out']
                for key, label, style in [('before', 'Original', 'o-'), ('common_support', 'Common support', 's--'),
                                           ('after', 'Reweighted', '^:')]:
                    ax.plot(np.arange(5), balance['standardized_mean_difference_'+key], style, label=label)
                ax.axhline(0, color='black', lw=.7)
                ax.set(title=variant, ylabel='Standardized mean difference', xticks=range(5),
                       xticklabels=['log p', 'st1 x', 'st1 y', 'true sx', 'true sy'])
                ax.legend()
            save(fig, 'held_out_truth_balance', 'Truth-input balance on held-out files', 'S06 spectrum/acceptance control',
                 'Common support and weighting are shown separately; truth variables only')

            payload, derived = {}, {}
            for dsid in DSIDS:
                with np.load(args.input/f'{dsid}_afterGeo_outcomes.npz') as archive:
                    p = {k: archive[k] for k in ['truth_features', 'fiducial', 'file_index']}
                payload[dsid] = {'truth': {'file_index': p['file_index']}}
                derived[dsid] = dict(features=p['truth_features'], fiducial=p['fiducial'])
            weights, repeated = truth_balance(payload, derived)
            assert repeated['balance'] == result['spectrum_balance']['afterGeo']['balance']
            edges = [np.linspace(np.log(5), np.log(5000), 61), np.linspace(-100, 100, 61),
                     np.linspace(-100, 100, 61), np.linspace(-.03, .03, 61), np.linspace(-.03, .03, 61)]
            fig, axes = plt.subplots(5, 2, figsize=(11, 16))
            for j, label in enumerate(['log truth p [GeV]', 'st1 truth x [mm]', 'st1 truth y [mm]', 'True sx', 'True sy']):
                for col, mode in enumerate(['Original', 'Reweighted']):
                    ax = axes[j, col]
                    for dsid, color in zip(DSIDS, COLORS):
                        eligible = derived[dsid]['fiducial'] & (payload[dsid]['truth']['file_index'] >= 3)
                        values = derived[dsid]['features'][eligible, j]
                        w = weights[dsid][eligible] if col else np.ones(len(values))
                        valid = np.isfinite(values)
                        in_range = valid & (values < edges[j][-1])
                        counts = np.histogram(values[in_range], edges[j], weights=w[in_range])[0]
                        sw2 = np.histogram(values[in_range], edges[j], weights=w[in_range]**2)[0]
                        total = w.sum()
                        ax.stairs(counts/total, edges[j], color=color, label=dsid)
                        bins[f'heldout_truth_{j}_{mode}_{dsid}'] = dict(edges=edges[j].tolist(),
                            counts=counts.tolist(), sumw2=sw2.tolist(), total_sumw=float(total),
                            below_sumw=float(w[valid & (values < edges[j][0])].sum()),
                            above_sumw=float(w[valid & (values >= edges[j][-1])].sum()), invalid=int((~valid).sum()))
                    ax.set(title=mode, xlabel=label, ylabel='Fraction / bin')
                    ax.legend()
            save(fig, 'held_out_truth_spectra', 'Original and balanced truth input spectra | afterGeo',
                 'S06 truth p/position/slope spectra', 'Held-out fiducial primary truth; joint-bin weights frozen on discovery; denominator includes flows')
    write_json(args.output/'bin_counts.json', bins)
    write_json(args.output/'figure_index.json', index)
    (args.output/'figure_index.md').write_text('\n'.join(
        f"- [{r['name']}]({r['png']}) · [PDF]({r['pdf']}) · {r['reference']}" for r in index)+'\n')


if __name__ == '__main__':
    main()
