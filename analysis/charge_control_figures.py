#!/usr/bin/env python3
"""Report production balance, native mode controls and actual charge crossings."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from charge_diagnostic_figures import weighted_histogram
from charge_modes import MODES
from slope_tables import histogram, write_json


DSIDS = ['100116', '100117']
VARIANTS = ['beforeGeo', 'afterGeo']
COLORS = ['#c62828', '#1565c0']


def intervals(ax, positions, rows, **options):
    y = 100*np.array([r['value'] for r in rows], dtype=float)
    ci = 100*np.array([r['interval95'] if 'interval95' in r else r['intervals']['95'] for r in rows], dtype=float)
    ax.errorbar(positions, y, yerr=[y-ci[:, 0], ci[:, 1]-y], **options)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    production = json.loads((args.root/'production_overlap_v1/summary.json').read_text())
    modes = json.loads((args.root/'mode_analysis_v1/summary.json').read_text())
    replay = json.loads((args.root/'replay_analysis_v1/events.json').read_text())
    main_result = json.loads((args.root/'main_analysis_v1/summary.json').read_text())
    assert production['balance']['held_out']['passes_frozen_balance_threshold']
    bins, index = {}, []
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2, 'legend.fontsize': 8})
    with PdfPages(args.output/'charge_controls_report.pdf') as pdf:
        def save(fig, name, title, scope, footer):
            fig.suptitle(title+'\n'+scope, fontsize=12)
            fig.text(.5, .012, footer, ha='center', fontsize=8)
            fig.tight_layout(rect=(0, .04, 1, .91))
            for suffix in ['png', 'pdf']:
                fig.savefig(args.output/(name+'.'+suffix), dpi=180)
            pdf.savefig(fig)
            plt.close(fig)
            index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
                reference='S06 Q3/Q4/Q5 controls', definition=scope+'; '+footer,
                bin_counts='bin_counts.json'))

        scope = 'Production-only truth | 500,000 primary muons per charge | No entrance-fiducial conditioning in training'
        balance = production['balance']['held_out']
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.3))
        x = np.arange(5)
        for name, label in [('original_SMD', 'Original'), ('overlap_SMD', 'Overlap weighted')]:
            axes[0].plot(x, balance[name], 'o-', label=label)
        for y in [-.05, .05]:
            axes[0].axhline(y, color='black', ls='--', lw=.8)
        axes[1].plot(x, balance['marginal_weighted_KS_distance'], 'o-')
        axes[1].axhline(.05, color='black', ls='--', lw=.8)
        for ax, ylabel in zip(axes, ['Standardized mean difference', 'Weighted marginal KS distance']):
            ax.set(xticks=x, xticklabels=['log p', 'x', 'y', 'sx', 'sy'], ylabel=ylabel)
        axes[0].legend()
        bins['production_balance'] = balance
        save(fig, 'production_balance', 'Production-spectrum balance on held-out files 3/4', scope,
             'Both frozen criteria pass: all |SMD| < 0.05 and all weighted marginal KS distances < 0.05')

        with np.load(args.root/'production_features_v1/features.npz') as archive:
            features = {dsid: archive[dsid] for dsid in DSIDS}
        with np.load(args.root/'production_overlap_v1/truth_weights.npz') as archive:
            weights = {dsid: archive[dsid] for dsid in DSIDS}
        masks = {}
        for dsid in DSIDS:
            with np.load(args.root/f'main_analysis_v1/{dsid}_afterGeo_outcomes.npz') as archive:
                masks[dsid] = archive['file_index'] >= 3
        edges = [np.linspace(np.log(1), np.log(10000), 61), np.linspace(-150, 150, 61),
                 np.linspace(-150, 150, 61), np.linspace(-.05, .05, 61), np.linspace(-.05, .05, 61)]
        fig, axes = plt.subplots(5, 2, figsize=(11, 15))
        for j, label in enumerate(['log production p [GeV]', 'Production x [mm]', 'Production y [mm]',
                                   'Production sx', 'Production sy']):
            for col, title in enumerate(['Original', 'Overlap weighted']):
                for dsid, color in zip(DSIDS, COLORS):
                    mask = masks[dsid]
                    w = weights[dsid][mask] if col else np.ones(mask.sum())
                    h = weighted_histogram(features[dsid][mask, j], w, edges[j])
                    axes[j, col].stairs(np.array(h['counts'])/h['total_sumw'], h['edges'], label=dsid, color=color)
                    bins[f'production|{dsid}|{title}|{j}'] = h
                axes[j, col].set(title=title, xlabel=label, ylabel='Weighted fraction / bin')
                axes[j, col].legend()
        save(fig, 'production_spectra', 'Held-out production spectra before and after overlap weighting', scope,
             'Files 0/1/2 train weights; files 3/4 validate; normalization includes flows; sumw/sumw2 retained')

        fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
        for variant, ax in zip(VARIANTS, axes):
            for offset, kind, title, color in [(-.06, 'raw', 'Raw', COLORS[0]),
                                                (.06, 'weighted', 'Production overlap', COLORS[1])]:
                rows = [main_result['charge_differences'][variant][split]['raw'] if kind == 'raw'
                        else production['differences'][split][variant] for split in ['all', 'held_out']]
                y = 100*np.array([r['mu_minus_minus_mu_plus'] for r in rows])
                ci = 100*np.array([r['interval95'] for r in rows])
                ax.errorbar(np.arange(2)+offset, y, yerr=[y-ci[:, 0], ci[:, 1]-y], fmt='o', label=title, color=color)
                bins[variant+'|production_difference|'+kind] = rows
            ax.axhline(0, color='black', lw=.8)
            ax.set(title=variant, xticks=[0, 1], xticklabels=['All files', 'Held-out files 3/4'],
                   ylabel=r'$f_{mis}(\mu^-)-f_{mis}(\mu^+)$ [percentage points]')
            ax.legend()
        save(fig, 'production_charge_difference', 'Charge-rate difference and 95% intervals', scope,
             'Reconstruction uses isolated smoothing ON; weighted uncertainty conditions on discovery-frozen model')

        fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
        for variant, ax in zip(VARIANTS, axes):
            for dsid, offset, color in zip(DSIDS, [-.06, .06], COLORS):
                rows = [production['rates']['held_out'][variant][dsid][key] for key in
                        ['fiducial_acceptance', 'truth_correct_efficiency', 'full_population_correct_efficiency']]
                intervals(ax, np.arange(3)+offset, rows, fmt='o', color=color, label=dsid)
                bins[variant+'|production_efficiency|'+dsid] = rows
            ax.set(title=variant, xticks=range(3), xticklabels=['Truth fiducial / all', 'Correct / fiducial', 'Correct / all'],
                   ylabel='Overlap-weighted probability [%]')
            ax.legend()
        save(fig, 'production_acceptance_efficiency', 'Separate truth acceptance and reconstruction efficiency', scope,
             'Held-out files only; conditioning on entrance kinematics is reported separately')

        mode_scope = 'Actual Calypso | 100,000 source events per charge/version/mode | Isolated smoothing ON'
        for observable, ylabel, name in [('conditional_misid', 'Charge misidentification [%]', 'mode_misid'),
                                         ('truth_correct_efficiency', 'Correct matched / fiducial truth [%]', 'mode_efficiency')]:
            fig, axes = plt.subplots(2, 2, figsize=(12, 8))
            for i, dsid in enumerate(DSIDS):
                for j, variant in enumerate(VARIANTS):
                    ax = axes[i, j]
                    rows = [modes['groups']['|'.join([dsid, variant, mode])][observable] for mode in MODES]
                    intervals(ax, np.arange(len(MODES)), rows, fmt='o', color=COLORS[j], capsize=3)
                    ax.set(title=dsid+' | '+variant, ylabel=ylabel, xticks=range(len(MODES)), xticklabels=MODES)
                    ax.tick_params(axis='x', rotation=25)
                    if observable == 'conditional_misid':
                        ax.set_ylim(bottom=0)
                    bins[dsid+'|'+variant+'|'+name] = rows
            save(fig, name, 'Native direction and station controls | 95% intervals', mode_scope,
                 'First 20,000 entries of each of five files; same charge-independent primary match and downstream Good Track cuts')

        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        for i, dsid in enumerate(DSIDS):
            for j, variant in enumerate(VARIANTS):
                ax = axes[i, j]
                bottom = np.zeros(2)
                for origin, color in zip(['forward', 'backward', 'both'], ['#388e3c', '#ef6c00', '#757575']):
                    values = np.array([modes['merge_origins']['|'.join([dsid, variant, station+'merge'])]
                                       ['counts']['selected'].get(origin, 0) for station in ['3', '4']])
                    ax.bar(np.arange(2), values, bottom=bottom, label=origin, color=color)
                    bottom += values
                    bins[f'{dsid}|{variant}|merge_origin|{origin}'] = values.tolist()
                ax.set(title=dsid+' | '+variant, xticks=[0, 1], xticklabels=['3-station merge', '4-station merge'],
                       ylabel='Selected matched tracks')
                ax.legend()
        save(fig, 'native_merge_origins', 'Actual merged-output origins', mode_scope,
             'Every output closes against the full native source hit set, q/p and upstream parameters; no approximate matching')

        wrong = [r for r in replay if r['wrong']]
        normal = [r for r in replay if not r['wrong']]
        crossings = [r['last_update_into_wrong'] for r in wrong if r['last_wrong_entry'] == 3]
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.3))
        for category, color in [('no_SDO', '#ef6c00'), ('primary_only', '#1565c0')]:
            rows = [r for r in crossings if r['hit_truth'] == category]
            axes[0].scatter([r['residual_mm'][0] for r in rows],
                            [r['gain_qop_GeVinv_per_mm'][0] for r in rows],
                            s=25, alpha=.8, color=color, label=f'{category}: N={len(rows)}')
        axes[0].set(xlabel='Innovation residual [mm]', ylabel=r'$K_\kappa$ [GeV$^{-1}$/mm]', title='Last CKF crossing into final wrong sign')
        axes[0].legend()
        for values, label, color in [([r['innovation_chi2'] for r in crossings], 'Last wrong-sign crossings (56)', '#ef6c00'),
                                     ([u['innovation_chi2'] for r in normal for u in r['per_update']], 'Normal-control updates', '#1565c0')]:
            h = histogram(values, np.linspace(0, 25, 51))
            axes[1].stairs(np.array(h['counts'])/h['raw_entries'], h['edges'], label=label, color=color)
            bins['crossing_chi2|'+label] = h
        axes[1].set(xlabel=r'Innovation $\chi^2$', ylabel='Fraction / bin')
        axes[1].legend()
        save(fig, 'charge_crossing_innovations', 'Measured Kalman updates that change the final sign',
             'All 59 afterGeo wrong signs + 60 frozen normal controls | Exact 58-field replay',
             'Three wrong signs persist from seed; no SDO association is not by itself proof of an electronic-noise hit')
        bins['crossing_operations'] = crossings

        pure = [r for r in wrong if r['final_purity'] == 1.]
        assert len(pure) == 6
        fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        for row, ax in zip(pure, axes.flat):
            charge = 1 if row['dsid'] == '100116' else -1
            updates = row['per_update']
            z = [u['z_F_mm'] for u in updates]
            for key, label, color in [('predicted_qop_GeVinv', 'Predicted', '#388e3c'),
                                       ('filtered_qop_GeVinv', 'Filtered', '#c62828')]:
                ax.plot(z, [charge*u[key] for u in updates], 'o-', ms=3, label=label, color=color)
            ax.axhline(0, color='black', lw=.8)
            ax.axhline(1./row['final_truth_p_GeV'], color='#1565c0', ls='--', label='Production q/p reference')
            ax.set(title=f"{row['dsid']} | source entry {row['source_entry']} | file {row['file_index']}",
                   xlabel='FASER z [mm]', ylabel=r'$q_{truth}\,\kappa$ [GeV$^{-1}$]')
            ax.legend(fontsize=7)
        save(fig, 'pure_wrong_sign_paths', 'All six wrong-sign final tracks with purity = 1',
             'Actual native afterGeo CKF | Same selected candidate and full hit sequence as 500k baseline',
             'Production q/p is a reference, not sensor-local momentum; final purity alone does not certify each CKF state')
        bins['pure_wrong_sign_paths'] = pure
    write_json(args.output/'bin_counts.json', bins)
    write_json(args.output/'figure_index.json', index)
    (args.output/'figure_index.md').write_text('\n'.join(
        f"- [{r['name']}]({r['png']}) · [PDF]({r['pdf']})" for r in index)+'\n')
    print(json.dumps(dict(figures=len(index), output=str(args.output)), indent=2))


if __name__ == '__main__':
    main()
