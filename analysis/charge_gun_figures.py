#!/usr/bin/env python3
"""Plot independent full-Calypso gun controls from frozen measured summaries."""
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
CHARGES = [r'$\mu^+$', r'$\mu^-$']
MODES = ['beforeGeo', 'beforeGeo_regular_both', 'afterGeo']
LABELS = ['Before', 'Before + regular refit', 'After']
COLORS = ['#424242', '#ef6c00', '#1565c0']
MOMENTA = [100, 1000]


def rate(ax, x, row, **kwargs):
    value = 100*row['value']
    low, high = np.asarray(row['intervals']['95'])*100
    ax.errorbar(x, value, yerr=[[value-low], [high-value]], capsize=3, **kwargs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads((args.root/'gun_analysis_v2/summary.json').read_text())
    assert data['all_actual_generated_conditions_exact']
    assert data['all_versions_same_truth'] and data['all_actual_libraries_match']
    args.output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2})
    index, counts = [], {}
    scope = 'Full Calypso G4 + digitization | 20,000 events per charge and momentum | Isolated smoothing ON'

    def group(dsid, momentum, mode):
        return data['groups'][f'{dsid}|{momentum}|{mode}']

    with PdfPages(args.output/'charge_gun_report.pdf') as report:
        def save(fig, name, title, footer):
            fig.suptitle(title, y=.985, fontsize=14)
            fig.text(.5, .943, scope, ha='center', fontsize=9)
            fig.text(.5, .018, footer, ha='center', fontsize=8)
            fig.tight_layout(rect=(.02, .045, .98, .925))
            fig.savefig(args.output/(name+'.png'), dpi=150)
            fig.savefig(args.output/(name+'.pdf'))
            report.savefig(fig)
            plt.close(fig)
            index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
                reference='S06 Q3/Q5 independent full-G4 muon gun',
                definition=scope+'; '+footer, bin_counts='bin_counts.json'))

        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        for column, momentum in enumerate(MOMENTA):
            for dsid, charge, marker, offset in zip(DSIDS, CHARGES, ['o', 's'], [-.06, .06]):
                for row, metric in enumerate(['conditional_misid', 'truth_correct_efficiency']):
                    values = []
                    for i, mode in enumerate(MODES):
                        value = group(dsid, momentum, mode)[metric]
                        rate(axes[row, column], i+offset, value, marker=marker,
                             color=COLORS[i], linestyle='none', label=charge if i == 0 else None)
                        values.append(value)
                    counts[f'{dsid}|{momentum}|{metric}'] = values
            axes[0, column].set(title=f'{momentum} GeV', ylabel='Charge misidentification [%]')
            axes[1, column].set(ylabel='Correct / generated fiducial truth [%]')
            if momentum == 1000:
                axes[0, column].set_yscale('log')
            else:
                axes[0, column].set_ylim(bottom=0)
            for ax in axes[:, column]:
                ax.set_xticks(range(3), LABELS)
                ax.legend()
        save(fig, 'gun_rates_efficiency', 'Charge and efficiency in identical generated spectra',
             '95% exact binomial intervals, including zero-count upper limits; mode-specific selected samples; 1 TeV misID axis is logarithmic')

        fig, axes = plt.subplots(1, 2, figsize=(12, 5.7))
        for momentum, ax in zip(MOMENTA, axes):
            for i, mode in enumerate(MODES):
                row = data['charge_differences'][f'{momentum}|{mode}']
                value = row['value']*100
                low, high = np.asarray(row['interval95'])*100
                ax.errorbar(i, value, yerr=[[value-low], [high-value]], marker='o', capsize=4, color=COLORS[i])
            ax.axhline(0, color='black', lw=.8)
            ax.set(title=f'{momentum} GeV', ylabel=r'$f(\mu^-)-f(\mu^+)$ [percentage points]',
                   xticks=range(3), xticklabels=LABELS)
        counts['charge_differences'] = data['charge_differences']
        save(fig, 'gun_charge_difference', 'Residual charge differences with identical launch conditions',
             'Conservative 95% difference bounds from simultaneous 97.5% Clopper-Pearson intervals; actual G4 interactions can differ by charge')

        for observable, family, name, title, xlabel, scale, log in [
            ('signed_qop', 'histograms', 'gun_signed_qop', 'Signed final curvature in the real apparatus', r'Final q/p [GeV$^{-1}$]', 1., True),
            ('relative_p', 'histograms', 'gun_relative_momentum', 'Final momentum response to generated momentum', '(reco p - production truth p) / production truth p', 1., True),
            ('st1_sy_minus_launch_fine', 'angular_histograms', 'gun_angle_residual', 'Physical slope residual at the actual launch plane', r'$s_y^{reco}-s_y^{launch}$ [$\mu$rad]', 1e6, False),
            ('st1_sy_pull', 'angular_histograms', 'gun_angle_pull', 'Formal physical slope pulls at the launch plane', r'$(s_y^{reco}-s_y^{launch})/\sigma_{s_y}$', 1., True),
            ('qop_pull', 'histograms', 'gun_qop_pull', 'Formal final curvature pulls', r'$(\kappa^{reco}-\kappa^{production})/\sigma_\kappa$', 1., True)]:
            fig, axes = plt.subplots(2, 2, figsize=(12, 9))
            for row, (dsid, charge) in enumerate(zip(DSIDS, CHARGES)):
                for column, momentum in enumerate(MOMENTA):
                    ax = axes[row, column]
                    for mode, label, color in zip(MODES, LABELS, COLORS):
                        hist = group(dsid, momentum, mode)[family][observable]
                        ax.stairs(np.asarray(hist['counts'])/hist['raw_entries'], np.asarray(hist['edges'])*scale,
                                  color=color, label=label)
                        counts[f'{dsid}|{momentum}|{mode}|{observable}'] = hist
                    ax.set(title=charge+f' | {momentum} GeV', xlabel=xlabel, ylabel='Fraction of selected tracks / bin')
                    if log:
                        ax.set_yscale('log')
                    ax.legend()
            footer = 'Normalization retains flow and invalid entries; each mode uses its own selected tracks; formal fitted covariance is not externally calibrated'
            if family == 'angular_histograms':
                footer = 'Actual Calypso projection to z_F=0 mm and actual G4 launch truth at that plane; no upstream material in this slope reference; flows retained'
            save(fig, name, title, footer)

        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        for column, momentum in enumerate(MOMENTA):
            for dsid, charge, marker, offset in zip(DSIDS, CHARGES, ['o', 's'], [-.06, .06]):
                rows = [group(dsid, momentum, mode)['angular_histograms']['residual_summary'] for mode in MODES]
                counts[f'{dsid}|{momentum}|angular_summary'] = rows
                axes[0, column].plot(np.arange(3)+offset, [row['absolute_residual_68']*1e6 for row in rows],
                                     marker=marker, label=charge)
                axes[1, column].plot(np.arange(3)+offset,
                    [100*row['absolute_residual_above_1mrad']/row['finite'] for row in rows], marker=marker, label=charge)
            axes[0, column].set(title=f'{momentum} GeV', ylabel=r'68% quantile of $|\Delta s_y|$ [$\mu$rad]')
            axes[1, column].set(ylabel=r'$|\Delta s_y|>1$ mrad / finite selected [%]')
            for ax in axes[:, column]:
                ax.set_xticks(range(3), LABELS)
                ax.legend()
        save(fig, 'gun_angle_core_tail', 'Slope core and large tails both respond to the refit chart',
             'Unbinned descriptive quantiles and raw tail fractions; separate mode-specific selected samples; no uncertainty is assigned to the quantiles')

        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        for column, momentum in enumerate(MOMENTA):
            for dsid, charge, marker, offset in zip(DSIDS, CHARGES, ['o', 's'], [-.06, .06]):
                for row, label in enumerate(['all_final', 'selected']):
                    values = [group(dsid, momentum, mode)['refit_counts'] for mode in MODES]
                    axes[row, column].plot(np.arange(3)+offset,
                        [100*value[label+'_fallback']/value[label+'_total'] for value in values], marker=marker, label=charge)
                counts[f'{dsid}|{momentum}|refit_counts'] = values
            axes[0, column].set(title=f'{momentum} GeV', ylabel='Fallback / all final tracks [%]')
            axes[1, column].set(ylabel='Fallback / selected tracks [%]')
            for ax in axes[:, column]:
                ax.set_xticks(range(3), LABELS)
                ax.legend()
        save(fig, 'gun_refit_fallback', 'Refit fallback depends strongly on the output denominator',
             'Actual refit return status; unselected fallback tracks remain in the all-final denominator; failed G4/reconstruction jobs are checked separately')

        fig, axes = plt.subplots(1, 2, figsize=(12, 5.7))
        for momentum, ax in zip(MOMENTA, axes):
            for dsid, charge, marker in zip(DSIDS, CHARGES, ['o', 's']):
                values = [1000*group(dsid, momentum, mode)['reconstruction_wall_seconds_per_source_event'] for mode in MODES]
                ax.plot(range(3), values, marker=marker, label=charge)
                counts[f'{dsid}|{momentum}|observed_wall_ms'] = values
            ax.set(title=f'{momentum} GeV', ylabel='Observed reconstruction wall / source event [ms]',
                   xticks=range(3), xticklabels=LABELS)
            ax.legend()
        save(fig, 'gun_observed_wall', 'Observed worker reconstruction time',
             'Different HTCondor workers and startup costs are included; these measurements do not establish a controlled CPU speed ratio')

        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        regions = ['xminus_yminus', 'xminus_yplus', 'xplus_yminus', 'xplus_yplus']
        for row, (dsid, charge) in enumerate(zip(DSIDS, CHARGES)):
            for column, momentum in enumerate(MOMENTA):
                ax = axes[row, column]
                for offset, mode, label, color in zip([-.12, 0, .12], MODES, LABELS, COLORS):
                    rows = [group(dsid, momentum, mode)['production_regions'][region]['conditional_misid'] for region in regions]
                    for i, value in enumerate(rows):
                        rate(ax, i+offset, value, marker='o', color=color, linestyle='none', label=label if i == 0 else None)
                    counts[f'{dsid}|{momentum}|{mode}|production_regions'] = rows
                ax.set(title=charge+f' | {momentum} GeV', ylabel='Charge misidentification [%]',
                       xticks=range(4), xticklabels=['x-, y-', 'x-, y+', 'x+, y-', 'x+, y+'])
                if momentum == 1000:
                    ax.set_yscale('symlog', linthresh=.02)
                ax.legend()
        save(fig, 'gun_spatial_regions', 'Launch quadrants in the actual detector and field',
             '95% exact binomial intervals; no exact reflection symmetry is assumed for the real apparatus; 1 TeV axes use symlog with 0.02% threshold')

    counts['migrations'] = data['migrations']
    counts['nominal_launch_precision'] = data['nominal_launch_precision']
    write_json(args.output/'bin_counts.json', counts)
    write_json(args.output/'figure_index.json', index)
    (args.output/'figure_index.md').write_text('\n'.join(
        '- ['+row['name']+']('+row['png']+') · [PDF]('+row['pdf']+')' for row in index)+'\n')
    print(json.dumps(dict(figures=len(index), output=str(args.output)), indent=2))


if __name__ == '__main__':
    main()
