#!/usr/bin/env python3
"""Show actual native refit, primary-hit and q/B charge controls."""
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
MODES = ['native', 'smoother', 'updater', 'both']
COLORS = ['#c62828', '#1565c0']


def rates(ax, x, rows, **options):
    y = 100*np.asarray([row['value'] for row in rows], dtype=float)
    interval = 100*np.asarray([row['intervals']['95'] for row in rows], dtype=float)
    ax.errorbar(x, y, yerr=[y-interval[:, 0], interval[:, 1]-y], **options)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    refit = json.loads((args.root/'refit_analysis_v1/summary.json').read_text())
    hits = json.loads((args.root/'hitcontrol_analysis_v2/summary.json').read_text())
    filtered_events = {(row['dsid'], row['sample_entry']): row for row in
        json.loads((args.root/'hitcontrol_analysis_v2/events.json').read_text()) if row['mode'] == 'truth_filter'}
    events = json.loads((args.root/'replay_analysis_v1/events.json').read_text())
    assert refit['all_nonfinal_CKF_rows_exact'] and hits['all_238_jobs_succeeded']
    counts, index = {}, []
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2, 'legend.fontsize': 8})
    with PdfPages(args.output/'charge_mechanisms_report.pdf') as pdf:
        def save(fig, name, title, scope, footer):
            fig.suptitle(title+'\n'+scope, fontsize=12)
            fig.text(.5, .012, footer, ha='center', fontsize=8)
            fig.tight_layout(rect=(0, .05, 1, .91))
            for extension in ['png', 'pdf']:
                fig.savefig(args.output/(name+'.'+extension), dpi=180)
            pdf.savefig(fig)
            plt.close(fig)
            index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
                definition=scope+'; '+footer, bin_counts='bin_counts.json', reference='S06 Q4/Q5 actual Calypso mechanisms'))

        scope = 'beforeGeo | 100,000 source events per charge | Same native CKF records | Isolated smoothing ON'
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
        for dsid, color, offset in zip(DSIDS, COLORS, [-.04, .04]):
            groups = [refit['groups'][dsid+'|'+mode] for mode in MODES]
            for ax, observable in zip(axes, ['conditional_misid', 'truth_correct_efficiency']):
                rows = [group[observable] for group in groups]
                rates(ax, np.arange(4)+offset, rows, fmt='o-', capsize=3, color=color, label=dsid)
                counts[dsid+'|'+observable] = rows
        for ax, label in zip(axes, ['Charge misidentification [%]', 'Correct / fiducial truth [%]']):
            ax.set(xticks=range(4), xticklabels=MODES, ylabel=label)
            ax.legend()
        axes[0].set_ylim(bottom=0)
        save(fig, 'refit_chart_rates', 'Refit update controls the additional wrong signs', scope,
             '95% exact binomial intervals; each mode retains its own selected denominator; no change to CKF finding')

        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
        for dsid, ax in zip(DSIDS, axes):
            groups = [refit['refit_counts'][dsid+'|'+mode] for mode in MODES]
            for offset, key, color, label in [(-.18, 'selected_correct_to_wrong_in_refit', '#c62828', 'Correct to wrong'),
                                                (.18, 'selected_wrong_to_correct_in_refit', '#388e3c', 'Wrong to correct')]:
                values = [group[key] for group in groups]
                bars = ax.bar(np.arange(4)+offset, values, width=.35, color=color, label=label)
                ax.bar_label(bars, fontsize=9)
                counts[dsid+'|'+key] = values
            ax.set(title=dsid, xticks=range(4), xticklabels=MODES, ylabel='Selected final tracks')
            ax.legend()
        counts['refit_all_final_and_selected'] = refit['refit_counts']
        save(fig, 'refit_sign_migrations', 'Actual sign changes inside the native refit', scope,
             'Selected refits all succeed; all-output fallback counts use a separate denominator and are retained in bin_counts.json')

        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
        for dsid, ax in zip(DSIDS, axes):
            for mode, color in zip(MODES, ['#424242', '#ef6c00', '#1565c0', '#388e3c']):
                row = refit['groups'][dsid+'|'+mode]['histograms']['relative_p']
                ax.stairs(np.asarray(row['counts'])/max(1, row['raw_entries']), row['edges'], color=color, label=mode)
                counts[dsid+'|'+mode+'|relative_p'] = row
            ax.set(title=dsid, xlabel='(reco p - production truth p) / production truth p', ylabel='Fraction of selected tracks / bin', yscale='log')
            ax.legend()
        save(fig, 'refit_momentum_regression', 'Momentum response in the charge refit control', scope,
             'Production momentum reference; normalization includes underflow, overflow and invalid entries; mode-specific selected samples')

        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
        outcomes = ['correct', 'wrong', 'quality_failed', 'unreconstructed_or_unmatched', 'low_purity', 'invalid_curvature']
        colors = ['#388e3c', '#c62828', '#ef6c00', '#757575', '#8e24aa', '#212121']
        for cohort, ax in zip(['wrong', 'normal'], axes):
            bottom = np.zeros(2)
            for outcome, color in zip(outcomes, colors):
                values = [hits['filter_groups'][dsid+'|'+cohort].get(outcome, 0) for dsid in DSIDS]
                bars = ax.bar(np.arange(2), values, bottom=bottom, color=color, label=outcome.replace('_', ' '))
                for bar, value, base in zip(bars, values, bottom):
                    if value:
                        ax.text(bar.get_x()+bar.get_width()/2, base+value/2, str(value), ha='center', va='center', color='white')
                bottom += values
            ax.set(title='Originally '+cohort, xticks=[0, 1], xticklabels=DSIDS, ylabel='Frozen source events')
        axes[1].legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=8)
        counts['primary_hit_filter'] = hits
        save(fig, 'primary_hit_filter_outcomes', 'Primary-SDO source-link control in the native combinatorial CKF',
             'All 59 afterGeo wrong signs and 60 normal controls | Seed means and full covariance unchanged',
             '17 recovered correct, 6 remain wrong, 36 fail quality; losing a selected track is not counted as charge recovery')

        pure = [row for row in events if row['wrong'] and row['final_purity'] == 1.]
        assert len(pure) == 6
        for row in pure:
            control = filtered_events[row['dsid'], row['sample_entry']]
            assert control['outcome'] == 'wrong'
            assert control['final_qop_GeVinv'] == row['final_qop_GeVinv']
        fig, ax = plt.subplots(figsize=(12, 5.5))
        # PDG/charge is fixed by the DSID, independently of reconstructed sign.
        values = [(1 if row['dsid'] == '100116' else -1)*row['final_qop_GeVinv']/np.sqrt(row['final_variance_GeVinv2']) for row in pure]
        ax.bar(np.arange(6), values, color=[COLORS[DSIDS.index(row['dsid'])] for row in pure])
        ax.axhline(0, color='black', lw=.8)
        ax.set(xticks=range(6), xticklabels=[row['dsid']+'\n'+str(row['sample_entry']) for row in pure],
               ylabel=r'$q_{truth}\,\kappa / \sigma_\kappa$', xlabel='Sample and frozen sample entry')
        counts['pure_wrong_significance'] = [dict(dsid=row['dsid'], sample_entry=row['sample_entry'],
            qop_GeVinv=row['final_qop_GeVinv'], variance_GeVinv2=row['final_variance_GeVinv2'],
            truth_signed_significance=value) for row, value in zip(pure, values)]
        save(fig, 'pure_wrong_significance', 'The six pure wrong-sign tracks survive primary-hit filtering unchanged',
             'afterGeo | Actual final q/p and variance | All six CKF crossings use primary-supported hits',
             'These are formal fitted significances; covariance calibration and individual physical causes require separate controls')
    write_json(args.output/'bin_counts.json', counts)
    write_json(args.output/'figure_index.json', index)
    (args.output/'figure_index.md').write_text('\n'.join(
        '- ['+row['name']+']('+row['png']+') · [PDF]('+row['pdf']+')' for row in index)+'\n')
    print(json.dumps(dict(figures=len(index), output=str(args.output)), indent=2))


if __name__ == '__main__':
    main()
