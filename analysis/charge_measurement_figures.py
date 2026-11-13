#!/usr/bin/env python3
"""Compare actual Calypso measurement and seed controls on frozen events."""
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


MODES = ['native', 'truth_filter', 'centroid', 'centroid_primary_filter', 'production_qop']
LABELS = ['Native', 'Primary hits', 'SiHit centroid', 'Centroid + primary hits', 'Production seed q/p']
COLORS = ['#424242', '#ef6c00', '#1565c0', '#388e3c', '#8e24aa']
OUTCOMES = ['correct', 'wrong', 'quality_failed', 'low_purity']


def load_controls(root):
    targets = json.loads((root/'charge_replay_v1/targets.json').read_text())
    records = {mode: {} for mode in MODES}
    identity = lambda row: (row['dsid'], row['sample_entry'])
    for row in targets:
        records['native'][identity(row)] = dict(row,
            original_wrong=row['wrong'], outcome='wrong' if row['wrong'] else 'correct',
            final_qop_GeVinv=row['final_qop_GeVinv'])
    for filename, mode_filter in [
        ('hitcontrol_analysis_v2/events.json', {'truth_filter'}),
        ('centroid_analysis_v1/events.json', {'centroid', 'centroid_primary_filter'}),
        ('seed_qop_analysis_v1/events.json', {'production_qop'})]:
        for row in json.loads((root/filename).read_text()):
            mode = row.get('mode', 'production_qop')
            if mode in mode_filter:
                key = identity(row)
                assert key not in records[mode]
                records[mode][key] = row
    assert len(records['native']) == 119
    for mode in MODES:
        assert records[mode].keys() == records['native'].keys()
        for key, row in records[mode].items():
            assert row['original_wrong'] == records['native'][key]['original_wrong']
            assert row['outcome'] in OUTCOMES
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    records = load_controls(args.root)
    args.output.mkdir(parents=True, exist_ok=False)
    counts, index = {}, []
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2})
    with PdfPages(args.output/'charge_measurement_report.pdf') as pdf:
        def save(fig, name, title, footer):
            fig.suptitle(title, fontsize=14, y=.99)
            fig.text(.5, .925, 'afterGeo | 59 native wrong signs + 60 frozen normal controls | Isolated smoothing ON',
                     ha='center', fontsize=9)
            fig.text(.5, .018, footer, ha='center', fontsize=8)
            fig.tight_layout(rect=(.02, .055, .98, .89))
            fig.savefig(args.output/(name+'.png'), dpi=150)
            fig.savefig(args.output/(name+'.pdf'))
            pdf.savefig(fig)
            plt.close(fig)
            index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
                reference='S06 Q4 measurement and seed single-factor controls',
                definition=footer, bin_counts='bin_counts.json'))

        fig, axes = plt.subplots(2, 2, figsize=(13, 9))
        for column, (dsid, charge) in enumerate([('100116', r'$\mu^+$'), ('100117', r'$\mu^-$')]):
            for row_index, wrong in enumerate([True, False]):
                ax = axes[row_index, column]
                selected = [key for key, row in records['native'].items()
                            if key[0] == dsid and row['original_wrong'] == wrong]
                values = [Counter(records[mode][key]['outcome'] for key in selected) for mode in MODES]
                bottom = np.zeros(len(MODES))
                for outcome, color in zip(OUTCOMES, ['#43a047', '#e53935', '#bdbdbd', '#fdd835']):
                    heights = np.array([value[outcome] for value in values])
                    ax.bar(range(len(MODES)), heights, bottom=bottom, color=color, label=outcome.replace('_', ' '))
                    for i, n in enumerate(heights):
                        if n:
                            ax.text(i, bottom[i]+n/2, str(n), ha='center', va='center', fontsize=9)
                    bottom += heights
                ax.set(title=charge+f" | Native {'wrong' if wrong else 'correct'} cohort, N={len(selected)}",
                       ylabel='Frozen source events', ylim=(0, len(selected)*1.17))
                ax.set_xticks(range(len(MODES)), LABELS, rotation=20, ha='right', fontsize=8)
                if row_index == 0 and column == 0:
                    ax.legend(ncol=2, fontsize=8, loc='upper center')
                counts[dsid+('|wrong' if wrong else '|normal')] = {
                    mode: dict(value) for mode, value in zip(MODES, values)}
        save(fig, 'measurement_control_outcomes', 'Measurement changes and seed curvature affect different wrong-sign events',
             'Targeted cohorts, not a global misID measurement; quality/purity failures are separate from selected correct signs; actual Calypso controls')

        pure = sorted(key for key, row in records['native'].items()
                      if row['original_wrong'] and row['purity'] == 1.)
        assert len(pure) == 6
        fig, axes = plt.subplots(2, 3, figsize=(14, 8.5))
        for ax, key in zip(axes.flat, pure):
            charge = 1. if key[0] == '100116' else -1.
            case = []
            for i, (mode, color) in enumerate(zip(MODES, COLORS)):
                row = records[mode][key]
                qop = row['final_qop_GeVinv']
                y = charge*qop if qop is not None else np.nan
                failed = row['outcome'] not in ['correct', 'wrong']
                ax.scatter(i, y, color=color, marker='x' if failed else 'o', s=55)
                case.append(dict(mode=mode, outcome=row['outcome'], final_qop_GeVinv=qop,
                                 truth_signed_final_qop_GeVinv=y))
            truth_qop = charge*records['production_qop'][key]['replacement_qop_GeVinv']
            ax.axhline(0, color='black', linewidth=.8)
            ax.axhline(truth_qop, color='#00838f', linestyle=':', label='Production |q/p|')
            ax.set(title=f'{key[0]} | source entry {key[1]}', ylabel=r'$q_{truth}(q/p)_{final}$ [GeV$^{-1}$]')
            ax.set_xticks(range(len(MODES)), ['Native', 'Primary', 'Centroid', 'Both', 'Seed q/p'], rotation=25, ha='right')
            limit = 1.4*max([abs(item['truth_signed_final_qop_GeVinv']) for item in case
                             if item['final_qop_GeVinv'] is not None]+[truth_qop])
            ax.set_ylim(-limit, limit)
            ax.legend(fontsize=8, loc='center', bbox_to_anchor=(.5, .65))
            counts[f'pure|{key[0]}|{key[1]}'] = dict(controls=case, production_abs_qop_GeVinv=truth_qop)
        save(fig, 'measurement_pure_wrong_cases', 'Six native purity-one wrong signs: final curvature and quality status',
             'Positive means correct sign; x marks a quality failure; primary-hit filtering alone leaves all six wrong; production q/p is not local sensor truth')

        fig, axes = plt.subplots(1, 2, figsize=(12, 5.8))
        edges = np.linspace(0, .2, 41)
        for ax, (dsid, label) in zip(axes, [('100116', r'$\mu^+$'), ('100117', r'$\mu^-$')]):
            for wrong, color, name in [(True, '#e53935', 'Native wrong'), (False, '#43a047', 'Normal controls')]:
                rows = [row for key, row in records['centroid'].items()
                        if key[0] == dsid and row['original_wrong'] == wrong]
                values = [row['maximum_common_measurement_change_mm'] for row in rows
                          if row['maximum_common_measurement_change_mm'] is not None]
                hist, _ = np.histogram(values, bins=edges)
                ax.stairs(hist, edges, color=color, label=name)
                counts[f'{dsid}|{name}|maximum_common_measurement_shift'] = dict(
                    edges=edges.tolist(), counts=hist.tolist(), raw_entries=len(values),
                    underflow=int(sum(value < edges[0] for value in values)),
                    overflow=int(sum(value > edges[-1] for value in values)))
            ax.set(title=label, xlabel='Maximum |common audited measurement shift| [mm]', ylabel='Frozen source events')
            ax.legend()
        save(fig, 'measurement_centroid_shifts', 'Observed measurement shifts from the actual SiHit-centroid control',
             'Maximum over common audited measurement hits in each event; seed means can change; measurement widths/errors, noise and missing-truth fallback are retained')
    write_json(args.output/'bin_counts.json', counts)
    write_json(args.output/'figure_index.json', index)
    (args.output/'figure_index.md').write_text('\n'.join(
        '- ['+row['name']+']('+row['png']+') · [PDF]('+row['pdf']+')' for row in index)+'\n')
    print(json.dumps(dict(figures=len(index), output=str(args.output)), indent=2))


if __name__ == '__main__':
    main()
