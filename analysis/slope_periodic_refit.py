#!/usr/bin/env python3
"""Measure angular-chart and interface controls in the actual Calypso refit."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_centroids import collect
from slope_controls import archive, save
from slope_refit_covariance import paired_fraction_change
from slope_statistics import event_ratio_interval
from slope_tables import histogram, join_truth, slope, write_json


CONTROLS = {
    'periodic': ('Periodic phi difference', 'periodic phi difference'),
    'canonical': ('Canonical phi/theta chart', 'periodic phi difference and canonical phi/theta, including covariance theta sign transform on pole reflection'),
    'basis': ('Covariance basis control', 'surface/curvilinear/refit-plane local covariance Jacobians; original reference displacement retained'),
    'plane': ('Input reference plane', 'refit starts at its actual input point; native covariance basis retained'),
    'basis_plane': ('Covariance basis and input plane', 'local covariance Jacobians and refit reference at the actual input point'),
    'updater': ('Canonical update angles', 'canonicalize the actual native refit update mean and theta covariance before subsequent transport'),
    'update_canonical': ('Canonical update and smoothing', 'canonical update mean/covariance plus periodic phi and canonical smoothing angles'),
    'regular_chart': ('Beam-regular smoothing chart', 'actual ACTS smoother in rotated longitude/latitude, with full covariance and transport Jacobian transformations; native forward fit unchanged'),
    'regular_update': ('Regular direction update', 'actual ACTS measurement update in rotated longitude/latitude, converting filtered mean/covariance back before native transport; native smoother retained'),
    'regular_both': ('Regular direction update and smoothing', 'actual ACTS measurement update and smoother in rotated longitude/latitude, with full covariance and transport Jacobian transformations'),
}


def exact_tables(left, right, names):
    checks = {}
    for name in names:
        assert set(left[name]) == set(right[name]), name
        for key, values in left[name].items():
            other = right[name][key]
            equal = np.array_equal(values, other) if values.dtype.kind in 'US' else np.array_equal(values, other, equal_nan=True)
            assert equal, (name, key)
        checks[name] = len(left[name])
    return checks


def native_closure(left_paths, right_paths, forward_invariant=True):
    checks = []
    for left, right in zip(left_paths, right_paths):
        a, b = (archive(p/'native_states.npz') for p in [left, right])
        if forward_invariant:
            for key in ['entry', 'track_index', 'stage', 'station', 'layer']:
                assert np.array_equal(a[key], b[key]), (left, key)
        components = [([1, 3, 4, 5], 'seed_CKF')]
        if forward_invariant:
            components.append(([8, 9], 'refit_forward'))
        for stages, label in components:
            mask, other = (np.isin(data['stage'], stages) for data in [a, b])
            for key in a:
                assert np.array_equal(a[key][mask], b[key][other], equal_nan=True), (left, label, key)
            checks.append(dict(shard=left.name, component=label, exact_rows=int(mask.sum())))
    return checks


def refit_outcomes(paths, tracks, cohorts, count):
    lookup = {(int(e), int(o)): i for i, (e, o) in enumerate(zip(tracks['source_entry'], tracks['ordinal']))}
    success = np.full(len(lookup), -1, dtype=np.int8)
    for path in paths:
        data = archive(path/'refit_outcomes.npz')
        for e, o, value in zip(data['source_entry'], data['ordinal'], data['refit_success']):
            index = lookup[(int(e), int(o))]
            assert success[index] == -1
            success[index] = value
    assert np.isin(success, [0, 1]).all()
    return {name: dict(selected=int(selected.sum()),
                       succeeded=int((selected & (success == 1)).sum()),
                       fallback=event_ratio_interval(tracks['sample_entry'], selected & (success == 0), selected, (0, count)))
            for name, selected in cohorts.items()}


def angle_flags(paths, tracks):
    lookup = {(int(e), int(o)): i for i, (e, o) in enumerate(zip(tracks['source_entry'], tracks['ordinal']))}
    flags = {name: np.zeros(len(tracks['ordinal']), bool)
             for name in ['phi_outside', 'theta_outside', 'large_phi_difference']}
    for path in paths:
        data = archive(path/'native_angles.npz')
        indices = np.array([lookup[(int(e), int(o))] for e, o in zip(data['entry'], data['track_index'])])
        refit = data['writer'] == 1
        for name, bad in [
            ('phi_outside', np.abs(data['ePHI_smt']) > np.pi),
            ('theta_outside', (data['eTHETA_smt'] < 0) | (data['eTHETA_smt'] > np.pi)),
            ('large_phi_difference', np.abs(data['ePHI_smt']-data['ePHI_prt']) > np.pi),
        ]:
            flags[name][indices[refit & bad]] = True
    return flags


def analyse_pair(native_paths, periodic_paths, output, report, label, validation_entry=None, baseline_paths=None,
                 control='periodic', generated_st1_truth=False):
    control_label, control_definition = CONTROLS[control]
    native, count = collect(native_paths)
    periodic, other_count = collect(periodic_paths)
    assert count == other_count
    result = dict(source_events=count, label=label,
                  unchanged_source_fields=exact_tables(native, periodic, ['events', 'truth', 'segments']),
                  control=control,
                  definition=('Actual Calypso; ' + control_definition) +
                             '; fixed native candidate identity and momentum selection. Majority and complete final-hit sets checked; same-hit subsets reported separately.',
                  validation_source_entry_min=validation_entry,
                  figure_histograms={}, plane_truth_definitions={})
    if baseline_paths is not None:
        baseline, n = collect(baseline_paths)
        assert n == count
        result['off_matches_original_baseline'] = exact_tables(native, baseline, ['events', 'truth', 'tracks', 'segments'])
        for a, b in zip(native_paths, baseline_paths):
            x, y = (archive(p/'native_states.npz') for p in [a, b])
            for key in x:
                assert np.array_equal(x[key], y[key], equal_nan=True), (a, key)
        result['native_stage_closure'] = native_closure(native_paths, periodic_paths,
            forward_invariant=control in ['periodic', 'canonical', 'regular_chart'])
    rb, ra, truth = native['tracks'], periodic['tracks'], native['truth']
    for key in ['source_entry', 'ordinal']:
        assert np.array_equal(rb[key], ra[key]), key
    index, matched = join_truth(native)
    pure = matched & (np.abs(truth['pdg'][index]) == 13) & (rb['station_bits'] & 14 == 14) & (rb['purity'] >= .7)
    validation = np.ones(len(pure), bool) if validation_entry is None else rb['source_entry'] >= validation_entry
    momentum = np.linalg.norm(rb['upstream_F'][:, 3:], axis=1)
    truth_p = np.linalg.norm(truth['production_momentum_F_GeV'][index], axis=1)
    same_hits = rb['hit_signatures'] == ra['hit_signatures']
    result.update(native_candidates=len(pure), same_complete_hit_sets=int(same_hits.sum()),
                  majority_changes=int(np.count_nonzero(rb['barcode'] != ra['barcode'])), configurations={}, planes={})
    if control in ['basis', 'plane', 'basis_plane', 'updater', 'update_canonical', 'regular_chart', 'regular_update', 'regular_both']:
        fixed = pure & validation & (50 < momentum) & (momentum < 250)
        result['refit_outcomes'] = {}
        for mode, paths, r in [('native', native_paths, rb), (control, periodic_paths, ra)]:
            near = (r['projection_valid'][:, 0] == 1) & (np.abs(slope(r['projection'][:, 0, 3:6])) < .0001)
            result['refit_outcomes'][mode] = refit_outcomes(paths, r,
                dict(all_candidates=np.ones(len(r['ordinal']), bool), fixed_native_p=fixed,
                     fixed_native_p_near_zero=fixed & near), count)
    for mode, table in [('native', native), (control, periodic)]:
        r = table['tracks']
        ti, linked = join_truth(table)
        selected = linked & (np.abs(truth['pdg'][ti]) == 13) & (r['station_bits'] & 14 == 14) & (r['purity'] >= .7)
        efficient = np.bincount(ti[selected], minlength=len(truth['barcode'])) > 0
        population = np.abs(truth['pdg']) == 13
        expected_charge = np.sign(-truth['pdg'][ti])
        p = np.linalg.norm(r['upstream_F'][:, 3:], axis=1)
        result['configurations'][mode] = dict(matched_pure_long=int(selected.sum()),
            efficiency=event_ratio_interval(truth['sample_entry'], efficient & population, population, (0, count)),
            wrong_charge=event_ratio_interval(r['sample_entry'], selected & (r['charge'] != expected_charge), selected, (0, count)),
            momentum_tail=event_ratio_interval(r['sample_entry'], selected & (p < .5*np.linalg.norm(truth['production_momentum_F_GeV'][ti], axis=1)), selected, (0, count)),
            propagation={plane: dict(selected=int(selected.sum()), failed=int((selected & (r['projection_valid'][:, i] == 0)).sum()),
                                    missing=int((selected & (r['projection_valid'][:, i] < 0)).sum()))
                         for i, plane in enumerate(['st1', 'VetoNu', 'sensor'])})
    if baseline_paths is not None:
        result['raw_angle_flags'] = {}
        fixed = pure & validation & (50 < momentum) & (momentum < 250)
        for mode, paths, r in [('native', native_paths, rb), (control, periodic_paths, ra)]:
            flags = angle_flags(paths, r)
            near = np.abs(slope(r['projection'][:, 0, 3:6])) < .0001
            result['raw_angle_flags'][mode] = {name: dict(selected=int(fixed.sum()), flagged=int((fixed & bad).sum()),
                near_zero=int((fixed & near).sum()), near_zero_flagged=int((fixed & near & bad).sum()))
                for name, bad in flags.items()}
    for i, plane in enumerate(['st1', 'VetoNu', 'sensor']):
        before = np.where(rb['projection_valid'][:, i] == 1, slope(rb['projection'][:, i, 3:6]), np.nan)
        after = np.where(ra['projection_valid'][:, i] == 1, slope(ra['projection'][:, i, 3:6]), np.nan)
        reference = slope(truth['truth_model_momentum_F_GeV'][index, i]) if i < 2 else slope(truth['truth_sensor_direction_F'][index])
        reference_label = 'Same-plane truth model' if i < 2 else 'Same-sensor SiHit direction'
        if generated_st1_truth and i == 0:
            assert np.max(np.abs(truth['production_position_F_mm'][:, 2])) < 1.e-9
            reference = slope(truth['production_momentum_F_GeV'][index])
            reference_label = 'Exact generated direction at st1'
        result['plane_truth_definitions'][plane] = reference_label
        item = {}
        for selection, window in [('all', np.ones(len(pure), bool)), ('native_reco_p', (50 < momentum) & (momentum < 250)),
                                  ('truth_p', (50 < truth_p) & (truth_p < 250))]:
            selected = pure & validation & window
            for hits in [False, True]:
                finite = selected & np.isfinite(before) & np.isfinite(after)
                if hits:
                    finite &= same_hits
                key = selection + ('|same_hits' if hits else '')
                item[key] = dict(common_finite=int(finite.sum()),
                    near_zero_native=event_ratio_interval(rb['sample_entry'], finite & (np.abs(before) < .0001), finite, (0, count)),
                    **{'near_zero_'+control: event_ratio_interval(rb['sample_entry'], finite & (np.abs(after) < .0001), finite, (0, count))},
                    paired_near_zero_change=paired_fraction_change(rb['sample_entry'], np.abs(before) < .0001, np.abs(after) < .0001, finite, count),
                    native=histogram(before[finite], np.linspace(-.005, .005, 251)),
                    **{control: histogram(after[finite], np.linspace(-.005, .005, 251))})
                if control.startswith('regular_'):
                    common_truth = finite & np.isfinite(reference)
                    truth_near = np.abs(reference) < .0001
                    item[key]['truth_comparison'] = dict(
                        reference=reference_label, common_finite=int(common_truth.sum()),
                        missing_reference=int((finite & ~np.isfinite(reference)).sum()),
                        truth_near_zero=int((common_truth & truth_near).sum()),
                        native_near_zero=int((common_truth & (np.abs(before) < .0001)).sum()),
                        control_near_zero=int((common_truth & (np.abs(after) < .0001)).sum()),
                        native_minus_truth=paired_fraction_change(rb['sample_entry'], truth_near,
                            np.abs(before) < .0001, common_truth, count),
                        control_minus_truth=paired_fraction_change(rb['sample_entry'], truth_near,
                            np.abs(after) < .0001, common_truth, count))
        result['planes'][plane] = item
        selected = pure & validation & (50 < momentum) & (momentum < 250)
        finite = selected & np.isfinite(before) & np.isfinite(after)
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        panels = {}
        for ax, edges, residual in [(axes[0, 0], np.linspace(-.005, .005, 251), False),
                                   (axes[0, 1], np.linspace(-.0005, .0005, 101), False),
                                   (axes[1, 0], np.linspace(-.001, .001, 201), True)]:
            view = 'residual' if residual else ('full' if ax is axes[0, 0] else 'zoom')
            panels[view] = {}
            series = [(before, 'Native refit', 'native'), (after, control_label, control)]
            if control.startswith('regular_') and not residual:
                series.append((reference, reference_label, 'reference'))
            for values, name, mode in series:
                h = histogram((values-reference if residual else values)[finite], edges)
                ax.stairs(h['counts'], edges*1000, label=name, linestyle='--' if mode == 'reference' else '-')
                panels[view][mode] = (dict(reference_path=f'planes/{plane}/native_reco_p/{mode}')
                                     if view == 'full' and mode != 'reference' else h)
            ax.set(xlabel='Reco sy - reference sy [mrad]' if residual else 'Physical sy [mrad]', ylabel='Fixed native candidates / bin')
            ax.grid(alpha=.2)
            ax.legend(fontsize=8)
        result['figure_histograms'][plane] = dict(
            units='Dimensionless sy and delta sy; figure axes display mrad',
            reference=reference_label, panels=panels,
            fraction_panel_reference=f'planes/{plane}/{{all,native_reco_p,truth_p}}')
        keys = ['all', 'native_reco_p', 'truth_p']
        for offset, (name, field) in enumerate([('Native refit', 'near_zero_native'), (control_label, 'near_zero_'+control)]):
            values = [item[k][field] for k in keys]
            y = np.array([v['estimate'] for v in values])
            ci = np.array([v['interval_95'] for v in values]).T
            axes[1, 1].errorbar(np.arange(3)+offset*.12, y, yerr=np.maximum(0, [y-ci[0], ci[1]-y]), fmt='o', label=name)
        axes[1, 1].set(xticks=np.arange(3)+.06, xticklabels=['All p', 'Native p cut', 'Truth p cut'],
                       ylabel='Fraction within |sy| < 100 microradians')
        axes[1, 1].legend(fontsize=8)
        axes[1, 1].grid(alpha=.2)
        fig.suptitle(f'{label} | {plane} | Actual native Calypso | {count:,} source events\nFixed native primary-muon candidates, purity >= 0.7; p = 50-250 GeV; common projections')
        save(fig, output, label+'_'+plane+'_'+control+'_refit', report)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['rdo', 'gun'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--native-root', type=Path, help='Original off production root when the control is in a different frozen run')
    parser.add_argument('--control', choices=list(CONTROLS), default='periodic')
    parser.add_argument('--native-mode', default='native', help='Production subdirectory containing the native comparator')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.output/(args.control+'_refit_report.pdf')) as report:
        if args.mode == 'rdo':
            native = sorted(((args.native_root or args.root)/('production/'+args.native_mode+'/shards/beforeGeo')).glob('s*'))
            periodic = sorted((args.root/('production/'+args.control+'/shards/beforeGeo')).glob('s*'))
            baseline = sorted((args.baseline/'shards/beforeGeo').glob('s*'))
            assert len(native) == len(periodic) == len(baseline) == 2
            cases = {'production_muplus': analyse_pair(native, periodic, args.output, report, 'production_muplus', 250000, baseline, args.control)}
        else:
            sources = json.loads((args.root/'manifest.json').read_text())['sources']
            cases = {}
            for source in sources:
                p = Path(source)
                task = json.loads((p/'source_task.json').read_text())
                label = 'muplus' if task['pdg'] == -13 else 'muminus'
                cases[str(task['pdg'])] = analyse_pair([p/'beforeGeo_native'], [p/('beforeGeo_'+args.control)], args.output, report, label, control=args.control,
                    generated_st1_truth=args.control.startswith('regular_') and task.get('gun_z_F_mm') == 0.)
    write_json(args.output/'summary.json', dict(mode=args.mode, control=args.control, cases=cases))
    print(args.output/'summary.json')


if __name__ == '__main__':
    main()
