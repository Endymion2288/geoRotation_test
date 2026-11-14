#!/usr/bin/env python3
"""Measure an encoding-only material change in actual isolated Calypso output."""
import argparse
import json
from pathlib import Path

import numpy as np

from charge_analysis import outcomes
from charge_joint_regression import compare, paired_rate, selected_momentum, summarize
from charge_refit_regression import angles
from charge_resource_analysis import resources
from charge_tables import derive as charge_decision, weighted_rate
from slope_analysis import EDGES, derive as slope_decision
from slope_statistics import peak_statistics
from slope_tables import histogram, merge_tables, slope, write_json


CONTROLS = ['legacy', 'corrected']
PLANES = ['st1', 'VetoNu', 'sensor']


def curvature_definition(tracks):
    """Read native curvature, or use the saved same-parameter Cartesian mean."""
    with np.errstate(divide='ignore', invalid='ignore'):
        calculated = tracks['charge']/np.linalg.norm(tracks['upstream_F'][:, 3:], axis=1)
    if 'qop_GeVinv' not in tracks:
        tracks['qop_GeVinv'] = calculated
        return dict(source='Saved upstream charge divided by the norm of its Cartesian momentum; '
                           'the S05 observer did not store a separate q/p branch.',
                    native_qop_variance='Not measured by the S05 observer')
    good = np.isfinite(calculated) & np.isfinite(tracks['qop_GeVinv']) & (tracks['qop_GeVinv'] != 0)
    error = float(np.max(np.abs(calculated[good]/tracks['qop_GeVinv'][good]-1.), initial=0.))
    assert error < 1e-10
    return dict(source='Native saved upstream q/p', Cartesian_norm_maximum_relative_difference=error)


def model_change(before, after):
    result = {}
    for target, name in enumerate(PLANES[:2]):
        a, b = [t['truth_model_valid'][:, target] for t in [before, after]]
        good = (a == 1) & (b == 1)
        momenta = [t['truth_model_momentum_F_GeV'][good, target] for t in [before, after]]
        difference = slope(momenta[1])-slope(momenta[0])
        relative_p = np.linalg.norm(momenta[1], axis=1)/np.linalg.norm(momenta[0], axis=1)-1.
        result[name] = dict(common_valid=int(good.sum()),
            validity_migration=np.bincount(3*(a+1)+b+1, minlength=9).reshape(3, 3).tolist(),
            slope_maximum_absolute_difference=float(np.max(np.abs(difference), initial=0.)),
            momentum_maximum_absolute_relative_difference=float(np.max(np.abs(relative_p), initial=0.)),
            slope_difference=histogram(difference, np.linspace(-1e-6, 1e-6, 201)))
    return result


def slope_summary(tables, dsid, windows, split):
    d = slope_decision(tables, 'afterGeo', dsid)
    tracks, truth = tables['tracks'], tables['truth']
    histograms, peaks, propagation = {}, {}, {}
    for i, plane in enumerate(PLANES):
        for selection in ['all', 'good95', 'matched_pure', 'reco_p_window', 'truth_p_window']:
            mask = d['selections'][selection]
            key = plane+'|'+selection
            histograms[key+'|physical'] = histogram(d['reconstructed'][mask, i], EDGES)
            histograms[key+'|residual'] = histogram(d['residuals'][mask, i], np.linspace(-.005, .005, 201))
            if windows and selection != 'good95':
                end = int(tables['events']['sample_entry'].max())+1
                peaks[key] = peak_statistics(d['reconstructed'][:, i], tracks['sample_entry'], mask,
                    windows.get(key, windows[plane]), (split, end))
        histograms[plane+'|truth_fiducial'] = histogram(d['truth_slopes'][d['fiducial_truth'], i], EDGES)
        propagation[plane] = dict(
            selected_long=int(d['long'].sum()),
            success=int(np.count_nonzero(d['long'] & (tracks['projection_valid'][:, i] == 1))),
            failed=int(np.count_nonzero(d['long'] & (tracks['projection_valid'][:, i] == 0))),
            missing=int(np.count_nonzero(d['long'] & (tracks['projection_valid'][:, i] < 0))))
    return dict(histograms=histograms, held_out_peaks=peaks, propagation=propagation,
        selections='Unchanged S05 definitions, including fiducial+H18 and purity>=0.7 for matched_pure; '
                   'all/reco-p/truth-p selections retain native long-track multiplicity.',
        fiducial_H18_truth=int(d['fiducial_truth'].sum()),
        native_long_tracks=int(d['long'].sum()), physical_definition='p_y,F/p_z,F = p_y,A/p_x,A')


def reduce(plan, root, dsid, scope, output):
    output.mkdir(parents=True, exist_ok=False)
    tables, decisions, coverage, groups, cpu, curvature = {}, {}, {}, {}, {}, {}
    for control in CONTROLS:
        source = Path(plan['sources'][scope][dsid][control])
        tables[control], coverage[control] = merge_tables(source, 'afterGeo')
        curvature[control] = curvature_definition(tables[control]['tracks'])
        # 100044 and 100116 are both mu+ (PDG=-13); the selection function's
        # DSID dispatch only determines this PDG/charge convention.
        decisions[control] = charge_decision(tables[control], '100116' if dsid == '100044' else dsid)
        jobs = json.loads((source/'afterGeo_tasks.json').read_text())
        cpu[control] = resources(jobs)
    old, new = tables['legacy'], tables['corrected']
    for name in ['truth', 'segments']:
        assert old[name].keys() == new[name].keys()
        for key in old[name]:
            if not (name == 'truth' and key.startswith('truth_model_')):
                np.testing.assert_array_equal(old[name][key], new[name][key], err_msg=name+'/'+key)
    for key in ['source_entry', 'sample_entry', 'file_index', 'run', 'event']:
        np.testing.assert_array_equal(old['events'][key], new['events'][key], err_msg='events/'+key)
    reference_change = model_change(old['truth'], new['truth'])
    native_reference_angles = angles(new, decisions['corrected'])
    # Isolate reconstruction changes from changes in the truth-model transport.
    # The native corrected model is measured above, and remains in source files.
    common = dict(new, truth=dict(new['truth']))
    for key in old['truth']:
        if key.startswith('truth_model_'):
            common['truth'][key] = old['truth'][key]
    tables['corrected'] = common
    for control in CONTROLS:
        groups[control] = summarize(tables[control], decisions[control])
    paired = compare(old, common, decisions['legacy'], decisions['corrected'])
    for field in ['common_same_candidate_angles', 'common_same_candidate_momentum']:
        paired[field] = {control: paired[field][key] for control, key in zip(CONTROLS, ['native', 'both'])}
    paired['identity_scope'] = 'Same final ordinal and identical ordered hit identifiers; no independent seed-lineage identity is asserted.'
    event_index = np.searchsorted(old['events']['sample_entry'], old['truth']['sample_entry'])
    count = len(old['events']['sample_entry'])
    relative = {c: selected_momentum(tables[c], decisions[c]) for c in CONTROLS}
    paired['momentum_tail_difference'] = paired_rate(event_index,
        relative['legacy'] < -.5, decisions['legacy']['selected'],
        relative['corrected'] < -.5, decisions['corrected']['selected'], count)
    replay = []
    for i in np.flatnonzero(decisions['legacy']['wrong'] | decisions['corrected']['wrong']):
        row = {key: old['truth'][key][i].item() for key in
               ['sample_entry', 'source_entry', 'file_index', 'run', 'event', 'barcode']}
        row['truth_p_GeV'] = float(decisions['legacy']['truth_p'][i])
        for control in CONTROLS:
            d, t = decisions[control], tables[control]['tracks']
            chosen = d['chosen'][i]
            row[control] = dict(status=int(d['status'][i]), ordinal=int(t['ordinal'][chosen]) if chosen >= 0 else None,
                qop_GeVinv=float(t['qop_GeVinv'][chosen]) if chosen >= 0 and np.isfinite(t['qop_GeVinv'][chosen]) else None)
        replay.append(row)
    write_json(output/'replay_events.json', replay)
    slopes = {}
    windows = json.loads(Path(plan['peak_windows']).read_text())
    for control in CONTROLS:
        slopes[control] = slope_summary(tables[control], dsid,
            windows['windows'] if scope == 'full' and dsid == '100044' else None,
            windows['held_out_boundary'])
    if scope == 'full' and dsid in ['100116', '100117']:
        baseline = json.loads((root/'main_analysis_v1/summary.json').read_text())
        expected = baseline['groups'][dsid+'|afterGeo']['r95']['conditional_misid']
        measured = groups['legacy']['conditional_misid']
        assert measured['numerator'] == expected['wrong']
        assert measured['denominator'] == expected['total']
        with np.load(root/'main_analysis_v1'/f'{dsid}_afterGeo_outcomes.npz') as archived:
            for key in ['sample_entry', 'file_index', 'source_entry', 'barcode']:
                np.testing.assert_array_equal(old['truth'][key], archived[key])
        with np.load(root/'production_overlap_v1/truth_weights.npz') as archive:
            weights = archive[dsid]
        for control in CONTROLS:
            d = decisions[control]
            groups[control]['production_overlap'] = {}
            for label, mask in [('all', np.ones(len(weights), dtype=bool)),
                                ('held_out', old['truth']['file_index'] >= 3)]:
                w = weights*mask
                groups[control]['production_overlap'][label] = dict(
                    conditional_misid=weighted_rate(d['wrong'], d['selected'], w),
                    correct_efficiency=weighted_rate(d['status'] == 1, d['fiducial'], w))
        outcomes(output, dsid, 'corrected', common, decisions['corrected'])
    write_json(output/'summary.json', dict(dsid=dsid, scope=scope, groups=groups,
        coverage=coverage, resources=cpu, paired=paired, slopes=slopes, curvature_definition=curvature,
        corrected_native_truth_model_angles=native_reference_angles,
        truth_model_change=reference_change, raw_truth_segments_and_event_identity_exact=True,
        reference='Both controls use the frozen legacy production-p/first-ST1-chord Calypso truth model on the same physical planes. '
                  'The corrected native truth-model change is measured separately; neither is local G4 momentum.',
        material_control='AfterGeo only, isolated smoothing ON; transpose 12 serialized rotation arrays, keep all other map entries and reconstruction settings.',
        source_plan=plan))
    print(dsid, scope, {c: groups[c]['conditional_misid'] for c in CONTROLS}, flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--dsid', choices=['100044', '100116', '100117'], required=True)
    parser.add_argument('--scope', choices=['smoke', 'full'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    reduce(json.loads(args.plan.read_text()), args.root, args.dsid, args.scope, args.output)


if __name__ == '__main__':
    main()
