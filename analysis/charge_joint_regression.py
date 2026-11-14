#!/usr/bin/env python3
"""Joint charge, momentum, angle, efficiency and fallback regression in Calypso."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from charge_history import PBINS, rate
from charge_mode_counter_analysis import close_counter_tracks, counts, merged_counts, selected_rows
from charge_modes import merge_origins, signatures
from charge_refit_regression import angles
from charge_resource_analysis import resources
from charge_tables import STATUS, derive, rank_candidates
from slope_tables import histogram, join_truth, merge_tables, slope, write_json


def cluster_rate(event_index, numerator, denominator, event_count):
    """Use source events as independent units, including both A-prime daughters."""
    den = np.bincount(event_index[denominator], minlength=event_count)
    num = np.bincount(event_index[numerator & denominator], minlength=event_count)
    total = int(den.sum())
    result = dict(numerator=int(num.sum()), denominator=total, source_events=event_count,
                  value=None, standard_error=None, interval95=[None, None],
                  method='Source-event cluster delta method')
    if not total:
        return result
    value = float(num.sum()/total)
    error = float(np.sqrt(event_count/max(1, event_count-1)*np.square(num-value*den).sum())/total)
    interval = [max(0., value-1.95996398454*error), min(1., value+1.95996398454*error)]
    if value in [0., 1.]:
        # With no successes, event resampling cannot provide a nonzero limit.
        # Bound any-success events, then allow every exposed particle in such
        # an event to be a success. The same construction bounds failures at 1.
        upper = min(1., rate(0, event_count)['intervals']['95'][1]*event_count*int(den.max())/total)
        interval = [0., upper] if value == 0. else [1.-upper, 1.]
        result['method'] += '; conservative any-event boundary bound'
    result.update(value=value, standard_error=error, interval95=interval)
    return result


def paired_rate(event_index, before_num, before_den, after_num, after_den, event_count):
    terms, values = [], []
    for num, den in [(before_num, before_den), (after_num, after_den)]:
        if not den.any():
            return dict(value=None, standard_error=None, interval95=[None, None])
        value = np.count_nonzero(num & den)/den.sum()
        values.append(value)
        terms.append((num.astype(float)*den-value*den)/den.sum())
    clustered = np.bincount(event_index, weights=terms[1]-terms[0], minlength=event_count)
    error = float(np.sqrt(event_count/max(1, event_count-1)*np.square(clustered).sum()))
    value = float(values[1]-values[0])
    interval = [value-1.95996398454*error, value+1.95996398454*error]
    method = 'Paired source-event cluster delta method'
    if error == 0.:
        # No observed paired fluctuations do not prove an exactly zero change.
        # Bound unobserved event changes using the largest observed exposure in
        # each arm, conditional on the observed denominator exposures.
        exposure = sum(np.bincount(event_index[den], minlength=event_count).max()/den.sum()
                       for den in [before_den, after_den])
        upper = min(1., rate(0, event_count)['intervals']['95'][1]*event_count*exposure)
        interval = [max(-1., value-upper), min(1., value+upper)]
        method += '; conservative any-event boundary bound conditional on observed exposures'
    return dict(value=value, standard_error=error, interval95=interval, method=method)


def derive_species(tables, dsid):
    if dsid in ['100116', '100117']:
        return derive(tables, dsid)
    truth, tracks, events = [tables[k] for k in ['truth', 'tracks', 'events']]
    primary = (np.abs(truth['pdg']) == 11) & (truth['status'] == 1) & \
        (truth['barcode'] > 0) & (truth['barcode'] < 200000)
    fiducial = primary & (np.linalg.norm(truth['production_position_F_mm'][:, :2], axis=1) < 100) & \
        np.all(np.linalg.norm(truth['truth_station_position_F_mm'][:, 1:, :2], axis=2) < 100, axis=1)
    index, matched = join_truth(tables)
    eligible = np.zeros(len(index), dtype=bool)
    eligible[matched] = primary[index[matched]]
    chosen, duplicates = rank_candidates(index, eligible, tracks, len(primary))
    present = chosen >= 0
    safe = np.maximum(chosen, 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        quality = (tracks['station_bits'] & 14 == 14) & (tracks['ndof'] > 0) & \
            (tracks['chi2']/tracks['ndof'] < 25) & (tracks['n_layers'] >= 7) & \
            (tracks['ndof']+5 >= 12) & (tracks['rmax_mm'] < 95)
    qop = tracks['qop_GeVinv']
    valid = np.isfinite(qop) & (qop != 0) & np.isin(tracks['charge'], [-1, 1])
    np.testing.assert_array_equal(np.sign(qop[valid]), tracks['charge'][valid])
    pure = tracks['purity'] > .7
    status = np.zeros(len(primary), dtype=np.int8)
    if len(qop):
        status[present & ~quality[safe]] = 5
        status[present & quality[safe] & ~pure[safe]] = 4
        status[present & quality[safe] & pure[safe] & ~valid[safe]] = 3
        selected = present & quality[safe] & pure[safe] & valid[safe]
        status[selected] = np.where(tracks['charge'][safe[selected]]*(-np.sign(truth['pdg'][selected])) < 0, 2, 1)
    status[~fiducial] = 6
    return dict(primary=primary, fiducial=fiducial, chosen=chosen, duplicates=duplicates,
        truth_p=np.linalg.norm(truth['production_momentum_F_GeV'], axis=1),
        track_quality=quality, track_valid=valid, track_matched=eligible,
        status=status, selected=(status == 1) | (status == 2), wrong=status == 2,
        missing_track_parameters=int(events['n_tracks'].sum()-len(qop)))


def selected_momentum(tables, decision):
    values = np.full(len(decision['selected']), np.nan)
    selected = decision['selected']
    qop = tables['tracks']['qop_GeVinv'][decision['chosen'][selected]]
    values[selected] = 1./np.abs(qop)/decision['truth_p'][selected]-1.
    return values


def summarize(tables, decision):
    truth, tracks, events = [tables[k] for k in ['truth', 'tracks', 'events']]
    index = np.searchsorted(events['sample_entry'], truth['sample_entry'])
    np.testing.assert_array_equal(events['sample_entry'][index], truth['sample_entry'])
    n = len(events['sample_entry'])
    selected, wrong, fid = [decision[k] for k in ['selected', 'wrong', 'fiducial']]
    relative = selected_momentum(tables, decision)
    rows = decision['chosen'][selected]
    qop = tracks['qop_GeVinv'][rows]
    quality = decision['track_quality']
    result = dict(source_events=n, primary_truth=int(decision['primary'].sum()),
        conditional_misid=cluster_rate(index, wrong, selected, n),
        correct_efficiency=cluster_rate(index, decision['status'] == 1, fid, n),
        wrong_efficiency=cluster_rate(index, wrong, fid, n),
        selected_efficiency=cluster_rate(index, selected, fid, n),
        momentum_below_half_production=cluster_rate(index, relative < -.5, selected, n),
        statuses={name: int(np.count_nonzero(decision['status'] == i)) for i, name in enumerate(STATUS)},
        duplicates=cluster_rate(index, decision['duplicates'] > 1, fid, n),
        quality_track_counts=dict(total=int(quality.sum()),
            unmatched_primary=int(np.count_nonzero(quality & ~decision['track_matched'])),
            low_purity=int(np.count_nonzero(quality & (tracks['purity'] <= .7))),
            invalid_curvature=int(np.count_nonzero(quality & ~decision['track_valid']))),
        track_parameters_missing=decision['missing_track_parameters'],
        angle=angles(tables, decision),
        momentum_residual_quantiles=np.quantile(relative[selected], [.025, .16, .5, .84, .975]).tolist() if selected.any() else None,
        histograms=dict(relative_p=histogram(relative[selected], np.linspace(-1.2, 1.2, 101)),
            signed_qop=histogram(qop, np.linspace(-.03, .03, 301)),
            wrong_signed_qop=histogram(qop[wrong[selected]], np.linspace(-.03, .03, 301)),
            truth_p=histogram(decision['truth_p'][selected], np.linspace(0, 5000, 201)),
            reco_p=histogram(1./np.abs(qop), np.linspace(0, 5000, 201))))
    result['truth_p_bins'] = []
    angle_valid = np.zeros(len(selected), dtype=bool)
    angle_residual = np.full(len(selected), np.nan)
    chosen_truth = np.flatnonzero(selected)
    valid = (tracks['projection_valid'][rows, 0] == 1) & (truth['truth_model_valid'][chosen_truth, 0] == 1)
    residual = slope(tracks['projection'][rows, 0, 3:6])-slope(truth['truth_model_momentum_F_GeV'][chosen_truth, 0])
    valid &= np.isfinite(residual)
    angle_valid[chosen_truth] = valid
    angle_residual[chosen_truth] = residual
    result['angle_above_1mrad'] = cluster_rate(index, np.abs(angle_residual) > .001, angle_valid, n)
    result['truth_charge_groups'] = {str(charge): dict(
        conditional_misid=cluster_rate(index, wrong, selected & (-np.sign(truth['pdg']) == charge), n),
        correct_efficiency=cluster_rate(index, decision['status'] == 1, fid & (-np.sign(truth['pdg']) == charge), n))
        for charge in [-1, 1]}
    for lo, hi in zip(PBINS[:-1], PBINS[1:]):
        mask = (decision['truth_p'] >= lo) & (decision['truth_p'] < hi)
        result['truth_p_bins'].append(dict(low_GeV=float(lo), high_GeV=float(hi),
            conditional_misid=cluster_rate(index, wrong, selected & mask, n),
            correct_efficiency=cluster_rate(index, decision['status'] == 1, fid & mask, n)))
    positions = truth['truth_station_position_F_mm'][:, 1]
    result['st1_truth_regions'] = {}
    for xsign in [-1, 1]:
        for ysign in [-1, 1]:
            mask = (positions[:, 0]*xsign >= 0) & (positions[:, 1]*ysign >= 0)
            result['st1_truth_regions'][f'x{xsign:+d}_y{ysign:+d}'] = dict(
                conditional_misid=cluster_rate(index, wrong, selected & mask, n),
                correct_efficiency=cluster_rate(index, decision['status'] == 1, fid & mask, n))
    return result


def compare(before, after, old, new):
    for name in ['truth', 'segments']:
        assert set(before[name]) == set(after[name])
        for field in before[name]:
            np.testing.assert_array_equal(before[name][field], after[name][field], err_msg=name+'/'+field)
    np.testing.assert_array_equal(before['events']['sample_entry'], after['events']['sample_entry'])
    common = old['selected'] & new['selected']
    same = np.zeros(len(common), dtype=bool)
    for i in np.flatnonzero(common):
        a, b = old['chosen'][i], new['chosen'][i]
        left, right = before['tracks'], after['tracks']
        if left['ordinal'][a] == right['ordinal'][b]:
            x = left['hit_identifiers'][left['hit_offsets'][a]:left['hit_offsets'][a+1]]
            y = right['hit_identifiers'][right['hit_offsets'][b]:right['hit_offsets'][b+1]]
            same[i] = np.array_equal(x, y)
    e = np.searchsorted(before['events']['sample_entry'], before['truth']['sample_entry'])
    n = len(before['events']['sample_entry'])
    difference = paired_rate(e, old['wrong'], old['selected'], new['wrong'], new['selected'], n)
    efficiency = paired_rate(e, old['status'] == 1, old['fiducial'], new['status'] == 1, new['fiducial'], n)
    return dict(truth_and_segments_exact=True, common_selected_truth=int(common.sum()),
        same_candidate_and_ordered_hits=int(same.sum()),
        migration_states=STATUS, migration=np.bincount(7*old['status']+new['status'], minlength=49).reshape(7, 7).tolist(),
        conditional_misid_difference=difference, correct_efficiency_difference=efficiency,
        common_same_candidate_angles=dict(native=angles(before, old, same), both=angles(after, new, same)),
        common_same_candidate_momentum=dict(native=histogram(selected_momentum(before, old)[same], np.linspace(-1.2, 1.2, 101)),
            both=histogram(selected_momentum(after, new)[same], np.linspace(-1.2, 1.2, 101))))


def collect_counters(root, group, records, direction):
    parts = defaultdict(list)
    for job in json.loads((Path(group['group'])/'beforeGeo_tasks.json').read_text()):
        folder = Path(job['output'])
        source = json.loads((folder/'source_task.json').read_text())
        if group['reused_native'] and group['mode'] == '3forward':
            with np.load(folder/'charge_lineage.npz') as data:
                mask = data['stage'] == 7
                parts['sample_entry'].append(data['source_entry'][mask]-source['skip']+source['sample_entry'])
                for field in ['ordinal', 'qop_GeVinv', 'status']:
                    parts[field].append(data[field][mask])
        else:
            if group['reused_native']:
                key = (group['dsid'], group['mode'], source['file_index'], source['skip'])
                record = records[key]
                assert json.loads((Path(record['validation_output'])/'counter_validation.json').read_text())['all_58_fields_exact']
                folder = Path(record['counter_source'])
            label = direction if group['mode'].endswith('merge') else 'forward'
            with np.load(folder/('mode_refit_'+label+'.npz')) as data:
                for field in ['sample_entry', 'ordinal', 'qop_GeVinv', 'status']:
                    parts[field].append(data[field])
    return {key: np.concatenate(value) for key, value in parts.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    manifest = json.loads((args.input/'manifest.json').read_text())
    groups = {'|'.join(g[k] for k in ['dsid', 'mode', 'control']): g for g in manifest['groups']}
    records = {(r['dsid'], r['mode'], r['task']['file_index'], r['task']['skip']): r for r in
               json.loads((args.root/'mode_counter_recovery_v2/effective_sources.json').read_text()) if r['variant'] == 'beforeGeo'}
    summaries, comparisons, coverage, counter_summary, resource_summary = {}, {}, {}, {}, {}
    for dsid in manifest['source_events_per_species']:
        for station in ['3', '4']:
            native_paths = {control: {} for control in ['native', 'both']}
            for direction in ['forward', 'backward', 'merge']:
                mode = station+direction
                data = {}
                for control in ['native', 'both']:
                    key = '|'.join([dsid, mode, control])
                    group = groups[key]
                    tables, cov = merge_tables(Path(group['group']), 'beforeGeo')
                    decision = derive_species(tables, dsid)
                    summaries[key] = summarize(tables, decision)
                    coverage[key] = cov
                    jobs = json.loads((Path(group['group'])/'beforeGeo_tasks.json').read_text())
                    resource_summary[key] = resources(jobs)
                    if direction != 'merge':
                        counters = collect_counters(args.root, group, records, direction)
                        status = close_counter_tracks(counters, tables['tracks'])
                        counter_summary[key] = counts(status, selected_rows(tables, decision))
                        native_paths[control][direction] = signatures(tables['tracks']), status, tables['tracks']
                    else:
                        streams = {}
                        for stream in ['forward', 'backward']:
                            counters = collect_counters(args.root, group, records, stream)
                            status = close_counter_tracks(counters, native_paths[control][stream][2])
                            np.testing.assert_array_equal(status, native_paths[control][stream][1])
                            streams[stream] = counts(counters['status'], np.zeros(len(counters['status']), dtype=bool))['all']
                        counter_summary[key] = dict(native_direction_calls=streams,
                            merged=merged_counts(native_paths[control]['forward'][:2], native_paths[control]['backward'][:2], tables['tracks'], decision))
                        selected_set, wrong_set = [set(decision['chosen'][decision[k]]) for k in ['selected', 'wrong']]
                        summaries[key]['merge_origins'] = merge_origins(native_paths[control]['forward'][0],
                            native_paths[control]['backward'][0], signatures(tables['tracks']), selected_set, wrong_set)
                    keep = {k: tables['truth'][k] for k in ['sample_entry', 'barcode', 'pdg']}
                    keep.update(status=decision['status'], selected=decision['selected'], wrong=decision['wrong'])
                    np.savez_compressed(args.output/(key.replace('|', '_')+'_outcomes.npz'), **keep)
                    data[control] = tables, decision
                    print(key, summaries[key]['conditional_misid'], flush=True)
                comparisons[dsid+'|'+mode] = compare(data['native'][0], data['both'][0], data['native'][1], data['both'][1])
                del data
            del native_paths
    write_json(args.output/'summary.json', dict(groups=summaries, comparisons=comparisons,
        coverage=coverage, refit_counters=counter_summary, resources=resource_summary,
        source_manifest=str(args.input/'manifest.json'), selection=manifest['holdout'],
        uncertainties='Source-event cluster delta; two electron daughters are not independent; paired changes use the same source events',
        reference='Relative momentum uses HepMC production p. Angle uses actual Calypso propagation of the production-p/first-st1-SiHit-chord model to z_F=0 mm; not local G4 p.',
        control='Only the established joint regular chart in native beforeGeo refit; isolated CKF smoothing ON in all paths'))


if __name__ == '__main__':
    main()
