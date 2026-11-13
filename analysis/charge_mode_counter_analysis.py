#!/usr/bin/env python3
"""Close actual native refit counters and resolve merged-output outcomes."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from charge_modes import signatures
from charge_tables import derive
from slope_tables import merge_tables, write_json


def collect_counters(root, records, dsid, variant, mode, direction='forward'):
    arrays = defaultdict(list)
    if mode == '3forward':
        jobs = json.loads((root/'mode_analysis_v1/forward_references'/dsid/(variant+'_tasks.json')).read_text())
        for job in jobs:
            folder = Path(job['output'])
            task = json.loads((folder/'source_task.json').read_text())
            with np.load(folder/'charge_lineage.npz') as data:
                final = data['stage'] == 7
                arrays['sample_entry'].append(data['source_entry'][final]-task['skip']+task['sample_entry'])
                for key in ['ordinal', 'qop_GeVinv', 'status']:
                    arrays[key].append(data[key][final])
    else:
        group = [r for r in records if (r['dsid'], r['variant'], r['mode']) == (dsid, variant, mode)]
        group.sort(key=lambda r: r['task']['sample_entry'])
        assert len(group) == 20
        for record in group:
            check = Path(record['validation_output'])/'counter_validation.json'
            assert json.loads(check.read_text())['all_58_fields_exact']
            with np.load(Path(record['counter_source'])/('mode_refit_'+direction+'.npz')) as data:
                for key in ['sample_entry', 'ordinal', 'qop_GeVinv', 'status']:
                    arrays[key].append(data[key])
    return {key: np.concatenate(parts) for key, parts in arrays.items()}


def close_counter_tracks(counters, tracks):
    valid = counters['status'] >= 0
    for key in ['sample_entry', 'ordinal', 'qop_GeVinv']:
        np.testing.assert_array_equal(counters[key][valid], tracks[key], err_msg='Native direction counter: '+key)
    return counters['status'][valid]


def counts(status, selected):
    result = {}
    for label, mask in [('all', np.ones(len(status), dtype=bool)), ('selected', selected)]:
        row = status[mask]
        result[label] = dict(total=len(row), success=int((row == 1).sum()),
            fallback=int((row == 0).sum()), invalid=int((row < 0).sum()))
    return result


def selected_rows(tables, decision):
    mask = np.zeros(len(tables['tracks']['sample_entry']), dtype=bool)
    mask[decision['chosen'][decision['selected']]] = True
    return mask


def merged_counts(forward, backward, tracks, decision):
    lookup = defaultdict(set)
    for signatures_, status in [forward, backward]:
        for signature, value in zip(signatures_, status):
            lookup[signature].add(int(value))
    masks = selected_rows(dict(tracks=tracks), decision)
    result = {key: Counter() for key in ['all_merged_output', 'selected_merged_output']}
    for i, signature in enumerate(signatures(tracks)):
        statuses = lookup[signature]
        assert statuses, 'Merged track has no fully identical native direction track'
        label = ('ambiguous_origin_status' if len(statuses) > 1 else
                 'success' if 1 in statuses else 'fallback' if 0 in statuses else 'invalid')
        result['all_merged_output']['total'] += 1
        result['all_merged_output'][label] += 1
        if masks[i]:
            result['selected_merged_output']['total'] += 1
            result['selected_merged_output'][label] += 1
    return {key: dict(value) for key, value in result.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    recovery = args.root/'mode_counter_recovery_v2'
    records = json.loads((recovery/'effective_sources.json').read_text())
    assert len(records) == 400
    for record in records:
        check = json.loads((Path(record['validation_output'])/'counter_validation.json').read_text())
        assert check['all_58_fields_exact'] and check['fields'] == 58
        if record['action'] in ['retry_failed_reconstruction', 'validate_completed_reconstruction']:
            assert (Path(record['validation_output'])/'exit_code.txt').read_text().strip() == '0'
    expected = json.loads((args.root/'mode_analysis_v1/summary.json').read_text())['groups']
    groups, coverage = {}, {}
    for dsid in ['100116', '100117']:
        for variant in ['beforeGeo', 'afterGeo']:
            for station in ['3', '4']:
                native = {}
                for direction in ['forward', 'backward']:
                    mode = station+direction
                    source = (args.root/'mode_analysis_v1/forward_references'/dsid if mode == '3forward'
                              else args.root/'modes_100k_v1'/dsid/mode)
                    tables, cov = merge_tables(source, variant)
                    decision = derive(tables, dsid)
                    key = '|'.join([dsid, variant, mode])
                    assert int(decision['selected'].sum()) == expected[key]['conditional_misid']['total']
                    counters = collect_counters(args.root, records, dsid, variant, mode)
                    status = close_counter_tracks(counters, tables['tracks'])
                    groups[key] = dict(native_direction=counts(status, selected_rows(tables, decision)),
                        all_refit_calls=counts(counters['status'], np.zeros(len(counters['status']), dtype=bool))['all'],
                        exact_source_event_ordinal_qop_closure=True)
                    coverage[key] = cov
                    native[direction] = (signatures(tables['tracks']), status.copy(), tables['tracks'])
                    del tables, decision
                mode = station+'merge'
                tables, cov = merge_tables(args.root/'modes_100k_v1'/dsid/mode, variant)
                decision = derive(tables, dsid)
                key = '|'.join([dsid, variant, mode])
                assert int(decision['selected'].sum()) == expected[key]['conditional_misid']['total']
                streams = {}
                for direction in ['forward', 'backward']:
                    counters = collect_counters(args.root, records, dsid, variant, mode, direction)
                    status = close_counter_tracks(counters, native[direction][2])
                    np.testing.assert_array_equal(status, native[direction][1], err_msg='Merge vs independent native refit status')
                    streams[direction] = counts(counters['status'], np.zeros(len(counters['status']), dtype=bool))['all']
                groups[key] = dict(native_direction_calls=streams,
                    merged=merged_counts(native['forward'][:2], native['backward'][:2], tables['tracks'], decision),
                    exact_internal_and_independent_direction_closure=True)
                coverage[key] = cov
                print(key, groups[key], flush=True)
                del tables, decision, native
    write_json(args.output/'summary.json', dict(groups=groups, coverage=coverage,
        all_400_counter_runs_58_fields_exact=True, source_events_per_charge_mode=100000,
        scope='Actual Calypso refit outcomes, both all calls and selected tracks; merged tracks traced by complete hits, q/p and upstream parameters',
        recovery_actions=json.loads((recovery/'manifest.json').read_text())['actions']))


if __name__ == '__main__':
    main()
