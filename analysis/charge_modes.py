#!/usr/bin/env python3
"""Compare actual Calypso directions/stations and close native merge origins."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from charge_analysis import paired_difference
from charge_history import PBINS
from charge_tables import STATUS, derive, measurement
from slope_tables import merge_tables, write_json


MODES = ['3forward', '3backward', '3merge', '4forward', '4backward', '4merge']


def forward_reference(root, output, dsid):
    source = root/f'main_{dsid}_v1'
    original = json.loads((source/'manifest.json').read_text())
    reference = output/'forward_references'/dsid
    reference.mkdir(parents=True)
    for variant in ['beforeGeo', 'afterGeo']:
        tasks = json.loads((source/(variant+'_tasks.json')).read_text())
        tasks = [task for task in tasks if
                 json.loads((Path(task['output'])/'source_task.json').read_text())['skip'] < 20000]
        assert len(tasks) == 20
        write_json(reference/(variant+'_tasks.json'), tasks)
    write_json(reference/'manifest.json', dict(original, source_events_per_variant=100000,
        ranges=[r for r in original['ranges'] if r['skip'] < 20000],
        source_rule='Existing main shards: first 20,000 entries of every source file'))
    return reference


def signatures(tracks):
    """Full hits and all saved upstream parameters, without a lossy hash."""
    result = []
    for i, entry in enumerate(tracks['sample_entry']):
        hits = tracks['hit_identifiers'][tracks['hit_offsets'][i]:tracks['hit_offsets'][i+1]]
        result.append((int(entry), tuple(sorted(map(int, hits))),
                       float(tracks['qop_GeVinv'][i]), tuple(tracks['upstream_raw'][i])))
    return result


def merge_origins(forward, backward, merged, selected_rows, wrong_rows):
    sources = [Counter(forward), Counter(backward)]
    outputs = Counter(merged)
    assert all(n <= sources[0][key]+sources[1][key] for key, n in outputs.items())
    counts = {scope: Counter() for scope in ['all_output', 'selected', 'wrong']}
    for i, key in enumerate(merged):
        a, b = key in sources[0], key in sources[1]
        assert a or b, 'Actual merge track has no exact native direction source'
        origin = 'both' if a and b else 'forward' if a else 'backward'
        counts['all_output'][origin] += 1
        if i in selected_rows:
            counts['selected'][origin] += 1
        if i in wrong_rows:
            counts['wrong'][origin] += 1
    return dict(forward_native=len(forward), backward_native=len(backward),
        merged_output=len(merged), exact_full_hit_qop_upstream_closure=True,
        counts={scope: dict(row) for scope, row in counts.items()})


def compact(tables, d):
    primary = d['primary']
    truth, tracks = tables['truth'], tables['tracks']
    rows = d['chosen'][primary]
    valid = rows >= 0
    result = {key: truth[key][primary] for key in ['sample_entry', 'file_index', 'barcode']}
    for key in ['status', 'selected', 'wrong', 'fiducial']:
        result[key] = d[key][primary]
    result['truth_p_GeV'] = d['truth_p'][primary]
    for key in ['qop_GeVinv', 'purity']:
        values = np.full(len(rows), np.nan)
        values[valid] = tracks[key][rows[valid]]
        result[key] = values
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    groups, migrations, origins, coverage = {}, {}, {}, {}
    for dsid in ['100116', '100117']:
        reference = forward_reference(args.root, args.output, dsid)
        for variant in ['beforeGeo', 'afterGeo']:
            decisions, direction_signatures = {}, {}
            for mode in MODES:
                source = reference if mode == '3forward' else args.root/'modes_100k_v1'/dsid/mode
                tables, cov = merge_tables(source, variant)
                d = derive(tables, dsid)
                key = '|'.join([dsid, variant, mode])
                groups[key] = measurement(tables, d)
                groups[key]['truth_p_bins'] = [dict(low_GeV=float(lo), high_GeV=float(hi),
                    **measurement(tables, d, (d['truth_p'] >= lo) & (d['truth_p'] < hi)))
                    for lo, hi in zip(PBINS[:-1], PBINS[1:])]
                groups[key]['refit_fallback'] = dict(value=None,
                    reason='Mode production did not enable refit counters; MessageSvc warning suppression prevents exact log counts')
                coverage[key] = cov
                row = compact(tables, d)
                np.savez_compressed(args.output/(key.replace('|', '_')+'_outcomes.npz'), **row)
                decisions[mode] = row
                current = signatures(tables['tracks'])
                if mode.endswith('merge'):
                    selected_rows = set(d['chosen'][d['selected']])
                    wrong_rows = set(d['chosen'][d['wrong']])
                    origins[key] = merge_origins(direction_signatures.pop(mode[0]+'forward'),
                        direction_signatures.pop(mode[0]+'backward'), current, selected_rows, wrong_rows)
                else:
                    direction_signatures[mode] = current
                print(key, groups[key]['conditional_misid'], flush=True)
                del tables, d, current
            baseline = decisions['3forward']
            for mode, current in decisions.items():
                for field in ['sample_entry', 'barcode', 'fiducial', 'truth_p_GeV']:
                    np.testing.assert_array_equal(baseline[field], current[field], err_msg=mode+'/'+field)
                matrix = np.bincount(7*baseline['status']+current['status'], minlength=49).reshape(7, 7)
                migrations['|'.join([dsid, variant, mode])] = dict(states=STATUS,
                    counts=matrix.tolist(), conditional_rate_difference=paired_difference(
                        baseline, current, np.ones(len(current['status']), dtype=bool)))
    result = dict(groups=groups, coverage=coverage, migrations=migrations, merge_origins=origins,
        source_events_per_charge_variant_mode=100000, source_rule='First 20,000 of each of five RDO files',
        runtime='Actual native Calypso before/after CKF and MergeTracksAlg; isolated smoothing ON',
        selection='Same primary/SDO matching and r95 selection as S06 main; require downstream stations 1-3 in every mode')
    write_json(args.output/'summary.json', result)


if __name__ == '__main__':
    main()
