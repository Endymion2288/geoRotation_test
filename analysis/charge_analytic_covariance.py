#!/usr/bin/env python3
"""Measure the isolated seed-covariance differential in actual Calypso outputs."""
import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path

from charge_component_refit import compare, components, metrics
from charge_refit_input import input_difference
from charge_ckf_chart_scan import FIELDS, STAGES
from slope_tables import write_json


def read_rows(path):
    with gzip.open(path, 'rt') as stream:
        return [json.loads(line) for line in stream]


def key(row):
    return row['series'], row['dsid'], row['entry'], row['component'], row['hit']


def analyze(root, analysis, output):
    manifest = json.loads((root/'manifest.json').read_text())
    reference = root.parent/'width_scan_v1/analysis_v2'
    old_summary = json.loads((reference/'summary.json').read_text())
    new_summary = json.loads((analysis/'summary.json').read_text())
    old_states = {key(r): r for r in read_rows(reference/'paired_state_metrics.jsonl.gz')}
    new_states = {key(r): r for r in read_rows(analysis/'paired_state_metrics.jsonl.gz')}
    common = old_states.keys() & new_states.keys()
    grouped = defaultdict(list)
    for k in common:
        a, b = old_states[k], new_states[k]
        if a['same_ordered_CKF_hits'] and b['same_ordered_CKF_hits']:
            grouped[(k[0], k[1])].append((a, b))
    output.mkdir(parents=True, exist_ok=False)
    groups, seed_rows, changes, checks = {}, [], [], []
    for point in manifest['points']:
        series, dsid = point['series'], point['dsid']
        tag = dict(series=series, dsid=dsid)
        name = series+'|'+dsid
        before, old, new = (components(Path(point[f])) for f in
                            ['beforeGeo', 'afterGeo_reference', 'afterGeo'])
        assert before.keys() == old.keys() == new.keys()
        seed_metrics = []
        for identity in before:
            a, b, c = before[identity]['seed'], old[identity]['seed'], new[identity]['seed']
            assert a['seed'] == b['seed'] == c['seed']
            # variance_qop is a redundant scalar copied from C, not a mean.
            assert b['variance_qop'] == b['cov'][4][4] and c['variance_qop'] == c['cov'][4][4]
            excluded = ('cov', 'variance_qop')
            assert {k:v for k,v in b.items() if k not in excluded} == {k:v for k,v in c.items() if k not in excluded}
            item = dict(**tag, entry=identity[0], component=identity[1],
                         finite_difference=input_difference(a, b), analytic=input_difference(a, c))
            seed_metrics.append(item)
            seed_rows.append(item)
        control, control_rows = compare(old, new)
        changes.extend(dict(**tag, **r) for r in control_rows)
        pairs = grouped[(series, dsid)]
        common_metrics = {stage: {field: {mode: metrics([p[i][stage][field] for p in pairs])
            for i,mode in enumerate(['finite_difference', 'analytic'])} for field in FIELDS}
            for stage in STAGES}
        groups[name] = dict(source_events=20, input_components=len(before),
            seed_input_metrics={mode: {field: metrics([r[mode][field] for r in seed_metrics])
                for field in ['position_difference_mm', 'direction_difference_norm',
                              'qop_difference_GeVinv', 'physical_covariance_scaled_max',
                              'direction_covariance_scaled_max']} for mode in ['finite_difference', 'analytic']},
            after_only_control=control,
            version_counts={mode: summary['groups'][name]['counts'] for mode,summary in
                            [('finite_difference', old_summary), ('analytic', new_summary)]},
            common_state_cohort=dict(states=len(pairs),
                components=len({(p[0]['entry'], p[0]['component']) for p in pairs}),
                same_ordered_hits_in_both_controls=True, metrics=common_metrics))
        checks.append(dict(**tag, input_components=len(before),
                           actual_after_seed_fields_except_covariance_exact=True))
    write_json(output/'summary.json', dict(groups=groups, source_events=60,
        scope='Fixed centers, before seeds/C, surfaces and actual Calypso settings; '
              'analytic versus historical central-difference after seed covariance only. '
              'Common state cohort is the intersection of same-hit paths in both controls. '
              'Finite Gaussian pushforward and likelihood-weighted posterior convergence are not tested.'))
    write_json(output/'after_control_paired_components.json', changes)
    with gzip.open(output/'seed_input_metrics.jsonl.gz', 'wt', compresslevel=4) as stream:
        for row in seed_rows:
            stream.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
    write_json(output/'validation.json', dict(checks=checks,
        actual_after_seed_non_C_fields_exact=True,
        all_seed_inputs_included=True, common_state_cohorts_explicit=True,
        reference_analysis=str(reference), analytic_analysis=str(analysis)))
    print(json.dumps(dict(groups=len(groups), input_components=len(seed_rows),
                          after_control_components=len(changes)), indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    analyze(args.input, args.analysis, args.output)


if __name__ == '__main__':
    main()
