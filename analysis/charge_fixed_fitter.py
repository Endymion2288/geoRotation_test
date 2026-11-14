#!/usr/bin/env python3
"""Compare actual native finder/fitter coverage for all frozen components."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path

from charge_component_refit import audit, components
from charge_forward_covariance import forward_states
from slope_tables import write_json


def terminal_records(folder):
    event, records = -1, {}
    for row in audit(folder/'prior_audit.jsonl.gz'):
        if row['stage'] == 'event':
            event += 1
        if row['stage'] in ('component_end', 'component_finder_error'):
            key = event, row['input_component']
            assert key not in records
            records[key] = row
    return records


def measurements(item):
    if len(item['selected']) != 1:
        return []
    candidate = next(iter(item['selected']))
    return [r for r in forward_states(item, candidate) if r['measurement_flag']]


def analyze(root, analysis):
    assert json.loads((analysis/'validation.json').read_text())['all_input_components_retained']
    manifest = json.loads((root/'manifest.json').read_text())
    assert manifest['component_solver'] == 'native_fixed_source_link_fitter'
    observations = json.loads((root/'observations.json').read_text())
    points = {(p['dsid'], p['variant'], p['mode']):Path(p['output']) for p in manifest['points']}
    results = {(r['dsid'], r['variant'], r['entry'], r['component']):r
               for r in audit(analysis/'components.jsonl.gz')}
    groups, rows = {}, []
    for dsid in ('100044', '100116', '100117'):
        for variant in ('beforeGeo', 'afterGeo'):
            old = components(points[dsid, variant, 'default'])
            new = components(points[dsid, variant, 'fixed_observation'])
            assert old.keys() == new.keys()
            terminals = {mode:terminal_records(points[dsid, variant, mode])
                         for mode in ('default', 'fixed_observation')}
            assert all(records.keys() == old.keys() for records in terminals.values())
            group = Counter(input_components=len(old))
            for key, a in old.items():
                b = new[key]
                assert a['seed'] == b['seed'], (dsid, variant, key, 'input seed differs')
                expected = observations[dsid]['events'][key[0]]['forward_ordered_hits']
                paths = {mode:measurements(item) for mode, item in [('default', a), ('fixed_observation', b)]}
                hit_lists = {mode:[r['hit'] for r in path] for mode, path in paths.items()}
                complete = {mode:len(item['selected']) == 1 and hit_lists[mode] == expected and
                            not any(r['outlier_flag'] for r in path)
                            for mode, item, path in [('default', a, paths['default']),
                                                     ('fixed_observation', b, paths['fixed_observation'])]}
                result_key = (dsid, variant, *key)
                assert not results[result_key]['complete'] or complete['fixed_observation']
                for mode, passed in complete.items():
                    group[mode+'_complete_components'] += passed
                transition = ('complete' if complete['default'] else 'incomplete')+'->'+(
                    'complete' if complete['fixed_observation'] else 'incomplete')
                group[transition] += 1
                row = dict(dsid=dsid, variant=variant, entry=key[0], component=key[1],
                    transition=transition, expected_hits=expected, actual_hits=hit_lists,
                    complete_observations=complete,
                    fixed_fitter_complete_likelihood=results[result_key]['complete'],
                    native_terminal_records={mode:records[key] for mode, records in terminals.items()})
                old_hits, new_hits = (hit_lists[m] for m in ('default', 'fixed_observation'))
                prefix = 0
                for ha, hb in zip(old_hits, new_hits):
                    if ha != hb:
                        break
                    prefix += 1
                row['common_measurement_prefix'] = prefix
                row['prefix_PF_parameters_exact'] = all(
                    ra[stage] == rb[stage]
                    for ra, rb in zip(paths['default'][:prefix], paths['fixed_observation'][:prefix])
                    for stage in ('predicted', 'filtered'))
                row['prefix_max_abs_qop_difference_GeVinv'] = max((
                    abs(ra[stage][4]-rb[stage][4])
                    for ra, rb in zip(paths['default'][:prefix], paths['fixed_observation'][:prefix])
                    for stage in ('predicted', 'filtered')), default=0.)
                if len(a['selected']) == 1:
                    missing = [hit for hit in expected if hit not in set(old_hits)]
                    row['first_missing_expected_hit'] = missing[0] if missing else None
                    if missing:
                        hit = missing[0]
                        matched = next((r for r in paths['fixed_observation'] if r['hit'] == hit), None)
                        row['first_missing_expected_ordinal'] = expected.index(hit)
                        row['first_missing_geometry_id'] = matched['geometry_id'] if matched else None
                        row['fixed_fitter_prediction_at_first_missing'] = (
                            dict(bound=matched['predicted'], surface=matched['surface']) if matched else None)
                else:
                    row['first_missing_expected_hit'] = None
                for mode, records in terminals.items():
                    terminal = records[key]
                    if terminal['stage'] == 'component_finder_error':
                        group[mode+'_raw_error_'+str(terminal['error'])] += 1
                rows.append(row)
            for mode in ('default', 'fixed_observation'):
                complete_events = 0
                for event in range(20):
                    relevant = [r for r in rows if r['dsid'] == dsid and r['variant'] == variant and r['entry'] == event]
                    if relevant:
                        complete_events += all(r['complete_observations'][mode] for r in relevant)
                group[mode+'_complete_observation_events'] = complete_events
            groups[dsid+'|'+variant] = dict(group)
    with gzip.open(analysis/'solver_coverage_rows.jsonl.gz', 'wt', compresslevel=4) as stream:
        for row in rows:
            stream.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
    write_json(analysis/'solver_coverage.json', dict(groups=groups, input_component_records=len(rows),
        common_inputs_exact=True, isolated_smoothing=True,
        source_events=60, native_geometry_navigator=True,
        fitter_external_surfaces='The existing ACTS KalmanFitter inserts measurement surfaces as external surfaces, ignoring boundary checks',
        error_scope='Native numeric errors retained together with recorded categories/messages; no category is inferred from an integer.',
        missing_scope='First missing observation identified only for a surviving native selected path. Failed paths have no saved state sequence; their first missing surface is unmeasured.'))
    print(json.dumps(groups))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--analysis', type=Path, required=True)
    args = parser.parse_args()
    analyze(args.input, args.analysis)


if __name__ == '__main__':
    main()
