#!/usr/bin/env python3
"""Compare actual native CKF chart controls in a common physical frame."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from charge_component_refit import components, compare, first, metrics
from charge_refit_input import input_contract, observer_gate, physical, refits
from slope_tables import write_json
from validate_charge_components import validate


SAMPLES = ('100044', '100116', '100117')
MODES = ('native', 'update', 'smoothing', 'both')
SERIES = 'matched_n32_f256_v1'


def source_paths(root, dsid):
    control = root/'common_prior_ckf_chart_v1'/SERIES/dsid
    reference = root/'common_prior_refit_input_v1/scale1_v1'/SERIES/dsid
    return control, reference


def gate(root, output):
    checks = []
    for dsid in SAMPLES:
        control, reference = source_paths(root, dsid)
        old = components(reference/'beforeGeo_regular_both')
        new = components(control/'native')
        check = observer_gate(old, new)
        assert refits(reference/'beforeGeo_regular_both') == refits(control/'native')
        checks.append(dict(dsid=dsid, runtime=validate(control/'native'), **check,
            refit_input_and_state_records_exact=True,
            contract=input_contract(control/'native', new, refits(control/'native'))))
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/'validation.json', dict(all_exact=True, checks=checks,
        scope='Default-off compiled shared chart adapter extraction and CKF plumbing; complete native chain replay.'))
    print(json.dumps(dict(exact_replays=len(checks), output=str(output)), indent=2))


def state_metrics(a, b, stage):
    ma, ca = physical(a, 'beforeGeo', stage)
    mb, cb = physical(b, 'afterGeo', stage)
    diagonal = np.maximum((np.abs(np.diag(ca))+np.abs(np.diag(cb)))/2., 1e-30)
    scaled = np.abs(ca-cb)/np.sqrt(np.outer(diagonal, diagonal))
    return dict(position_difference_mm=float(np.linalg.norm(ma[:3]-mb[:3])),
        direction_difference_norm=float(np.linalg.norm(ma[3:6]-mb[3:6])),
        qop_difference_GeVinv=float(abs(ma[6]-mb[6])),
        physical_covariance_scaled_max=float(np.max(scaled)),
        direction_covariance_scaled_max=float(np.max(scaled[3:6, 3:6])),
        before_qop_GeVinv=float(ma[6]), after_qop_GeVinv=float(mb[6]))


def measurement_pairs(before, after):
    output, examples = [], []
    for key, a in before.items():
        b = after[key]
        ia, ib = first(a), first(b)
        if ia is None or ib is None:
            continue
        same_hits = a['selected'][ia]['hits'] == b['selected'][ib]['hits']
        # Pair the actual measurement IDs on either path; the complete-path
        # cohort is recorded separately and is used for conditional summaries.
        def keyed(item, candidate):
            selected = [r for r in item['states'][candidate] if r['measurement_flag']]
            rows = {(r['hit'], r['outlier_flag'], r['hole_flag']): r for r in selected}
            assert len(rows) == len(selected)
            return rows
        sa, sb = keyed(a, ia), keyed(b, ib)
        common = sorted(sa.keys() & sb.keys(), key=lambda k: sa[k]['path_length_mm'])
        for hit in common:
            ra, rb = sa[hit], sb[hit]
            record = dict(entry=key[0], component=key[1], hit=hit[0],
                geometry_id_before=ra['geometry_id'], geometry_id_after=rb['geometry_id'],
                same_ordered_CKF_hits=same_hits, path_length_before_mm=ra['path_length_mm'],
                path_length_after_mm=rb['path_length_mm'])
            for stage in ('predicted', 'filtered', 'smoothed'):
                record[stage] = state_metrics(ra, rb, stage)
            output.append(record)
            if key in [(1, 23), (3, 19), (12, 2), (0, 25)]:
                examples.append(record)
    return output, examples


def analyze(root, output):
    control_root = root/'common_prior_ckf_chart_v1'
    assert json.loads((control_root/'gate_v1/validation.json').read_text())['all_exact']
    output.mkdir(parents=True, exist_ok=False)
    groups, rows, states, examples, checks = {}, [], [], [], []
    for dsid in SAMPLES:
        control, reference = source_paths(root, dsid)
        after = components(reference/'afterGeo_native')
        baseline = components(control/'native')
        for mode in MODES:
            before = baseline if mode == 'native' else components(control/mode)
            checks.append(dict(dsid=dsid, mode=mode, runtime=validate(control/mode)))
            assert baseline.keys() == before.keys() == after.keys()
            for key in baseline:
                assert baseline[key]['seed'] == before[key]['seed']
                assert before[key]['selected'].keys() == before[key]['converted'].keys() == before[key]['final'].keys()
            group, component_rows = compare(before, after)
            paired_states, case_states = measurement_pairs(before, after)
            group['all_selected_candidates'] = dict(before=sum(len(r['selected']) for r in before.values()),
                                                    after=sum(len(r['selected']) for r in after.values()))
            group['physical_state_metrics'] = {}
            for cohort, predicate in [('common_measurements', lambda r: True),
                                      ('same_ordered_CKF_hits', lambda r: r['same_ordered_CKF_hits'])]:
                available = [r for r in paired_states if predicate(r)]
                group['physical_state_metrics'][cohort] = dict(n_states=len(available),
                    components=len({(r['entry'], r['component']) for r in available}),
                    metrics={stage: {field: metrics([r[stage][field] for r in available])
                        for field in ['position_difference_mm', 'direction_difference_norm', 'qop_difference_GeVinv',
                                      'physical_covariance_scaled_max', 'direction_covariance_scaled_max']}
                        for stage in ('predicted', 'filtered', 'smoothed')})
            # Freeze a common same-hit cohort across controls below; individual
            # same-hit samples can change with CKF update or smoothing.
            name = dsid+'|'+mode
            groups[name] = group
            rows.extend(dict(dsid=dsid, mode=mode, **r) for r in component_rows)
            states.extend(dict(dsid=dsid, mode=mode, **r) for r in paired_states)
            examples.extend(dict(dsid=dsid, mode=mode, **r) for r in case_states)
    common_cohorts = {}
    for dsid in SAMPLES:
        subsets = []
        for mode in MODES:
            subsets.append({(r['entry'], r['component'], r['hit']) for r in states
                if r['dsid'] == dsid and r['mode'] == mode and r['same_ordered_CKF_hits']})
        common = set.intersection(*subsets)
        common_cohorts[dsid] = dict(n_states=len(common),
            components=len({(e, c) for e, c, _ in common}), modes={})
        for mode in MODES:
            available = [r for r in states if r['dsid'] == dsid and r['mode'] == mode
                         and (r['entry'], r['component'], r['hit']) in common]
            assert len(available) == len(common)
            common_cohorts[dsid]['modes'][mode] = {stage: {field: metrics([r[stage][field] for r in available])
                for field in ['direction_difference_norm', 'qop_difference_GeVinv', 'direction_covariance_scaled_max']}
                for stage in ('predicted', 'filtered', 'smoothed')}
    write_json(output/'summary.json', dict(groups=groups, common_cohorts=common_cohorts,
        scope='Actual CKF chart controls, identical finite-component seed inputs; correlated diagnostics, no aggregated posterior.'))
    write_json(output/'paired_components.json', rows)
    # Compact state metrics are sufficient for plots; retain full traces in
    # the original compressed native audits rather than duplicating matrices.
    import gzip
    with gzip.open(output/'paired_state_metrics.jsonl.gz', 'wt', compresslevel=4) as stream:
        for row in states:
            stream.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
    write_json(output/'case_state_metrics.json', examples)
    write_json(output/'validation.json', dict(checks=checks, complete_seed_records_exact_across_controls=True,
        source_default_gate=str(control_root/'gate_v1/validation.json'),
        all_component_candidates_retained_in_native_audits=True,
        first_candidate_selection='measurements descending, chi2 ascending, candidate ID ascending'))
    print(json.dumps(dict(groups=len(groups), paired_components=len(rows), paired_measurement_states=len(states)), indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--phase', choices=['gate', 'analyze'], required=True)
    args = parser.parse_args()
    if args.phase == 'gate':
        gate(args.root, args.output)
    else:
        analyze(args.root, args.output)


if __name__ == '__main__':
    main()
