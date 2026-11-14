#!/usr/bin/env python3
"""Check width sensitivity and locate residuals in the actual native CKF chain."""
import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path

from charge_ckf_chart import measurement_pairs
from charge_component_refit import compare, components, first, metrics
from charge_refit_input import input_difference, refits
from slope_tables import write_json
from validate_charge_components import validate


STAGES = ('predicted', 'filtered', 'smoothed')
FIELDS = ('position_difference_mm', 'direction_difference_norm',
          'qop_difference_GeVinv', 'physical_covariance_scaled_max',
          'direction_covariance_scaled_max')
# Descriptive thresholds locate the first measurement along a path with a
# residual. P/F/S are tested at that location, not in algorithm execution
# order: a smoothed value can include a later downstream cause.
# These thresholds are not acceptance cuts; all states enter the summaries.
DIAGNOSTIC_THRESHOLDS = dict(direction_difference_norm=1e-8,
                             qop_difference_GeVinv=1e-9,
                             direction_covariance_scaled_max=1e-5)


def source_states(item):
    candidate = first(item)
    if candidate is None:
        return {}
    return {r['hit']: r for r in item['states'][candidate] if r['measurement_flag']}


def detail(row, before, after):
    key = row['entry'], row['component']
    return dict(metrics=row, before_seed=before[key]['seed'], after_seed=after[key]['seed'],
                before_state=source_states(before[key])[row['hit']],
                after_state=source_states(after[key])[row['hit']])


def analyze(root, output):
    manifest = json.loads((root/'manifest.json').read_text())
    output.mkdir(parents=True, exist_ok=False)
    groups, checks, rows, states, extrema, growth, refit_rows = {}, [], [], [], [], [], []
    for point in manifest['points']:
        series, dsid = point['series'], point['dsid']
        tag = dict(series=series, dsid=dsid)
        folders = {variant: Path(point[variant]) for variant in ('beforeGeo', 'afterGeo')}
        for variant, folder in folders.items():
            checks.append(dict(**tag, variant=variant, runtime=validate(folder)))
        before, after = (components(folders[v]) for v in ('beforeGeo', 'afterGeo'))
        group, paired = compare(before, after)
        state_rows, _ = measurement_pairs(before, after)
        same = [r for r in state_rows if r['same_ordered_CKF_hits']]
        group['physical_state_metrics'] = dict(n_states=len(same),
            components=len({(r['entry'], r['component']) for r in same}),
            same_ordered_hits_only=True,
            metrics={stage: {field: metrics([r[stage][field] for r in same])
                for field in FIELDS} for stage in STAGES})
        # Preserve the actual bound states, C, surfaces, measurements and
        # transport Jacobians for the largest residual of every field/stage.
        for stage in STAGES:
            for field in FIELDS:
                if same:
                    worst = max(same, key=lambda r:r[stage][field])
                    extrema.append(dict(**tag, stage=stage, field=field,
                        value=worst[stage][field], **detail(worst, before, after)))
        by_component = defaultdict(list)
        for row in same:
            by_component[(row['entry'], row['component'])].append(row)
        for key, path in by_component.items():
            path.sort(key=lambda r:r['path_length_before_mm'])
            for field, threshold in DIAGNOSTIC_THRESHOLDS.items():
                found = False
                for row in path:
                    for stage in STAGES:
                        if row[stage][field] > threshold:
                            growth.append(dict(**tag, field=field, threshold=threshold,
                                first_stage=stage, **detail(row, before, after)))
                            found = True
                            break
                    if found:
                        break
        native_refits = {variant: refits(folder) for variant,folder in folders.items()}
        for pair in paired:
            if not pair['before_selected'] or not pair['after_selected']:
                continue
            key = pair['entry'], pair['component']
            a = native_refits['beforeGeo'].get((*key, first(before[key])), {})
            b = native_refits['afterGeo'].get((*key, first(after[key])), {})
            if 'refit_input' in a and 'refit_input' in b:
                refit_rows.append(dict(**tag, entry=key[0], component=key[1],
                    same_ordered_CKF_hits=pair['same_selected_hit_sequence'],
                    before_refit_success=pair['before_refit_success'],
                    after_refit_success=pair['after_refit_success'],
                    **input_difference(a['refit_input'], b['refit_input'])))
        groups[series+'|'+dsid] = group
        rows.extend(dict(**tag, **r) for r in paired)
        states.extend(dict(**tag, **r) for r in state_rows)
    write_json(output/'summary.json', dict(groups=groups,
        source_events=60, series=manifest['series'], diagnostic_thresholds=DIAGNOSTIC_THRESHOLDS,
        scope='Actual native regular CKF/refit versus rotated native, diagonal scale 1; '
              'correlated finite-prior components, no aggregated mixture posterior or production rate.'))
    write_json(output/'paired_components.json', rows)
    write_json(output/'refit_input_differences.json', refit_rows)
    write_json(output/'extreme_native_states.json', extrema)
    with gzip.open(output/'first_growth_native_states.jsonl.gz', 'wt', compresslevel=4) as stream:
        for row in growth:
            stream.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
    with gzip.open(output/'paired_state_metrics.jsonl.gz', 'wt', compresslevel=4) as stream:
        for row in states:
            stream.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
    write_json(output/'validation.json', dict(checks=checks,
        actual_seed_inputs_surfaces_C_event_coverage_and_loaded_libraries_checked=True,
        all_measurement_states_retained=True,
        extreme_states_include_full_native_mean_C_measurement_and_jacobian=True,
        growth_thresholds_are_descriptive_not_acceptance_cuts=True))
    print(json.dumps(dict(groups=len(groups), paired_components=len(rows),
        paired_measurement_states=len(states), extreme_states=len(extrema),
        first_growth_records=len(growth)), indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    analyze(args.input, args.output)


if __name__ == '__main__':
    main()
