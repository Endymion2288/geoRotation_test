#!/usr/bin/env python3
"""Audit explicit step limits inside the actual Calypso CKF chain."""
import argparse
from collections import Counter
import gzip
from itertools import zip_longest
import json
from pathlib import Path

import numpy as np

from bound_surface_coordinates import direction_f
from charge_component_refit import audit, components, first, metrics, qop, sign
from charge_forward_covariance import common, forward_states, matrix_metrics, signature
from slope_tables import write_json
from validate_charge_components import validate


MODES = ('default', 'step100mm', 'step10mm')
SAMPLES = ('100044', '100116', '100117')


def gate(root, output):
    checks = []
    for point in json.loads((root/'manifest.json').read_text())['points']:
        old, new = Path(point['reference']), Path(point['default'])
        assert (new/'exit_code.txt').read_text().strip() == '0'
        count = 0
        for a, b in zip_longest(audit(old/'prior_audit.jsonl.gz'), audit(new/'prior_audit.jsonl.gz')):
            assert a == b, (point['dsid'], point['variant'], count)
            count += 1
        smoothing = 0
        for a, b in zip_longest(audit(old/'smoothing_operations.jsonl.gz'), audit(new/'smoothing_operations.jsonl.gz')):
            assert a == b, (point['dsid'], point['variant'], 'smoothing', smoothing)
            smoothing += 1
        fields = 0
        for name in ('events', 'truth', 'tracks', 'segments', 'refit_outcomes'):
            if not (old/(name+'.npz')).exists() and not (new/(name+'.npz')).exists():
                continue
            with np.load(old/(name+'.npz')) as a, np.load(new/(name+'.npz')) as b:
                assert a.files == b.files
                for field in a.files:
                    np.testing.assert_array_equal(a[field], b[field])
                    fields += 1
        run = json.loads((new/'run.json').read_text())
        assert run['success'] and run['isolate_smoothing'] and run['ckf_max_step_mm'] is None
        assert run['actual_ckf_max_step_mm'] == 0.
        checks.append(dict(dsid=point['dsid'], variant=point['variant'], audit_rows_exact=count,
            smoothing_rows_exact=smoothing, array_fields_exact=fields,
            loaded_libraries=run['loaded_acts_libraries'], runtime=validate(new)))
    output.mkdir(exist_ok=False)
    write_json(output/'validation.json', dict(default_replay_exact=True, checks=checks))


def state_difference(a, b, variant_a, variant_b, stage):
    _, ca, _ = common(a, variant_a, stage)
    _, cb, _ = common(b, variant_b, stage)
    da, _ = direction_f(a[stage], variant_a)
    db, _ = direction_f(b[stage], variant_b)
    return dict(**matrix_metrics(ca, cb), direction_difference_norm=float(np.linalg.norm(da-db)),
        qop_difference_GeVinv=float(abs(a[stage][4]-b[stage][4])))


def compare_path(a, b, va, vb, key, tag, rows):
    ia, ib = first(a), first(b)
    counts = Counter(inputs=1)
    counts['before_selected'] += ia is not None
    counts['after_selected'] += ib is not None
    counts['one_sided_selected'] += (ia is None) != (ib is None)
    for label, item, index in [('before', a, ia), ('after', b, ib)]:
        if index is None:
            continue
        selected_sign = sign(qop(item['selected'][index], 'selected'))
        counts[label+'_selected_wrong'] += selected_sign == -item['truth_charge']
        counts[label+'_selected_invalid'] += selected_sign == 0
        final = item['final'].get(index)
        if final:
            counts[label+'_finals'] += 1
            final_sign = sign(qop(final, 'final'))
            counts[label+'_final_wrong'] += final_sign == -item['truth_charge']
            counts[label+'_final_invalid'] += final_sign == 0
            counts[label+'_refit_success'] += final['refit'] == 1
    if ia is None or ib is None:
        return counts
    counts['both_selected'] += 1
    counts['selected_sign_difference'] += sign(qop(a['selected'][ia], 'selected')) != sign(qop(b['selected'][ib], 'selected'))
    fa, fb = a['final'].get(ia), b['final'].get(ib)
    if fa and fb:
        counts['both_final'] += 1
        counts['final_sign_difference'] += sign(qop(fa, 'final')) != sign(qop(fb, 'final'))
    if a['selected'][ia]['hits'] != b['selected'][ib]['hits']:
        counts['different_ordered_hits'] += 1
        return counts
    counts['same_ordered_hits'] += 1
    sa, sb = forward_states(a, ia), forward_states(b, ib)
    if list(map(signature, sa)) != list(map(signature, sb)):
        counts['different_full_state_signatures'] += 1
        return counts
    counts['same_full_state_signatures'] += 1
    counts['terminal_only'] += len(sa) == 1
    for ordinal, (ra, rb) in enumerate(zip(sa, sb)):
        item = dict(tag, entry=key[0], component=key[1], forward_index=ordinal,
            geometry_id=ra['geometry_id'], path_length_mm=ra['path_length_mm'],
            measurement_flag=ra['measurement_flag'])
        for stage in ('predicted', 'filtered', 'smoothed'):
            item[stage] = state_difference(ra, rb, va, vb, stage)
        rows.append(item)
        counts['paired_states'] += 1
    return counts


def analyze(root, output):
    assert json.loads((root/'default_gate_v1/validation.json').read_text())['default_replay_exact']
    manifest = json.loads((root/'manifest.json').read_text())
    points = {(p['dsid'], p['variant']): p for p in manifest['points']}
    groups, all_rows, runtime_checks = {}, [], []
    for dsid in SAMPLES:
        data = {}
        for variant in ('beforeGeo', 'afterGeo'):
            point = points[dsid, variant]
            for mode in MODES:
                folder = Path(point[mode])
                assert (folder/'exit_code.txt').read_text().strip() == '0'
                run = json.loads((folder/'run.json').read_text())
                assert run['success'] and run['isolate_smoothing']
                assert run['ckf_max_step_mm'] == dict(default=None, step100mm=100., step10mm=10.)[mode]
                assert run['actual_ckf_max_step_mm'] == dict(default=0., step100mm=100., step10mm=10.)[mode]
                for library in run['loaded_acts_libraries']:
                    assert Path(library).resolve() == (Path(run['library'])/Path(library).name).resolve()
                runtime_checks.append(dict(mode=mode, **validate(folder)))
                data[variant, mode] = components(folder)
        reference_keys = data['beforeGeo', 'default'].keys()
        assert all(item.keys() == reference_keys for item in data.values())
        for (variant, mode), records in data.items():
            assert all(records[k]['seed'] == data[variant, 'default'][k]['seed'] for k in records)
        comparisons = [(mode, 'between_geometry', ('beforeGeo', mode), ('afterGeo', mode)) for mode in MODES]
        comparisons += [(mode, variant+'_step_change', (variant, 'default'), (variant, mode))
                        for variant in ('beforeGeo', 'afterGeo') for mode in MODES[1:]]
        for mode, kind, a_key, b_key in comparisons:
            rows, counts = [], Counter()
            for key in reference_keys:
                counts.update(compare_path(data[a_key][key], data[b_key][key], a_key[0], b_key[0], key,
                    dict(dsid=dsid, mode=mode, comparison=kind), rows))
            fields = ('direction_scaled_max', 'direction_absolute_max_rad2',
                      'direction_difference_norm', 'qop_difference_GeVinv')
            groups[dsid+'|'+kind+'|'+mode] = dict(counts=counts,
                stages={stage:{field:metrics([r[stage][field] for r in rows]) for field in fields}
                        for stage in ('predicted', 'filtered', 'smoothed')})
            all_rows.extend(rows)
    output.mkdir(exist_ok=False)
    with gzip.open(output/'paired_states.jsonl.gz', 'wt', compresslevel=4) as stream:
        for row in all_rows:
            stream.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
    write_json(output/'summary.json', dict(groups=groups, source_events=60,
        runtime_checks=runtime_checks,
        scope='Complete native CKF/conversion/refit under explicit CKF step limits. Same physical component inputs; '
              'first-ranked selected candidate per component. Selection, hit changes, invalid signs, final sign '
              'and refit failures retained. All native candidate totals and errors are in runtime_checks. '
              'No finite-mixture posterior aggregation.'))
    write_json(output/'validation.json', dict(default_replay_exact=True, complete_seed_inputs_per_variant_exact=True,
        input_keys_exact=True, paired_states=len(all_rows), comparisons=len(groups),
        runtime_checks=runtime_checks))
    print(json.dumps(dict(comparisons=len(groups), paired_states=len(all_rows))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--phase', choices=['gate', 'analyze'], required=True)
    args = parser.parse_args()
    dict(gate=gate, analyze=analyze)[args.phase](args.input, args.output)


if __name__ == '__main__':
    main()
