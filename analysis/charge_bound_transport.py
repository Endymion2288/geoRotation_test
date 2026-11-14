#!/usr/bin/env python3
"""Compare native Calypso full-C transport with local physical input closure."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from charge_component_refit import metrics
from charge_forward_covariance import common, matrix_metrics
from bound_surface_coordinates import direction_f, surface_position, ROTATE_F_TO_A
from slope_tables import write_json


def loaded_libraries(root, variant, run, revision):
    overlay = root/('overlay_'+variant+('_v2' if revision == 2 else ''))/'lib'
    expected = {'libCalypsoSlopeDiagnostics.so', 'libFaserActsGeometry.so', 'libFaserActsGeometryLib.so'}
    assert {Path(p).name for p in run['loaded_acts_libraries']} == expected
    assert all(Path(p).resolve() == (overlay/Path(p).name).resolve() for p in run['loaded_acts_libraries'])


def gate(root, output, revision=1):
    checks = []
    for variant in ('beforeGeo', 'afterGeo'):
        folder = root/('legacy_replay_v2' if revision == 2 else 'legacy_replay')/variant
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        task = json.loads((folder/'source_task.json').read_text())
        reference = Path(task['reference'])
        with np.load(reference/'transport.npz') as old, np.load(folder/'transport.npz') as new:
            assert set(old.files).issubset(new.files)
            for field in old.files:
                np.testing.assert_array_equal(old[field], new[field], err_msg=variant+':'+field)
            fields, rows = len(old.files), len(old['source_entry'])
        for name in ('material_trace.csv', 'material_trace.csv.calls.csv'):
            assert (reference/name).read_bytes() == (folder/name).read_bytes()
        run = json.loads((folder/'run.json').read_text())
        assert run['success'] and run['max_step_size_m'] == 10.
        loaded_libraries(root, variant, run, revision)
        checks.append(dict(variant=variant, reference=str(reference), new=str(folder),
            original_array_fields_exact=fields, original_inputs=rows, material_csv_bytes_exact=True,
            loaded_libraries=run['loaded_acts_libraries']))
    output.mkdir(exist_ok=False)
    write_json(output/'validation.json', dict(all_legacy_fields_and_material_calls_exact=True,
        checks=checks, scope='No-C legacy inputs; full-C inputs are a separately validated new opt-in mode.'))


def load(folder):
    assert (folder/'exit_code.txt').read_text().strip() == '0', folder
    task = json.loads((folder/'source_task.json').read_text())
    run = json.loads((folder/'run.json').read_text())
    assert run['success']
    with np.load(folder/'transport.npz') as archive:
        arrays = {name: archive[name] for name in archive.files}
    rows = {}
    for i in range(len(arrays['source_entry'])):
        key = int(arrays['source_entry'][i]), int(arrays['candidate'][i])
        assert key not in rows
        rows[key] = {name: value[i] for name, value in arrays.items()}
    inputs = {}
    for line in Path(task['states']).read_text().splitlines():
        fields = line.split()
        assert len(fields) == 52
        entry, seed, candidate = map(int, fields[:3])
        assert int(fields[6]) == 2
        mean = np.asarray(list(map(float, fields[9:15])))
        covariance = np.asarray(list(map(float, fields[15:51]))).reshape(6, 6)
        row = rows[entry, candidate]
        assert row['input_has_covariance'] == 1 and row['seed'] == seed
        np.testing.assert_array_equal(row['input_covariance_native'], covariance)
        np.testing.assert_allclose(row['input_angles_native'], mean[2:4], rtol=0, atol=2e-15)
        assert row['input_qop_GeVinv'] == mean[4]
        inputs[entry, candidate] = dict(bound=mean, cov=covariance,
            surface=row['input_surface_native'], surface_type=int(row['input_surface_type']),
            position_F_mm=row['input_position_F_mm'], momentum_F_GeV=row['input_momentum_F_GeV'])
    assert rows.keys() == inputs.keys()
    return rows, inputs, run


def output_state(row):
    return dict(bound=row['projection_bound_native'][0],
        cov=row['projection_covariance_native'][0], surface=row['projection_surface_native'][0],
        surface_type=int(row['projection_surface_type'][0]),
        position_F_mm=row['projection_F'][0, :3], momentum_F_GeV=row['projection_F'][0, 3:])


def native_checks(row, variant):
    """Check typed local-to-global coordinates against native Calypso fields."""
    position, _ = surface_position(row['bound'], row['surface'], row['surface_type'])
    if variant == 'afterGeo':
        position = ROTATE_F_TO_A.T @ position
    np.testing.assert_allclose(position, row['position_F_mm'], rtol=0, atol=2e-9)
    direction, _ = direction_f(row['bound'], variant)
    momentum = row['momentum_F_GeV']
    np.testing.assert_allclose(direction, momentum/np.linalg.norm(momentum), rtol=0, atol=3e-15)
    return direction


def difference(a, b):
    assert a['surface_type'] == b['surface_type']
    da, db = native_checks(a, 'beforeGeo'), native_checks(b, 'afterGeo')
    _, common_a, _ = common(a, 'beforeGeo')
    _, common_b, _ = common(b, 'afterGeo')
    return dict(**matrix_metrics(common_a, common_b),
        position_difference_mm=float(np.linalg.norm(a['position_F_mm']-b['position_F_mm'])),
        direction_difference_norm=float(np.linalg.norm(da-db)),
        qop_difference_GeVinv=float(abs(a['bound'][4]-b['bound'][4])),
        before_qop_GeVinv=float(a['bound'][4]), after_qop_GeVinv=float(b['bound'][4]))


def analyze(root, output, revision=1):
    assert json.loads((root/('legacy_gate_v3' if revision == 2 else 'legacy_gate_v1')/'validation.json').read_text())['all_legacy_fields_and_material_calls_exact']
    manifest = json.loads((root/'manifest.json').read_text())
    cases = json.loads((root/'cases.json').read_text())
    groups, all_rows, checks, surfaces = {}, [], [], {}
    rotate = np.eye(4)
    rotate[:3, :3] = ROTATE_F_TO_A
    for dsid in ('100044', '100116', '100117'):
        selected_cases = [c for c in cases if c['dsid'] == dsid]
        for control in manifest['controls']:
            data, inputs, runs = {}, {}, {}
            for variant in ('beforeGeo', 'afterGeo'):
                folder = root/('runs_v2' if revision == 2 else 'runs')/dsid/variant/control
                data[variant], inputs[variant], runs[variant] = load(folder)
                assert len(data[variant]) == len(selected_cases)
                loaded_libraries(root, variant, runs[variant], revision)
            grouped = defaultdict(list)
            counts = Counter()
            for case in selected_cases:
                key = case['entry'], case['candidate']
                a, b = data['beforeGeo'][key], data['afterGeo'][key]
                assert a['geometry_id'] == b['geometry_id'] == case['source_geometry_id']
                ia, ib = inputs['beforeGeo'][key], inputs['afterGeo'][key]
                for v, state in [('beforeGeo', ia), ('afterGeo', ib)]:
                    catalog_key = v+'|'+str(case['source_geometry_id'])
                    metadata = dict(variant=v, geometry_id=case['source_geometry_id'],
                        surface_type=state['surface_type'], surface=state['surface'].tolist())
                    assert catalog_key not in surfaces or surfaces[catalog_key] == metadata
                    surfaces[catalog_key] = metadata
                surface_delta = np.abs(rotate @ ia['surface']-ib['surface']).max()
                assert surface_delta < 2e-9
                item = dict(case, control=control, input=difference(ia, ib),
                    source_surface_type=ia['surface_type'],
                    before_valid=int(a['projection_valid'][0]), after_valid=int(b['projection_valid'][0]))
                if case['mode'] == 'fixed_before_input':
                    assert item['input']['position_difference_mm'] < 2e-9
                    assert item['input']['direction_difference_norm'] < 3e-15
                    assert item['input']['qop_difference_GeVinv'] == 0.
                    assert item['input']['all_scaled_max'] < 5e-10
                    assert item['input']['direction_scaled_max'] < 5e-10
                counts[case['mode']+'_inputs'] += 1
                counts[case['mode']+'_both_valid'] += item['before_valid'] == item['after_valid'] == 1
                counts[case['mode']+'_before_failed'] += item['before_valid'] != 1
                counts[case['mode']+'_after_failed'] += item['after_valid'] != 1
                if item['before_valid'] == item['after_valid'] == 1:
                    sa, sb = output_state(a), output_state(b)
                    assert np.abs(rotate @ sa['surface']-sb['surface']).max() < 2e-9
                    item['output'] = difference(sa, sb)
                    item['target_surface_type'] = sa['surface_type']
                    for v, state in [('beforeGeo', sa), ('afterGeo', sb)]:
                        catalog_key = v+'|'+str(case['target_geometry_id'])
                        metadata = dict(variant=v, geometry_id=case['target_geometry_id'],
                            surface_type=state['surface_type'], surface=state['surface'].tolist())
                        assert catalog_key not in surfaces or surfaces[catalog_key] == metadata
                        surfaces[catalog_key] = metadata
                    grouped[case['mode']].append(item['output'])
                all_rows.append(item)
            groups[dsid+'|'+control] = dict(counts=counts,
                modes={mode: {field: metrics([r[field] for r in rows]) for field in
                    ('direction_scaled_max', 'direction_absolute_max_rad2', 'position_difference_mm',
                     'direction_difference_norm', 'qop_difference_GeVinv')} for mode, rows in grouped.items()})
            checks.append(dict(dsid=dsid, control=control, inputs_per_version=len(selected_cases),
                input_covariance_exactly_loaded=True, physical_fixed_input_checked=True,
                actual_loaded_libraries={v: runs[v]['loaded_acts_libraries'] for v in runs}))
    output.mkdir(exist_ok=False)
    assert len(all_rows) == len(cases)*len(manifest['controls'])
    write_json(output/'paired_cases.json', all_rows)
    write_json(output/'surface_catalog.json', surfaces)
    write_json(output/'summary.json', dict(groups=groups,
        scope='Actual Calypso bound-to-bound transport with full covariance. Surface material splitting can '
              'differ from an internal CKF segment; this is a controlled local propagation experiment, '
              'not exact CKF replay or finite-width posterior equivalence. No failed case removed.'))
    write_json(output/'validation.json', dict(checks=checks, legacy_replay_exact=True,
        complete_mean_covariance_surface_input_verified=True, paired_cases=len(all_rows),
        failures_retained=True, native_output_representations_cross_checked=True))
    print(json.dumps(dict(paired_cases=len(all_rows), groups=len(groups))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--phase', choices=['gate', 'analyze'], required=True)
    parser.add_argument('--revision', type=int, choices=[1, 2], default=1)
    args = parser.parse_args()
    dict(gate=gate, analyze=analyze)[args.phase](args.input, args.output, args.revision)


if __name__ == '__main__':
    main()
