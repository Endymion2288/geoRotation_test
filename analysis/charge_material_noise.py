#!/usr/bin/env python3
"""Audit material variances recorded by the actual Calypso extrapolator."""
import argparse
from collections import Counter, defaultdict
import csv
from itertools import zip_longest
import json
from pathlib import Path

import numpy as np

from bound_surface_coordinates import ROTATE_F_TO_A
from charge_component_refit import metrics
from slope_tables import write_json


FIELDS = ('native_sigma_phi2_rad2', 'native_sigma_theta2_rad2', 'native_sigma_qop2_GeVinv2')


def trace_gate(old, new, enabled):
    names = ('material_trace.csv', 'material_trace.csv.calls.csv')
    if not enabled:
        assert all(not (folder/name).exists() for folder in (old, new) for name in names)
        return dict(trace_recorded=False, original_trace_rows_exact=0, native_calls_exact=None)
    count = 0
    with (old/names[0]).open() as a, (new/names[0]).open() as b:
        left, right = csv.DictReader(a), csv.DictReader(b)
        assert right.fieldnames == left.fieldnames+list(FIELDS)
        for ra, rb in zip_longest(left, right):
            assert ra is not None and rb is not None
            assert ra == {name: rb[name] for name in left.fieldnames}
            assert all(np.isfinite(float(rb[name])) and float(rb[name]) >= 0. for name in FIELDS)
            count += 1
    assert (old/names[1]).read_bytes() == (new/names[1]).read_bytes()
    return dict(trace_recorded=True, original_trace_rows_exact=count, native_calls_exact=True)


def gate(root, output):
    manifest = json.loads((root/'manifest.json').read_text())
    checks = []
    for point in manifest['points']:
        old, new = Path(point['reference']), Path(point['output'])
        assert (new/'exit_code.txt').read_text().strip() == '0'
        with np.load(old/'transport.npz') as a, np.load(new/'transport.npz') as b:
            assert a.files == b.files
            for field in a.files:
                np.testing.assert_array_equal(a[field], b[field])
            fields, inputs = len(a.files), len(a['source_entry'])
        a, b = [json.loads((p/'source_task.json').read_text()) for p in (old, new)]
        assert all(a[key] == b[key] for key in ('input', 'states', 'material', 'variant', 'options'))
        trace = trace_gate(old, new, '--trace-material' in a['options'])
        run = json.loads((new/'run.json').read_text())
        assert run['success']
        overlay = root/('overlay_'+point['variant'])/'lib'
        assert {Path(p).name for p in run['loaded_acts_libraries']} == {
            'libCalypsoSlopeDiagnostics.so', 'libFaserActsGeometry.so', 'libFaserActsGeometryLib.so'}
        assert all(Path(p).resolve() == (overlay/Path(p).name).resolve() for p in run['loaded_acts_libraries'])
        checks.append(dict(point, transport_array_fields_exact=fields, inputs=inputs,
            **trace, loaded_libraries=run['loaded_acts_libraries']))
    assert (root/'software_probe/CalypsoSlopeDiagnosticsAlg.cxx').read_bytes() == (root/'CalypsoSlopeDiagnosticsAlg.cxx').read_bytes()
    output.mkdir(exist_ok=False)
    write_json(output/'validation.json', dict(all_transport_arrays_and_original_trace_columns_exact=True,
        source_compiled_copy_exact=True, recorded_trace_jobs=sum(c['trace_recorded'] for c in checks),
        unrecorded_trace_jobs=sum(not c['trace_recorded'] for c in checks), checks=checks))


def trace_rows(folder):
    if not (folder/'material_trace.csv').exists():
        return None
    groups = defaultdict(list)
    with (folder/'material_trace.csv').open() as stream:
        for row in csv.DictReader(stream):
            assert int(row['target']) == 0
            groups[int(row['source_entry']), int(row['candidate'])].append(row)
    return groups


def noise(row, variant):
    direction_f = np.asarray([float(row['dirF'+axis]) for axis in 'xyz'])
    np.testing.assert_allclose(direction_f @ direction_f, 1., rtol=0., atol=5e-15)
    rotation = np.eye(3) if variant == 'beforeGeo' else ROTATE_F_TO_A
    native = rotation @ direction_f
    native /= np.linalg.norm(native)
    x, y, z = native
    radius = np.hypot(x, y)
    assert radius > 0.
    derivative = np.array([[-y, x*z/radius], [x, y*z/radius], [0., -radius]])
    variances = np.asarray([float(row[name]) for name in FIELDS])
    physical = rotation.T @ derivative @ np.diag(variances[:2]) @ derivative.T @ rotation
    result = dict(direction_F=direction_f.tolist(), variance_phi_rad2=float(variances[0]),
        variance_theta_rad2=float(variances[1]), variance_qop_GeVinv2=float(variances[2]),
        physical_direction_Q=physical.tolist(), scattering=bool(variances[1] > 0.))
    if not result['scattering']:
        result.update(isotropy_error=None, float_square_relation_error=None,
            double_square_relation_error=None, expected_isotropy_from_float_square=None,
            same_double_square_direction_Q=physical.tolist())
        return result
    # Recover the native float theta0 from its recorded phi variance. This
    # checks a runtime arithmetic relation, without recalculating scattering.
    native_square = variances[0]*radius**2
    theta0 = np.float32(np.sqrt(native_square))
    square_double = float(theta0)**2
    square_float = float(np.float32(theta0*theta0))
    corrected = rotation.T @ derivative @ np.diag([variances[0], square_double]) @ derivative.T @ rotation
    result.update(isotropy_error=float(native_square/variances[1]-1.),
        float_square_relation_error=float((variances[1]-square_float)/variances[1]),
        double_square_relation_error=float((native_square-square_double)/square_double),
        expected_isotropy_from_float_square=float(square_double/square_float-1.),
        same_double_square_direction_Q=corrected.tolist())
    return result


def noise_difference(a, b):
    qa, qb = [np.asarray(r['physical_direction_Q']) for r in (a, b)]
    ca, cb = [np.asarray(r['same_double_square_direction_Q']) for r in (a, b)]
    scale = max((a['variance_theta_rad2']+b['variance_theta_rad2'])/2., 1e-300)
    return dict(physical_direction_Q_scaled_max=float(np.abs(qa-qb).max()/scale),
        physical_direction_Q_absolute_max=float(np.abs(qa-qb).max()),
        same_double_square_Q_scaled_max=float(np.abs(ca-cb).max()/scale),
        direction_difference_norm=float(np.linalg.norm(np.asarray(a['direction_F'])-b['direction_F'])),
        theta_variance_relative_difference=float(abs(a['variance_theta_rad2']-b['variance_theta_rad2'])/scale),
        qop_variance_absolute_difference_GeVinv2=abs(a['variance_qop_GeVinv2']-b['variance_qop_GeVinv2']))


def analyze(root, output, gate_validation=None):
    gate_validation = gate_validation or root/'legacy_gate_v1/validation.json'
    assert json.loads(gate_validation.read_text())['all_transport_arrays_and_original_trace_columns_exact']
    manifest = json.loads((root/'manifest.json').read_text())
    cases = json.loads(Path(manifest['cases']).read_text())
    paths = {(Path(p['output']).parts[-3], p['variant'], Path(p['output']).name): Path(p['output'])
             for p in manifest['points']}
    groups, paired, sequences = {}, [], []
    controls = sorted({k[2] for k in paths})
    for dsid in ('100044', '100116', '100117'):
        for control in controls:
            data = {v:trace_rows(paths[dsid, v, control]) for v in ('beforeGeo', 'afterGeo')}
            for mode in ('native_inputs', 'fixed_before_input'):
                counts, rows = Counter(), []
                for case in (c for c in cases if c['dsid'] == dsid and c['mode'] == mode):
                    key = case['entry'], case['candidate']
                    counts['input_links'] += 1
                    if any(data[v] is None for v in data):
                        assert all(data[v] is None for v in data)
                        counts['unrecorded_links'] += 1
                        sequences.append(dict(case, control=control, noise_recorded=False,
                            sequence_equal=None, before_interactions=None, after_interactions=None))
                        continue
                    a, b = data['beforeGeo'].get(key, []), data['afterGeo'].get(key, [])
                    signature = lambda r: tuple(r[k] for k in ('interaction', 'geometry_id', 'surface', 'updated_volume_step'))
                    equal = list(map(signature, a)) == list(map(signature, b))
                    counts['same_trace_sequence'] += equal
                    counts['different_trace_sequence'] += not equal
                    sequences.append(dict(case, control=control, noise_recorded=True, sequence_equal=equal,
                        before_interactions=len(a), after_interactions=len(b)))
                    if not equal:
                        continue
                    for ra, rb in zip(a, b):
                        na, nb = noise(ra, 'beforeGeo'), noise(rb, 'afterGeo')
                        counts['paired_interactions'] += 1
                        counts['both_scattering'] += na['scattering'] and nb['scattering']
                        counts['one_sided_scattering'] += na['scattering'] != nb['scattering']
                        row = dict(case, control=control, interaction=int(ra['interaction']),
                            geometry_id=int(ra['geometry_id']), surface=bool(int(ra['surface'])),
                            before=na, after=nb, difference=noise_difference(na, nb))
                        rows.append(row)
                groups[dsid+'|'+control+'|'+mode] = dict(counts=counts,
                    differences={field:metrics([r['difference'][field] for r in rows]) for field in
                        ('physical_direction_Q_scaled_max', 'same_double_square_Q_scaled_max',
                         'direction_difference_norm', 'theta_variance_relative_difference',
                         'qop_variance_absolute_difference_GeVinv2')},
                    native_relations={v:{field:metrics([r[v][field] for r in rows if r[v]['scattering']])
                        for field in ('isotropy_error', 'float_square_relation_error', 'double_square_relation_error')}
                        for v in ('before', 'after')})
                paired.extend(rows)
    output.mkdir(exist_ok=False)
    import gzip
    with gzip.open(output/'paired_interactions.jsonl.gz', 'wt', compresslevel=4) as stream:
        for row in paired:
            stream.write(json.dumps(row, separators=(',', ':'), allow_nan=False)+'\n')
    write_json(output/'trace_sequences.json', sequences)
    write_json(output/'summary.json', dict(groups=groups,
        scope='Native RecordedMaterial noise variances in Calypso propagationSteps. '
              'The same-double-square Q is an analytic arithmetic diagnostic using these recorded terms; '
              'it is not a rerun or a corrected track covariance. Independent links restart from frozen states. '
              'Full bound transport arrays and original trace columns reproduce the reference exactly.'))
    write_json(output/'validation.json', dict(transport_and_trace_gate_passed=True,
        paired_interactions=len(paired), complete_trace_sequence_counts=len(sequences),
        different_trace_sequences=sum(r['sequence_equal'] is False for r in sequences),
        unrecorded_links=sum(not r['noise_recorded'] for r in sequences),
        failed_or_mismatched_sequences_retained=True))
    print(json.dumps(dict(paired_interactions=len(paired), links=len(sequences))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--phase', choices=['gate', 'analyze'], required=True)
    parser.add_argument('--gate-validation', type=Path)
    args = parser.parse_args()
    if args.phase == 'gate':
        gate(args.input.resolve(), args.output.resolve())
    else:
        analyze(args.input.resolve(), args.output.resolve(), args.gate_validation)


if __name__ == '__main__':
    main()
