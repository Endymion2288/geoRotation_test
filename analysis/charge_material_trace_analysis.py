#!/usr/bin/env python3
"""Audit native material records against the unchanged Calypso endpoints."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from slope_tables import write_json


INTEGER = {'source_entry', 'seed', 'candidate', 'target', 'interaction', 'geometry_id',
           'surface', 'updated_volume_step', 'steps', 'interactions'}


def read(path):
    with path.open() as stream:
        return [{k: int(v) if k in INTEGER else float(v) for k, v in row.items()}
                for row in csv.DictReader(stream)]


def key(row):
    return tuple(row[k] for k in ['source_entry', 'candidate', 'target'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = args.root/'physical_prior_material_trace_v3'
    manifest = json.loads((source/'manifest.json').read_text())
    validation, cases = {}, {}
    for dsid in manifest['samples']:
        calls, interactions = {}, {}
        for variant in ['beforeGeo', 'afterGeo']:
            folder = source/dsid/variant
            assert (folder/'exit_code.txt').read_text().strip() == '0'
            run = json.loads((folder/'run.json').read_text())
            assert run['success'] and run['material_trace']
            with np.load(folder/'transport.npz') as file:
                current = dict(file)
            with np.load(args.root/'physical_prior_transport_v2'/dsid/'material_on'/variant/'transport.npz') as old:
                index = np.searchsorted(old['candidate'], current['candidate'])
                assert np.array_equal(old['candidate'][index], current['candidate'])
                for field in old.files:
                    np.testing.assert_array_equal(old[field][index], current[field], err_msg=dsid+'/'+variant+'/'+field)
                fields = len(old.files)
            summary = read(folder/'material_trace.csv.calls.csv')
            calls[variant] = {key(row): row for row in summary}
            assert len(calls[variant]) == len(summary) == 3*len(current['candidate'])
            rows = read(folder/'material_trace.csv')
            grouped = {}
            for row in rows:
                grouped.setdefault(key(row), []).append(row)
            interactions[variant] = grouped
            closure = []
            for row in summary:
                i = int(np.flatnonzero(current['candidate'] == row['candidate'])[0])
                expected = np.linalg.norm(current['projection_F'][i, row['target'], 3:])
                np.testing.assert_allclose(row['normal_end_p_GeV'], expected, rtol=1e-14, atol=1e-12)
                closure.append(row['sum_deltaP_GeV']-(row['normal_end_p_GeV']-row['input_p_GeV']))
            validation[dsid+'|'+variant] = dict(original_v2_fields_exact=fields,
                normal_endpoint_matches_original=True, calls=len(summary), interactions=len(rows),
                record_overload_energy_closure_max_abs_GeV=float(np.max(np.abs(closure))),
                interpretation='Material records use the native propagationSteps overload with SteppingLogger; endpoint p remains from native propagate. Their measured closure is reported, not assumed exact.')
        assert set(calls['beforeGeo']) == set(calls['afterGeo'])
        for identity in sorted(calls['beforeGeo']):
            left, right = [calls[v][identity] for v in ['beforeGeo', 'afterGeo']]
            a, b = [interactions[v].get(identity, []) for v in ['beforeGeo', 'afterGeo']]
            signature = lambda rows: [(r['geometry_id'], r['surface']) for r in rows]
            sequence = signature(a) == signature(b)
            differences = []
            if sequence:
                for before, after in zip(a, b):
                    material_fields = ['thickness_mm']+[f'material_{i}' for i in range(5)]
                    scale = max(abs(after[k]-before[k])/max(abs(before[k]), 1e-20) for k in material_fields)
                    differences.append(dict(interaction=before['interaction'], geometry_id=before['geometry_id'],
                        position_difference_mm=float(np.linalg.norm([after[k]-before[k] for k in ['xF_mm', 'yF_mm', 'zF_mm']])),
                        path_correction_difference=after['path_correction']-before['path_correction'],
                        max_material_relative_difference=scale,
                        before=before, after=after))
            case = dict(dsid=dsid, source_entry=identity[0], candidate=identity[1], target=identity[2],
                identical_geometry_sequence=sequence, before_call=left, after_call=right,
                native_endpoint_relative_p=right['normal_end_p_GeV']/left['normal_end_p_GeV']-1.,
                material_X0_difference=right['material_X0']-left['material_X0'],
                material_rows=differences)
            changed = [row for row in differences if row['max_material_relative_difference'] > 1e-5]
            case['first_material_difference_above_1e_minus_5'] = changed[0] if changed else None
            cases[dsid+'|'+'|'.join(map(str, identity))] = case
    args.output.mkdir(exist_ok=False)
    values = list(cases.values())
    write_json(args.output/'summary.json', dict(validation=validation, cases=cases,
        paired_calls=len(cases), selected_physical_points=sum(s['input_points'] for s in manifest['samples'].values()),
        identical_geometry_sequences=sum(r['identical_geometry_sequence'] for r in values),
        calls_with_material_value_difference=sum(r['first_material_difference_above_1e_minus_5'] is not None for r in values),
        max_paired_material_position_difference_mm=max((r['position_difference_mm'] for v in values for r in v['material_rows']), default=0.),
        source_manifest=str(source/'manifest.json'),
        scope='Selected extrema plus fixed controls, not an unbiased incidence estimate. The unchanged native endpoint reproduces v2 exactly; material records diagnose actual slab values and surface sequence.'))
    print('Paired calls', len(cases), '; identical geometry sequences', sum(r['identical_geometry_sequence'] for r in values),
          '; calls with differing material values', sum(r['first_material_difference_above_1e_minus_5'] is not None for r in values))


if __name__ == '__main__':
    main()
