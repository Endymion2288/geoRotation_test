#!/usr/bin/env python3
"""Validate native Calypso bin queries and the corrected finite-prior transport."""
import argparse
import csv
from functools import lru_cache
import json
from pathlib import Path

import numpy as np

from charge_prior_transport_analysis import compare
from slope_tables import write_json


def arrays(path):
    assert (path/'exit_code.txt').read_text().strip() == '0'
    with np.load(path/'transport.npz') as file:
        order = np.argsort(file['candidate'])
        return {key: file[key][order] for key in file.files}


@lru_cache(maxsize=None)
def surface_materials(path):
    entries = json.loads(path.read_text())['Surfaces']['entries']
    return {((e.get('volume', 0) << 56) | (e.get('boundary', 0) << 48) |
             (e.get('layer', 0) << 36) | (e.get('approach', 0) << 28) |
             (e.get('sensitive', 0) << 8)): e['value']['material'] for e in entries}


def query_table(path, materials):
    with (path/'material_trace.csv').open() as file:
        rows = list(csv.DictReader(file))
    result = {}
    for row in rows:
        key = tuple(int(row[name]) for name in
            ['source_entry', 'seed', 'candidate', 'target', 'interaction'])
        assert key not in result
        material = materials[int(row['geometry_id'])]
        assert int(row['binned']) == 1
        encoded = material['binUtility'].get('transform', {})
        # Independent ACTS from_json rule, not the Python rotation helper.
        expected = np.eye(4)
        if encoded.get('rotation') is not None:
            expected[:3, :3] = np.asarray(encoded['rotation']).reshape(3, 3)
        expected[:3, 3] = encoded.get('translation') or [0., 0., 0.]
        actual = np.array([[float(row[f'bin_transform_{r}{c}']) for c in range(4)] for r in range(3)])
        np.testing.assert_array_equal(actual, expected[:3])
        bins = [int(row[f'bin_{i}']) for i in range(3)]
        cell = material['data'][bins[1]][bins[0]]
        raw = np.array([float(row['raw_thickness_mm'])]+
                       [float(row[f'raw_material_{i}']) for i in range(5)])
        np.testing.assert_array_equal(raw, [cell['thickness']]+cell['material'])
        result[key] = dict(geometry_id=int(row['geometry_id']), bins=bins,
            local_mm=[float(row['bin_local_'+axis]) for axis in ['x', 'y', 'z']],
            position_F_mm=[float(row[axis+'F_mm']) for axis in ['x', 'y', 'z']],
            native_transform=actual.tolist(), raw_material=raw.tolist())
    return result


def paired_queries(before, after):
    assert before.keys() == after.keys()
    changes, sequence_changes, local_error, position_error = 0, 0, 0., 0.
    examples = []
    for key in before:
        a, b = before[key], after[key]
        sequence_changes += a['geometry_id'] != b['geometry_id']
        changes += a['raw_material'] != b['raw_material']
        local_error = max(local_error, np.max(np.abs(np.array(a['local_mm'])-b['local_mm'])))
        position_error = max(position_error, np.max(np.abs(np.array(a['position_F_mm'])-b['position_F_mm'])))
        if len(examples) < 3 and a['raw_material'] != b['raw_material']:
            examples.append(dict(key=key, before=a, after=b))
    return dict(interactions=len(before), changed_geometry_id=sequence_changes,
        changed_raw_material=changes, changed_bin_indices=sum(before[k]['bins'] != after[k]['bins'] for k in before),
        maximum_local_coordinate_difference_mm=float(local_error),
        maximum_physical_position_difference_mm=float(position_error), examples=examples)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    source = root/'physical_prior_material_bins_v4'
    manifest = json.loads((source/'manifest.json').read_text())
    previous = root/'physical_prior_material_trace_v3'
    prior = root/'physical_prior_transport_v2'
    old_manifest = json.loads((prior/'manifest.json').read_text())
    cloud_root = Path(old_manifest['previous'])
    legacy = json.loads((root/'physical_prior_transport_analysis_v2/summary.json').read_text())
    validation, queries, groups = {}, {}, {}
    for dsid in ['100044', '100116', '100117']:
        tables = {}
        for case in ['beforeGeo_trace', 'afterGeo_trace', 'afterGeo_corrected_trace']:
            folder = source/dsid/case
            task = json.loads(Path(manifest['source_tasks'][dsid+'|'+case]).read_text())
            run = json.loads((folder/'run.json').read_text())
            assert run['success'] and run['args']['material'] == task['material']
            actual = [Path(p).resolve() for p in run['loaded_acts_libraries']]
            assert (Path(manifest['observer']).parent/'lib/libCalypsoSlopeDiagnostics.so').resolve() in actual
            for name in ['libFaserActsGeometry.so', 'libFaserActsGeometryLib.so']:
                assert (Path(manifest['libraries'][task['variant']])/name).resolve() in actual
            current = arrays(folder)
            if 'corrected' not in case:
                old = arrays(previous/dsid/task['variant'])
                for key in old:
                    np.testing.assert_array_equal(current[key], old[key], err_msg=dsid+'/'+case+'/'+key)
                validation[dsid+'|'+case] = dict(original_endpoint_fields_exact=len(old))
            tables[case] = query_table(folder, surface_materials(Path(task['material'])))
        queries[dsid+'|legacy'] = paired_queries(tables['beforeGeo_trace'], tables['afterGeo_trace'])
        queries[dsid+'|corrected'] = paired_queries(tables['beforeGeo_trace'], tables['afterGeo_corrected_trace'])
        with np.load(cloud_root/dsid/'physical_cloud.npz') as file:
            cloud = dict(file)
        before = arrays(prior/dsid/'material_on/beforeGeo')
        after_folder = source/dsid/'afterGeo_corrected_full'
        after = arrays(after_folder)
        run = json.loads((after_folder/'run.json').read_text())
        assert run['success'] and run['args']['material'] == manifest['corrected_map']
        assert run['interaction_energy_loss'] and run['interaction_scattering']
        actual = [Path(p).resolve() for p in run['loaded_acts_libraries']]
        for name in ['libCalypsoSlopeDiagnostics.so', 'libFaserActsGeometry.so', 'libFaserActsGeometryLib.so']:
            assert (Path(manifest['libraries']['afterGeo'])/name).resolve() in actual
        groups[dsid+'|row_major_map'] = compare(before, after, cloud)
        groups[dsid+'|legacy_map'] = legacy['groups'][dsid+'|material_on']
    args.output.mkdir(exist_ok=False)
    write_json(args.output/'summary.json', dict(validation=validation, queries=queries, groups=groups,
        source_manifest=str(source/'manifest.json'),
        plot_controls=[dict(key='legacy_map', label='Legacy rotated map', color='#c43c39'),
                       dict(key='row_major_map', label='Row-major rotated map', color='#205fb2')],
        plot_title='native Calypso | legacy/corrected material maps',
        scope='Frozen 44 extreme/control points for native bin queries; full 7680 physical prior points for corrected transport. '
              'The matrix and material-bin values are read from actual Calypso surfaces. '
              'Transport sign preservation does not imply complete nonlinear CKF posterior closure.'))
    print('Native material-bin validation and full corrected transport reduction complete')


if __name__ == '__main__':
    main()
