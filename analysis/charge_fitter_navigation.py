#!/usr/bin/env python3
"""Reduce real Calypso navigation logs and verify isolated component replays."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import re

import numpy as np

from bound_surface_coordinates import ROTATE_F_TO_A
from charge_component_refit import audit, components
from charge_forward_covariance import forward_states
from slope_tables import write_json
from validate_charge_components import validate


def component_rows(folder, key):
    entry, rows = -1, []
    for row in audit(folder/'prior_audit.jsonl.gz'):
        if row['stage'] == 'event':
            entry += 1
        if (entry, row.get('input_component', -1)) == key:
            rows.append(row)
    return rows


def log_records(path, variant):
    forward, backward = [], []
    prepared, geometry, markers, limits = [], [], [], []
    reversed_direction = False
    with gzip.open(path, 'rt') as stream:
        for number, line in enumerate(stream, 1):
            if 'KalmanFitter' not in line:
                continue
            if 'Reverse navigation direction after smoothing' in line:
                reversed_direction = True
            if 'Preparing ' in line and ' input measurements' in line:
                prepared.append(int(re.search(r'Preparing (\d+) input measurements', line)[1]))
            if 'KalmanFitter step at pos:' in line:
                pos, rest = line.split('step at pos:', 1)[1].split(' dir:', 1)
                direction, momentum = rest.split(' momentum:', 1)
                pos, direction = np.fromstring(pos, sep=' '), np.fromstring(direction, sep=' ')
                assert pos.shape == direction.shape == (3,)
                if variant == 'afterGeo':
                    pos, direction = ROTATE_F_TO_A.T@pos, ROTATE_F_TO_A.T@direction
                (backward if reversed_direction else forward).append([
                    *pos.tolist(), *direction.tolist(), float(momentum)])
            match = re.search(r'Measurement surface vol=(\d+)\|lay=(\d+)\|sen=(\d+) detected', line)
            if match:
                assert not reversed_direction
                volume, layer, sensitive = map(int, match.groups())
                geometry.append((volume << 56) | (layer << 36) | (sensitive << 8))
            match = re.search(r'step count limit of (\d+) \(did (\d+) steps\)', line)
            if match:
                limits.append(dict(configured_steps=int(match[1]), actual_steps=int(match[2])))
            if any(token in line for token in ('Last surface on layer', 'Reverse navigation direction',
                    'Smoothing successful', 'Target surface reached', 'ERROR')):
                markers.append(dict(line=number, text=line.rstrip()))
    assert len(prepared) == 1 and forward
    return dict(prepared_measurements=prepared[0], detected_geometry_ids=geometry,
        step_limits=limits, markers=markers, forward=np.asarray(forward),
        backward=np.asarray(backward).reshape(-1, 7))


def observation_shape(expected, actual):
    if actual == expected:
        return 'complete'
    if not actual:
        return 'empty'
    if actual == expected[-len(actual):]:
        return 'missing_prefix_only'
    if actual == expected[:len(actual)]:
        return 'missing_suffix_only'
    indices = [expected.index(hit) for hit in actual]
    if indices == sorted(set(indices)):
        return 'missing_both_ends' if indices == list(range(indices[0], indices[-1]+1)) else 'internal_gaps'
    return 'reordered_or_repeated'


def analyze(root, output):
    original = root.parent
    manifest = json.loads((root/'manifest.json').read_text())
    observations = json.loads((original/'observations.json').read_text())
    groups = defaultdict(Counter)
    for row in audit(original/'analysis_v1/solver_coverage_rows.jsonl.gz'):
        group = groups[row['dsid']+'|'+row['variant']]
        actual = row['actual_hits']['fixed_observation']
        group[observation_shape(row['expected_hits'], actual)] += 1
        terminal = row['native_terminal_records']['fixed_observation']
        if terminal['stage'] == 'component_finder_error':
            group[terminal['error_category']+': '+terminal['error_message']] += 1
    points, paths = [], {}
    for point in manifest['points']:
        folder = Path(point['output'])
        runtime = validate(folder)
        key = point['entry'], point['component']
        fresh, archived = component_rows(folder, key), component_rows(Path(point['reference']), key)
        assert fresh == archived, (point['case'], point['variant'], 'component output changed')
        run = json.loads((folder/'run.json').read_text())
        assert run['actual_ckf_acts_logging'] == run['ckf_acts_logging'] == 'VERBOSE'
        assert Path(run['library']).resolve() == (original/('overlay_'+point['variant'])/'lib').resolve()
        item = components(folder)[key]
        assert len(item['selected']) <= 1
        states = forward_states(item, next(iter(item['selected']))) if item['selected'] else []
        measurements = [r for r in states if r['measurement_flag']]
        record = log_records(folder/'tracking.log.gz', point['variant'])
        assert record['detected_geometry_ids'] == [r['geometry_id'] for r in measurements]
        expected = observations['100044']['events'][key[0]]['forward_ordered_hits']
        assert record['prepared_measurements'] == len(expected)
        name = point['case']+'_'+point['variant']
        paths[name+'_forward'] = record.pop('forward')
        paths[name+'_backward'] = record.pop('backward')
        forward = paths[name+'_forward']
        terminal = [r for r in fresh if r['stage'] in ('component_end', 'component_finder_error')]
        assert len(terminal) == 1
        points.append(dict(point, runtime=runtime, component_audit_rows_exact=len(fresh),
            native_terminal=terminal[0], expected_hits=expected,
            actual_hits=[r['hit'] for r in measurements],
            missing_hits=[hit for hit in expected if hit not in {r['hit'] for r in measurements}],
            native_log=record, forward_log_records=len(forward),
            native_forward_first_position_F_mm=forward[0, :3].tolist(),
            native_forward_last_position_F_mm=forward[-1, :3].tolist(),
            native_forward_z_range_F_mm=[float(forward[:, 2].min()), float(forward[:, 2].max())]))
    output.mkdir(exist_ok=False)
    np.savez_compressed(output/'native_log_paths.npz', **paths)
    write_json(output/'summary.json', dict(groups={k:dict(v) for k, v in groups.items()}, points=points,
        all_eight_component_replays_exact=True,
        scope='All 3520 solver-coverage records retained; four selected native navigation examples per version. '
              'Log positions are printed by the native ACTS logger at its limited decimal precision, '
              'not full-precision numerical transport comparisons. No propagation or fit is implemented here.'))
    print(json.dumps(dict(groups={k:dict(v) for k, v in groups.items()},
        replayed_components=len(points), component_rows_exact=sum(p['component_audit_rows_exact'] for p in points))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    analyze(args.input, args.output)
