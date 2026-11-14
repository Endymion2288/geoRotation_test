#!/usr/bin/env python3
"""Check native compatibility and both actual refit paths before expansion."""
import argparse
import json
from pathlib import Path

import numpy as np

from validate_charge_alignment import equal_tables


def full(input_path):
    manifest = json.loads((input_path/'manifest.json').read_text())
    checks = []
    for group in manifest['groups']:
        jobs = json.loads((Path(group['group'])/'beforeGeo_tasks.json').read_text())
        ranges = []
        for job in jobs:
            folder = Path(job['output'])
            assert (folder/'exit_code.txt').read_text().strip() == '0'
            source = json.loads((folder/'source_task.json').read_text())
            run = json.loads((folder/'run.json').read_text())
            assert run['success'] and run['isolate_smoothing']
            active = group['control'] == 'both'
            assert run['refit_smoothing_chart_regular'] == active
            assert run['refit_update_chart_regular'] == active
            expected = (['forward', 'backward'] if group['mode'].endswith('merge') else
                ['backward'] if group['mode'].endswith('backward') else ['forward'])
            if not group['reused_native']:
                actual = run['actual_refit_controls_by_direction']
                assert set(actual) == set(expected)
                for direction in expected:
                    for name in ['RegularRefitSmoothingChart', 'RegularRefitUpdateChart']:
                        assert actual[direction][name] == active
            ranges.append((source['file_index'], source['skip'], source['events']))
        expected_ranges = [(r['file_index'], r['skip'], r['events']) for r in manifest['ranges'][group['dsid']]]
        assert sorted(ranges) == sorted(expected_ranges)
        checks.append(dict(dsid=group['dsid'], mode=group['mode'], control=group['control'],
            jobs=len(jobs), source_events=sum(r[2] for r in ranges),
            all_requested_directions_regular=True, reused_native=group['reused_native']))
    result = dict(groups=len(checks), tasks=sum(r['jobs'] for r in checks),
        all_frozen_ranges_exact=True, all_actual_controls_exact=True, checks=checks)
    (input_path/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'checks'}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['smoke', 'full'], default='smoke')
    parser.add_argument('--input', type=Path)
    args = parser.parse_args()
    if args.phase == 'full':
        full(args.input or args.root/'joint_regression_effective_v1')
        return
    output = args.root/'joint_regression_smoke_v1'
    manifest = json.loads((output/'manifest.json').read_text())
    checks = {}
    for pair in manifest['compatibility']:
        current, original = Path(pair['current']), Path(pair['original'])
        assert (current/'exit_code.txt').read_text().strip() == '0'
        assert (original/'exit_code.txt').read_text().strip() == '0'
        checks[pair['dsid']+'|'+pair['mode']] = dict(fields_exact=len(equal_tables(current, original)))
    directions = {}
    for group in manifest['groups']:
        jobs = json.loads((Path(group['group'])/'beforeGeo_tasks.json').read_text())
        assert len(jobs) == 1
        folder = Path(jobs[0]['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        run = json.loads((folder/'run.json').read_text())
        assert run['success'] and run['isolate_smoothing']
        expected = (['forward', 'backward'] if group['mode'].endswith('merge') else
                    ['backward'] if group['mode'].endswith('backward') else ['forward'])
        actual = run['actual_refit_controls_by_direction']
        assert set(actual) == set(expected)
        active = group['control'] == 'both'
        for direction in expected:
            for name in ['RegularRefitSmoothingChart', 'RegularRefitUpdateChart']:
                assert actual[direction][name] == active
            counters = folder/('mode_refit_'+direction+'.npz')
            # The primary backward path uses the historical "forward" file label.
            if len(expected) == 1:
                counters = folder/'mode_refit_forward.npz'
            with np.load(counters) as data:
                assert np.isin(data['status'], [-1, 0, 1]).all()
        key = '|'.join(group[k] for k in ['dsid', 'mode', 'control'])
        directions[key] = actual
    result = dict(controls=len(directions), compatibility=checks,
        all_native_58_fields_exact=all(v['fields_exact'] == 58 for v in checks.values()),
        all_requested_directions_regular=True, actual_refit_controls_by_direction=directions,
        scope='Four species, six native modes, twenty events each; joint refit control applied to every actual CKF direction')
    (output/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'actual_refit_controls_by_direction'}, indent=2))


if __name__ == '__main__':
    main()
