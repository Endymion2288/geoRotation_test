#!/usr/bin/env python3
"""Check encoding-only configuration and exact native 20-event replay."""
import argparse
import copy
import json
from pathlib import Path

import numpy as np


def check_maps(manifest):
    legacy, corrected = [json.loads(Path(manifest['maps'][key]).read_text())
                         for key in ['legacy', 'corrected']]
    restored = copy.deepcopy(corrected)
    frame = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
    count = 0
    for section in ['Surfaces', 'Volumes']:
        for a, b, c in zip(legacy[section]['entries'], corrected[section]['entries'], restored[section]['entries']):
            x, y, z = [entry['value']['material']['binUtility']['transform'] for entry in [a, b, c]]
            np.testing.assert_array_equal(np.array(x['rotation']).reshape(3, 3), frame.T)
            np.testing.assert_array_equal(np.array(y['rotation']).reshape(3, 3), frame)
            z['rotation'] = x['rotation']
            count += 1
    assert legacy == restored and count == manifest['transformed_entries'] == 12
    return dict(transposed_rotation_arrays=count, every_other_JSON_value_identical=True,
                decoding='Independent logical rows, matching actual ACTS decoder')


def read_table(path, cutoff=None):
    with np.load(path) as file:
        data = dict(file)
    if cutoff is None:
        return data
    mask = data['source_entry'] < cutoff
    count = int(mask.sum())
    assert mask[:count].all() and not mask[count:].any()
    end = int(data['hit_offsets'][count]) if 'hit_offsets' in data else None
    return {key: value[:end] if key == 'hit_identifiers' else
            value[:count+1] if key == 'hit_offsets' else value[:count]
            for key, value in data.items()}


def check_run(folder, task, expected):
    assert (folder/'exit_code.txt').read_text().strip() == '0'
    run = json.loads((folder/'run.json').read_text())
    assert run['success'] and run['args']['material'] == task['material']
    assert run['args']['isolate_smoothing']
    assert sorted(str(Path(p).resolve()) for p in run['loaded_acts_libraries']) == \
           sorted(str(Path(p).resolve()) for p in expected)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['maps', 'smoke', 'full'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()/'material_encoding_control_v1'
    manifest = json.loads((root/'manifest.json').read_text())
    if args.phase == 'maps':
        result = check_maps(manifest)
        target = root/'maps_validation.json'
        target.write_text(json.dumps(result, indent=2)+'\n')
        print(target, result)
        return
    records = {}
    for dsid, group in manifest['groups'].items():
        if args.phase == 'smoke':
            jobs = json.loads(Path(group['smoke']).read_text())
            tasks = [json.loads(Path(job['command'][-1]).read_text()) for job in jobs]
            legacy, corrected = [Path(t['output']) for t in tasks]
            for folder, task in zip([legacy, corrected], tasks):
                check_run(folder, task, group['original_loaded_libraries'])
            first = json.loads((Path(group['source'])/'afterGeo_tasks.json').read_text())[0]
            original = Path(first['output'])
            fields = 0
            for table in ['events', 'truth', 'tracks', 'segments']:
                old = read_table(original/(table+'.npz'), 20)
                repeated = read_table(legacy/(table+'.npz'))
                assert old.keys() == repeated.keys()
                for key in old:
                    np.testing.assert_array_equal(old[key], repeated[key], err_msg=dsid+'/'+table+'/'+key)
                    fields += 1
                if table in ['truth', 'segments']:
                    changed = read_table(corrected/(table+'.npz'))
                    for key in old:
                        if not key.startswith('truth_model_'):
                            np.testing.assert_array_equal(old[key], changed[key], err_msg='input '+dsid+'/'+table+'/'+key)
            records[dsid] = dict(legacy_replay_exact=True, fields=fields, source_events=20,
                corrected_raw_truth_and_segments_exact=True)
        else:
            assert json.loads((root/'smoke_validation.json').read_text())['all_exact']
            old_jobs = json.loads((Path(group['source'])/'afterGeo_tasks.json').read_text())
            jobs = json.loads((Path(group['full'])/'afterGeo_tasks.json').read_text())
            assert len(old_jobs) == len(jobs) == group['tasks']
            events = 0
            for old_job, job in zip(old_jobs, jobs):
                old = json.loads(Path(old_job['command'][-1]).read_text())
                task = json.loads(Path(job['command'][-1]).read_text())
                assert {k: v for k, v in old.items() if k not in ['material', 'output']} == \
                       {k: v for k, v in task.items() if k not in ['material', 'output']}
                assert task['material'] == manifest['maps']['corrected']
                check_run(Path(task['output']), task, group['original_loaded_libraries'])
                actual = json.loads((Path(task['output'])/'summary.json').read_text())
                assert actual['source_events'] == task['events']
                events += actual['source_events']
            assert events == group['source_events']
            records[dsid] = dict(tasks=len(jobs), source_events=events, frozen_ranges_exact=True)
    target = root/(args.phase+'_validation.json')
    target.write_text(json.dumps(dict(all_exact=True, groups=records), indent=2)+'\n')
    print(target, records)


if __name__ == '__main__':
    main()
