#!/usr/bin/env python3
"""Freeze an encoding-only material control on the existing S05/S06 inputs."""
import argparse
import copy
import json
from pathlib import Path

import numpy as np


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output = root/'material_encoding_control_v1'
    output.mkdir(exist_ok=False)
    material_path = Path('/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json')
    original = json.loads(material_path.read_text())
    corrected = copy.deepcopy(original)
    frame = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
    changed = 0
    for section in ['Surfaces', 'Volumes']:
        for entry in corrected[section]['entries']:
            transform = entry['value']['material']['binUtility']['transform']
            rotation = np.array(transform['rotation']).reshape(3, 3)
            np.testing.assert_array_equal(rotation, frame.T)
            transform['rotation'] = rotation.T.flatten().tolist()
            changed += 1
    maps = {}
    for name, data in [('legacy', original), ('corrected', corrected)]:
        path = output/(name+'_afs_map.json')
        path.write_text(json.dumps(data, separators=(',', ':'))+'\n')
        maps[name] = str(path)
    sources = {dsid: root/('main_'+dsid+'_v1') for dsid in ['100116', '100117']}
    sources['100044'] = root.parents[1]/'stage05_slopes/isolated_100044_v1'
    groups = {}
    for dsid, source in sources.items():
        old_jobs = json.loads((source/'afterGeo_tasks.json').read_text())
        first = json.loads(Path(old_jobs[0]['command'][-1]).read_text())
        original_run = json.loads((Path(first['output'])/'run.json').read_text())
        library = original_run['library']
        assert '--isolate-smoothing' in first['options']
        group = output/dsid
        group.mkdir()
        smoke = []
        for control in maps:
            task = dict(first, events=20, material=maps[control],
                        output=str(group/'smoke'/control))
            path = group/('smoke_'+control+'_source.json')
            dump(path, task)
            command = old_jobs[0]['command'][:-1]+[str(path)]
            smoke.append(dict(kind='analysis', command=command, output=task['output']))
        dump(group/'smoke_tasks.json', smoke)
        full = group/'full'
        full.mkdir()
        jobs = []
        event_count = 0
        for i, old_job in enumerate(old_jobs):
            old_task = json.loads(Path(old_job['command'][-1]).read_text())
            assert old_task['material'] == str(material_path)
            task = dict(old_task, material=maps['corrected'],
                        output=str(full/'shards/afterGeo'/f's{i:04d}'))
            path = full/f'afterGeo_source_{i:04d}.json'
            dump(path, task)
            jobs.append(dict(kind='analysis', command=old_job['command'][:-1]+[str(path)],
                             output=task['output']))
            event_count += task['events']
        dump(full/'afterGeo_tasks.json', jobs)
        manifest = json.loads((source/'manifest.json').read_text())
        manifest['material_encoding_control'] = dict(source=str(source), maps=maps,
            actual_library=library, unchanged='Existing frozen runner, worker analysis, inputs, '
            'ranges, flags and Calypso libraries; only material rotation serialization changes')
        dump(full/'manifest.json', manifest)
        groups[dsid] = dict(source=str(source), library=library, smoke=str(group/'smoke_tasks.json'),
            full=str(full), tasks=len(jobs), source_events=event_count,
            original_loaded_libraries=original_run['loaded_acts_libraries'])
    dump(output/'manifest.json', dict(groups=groups, maps=maps, source_material=str(material_path),
        source_bytes=material_path.stat().st_size, transformed_entries=changed,
        change='Transpose only the 12 existing serialized rotation arrays. Geometry IDs, '
               'translations, bin axes, bin data, metadata and surface presence are identical.',
        gate='Validate the legacy 20-event replay against the archived prefix before full submission.'))
    print(output, groups)


if __name__ == '__main__':
    main()
