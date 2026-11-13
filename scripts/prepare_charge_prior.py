#!/usr/bin/env python3
"""Expand the validated first-order matched seed prior to all five RDO files."""
import argparse
import json
from pathlib import Path


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    validation = json.loads((root/'prior_validation_v1/summary.json').read_text())
    for dsid in ['100044', '100116', '100117']:
        row = validation['samples'][dsid]
        assert row['default_reconstruction_unchanged'] and row['audit_read_only'] and row['seed_means_identical']
        assert row['maximum']['physical_covariance_relative'] < 1e-11
    output = root/'matched_prior_500k_v1'
    output.mkdir(exist_ok=False)
    combined = []
    for dsid in ['100116', '100117']:
        main_root = root/f'main_{dsid}_v1'
        manifest = json.loads((main_root/'manifest.json').read_text())
        group = output/dsid
        group.mkdir()
        tasks = []
        original_tasks = json.loads((main_root/'afterGeo_tasks.json').read_text())
        smoke = json.loads((root/'prior_smoke_v1'/dsid/'matched/source_task.json').read_text())
        for shard, job in enumerate(original_tasks):
            original = json.loads((Path(job['output'])/'source_task.json').read_text())
            task = dict(original, options=smoke['options'], runner=smoke['runner'],
                output=str(group/'shards/afterGeo'/f's{shard:04d}'))
            path = group/f'afterGeo_source_{shard:04d}.json'
            dump(path, task)
            tasks.append(dict(kind='analysis', command=['python3',
                str(root/'prior_smoke_v1_code/scripts/run_slope_task.py'), '--task', str(path)], output=task['output']))
        assert len(tasks) == 100 and manifest['source_events_per_variant'] == 500000
        dump(group/'afterGeo_tasks.json', tasks)
        dump(group/'manifest.json', dict(manifest, stage='S06-matched-prior-control',
            native_baseline=str(main_root), validation=str(root/'prior_validation_v1/summary.json'),
            prior_scope='Full first-order J C_F J^T at each native seed; finite-width nonlinear spherical prior equality is not assumed'))
        combined.extend(tasks)
    dump(output/'tasks.json', combined)
    print(output/'tasks.json', len(combined), 'jobs')


if __name__ == '__main__':
    main()
