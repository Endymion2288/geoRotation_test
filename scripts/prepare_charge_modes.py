#!/usr/bin/env python3
"""Freeze balanced file coverage for native backward, IFT and merge controls."""
import argparse
import json
from pathlib import Path


MODES = {'3backward': ['--backward'], '4forward': ['--four-station'],
         '4backward': ['--four-station', '--backward'],
         '3merge': ['--merge-directions'], '4merge': ['--four-station', '--merge-directions']}


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    validation = json.loads((root/'mode_smoke_v1/validation.json').read_text())
    assert len(validation['default_agreement']) == 4 and len(validation['actual_merge']) == 8
    output = root/'modes_100k_v1'
    output.mkdir(exist_ok=False)
    # Reuse the exact entry points already validated in all 24 mode jobs.
    code = root/'mode_smoke_v1_code'
    combined = {v: [] for v in ['beforeGeo', 'afterGeo']}
    for dsid in ['100116', '100117']:
        source_root = root/f'main_{dsid}_v1'
        original = json.loads((source_root/'manifest.json').read_text())
        ranges = [r for r in original['ranges'] if r['skip'] < 20000]
        assert len(ranges) == 20 and sum(r['events'] for r in ranges) == 100000
        assert {r['file_index'] for r in ranges} == set(range(5))
        for mode, extra in MODES.items():
            group = output/dsid/mode
            group.mkdir(parents=True)
            for variant in combined:
                first = json.loads((source_root/f'{variant}_source_0000.json').read_text())
                tasks = []
                for shard, source in enumerate(ranges):
                    task = dict(source, dsid=dsid, variant=variant, material=first['material'],
                        runner=str(code/'configs/slope_tracking.py'),
                        options=['--no-seed-writer', '--actual-covariance', '--isolate-smoothing', *extra],
                        output=str(group/'shards'/variant/f's{shard:04d}'))
                    path = group/f'{variant}_source_{shard:04d}.json'
                    dump(path, task)
                    job = dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                                                        '--task', str(path)], output=task['output'])
                    tasks.append(job)
                dump(group/(variant+'_tasks.json'), tasks)
                combined[variant].extend(tasks)
            manifest = dict(original, stage='S06-mode-control', ranges=ranges,
                source_events_per_variant=100000,
                runtime_invariants=dict(audit=False, reset_event_maps=False, isolate_smoothing=True),
                native_mode=mode, source_rule='First 20,000 entries of each of five source files',
                main_forward_reference=str(source_root), default_and_merge_smoke=str(root/'mode_smoke_v1/validation.json'))
            dump(group/'manifest.json', manifest)
    for variant, tasks in combined.items():
        dump(output/(variant+'_tasks.json'), tasks)
        print(variant, len(tasks), 'jobs', sum(r['events'] for r in ranges), 'events per charge/mode')


if __name__ == '__main__':
    main()
