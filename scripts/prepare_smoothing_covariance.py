#!/usr/bin/env python3
"""Freeze full 110051 isolated-smoothing geometry and seed-variance comparisons."""
import argparse
import datetime
import json
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True,
        help='Completed native geometry comparison, including frozen macro and inputs')
    parser.add_argument('--events-per-shard', type=int, default=1000)
    args = parser.parse_args()
    assert args.events_per_shard > 0
    repo = Path(__file__).resolve().parents[1]
    root = args.output.resolve()
    reference = args.reference.resolve()
    root.mkdir(parents=True, exist_ok=True)
    frozen = json.loads((reference/'inputs.json').read_text())
    sample, = frozen['samples']
    files = sample['files']
    actual = sorted(Path(files[0]['path']).parent.glob('*.root'))
    assert [str(p) for p in actual] == [f['path'] for f in files]
    assert all(p.stat().st_size == f['size_bytes'] for p, f in zip(actual, files))
    for relative in ['configs/reference_tracking.py', 'scripts/reference_truth_observer.py']:
        target = root/'code'/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, target)
    shutil.copytree(reference/'reference', root/'reference', dirs_exist_ok=True)
    if not (root/'software_tools').exists():
        (root/'software_tools').symlink_to(reference/'software_tools', target_is_directory=True)
    cases = {}
    for case, origin, inflation in [
            ('beforeGeo', 'beforeGeo', [100, 100, 100, 100, 1000]),
            ('afterGeo', 'afterGeo', [100, 100, 100, 10000, 100000]),
            ('afterGeo_originalVariance', 'afterGeo', [100, 100, 100, 100, 1000])]:
        old = json.loads((reference/(origin+'_tasks.json')).read_text())
        tasks = []
        for task in old:
            for skip in range(0, task['events'], args.events_per_shard):
                tasks.append(dict(task, case=case, isolate_smoothing=True, skip=skip,
                    events=min(args.events_per_shard, task['events']-skip),
                    input_entries=files[task['file_index']]['entries'],
                    initial_variance_inflation=inflation,
                    runner=str(root/'code/configs/reference_tracking.py'),
                    macro=str(root/'reference/adapted/ckf_residual_perLayer.C'),
                    output=str(root/'shards'/case/f"f{task['file_index']:02d}_e{skip:04d}")))
        (root/(case+'_tasks.json')).write_text(json.dumps(tasks, indent=2)+'\n')
        smoke = dict(tasks[0], events=20, case=case+'_smoke',
            output=str(root/'smoke'/case))
        (root/(case+'_smoke_tasks.json')).write_text(json.dumps([smoke], indent=2)+'\n')
        shared = dict(smoke, isolate_smoothing=False,
            output=str(root/'smoke'/(case+'_shared')))
        if case in ['beforeGeo', 'afterGeo']:
            (root/(case+'_shared_smoke_tasks.json')).write_text(json.dumps([shared], indent=2)+'\n')
        cases[case] = dict(variant=tasks[0]['variant'], material=tasks[0]['material'],
            initial_variance_inflation=inflation,
            software='software_before' if origin == 'beforeGeo' else 'software_after')
    frozen.update(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        previous_native_comparison=str(reference), cases=cases,
        isolate_smoothing=True, reset_event_maps=False, audit=False,
        directory_and_file_sizes_rechecked=True, events_per_shard=args.events_per_shard,
        scope='Full 20-file 110051 sample; geometry comparison and afterGeo variance control')
    (root/'inputs.json').write_text(json.dumps(frozen, indent=2)+'\n')
    print(json.dumps(dict(cases=cases, source_events=frozen['source_events']), indent=2))


if __name__ == '__main__':
    main()
