#!/usr/bin/env python3
"""Freeze source ranges and native Calypso entry points for S05 Condor jobs."""
import argparse
import json
from pathlib import Path
import shutil


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--name', required=True)
    parser.add_argument('--events', type=int, required=True)
    parser.add_argument('--skip', type=int, default=0)
    parser.add_argument('--chunk', type=int, default=5000)
    parser.add_argument('--dsid', default='100044')
    parser.add_argument('--samples', type=Path, help='Frozen input manifest; defaults to configs/samples.json')
    parser.add_argument('--ntuple-check', action='store_true')
    parser.add_argument('--actual-covariance', action='store_true')
    parser.add_argument('--no-observer', action='store_true')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.root.resolve()/args.name
    root.mkdir(parents=True, exist_ok=False)
    code = root/'code'
    for folder in ['scripts', 'configs']:
        (code/folder).mkdir(parents=True)
    for name in ['run_slope_task.py', 'ntuple_tracking_check.py']:
        shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
    shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
    samples = json.loads((args.samples or repo/'configs/samples.json').read_text())['samples']
    sample = next(s for s in samples if s['dsid'] == args.dsid)
    if not 0 <= args.skip < args.skip+args.events <= sample['total_entries']:
        raise ValueError('Requested range is outside the frozen sample')
    files = sample['files']
    inputs, offset = [], 0
    for i, f in enumerate(files):
        begin, end = max(args.skip, offset), min(args.skip+args.events, offset+f['entries'])
        if end > begin:
            for start in range(begin, end, args.chunk):
                inputs.append(dict(input=f['path'], file_index=i, skip=start-offset,
                                   events=min(args.chunk, end-start), sample_entry=start))
        offset += f['entries']
    options = ['--no-seed-writer']
    if args.ntuple_check:
        options.append('--ntuple-check')
    if args.actual_covariance:
        options.append('--actual-covariance')
    if args.no_observer:
        if not args.ntuple_check:
            parser.error('--no-observer requires --ntuple-check')
        options.append('--no-observer')
    materials = {
        'beforeGeo': '/eos/home-x/xcheng/FASER/geoRotation_test_results/stage02_closure/material-current-source.json',
        'afterGeo': '/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json',
    }
    for variant, material in materials.items():
        tasks = []
        for shard, source in enumerate(inputs):
            task = dict(source, variant=variant, material=material, dsid=args.dsid,
                        runner=str(code/'configs/slope_tracking.py'), options=options,
                        output=str(root/'shards'/variant/f's{shard:04d}'))
            source_file = root/f'{variant}_source_{shard:04d}.json'
            dump(source_file, task)
            tasks.append(dict(kind='analysis',
                command=['python3', str(code/'scripts/run_slope_task.py'), '--task', str(source_file)],
                output=task['output']))
        dump(root/(variant+'_tasks.json'), tasks)
    recorded_args = {key: str(value) if isinstance(value, Path) else value for key,value in vars(args).items()}
    dump(root/'manifest.json', dict(stage='S05', sample=sample, args=recorded_args,
        source_events_per_variant=sum(s['events'] for s in inputs), ranges=inputs,
        native_tracking=True, coordinate_only_observer=True,
        held_out_split_sample_entry=250000 if args.dsid == '100044' else args.skip+args.events//2))
    print(root)


if __name__ == '__main__':
    main()
