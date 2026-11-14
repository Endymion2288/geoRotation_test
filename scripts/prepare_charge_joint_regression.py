#!/usr/bin/env python3
"""Freeze native/regular-refit validation across species and actual CKF modes."""
import argparse
import json
from pathlib import Path
import shutil


MODES = {
    '3forward': [], '3backward': ['--backward'],
    '3merge': ['--merge-directions'], '4forward': ['--four-station'],
    '4backward': ['--four-station', '--backward'],
    '4merge': ['--four-station', '--merge-directions'],
}
CONTROLS = {'native': [], 'both': ['--regular-refit-smoothing-chart',
                                  '--regular-refit-update-chart']}
SPECIES = ['100116', '100117', '110033', '110051']


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def old_muon_group(root, dsid, mode):
    return (root/'mode_analysis_v1/forward_references'/dsid if mode == '3forward'
            else root/'modes_100k_v1'/dsid/mode)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['smoke', 'full'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    code = root/'joint_regression_v1_code'
    output = root/('joint_regression_smoke_v1' if args.phase == 'smoke'
                   else 'joint_regression_full_v1')
    if args.phase == 'smoke':
        for folder in ['configs', 'scripts']:
            (code/folder).mkdir(parents=True, exist_ok=False)
        shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
        for name in ['run_slope_task.py', 'ntuple_tracking_check.py', 'charge_diagnostics.py']:
            shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
    else:
        validation = json.loads((root/'joint_regression_smoke_v1/validation.json').read_text())
        assert validation['all_native_58_fields_exact']
        assert validation['all_requested_directions_regular']
        assert validation['controls'] == 48
    output.mkdir(exist_ok=False)
    sample_definitions = {s['dsid']: s for s in json.loads((repo/'configs/samples.json').read_text())['samples']}
    samples, ranges = {}, {}
    for dsid in SPECIES:
        if dsid in ['100116', '100117']:
            original = json.loads((root/f'main_{dsid}_v1/manifest.json').read_text())
            sample = original['sample']
            current = ([dict(original['ranges'][0], events=20)] if args.phase == 'smoke'
                       else [r for r in original['ranges'] if r['file_index'] in [3, 4] and r['skip'] < 20000])
            assert sum(r['events'] for r in current) == (20 if args.phase == 'smoke' else 40000)
        else:
            # Read only ROOT metadata; all reconstruction goes through Condor.
            import uproot
            sample = dict(sample_definitions[dsid])
            paths = sorted(Path(sample['files'][0]['path']).parent.glob('*RDO.root'))
            assert len(paths) == (10 if dsid == '110033' else 20)
            files = []
            for path in paths:
                with uproot.open(path) as source:
                    entries = source['CollectionTree'].num_entries
                files.append(dict(path=str(path), entries=entries, size_bytes=path.stat().st_size))
            sample.update(files=files, total_entries=sum(f['entries'] for f in files))
            offset = 0
            current = []
            for i, f in enumerate(files):
                if (args.phase == 'smoke' and i == 0) or (args.phase == 'full' and i >= len(files)//2):
                    for skip in range(0, 20 if args.phase == 'smoke' else f['entries'], 1000):
                        current.append(dict(input=f['path'], file_index=i, skip=skip,
                            events=20 if args.phase == 'smoke' else min(1000, f['entries']-skip),
                            sample_entry=offset+skip))
                offset += f['entries']
            assert sum(r['events'] for r in current) == (20 if args.phase == 'smoke' else 10000)
        samples[dsid], ranges[dsid] = sample, current
    combined, groups, compatibility = [], [], []
    for dsid in SPECIES:
        for mode, mode_options in MODES.items():
            for control, extra in CONTROLS.items():
                group = output/dsid/mode/control
                group.mkdir(parents=True)
                tasks = []
                sources = ranges[dsid]
                reused = args.phase == 'full' and dsid in ['100116', '100117'] and control == 'native'
                if reused:
                    original_jobs = json.loads((old_muon_group(root, dsid, mode)/'beforeGeo_tasks.json').read_text())
                    for job in original_jobs:
                        source = json.loads((Path(job['output'])/'source_task.json').read_text())
                        if source['file_index'] in [3, 4] and source['skip'] < 20000:
                            tasks.append(job)
                    assert len(tasks) == len(sources)
                else:
                    for i, source in enumerate(sources):
                        task = dict(source, dsid=dsid, variant='beforeGeo',
                            material=str(root.parent.parent/'stage02_closure/material-current-source.json'),
                            runner=str(code/'configs/slope_tracking.py'),
                            options=['--no-seed-writer', '--actual-covariance', '--isolate-smoothing',
                                     '--mode-refit-counters', *mode_options, *extra],
                            output=str(group/'shards/beforeGeo'/f's{i:04d}'))
                        path = group/f'beforeGeo_source_{i:04d}.json'
                        dump(path, task)
                        job = dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                                                            '--task', str(path)], output=task['output'])
                        tasks.append(job)
                        combined.append(job)
                dump(group/'beforeGeo_tasks.json', tasks)
                dump(group/'manifest.json', dict(stage='S06-joint-refit-regression', sample=samples[dsid],
                    source_events_per_variant=sum(r['events'] for r in sources), ranges=sources,
                    runtime_invariants=dict(audit=False, reset_event_maps=False, isolate_smoothing=True),
                    control=control, mode=mode, reused_native=reused))
                row = dict(dsid=dsid, mode=mode, control=control, group=str(group), reused_native=reused)
                groups.append(row)
                if args.phase == 'smoke' and control == 'native':
                    current = Path(tasks[0]['output'])
                    if dsid in ['100116', '100117']:
                        original = (root/'charge_smoke/beforeGeo'/dsid/'off' if mode == '3forward'
                                    else root/'mode_smoke_v1/beforeGeo'/dsid/mode)
                    else:
                        task = json.loads((group/'beforeGeo_source_0000.json').read_text())
                        original = group/'legacy_native'
                        task.update(runner=str(root/'mode_counter_v1_code/configs/slope_tracking.py'),
                                    output=str(original))
                        path = group/'legacy_native_source.json'
                        dump(path, task)
                        combined.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                            '--task', str(path)], output=str(original)))
                    compatibility.append(dict(dsid=dsid, mode=mode, current=str(current), original=str(original)))
    dump(output/'tasks.json', combined)
    dump(output/'manifest.json', dict(phase=args.phase, groups=groups, jobs=len(combined),
        compatibility=compatibility, samples=samples, ranges=ranges,
        library=str(root/'overlay_refit_before/lib'),
        source_events_per_species={k: sum(r['events'] for r in v) for k, v in ranges.items()},
        selection_frozen_before_new_results=True,
        scope='Native beforeGeo versus the existing S05 joint regular refit chart; all actual forward/backward and three/four-station merge paths; isolated smoothing ON',
        holdout='Muon files 3/4, entries [0,20000): reserved charge-validation files, with prior 3forward validation already reported. A-prime second half of input files: no joint-refit control previously run; 110033 files 5-9 were outside S05 input. Baseline results have been inspected; this is frozen patch validation, not a new blind discovery sample.'))
    print(output/'tasks.json', len(combined), 'new jobs')


if __name__ == '__main__':
    main()
