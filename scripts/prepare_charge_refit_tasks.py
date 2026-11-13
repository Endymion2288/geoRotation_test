#!/usr/bin/env python3
"""Freeze native beforeGeo refit-chart controls after the S06 build."""
import argparse
import json
from pathlib import Path
import shutil


MODES = {'native': [], 'smoother': ['--regular-refit-smoothing-chart'],
         'updater': ['--regular-refit-update-chart'],
         'both': ['--regular-refit-smoothing-chart', '--regular-refit-update-chart']}


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['smoke', 'full'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    code = root/'refit_control_v1_code'
    if args.phase == 'smoke':
        assert (root/'refit_control_build_v1_job/exit_code.txt').read_text().strip() == '0'
        overlay = root/'overlay_refit_before'
        overlay.mkdir(exist_ok=False)
        sources = [root/'software_refit_before', root/'observer']
        provenance = {}
        for folder in ['lib', 'python']:
            (overlay/folder).mkdir()
            for source in sources:
                for item in (source/folder).iterdir():
                    target = overlay/folder/item.name
                    assert not target.exists()
                    target.symlink_to(item.resolve(), target_is_directory=item.is_dir())
                    provenance[str(target.relative_to(overlay))] = str(item.resolve())
        dump(overlay/'software.json', dict(sources=list(map(str, sources)), files=provenance,
            control='Native beforeGeo S06 observers plus default-OFF S05 regular refit charts'))
        for folder in ['configs', 'scripts']:
            (code/folder).mkdir(parents=True, exist_ok=False)
        shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
        for name in ['run_slope_task.py', 'ntuple_tracking_check.py', 'charge_diagnostics.py']:
            shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
        output = root/'refit_control_smoke_v1'
        output.mkdir(exist_ok=False)
    else:
        valid = json.loads((root/'refit_control_smoke_v1/validation.json').read_text())
        assert valid['native_58_fields_exact'] and valid['all_CKF_rows_exact']
        output = root/'refit_control_100k_v1'
        output.mkdir(exist_ok=False)
    combined = []
    for dsid in ['100116', '100117']:
        main_root = root/f'main_{dsid}_v1'
        original_manifest = json.loads((main_root/'manifest.json').read_text())
        if args.phase == 'smoke':
            original = [json.loads((root/'charge_smoke/beforeGeo'/dsid/'on/source_task.json').read_text())]
        else:
            original = [json.loads((Path(job['output'])/'source_task.json').read_text()) for job in
                        json.loads((main_root/'beforeGeo_tasks.json').read_text())]
            original = [t for t in original if t['skip'] < 20000]
            assert len(original) == 20 and sum(t['events'] for t in original) == 100000
        for mode, options in MODES.items():
            if args.phase == 'full' and mode == 'native':
                continue  # Reuse the original, exactly validated native shards.
            group = output/dsid/mode
            group.mkdir(parents=True)
            tasks = []
            for i, original_task in enumerate(original):
                task = dict(original_task, runner=str(code/'configs/slope_tracking.py'),
                    options=['--no-seed-writer', '--actual-covariance', '--isolate-smoothing',
                             '--charge-diagnostics', *options],
                    original_output=original_task['output'],
                    output=str(group/'shards/beforeGeo'/f's{i:04d}'))
                path = group/f'beforeGeo_source_{i:04d}.json'
                dump(path, task)
                tasks.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                                                           '--task', str(path)], output=task['output']))
            dump(group/'beforeGeo_tasks.json', tasks)
            dump(group/'manifest.json', dict(original_manifest, stage='S06-refit-chart-control',
                source_events_per_variant=sum(t['events'] for t in original),
                ranges=[{k: t[k] for k in ['file_index', 'skip', 'events', 'sample_entry', 'input']}
                        for t in original], refit_control=mode,
                runtime_invariants=dict(audit=False, reset_event_maps=False, isolate_smoothing=True),
                baseline=str(main_root)))
            combined.extend(tasks)
    dump(output/'tasks.json', combined)
    print(output/'tasks.json', len(combined), 'jobs')


if __name__ == '__main__':
    main()
