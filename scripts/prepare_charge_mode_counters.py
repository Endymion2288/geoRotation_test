#!/usr/bin/env python3
"""Record explicit native per-direction refit outcomes with unchanged tracks."""
import argparse
import json
from pathlib import Path
import shutil


MODES = ['3forward', '3backward', '3merge', '4forward', '4backward', '4merge']


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['smoke', 'full'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    code = root/'mode_counter_v1_code'
    if args.phase == 'smoke':
        for folder in ['scripts', 'configs']:
            (code/folder).mkdir(parents=True, exist_ok=False)
        for name in ['run_charge_mode_counter.py', 'run_slope_task.py', 'ntuple_tracking_check.py',
                     'charge_diagnostics.py', 'validate_charge_audit.py']:
            shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
        shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
        output = root/'mode_counter_smoke_v1'
    else:
        smoke = root/'mode_counter_smoke_v1'
        for variant in ['beforeGeo', 'afterGeo']:
            for job in json.loads((smoke/(variant+'_tasks.json')).read_text()):
                folder = Path(job['output'])
                assert (folder/'exit_code.txt').read_text().strip() == '0'
                assert json.loads((folder/'counter_validation.json').read_text())['all_58_fields_exact']
        output = root/'mode_counter_100k_v1'
    output.mkdir(exist_ok=False)
    for variant in ['beforeGeo', 'afterGeo']:
        combined = []
        for dsid in ['100116', '100117']:
            for mode in MODES:
                if args.phase == 'smoke':
                    original = (root/'charge_smoke'/variant/dsid/'off' if mode == '3forward'
                                else root/'mode_smoke_v1'/variant/dsid/mode)
                    originals = [json.loads((original/'source_task.json').read_text())]
                else:
                    if mode == '3forward':
                        continue  # Already covered by exact original main charge lineage.
                    jobs = json.loads((root/'modes_100k_v1'/dsid/mode/(variant+'_tasks.json')).read_text())
                    originals = [json.loads((Path(job['output'])/'source_task.json').read_text()) for job in jobs]
                    assert len(originals) == 20 and sum(t['events'] for t in originals) == 100000
                for i, original in enumerate(originals):
                    options = [x for x in original['options'] if x != '--charge-diagnostics']
                    task = dict(original, options=[*options, '--mode-refit-counters'],
                        runner=str(code/'configs/slope_tracking.py'), original_output=original['output'],
                        output=str(output/variant/dsid/mode/f's{i:04d}'))
                    path = output/f'{variant}_{dsid}_{mode}_{i:04d}.json'
                    dump(path, task)
                    combined.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_charge_mode_counter.py'),
                                                                 '--task', str(path)], output=task['output']))
        dump(output/(variant+'_tasks.json'), combined)
        print(variant, len(combined), 'jobs')
    dump(output/'manifest.json', dict(phase=args.phase,
        scope='Read-only ChargeDiagnosticFile per native CKF direction; exact 58-field comparison to original output before compact duplicate tables are removed',
        forward_reference='Existing main lineage provides all 3forward counts'))


if __name__ == '__main__':
    main()
