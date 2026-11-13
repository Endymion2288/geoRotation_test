#!/usr/bin/env python3
"""Freeze every rotated wrong-sign primary and a file-stratified normal control."""
import argparse
from collections import Counter
import json
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    smoke = json.loads((root/'full_audit_smoke_v1/validation.json').read_text())
    assert all(v['all_reconstruction_arrays_identical'] and v['fields_exact'] == 58 for v in smoke.values())
    repo = Path(__file__).resolve().parents[1]
    output = root/'charge_replay_v1'
    output.mkdir(exist_ok=False)
    code = root/'charge_replay_v1_code'
    for folder in ['scripts', 'configs']:
        (code/folder).mkdir(parents=True)
    shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
    for name in ['run_slope_task.py', 'run_charge_replay.py', 'validate_charge_audit.py',
                 'slope_stage_extract.py', 'audit_reduce.py', 'charge_audit.py']:
        shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
    events = json.loads((root/'lineage_summary_v1/replay_events.json').read_text())
    targets, tasks, counts = [], [], Counter()
    for dsid, expected in [('100116', 27), ('100117', 32)]:
        rows = events[dsid+'|afterGeo']
        wrong = [r for r in rows if r['wrong']]
        assert len(wrong) == expected
        normal = []
        for file_index in range(5):
            choices = sorted([r for r in rows if not r['wrong'] and r['file_index'] == file_index],
                             key=lambda r: r['sample_entry'])
            assert len(choices) >= 6
            # This rule was frozen before inspecting audit residuals.
            normal.extend(choices[:6])
        for row in [*wrong, *normal]:
            name = f"f{row['file_index']:02d}_e{row['source_entry']:08d}"
            task = dict(input=row['input'], dsid=dsid, file_index=row['file_index'],
                skip=row['source_entry'], events=1, sample_entry=row['sample_entry'], variant='afterGeo',
                material='/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json',
                runner=str(code/'configs/slope_tracking.py'),
                options=['--no-seed-writer', '--actual-covariance', '--isolate-smoothing', '--stage-audit', '--charge-audit'],
                output=str(output/dsid/name), replay_targets=[row],
                original_output=str(root/f"main_{dsid}_v1/shards/afterGeo/s{row['sample_entry']//5000:04d}"))
            path = output/(dsid+'_'+name+'.json')
            path.write_text(json.dumps(task, indent=2)+'\n')
            tasks.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_charge_replay.py'),
                                                        '--task', str(path)], output=task['output']))
            targets.append(dict(dsid=dsid, **row))
            counts[dsid+'|'+('wrong' if row['wrong'] else 'normal')] += 1
    assert len({(r['dsid'], r['file_index'], r['source_entry']) for r in targets}) == len(targets)
    (output/'tasks.json').write_text(json.dumps(tasks, indent=2)+'\n')
    (output/'targets.json').write_text(json.dumps(targets, indent=2)+'\n')
    (output/'manifest.json').write_text(json.dumps(dict(counts=dict(counts),
        selection='All 59 selected wrong signs; first six saved modulus-997 normal controls per source file and charge',
        freeze='Targets frozen before inspection of per-update audit r/S/K/C',
        same_candidate_requirement='All 58 reconstruction arrays must reproduce the original source entry exactly',
        diagnostic_smoke=str(root/'full_audit_smoke_v1/validation.json')), indent=2)+'\n')
    print(json.dumps(dict(counts=counts, jobs=len(tasks)), indent=2))


if __name__ == '__main__':
    main()
