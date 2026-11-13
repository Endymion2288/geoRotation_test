#!/usr/bin/env python3
"""Prepare q/p-only seed controls through the existing actual Calypso CKF."""
import argparse
import gzip
import json
from pathlib import Path
import shutil

import numpy as np


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def override(folder, mode, output):
    with gzip.open(folder/'charge_audit.jsonl.gz', 'rt') as source:
        rows = [json.loads(line) for line in source]
    seeds = [row for row in rows if row['stage'] == 'seed']
    original = json.loads((folder/'source_task.json').read_text())
    pdg = 13 if original['dsid'] == '100117' else -13
    primary = {row['source_entry']: row for row in rows if row['stage'] == 'truth'
               and row['pdg'] == pdg and row['status'] == 1 and 0 < row['barcode'] < 200000}
    lines = []
    for seed in seeds:
        mean = seed['bound'].copy()
        if mode == 'production_qop':
            particle = primary[seed['source_entry']]
            p = np.linalg.norm(particle['momentum_F_MeV'])/1000.
            mean[4] = (-1 if pdg > 0 else 1)/p
        lines.append(' '.join([str(seed['source_entry']-original['skip']), str(seed['seed']),
                               *[format(value, '.17g') for value in mean]]))
    output.write_text('\n'.join(lines)+'\n')
    return len(seeds)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['smoke', 'replay', 'mean-smoke', 'centroid-fixed'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    code = root/'seed_qop_v1_code'
    if args.phase == 'smoke':
        for folder in ['configs', 'scripts']:
            (code/folder).mkdir(parents=True, exist_ok=False)
        shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
        for name in ['run_slope_task.py', 'ntuple_tracking_check.py', 'slope_stage_extract.py',
                     'audit_reduce.py', 'charge_audit.py']:
            shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
        originals = [root/'hitcontrol_smoke_v2'/dsid/'default/f00_e00000000'
                     for dsid in ['100044', '100116', '100117']]
        modes = ['identical', 'production_qop']
        output = root/'seed_qop_smoke_v1'
    elif args.phase == 'replay':
        valid = json.loads((root/'seed_qop_smoke_v1/validation.json').read_text())
        assert valid['all_identical_58_fields_exact'] and valid['all_changes_qop_only_with_identical_covariance']
        originals = [Path(job['output']) for job in json.loads((root/'charge_replay_v1/tasks.json').read_text())]
        modes = ['production_qop']
        output = root/'seed_qop_replay_v1'
    elif args.phase == 'mean-smoke':
        originals = [root/'hitcontrol_smoke_v2'/dsid/'default/f00_e00000000'
                     for dsid in ['100044', '100116', '100117']]
        modes = ['identical_mean']
        output = root/'seed_mean_smoke_v1'
    else:
        centroid = json.loads((root/'centroid_analysis_v1/events.json').read_text())
        retained = [row for row in centroid if row['mode'] == 'centroid'
                    and (not row['original_wrong'] or row['original_final_purity'] == 1.)]
        assert len(retained) == 66
        for row in retained:
            response = row['seed_response']
            assert response['original_seeds'] == response['centroid_seeds'] == response['same_id_hits_surface']
            assert response['removed_seed_hit_sets'] == response['added_seed_hit_sets'] == response['covariance_changed'] == 0
        lookup = {row['dsid']+'_'+str(row['sample_entry']) for row in retained}
        originals = []
        for job in json.loads((root/'charge_replay_v1/tasks.json').read_text()):
            folder = Path(job['output'])
            task = json.loads((folder/'source_task.json').read_text())
            if task['dsid']+'_'+str(task['sample_entry']) in lookup:
                originals.append(folder)
        assert len(originals) == 66
        modes = ['centroid_fixed_seed']
        output = root/'centroid_fixed_seed_v1'
    output.mkdir(exist_ok=False)
    jobs, sources = [], []
    for original in originals:
        task = json.loads((original/'source_task.json').read_text())
        for mode in modes:
            name = task['dsid']+'_'+mode+f"_f{task['file_index']:02d}_e{task['skip']:08d}"
            filename = output/(name+'_seeds.txt')
            seeds = override(original, mode, filename)
            destination = output/task['dsid']/mode/f"f{task['file_index']:02d}_e{task['skip']:08d}"
            fixed = mode == 'centroid_fixed_seed'
            options = ['--no-seed-writer', '--actual-covariance', '--isolate-smoothing',
                '--stage-audit', '--charge-audit', '--seed-override-file', str(filename),
                '--seed-override-mode', 'mean' if fixed or mode == 'identical_mean' else 'qop']
            if fixed:
                options.append('--truth-cluster-position')
            source = dict(task, runner=str(code/'configs/slope_tracking.py'), options=options,
                output=str(destination), original_output=str(original), override_seed_rows=seeds)
            sourcefile = output/(name+'.json')
            dump(sourcefile, source)
            jobs.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                                                       '--task', str(sourcefile)], output=str(destination)))
            sources.append(source)
    dump(output/'tasks.json', jobs)
    dump(output/'manifest.json', dict(phase=args.phase, jobs=len(jobs), source_tasks=sources,
        software=str(root/('overlay_centroid_v1' if args.phase == 'centroid-fixed' else 'overlay_hitcontrol_v2')),
        control=('Actual CKF SeedOverrideMode=mean; identical six native seed mean parameters'
                 if args.phase == 'mean-smoke' else
                 'Actual CKF SeedOverrideMode=mean. SiHit-centroid measurements with the original native six seed means and unchanged full covariance, surfaces and seed hits'
                 if args.phase == 'centroid-fixed' else
                 'Actual CKF SeedOverrideMode=qop. Replace each native seed q/p by charge / stored production p; all other mean elements and full covariance retained'),
        interpretation=('Six native purity-one wrong signs plus the original 60 normal controls; choose only events whose centroid seed IDs/hits/surfaces/full covariance already exactly match native. Native CKF can still branch and change selected hits.'
                        if args.phase == 'centroid-fixed' else
                        'Production momentum is not local sensor G4 momentum. This diagnoses mean-curvature sensitivity and does not define a perfect truth seed or fix covariance.')))
    print(output/'tasks.json', len(jobs), 'jobs')


if __name__ == '__main__':
    main()
