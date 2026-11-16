#!/usr/bin/env python3
"""Freeze full 110051 jobs and seed-mean controls for the actual latest Calypso."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['baseline', 'controls', 'transport', 'seed'], required=True)
    parser.add_argument('--tag', help='Unique code and task manifest name')
    parser.add_argument('--files', type=int, nargs='+', help='Explicit input file indices')
    parser.add_argument('--controls', nargs='+', help='Restrict the frozen control matrix')
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    tag=args.tag or args.phase
    code = root/('code_'+tag)
    code.mkdir(exist_ok=False)
    for relative in ['configs/reference_tracking.py', 'scripts/reference_truth_observer.py',
                     'scripts/run_latest_momentum.py', 'analysis/latest_momentum.py']:
        target = code/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, target)
    reference = root/'reference'
    if not reference.exists():
        old = repo.parent/'geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/reference/adapted'
        reference.mkdir()
        for name in ['ckf_residual_perLayer.C', 'CommonUtils.h', 'SetStyle.hpp']:
            shutil.copy2(old/name, reference/name)
    inputs = sorted(Path('/eos/experiment/faser/data0/sim/mc24/foresee/110051/rdo/master_dev').glob('*RDO.root'))
    assert len(inputs) == 20
    options = {'baseline': []} if args.phase == 'baseline' else {
        'truth_qop': ['--seed-override-mode', 'qop'],
        'truth_mean': ['--seed-override-mode', 'full'],
        'no_energy_loss': ['--ckf-energy-loss', '0'],
        'no_scattering': ['--ckf-scattering', '0'],
        'step10': ['--ckf-max-step-mm', '10'],
        'qop_prior_100': ['--initial-variance-inflation', '100', '100', '100', '10000', '100'],
        'qop_prior_1': ['--initial-variance-inflation', '100', '100', '100', '10000', '1'],
        'truth_qop_prior1': ['--seed-override-mode', 'qop',
                           '--initial-variance-inflation', '100', '100', '100', '10000', '1'],
    }
    if args.phase == 'transport':
        options={k:v for k,v in options.items() if not k.startswith('truth_')}
    elif args.phase == 'seed':
        options={k:v for k,v in options.items() if k.startswith('truth_')}
    if args.controls:
        options={k:v for k,v in options.items() if k in args.controls}
    jobs, source_tasks = [], []
    # Start controls on file 0. Expand supported mechanisms after paired checks.
    selected_inputs = [(i,p) for i,p in enumerate(inputs) if i in
        (args.files if args.files is not None else range(20) if args.phase=='baseline' else [0])]
    for group, base_options in options.items():
        for file_index, path in selected_inputs:
            run_options = list(base_options)
            override_counts = None
            if group.startswith('truth_'):
                baseline = root/'runs/baseline'/f'file{file_index:02d}'
                if (root/'baseline_replacements.json').exists():
                    replacements=json.loads((root/'baseline_replacements.json').read_text())
                    baseline=Path(replacements.get(baseline.name,str(baseline)))
                seeds = np.load(baseline/'writer.npz')['seeds']
                override = root/f'{group}_file{file_index:02d}.txt'
                valid = np.all(np.isfinite(seeds[:, 9:15]), axis=1)
                # writer seeds: entry, seedIndex, barcode, reco[6], truth[6].
                np.savetxt(override, np.column_stack([seeds[valid, :2], seeds[valid, 9:15]]),
                           fmt=['%d', '%d', *(['%.17g']*6)])
                run_options += ['--seed-override-file', str(override)]
                override_counts=dict(total_seeds=len(seeds),valid_propagated_truth=int(valid.sum()),
                    missing_truth_unchanged=int((~valid).sum()),
                    covariance='Original full seed covariance unchanged unless explicitly scanned')
            destination = root/'runs'/group/f'file{file_index:02d}'
            task = dict(input=str(path), file_index=file_index, events=1000, skip=0,
                output=str(destination), material=str(root/'material-latest-rowmajor.json'),
                audit=True, retain_audit=False, options=run_options,
                macro=str(reference/'ckf_residual_perLayer.C'))
            if override_counts is not None:
                task['seed_override_coverage']=override_counts
            source = root/'tasks'/group/f'file{file_index:02d}.json'
            source.parent.mkdir(parents=True, exist_ok=True)
            source.write_text(json.dumps(task, indent=2)+'\n')
            source_tasks.append(task)
            jobs.append(dict(kind='analysis', command=['python3',str(code/'scripts/run_latest_momentum.py'),
                '--task',str(source)], output=str(destination)))
    (root/(tag+'_tasks.json')).write_text(json.dumps(jobs, indent=2)+'\n')
    (root/(tag+'_manifest.json')).write_text(json.dumps(dict(tasks=source_tasks,
        source_events=sum(t['events'] for t in source_tasks),
        macro_source_commit='33ca2281b961935456dfd4dc9fc06235e24b551e',
        cuts='GitHub chi2Sum/NDF<20, rMax_smt<90 mm, NDF>9; final SDO purity>0.7; original good-seed gate',
        truth='Native writer t_eQOP (SiHit-linked HepMC p); operation diagnostics also explicit production truth',
        isolate_smoothing=True), indent=2)+'\n')


if __name__ == '__main__':
    main()
