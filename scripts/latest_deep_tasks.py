#!/usr/bin/env python3
"""Freeze full 20-file task matrices for latest-Calypso seed, covariance, and B-field controls."""
import argparse
import json
from pathlib import Path
import shutil


# Effective B from file 0 [0,100) deterministic native full-field propagation + CircleFit closure:
# file 0 must be excluded when evaluating selection windows or calibrated seed gates!
CALIBRATED_SEED_B_TESLA = '0.473662'

GROUP_OPTIONS = {
    'baseline': [],
    'fix_cov_loc_swap': ['--swap-cov-loc'],
    'qop_1000_only': ['--initial-variance-inflation', '100', '100', '100', '10000', '1000'],
    'qop_1000_loc_fixed': ['--swap-cov-loc', '--initial-variance-inflation', '100', '100', '100', '10000', '1000'],
    'seed_b_04737': ['--seed-b-tesla', CALIBRATED_SEED_B_TESLA],
    'combined_phys_fix': [
        '--seed-b-tesla', CALIBRATED_SEED_B_TESLA,
        '--swap-cov-loc',
        '--initial-variance-inflation', '100', '100', '100', '10000', '1000',
    ],
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--tag', required=True)
    parser.add_argument('--groups', nargs='+', default=list(GROUP_OPTIONS.keys()))
    parser.add_argument('--files', type=int, nargs='+')
    args = parser.parse_args()

    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    code = root / f'code_{args.tag}'
    code.mkdir(exist_ok=False)
    for relative in (
        'configs/latest_deep_tracking.py',
        'scripts/reference_truth_observer.py',
        'scripts/run_latest_deep.py',
        'analysis/latest_deep_reduce.py',
    ):
        target = code / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / relative, target)

    reference = root / 'reference'
    assert reference.exists(), reference
    inputs = sorted(Path('/eos/experiment/faser/data0/sim/mc24/foresee/110051/rdo/master_dev').glob('*RDO.root'))
    assert len(inputs) == 20
    file_set = set(args.files) if args.files is not None else set(range(20))

    jobs, source_tasks = [], []
    for group in args.groups:
        assert group in GROUP_OPTIONS, group
        base_options = GROUP_OPTIONS[group]
        for file_index, path in enumerate(inputs):
            if file_index not in file_set:
                continue
            destination = root / 'runs' / group / f'file{file_index:02d}'
            task = dict(
                group=group,
                input=str(path),
                file_index=file_index,
                events=1000,
                skip=0,
                output=str(destination),
                material=str(root / 'material-latest-rowmajor.json'),
                audit=True,
                retain_audit=False,
                options=list(base_options),
                macro=str(reference / 'ckf_residual_perLayer.C'),
                macro_label=f'latest_{group}',
            )
            source = root / 'tasks' / group / f'file{file_index:02d}.json'
            source.parent.mkdir(parents=True, exist_ok=True)
            source.write_text(json.dumps(task, indent=2) + '\n')
            source_tasks.append(task)
            jobs.append(dict(
                kind='analysis',
                command=['python3', str(code / 'scripts/run_latest_deep.py'), '--task', str(source)],
                output=str(destination),
            ))

    (root / f'{args.tag}_tasks.json').write_text(json.dumps(jobs, indent=2) + '\n')
    (root / f'{args.tag}_manifest.json').write_text(json.dumps(dict(
        groups=args.groups,
        tasks_count=len(jobs),
        source_events_per_group=sum(1000 for i in range(20) if i in file_set),
        calibrated_seed_b_tesla=float(CALIBRATED_SEED_B_TESLA),
        calibration_note='0.473662 T determined on file 0 [0,100) only; independent evaluations must exclude file 0',
        macro_source_commit='33ca2281b961935456dfd4dc9fc06235e24b551e',
        isolate_smoothing=True,
    ), indent=2) + '\n')
    print(f'Prepared {len(jobs)} jobs in {root / (args.tag + "_tasks.json")}')


if __name__ == '__main__':
    main()
