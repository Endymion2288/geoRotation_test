#!/usr/bin/env python3
"""Freeze multi-sample step-by-step ablation, particle-hypothesis check, and seeding-improvement tasks."""
import argparse
import json
from pathlib import Path
import shutil


GROUPS = {
    'C0_baseline': [
        '--legacy-swapped-cov-loc',
        '--initial-variance-inflation', '100', '100', '100', '10000', '100000',
        '--seed-b-tesla', '0.55',
        '--fix-gety-nan', '0',
        '--use-local-tangent', '0',
        '--abs-pdg-hypothesis', '13',
    ],
    'C1_cov_loc': [
        '--initial-variance-inflation', '100', '100', '100', '10000', '100000',
        '--seed-b-tesla', '0.55',
        '--fix-gety-nan', '0',
        '--use-local-tangent', '0',
        '--abs-pdg-hypothesis', '13',
    ],
    'C2_orig_var': [
        '--initial-variance-inflation', '100', '100', '100', '100', '1000',
        '--seed-b-tesla', '0.55',
        '--fix-gety-nan', '0',
        '--use-local-tangent', '0',
        '--abs-pdg-hypothesis', '13',
    ],
    'C3_seed_b': [
        '--initial-variance-inflation', '100', '100', '100', '100', '1000',
        '--seed-b-tesla', '0.473662',
        '--fix-gety-nan', '0',
        '--use-local-tangent', '0',
        '--abs-pdg-hypothesis', '13',
    ],
    'C4_no_nan': [
        '--initial-variance-inflation', '100', '100', '100', '100', '1000',
        '--seed-b-tesla', '0.473662',
        '--fix-gety-nan', '1',
        '--use-local-tangent', '0',
        '--abs-pdg-hypothesis', '13',
    ],
    'C5_local_tan': [
        '--initial-variance-inflation', '100', '100', '100', '100', '1000',
        '--seed-b-tesla', '0.473662',
        '--fix-gety-nan', '1',
        '--use-local-tangent', '1',
        '--abs-pdg-hypothesis', '13',
    ],
    'H1_electron_hypo': [
        '--initial-variance-inflation', '100', '100', '100', '100', '1000',
        '--seed-b-tesla', '0.473662',
        '--fix-gety-nan', '1',
        '--use-local-tangent', '1',
        '--abs-pdg-hypothesis', '11',
    ],
    'S1_dedup_prefilter': [
        '--initial-variance-inflation', '100', '100', '100', '100', '1000',
        '--seed-b-tesla', '0.473662',
        '--fix-gety-nan', '1',
        '--use-local-tangent', '1',
        '--abs-pdg-hypothesis', '13',
        '--seed-prefilter', '1',
        '--seed-improved-dedup', '1',
        '--seed-max-reduced-chi2', '35.0',
    ],
    'S2_strip_field_fit': [
        '--initial-variance-inflation', '100', '100', '100', '100', '1000',
        '--seed-b-tesla', '0.473662',
        '--fix-gety-nan', '1',
        '--use-local-tangent', '1',
        '--abs-pdg-hypothesis', '13',
        '--seed-prefilter', '1',
        '--seed-improved-dedup', '1',
        '--seed-strip-field-fit', '1',
        '--seed-max-reduced-chi2', '35.0',
    ],
    'S3_strip_only': [
        '--initial-variance-inflation', '100', '100', '100', '100', '1000',
        '--seed-b-tesla', '0.473662',
        '--fix-gety-nan', '1',
        '--use-local-tangent', '1',
        '--abs-pdg-hypothesis', '13',
        '--seed-prefilter', '0',
        '--seed-improved-dedup', '0',
        '--seed-strip-field-fit', '1',
        '--seed-max-reduced-chi2', '1000000000.0',
    ],
}


def sample_shards():
    shards = []
    # 1. 110051: A' -> e+e- (eps=1e-4), all 20 files = 20,000 events (20 shards, each 2x500 parallel)
    files_110051 = sorted(Path('/eos/experiment/faser/data0/sim/mc24/foresee/110051/rdo/master_dev').glob('*-RDO.root'))
    assert len(files_110051) == 20
    for fi, fp in enumerate(files_110051):
        shards.append({
            'sample': '110051',
            'shard_id': f'file{fi:02d}',
            'run_macro': True,
            'sub_ranges': [
                {'input': str(fp), 'file_index': fi, 'skip': 0, 'events': 500},
                {'input': str(fp), 'file_index': fi, 'skip': 500, 'events': 500},
            ],
        })

    # 2. 110033: A' -> e+e- (eps=1e-5), first 5 files = 10,000 events (5 shards, each 2x1000 parallel)
    files_110033 = sorted(Path('/eos/experiment/faser/data0/sim/mc24/foresee/110033/rdo/master_dev').glob('*-RDO.root'))[:5]
    assert len(files_110033) == 5
    for fi, fp in enumerate(files_110033):
        shards.append({
            'sample': '110033',
            'shard_id': f'file{fi:02d}',
            'run_macro': False,
            'sub_ranges': [
                {'input': str(fp), 'file_index': fi, 'skip': 0, 'events': 1000},
                {'input': str(fp), 'file_index': fi, 'skip': 1000, 'events': 1000},
            ],
        })

    # 3. 100116: PG mu+ 2D FLUKA, [0, 10000) = 10,000 events (2 shards, each 2x2500 parallel)
    fp_100116 = '/eos/experiment/faser/data0/sim/mc24/particle_gun/100116/rdo/s0012/FaserMC-MC24_PG_mupl_fasernu_2D_fluka-100116-00000-00009-s0012-RDO.root'
    for si in range(2):
        base = si * 5000
        shards.append({
            'sample': '100116',
            'shard_id': f'shard{si:02d}',
            'run_macro': False,
            'sub_ranges': [
                {'input': fp_100116, 'file_index': 0, 'skip': base, 'events': 2500},
                {'input': fp_100116, 'file_index': 0, 'skip': base + 2500, 'events': 2500},
            ],
        })

    # 4. 100117: PG mu- 2D FLUKA, [0, 10000) = 10,000 events (2 shards, each 2x2500 parallel)
    fp_100117 = '/eos/experiment/faser/data0/sim/mc24/particle_gun/100117/rdo/s0012/FaserMC-MC24_PG_mumi_fasernu_2D_fluka-100117-00000-00009-s0012-RDO.root'
    for si in range(2):
        base = si * 5000
        shards.append({
            'sample': '100117',
            'shard_id': f'shard{si:02d}',
            'run_macro': False,
            'sub_ranges': [
                {'input': fp_100117, 'file_index': 0, 'skip': base, 'events': 2500},
                {'input': fp_100117, 'file_index': 0, 'skip': base + 2500, 'events': 2500},
            ],
        })

    # 5. 100130: PG K_S -> pi+ pi- at end_fasernu, [0, 10000) = 10,000 events (4 shards, each 2x1250 parallel)
    fp_100130 = '/eos/experiment/faser/data0/sim/mc24/particle_gun/100130/rdo/s0013/FaserMC-MC24_PG_kshort_end_fasernu_xin-100130-00000-00009-s0013-RDO.root'
    for si in range(4):
        base = si * 2500
        shards.append({
            'sample': '100130',
            'shard_id': f'shard{si:02d}',
            'run_macro': False,
            'sub_ranges': [
                {'input': fp_100130, 'file_index': 0, 'skip': base, 'events': 1250},
                {'input': fp_100130, 'file_index': 0, 'skip': base + 1250, 'events': 1250},
            ],
        })
    return shards


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--tag', default='validation_campaign_v1')
    parser.add_argument('--groups', nargs='+', default=list(GROUPS.keys()))
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    root = args.root.resolve()
    material = root / 'material-latest-rowmajor.json'
    assert material.is_file(), material

    reference = root / 'reference'
    if not reference.exists():
        old_ref = repo.parent / 'geoRotation_test_results/stage04_seed_purity/20261113_seed_good_v2/reference/adapted'
        reference.mkdir(parents=True)
        for name in ['ckf_residual_perLayer.C', 'CommonUtils.h', 'SetStyle.hpp']:
            shutil.copy2(old_ref / name, reference / name)

    code_dir = root / f'code_{args.tag}'
    if code_dir.exists():
        shutil.rmtree(code_dir)
    for rel in [
        'configs/latest_validation_tracking.py',
        'scripts/run_validation_shard.py',
    ]:
        dst = code_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / rel, dst)

    task_sources = root / f'tasks_{args.tag}'
    task_sources.mkdir(parents=True, exist_ok=True)
    shards = sample_shards()

    condor_tasks = []
    for gname in args.groups:
        opts = GROUPS[gname]
        for sh in shards:
            sample = sh['sample']
            sid = sh['shard_id']
            out_dir = root / 'runs' / gname / sample / sid
            src_payload = {
                'group': gname,
                'sample': sample,
                'shard_id': sid,
                'workers': 2,
                'sub_ranges': sh['sub_ranges'],
                'material': str(material),
                'options': opts,
                'macro': str(reference / 'ckf_residual_perLayer.C') if sh['run_macro'] else None,
                'macro_label': gname,
                'output': str(out_dir),
            }
            src_file = task_sources / f'{gname}_{sample}_{sid}.json'
            src_file.write_text(json.dumps(src_payload, indent=2) + '\n')
            condor_tasks.append({
                'kind': 'analysis',
                'command': ['python3', str(code_dir / 'scripts/run_validation_shard.py'), '--task', str(src_file)],
                'output': str(out_dir),
            })

    tasks_json = root / f'{args.tag}_tasks.json'
    tasks_json.write_text(json.dumps(condor_tasks, indent=2) + '\n')
    manifest = {
        'groups': list(args.groups),
        'samples': {
            '110051': '20 files, 20,000 events (Aprime->e+e-, eps=1e-4; held-out 19 files exclude file00)',
            '110033': '5 files, 10,000 events (Aprime->e+e-, eps=1e-5)',
            '100116': '10,000 events (PG mu+ 2D FLUKA)',
            '100117': '10,000 events (PG mu- 2D FLUKA)',
            '100130': '10,000 events (PG K_S hadrons at end_fasernu)',
        },
        'shards_per_group': len(shards),
        'total_condor_jobs': len(condor_tasks),
        'parallel_workers_per_job': 2,
        'total_reconstructed_events': len(args.groups) * 60000,
    }
    (root / f'{args.tag}_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
