#!/usr/bin/env python3
"""Multi-core parallel worker executing concurrent Athena sub-shards on local scratch and merging binary tables."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
import uproot


EVENT_COLS = 16
SIGNAL_COLS = 38
SEED_COLS = 34
TRACK_COLS = 78


def _run_one_subshard(args_tuple):
    sub_idx, sub_spec, task, code_dir, scratch_root = args_tuple
    sub_dir = scratch_root / f'sub_{sub_idx:02d}'
    sub_dir.mkdir(parents=True, exist_ok=True)

    env = dict(os.environ)
    env['OMP_NUM_THREADS'] = '1'
    env['OPENBLAS_NUM_THREADS'] = '1'
    env['MKL_NUM_THREADS'] = '1'

    write_states = bool(task.get('macro'))
    cmd = [
        'python3', str(code_dir / 'configs/latest_validation_tracking.py'),
        '--input', sub_spec['input'],
        '--events', str(sub_spec['events']),
        '--skip', str(sub_spec['skip']),
        '--material', task['material'],
        '--audit-prefix', 'fast_audit',
        *task.get('options', []),
    ]
    if write_states:
        cmd.append('--write-states-root')

    t0 = time.monotonic()
    log_path = sub_dir / 'tracking.log'
    with log_path.open('w') as logf:
        rc = subprocess.run(cmd, cwd=sub_dir, env=env, stdout=logf, stderr=subprocess.STDOUT).returncode
    wall_s = time.monotonic() - t0
    if rc != 0:
        tail = log_path.read_text(errors='replace')[-4000:] if log_path.exists() else ''
        raise RuntimeError(f'Sub-shard {sub_idx} failed (rc={rc}):\n{tail}')

    def load_bin(suffix, ncols):
        p = sub_dir / f'fast_audit_{suffix}.bin'
        if not p.exists() or p.stat().st_size == 0:
            return np.zeros((0, ncols), dtype=np.float64)
        raw = np.fromfile(p, dtype=np.float64)
        assert raw.size % ncols == 0, f'{p}: {raw.size} not divisible by {ncols}'
        arr = raw.reshape(-1, ncols).copy()
        global_offset = float(sub_spec['file_index'] * 100000 + sub_spec['skip'])
        arr[:, 0] += global_offset
        return arr

    ev = load_bin('events', EVENT_COLS)
    sig = load_bin('signals', SIGNAL_COLS)
    seeds = load_bin('seeds', SEED_COLS)
    trks = load_bin('tracks', TRACK_COLS)
    assert len(ev) == sub_spec['events'], f'Expected {sub_spec["events"]} events, got {len(ev)}'

    counts_all_path = None
    counts_good_path = None
    if write_states:
        states_root = sub_dir / 'track_states.root'
        assert states_root.exists(), states_root
        link = sub_dir / 'reference_3station_forward.root'
        if not link.exists():
            link.symlink_to('track_states.root')

        with uproot.open(states_root) as f:
            n_tree = f['tree'].num_entries
        assert len(trks) == n_tree, f'tracks.bin ({len(trks)}) != tree ({n_tree})'

        # Compute seed_good flag per track from seeds.bin (exact ancestor seed in same sub-shard)
        # seed_cols: [0:entry, 1:iseed, ..., 11:p_seed_GeV, 12:q_seed, ..., 15:phi_A, 16:theta_A, ..., 25:truth_q, 26:p_prod_GeV, 27:p_st1_GeV, ..., 29:seed_truth_valid, ..., 33:t_theta_A]
        seed_good_map = {}
        for srow in seeds:
            ent_i = int(round(srow[0]))
            iseed_i = int(round(srow[1]))
            p_s = srow[11]
            q_s = srow[12]
            th_s = srow[16]
            tq = srow[25]
            tp = srow[26]
            tvalid = srow[29] > 0.5
            tth = srow[33]
            is_good = (
                tvalid and np.isfinite(p_s) and np.isfinite(tp) and tp > 0.0
                and abs(p_s / tp - 1.0) < 0.1
                and (q_s * tq > 0.0)
                and np.isfinite(th_s) and np.isfinite(tth)
                and abs(th_s - tth) < 1e-4
            )
            seed_good_map[(ent_i, iseed_i)] = 1 if is_good else 0

        purities = trks[:, 15] if len(trks) else np.zeros(0, dtype=np.float64)
        good_flags = np.asarray(
            [seed_good_map.get((int(round(r[0])), int(round(r[2]))), 0) for r in trks],
            dtype=np.int32,
        ) if len(trks) else np.zeros(0, dtype=np.int32)

        label = task.get('macro_label', 'latest_val')
        macro = task['macro']

        # 1) All GitHub + purity > 0.7 (set seed_good = 1 everywhere)
        with uproot.recreate(sub_dir / 'production.root') as friend:
            friend.mktree('production', {'final_purity': 'float64', 'seed_good': 'int32'})
            if n_tree:
                friend['production'].extend({
                    'final_purity': purities,
                    'seed_good': np.ones(n_tree, dtype=np.int32),
                })
        inv_all = f'{macro}("reference_3station_forward.root","counts_all","{label}","t_eQOP",false)'
        subprocess.run(['root', '-l', '-b', '-q', inv_all], check=True, cwd=sub_dir, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        counts_all_path = sub_dir / 'counts_all.root'

        # 2) GitHub + purity > 0.7 + original good-seed gate
        with uproot.recreate(sub_dir / 'production.root') as friend:
            friend.mktree('production', {'final_purity': 'float64', 'seed_good': 'int32'})
            if n_tree:
                friend['production'].extend({
                    'final_purity': purities,
                    'seed_good': good_flags,
                })
        inv_good = f'{macro}("reference_3station_forward.root","counts_good","{label}","t_eQOP",false)'
        subprocess.run(['root', '-l', '-b', '-q', inv_good], check=True, cwd=sub_dir, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        counts_good_path = sub_dir / 'counts_good.root'

        # Remove bulky track_states.root immediately
        states_root.unlink(missing_ok=True)

    schema_json = json.loads((sub_dir / 'fast_audit_schema.json').read_text())
    run_meta = json.loads((sub_dir / 'run.json').read_text())
    return {
        'sub_idx': sub_idx,
        'wall_s': wall_s,
        'events': ev,
        'signals': sig,
        'seeds': seeds,
        'tracks': trks,
        'schema': schema_json,
        'run_meta': run_meta,
        'log_path': str(log_path),
        'counts_all': str(counts_all_path) if counts_all_path else None,
        'counts_good': str(counts_good_path) if counts_good_path else None,
    }


def _merge_root_histograms(src_paths, dst_path):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    names = [f'h_{c}_{s}_{l}' for c in range(3) for s in (1, 2, 3) for l in (0, 1, 2)]
    merged = {}
    for p in src_paths:
        f = ROOT.TFile.Open(str(p))
        for name in names:
            h = f.Get(name)
            assert h, f'Missing {name} in {p}'
            if name not in merged:
                clone = h.Clone(name)
                clone.SetDirectory(0)
                merged[name] = clone
            else:
                merged[name].Add(h)
        f.Close()
    out = ROOT.TFile.Open(str(dst_path), 'RECREATE')
    for h in merged.values():
        h.Write()
    out.Close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    task = json.loads(args.task.read_text())
    destination = Path(task['output'])
    destination.mkdir(parents=True, exist_ok=True)
    code_dir = Path(__file__).resolve().parents[1]

    scratch_base = Path(os.environ.get('TMPDIR', '/tmp'))
    scratch_root = scratch_base / f'faser_val_shard_{os.getpid()}'
    scratch_root.mkdir(parents=True, exist_ok=True)

    sub_ranges = task['sub_ranges']
    workers = min(len(sub_ranges), int(task.get('workers', len(sub_ranges))))

    t0 = time.monotonic()
    try:
        work_items = [(i, sr, task, code_dir, scratch_root) for i, sr in enumerate(sub_ranges)]
        if workers > 1:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                results = list(pool.map(_run_one_subshard, work_items))
        else:
            results = [_run_one_subshard(w) for w in work_items]

        results.sort(key=lambda r: r['sub_idx'])
        ev = np.vstack([r['events'] for r in results])
        sig = np.vstack([r['signals'] for r in results]) if any(len(r['signals']) for r in results) else np.zeros((0, SIGNAL_COLS))
        seeds = np.vstack([r['seeds'] for r in results]) if any(len(r['seeds']) for r in results) else np.zeros((0, SEED_COLS))
        trks = np.vstack([r['tracks'] for r in results]) if any(len(r['tracks']) for r in results) else np.zeros((0, TRACK_COLS))
        schema = results[0]['schema']

        np.savez_compressed(
            destination / 'shard_data.npz',
            events=ev,
            signals=sig,
            seeds=seeds,
            tracks=trks,
            event_cols=np.asarray(schema['event_cols']),
            signal_cols=np.asarray(schema['signal_cols']),
            seed_cols=np.asarray(schema['seed_cols']),
            track_cols=np.asarray(schema['track_cols']),
        )

        if results[0]['counts_all']:
            _merge_root_histograms([r['counts_all'] for r in results], destination / 'native_counts_all.root')
            _merge_root_histograms([r['counts_good'] for r in results], destination / 'native_counts_good.root')

        with gzip.open(destination / 'tracking.log.gz', 'wb') as gz:
            for r in results:
                lp = Path(r['log_path'])
                if lp.exists():
                    gz.write(f'=== sub_shard {r["sub_idx"]} ===\n'.encode())
                    with lp.open('rb') as src:
                        shutil.copyfileobj(src, gz)

        total_wall = time.monotonic() - t0
        summary = {
            'sample': task['sample'],
            'group': task['group'],
            'shard_id': task['shard_id'],
            'parallel_workers': workers,
            'sub_shards': len(sub_ranges),
            'total_events': int(len(ev)),
            'total_signals': int(len(sig)),
            'total_seeds': int(len(seeds)),
            'total_tracks': int(len(trks)),
            'wall_s': total_wall,
            'sub_walls_s': [r['wall_s'] for r in results],
            'run_meta_0': results[0]['run_meta'],
        }
        (destination / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    finally:
        shutil.rmtree(scratch_root, ignore_errors=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
