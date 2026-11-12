"""Physical definitions and compact native Calypso table handling for S05."""
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np


def write_json(path,value):
    """Close a complete local file before sequential copying to EOS."""
    path=Path(path)
    with tempfile.TemporaryDirectory(prefix='slope-json-') as directory:
        local=Path(directory)/path.name
        local.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
        with local.open('rb') as source,path.open('wb') as target:
            shutil.copyfileobj(source,target,length=1024*1024)
        assert path.stat().st_size==local.stat().st_size


def physical(vectors, rotated):
    values = np.asarray(vectors)
    return np.stack([-values[..., 2], values[..., 1], values[..., 0]], axis=-1) if rotated else values.copy()


def slope(momentum):
    momentum = np.asarray(momentum)
    with np.errstate(divide='ignore', invalid='ignore'):
        return momentum[..., 1]/momentum[..., 2]


def angle(momentum):
    momentum = np.asarray(momentum)
    return np.arctan2(momentum[..., 1], momentum[..., 2])


def histogram(values, edges):
    values = np.asarray(values)
    valid = np.isfinite(values)
    counts, _ = np.histogram(values[valid], edges)
    # ROOT treats the upper edge as overflow; numpy.histogram includes it.
    upper_equal = np.count_nonzero(values[valid] == edges[-1])
    counts[-1] -= upper_equal
    return dict(edges=np.asarray(edges).tolist(), counts=counts.tolist(),
        sumw=counts.tolist(), sumw2=counts.tolist(), raw_entries=len(values),
        finite=int(valid.sum()), invalid=int((~valid).sum()),
        underflow=int(np.count_nonzero(values[valid] < edges[0])),
        overflow=int(np.count_nonzero(values[valid] >= edges[-1])))


def event_first_mask(events, raw_pz_GeV, eligible, lower=50., upper=250.):
    """PDF cut on first stored *long* track; keep all tracks of that event."""
    events = np.asarray(events)
    rows = np.flatnonzero(eligible)
    if not len(rows):
        return np.zeros(len(events), dtype=bool)
    keys, first = np.unique(events[rows], return_index=True)
    pz = raw_pz_GeV[rows[first]]
    passing = keys[(pz > lower) & (pz < upper)]
    return np.isin(events, passing) & eligible


def merge_tables(root, variant):
    """Check exact per-file coverage and merge only successful frozen tasks."""
    root = Path(root)
    tasks = json.loads((root/(variant+'_tasks.json')).read_text())
    manifest = json.loads((root/'manifest.json').read_text())
    rows = {name: [] for name in ['events', 'truth', 'tracks', 'segments']}
    offsets, hit_arrays = [0], []
    seen = set()
    resources = []
    for task in tasks:
        output = Path(task['output'])
        if (output/'exit_code.txt').read_text().strip() != '0':
            raise RuntimeError(f'Unsuccessful task: {output}')
        source = json.loads((output/'source_task.json').read_text())
        summary = json.loads((output/'summary.json').read_text())
        run = json.loads((output/'run.json').read_text())
        assert source['output'] == str(output) and source['variant'] == variant
        assert source['events'] == summary['source_events'] and run['success']
        assert not any(run[k] for k in ['audit', 'reset_event_maps', 'isolate_smoothing'])
        resources.append(dict(output=str(output),wall_s=run['wall_s']))
        for name in rows:
            with np.load(output/(name+'.npz'), allow_pickle=False) as archive:
                data = {key:archive[key] for key in archive.files}
            count = len(data['source_entry'])
            data['file_index'] = np.full(count, source['file_index'], dtype=np.int32)
            data['sample_entry'] = data['source_entry']-source['skip']+source['sample_entry']
            if name == 'events':
                expected = np.arange(source['skip'], source['skip']+source['events'])
                assert np.array_equal(data['source_entry'], expected)
                identities = {(source['file_index'], int(i)) for i in expected}
                assert not seen.intersection(identities)
                seen.update(identities)
            if name == 'tracks':
                lengths = np.diff(data.pop('hit_offsets'))
                offsets.extend((offsets[-1]+np.cumsum(lengths)).tolist())
                hit_arrays.append(data.pop('hit_identifiers'))
            rows[name].append(data)
    assert len(seen) == manifest['source_events_per_variant']
    merged = {name:{key:np.concatenate([p[key] for p in parts]) for key in parts[0]} for name,parts in rows.items()}
    merged['tracks']['hit_offsets'] = np.array(offsets, dtype=np.int64)
    merged['tracks']['hit_identifiers'] = np.concatenate(hit_arrays)
    return merged, dict(variant=variant, source_events=len(seen), tasks=len(tasks),
                        complete_exact_coverage=True, resources=resources)


def join_truth(tables):
    tracks, truth = tables['tracks'], tables['truth']
    keys = ['sample_entry', 'barcode']
    dtype = np.dtype([(key,np.int64) for key in keys])
    def record(data):
        result = np.empty(len(data['barcode']), dtype=dtype)
        for key in keys:
            result[key] = data[key]
        return result
    tkeys, rkeys = record(truth), record(tracks)
    order = np.argsort(tkeys)
    sorted_keys = tkeys[order]
    assert len(np.unique(sorted_keys)) == len(sorted_keys)
    insert = np.searchsorted(sorted_keys, rkeys)
    valid = insert < len(sorted_keys)
    clipped = np.minimum(insert, max(0,len(sorted_keys)-1))
    if len(sorted_keys):
        valid &= sorted_keys[clipped] == rkeys
        indices = order[clipped]
    else:
        indices = np.zeros(len(rkeys), dtype=int)
    return indices, valid
