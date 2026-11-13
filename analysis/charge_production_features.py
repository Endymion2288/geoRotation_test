#!/usr/bin/env python3
"""Read production-only muon covariates without conditioning on tracker entrance."""
import argparse
import json
from pathlib import Path

import numpy as np

from slope_tables import write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--outcomes', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    payload, summary = {}, {}
    for dsid in ['100116', '100117']:
        tasks = json.loads((args.root/f'main_{dsid}_v1/afterGeo_tasks.json').read_text())
        blocks, entries = [], []
        for task in tasks:
            output = Path(task['output'])
            assert (output/'exit_code.txt').read_text().strip() == '0'
            source = json.loads((output/'source_task.json').read_text())
            with np.load(output/'truth.npz') as truth:
                p = truth['production_momentum_F_GeV']
                position = truth['production_position_F_mm']
                assert np.all(truth['pdg'] == (-13 if dsid == '100116' else 13))
                assert np.all(truth['status'] == 1)
                assert np.all((truth['barcode'] > 0) & (truth['barcode'] < 200000))
                with np.errstate(divide='ignore', invalid='ignore'):
                    blocks.append(np.column_stack([np.log(np.linalg.norm(p, axis=1)),
                        position[:, 0], position[:, 1], p[:, 0]/p[:, 2], p[:, 1]/p[:, 2]]))
                entries.append(truth['source_entry']-source['skip']+source['sample_entry'])
        features, sample_entry = np.concatenate(blocks), np.concatenate(entries)
        with np.load(args.outcomes/f'{dsid}_afterGeo_outcomes.npz') as outcomes:
            np.testing.assert_array_equal(sample_entry, outcomes['sample_entry'])
            assert len(np.unique(sample_entry)) == len(sample_entry) == 500000
            np.testing.assert_allclose(features[:, 0], np.log(outcomes['truth_p_GeV']), rtol=0, atol=0)
        payload[dsid] = features
        payload[dsid+'_sample_entry'] = sample_entry
        summary[dsid] = dict(source_primary_muons=len(features), finite=int(np.isfinite(features).all(axis=1).sum()),
            mean=np.mean(features, axis=0).tolist(), std=np.std(features, axis=0).tolist())
    np.savez_compressed(args.output/'features.npz', **payload)
    write_json(args.output/'summary.json', dict(samples=summary,
        features=['log production p', 'production x', 'production y', 'production sx', 'production sy'],
        reference='FASER frame at production; no station/hit/direction-after-field inputs',
        population='All 500,000 primary muons per charge, before the truth-fiducial requirement'))


if __name__ == '__main__':
    main()
