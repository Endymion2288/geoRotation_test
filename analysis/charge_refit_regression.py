#!/usr/bin/env python3
"""Add same-plane angle and observed wall regression to measured refit controls."""
import argparse
import json
from pathlib import Path

import numpy as np

from charge_tables import derive
from slope_tables import histogram, merge_tables, slope, write_json


def angles(tables, decision, truth_mask=None):
    tracks, truth = tables['tracks'], tables['truth']
    selected = decision['selected'] if truth_mask is None else decision['selected'] & truth_mask
    chosen = decision['chosen'][selected]
    truth_rows = np.flatnonzero(selected)
    projection = tracks['projection'][chosen, 0]
    valid = tracks['projection_valid'][chosen, 0] == 1
    reference_valid = truth['truth_model_valid'][truth_rows, 0] == 1
    reference = slope(truth['truth_model_momentum_F_GeV'][truth_rows, 0])
    residual = slope(projection[:, 3:6])-reference
    finite = valid & reference_valid & np.isfinite(residual)
    pull_valid = finite & np.isfinite(projection[:, 6]) & (projection[:, 6] > 0)
    return dict(
        residual_fine=histogram(residual[finite], np.linspace(-.0005, .0005, 201)),
        residual_wide=histogram(residual[finite], np.linspace(-.005, .005, 201)),
        pull=histogram(residual[pull_valid]/np.sqrt(projection[pull_valid, 6]), np.linspace(-20, 20, 201)),
        selected_tracks=len(chosen), finite=int(finite.sum()),
        failed_track_projections=int((~valid).sum()), failed_truth_model=int((~reference_valid).sum()),
        invalid_pull_variance=int(np.count_nonzero(finite & ~pull_valid)),
        quantile_probabilities=[.025, .16, .5, .84, .975],
        residual_quantiles=np.quantile(residual[finite], [.025, .16, .5, .84, .975]).tolist() if finite.any() else None,
        absolute_residual_68=float(np.quantile(np.abs(residual[finite]), .68)) if finite.any() else None,
        absolute_residual_above_1mrad=int(np.count_nonzero(finite & (np.abs(residual) > .001))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    baseline = json.loads((args.root/'refit_analysis_v1/summary.json').read_text())
    groups = {}
    for dsid in ['100116', '100117']:
        native = None
        for mode in ['native', 'smoother', 'updater', 'both']:
            source = (args.root/'mode_analysis_v1/forward_references'/dsid if mode == 'native'
                      else args.root/'refit_control_100k_v1'/dsid/mode)
            tables, coverage = merge_tables(source, 'beforeGeo')
            decision = derive(tables, dsid)
            key = dsid+'|'+mode
            assert int(decision['selected'].sum()) == baseline['groups'][key]['conditional_misid']['total']
            assert int(decision['wrong'].sum()) == baseline['groups'][key]['conditional_misid']['wrong']
            current = angles(tables, decision)
            with np.load(args.root/'refit_analysis_v1'/(key.replace('|', '_')+'_outcomes.npz')) as outcome:
                np.testing.assert_array_equal(outcome['status'], decision['status'])
            groups[key] = dict(angles=current,
                observed_reconstruction_wall_seconds_per_source_event=sum(row['wall_s'] for row in coverage['resources'])/coverage['source_events'],
                coverage=coverage,
                reference='Actual Calypso track and truth-model propagation to z_F=0 mm. Truth model attaches production-p magnitude to the first st1 SiHit chord; this is not local G4 momentum or independently calibrated truth uncertainty.')
            if native is None:
                native = tables, decision
            else:
                old, old_d = native
                for field in ['sample_entry', 'barcode']:
                    np.testing.assert_array_equal(tables['truth'][field], old['truth'][field])
                common = decision['selected'] & old_d['selected']
                same = np.zeros(len(common), dtype=bool)
                for i in np.flatnonzero(common):
                    a, b = old_d['chosen'][i], decision['chosen'][i]
                    left, right = old['tracks'], tables['tracks']
                    if left['ordinal'][a] != right['ordinal'][b]:
                        continue
                    ah = left['hit_identifiers'][left['hit_offsets'][a]:left['hit_offsets'][a+1]]
                    bh = right['hit_identifiers'][right['hit_offsets'][b]:right['hit_offsets'][b+1]]
                    same[i] = np.array_equal(ah, bh)
                groups[key]['common_selected_truth'] = int(common.sum())
                groups[key]['same_candidate_and_ordered_hits'] = int(same.sum())
                groups[key]['common_same_candidate_angles'] = dict(
                    native=angles(old, old_d, same), control=angles(tables, decision, same))
            print(key, current['finite'], current['absolute_residual_68'], flush=True)
            del tables, decision
    write_json(args.output/'summary.json', dict(groups=groups,
        all_source_outcome_status_arrays_exact=True, source_events_per_charge_mode=100000,
        scope='Same 100k per charge and four beforeGeo refit modes as measured charge analysis; isolated smoothing ON',
        wall_interpretation='Different native/control worker batches and startup costs; descriptive observed wall, not a controlled CPU ratio.'))


if __name__ == '__main__':
    main()
