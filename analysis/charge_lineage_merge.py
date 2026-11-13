#!/usr/bin/env python3
"""Close full Calypso charge-path reductions against truth-level measurements."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from charge_lineage_analysis import ORIGINS
from slope_tables import histogram, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    baseline = json.loads((args.root/'main_analysis_v1/summary.json').read_text())
    result, examples = {}, {}
    for dsid in ['100116', '100117']:
        charge = 1 if dsid == '100116' else -1
        for variant in ['beforeGeo', 'afterGeo']:
            counts, parts, candidates = Counter(), [], []
            cases, max_conversion, max_barcode = [], 0., 0
            for shard in range(100):
                folder = args.root/'lineage_v1'/dsid/variant/f's{shard:04d}'
                assert (folder/'job/exit_code.txt').read_text().strip() == '0'
                summary = json.loads((folder/'data/summary.json').read_text())
                assert summary['all_final_qop_and_identity_closed']
                counts.update(summary['counts'])
                counts['source_events'] += summary['source_events']
                counts['native_candidates'] += summary['native_candidates']
                counts['candidate_primary_majority_ties'] += summary['candidate_primary_majority_ties']
                max_conversion = max(max_conversion, summary['conversion_qop_max_relative_error'])
                max_barcode = max(max_barcode, summary['barcode_string_max_length'])
                cases.extend(json.loads((folder/'data/replay_events.json').read_text()))
                for name, target in [('finals', parts), ('candidates', candidates)]:
                    with np.load(folder/f'data/{name}.npz') as archive:
                        target.append({k: archive[k] for k in archive.files})
            merged = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
            candidate = {k: np.concatenate([p[k] for p in candidates]) for k in candidates[0]}
            key = dsid+'|'+variant
            expected = baseline['groups'][key]['r95']['conditional_misid']
            assert counts['selected'] == expected['total'] and counts['wrong'] == expected['wrong']
            assert counts['source_events'] == 500000
            with np.load(args.root/f'main_analysis_v1/{dsid}_{variant}_outcomes.npz') as archive:
                outcome = {k: archive[k] for k in archive.files}
            lookup = np.searchsorted(outcome['sample_entry'], merged['sample_entry'])
            np.testing.assert_array_equal(outcome['sample_entry'][lookup], merged['sample_entry'])
            selected, wrong = merged['selected'], merged['wrong']
            np.testing.assert_array_equal(merged['ordinal'][selected], outcome['ordinal'][lookup[selected]])
            np.testing.assert_array_equal(merged['final_qop_GeVinv'][selected], outcome['qop_GeVinv'][lookup[selected]])
            with np.errstate(divide='ignore', invalid='ignore'):
                z = charge*merged['final_qop_GeVinv']/np.sqrt(merged['final_variance_GeVinv2'])
                pull = (merged['final_qop_GeVinv']-charge/outcome['truth_p_GeV'][lookup]) / \
                    np.sqrt(merged['final_variance_GeVinv2'])
            group = dict(counts=dict(counts), max_conversion_relative_error=max_conversion,
                barcode_string_max_length=max_barcode,
                first_error={label: int(np.count_nonzero(selected & (merged['first_wrong'] == i)))
                             for i, label in enumerate(ORIGINS)},
                final_wrong_entry={label: int(np.count_nonzero(wrong & (merged['final_wrong_entry'] == i)))
                                   for i, label in enumerate(ORIGINS)},
                selected_wrong_curvature_significance_gt3=int(np.count_nonzero(wrong & (z < -3))),
                selected_wrong_curvature_significance_gt5=int(np.count_nonzero(wrong & (z < -5))),
                selected_invalid_curvature_variance=int(np.count_nonzero(selected & ~np.isfinite(z))),
                selected_invalid_or_nonpositive_variance=int(np.count_nonzero(selected &
                    (~np.isfinite(merged['final_variance_GeVinv2']) | (merged['final_variance_GeVinv2'] <= 0)))),
                histograms=dict(curvature_significance=histogram(z[selected], np.linspace(-20, 40, 241)),
                    wrong_curvature_significance=histogram(z[wrong], np.linspace(-20, 0, 101)),
                    production_qop_pull=histogram(pull[selected], np.linspace(-20, 20, 201))),
                all_candidate_primary_supported=int(candidate['primary_supported'].sum()),
                all_candidate_primary_supported_ties=int(np.count_nonzero(
                    candidate['primary_supported'] & candidate['primary_majority_tied'])),
                all_selected_decisions_closed_against_full_measurement=True)
            for split in ['discovery', 'held_out']:
                mask = merged['sample_entry'] < 300000 if split == 'discovery' else merged['sample_entry'] >= 300000
                group[split] = dict(selected=int(np.count_nonzero(selected & mask)),
                    wrong=int(np.count_nonzero(wrong & mask)),
                    final_wrong_entry={label: int(np.count_nonzero(wrong & mask & (merged['final_wrong_entry'] == i)))
                                       for i, label in enumerate(ORIGINS)})
            result[key] = group
            examples[key] = cases
            np.savez_compressed(args.output/(dsid+'_'+variant+'_finals.npz'), **merged)
            np.savez_compressed(args.output/(dsid+'_'+variant+'_candidates.npz'), **candidate)
            print(key, group['final_wrong_entry'], flush=True)
    write_json(args.output/'summary.json', dict(groups=result, origin_codes=ORIGINS,
        native_conversion_signed_qop_preserved=True,
        interpretation='First recorded sign error and final entry into wrong sign are distinct; native reference fit includes its own smoothing and target propagation',
        isolated_smoothing='Selected isolated S is observed separately and does not replace the cached native reference mean'))
    write_json(args.output/'replay_events.json', examples)


if __name__ == '__main__':
    main()
