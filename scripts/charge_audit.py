"""Retain actual Calypso r/S/K/C on a small, frozen charge replay sample."""
from collections import defaultdict
import gzip
import json

import numpy as np

from audit_reduce import trace_record
from slope_stage_extract import events


def extract(path, output, task):
    entries = 0
    with gzip.open(output/'charge_audit.jsonl.gz', 'wt') as destination:
        for ordinal, records in enumerate(events(path)):
            groups = defaultdict(list)
            summary = []
            for row in records:
                if row['stage'] in ['selected_state', 'selected_shared_state', 'refit_state']:
                    groups[row['stage'], row['seed'], row['candidate']].append(row)
                elif row['stage'] != 'truth_hit':
                    summary.append(row)
            for states in groups.values():
                trace = trace_record(states)
                ordered = list(reversed(states))
                for i, (step, state) in enumerate(zip(trace['steps'], ordered)):
                    step['z_F_mm'] = state['surface'][0 if task['variant'] != 'beforeGeo' else 2][3]
                    for key in ['predicted_cov', 'filtered_cov', 'smoothed_cov', 'jacobian',
                                'projector', 'measurement', 'measurement_cov', 'surface']:
                        if key in state:
                            step[key] = state[key]
                    if 'update_residual_mm' in step:
                        c, h = [np.array(state[k], dtype=float) for k in ['predicted_cov', 'projector']]
                        innovation = h@c@h.T+np.array(state['measurement_cov'])
                        step['update_gain'] = np.linalg.solve(innovation, h@c).T.tolist()
                    if 'smoothing_delta_qop' in step:
                        nxt = ordered[i+1]
                        cf, cp, jac = [np.array(x, dtype=float) for x in
                                       [state['filtered_cov'], nxt['predicted_cov'], nxt['jacobian']]]
                        step['smoothing_gain'] = np.linalg.solve(cp, jac@cf).T.tolist()
                summary.append(trace)
            for row in summary:
                row['source_entry'] = task['skip']+ordinal
                row['sample_entry'] = task['sample_entry']+ordinal
                destination.write(json.dumps(row, allow_nan=False, separators=(',', ':'))+'\n')
            entries += 1
    assert entries == task['events']
    return dict(source_events=entries, artifact='charge_audit.jsonl.gz',
                source='Actual Calypso audit states; r/S/K recomputed from recorded measurement and covariance',
                isolated_states_are_separate_from_native_cached_reference=True)
