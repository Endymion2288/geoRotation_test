#!/usr/bin/env python3
"""Verify the scope of the typed-surface correction on native audit metrics."""
import argparse
import json
from pathlib import Path

from charge_component_refit import audit
from charge_forward_covariance import STAGES
from slope_tables import write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--corrected', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    old = list(audit(args.previous/'paired_forward_states.jsonl.gz'))
    new = list(audit(args.corrected/'paired_forward_states.jsonl.gz'))
    assert len(old) == len(new) == 14157
    fields = ('direction_scaled_max', 'direction_absolute_max_rad2', 'direction_frobenius_rad2',
              'physical_direction_C_scaled_max', 'direction_difference_norm', 'angular_difference_rad',
              'qop_difference_GeVinv')
    maxima = {field: 0. for field in fields}
    exact = {field: True for field in fields}
    measured = 0
    for a, b in zip(old, new):
        for field in ('dsid', 'entry', 'component', 'forward_index', 'geometry_id', 'hit',
                      'measurement_flag', 'outlier_flag', 'hole_flag', 'material'):
            assert a[field] == b[field], field
        for stage in STAGES:
            for field in fields:
                delta = abs(a[stage][field]-b[stage][field])
                maxima[field] = max(maxima[field], delta)
                exact[field] &= a[stage][field] == b[stage][field]
        measured += b['surface_type'] is not None
    # Coordinate-diagnostic roundoff can depend on worker CPU/BLAS dispatch.
    # The tolerance is on the change of the scaled comparison metric itself.
    assert maxima['direction_scaled_max'] < 1e-12
    assert maxima['physical_direction_C_scaled_max'] < 1e-10
    assert maxima['direction_difference_norm'] < 1e-14
    assert maxima['qop_difference_GeVinv'] == 0.
    write_json(args.output, dict(paired_states=len(new), identifiers_material_and_flags_exact=True,
        position_measured_states=measured, position_unmeasured_states=len(new)-measured,
        direction_metric_absolute_changes=maxima, direction_metrics_bitwise_equal=exact,
        previous=str(args.previous), corrected=str(args.corrected)))
    print(json.dumps(maxima))


if __name__ == '__main__':
    main()
