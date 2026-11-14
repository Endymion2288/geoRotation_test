#!/usr/bin/env python3
"""Tabulate the frozen Calypso measurement controls and their actual resources."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np


def quantiles(values):
    data = np.asarray(values, dtype=float)
    assert np.isfinite(data).all()
    return dict(n=len(data), **{name: float(value) for name, value in
        zip(['min', 'median', 'p90', 'max'], np.quantile(data, [0, .5, .9, 1]))}) if len(data) else dict(n=0)


def resources(control):
    rows = []
    names = {'User time (seconds)': 'user_s', 'System time (seconds)': 'system_s',
             'Maximum resident set size (kbytes)': 'rss_KiB'}
    for job in json.loads((control/'all_tasks.json').read_text()):
        folder = Path(job['output'])
        task = json.loads((folder/'source_task.json').read_text())
        row = dict(dsid=task['dsid'], sample_entry=task['sample_entry'], mode=task['coordinate_mode'])
        for line in (folder/'tracking_resource.txt').read_text().splitlines():
            for prefix, key in names.items():
                if line.strip().startswith(prefix+':'):
                    row[key] = float(line.split(':', 1)[1])
        assert all(key in row for key in names.values()), str(folder)
        rows.append(row)
    return dict(jobs=len(rows), total_cpu_s=sum(row['user_s']+row['system_s'] for row in rows),
        per_job_cpu_s=quantiles([row['user_s']+row['system_s'] for row in rows]),
        per_job_rss_KiB=quantiles([row['rss_KiB'] for row in rows]), rows=rows)


def details(control, source, output):
    summary = json.loads((source/'summary.json').read_text())
    events = json.loads((source/'events.json').read_text())
    coordinates = json.loads((source/'coordinates.json').read_text())
    groups = defaultdict(Counter)
    for key, counts in summary['groups'].items():
        dsid, cohort, mode = key.split('|')
        groups[cohort+'|'+mode].update(counts)
    spatial = {}
    for cohort in ['pure_wrong', 'normal']:
        rows = [r for r in coordinates if r['cohort'] == cohort]
        measured = [r for r in rows if 'measured_minus_primary_centroid_mm' in r]
        spatial[cohort] = dict(states=len(rows), primary_sihit_available=len(measured),
            strip_multiplicity=dict(Counter(str(r['n_strips']) for r in rows)),
            channels_without_sdo=sum(r['channels_without_sdo'] for r in rows),
            states_with_channels_without_sdo=sum(r['channels_without_sdo'] > 0 for r in rows),
            states_with_other_sdo=sum(r['channels_with_other'] > 0 for r in rows),
            states_without_primary_sdo=sum(r['channels_with_primary'] == 0 for r in rows),
            absolute_cluster_centroid_mm=quantiles([abs(r['measured_minus_primary_centroid_mm']) for r in measured]),
            cluster_sigma_mm=quantiles([r['cluster_sigma_mm'] for r in rows]),
            primary_endpoint_span_mm=quantiles([r['primary_endpoints_x_range_mm'] for r in measured]),
            connected_component_counts=dict(Counter(str(r['primary_endpoint_connected_components_10nm']) for r in measured)))
    pure = []
    for row in events:
        if row['cohort'] != 'pure_wrong':
            continue
        hits = [r for r in coordinates if r['dsid'] == row['dsid'] and r['sample_entry'] == row['sample_entry']
                and 'measured_minus_primary_centroid_mm' in r]
        result = dict(dsid=row['dsid'], sample_entry=row['sample_entry'], modes={})
        result['largest_coordinate_offset'] = max(hits, key=lambda r: abs(r['measured_minus_primary_centroid_mm']))
        for mode, state in row['modes'].items():
            result['modes'][mode] = {key: state[key] for key in ['status', 'final_qop_GeVinv',
                'same_ordered_hits', 'failed_quality_conditions', 'native_reference_qop_GeVinv',
                'last_update_into_wrong']}
        pure.append(result)
    result = dict(source=str(source), groups={key: dict(v) for key, v in groups.items()},
        spatial=spatial, pure_events=pure,
        actual_measurement_updates_checked=sum(state['updates_checked'] for row in events for state in row['modes'].values()),
        resources=resources(control))
    output.mkdir(parents=True, exist_ok=False)
    (output/'details.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ['pure_events', 'resources']}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--control', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    details(args.control, args.source, args.output)


if __name__ == '__main__':
    main()
