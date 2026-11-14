#!/usr/bin/env python3
"""Record runtime resources and fixed-observation provenance for native fits."""
import argparse
import json
from pathlib import Path
import re

import numpy as np

from charge_component_refit import audit, metrics
from charge_forward_covariance import matrix_metrics
from slope_tables import write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.input, args.output
    manifest = json.loads((root/'manifest.json').read_text())
    jobs, comparisons = [], {}
    for point in manifest['points']:
        folder = Path(point['output'])
        fields = {}
        for prefix, filename in [('native', 'tracking_resource.txt'), ('worker', 'resource.txt')]:
            resource = (folder/filename).read_text()
            for label, key in [('User time (seconds)', 'user_cpu_s'), ('System time (seconds)', 'system_cpu_s'),
                               ('Maximum resident set size (kbytes)', 'max_RSS_KiB')]:
                fields[prefix+'_'+key] = float(re.search(re.escape(label)+r':\s*([\d.]+)', resource)[1])
        jobs.append(dict(dsid=point['dsid'], variant=point['variant'], mode=point['mode'], **fields))
    rows = list(audit(args.analysis/'components.jsonl.gz'))
    for dsid in ('100044', '100116', '100117'):
        a = {(r['entry'], r['component']):r for r in rows if r['dsid'] == dsid and r['variant'] == 'beforeGeo'}
        b = {(r['entry'], r['component']):r for r in rows if r['dsid'] == dsid and r['variant'] == 'afterGeo'}
        assert a.keys() == b.keys()
        common = [key for key in a if a[key]['complete'] and b[key]['complete']]
        comparisons[dsid] = dict(input_components=len(a), complete_both=len(common),
            log_likelihood_difference=metrics([b[k]['log_likelihood']-a[k]['log_likelihood'] for k in common]),
            native_float_chi2_sum_difference=metrics([a[k]['innovation_chi2']-a[k]['native_chi2'] for k in common]
                +[b[k]['innovation_chi2']-b[k]['native_chi2'] for k in common]))
    output.mkdir(exist_ok=False)
    write_json(output/'resource_summary.json', dict(jobs=jobs,
        native_cpu_seconds=sum(r['native_user_cpu_s']+r['native_system_cpu_s'] for r in jobs),
        worker_cpu_seconds=sum(r['worker_user_cpu_s']+r['worker_system_cpu_s'] for r in jobs),
        max_RSS_KiB=max(r['worker_max_RSS_KiB'] for r in jobs)))
    write_json(output/'component_comparisons.json', comparisons)
    event_rows = list(audit(args.analysis/'events.jsonl.gz'))
    full_comparisons = []
    for dsid in ('100044', '100116', '100117'):
        indexed = {(r['variant'], r['entry']):r for r in event_rows if r['dsid'] == dsid}
        for entry in range(20):
            a, b = (indexed[v, entry] for v in ('beforeGeo', 'afterGeo'))
            if not (a['complete'] and b['complete']):
                continue
            a, b = a['posterior'], b['posterior']
            ca, cb = np.asarray(a['cov']), np.asarray(b['cov'])
            denominator = np.sqrt(np.maximum((np.diag(ca)+np.diag(cb))/2., 1e-30))
            delta = np.asarray(b['mean'])-np.asarray(a['mean'])
            delta[2] = np.arctan2(np.sin(delta[2]), np.cos(delta[2]))
            full_comparisons.append(dict(dsid=dsid, entry=entry,
                **matrix_metrics(ca, cb), mean_difference=delta.tolist(),
                mean_difference_over_sigma=(delta/denominator).tolist(),
                max_mean_difference_over_sigma=float(np.max(np.abs(delta/denominator)))))
    write_json(output/'full_posterior_comparisons.json', dict(events=full_comparisons,
        scope='Full six-dimensional first-sensor common-chart moments; within and between C retained. '
              'Higher mixture shape differences remain available in all weighted native components. '
              'A moment comparison alone does not prove finite-distribution or sampling convergence.'))
    observations = json.loads((root/'observations.json').read_text())
    proof = {}
    for dsid, record in observations.items():
        path = root/'observations'/(dsid+'_hits.txt')
        lines = [list(map(int, line.split())) for line in path.read_text().splitlines()]
        assert len(lines) == len(record['events']) == 20
        for row, line in zip(record['events'], lines):
            assert line == [row['entry'], len(row['reverse_ordered_hits']), *row['reverse_ordered_hits']]
            assert row['forward_ordered_hits'] == list(reversed(row['reverse_ordered_hits']))
        proof[dsid] = dict(source=record['source'], source_events=20,
            frozen_sequences=sum(bool(r['forward_ordered_hits']) for r in record['events']),
            events_without_observations=[r['entry'] for r in record['events'] if not r['forward_ordered_hits']],
            selector_file_exact=True)
    write_json(output/'observation_validation.json', dict(samples=proof,
        observation_selection='One reference native afterGeo selected cluster sequence per event, frozen before component fits',
        unseeded_events_retained=True, no_component_outcome_used_to_choose_observations=True,
        full_runtime_input_and_library_gates=str(args.analysis/'validation.json')))


if __name__ == '__main__':
    main()
