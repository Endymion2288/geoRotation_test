#!/usr/bin/env python3
"""Read actual Calypso process CPU/RSS separately from event-loop wall time."""
import argparse
import json
from pathlib import Path
import re

from slope_tables import write_json


def resources(jobs):
    records = []
    for job in jobs:
        folder = Path(job['output'])
        task = json.loads((folder/'source_task.json').read_text())
        run = json.loads((folder/'run.json').read_text())
        assert run['success'] and run['isolate_smoothing']
        text = (folder/'tracking_resource.txt').read_text()
        def field(name):
            match = re.search(r'^\s*'+re.escape(name)+r':\s*(\S+)\s*$', text, re.MULTILINE)
            assert match, name
            return float(match[1])
        assert field('Exit status') == 0
        records.append(dict(output=str(folder), source_events=task['events'],
            process_cpu_seconds=field('User time (seconds)')+field('System time (seconds)'),
            process_maximum_RSS_KiB=field('Maximum resident set size (kbytes)'),
            event_loop_wall_seconds=run['wall_s']))
    count = sum(row['source_events'] for row in records)
    return dict(source_events=count, processes=len(records),
        process_cpu_ms_per_source_event=1000*sum(row['process_cpu_seconds'] for row in records)/count,
        event_loop_wall_ms_per_source_event=1000*sum(row['event_loop_wall_seconds'] for row in records)/count,
        maximum_process_RSS_KiB=max(row['process_maximum_RSS_KiB'] for row in records), records=records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    refit, gun = {}, {}
    for dsid in ['100116', '100117']:
        for mode in ['native', 'smoother', 'updater', 'both']:
            source = (args.root/'mode_analysis_v1/forward_references'/dsid if mode == 'native'
                      else args.root/'refit_control_100k_v1'/dsid/mode)
            jobs = json.loads((source/'beforeGeo_tasks.json').read_text())
            refit[dsid+'|'+mode] = resources(jobs)
    manifest = json.loads((args.root/'symmetric_gun_20k_v1/manifest.json').read_text())
    for dsid in ['100116', '100117']:
        for momentum in [100, 1000]:
            sources = [task for task in manifest['source_tasks'] if task['dsid'] == dsid and task['p_GeV'] == momentum]
            for mode in ['beforeGeo', 'beforeGeo_regular_both', 'afterGeo']:
                gun[f'{dsid}|{momentum}|{mode}'] = resources([dict(output=str(Path(task['output'])/mode)) for task in sources])
    write_json(args.output/'summary.json', dict(refit=refit, gun=gun,
        definitions=dict(process_cpu='GNU time user+system for complete tracking process, including initialization',
            wall='acc.run wall including Athena initialization/finalization, excluding prior Python configuration; separate from whole-process elapsed time',
            RSS='Maximum resident set size observed by GNU time; KiB'),
        interpretation='Recorded worker measurements. Batches ran on different workers and at different times; CPU/event is descriptive and does not establish a controlled speed ratio.'))
    print(json.dumps({group: {key: value['process_cpu_ms_per_source_event'] for key, value in rows.items()}
                      for group, rows in [('refit', refit), ('gun', gun)]}, indent=2))


if __name__ == '__main__':
    main()
