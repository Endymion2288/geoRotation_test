#!/usr/bin/env python3
"""Aggregate measured GNU time records for the completed S04 matrices."""
import argparse
import json
from pathlib import Path


def wall_seconds(value):
    parts=[float(p) for p in value.split(':')]
    return sum(p*60**i for i,p in enumerate(reversed(parts)))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    groups={
        'discovery_strict':'initial_tasks.json',
        'independent_current_and_strict':'independent_tasks.json',
        'operation_replays_current_and_strict':'operations_tasks.json',
        'diagnostic_simulation_final':'simulation_final_tasks.json',
        'mode_matrix':'modes_tasks.json',
        'case_controls':'case_controls_tasks.json',
        'toy':'toy_tasks.json',
    }
    result={}
    for label,name in groups.items():
        tasks=json.loads((args.root/name).read_text());measurements=[]
        for task in tasks:
            root=Path(task['output'])
            assert (root/'exit_code.txt').read_text().strip()=='0',root
            fields={}
            for line in (root/'resource.txt').read_text().splitlines():
                for key,kind in [('User time (seconds)','user_seconds'),('System time (seconds)','system_seconds'),
                                 ('Elapsed (wall clock) time','wall_seconds'),('Maximum resident set size','max_rss_kB')]:
                    if key in line:
                        value=line.split(': ',1)[1]
                        fields[kind]=wall_seconds(value) if kind=='wall_seconds' else float(value)
            assert len(fields)==4,(root,fields)
            measurements.append(dict(output=str(root),**fields))
        result[label]=dict(tasks=len(tasks),source_events=sum(t.get('events',0) for t in tasks),
            summed_wall_hours=sum(m['wall_seconds'] for m in measurements)/3600,
            summed_cpu_hours=sum(m['user_seconds']+m['system_seconds'] for m in measurements)/3600,
            max_rss_MiB=max(m['max_rss_kB'] for m in measurements)/1024,measurements=measurements)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
