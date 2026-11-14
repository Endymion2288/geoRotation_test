#!/usr/bin/env python3
"""Collect step-control final q/p, resources and source provenance."""
import argparse
import gzip
import json
from pathlib import Path
import re
import subprocess

from charge_component_refit import audit, metrics, qop
from slope_tables import write_json


def compact_tracks(folder):
    event, components = -1, {}
    with gzip.open(folder/'prior_audit.jsonl.gz', 'rt') as stream:
        for line in stream:
            stage = re.search(r'"stage"\s*:\s*"([^"]+)"', line)[1]
            if stage not in ('event', 'selected', 'final'):
                continue
            row = json.loads(line)
            if stage == 'event':
                event += 1
            elif row.get('input_component', -1) >= 0:
                item = components.setdefault((event, row['input_component']), dict(selected={}, final={}))
                item[stage][row['candidate']] = row
    result = {}
    for key, item in components.items():
        selected = min(item['selected'].values(), key=lambda r:(-r['measurements'], r['chi2'], r['candidate']))
        result[key] = dict(selected=qop(selected, 'selected'),
            final=qop(item['final'][selected['candidate']], 'final'))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.input, args.output
    manifest = json.loads((root/'manifest.json').read_text())
    output.mkdir(exist_ok=False)
    tracks, resources = {}, []
    for point in manifest['points']:
        for mode in manifest['modes']:
            folder = Path(point[mode])
            tracks[point['dsid'], point['variant'], mode] = compact_tracks(folder)
            fields = {}
            for prefix, filename in [('native', 'tracking_resource.txt'), ('worker', 'resource.txt')]:
                text = (folder/filename).read_text()
                for label, key in [('User time (seconds)', 'user_cpu_s'), ('System time (seconds)', 'system_cpu_s'),
                                   ('Maximum resident set size (kbytes)', 'max_rss_KiB')]:
                    fields[prefix+'_'+key] = float(re.search(re.escape(label)+r':\s*([\d.]+)', text)[1])
            resources.append(dict(dsid=point['dsid'], variant=point['variant'], mode=mode, **fields))
    comparisons = {}
    for dsid in ('100044', '100116', '100117'):
        pairs = [(mode, 'between_geometry', ('beforeGeo', mode), ('afterGeo', mode)) for mode in manifest['modes']]
        pairs += [(mode, variant+'_step_change', (variant, 'default'), (variant, mode))
                  for variant in ('beforeGeo', 'afterGeo') for mode in ('step100mm', 'step10mm')]
        for mode, kind, a_key, b_key in pairs:
            a, b = tracks[(dsid, *a_key)], tracks[(dsid, *b_key)]
            assert a.keys() == b.keys()
            comparisons[dsid+'|'+kind+'|'+mode] = {stage:metrics([b[k][stage]-a[k][stage] for k in a])
                for stage in ('selected', 'final')}
    extremes = {}
    for row in audit(root/'analysis_v1/paired_states.jsonl.gz'):
        if row['comparison'] == 'between_geometry' and row['dsid'] == '100117' and row['entry'] == 1 and row['component'] == 24:
            extremes.setdefault(row['mode'], []).append(row)
    proof = {}
    for variant, info in manifest['variants'].items():
        source = Path(info['source'])
        tracked = subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True)
        assert tracked == Path(info['tracked_patch']).read_text()
        saved = root/('source_additions_'+variant)
        names = [p.relative_to(saved) for p in saved.rglob('*') if p.is_file()]
        assert all((saved/p).read_bytes() == (source/p).read_bytes() for p in names)
        subprocess.run(['git', '-C', info['source_base'], 'apply', '--check', info['increment']], check=True)
        subprocess.run(['git', '-C', str(source), 'apply', '--reverse', '--check',
                        str(root/(variant+'_complete_source.patch'))], check=True)
        proof[variant] = dict(tracked_diff_exact=True, added_sources_byte_equal=True,
            complete_source_patch_reverse_check=True, step_increment_forward_check=True,
            compiled_manifest=str(root/('software_'+variant)/'software.json'))
    write_json(output/'resource_summary.json', dict(jobs=resources,
        native_cpu_seconds=sum(r['native_user_cpu_s']+r['native_system_cpu_s'] for r in resources),
        worker_cpu_seconds=sum(r['worker_user_cpu_s']+r['worker_system_cpu_s'] for r in resources),
        max_RSS_KiB=max(r['worker_max_rss_KiB'] for r in resources)))
    write_json(output/'qop_comparisons.json', dict(absolute_qop_differences_GeVinv=comparisons,
        scope='First-ranked actual native selected and final tracks; complete common input and selection counts are in analysis_v1.'))
    write_json(output/'extreme_path.json', extremes)
    write_json(output/'source_validation.json', proof)
    print(json.dumps(dict(native_jobs=len(resources), comparisons=len(comparisons))))


if __name__ == '__main__':
    main()
