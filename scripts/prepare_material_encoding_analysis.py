#!/usr/bin/env python3
"""Freeze small-replay and complete-sample analysis of one material encoding."""
import argparse
import json
from pathlib import Path
import shutil


def dump(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--revision', default='analysis_v1')
    args = parser.parse_args()
    root = args.root.resolve()
    control = root/'material_encoding_control_v1'
    assert json.loads((control/'smoke_validation.json').read_text())['all_exact']
    manifest = json.loads((control/'manifest.json').read_text())
    output = control/args.revision
    output.mkdir(exist_ok=False)
    code = output/'code'
    shutil.copytree(Path(__file__).resolve().parents[1]/'analysis', code/'analysis',
                    ignore=shutil.ignore_patterns('__pycache__'))
    (code/'scripts').mkdir()
    shutil.copy2(Path(__file__).resolve().parent/'validate_material_encoding_control.py',
                 code/'scripts/validate_material_encoding_control.py')
    sources = {scope: {} for scope in ['smoke', 'full']}
    for dsid, group in manifest['groups'].items():
        original = json.loads((Path(group['source'])/'manifest.json').read_text())
        sources['smoke'][dsid] = {}
        for task in json.loads(Path(group['smoke']).read_text()):
            name = Path(task['output']).name
            view = output/'smoke_views'/dsid/name
            view.mkdir(parents=True)
            dump(view/'afterGeo_tasks.json', [task])
            dump(view/'manifest.json', dict(original, source_events_per_variant=20))
            sources['smoke'][dsid][name] = str(view)
        sources['full'][dsid] = dict(legacy=group['source'], corrected=group['full'])
    plan = dict(sources=sources, material_manifest=str(control/'manifest.json'),
        peak_windows=str(root.parents[1]/'stage05_slopes/main_analysis_v4/results/peak_windows.json'),
        truth_reference='Fixed legacy model in the primary comparison; corrected model change measured separately',
        production_weights=str(root/'production_overlap_v1/truth_weights.npz'))
    dump(output/'plan.json', plan)
    for scope in sources:
        jobs = []
        for dsid in sources[scope]:
            destination = output/scope/dsid
            jobs.append(dict(kind='analysis', command=['python3', str(code/'analysis/material_encoding_reconstruction.py'),
                '--root', str(root), '--plan', str(output/'plan.json'), '--scope', scope,
                '--dsid', dsid, '--output', str(destination)], output=str(output/'jobs'/scope/dsid)))
        dump(output/(scope+'_tasks.json'), jobs)
    dump(output/'full_validation_tasks.json', [dict(kind='analysis', command=['python3',
        str(code/'scripts/validate_material_encoding_control.py'),
        '--root', str(root), '--phase', 'full'], output=str(output/'jobs/full_validation'))])
    print(output)


if __name__ == '__main__':
    main()
