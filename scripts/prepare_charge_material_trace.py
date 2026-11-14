#!/usr/bin/env python3
"""Select measured transport differences for native material-interaction traces."""
import argparse
import json
from pathlib import Path
import shutil


def dump(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    baseline = root/'physical_prior_transport_v2'
    common = root/'software_material_trace_v3'
    assert (common/'software.json').exists()
    data = json.loads((root/'physical_prior_transport_analysis_v2/summary.json').read_text())
    old = json.loads((baseline/'manifest.json').read_text())
    output = root/'physical_prior_material_trace_v3'
    output.mkdir(exist_ok=False)
    code = output/'code'
    for folder in ['scripts', 'configs']:
        (code/folder).mkdir(parents=True)
    shutil.copy2(repo/'scripts/run_slope_transport.py', code/'scripts/run_slope_transport.py')
    shutil.copy2(repo/'configs/charge_prior_transport.py', code/'configs/charge_prior_transport.py')
    libraries, jobs, samples = {}, {}, {}
    for variant in ['beforeGeo', 'afterGeo']:
        previous = Path(old['libraries'][variant]).parent
        overlay = output/('software_'+variant)
        (overlay/'lib').mkdir(parents=True)
        (overlay/'python').mkdir()
        for p in (previous/'lib').iterdir():
            (overlay/'lib'/p.name).symlink_to(common/'lib'/p.name if 'CalypsoSlopeDiagnostics' in p.name else p.resolve())
        for p in (previous/'python').iterdir():
            (overlay/'python'/p.name).symlink_to(common/'python'/p.name if p.name == 'CalypsoSlopeDiagnostics' else p.resolve())
        libraries[variant] = str(overlay/'lib')
        jobs[variant] = []
    for dsid in old['samples']:
        rows = data['groups'][dsid+'|material_on']
        # Include the largest measured momentum differences at both nontrivial
        # targets plus one deterministic draw from each of the five priors.
        candidates = sorted(set([i*512 for i in range(5)]+[
            r['candidate'] for target in rows[1:] for r in target['largest_relative_p_differences']]))
        group = output/dsid
        group.mkdir()
        for variant in libraries:
            task = json.loads((baseline/dsid/'material_on'/(variant+'_source.json')).read_text())
            selected = [line for line in Path(task['states']).read_text().splitlines(True)
                        if int(line.split()[2]) in candidates]
            assert len(selected) == len(candidates)
            states = group/(variant+'_states.txt')
            states.write_text(''.join(selected))
            task.update(states=str(states), options=['--trace-material'],
                runner=str(code/'configs/charge_prior_transport.py'), output=str(group/variant))
            path = group/(variant+'_source.json')
            dump(path, task)
            jobs[variant].append(dict(kind='analysis', command=['python3',
                str(code/'scripts/run_slope_transport.py'), '--task', str(path)], output=task['output']))
        samples[dsid] = dict(candidate_ids=candidates, input_points=len(candidates),
            selection='Five largest abs(relative-p) differences at VetoNu and first-st1 target, plus draw zero of every seed prior')
    for variant, tasks in jobs.items():
        dump(output/(variant+'_tasks.json'), tasks)
    dump(output/'manifest.json', dict(samples=samples, libraries=libraries,
        source_manifest=str(baseline/'manifest.json'), endpoint_observer=str(common/'software.json'),
        native_overload='Existing IFaserActsExtrapolationTool::propagationSteps(start,target,direction)',
        validation='Normal endpoint output must equal the v2 material-on subset; compare accumulated deltaP to actual endpoint momentum change'))
    print(output, samples)


if __name__ == '__main__':
    main()
