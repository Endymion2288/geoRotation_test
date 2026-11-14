#!/usr/bin/env python3
"""Freeze common actual measurements for finite-prior native Calypso fits."""
import argparse
import gzip
import json
from pathlib import Path
import shutil

from prepare_cluster_measurements import dump


SAMPLES = ('100044', '100116', '100117')


def observation_rows(folder):
    entry, events = -1, {}
    with gzip.open(folder/'prior_audit.jsonl.gz', 'rt') as stream:
        for line in stream:
            row = json.loads(line)
            if row['stage'] == 'event':
                entry += 1
                events[entry] = dict(native_seeds=row['seeds'], selected=[])
            elif row['stage'] == 'selected':
                events[entry]['selected'].append(row)
    assert entry == 19
    rows = []
    for entry, event in events.items():
        # This defines one common observation for the diagnostic. No finite
        # component is chosen, ranked or removed based on its fit outcome.
        chosen = min(event['selected'], key=lambda r:(-r['measurements'], r['chi2'], r['candidate'])) if event['selected'] else None
        assert (event['native_seeds'] == 0) == (chosen is None)
        hits = chosen['hits'] if chosen else []
        assert len(hits) == len(set(hits))
        rows.append(dict(entry=entry, native_seeds=event['native_seeds'],
            selected_seed=chosen['seed'] if chosen else None,
            reference_candidate=chosen['candidate'] if chosen else None,
            reverse_ordered_hits=hits, forward_ordered_hits=list(reversed(hits))))
    return rows


def prepare(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_posterior_v1'
    output.mkdir(exist_ok=False)
    (output/'tasks').mkdir()
    code = output/'run_code'
    for relative in ('configs/slope_tracking.py', 'scripts/run_slope_task.py'):
        target = code/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, target)
    observations, points, jobs = {}, [], {v:[] for v in ('beforeGeo', 'afterGeo')}
    step = root/'common_prior_ckf_chart_v1/ckf_step_control_v1'
    for dsid in SAMPLES:
        reference = root/'common_prior_ckf_v1/smoke/afterGeo'/dsid/'audit'
        observations[dsid] = dict(source=str(reference), events=observation_rows(reference))
        hits_file = output/'observations'/(dsid+'_hits.txt')
        hits_file.parent.mkdir(exist_ok=True)
        hits_file.write_text(''.join(' '.join(map(str, [r['entry'], len(r['reverse_ordered_hits']),
            *r['reverse_ordered_hits']]))+'\n' for r in observations[dsid]['events']))
        for variant in ('beforeGeo', 'afterGeo'):
            reference = step/'runs'/dsid/variant/'step10mm'
            old = json.loads((reference/'source_task.json').read_text())
            for mode in ('default', 'fixed_observation'):
                folder = output/'smoke'/dsid/variant/mode
                options = old['options'].copy()
                if mode == 'fixed_observation':
                    options += ['--source-link-selection-file', str(hits_file),
                                '--ckf-chi2-max', '1e30', '--ckf-nmax', '1']
                task = dict(old, runner=str(code/'configs/slope_tracking.py'), output=str(folder), options=options,
                    posterior_control=mode, reference=str(reference), observation_file=str(hits_file))
                task_path = output/'tasks'/(dsid+'_'+variant+'_'+mode+'.json')
                dump(task_path, task)
                jobs[variant].append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                    '--task', str(task_path)], output=str(folder)))
                points.append(dict(dsid=dsid, variant=variant, mode=mode, output=str(folder), reference=str(reference)))
    for variant, tasks in jobs.items():
        dump(output/(variant+'_smoke_tasks.json'), tasks)
    dump(output/'observations.json', observations)
    dump(output/'manifest.json', dict(points=points, reference_step_control=str(step),
        before_library=str(step/'overlay_beforeGeo/lib'),
        after_library=str(root/'common_prior_ckf_chart_v1/smoothing_operations_v1/overlay_afterGeo/lib'),
        source_events=60, isolated_smoothing=True, component_series='matched_n32_f256_v1',
        scope='Same existing finite physical-prior components. One frozen actual cluster measurement sequence '
              'per source event, defined before component fits. Native selector chi2Max=1e30 and nMax=1 '
              'are explicit diagnostics, not production settings. All inputs, failures, incomplete '
              'observation paths and refit fallback remain in coverage; no posterior is claimed if '
              'required likelihoods are missing. Six default replays must reproduce archived outputs.'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['runs', 'analysis', 'figures', 'validation'], default='runs')
    parser.add_argument('--analysis-version', default='v1')
    parser.add_argument('--figure-version', default='v1')
    parser.add_argument('--validation-version', default='v1')
    args = parser.parse_args()
    if args.phase == 'runs':
        prepare(args.root.resolve())
    elif args.phase == 'analysis':
        prepare_analysis(args.root.resolve(), args.analysis_version)
    elif args.phase == 'figures':
        prepare_figures(args.root.resolve(), args.analysis_version, args.figure_version)
    else:
        prepare_validation(args.root.resolve(), args.analysis_version, args.figure_version, args.validation_version)


def prepare_analysis(root, version, study_name='common_prior_posterior_v1'):
    repo = Path(__file__).resolve().parents[1]
    output = root/study_name
    code = output/('analysis_code_'+version)
    code.mkdir(exist_ok=False)
    files = {
        'analysis': ('charge_posterior.py', 'bound_surface_coordinates.py',
            'charge_component_refit.py', 'charge_forward_covariance.py', 'charge_smoothing_operations.py',
            'charge_refit_input.py', 'slope_tables.py'),
        'scripts': ('validate_charge_components.py',)}
    for folder, names in files.items():
        (code/folder).mkdir()
        for name in names:
            shutil.copy2(repo/folder/name, code/folder/name)
    runner = code/'analysis/charge_posterior.py'
    for phase, destination in [('gate', 'default_gate_'+version), ('analyze', 'analysis_'+version)]:
        command = ['python3', str(runner), '--phase', phase, '--input', str(output),
            '--output', str(output/destination)]
        if phase == 'analyze':
            command += ['--gate-validation', str(output/('default_gate_'+version+'/validation.json'))]
        dump(output/(phase+'_tasks_'+version+'.json'), [dict(kind='analysis', command=command,
            output=str(output/('jobs/'+phase+'_'+version)))])


def prepare_figures(root, version, figure_version, study_name='common_prior_posterior_v1',
                    report_name='report_v21', previous_report='report_v20'):
    repo = Path(__file__).resolve().parents[1]
    output = root/study_name
    code = output/('figure_code_'+figure_version)
    shutil.copytree(output/('analysis_code_'+version), code,
        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for name in ('charge_posterior_figures.py', 'charge_posterior_artifacts.py', 'validate_charge_gun_figures.py'):
        shutil.copy2(repo/'analysis'/name, code/'analysis'/name)
    for phase in ('plot', 'validate'):
        command = ['python3', str(code/'analysis/charge_posterior_figures.py'), '--phase', phase,
            '--input', str(output/('analysis_'+version)), '--output', str(output/('figures_'+figure_version+'/results'))]
        if phase == 'validate':
            command += ['--pdf-tools', str(root.parent.parent/'stage05_geometry_comparison/software_tools')]
        dump(output/(phase+'_tasks_'+figure_version+'.json'), [dict(kind='analysis', command=command,
            output=str(output/('jobs/'+phase+'_'+figure_version)))])
    dump(output/('artifact_tasks_'+figure_version+'.json'), [dict(kind='analysis', command=['python3',
        str(code/'analysis/charge_posterior_artifacts.py'), '--input', str(output),
        '--analysis', str(output/('analysis_'+version)), '--output', str(output/('auxiliary_'+figure_version))],
        output=str(output/('jobs/artifacts_'+figure_version)))])
    report_task = (root/'common_prior_posterior_v1/report_tasks_v3.json'
        if study_name.startswith('common_prior_fixed_fitter_') else
        root/'common_prior_ckf_chart_v1/material_noise_v1/report_tasks_v1.json')
    old = json.loads(report_task.read_text())[0]
    command = old['command'].copy()
    command.insert(command.index('--output'), str(output/('figures_'+figure_version+'/results')))
    command[command.index('--output')+1] = str(root/report_name)
    command[command.index('--previous-report')+1] = str(root/previous_report)
    dump(output/('report_tasks_'+figure_version+'.json'), [dict(kind='analysis', command=command,
        output=str(output/('jobs/'+report_name)))])


def prepare_validation(root, analysis_version, figure_version, validation_version,
                       study_name='common_prior_posterior_v1'):
    repo = Path(__file__).resolve().parents[1]
    output = root/study_name
    code = output/('validation_code_'+validation_version)
    shutil.copytree(output/('figure_code_'+figure_version), code,
        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    shutil.copy2(repo/'analysis/charge_posterior_figures.py', code/'analysis/charge_posterior_figures.py')
    dump(output/('validate_tasks_'+validation_version+'.json'), [dict(kind='analysis', command=['python3',
        str(code/'analysis/charge_posterior_figures.py'), '--phase', 'validate',
        '--input', str(output/('analysis_'+analysis_version)),
        '--output', str(output/('figures_'+figure_version+'/results')),
        '--pdf-tools', str(root.parent.parent/'stage05_geometry_comparison/software_tools')],
        output=str(output/('jobs/validate_'+validation_version)))])


if __name__ == '__main__':
    main()
