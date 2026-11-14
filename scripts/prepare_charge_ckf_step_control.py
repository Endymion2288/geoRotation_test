#!/usr/bin/env python3
"""Prepare explicit step limits inside the actual native Calypso CKF."""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_charge_component_refit import PACKAGE, provenance, replace_once
from prepare_cluster_measurements import dump


def source(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/ckf_step_control_v1'
    output.mkdir(exist_ok=False)
    old = repo.parent/'geoRotation_calypso_smoothing_audit_beforeGeo'
    target = repo.parent/'geoRotation_calypso_ckf_step_beforeGeo'
    commit = subprocess.check_output(['git', '-C', str(old), 'rev-parse', 'HEAD'], text=True).strip()
    subprocess.run(['git', '-C', str(old), 'worktree', 'add', '--detach', str(target), commit], check=True)
    for name in ('FaserActsGeometry', 'FaserActsGeometryInterfaces', 'FaserActsKalmanFilter'):
        relative = Path('Tracking/Acts')/name
        shutil.copytree(old/relative, target/relative, dirs_exist_ok=True,
            ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    changes = {
        'CKF2.h': replace_once((old/PACKAGE/'CKF2.h').read_text(),
            '  Gaudi::Property<double> m_maxSteps',
            '  Gaudi::Property<double> m_maxStepSize {this,"MaxStepSize",0.,\n'
            '      "Optional explicit CKF propagation step limit in ACTS mm; zero keeps native default"};\n'
            '  Gaudi::Property<double> m_maxSteps'),
        'CKF2.cxx': replace_once((old/PACKAGE/'CKF2.cxx').read_text(),
            '  pOptions.maxSteps = m_maxSteps;',
            '  pOptions.maxSteps = m_maxSteps;\n  if(m_maxStepSize>0.) pOptions.maxStepSize=m_maxStepSize;')}
    patch = ''
    for name, text in changes.items():
        relative = str(PACKAGE/name)
        (target/relative).write_text(text)
        patch += 'diff --git a/'+relative+' b/'+relative+'\n'+''.join(difflib.unified_diff(
            (old/relative).read_text().splitlines(True), text.splitlines(True),
            fromfile='a/'+relative, tofile='b/'+relative))
    increment = repo/'patches/0043_before_native_ckf_step_control.patch'
    increment.write_text(patch)
    subprocess.run(['git', '-C', str(old), 'apply', '--check', str(increment)], check=True)
    subprocess.run(['git', '-C', str(target), 'apply', '--reverse', '--check', str(increment)], check=True)
    tracked = output/'beforeGeo_source.patch'
    tracked.write_text(subprocess.check_output(['git', '-C', str(target), 'diff', 'HEAD'], text=True))
    previous = root/'common_prior_ckf_chart_v1/smoothing_operations_v1/software_beforeGeo'
    manifest = dict(variants=dict(beforeGeo=dict(source=str(target), source_base=str(old), commit=commit,
        tracked_patch=str(tracked), increment=str(increment), geometry_reused_from=str(previous))),
        afterGeo_unchanged_overlay=str(root/'common_prior_ckf_chart_v1/smoothing_operations_v1/overlay_afterGeo'),
        isolated_smoothing=True,
        scope='Explicit step upper limit in the real CKF, default zero preserves the original options. '
              'All existing update/smoothing/material/refit operations retained. Exact default replay required.')
    dump(output/'manifest.json', manifest)
    code = output/'build_code'
    code.mkdir()
    shutil.copy2(repo/'scripts/build_native_calypso.py', code/'build_native_calypso.py')
    dump(output/'build_tasks.json', [dict(kind='analysis', command=['python3', str(code/'build_native_calypso.py'),
        '--source', str(target), '--output', str(output/'software_beforeGeo'), '--jobs', '4',
        '--source-patch', str(tracked), '--reuse-geometry', str(previous)],
        output=str(output/'jobs/build_beforeGeo'))])
    provenance(root, output)


def tasks(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/ckf_step_control_v1'
    original = root/'common_prior_ckf_chart_v1/smoothing_operations_v1'
    manifest = json.loads((output/'manifest.json').read_text())
    code = output/'run_code'
    for relative in ('configs/slope_tracking.py', 'scripts/run_slope_task.py'):
        target = code/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, target)
    (output/'tasks').mkdir()
    modes = dict(default=[], step100mm=['--ckf-max-step-mm', '100'], step10mm=['--ckf-max-step-mm', '10'])
    jobs = {(variant, mode): [] for variant in ('beforeGeo', 'afterGeo') for mode in modes}
    points = []
    for point in json.loads((original/'manifest.json').read_text())['points']:
        dsid, variant = point['dsid'], point['variant']
        reference = Path(point['on'])
        old = json.loads((reference/'source_task.json').read_text())
        saved = dict(dsid=dsid, variant=variant, reference=str(reference))
        for mode, options in modes.items():
            folder = output/'runs'/dsid/variant/mode
            task = dict(old, runner=str(code/'configs/slope_tracking.py'), output=str(folder),
                options=old['options']+options, ckf_step_control=mode, reference=str(reference))
            path = output/'tasks'/(dsid+'_'+variant+'_'+mode+'.json')
            dump(path, task)
            jobs[variant, mode].append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                '--task', str(path)], output=str(folder)))
            saved[mode] = str(folder)
        points.append(saved)
    for (variant, mode), rows in jobs.items():
        dump(output/(variant+'_'+mode+'_tasks.json'), rows)
    for variant in ('beforeGeo', 'afterGeo'):
        dump(output/(variant+'_all_tasks.json'), [row for mode in modes for row in jobs[variant, mode]])
    manifest.update(points=points, modes=modes,
        input_definition='Same complete n32 C/256 physical seed components and analytic after seed-C as the native smoother audit; '
                         'before regular CKF/refit charts, isolated smoothing ON in both versions')
    dump(output/'manifest.json', manifest)


def overlay(root):
    output = root/'common_prior_ckf_chart_v1/ckf_step_control_v1'
    previous = root/'common_prior_ckf_chart_v1/smoothing_operations_v1/overlay_beforeGeo'
    target = output/'overlay_beforeGeo'
    shutil.copytree(previous, target)
    built = output/'software_beforeGeo'
    assert (built/'software.json').exists()
    for folder in ('lib', 'python'):
        shutil.copytree(built/folder, target/folder, dirs_exist_ok=True)


def analysis_tasks(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/ckf_step_control_v1'
    code = output/'analysis_code_v1'
    code.mkdir(exist_ok=False)
    files = {
        'analysis': ('charge_ckf_step_control.py', 'charge_ckf_step_figures.py', 'bound_surface_coordinates.py',
            'charge_component_refit.py', 'charge_forward_covariance.py', 'charge_smoothing_operations.py',
            'charge_refit_input.py', 'slope_tables.py', 'validate_charge_gun_figures.py'),
        'scripts': ('validate_charge_components.py',)}
    for folder, names in files.items():
        (code/folder).mkdir()
        for name in names:
            shutil.copy2(repo/folder/name, code/folder/name)
    runner = code/'analysis/charge_ckf_step_control.py'
    for phase, destination in [('gate', 'default_gate_v1'), ('analyze', 'analysis_v1')]:
        dump(output/(phase+'_tasks_v1.json'), [dict(kind='analysis', command=['python3', str(runner),
            '--phase', phase, '--input', str(output), '--output', str(output/destination)],
            output=str(output/('jobs/'+phase+'_v1')))])
    figures = code/'analysis/charge_ckf_step_figures.py'
    for phase in ('plot', 'validate'):
        command = ['python3', str(figures), '--phase', phase, '--input', str(output/'analysis_v1'),
            '--output', str(output/'figures_v1/results')]
        if phase == 'validate':
            command += ['--pdf-tools', str(root.parent.parent/'stage05_geometry_comparison/software_tools')]
        dump(output/(phase+'_tasks_v1.json'), [dict(kind='analysis', command=command,
            output=str(output/('jobs/'+phase+'_v1')))])


def figure_tasks(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/ckf_step_control_v1'
    code = output/'figure_code_v2'
    shutil.copytree(output/'analysis_code_v1', code,
        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    shutil.copy2(repo/'analysis/charge_ckf_step_figures.py', code/'analysis/charge_ckf_step_figures.py')
    for phase in ('plot', 'validate'):
        command = ['python3', str(code/'analysis/charge_ckf_step_figures.py'), '--phase', phase,
            '--input', str(output/'analysis_v1'), '--output', str(output/'figures_v2/results')]
        if phase == 'validate':
            command += ['--pdf-tools', str(root.parent.parent/'stage05_geometry_comparison/software_tools')]
        dump(output/(phase+'_tasks_v2.json'), [dict(kind='analysis', command=command,
            output=str(output/('jobs/'+phase+'_v2')))])


def artifact_tasks(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/ckf_step_control_v1'
    code = output/'artifact_code_v1'
    shutil.copytree(output/'analysis_code_v1', code,
        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    shutil.copy2(repo/'analysis/charge_ckf_step_artifacts.py', code/'analysis/charge_ckf_step_artifacts.py')
    dump(output/'artifact_tasks_v1.json', [dict(kind='analysis', command=['python3',
        str(code/'analysis/charge_ckf_step_artifacts.py'), '--input', str(output),
        '--output', str(output/'auxiliary_v1')], output=str(output/'jobs/artifacts_v1'))])
    previous = root/'common_prior_ckf_chart_v1/bound_transport_v1/report_tasks_v1.json'
    task = json.loads(previous.read_text())[0]
    command = task['command'].copy()
    command.insert(command.index('--output'), str(output/'figures_v2/results'))
    command[command.index('--output')+1] = str(root/'report_v19')
    command[command.index('--previous-report')+1] = str(root/'report_v18')
    dump(output/'report_tasks_v1.json', [dict(kind='analysis', command=command,
        output=str(output/'jobs/report_v19'))])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['source', 'tasks', 'overlay', 'analysis-tasks', 'figure-tasks', 'artifact-tasks'], required=True)
    args = parser.parse_args()
    {'source': source, 'tasks': tasks, 'overlay': overlay, 'analysis-tasks': analysis_tasks,
     'figure-tasks': figure_tasks, 'artifact-tasks': artifact_tasks}[args.phase](args.root.resolve())


if __name__ == '__main__':
    main()
