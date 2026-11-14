#!/usr/bin/env python3
"""Record the native ACTS material-noise terms through the Calypso probe."""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_charge_component_refit import replace_once
from prepare_cluster_measurements import dump


UNIT = Path('Tracking/Acts/FaserActsKalmanFilter/src/CalypsoSlopeDiagnosticsAlg.cxx')


def prepare(root):
    repo = Path(__file__).resolve().parents[1]
    reference = root/'common_prior_ckf_chart_v1/bound_transport_v1'
    output = root/'common_prior_ckf_chart_v1/material_noise_v1'
    output.mkdir(exist_ok=False)
    previous = json.loads((reference/'manifest.json').read_text())
    source = Path(previous['source'])
    old = (source/UNIT).read_text()
    assert old == (reference/'software_probe_v3'/UNIT.name).read_text()
    text = replace_once(old,
        'raw_material_0,raw_material_1,raw_material_2,raw_material_3,raw_material_4\\n',
        'raw_material_0,raw_material_1,raw_material_2,raw_material_3,raw_material_4,'
        'native_sigma_phi2_rad2,native_sigma_theta2_rad2,native_sigma_qop2_GeVinv2\\n')
    text = replace_once(text, "        m_trace<<'\\n';deltaP+=hit.deltaP;",
        "        m_trace<<','<<hit.sigmaPhi2<<','<<hit.sigmaTheta2<<','<<hit.sigmaQoP2;\n"
        "        m_trace<<'\\n';deltaP+=hit.deltaP;")
    increment = repo/'patches/0044_native_recorded_material_noise.patch'
    relative = str(UNIT)
    increment.write_text('diff --git a/'+relative+' b/'+relative+'\n'+''.join(difflib.unified_diff(
        old.splitlines(True), text.splitlines(True), fromfile='a/'+relative, tofile='b/'+relative)))
    (output/'CalypsoSlopeDiagnosticsAlg_before.cxx').write_text(old)
    (source/UNIT).write_text(text)
    (output/UNIT.name).write_text(text)
    subprocess.run(['git', '-C', str(source), 'apply', '--reverse', '--check', str(increment)], check=True)
    code = output/'code'
    for folder, names in {'scripts': ['build_slope_diagnostics.py', 'run_slope_transport.py'],
                          'configs': ['charge_prior_transport.py']}.items():
        (code/folder).mkdir(parents=True)
        for name in names:
            shutil.copy2(repo/folder/name, code/folder/name)
    dump(output/'build_tasks.json', [dict(kind='analysis', command=['python3',
        str(code/'scripts/build_slope_diagnostics.py'), '--source', str(source),
        '--output', str(output/'software_probe')], output=str(output/'jobs/build_probe'))])
    points = []
    for variant in ('beforeGeo', 'afterGeo'):
        jobs = []
        for index, old_job in enumerate(json.loads((reference/(variant+'_tasks_v2.json')).read_text())):
            old_path = Path(old_job['command'][-1])
            old_task = json.loads(old_path.read_text())
            relative_output = Path(old_task['output']).relative_to(reference/'runs_v2')
            destination = output/'runs'/relative_output
            task = dict(old_task, runner=str(code/'configs/charge_prior_transport.py'),
                output=str(destination), reference=str(old_task['output']))
            task_path = output/'tasks'/(variant+'_'+str(index)+'.json')
            task_path.parent.mkdir(exist_ok=True)
            dump(task_path, task)
            jobs.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_transport.py'),
                '--task', str(task_path)], output=str(destination)))
            points.append(dict(variant=variant, output=str(destination), reference=old_task['output']))
        dump(output/(variant+'_tasks.json'), jobs)
    dump(output/'manifest.json', dict(source=str(source), source_base=str(reference/'software_probe_v3'),
        source_commit=previous['source_commit'], increment=str(increment), points=points,
        inputs_unchanged=str(reference/'inputs_v2'), cases=str(reference/'cases.json'),
        native_geometry_and_extrapolation_libraries_unchanged=previous['reference_libraries'],
        scope='Append the actual RecordedMaterial sigmaPhi2/sigmaTheta2/sigmaQoP2 values to the existing '
              'Calypso propagationSteps trace. No tracking or material formula changes. '
              'This is the trace propagation call, not a new observation inside bound propagate. '
              'Complete native transport-array replay and original CSV-column replay are required.'))


def overlay(root):
    output = root/'common_prior_ckf_chart_v1/material_noise_v1'
    reference = root/'common_prior_ckf_chart_v1/bound_transport_v1'
    built = output/'software_probe'
    assert (built/'CalypsoSlopeDiagnosticsAlg.cxx').read_bytes() == (output/UNIT.name).read_bytes()
    for variant in ('beforeGeo', 'afterGeo'):
        target = output/('overlay_'+variant)
        for folder in ('lib', 'python'):
            (target/folder).mkdir(parents=True, exist_ok=False)
            for path in (reference/('overlay_'+variant+'_v2')/folder).iterdir():
                origin = built/folder/path.name if 'CalypsoSlopeDiagnostics' in path.name else path.resolve()
                (target/folder/path.name).symlink_to(origin)


def analysis_tasks(root, version='v1'):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/material_noise_v1'
    code = output/('analysis_code_'+version)
    for relative in ['scripts/validate_charge_components.py']+['analysis/'+name for name in
            ('charge_material_noise.py', 'bound_surface_coordinates.py', 'charge_component_refit.py', 'slope_tables.py')]:
        target = code/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        assert not target.exists()
        shutil.copy2(repo/relative, target)
    for phase, destination in [('gate', 'legacy_gate_'+version), ('analyze', 'analysis_'+version)]:
        command = ['python3',
            str(code/'analysis/charge_material_noise.py'), '--phase', phase, '--input', str(output),
            '--output', str(output/destination)]
        if phase == 'analyze':
            command += ['--gate-validation', str(output/('legacy_gate_'+version)/'validation.json')]
        dump(output/(phase+'_tasks_'+version+'.json'), [dict(kind='analysis', command=command,
            output=str(output/('jobs/'+phase+'_'+version)))])


def figure_tasks(root, version='v2'):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/material_noise_v1'
    code = output/'figure_code_v1'
    shutil.copytree(output/('analysis_code_'+version), code,
        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for name in ('charge_material_noise_figures.py', 'charge_material_noise_artifacts.py',
                 'validate_charge_gun_figures.py'):
        shutil.copy2(repo/'analysis'/name, code/'analysis'/name)
    for phase in ('plot', 'validate'):
        command = ['python3', str(code/'analysis/charge_material_noise_figures.py'), '--phase', phase,
            '--input', str(output/('analysis_'+version)), '--output', str(output/'figures_v1/results')]
        if phase == 'validate':
            command += ['--pdf-tools', str(root.parent.parent/'stage05_geometry_comparison/software_tools')]
        dump(output/(phase+'_tasks_v1.json'), [dict(kind='analysis', command=command,
            output=str(output/('jobs/'+phase+'_v1')))])
    dump(output/'artifact_tasks_v1.json', [dict(kind='analysis', command=['python3',
        str(code/'analysis/charge_material_noise_artifacts.py'), '--input', str(output),
        '--output', str(output/'auxiliary_v1')], output=str(output/'jobs/artifacts_v1'))])
    old = json.loads((root/'common_prior_ckf_chart_v1/ckf_step_control_v1/report_tasks_v1.json').read_text())[0]
    command = old['command'].copy()
    command.insert(command.index('--output'), str(output/'figures_v1/results'))
    command[command.index('--output')+1] = str(root/'report_v20')
    command[command.index('--previous-report')+1] = str(root/'report_v19')
    dump(output/'report_tasks_v1.json', [dict(kind='analysis', command=command,
        output=str(output/'jobs/report_v20'))])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['prepare', 'overlay', 'analysis-tasks', 'figure-tasks'], required=True)
    parser.add_argument('--analysis-version', default='v1')
    args = parser.parse_args()
    if args.phase == 'analysis-tasks':
        analysis_tasks(args.root.resolve(), args.analysis_version)
    elif args.phase == 'figure-tasks':
        figure_tasks(args.root.resolve(), args.analysis_version)
    else:
        {'prepare': prepare, 'overlay': overlay}[args.phase](args.root.resolve())


if __name__ == '__main__':
    main()
