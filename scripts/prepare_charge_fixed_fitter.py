#!/usr/bin/env python3
"""Prepare a real Calypso fixed-source-link fitter control for missing paths."""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_charge_component_refit import PACKAGE, provenance, replace_once
from prepare_cluster_measurements import dump


def source(root, version, resume=False):
    repo = Path(__file__).resolve().parents[1]
    output = root/('common_prior_fixed_fitter_'+version)
    output.mkdir(exist_ok=resume)
    if resume:
        assert not list(output.glob('*_worker/submission.json'))
        assert not list(output.glob('software_*/software.json'))
    (output/'build_code').mkdir(exist_ok=resume)
    shutil.copy2(repo/'scripts/build_native_calypso.py', output/'build_code/build_native_calypso.py')
    chart = root/'common_prior_ckf_chart_v1'
    smoother = json.loads((chart/'smoothing_operations_v1/manifest.json').read_text())
    step = json.loads((chart/'ckf_step_control_v1/manifest.json').read_text())
    references = dict(beforeGeo=step['variants']['beforeGeo'], afterGeo=smoother['variants']['afterGeo'])
    variants, jobs = {}, []
    for variant, info in references.items():
        old = Path(info['source'])
        assert subprocess.check_output(['git', '-C', str(old), 'diff', 'HEAD'], text=True) == Path(info['tracked_patch']).read_text()
        target = repo.parent/('geoRotation_calypso_fixed_fitter_'+version+'_'+variant)
        commit = subprocess.check_output(['git', '-C', str(old), 'rev-parse', 'HEAD'], text=True).strip()
        if resume and (target/'.git').exists():
            assert subprocess.check_output(['git', '-C', str(target), 'rev-parse', 'HEAD'], text=True).strip() == commit
        else:
            subprocess.run(['git', '-C', str(old), 'worktree', 'add', '--detach', str(target), commit], check=True)
        for name in ('FaserActsGeometry', 'FaserActsGeometryInterfaces', 'FaserActsKalmanFilter'):
            relative = Path('Tracking/Acts')/name
            shutil.copytree(old/relative, target/relative, dirs_exist_ok=True,
                ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        changes = {}
        changes['CKF2.h'] = replace_once((old/PACKAGE/'CKF2.h').read_text(),
            '  std::map<size_t,std::vector<SeedComponentInput>> m_seedComponentInputs;',
            '  std::map<size_t,std::vector<SeedComponentInput>> m_seedComponentInputs;\n'
            '  Gaudi::Property<bool> m_componentFixedHitFitter {this,"ComponentFixedHitFitter",false,\n'
            '      "Diagnostic: use the native fixed-source-link KalmanFitter per prior component"};\n'
            '  std::shared_ptr<TrackFitterFunction> m_componentFitter;')
        text = (old/PACKAGE/'CKF2.cxx').read_text()
        anchor = ('  m_finder = makeTrackFinderFunction(trackingGeometry, magneticField,\n'
            '                                  m_resolvePassive, m_resolveMaterial, m_resolveSensitive, *m_logger,\n'
            '                                  m_regularCkfUpdateChart, m_regularCkfSmoothingChart);')
        fitter_call = ('true,true,0.,Acts::FreeToBoundCorrection(false),*m_logger,\n'
            '        false,false,false,m_regularCkfSmoothingChart,m_regularCkfUpdateChart);'
            if variant == 'beforeGeo' else
            'm_scattering,m_energyLoss,0.,Acts::FreeToBoundCorrection(false),*m_logger);')
        if variant == 'afterGeo':
            anchor = ('  m_finder = makeTrackFinderFunction(trackingGeometry, magneticField,\n'
                '                                  m_resolvePassive, m_resolveMaterial, m_resolveSensitive, *m_logger, m_scattering, m_energyLoss);')
        text = replace_once(text, anchor, anchor+'\n'
            '  if(m_componentFixedHitFitter) {\n'
            '    if(m_seedComponentInputFile.value().empty() || m_sourceLinkSelectionFile.value().empty() ||\n'
            '       !m_isolateSmoothing) return StatusCode::FAILURE;\n'
            '    m_componentFitter=makeTrackFitterFunction(trackingGeometry,magneticField,\n'
            '        '+fitter_call+'\n'
            '  }')
        anchor = ('      auto result=(*m_finder)(initial,options,calibrator,slAccessor,\n'
                  '                             Acts::MeasurementSelector(measurementSelectorCfg),tracksTemp);')
        replacement = '''      // This opt-in control calls Calypso's existing fixed-source-link fitter.
      // Its native geometry navigator, material, updates and smoother remain active.
      // It is a measured observation-coverage control, not assumed full coverage.
      auto result=[&]() -> TrackFinderResult {
        if(!m_componentFixedHitFitter)
          return (*m_finder)(initial,options,calibrator,slAccessor,
                            Acts::MeasurementSelector(measurementSelectorCfg),tracksTemp);
        std::vector<Acts::SourceLink> fixedLinks;
        fixedLinks.reserve(slAccessor.container->size());
        m_audit->begin("component_fixed_observations",component.seed);
        m_audit->out << ",\\\"count\\\":" << slAccessor.container->size() << ",\\\"hit_ids\\\":[";
        bool comma=false;
        for(const auto& link:*slAccessor.container) {
          fixedLinks.emplace_back(link);
          if(comma) m_audit->out << ',';
          m_audit->out << link.hit()->identify().get_compact(); comma=true;
        }
        m_audit->out << ']'; m_audit->end();
        auto fitted=(*m_componentFitter)(fixedLinks,initial,options,calibrator,tracksTemp);
        if(!fitted.ok()) return fitted.error();
        return std::vector<FaserActsTrackContainer::TrackProxy>{fitted.value()};
      }();'''
        text = replace_once(text, anchor, replacement)
        anchor = '        m_audit->out << ",\\\"error\\\":" << result.error().value(); m_audit->end();'
        changes['CKF2.cxx'] = replace_once(text, anchor,
            '        m_audit->out << ",\\\"error\\\":" << result.error().value();\n'
            '        if(m_componentFixedHitFitter) m_audit->out << ",\\\"error_category\\\":"\n'
            '            << std::quoted(result.error().category().name()) << ",\\\"error_message\\\":"\n'
            '            << std::quoted(result.error().message());\n'
            '        m_audit->end();')
        increment = ''
        for name, text in changes.items():
            relative = str(PACKAGE/name)
            (target/relative).write_text(text)
            increment += 'diff --git a/'+relative+' b/'+relative+'\n'+''.join(difflib.unified_diff(
                (old/relative).read_text().splitlines(True), text.splitlines(True),
                fromfile='a/'+relative, tofile='b/'+relative))
        patch = repo/'patches'/('0046_native_component_fixed_fitter_'+variant+'.patch')
        patch.write_text(increment)
        subprocess.run(['git', '-C', str(old), 'apply', '--check', str(patch)], check=True)
        subprocess.run(['git', '-C', str(target), 'apply', '--reverse', '--check', str(patch)], check=True)
        tracked = output/(variant+'_tracked_source.patch')
        tracked.write_text(subprocess.check_output(['git', '-C', str(target), 'diff', 'HEAD'], text=True))
        previous = chart/('ckf_step_control_v1' if variant == 'beforeGeo' else 'smoothing_operations_v1')/('software_'+variant)
        variants[variant] = dict(source=str(target), source_base=str(old), commit=commit,
            tracked_patch=str(tracked), increment=str(patch), geometry_reused_from=str(previous))
        jobs.append(dict(kind='analysis', command=['python3', str(output/'build_code/build_native_calypso.py'),
            '--source', str(target), '--output', str(output/('software_'+variant)), '--jobs', '4',
            '--source-patch', str(tracked), '--reuse-geometry', str(previous)],
            output=str(output/('jobs/build_'+variant))))
    dump(output/'manifest.json', dict(variants=variants, isolated_smoothing=True,
        scope='Default-false native fixed-source-link KalmanFitter component control. Actual Calypso '
              'geometry navigator/material/calibrator/update/smoother and conversion/refit retained. '
              'Regular before direction charts; reverse filtering disabled for this diagnostic. '
              'Coverage and exact default replays must be measured before claiming posterior completion.'))
    dump(output/'build_tasks.json', jobs)
    provenance(root, output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--version', default='v2', choices=['v2'])
    parser.add_argument('--phase', choices=['source', 'runs', 'analysis', 'coverage', 'pipeline', 'figures', 'validation'], default='source')
    parser.add_argument('--resume', action='store_true', help='Resume an incomplete preparation before its corresponding jobs are submitted')
    parser.add_argument('--analysis-version', default='v1')
    parser.add_argument('--figure-version', default='v1')
    parser.add_argument('--validation-version', default='v1')
    args = parser.parse_args()
    if args.phase == 'source':
        source(args.root.resolve(), args.version, args.resume)
    elif args.phase == 'runs':
        runs(args.root.resolve(), args.version, args.resume)
    else:
        from prepare_charge_posterior import prepare_analysis, prepare_figures, prepare_validation
        root = args.root.resolve()
        study = 'common_prior_fixed_fitter_'+args.version
        if args.phase == 'analysis':
            prepare_analysis(root, args.analysis_version, study)
        elif args.phase == 'pipeline':
            output = root/study
            code = output/('pipeline_code_'+args.analysis_version)
            shutil.copytree(output/('analysis_code_'+args.analysis_version), code,
                ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
            for name in ('charge_fixed_fitter.py', 'charge_fixed_fitter_pipeline.py'):
                shutil.copy2(Path(__file__).resolve().parents[1]/'analysis'/name, code/'analysis'/name)
            dump(output/('pipeline_tasks_'+args.analysis_version+'.json'), [dict(kind='analysis', command=[
                'python3', str(code/'analysis/charge_fixed_fitter_pipeline.py'), '--input', str(output),
                '--analysis-version', args.analysis_version],
                output=str(output/('jobs/pipeline_'+args.analysis_version)))])
        elif args.phase == 'coverage':
            output = root/study
            code = output/('coverage_code_'+args.analysis_version)
            shutil.copytree(output/('analysis_code_'+args.analysis_version), code,
                ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
            shutil.copy2(Path(__file__).resolve().parents[1]/'analysis/charge_fixed_fitter.py',
                         code/'analysis/charge_fixed_fitter.py')
            dump(output/('coverage_tasks_'+args.analysis_version+'.json'), [dict(kind='analysis', command=[
                'python3', str(code/'analysis/charge_fixed_fitter.py'), '--input', str(output),
                '--analysis', str(output/('analysis_'+args.analysis_version))],
                output=str(output/('jobs/coverage_'+args.analysis_version)))])
        elif args.phase == 'figures':
            prepare_figures(root, args.analysis_version, args.figure_version, study, 'report_v22', 'report_v21')
        else:
            prepare_validation(root, args.analysis_version, args.figure_version, args.validation_version, study)


def runs(root, version, resume=False):
    repo = Path(__file__).resolve().parents[1]
    output = root/('common_prior_fixed_fitter_'+version)
    previous = root/'common_prior_posterior_v1'
    manifest = json.loads((output/'manifest.json').read_text())
    if resume:
        assert all(Path(json.loads(path.read_text())['task_file']).name == 'build_tasks.json'
            for path in output.glob('*_worker/submission.json'))
    (output/'tasks').mkdir(exist_ok=resume)
    code = output/'run_code'
    for relative in ('configs/slope_tracking.py', 'scripts/run_slope_task.py'):
        path = code/relative
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, path)
    assert '--component-fixed-hit-fitter' in (code/'configs/slope_tracking.py').read_text()
    shutil.copytree(previous/'observations', output/'observations', dirs_exist_ok=resume)
    shutil.copy2(previous/'observations.json', output/'observations.json')
    before = root/'common_prior_ckf_chart_v1/ckf_step_control_v1/overlay_beforeGeo'
    after = root/'common_prior_ckf_chart_v1/smoothing_operations_v1/overlay_afterGeo'
    for variant, reference in [('beforeGeo', before), ('afterGeo', after)]:
        built = output/('software_'+variant)
        assert (output/('jobs/build_'+variant)/'exit_code.txt').read_text().strip() == '0'
        info = json.loads((built/'software.json').read_text())
        assert info['commit'] == manifest['variants'][variant]['commit']
        assert info['source_patch'] == Path(manifest['variants'][variant]['tracked_patch']).read_text()
        overlay = output/('overlay_'+variant)
        shutil.copytree(reference, overlay, dirs_exist_ok=resume)
        for folder in ('lib', 'python'):
            shutil.copytree(built/folder, overlay/folder, dirs_exist_ok=True)
    points, tasks = [], {v:[] for v in ('beforeGeo', 'afterGeo')}
    for point in json.loads((previous/'manifest.json').read_text())['points']:
        if point['mode'] != 'fixed_observation':
            continue
        dsid, variant = point['dsid'], point['variant']
        reference = Path(point['output'])
        old = json.loads((reference/'source_task.json').read_text())
        for mode, addition in [('default', []), ('fixed_observation', ['--component-fixed-hit-fitter'])]:
            folder = output/'runs'/dsid/variant/mode
            options = old['options'].copy()
            observation = output/'observations'/(dsid+'_hits.txt')
            assert observation.read_bytes() == Path(old['observation_file']).read_bytes()
            options[options.index('--source-link-selection-file')+1] = str(observation)
            task = dict(old, runner=str(code/'configs/slope_tracking.py'), options=options+addition,
                output=str(folder), reference=str(reference), observation_file=str(observation),
                component_solver='native_fixed_source_link_fitter' if addition else 'native_CKF')
            path = output/'tasks'/(dsid+'_'+variant+'_'+mode+'.json')
            dump(path, task)
            points.append(dict(dsid=dsid, variant=variant, mode=mode, output=str(folder), reference=str(reference)))
            tasks[variant].append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                '--task', str(path)], output=str(folder)))
    for variant, jobs in tasks.items():
        dump(output/(variant+'_tasks.json'), jobs)
    manifest.update(points=points, source_events=60, component_series='matched_n32_f256_v1',
        before_library=str(output/'overlay_beforeGeo/lib'), after_library=str(output/'overlay_afterGeo/lib'),
        component_solver='native_fixed_source_link_fitter')
    dump(output/'manifest.json', manifest)


if __name__ == '__main__':
    main()
