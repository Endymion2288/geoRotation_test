#!/usr/bin/env python3
"""Add full-C surface-to-surface inputs to the actual Calypso transport probe."""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

from prepare_charge_component_refit import replace_once
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from bound_surface_coordinates import rotate_tangent


UNIT = Path('Tracking/Acts/FaserActsKalmanFilter/src/CalypsoSlopeDiagnosticsAlg.cxx')


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def source(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/bound_transport_v1'
    output.mkdir(exist_ok=False)
    old = repo.parent/'geoRotation_calypso_material_bins_v4'
    target = repo.parent/'geoRotation_calypso_bound_transport_v1'
    commit = subprocess.check_output(['git', '-C', str(old), 'rev-parse', 'HEAD'], text=True).strip()
    subprocess.run(['git', '-C', str(old), 'worktree', 'add', '--detach', str(target), commit], check=True)
    previous = (old/UNIT).read_text()
    text = replace_once(previous,
        'double m_outputParameters[3][6]{},m_outputCovariance[3][6][6]{};',
        'double m_outputParameters[3][6]{},m_outputCovariance[3][6][6]{};\n'
        '  double m_outputBound[3][6]{},m_inputCovariance[6][6]{},m_inputSurface[4][4]{},m_outputSurface[3][4][4]{};')
    text = replace_once(text,
        'Acts::BoundVector values;};',
        'Acts::BoundVector values;Acts::BoundSquareMatrix covariance;unsigned long long targetGeometry=0;};')
    text = replace_once(text,
        'const Acts::BoundTrackParameters&,const std::array<double,3>&);',
        'const Acts::BoundTrackParameters&,const std::array<double,3>&,const Acts::Surface* exactTarget=nullptr);')
    text = replace_once(text,
        '      m_inputs[entry].push_back(row);',
        '      if(row.kind==2) {\n'
        '        for(int r=0;r<6;++r) for(int c=0;c<6;++c) if(!(file>>row.covariance(r,c))) return StatusCode::FAILURE;\n'
        '        if(!(file>>row.targetGeometry)) return StatusCode::FAILURE;\n'
        '      }\n'
        '      m_inputs[entry].push_back(row);')
    text = replace_once(text,
        '  m_tree->Branch("projection_valid",m_valid,"projection_valid[3]/I");',
        '  m_tree->Branch("projection_valid",m_valid,"projection_valid[3]/I");\n'
        '  m_tree->Branch("input_covariance_native",m_inputCovariance,"input_covariance_native[6][6]/D");\n'
        '  m_tree->Branch("input_surface_native",m_inputSurface,"input_surface_native[4][4]/D");\n'
        '  m_tree->Branch("projection_bound_native",m_outputBound,"projection_bound_native[3][6]/D");\n'
        '  m_tree->Branch("projection_surface_native",m_outputSurface,"projection_surface_native[3][4][4]/D");')
    text = replace_once(text,
        'const Acts::BoundTrackParameters& start,const std::array<double,3>& targets) {',
        'const Acts::BoundTrackParameters& start,const std::array<double,3>& targets,const Acts::Surface* exactTarget) {')
    text = replace_once(text,
        '  m_hasCov=start.covariance().has_value();m_charge=start.charge();m_inputQop=start.parameters()[4];',
        '  m_hasCov=start.covariance().has_value();m_charge=start.charge();m_inputQop=start.parameters()[4];\n'
        '  for(int r=0;r<6;++r) for(int c=0;c<6;++c) m_inputCovariance[r][c]=m_hasCov ? (*start.covariance())(r,c) : NAN;\n'
        '  const auto& inputTransform=start.referenceSurface().transform(gctx).matrix();\n'
        '  for(int r=0;r<4;++r) for(int c=0;c<4;++c) m_inputSurface[r][c]=inputTransform(r,c);')
    text = replace_once(text,
        '    if(!std::isfinite(targets[i])) {m_valid[i]=-1;continue;}\n'
        '    auto target=plane(targets[i]);',
        '    for(double& v:m_outputBound[i]) v=NAN;\n'
        '    for(auto& row:m_outputSurface[i]) for(double& v:row) v=NAN;\n'
        '    if(!std::isfinite(targets[i])) {m_valid[i]=-1;continue;}\n'
        '    std::shared_ptr<const Acts::Surface> target=(i==0 && exactTarget) ? exactTarget->getSharedPtr() : plane(targets[i]);')
    text = replace_once(text,
        '    m_outputQop[i]=result->parameters()[4];m_outputCharge[i]=result->charge();',
        '    m_outputQop[i]=result->parameters()[4];m_outputCharge[i]=result->charge();\n'
        '    for(int j=0;j<6;++j) m_outputBound[i][j]=result->parameters()[j];\n'
        '    const auto& outputTransform=result->referenceSurface().transform(gctx).matrix();\n'
        '    for(int r=0;r<4;++r) for(int c=0;c<4;++c) m_outputSurface[i][r][c]=outputTransform(r,c);')
    text = replace_once(text, '      if(row.kind==1) {', '      if(row.kind==1 || row.kind==2) {')
    text = replace_once(text,
        '      Acts::BoundTrackParameters start(std::move(surface),pars,std::nullopt,Acts::ParticleHypothesis::muon());\n'
        '      propagate(ctx,gctx,start,{0.,-3112.,first.count(row.barcode) ? first.at(row.barcode) : NAN});',
        '      std::optional<Acts::BoundSquareMatrix> covariance;\n'
        '      if(row.kind==2) covariance=row.covariance;\n'
        '      Acts::BoundTrackParameters start(std::move(surface),pars,covariance,Acts::ParticleHypothesis::muon());\n'
        '      if(row.kind==2) {\n'
        '        const auto* exact=tracking->findSurface(Acts::GeometryIdentifier(row.targetGeometry));\n'
        '        if(!exact) return StatusCode::FAILURE;\n'
        '        propagate(ctx,gctx,start,{toF(exact->center(gctx)).z(),NAN,NAN},exact);\n'
        '      } else propagate(ctx,gctx,start,{0.,-3112.,first.count(row.barcode) ? first.at(row.barcode) : NAN});')
    (target/UNIT).write_text(text)
    relative = str(UNIT)
    patch = 'diff --git a/'+relative+' b/'+relative+'\n'+''.join(difflib.unified_diff(
        previous.splitlines(True), text.splitlines(True), fromfile='a/'+relative, tofile='b/'+relative))
    increment = repo/'patches/0040_full_covariance_native_transport_probe.patch'
    increment.write_text(patch)
    subprocess.run(['git', '-C', str(old), 'apply', '--check', str(increment)], check=True)
    subprocess.run(['git', '-C', str(target), 'apply', '--reverse', '--check', str(increment)], check=True)
    code = output/'code'
    (code/'scripts').mkdir(parents=True)
    shutil.copy2(repo/'scripts/build_slope_diagnostics.py', code/'scripts/build_slope_diagnostics.py')
    jobs = [dict(kind='analysis', command=['python3', str(code/'scripts/build_slope_diagnostics.py'),
        '--source', str(target), '--output', str(output/'software_probe')], output=str(output/'jobs/build_probe'))]
    dump(output/'build_tasks.json', jobs)
    dump(output/'manifest.json', dict(source=str(target), source_base=str(old), source_commit=commit,
        source_file=str(UNIT), source_increment=str(increment), algorithm_source=str(output/'CalypsoSlopeDiagnosticsAlg.cxx'),
        scope='Actual Calypso extrapolation tool, complete bound mean/C and actual native source/target surfaces. '
              'Probe extension only; native transport and CKF algorithms unchanged.'))
    shutil.copy2(target/UNIT, output/'CalypsoSlopeDiagnosticsAlg.cxx')
    print(output)


def tasks(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/bound_transport_v1'
    code = output/'run_code'
    for folder, names in {'scripts': ['run_slope_transport.py'],
                          'configs': ['charge_prior_transport.py']}.items():
        (code/folder).mkdir(parents=True, exist_ok=False)
        for name in names:
            shutil.copy2(repo/folder/name, code/folder/name)
    (output/'inputs').mkdir()
    (output/'tasks').mkdir()
    previous = json.loads((root/'physical_prior_material_bins_v4/manifest.json').read_text())
    audit = root/'common_prior_ckf_chart_v1/smoothing_operations_v1'
    paths = json.loads((root/'common_prior_ckf_chart_v1/forward_covariance_v1/analysis_v2/extreme_forward_paths.json').read_text())
    lines = {(dsid, variant): [] for dsid in ('100044', '100116', '100117')
             for variant in ('beforeGeo', 'afterGeo')}
    cases = []
    for path in paths:
        dsid = path['dsid']
        for ordinal in range(1, len(path['rows'])):
            a, b = path['before_states'][ordinal-1], path['after_states'][ordinal-1]
            targets = {v: path[prefix+'_states'][ordinal] for v, prefix in
                       [('beforeGeo', 'before'), ('afterGeo', 'after')]}
            mean_f, covariance_f = np.asarray(a['filtered']), np.asarray(a['filtered_cov'])
            mapped, mapped_covariance = rotate_tangent(mean_f, covariance_f,
                np.asarray(a['surface']), np.asarray(b['surface']))
            for mode in ('native_inputs', 'fixed_before_input'):
                candidate = len(lines[dsid, 'beforeGeo'])
                cases.append(dict(dsid=dsid, entry=path['entry'], component=path['component'],
                    candidate=candidate, mode=mode, from_forward_index=ordinal-1, to_forward_index=ordinal,
                    source_geometry_id=a['geometry_id'], target_geometry_id=targets['beforeGeo']['geometry_id'],
                    source_path_length_mm=a['path_length_mm'], target_path_length_mm=targets['beforeGeo']['path_length_mm']))
                for variant in ('beforeGeo', 'afterGeo'):
                    row = a if variant == 'beforeGeo' else b
                    mean, cov = np.asarray(row['filtered']), np.asarray(row['filtered_cov'])
                    if variant == 'afterGeo' and mode == 'fixed_before_input':
                        mean, cov = mapped, mapped_covariance
                    fields = [path['entry'], path['before_seed']['seed'], candidate, 0, 0, 0, 2,
                              row['geometry_id'], int(np.sign(mean[4]))]
                    values = list(mean)+list(cov.ravel())
                    lines[dsid, variant].append(' '.join(map(str, fields))+' '+
                        ' '.join(format(value, '.17g') for value in values)+' '+str(targets[variant]['geometry_id'])+'\n')
    dump(output/'cases.json', cases)
    jobs = {v: [] for v in ('beforeGeo', 'afterGeo')}
    replay_jobs = {v: [] for v in jobs}
    controls = dict(native=['--trace-material'], material_off=['--disable-material'],
                    step100mm=['--trace-material', '--max-step-size-m', '0.1'],
                    step10mm=['--trace-material', '--max-step-size-m', '0.01'])
    for dsid in ('100044', '100116', '100117'):
        for variant in jobs:
            state_file = output/'inputs'/(dsid+'_'+variant+'.txt')
            state_file.write_text(''.join(lines[dsid, variant]))
            for control, options in controls.items():
                source_task = json.loads((audit/'runs'/dsid/variant/'on/source_task.json').read_text())
                task = dict(dsid=dsid, variant=variant, input=source_task['input'], material=source_task['material'],
                    events=20, skip=0, states=str(state_file), runner=str(code/'configs/charge_prior_transport.py'),
                    options=options, output=str(output/'runs'/dsid/variant/control))
                path = output/'tasks'/(dsid+'_'+variant+'_'+control+'.json')
                dump(path, task)
                jobs[variant].append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_transport.py'),
                    '--task', str(path)], output=task['output']))
    for variant in jobs:
        dump(output/(variant+'_tasks.json'), jobs[variant])
        case = 'beforeGeo_trace' if variant == 'beforeGeo' else 'afterGeo_corrected_trace'
        old = json.loads(Path(previous['source_tasks']['100117|'+case]).read_text())
        task = dict(old, runner=str(code/'configs/charge_prior_transport.py'),
                    output=str(output/'legacy_replay'/variant), reference=old['output'])
        path = output/'tasks'/('legacy_'+variant+'.json')
        dump(path, task)
        replay_jobs[variant].append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_transport.py'),
            '--task', str(path)], output=task['output']))
        dump(output/(variant+'_replay_tasks.json'), replay_jobs[variant])
    manifest = json.loads((output/'manifest.json').read_text())
    manifest.update(cases=str(output/'cases.json'), controls=controls,
        source_paths=str(root/'common_prior_ckf_chart_v1/forward_covariance_v1/analysis_v2/extreme_forward_paths.json'),
        inputs_from_isolated_smoothing=True, complete_native_input_rows={dsid: len(lines[dsid, 'beforeGeo']) for dsid in ('100044', '100116', '100117')},
        reference_libraries=previous['libraries'],
        interpretation='Transport-only local closure experiment. Actual native source/target surfaces and full input covariance; '
                       'fixed-before maps the same complete tangent Gaussian with an analytic Jacobian. '
                       'Material controls and step limits are explicit; not a CKF posterior or finite-distribution equivalence claim.')
    dump(output/'manifest.json', manifest)
    print(json.dumps(manifest['complete_native_input_rows']))


def overlay(root):
    output = root/'common_prior_ckf_chart_v1/bound_transport_v1'
    manifest = json.loads((output/'manifest.json').read_text())
    built = output/'software_probe'
    assert (built/'software.json').is_file()
    for variant, library in manifest['reference_libraries'].items():
        previous = Path(library).parent
        target = output/('overlay_'+variant)
        (target/'lib').mkdir(parents=True, exist_ok=False)
        (target/'python').mkdir()
        for path in (previous/'lib').iterdir():
            origin = built/'lib'/path.name if 'CalypsoSlopeDiagnostics' in path.name else path.resolve()
            (target/'lib'/path.name).symlink_to(origin)
        for path in (previous/'python').iterdir():
            origin = built/'python'/path.name if path.name == 'CalypsoSlopeDiagnostics' else path.resolve()
            (target/'python'/path.name).symlink_to(origin)
    assert (built/'CalypsoSlopeDiagnosticsAlg.cxx').read_bytes() == Path(manifest['algorithm_source']).read_bytes()
    dump(output/'source_validation.json', dict(observer_source_exact=True,
        compiled_component=json.loads((built/'software.json').read_text()),
        before_after_geometry_and_field_libraries_reused=manifest['reference_libraries'],
        native_extrapolation_tool_unchanged=True, legacy_replay_required=True))


def passive_surfaces(root):
    """Find actual approach/boundary surfaces outside ACTS' sensitive index."""
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/bound_transport_v1'
    manifest = json.loads((output/'manifest.json').read_text())
    source = Path(manifest['source'])
    previous = (source/UNIT).read_text()
    text = replace_once(previous, '#include "Acts/Geometry/TrackingVolume.hpp"',
        '#include "Acts/Geometry/TrackingVolume.hpp"\n#include "Acts/Geometry/ApproachDescriptor.hpp"')
    finder = '''  const Acts::Surface* completeSurface(const Acts::TrackingVolume& volume,unsigned long long id) const {
    for(const auto& boundary:volume.boundarySurfaces()) {
      const auto& surface=boundary->surfaceRepresentation();
      if(surface.geometryId().value()==id) return &surface;
    }
    if(volume.confinedLayers()) for(const auto& layer:volume.confinedLayers()->arrayObjects()) {
      const auto& surface=layer->surfaceRepresentation();
      if(surface.geometryId().value()==id) return &surface;
      if(layer->approachDescriptor()) for(const auto* approach:layer->approachDescriptor()->containedSurfaces())
        if(approach->geometryId().value()==id) return approach;
      if(layer->surfaceArray()) for(const auto* sensitive:layer->surfaceArray()->surfaces())
        if(sensitive->geometryId().value()==id) return sensitive;
    }
    if(volume.confinedVolumes()) for(const auto& child:volume.confinedVolumes()->arrayObjects())
      if(const auto* found=completeSurface(*child,id)) return found;
    for(const auto& child:volume.denseVolumes())
      if(const auto* found=completeSurface(*child,id)) return found;
    return nullptr;
  }
'''
    text = replace_once(text,
        '  void propagate(const EventContext&,const Acts::GeometryContext&,const Acts::BoundTrackParameters&,const std::array<double,3>&,const Acts::Surface* exactTarget=nullptr);',
        finder+'  void propagate(const EventContext&,const Acts::GeometryContext&,const Acts::BoundTrackParameters&,const std::array<double,3>&,const Acts::Surface* exactTarget=nullptr);')
    text = replace_once(text, 'tracking->findSurface(Acts::GeometryIdentifier(row.geometry))',
                       'completeSurface(*tracking->highestTrackingVolume(),row.geometry)')
    text = replace_once(text, 'tracking->findSurface(Acts::GeometryIdentifier(row.targetGeometry))',
                       'completeSurface(*tracking->highestTrackingVolume(),row.targetGeometry)')
    (source/UNIT).write_text(text)
    relative = str(UNIT)
    patch = 'diff --git a/'+relative+' b/'+relative+'\n'+''.join(difflib.unified_diff(
        previous.splitlines(True), text.splitlines(True), fromfile='a/'+relative, tofile='b/'+relative))
    increment = repo/'patches/0041_native_passive_surface_transport_lookup.patch'
    increment.write_text(patch)
    subprocess.run(['git', '-C', str(source), 'apply', '--reverse', '--check', str(increment)], check=True)
    shutil.copy2(source/UNIT, output/'CalypsoSlopeDiagnosticsAlg_v2.cxx')
    manifest['passive_lookup_increment'] = str(increment)
    manifest['probe_v2_source'] = str(output/'CalypsoSlopeDiagnosticsAlg_v2.cxx')
    dump(output/'manifest.json', manifest)
    dump(output/'build_tasks_v2.json', [dict(kind='analysis', command=['python3', str(output/'code/scripts/build_slope_diagnostics.py'),
        '--source', str(source), '--output', str(output/'software_probe_v2')], output=str(output/'jobs/build_probe_v2'))])
    for variant in ('beforeGeo', 'afterGeo'):
        for suffix, directory in [('_tasks', 'runs_v2'), ('_replay_tasks', 'legacy_replay_v2')]:
            old_tasks = json.loads((output/(variant+suffix+'.json')).read_text())
            new_tasks = []
            for old_task in old_tasks:
                old_source = Path(old_task['command'][-1])
                task = json.loads(old_source.read_text())
                old_output = Path(task['output'])
                relative_output = old_output.relative_to(output/('runs' if suffix == '_tasks' else 'legacy_replay'))
                task['output'] = str(output/directory/relative_output)
                new_source = old_source.with_stem(old_source.stem+'_v2')
                dump(new_source, task)
                new_tasks.append(dict(kind='analysis', command=old_task['command'][:-1]+[str(new_source)], output=task['output']))
            dump(output/(variant+suffix+'_v2.json'), new_tasks)


def overlay_v2(root):
    output = root/'common_prior_ckf_chart_v1/bound_transport_v1'
    manifest = json.loads((output/'manifest.json').read_text())
    built = output/manifest.get('latest_probe_software', 'software_probe_v2')
    assert (built/'software.json').is_file()
    for variant in manifest['reference_libraries']:
        previous = output/('overlay_'+variant)
        target = output/('overlay_'+variant+'_v2')
        (target/'lib').mkdir(parents=True, exist_ok=False)
        (target/'python').mkdir()
        for path in (previous/'lib').iterdir():
            origin = built/'lib'/path.name if 'CalypsoSlopeDiagnostics' in path.name else path.resolve()
            (target/'lib'/path.name).symlink_to(origin)
        for path in (previous/'python').iterdir():
            origin = built/'python'/path.name if path.name == 'CalypsoSlopeDiagnostics' else path.resolve()
            (target/'python'/path.name).symlink_to(origin)
    assert (built/'CalypsoSlopeDiagnosticsAlg.cxx').read_bytes() == Path(manifest['probe_v2_source']).read_bytes()
    dump(output/'source_validation_v2.json', dict(observer_source_exact=True,
        compiled_component=json.loads((built/'software.json').read_text()),
        before_after_geometry_and_field_libraries_reused=manifest['reference_libraries'],
        passive_surfaces_found_by_traversing_actual_native_volumes_layers_approaches=True,
        native_extrapolation_tool_unchanged=True, legacy_replay_required=True))


def surface_metadata(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/bound_transport_v1'
    manifest = json.loads((output/'manifest.json').read_text())
    source = Path(manifest['source'])
    previous = (source/UNIT).read_text()
    text = replace_once(previous, '  int m_valid[3]{};',
        '  int m_valid[3]{},m_inputSurfaceType=0,m_outputSurfaceType[3]{};')
    text = replace_once(text,
        '  m_tree->Branch("input_surface_native",m_inputSurface,"input_surface_native[4][4]/D");',
        '  m_tree->Branch("input_surface_native",m_inputSurface,"input_surface_native[4][4]/D");\n'
        '  m_tree->Branch("input_surface_type",&m_inputSurfaceType,"input_surface_type/I");\n'
        '  m_tree->Branch("projection_surface_type",m_outputSurfaceType,"projection_surface_type[3]/I");')
    text = replace_once(text, '  const auto& inputTransform=start.referenceSurface().transform(gctx).matrix();',
        '  m_inputSurfaceType=static_cast<int>(start.referenceSurface().type());\n'
        '  const auto& inputTransform=start.referenceSurface().transform(gctx).matrix();')
    text = replace_once(text, '    m_valid[i]=0;', '    m_valid[i]=0;m_outputSurfaceType[i]=-1;')
    text = replace_once(text, '    const auto& outputTransform=result->referenceSurface().transform(gctx).matrix();',
        '    m_outputSurfaceType[i]=static_cast<int>(result->referenceSurface().type());\n'
        '    const auto& outputTransform=result->referenceSurface().transform(gctx).matrix();')
    (source/UNIT).write_text(text)
    relative = str(UNIT)
    patch = 'diff --git a/'+relative+' b/'+relative+'\n'+''.join(difflib.unified_diff(
        previous.splitlines(True), text.splitlines(True), fromfile='a/'+relative, tofile='b/'+relative))
    increment = repo/'patches/0042_native_transport_surface_types.patch'
    increment.write_text(patch)
    subprocess.run(['git', '-C', str(source), 'apply', '--reverse', '--check', str(increment)], check=True)
    shutil.copy2(source/UNIT, output/'CalypsoSlopeDiagnosticsAlg_v3.cxx')
    manifest['surface_metadata_increment'] = str(increment)
    manifest['probe_v2_source'] = str(output/'CalypsoSlopeDiagnosticsAlg_v3.cxx')
    manifest['latest_probe_software'] = 'software_probe_v3'
    dump(output/'manifest.json', manifest)
    dump(output/'build_tasks_v3.json', [dict(kind='analysis', command=['python3', str(output/'code/scripts/build_slope_diagnostics.py'),
        '--source', str(source), '--output', str(output/'software_probe_v3')], output=str(output/'jobs/build_probe_v3'))])


def surface_inputs(root):
    """Preserve native local coordinates exactly for fully rotated surfaces."""
    output = root/'common_prior_ckf_chart_v1/bound_transport_v1'
    manifest = json.loads((output/'manifest.json').read_text())
    paths = json.loads(Path(manifest['source_paths']).read_text())
    paths = {(p['dsid'], p['entry'], p['component']): p for p in paths}
    cases = json.loads((output/'cases.json').read_text())
    target = output/'inputs_v2'
    target.mkdir(exist_ok=False)
    code = output/'run_code_v2'
    for folder, name in [('scripts', 'run_slope_transport.py'), ('configs', 'charge_prior_transport.py')]:
        (code/folder).mkdir(parents=True, exist_ok=True)
        shutil.copy2(Path(__file__).resolve().parents[1]/folder/name, code/folder/name)
    for dsid in ('100044', '100116', '100117'):
        sample = [c for c in cases if c['dsid'] == dsid]
        for variant in ('beforeGeo', 'afterGeo'):
            original = (output/'inputs'/(dsid+'_'+variant+'.txt')).read_text().splitlines()
            assert len(original) == len(sample)
            rows = []
            for case, line in zip(sample, original):
                fields = line.split()
                assert int(fields[2]) == case['candidate']
                if variant == 'afterGeo' and case['mode'] == 'fixed_before_input':
                    path = paths[dsid, case['entry'], case['component']]
                    index = case['from_forward_index']
                    a, b = path['before_states'][index], path['after_states'][index]
                    mean, covariance = rotate_tangent(np.asarray(a['filtered']), np.asarray(a['filtered_cov']),
                        np.asarray(a['surface']), np.asarray(b['surface']))
                    fields[9:51] = [format(v, '.17g') for v in np.r_[mean, covariance.ravel()]]
                rows.append(' '.join(fields)+'\n')
            (target/(dsid+'_'+variant+'.txt')).write_text(''.join(rows))
    for variant in ('beforeGeo', 'afterGeo'):
        for job in json.loads((output/(variant+'_tasks_v2.json')).read_text()):
            path = Path(job['command'][-1])
            task = json.loads(path.read_text())
            task['states'] = str(target/Path(task['states']).name)
            task['runner'] = str(code/'configs/charge_prior_transport.py')
            dump(path, task)
    manifest['fixed_input_mapping'] = dict(inputs=str(target), native_local_coordinates_preserved_exactly=True,
        analytic_global_angular_mean_and_complete_covariance_mapping=True,
        paired_full_surface_bases_checked=True, previous_inputs_preserved=str(output/'inputs'))
    dump(output/'manifest.json', manifest)


def analysis_tasks(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1/bound_transport_v1'
    code = output/'analysis_code_v4'
    files = ['scripts/validate_charge_components.py']+['analysis/'+name for name in
        ['charge_bound_transport.py', 'charge_bound_transport_figures.py', 'bound_surface_coordinates.py',
         'charge_forward_covariance.py', 'charge_component_refit.py', 'charge_refit_input.py',
         'charge_smoothing_operations.py', 'slope_tables.py', 'validate_charge_gun_figures.py']]
    code.mkdir(exist_ok=False)
    for relative in files:
        target = code/relative
        target.parent.mkdir(exist_ok=True, parents=True)
        shutil.copy2(repo/relative, target)
    dump(output/'analysis_tasks_v2.json', [dict(kind='analysis', command=['python3',
        str(code/'analysis/charge_bound_transport.py'), '--input', str(output),
        '--output', str(output/'analysis_v2'), '--phase', 'analyze', '--revision', '2'],
        output=str(output/'jobs/analysis_v2'))])
    manifest = json.loads((output/'manifest.json').read_text())
    manifest['analysis_code'] = str(code)
    manifest['analysis_files'] = files
    dump(output/'manifest.json', manifest)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['source', 'tasks', 'overlay', 'passive-surfaces', 'overlay-v2', 'surface-metadata', 'surface-inputs', 'analysis-tasks'], default='source')
    args = parser.parse_args()
    {'source': source, 'tasks': tasks, 'overlay': overlay, 'passive-surfaces': passive_surfaces,
     'overlay-v2': overlay_v2, 'surface-metadata': surface_metadata,
     'surface-inputs': surface_inputs, 'analysis-tasks': analysis_tasks}[args.phase](args.root.resolve())


if __name__ == '__main__':
    main()
