#!/usr/bin/env python3
"""Add a read-only native beforeGeo refit-input audit for the common prior."""
import argparse
import difflib
import json
from pathlib import Path
import re
import shutil
import subprocess

from prepare_charge_component_refit import provenance, replace_once
from prepare_cluster_measurements import dump


def prepare_source(root):
    output = root/'common_prior_refit_input_v1'
    output.mkdir(exist_ok=False)
    repo = Path(__file__).resolve().parents[1]
    old = Path('/eos/home-x/xcheng/FASER/geoRotation_calypso_component_refit_beforeGeo')
    source = old.parent/'geoRotation_calypso_component_refit_input_beforeGeo'
    commit = subprocess.check_output(['git', '-C', str(old), 'rev-parse', 'HEAD'], text=True).strip()
    subprocess.run(['git', '-C', str(old), 'worktree', 'add', '--detach', str(source), commit], check=True)
    for package in ('FaserActsGeometry', 'FaserActsGeometryInterfaces', 'FaserActsKalmanFilter'):
        relative = Path('Tracking/Acts')/package
        shutil.copytree(old/relative, source/relative, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    package = Path('Tracking/Acts/FaserActsKalmanFilter/src')
    changes = {}
    text = (old/package/'KalmanFitterTool.h').read_text()
    kind = re.search(r'\b(struct|class) TrackingAudit\b', (old/package/'TrackingAudit.h').read_text())[1]
    text = replace_once(text, 'class FaserSCT_ID;', 'class FaserSCT_ID;\n'+kind+' TrackingAudit;')
    declaration = '''  std::unique_ptr<Trk::Track> fitAudited(
      const EventContext &ctx, const Acts::GeometryContext &gctx, Trk::Track *inputTrack,
      const Acts::BoundVector &inputVector, bool isMC, TrackingAudit *audit,
      int auditSeed, int auditCandidate) const;

'''
    text = replace_once(text, '  virtual Acts::MagneticFieldContext getMagneticFieldContext',
                        declaration+'  virtual Acts::MagneticFieldContext getMagneticFieldContext')
    changes['KalmanFitterTool.h'] = text
    text = (old/package/'KalmanFitterTool.cxx').read_text()
    text = replace_once(text, '#include "KalmanFitterTool.h"', '#include "KalmanFitterTool.h"\n#include "TrackingAudit.h"')
    signature = '''KalmanFitterTool::fit(const EventContext &ctx, const Acts::GeometryContext &gctx,
                      Trk::Track* inputTrack, const Acts::BoundVector& inputVector,
                      bool isMC) const {'''
    replacement = '''KalmanFitterTool::fit(const EventContext &ctx, const Acts::GeometryContext &gctx,
                      Trk::Track* inputTrack, const Acts::BoundVector& inputVector,
                      bool isMC) const {
  return fitAudited(ctx,gctx,inputTrack,inputVector,isMC,nullptr,-1,-1);
}

std::unique_ptr<Trk::Track>
KalmanFitterTool::fitAudited(const EventContext &ctx, const Acts::GeometryContext &gctx,
                      Trk::Track* inputTrack, const Acts::BoundVector& inputVector,
                      bool isMC, TrackingAudit* audit, int auditSeed, int auditCandidate) const {'''
    text = replace_once(text, signature, replacement)
    prefix, separator, body = text.partition('KalmanFitterTool::fitAudited(')
    assert separator
    body = replace_once(body,
        '    trackParameters.covariance() = (trackParameters.covariance().value())*10; \n  }',
        '    trackParameters.covariance() = (trackParameters.covariance().value())*10; \n  }\n'
        '  if(audit) audit->parameters("refit_input",gctx,trackParameters,auditSeed,auditCandidate);')
    body = replace_once(body, '    const auto& track = result.value();',
        '    const auto& track = result.value();\n'
        '    if(audit) audit->track("refit",gctx,track,auditSeed,auditCandidate);')
    changes['KalmanFitterTool.cxx'] = prefix+separator+body
    text = (old/package/'CKF2.cxx').read_text()
    text = replace_once(text,
        'm_kalmanFitterTool1->fit(ctx,gctx,converted.get(),Acts::BoundVector::Zero(),m_isMC)',
        'm_kalmanFitterTool1->fitAudited(ctx,gctx,converted.get(),Acts::BoundVector::Zero(),'
        'm_isMC,m_audit.get(),component.seed,selected.index)')
    changes['CKF2.cxx'] = text
    patch = ''
    for name, text in changes.items():
        (source/package/name).write_text(text)
        filename = str(package/name)
        patch += 'diff --git a/'+filename+' b/'+filename+'\n'
        patch += ''.join(difflib.unified_diff((old/package/name).read_text().splitlines(True),
            text.splitlines(True), fromfile='a/'+filename, tofile='b/'+filename))
    increment = repo/'patches/0037_component_refit_input_beforeGeo.patch'
    increment.write_text(patch)
    tracked = output/'beforeGeo_source.patch'
    tracked.write_text(subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True))
    code = output/'build_code'
    code.mkdir()
    shutil.copy2(repo/'scripts/build_native_calypso.py', code/'build_native_calypso.py')
    previous = root/'common_prior_component_refit_v1/software_beforeGeo'
    dump(output/'manifest.json', dict(variants=dict(beforeGeo=dict(source=str(source),
        source_base=str(old), commit=commit, increment=str(increment), tracked_patch=str(tracked),
        geometry_reused_from=str(previous))), isolated_smoothing=True,
        scope='Read-only native refit seed/state audit; no covariance, surface, measurement or fitting formula changes.'))
    dump(output/'build_tasks.json', [dict(kind='analysis', command=['python3',str(code/'build_native_calypso.py'),
        '--source', str(source), '--output', str(output/'software_beforeGeo'), '--jobs', '4',
        '--source-patch', str(tracked), '--reuse-geometry', str(previous)],
        output=str(output/'jobs/build_beforeGeo'))])
    provenance(root, output)
    print(json.dumps(dict(output=str(output), native_source=str(source), increment=str(increment)), indent=2))


def prepare_tasks(root, run_name, series_list, variants, covariance_scale):
    previous = root/'common_prior_component_refit_v1'
    output = root/'common_prior_refit_input_v1'
    run = output/run_name
    run.mkdir(exist_ok=False)
    code = run/'code'
    shutil.copytree(previous/'run_code', code)
    (run/'tasks').mkdir()
    jobs = {variant: [] for variant in variants}
    for series in series_list:
        for variant in variants:
            modes = ('native', 'regular_both') if variant == 'beforeGeo' else ('native',)
            for dsid in ('100044', '100116', '100117'):
                for mode in modes:
                    source = previous/'tasks'/(series+'_'+dsid+'_'+variant+'_'+mode+'.json')
                    task = json.loads(source.read_text())
                    options = task['options'].copy()
                    if covariance_scale is not None:
                        if '--refit-covariance-scale' in options:
                            index = options.index('--refit-covariance-scale')
                            del options[index:index+2]
                        options += ['--refit-covariance-scale', str(covariance_scale)]
                    task.update(runner=str(code/'configs/slope_tracking.py'), options=options,
                        output=str(run/series/dsid/(variant+'_'+mode)),
                        refit_input_audit=True, reference_component_task=str(source))
                    path = run/'tasks'/source.name
                    dump(path, task)
                    jobs[variant].append(dict(kind='analysis', command=['python3',
                        str(code/'scripts/run_slope_task.py'), '--task', str(path)], output=task['output']))
    for variant, group in jobs.items():
        dump(run/(variant+'_tasks.json'), group)
    dump(run/'manifest.json', dict(series=series_list, variants=variants,
        source_full_chain=str(previous), isolated_smoothing=True,
        refit_diagonal_covariance_scale=covariance_scale,
        scope='Read-only input audit replay; optional explicitly recorded diagonal covariance scale control.'))
    print(json.dumps(dict(output=str(run), tasks={v: len(t) for v, t in jobs.items()}), indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['source', 'tasks'], default='source')
    parser.add_argument('--run-name', default='replay_v1')
    parser.add_argument('--series', nargs='+', default=['matched_n32_f256_v1'])
    parser.add_argument('--variants', nargs='+', choices=['beforeGeo', 'afterGeo'], default=['beforeGeo'])
    parser.add_argument('--refit-covariance-scale', type=float)
    args = parser.parse_args()
    if args.phase == 'source':
        prepare_source(args.root.resolve())
    else:
        prepare_tasks(args.root.resolve(), args.run_name, args.series, args.variants,
                      args.refit_covariance_scale)


if __name__ == '__main__':
    main()
