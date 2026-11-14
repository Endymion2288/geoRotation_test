#!/usr/bin/env python3
"""Extend the real Calypso component diagnostic through conversion and refit."""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_cluster_measurements import dump


VARIANTS = ('beforeGeo', 'afterGeo')
SAMPLES = ('100044', '100116', '100117')
PACKAGE = Path('Tracking/Acts/FaserActsKalmanFilter/src')

PROPERTY = '''  Gaudi::Property<bool> m_componentRefit {this, "ComponentRefit", false,
      "Diagnostic: run each finite-prior component through native conversion and refit"};
'''

CONVERSION = '''
        if (m_componentRefit) {
          std::optional<Acts::BoundTrackParameters> fittedParams=std::nullopt;
          if (dest.hasReferenceSurface() && m_addFittedParamsToTrack) {
            fittedParams=Acts::BoundTrackParameters(
                dest.referenceSurface().getSharedPtr(),dest.parameters(),dest.covariance(),
                Acts::ParticleHypothesis::muon());
          }
          auto converted=m_createTrkTrackTool->createTrack(
              gctx,dest,fittedParams,m_backwardPropagation);
          m_audit->trk("converted",converted.get(),component.seed,selected.index);
          auto refitted=REFIT_CALL;
          m_audit->trk("final",refitted ? refitted.get() : converted.get(),
                       component.seed,selected.index,refitted ? 1 : 0);
          if (refitted) outputTracks->push_back(std::move(refitted));
          else outputTracks->push_back(std::move(converted));
        }
'''


def replace_once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new, 1)


def source(root, output):
    output.mkdir(exist_ok=False)
    repo = Path(__file__).resolve().parents[1]
    code = output/'build_code'
    code.mkdir()
    shutil.copy2(repo/'scripts/build_native_calypso.py', code/'build_native_calypso.py')
    base = Path('/eos/home-x/xcheng/FASER/geoRotation_calypso')
    tasks, provenance = [], {}
    for variant in VARIANTS:
        old = base.parent/('geoRotation_calypso_components_'+variant)
        target = base.parent/('geoRotation_calypso_component_refit_'+variant)
        commit = subprocess.check_output(['git', '-C', str(old), 'rev-parse', 'HEAD'], text=True).strip()
        subprocess.run(['git', '-C', str(base), 'worktree', 'add', '--detach', str(target), commit], check=True)
        for name in ('FaserActsGeometry', 'FaserActsKalmanFilter', 'FaserActsGeometryInterfaces'):
            relative = Path('Tracking/Acts')/name
            shutil.copytree(old/relative, target/relative, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        changes = {}
        header = (old/PACKAGE/'CKF2.h').read_text()
        changes['CKF2.h'] = replace_once(header,
            '  std::map<size_t,std::vector<SeedComponentInput>> m_seedComponentInputs;',
            '  std::map<size_t,std::vector<SeedComponentInput>> m_seedComponentInputs;\n'+PROPERTY)
        text = (old/PACKAGE/'CKF2.cxx').read_text()
        call = ('m_kalmanFitterTool1->fit(ctx,gctx,converted.get(),Acts::BoundVector::Zero(),m_isMC)'
                if variant == 'beforeGeo' else
                'm_kalmanFitterTool1->fitAudited(ctx,gctx,converted.get(),Acts::BoundVector::Zero(),'
                'm_isMC,m_audit.get(),component.seed,selected.index)')
        text = replace_once(text,
            '        m_audit->track("selected",gctx,dest,component.seed,selected.index);',
            '        m_audit->track("selected",gctx,dest,component.seed,selected.index);'+
            CONVERSION.replace('REFIT_CALL', call))
        text = replace_once(text, '  // Each component invokes this version\'s actual finder, calibrator, transport,',
            '  // Optional ComponentRefit continues the same candidate through native conversion/refit.\n'
            '  // Each component invokes this version\'s actual finder, calibrator, transport,')
        text = replace_once(text,
            '  // compete in the event-level shared-hit selector. Only CKF posterior states\n'
            '  // are observed in this opt-in mode; conversion/refit are not invoked.',
            '  // compete in the event-level shared-hit selector. The default component\n'
            '  // mode observes CKF states; ComponentRefit adds native conversion/refit.')
        changes['CKF2.cxx'] = text
        patch = ''
        for name, text in changes.items():
            (target/PACKAGE/name).write_text(text)
            filename = str(PACKAGE/name)
            patch += 'diff --git a/'+filename+' b/'+filename+'\n'
            patch += ''.join(difflib.unified_diff((old/PACKAGE/name).read_text().splitlines(True),
                text.splitlines(True), fromfile='a/'+filename, tofile='b/'+filename))
        increment = repo/'patches'/('0036_component_native_refit_'+variant+'.patch')
        increment.write_text(patch)
        tracked = subprocess.check_output(['git', '-C', str(target), 'diff', 'HEAD'], text=True)
        tracked_patch = output/(variant+'_source.patch')
        tracked_patch.write_text(tracked)
        software = root/'common_prior_components_v4'/('software_'+variant)
        tasks.append(dict(kind='analysis', command=['python3', str(code/'build_native_calypso.py'),
            '--source', str(target), '--output', str(output/('software_'+variant)), '--jobs', '4',
            '--source-patch', str(tracked_patch), '--reuse-geometry', str(software)],
            output=str(output/'jobs'/('build_'+variant))))
        provenance[variant] = dict(source=str(target), source_base=str(old), commit=commit,
            increment=str(increment), tracked_patch=str(tracked_patch), geometry_reused_from=str(software))
    dump(output/'build_tasks.json', tasks)
    dump(output/'manifest.json', dict(variants=provenance, isolated_smoothing=True,
        scope='Actual Calypso CKF and per-component native conversion/refit; independent component calls. '
              'No likelihood aggregation or production mixture fit is introduced.'))


def tasks(root, output, series_list, task_prefix=''):
    repo = Path(__file__).resolve().parents[1]
    code = output/'run_code'
    (code/'configs').mkdir(parents=True, exist_ok=True)
    (code/'scripts').mkdir(exist_ok=True)
    (output/'tasks').mkdir(exist_ok=True)
    for relative in ('configs/slope_tracking.py', 'scripts/run_slope_task.py'):
        target = code/relative
        if target.exists():
            assert target.read_bytes() == (repo/relative).read_bytes()
        else:
            shutil.copy2(repo/relative, target)
    jobs = {variant: [] for variant in VARIANTS}
    for series in series_list:
        previous = root/'common_prior_components_v4'/series
        for variant in VARIANTS:
            modes = ('native', 'regular_both') if variant == 'beforeGeo' else ('native',)
            for dsid in SAMPLES:
                original = json.loads((previous/dsid/(variant+'_source.json')).read_text())
                for mode in modes:
                    options = [*original['options'], '--component-refit']
                    if mode == 'regular_both':
                        options += ['--regular-refit-update-chart', '--regular-refit-smoothing-chart']
                    task = dict(original, runner=str(code/'configs/slope_tracking.py'), options=options,
                        component_refit=True, refit_chart=mode, prior_series=series,
                        output=str(output/series/dsid/(variant+'_'+mode)))
                    path = output/'tasks'/(series+'_'+dsid+'_'+variant+'_'+mode+'.json')
                    dump(path, task)
                    jobs[variant].append(dict(kind='analysis', command=['python3',
                        str(code/'scripts/run_slope_task.py'), '--task', str(path)], output=task['output']))
    for variant in VARIANTS:
        dump(output/(task_prefix+variant+'_tasks.json'), jobs[variant])


def provenance(root, output):
    """Freeze added headers as well as tracked diffs used by the native build."""
    manifest = json.loads((output/'manifest.json').read_text())
    checks = {}
    for variant, info in manifest['variants'].items():
        source = Path(info['source'])
        tracked = subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True)
        assert tracked == Path(info['tracked_patch']).read_text()
        patch = tracked
        additions = subprocess.check_output(
            ['git', '-C', str(source), 'ls-files', '--others', '--exclude-standard'], text=True).splitlines()
        saved = output/('source_additions_'+variant)
        saved.mkdir(exist_ok=False)
        for name in additions:
            assert name.startswith('Tracking/Acts/') and name.endswith(('.h', '.cxx'))
            original = source/name
            relative = Path(name)
            target = saved/relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original, target)
            assert target.read_bytes() == original.read_bytes()
            patch += 'diff --git a/'+name+' b/'+name+'\nnew file mode 100644\n'
            patch += ''.join(difflib.unified_diff([], original.read_text().splitlines(True),
                fromfile='/dev/null', tofile='b/'+name))
        full = output/(variant+'_complete_source.patch')
        full.write_text(patch)
        # Applying in reverse must remove additions and restore every tracked
        # file to the actual branch commit, without changing the worktree.
        subprocess.run(['git', '-C', str(source), 'apply', '--reverse', '--check', str(full)], check=True)
        checks[variant] = dict(source=str(source), complete_patch=str(full),
                               added_source_files=additions, added_files_byte_equal=True,
                               tracked_diff_exact=True, complete_patch_reverse_check=True,
                               build_manifest=str(output/('software_'+variant)/'software.json'))
    dump(output/'source_validation.json', checks)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['source', 'tasks', 'provenance'], required=True)
    parser.add_argument('--output-name', default='common_prior_component_refit_v1')
    parser.add_argument('--series', nargs='+', default=['matched_n32_f256_v1', 'matched_n64_f1024_v1'])
    parser.add_argument('--task-prefix', default='')
    args = parser.parse_args()
    root = args.root.resolve()
    if args.phase == 'source':
        source(root, root/args.output_name)
    elif args.phase == 'tasks':
        tasks(root, root/args.output_name, args.series, args.task_prefix)
    else:
        provenance(root, root/args.output_name)


if __name__ == '__main__':
    main()
