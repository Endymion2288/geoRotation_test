#!/usr/bin/env python3
"""Observe actual Calypso backward smoothing, with default/observer replay gates."""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_charge_component_refit import PACKAGE, provenance, replace_once
from prepare_cluster_measurements import dump
from validate_charge_components import input_rows


NAME = 'smoothing_operations_v1'
SERIES = 'matched_n32_f256_v1'
VARIANTS = ('beforeGeo', 'afterGeo')

PROPERTIES = '''  Gaudi::Property<std::string> m_smoothingOperationFile {this,"SmoothingOperationFile","",
      "Read-only actual ACTS backward-operation audit on private component candidates"};
  Gaudi::Property<std::string> m_smoothingOperationSelectionFile {this,"SmoothingOperationSelectionFile","",
      "Local event/component pairs to observe, without restricting reconstruction"};
  std::ofstream m_smoothingOperations;
  std::set<std::pair<size_t,int>> m_smoothingOperationSelections;
'''

INITIALIZE = '''  if (!m_smoothingOperationFile.value().empty()) {
    if (m_seedComponentInputFile.value().empty() || !m_isolateSmoothing ||
        m_smoothingOperationSelectionFile.value().empty()) return StatusCode::FAILURE;
    std::ifstream input(m_smoothingOperationSelectionFile.value());
    if (!input) return StatusCode::FAILURE;
    size_t event; int component;
    while (input >> event >> component) {
      if (component<0 || !m_smoothingOperationSelections.emplace(event,component).second)
        return StatusCode::FAILURE;
    }
    if (!input.eof() || m_smoothingOperationSelections.empty()) return StatusCode::FAILURE;
    m_smoothingOperations.open(m_smoothingOperationFile.value());
    if (!m_smoothingOperations) return StatusCode::FAILURE;
  }
'''


def output_path(root):
    return root/'common_prior_ckf_chart_v1'/NAME


def source(root):
    repo = Path(__file__).resolve().parents[1]
    output = output_path(root)
    output.mkdir(exist_ok=False)
    code = output/'build_code'
    code.mkdir()
    shutil.copy2(repo/'scripts/build_native_calypso.py', code/'build_native_calypso.py')
    variants, jobs = {}, []
    for variant in VARIANTS:
        before = variant == 'beforeGeo'
        old = repo.parent/('geoRotation_calypso_ckf_chart_beforeGeo' if before else
                           'geoRotation_calypso_component_refit_afterGeo')
        target = repo.parent/('geoRotation_calypso_smoothing_audit_'+variant)
        commit = subprocess.check_output(['git', '-C', str(old), 'rev-parse', 'HEAD'], text=True).strip()
        subprocess.run(['git', '-C', str(old), 'worktree', 'add', '--detach', str(target), commit], check=True)
        for name in ('FaserActsGeometry', 'FaserActsGeometryInterfaces', 'FaserActsKalmanFilter'):
            relative = Path('Tracking/Acts')/name
            shutil.copytree(old/relative, target/relative, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        changes = {'CandidateSmoothingAudit.h': (repo/'patches/CandidateSmoothingAudit.h').read_text()}
        if not (old/PACKAGE/'RefitDirectionChart.h').exists():
            changes['RefitDirectionChart.h'] = (repo.parent/'geoRotation_calypso_ckf_chart_beforeGeo'/
                                               PACKAGE/'RefitDirectionChart.h').read_text()
        header = replace_once((old/PACKAGE/'CKF2.h').read_text(), '#include <set>',
                              '#include <set>\n#include <fstream>')
        changes['CKF2.h'] = replace_once(header, '  Gaudi::Property<bool> m_componentRefit',
                                        PROPERTIES+'  Gaudi::Property<bool> m_componentRefit')
        ckf = replace_once((old/PACKAGE/'CKF2.cxx').read_text(),
            '#include "Acts/TrackFitting/GainMatrixSmoother.hpp"',
            '#include "Acts/TrackFitting/GainMatrixSmoother.hpp"\n#include "CandidateSmoothingAudit.h"')
        ckf = replace_once(ckf, '  if (!m_fullSeedInputFile.value().empty()) {',
                          INITIALIZE+'  if (!m_fullSeedInputFile.value().empty()) {')
        original = ('smoothCkfCandidate(m_regularCkfSmoothingChart,' if before else
                    'Acts::GainMatrixSmoother()(')
        call = original+'\n            gctx,selectedTracks.trackStateContainer(),dest.tipIndex(),*m_logger)'
        audited = '''(m_smoothingOperations.is_open() &&
            m_smoothingOperationSelections.count({m_numberOfEvents-1,component.component}))
            ? CandidateSmoothingAudit::run(REGULAR,gctx,
                selectedTracks.trackStateContainer(),dest.tipIndex(),*m_logger,m_smoothingOperations,
                m_numberOfEvents-1,component.seed,component.component,selected.index)
            : CALL'''.replace('REGULAR', 'm_regularCkfSmoothingChart' if before else 'false').replace('CALL', call)
        ckf = replace_once(ckf, call, audited)
        changes['CKF2.cxx'] = ckf
        increment = ''
        for name, content in changes.items():
            relative = str(PACKAGE/name)
            (target/PACKAGE/name).write_text(content)
            previous = old/PACKAGE/name
            increment += 'diff --git a/'+relative+' b/'+relative+'\n'
            if not previous.exists():
                increment += 'new file mode 100644\n'
            increment += ''.join(difflib.unified_diff(
                previous.read_text().splitlines(True) if previous.exists() else [],
                content.splitlines(True), fromfile='a/'+relative if previous.exists() else '/dev/null',
                tofile='b/'+relative))
        patch = repo/'patches'/('0039_smoothing_operation_audit_'+variant+'.patch')
        patch.write_text(increment)
        subprocess.run(['git', '-C', str(old), 'apply', '--check', str(patch)], check=True)
        subprocess.run(['git', '-C', str(target), 'apply', '--reverse', '--check', str(patch)], check=True)
        tracked = output/(variant+'_source.patch')
        tracked.write_text(subprocess.check_output(['git', '-C', str(target), 'diff', 'HEAD'], text=True))
        previous = (root/'common_prior_ckf_chart_v1/software_beforeGeo' if before else
                    root/'common_prior_component_refit_v1/software_afterGeo')
        variants[variant] = dict(source=str(target), source_base=str(old), commit=commit,
            tracked_patch=str(tracked), increment=str(patch), geometry_reused_from=str(previous))
        jobs.append(dict(kind='analysis', command=['python3', str(code/'build_native_calypso.py'),
            '--source', str(target), '--output', str(output/('software_'+variant)), '--jobs', '4',
            '--source-patch', str(tracked), '--reuse-geometry', str(previous)],
            output=str(output/'jobs'/('build_'+variant))))
    dump(output/'build_tasks.json', jobs)
    dump(output/'manifest.json', dict(variants=variants, isolated_smoothing=True,
        scope='Read-only actual ACTS calculate() observation; default and enabled replay required.'))
    provenance(root, output)
    print(output)


def tasks(root):
    repo = Path(__file__).resolve().parents[1]
    output = output_path(root)
    reference = root/'common_prior_ckf_chart_v1/analytic_covariance_v1'
    code = output/'run_code'
    (output/'tasks').mkdir(exist_ok=True)
    for relative in ('configs/slope_tracking.py', 'scripts/run_slope_task.py'):
        target = code/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, target)
    manifest = json.loads((output/'manifest.json').read_text())
    manifest['points'] = []
    jobs = {(v, mode): [] for v in VARIANTS for mode in ('off', 'on')}
    for point in json.loads((reference/'manifest.json').read_text())['points']:
        if point['series'] != SERIES:
            continue
        dsid = point['dsid']
        for variant in VARIANTS:
            original_folder = Path(point[variant])
            original = json.loads((original_folder/'source_task.json').read_text())
            rows = input_rows(Path(original['component_file']))
            selection = output/'inputs'/(dsid+'_'+variant+'_selection.txt')
            selection.parent.mkdir(parents=True, exist_ok=True)
            selection.write_text(''.join(f'{entry} {component}\n' for entry,component in sorted(rows)))
            result = dict(dsid=dsid, variant=variant, reference=str(original_folder),
                          selection=str(selection), input_components=len(rows))
            for mode in ('off', 'on'):
                folder = output/'runs'/dsid/variant/mode
                options = original['options'].copy()
                if mode == 'on':
                    options += ['--smoothing-operation-selection-file', str(selection)]
                task = dict(original, output=str(folder), runner=str(code/'configs/slope_tracking.py'),
                            options=options, smoothing_operation_observer=mode)
                path = output/'tasks'/(dsid+'_'+variant+'_'+mode+'.json')
                dump(path, task)
                jobs[variant,mode].append(dict(kind='analysis', command=['python3',
                    str(code/'scripts/run_slope_task.py'), '--task', str(path)], output=str(folder)))
                result[mode] = str(folder)
            manifest['points'].append(result)
    for (variant,mode), items in jobs.items():
        dump(output/(variant+'_'+mode+'_tasks.json'), items)
    dump(output/'manifest.json', manifest)
    print(json.dumps(dict(points=len(manifest['points']), tasks=sum(map(len,jobs.values())))))


def overlay(root, variants=VARIANTS):
    output = output_path(root)
    for variant in variants:
        previous = (root/'common_prior_ckf_chart_v1/overlay_beforeGeo' if variant == 'beforeGeo' else
                    root/'common_prior_component_refit_v1/overlay_afterGeo')
        target = output/('overlay_'+variant)
        shutil.copytree(previous, target)
        software = output/('software_'+variant)
        for folder in ('lib', 'python'):
            shutil.copytree(software/folder, target/folder, dirs_exist_ok=True)
        dump(target/'smoothing_audit_overlay.json', dict(base=str(previous), replacement=str(software),
            source_validation=str(output/'source_validation.json')))
    print(output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['source', 'tasks', 'overlay'], required=True)
    parser.add_argument('--variant', choices=VARIANTS, nargs='+', help='Overlay only these completed builds')
    args = parser.parse_args()
    if args.phase == 'overlay':
        overlay(args.root.resolve(), args.variant or VARIANTS)
    else:
        {'source': source, 'tasks': tasks}[args.phase](args.root.resolve())


if __name__ == '__main__':
    main()
