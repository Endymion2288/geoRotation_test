#!/usr/bin/env python3
"""Reuse the actual refit chart adapters in native Calypso CKF controls."""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_charge_component_refit import provenance, replace_once
from prepare_cluster_measurements import dump


PACKAGE = Path('Tracking/Acts/FaserActsKalmanFilter/src')


def prepare_source(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1'
    output.mkdir(exist_ok=False)
    old = repo.parent/'geoRotation_calypso_component_refit_input_beforeGeo'
    source = repo.parent/'geoRotation_calypso_ckf_chart_beforeGeo'
    commit = subprocess.check_output(['git', '-C', str(old), 'rev-parse', 'HEAD'], text=True).strip()
    subprocess.run(['git', '-C', str(old), 'worktree', 'add', '--detach', str(source), commit], check=True)
    for name in ('FaserActsGeometry', 'FaserActsGeometryInterfaces', 'FaserActsKalmanFilter'):
        relative = Path('Tracking/Acts')/name
        shutil.copytree(old/relative, source/relative, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    changes = {}
    fitter = (old/PACKAGE/'TrackFitterFunction.cxx').read_text()
    start = fitter.index('// Direction-chart control inside the real Calypso fitter.')
    end = fitter.index('// Keep the native ACTS update and its chi2.')
    adapters = fitter[start:end]
    header = '''#ifndef FASER_REGULAR_DIRECTION_KALMAN_H
#define FASER_REGULAR_DIRECTION_KALMAN_H
#include "RefitDirectionChart.h"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
'''+adapters+'''
inline Acts::Result<void> smoothCkfCandidate(bool regular,
    const Acts::GeometryContext& gctx, Acts::VectorMultiTrajectory& states,
    std::size_t tip, const Acts::Logger& logger) {
  if (regular) return RegularDirectionGainMatrixSmoother()(gctx,states,tip,logger);
  return Acts::GainMatrixSmoother()(gctx,states,tip,logger);
}
#endif
'''
    changes['RegularDirectionKalman.h'] = header
    changes['TrackFitterFunction.cxx'] = fitter[:start]+'#include "RegularDirectionKalman.h"\n\n'+fitter[end:]
    finder = (old/PACKAGE/'TrackFinderFunction.cxx').read_text()
    finder = replace_once(finder, '#include "TrackFinderFunction.h"',
                          '#include "TrackFinderFunction.h"\n#include "RegularDirectionKalman.h"')
    finder = replace_once(finder, '  Acts::GainMatrixSmoother kfSmoother;',
        '  Acts::GainMatrixSmoother kfSmoother;\n'
        '  RegularDirectionGainMatrixUpdater regularUpdater;\n'
        '  RegularDirectionGainMatrixSmoother regularSmoother;\n'
        '  bool regularUpdate=false, regularSmoothing=false;')
    update = '''    extensions.updater.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
        &kfUpdater);'''
    finder = replace_once(finder, update,
        '    if (regularUpdate) {\n'
        '      extensions.updater.connect<&RegularDirectionGainMatrixUpdater::operator()>(&regularUpdater);\n'
        '    } else {\n'+update+'\n    }')
    smoothing = '''    extensions.smoother.connect<
        &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
        &kfSmoother);'''
    finder = replace_once(finder, smoothing,
        '    if (regularSmoothing) {\n'
        '      extensions.smoother.connect<&RegularDirectionGainMatrixSmoother::operator()>(&regularSmoother);\n'
        '    } else {\n'+smoothing+'\n    }')
    finder = replace_once(finder, '    const Acts::Logger& logger) {',
        '    const Acts::Logger& logger, bool regularUpdate, bool regularSmoothing) {')
    finder = replace_once(finder, '  return std::make_shared<TrackFinderFunctionImpl>(std::move(trackFinder));',
        '  auto result=std::make_shared<TrackFinderFunctionImpl>(std::move(trackFinder));\n'
        '  result->regularUpdate=regularUpdate;\n'
        '  result->regularSmoothing=regularSmoothing;\n'
        '  return result;')
    changes['TrackFinderFunction.cxx'] = finder
    changes['TrackFinderFunction.h'] = replace_once((old/PACKAGE/'TrackFinderFunction.h').read_text(),
        '  const Acts::Logger& logger);',
        '  const Acts::Logger& logger, bool regularUpdate=false, bool regularSmoothing=false);')
    changes['CKF2.h'] = replace_once((old/PACKAGE/'CKF2.h').read_text(),
        '  Gaudi::Property<bool> m_componentRefit',
        '  Gaudi::Property<bool> m_regularCkfUpdateChart {this,"RegularCkfUpdateChart",false,\n'
        '      "Diagnostic: native CKF update in the existing beam-regular direction chart"};\n'
        '  Gaudi::Property<bool> m_regularCkfSmoothingChart {this,"RegularCkfSmoothingChart",false,\n'
        '      "Diagnostic: native CKF and private smoothing in the existing beam-regular direction chart"};\n'
        '  Gaudi::Property<bool> m_componentRefit')
    ckf = (old/PACKAGE/'CKF2.cxx').read_text()
    ckf = replace_once(ckf, '#include "Acts/TrackFitting/GainMatrixSmoother.hpp"',
        '#include "Acts/TrackFitting/GainMatrixSmoother.hpp"\n#include "RegularDirectionKalman.h"')
    ckf = replace_once(ckf,
        'm_resolvePassive, m_resolveMaterial, m_resolveSensitive, *m_logger);',
        'm_resolvePassive, m_resolveMaterial, m_resolveSensitive, *m_logger,\n'
        '                                  m_regularCkfUpdateChart, m_regularCkfSmoothingChart);')
    assert ckf.count('Acts::GainMatrixSmoother()(') == 2
    changes['CKF2.cxx'] = ckf.replace('Acts::GainMatrixSmoother()(',
        'smoothCkfCandidate(m_regularCkfSmoothingChart,')
    increment = ''
    for name, content in changes.items():
        path = source/PACKAGE/name
        path.write_text(content)
        relative = str(PACKAGE/name)
        previous = old/PACKAGE/name
        increment += 'diff --git a/'+relative+' b/'+relative+'\n'
        if not previous.exists():
            increment += 'new file mode 100644\n'
        increment += ''.join(difflib.unified_diff(previous.read_text().splitlines(True) if previous.exists() else [],
            content.splitlines(True), fromfile='a/'+relative if previous.exists() else '/dev/null', tofile='b/'+relative))
    patch = repo/'patches/0038_native_ckf_direction_chart.patch'
    patch.write_text(increment)
    tracked = output/'beforeGeo_source.patch'
    tracked.write_text(subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True))
    code = output/'build_code';code.mkdir()
    shutil.copy2(repo/'scripts/build_native_calypso.py', code/'build_native_calypso.py')
    previous = root/'common_prior_refit_input_v1/software_beforeGeo'
    dump(output/'manifest.json', dict(variants=dict(beforeGeo=dict(source=str(source),
        source_base=str(old), commit=commit, tracked_patch=str(tracked), increment=str(patch),
        geometry_reused_from=str(previous))), isolated_smoothing=True,
        scope='Existing actual Calypso refit chart adapters reused in native CKF update/smoothing; '
              'all switches default off, geometry/transport/measurements unchanged.'))
    dump(output/'build_tasks.json', [dict(kind='analysis', command=['python3',
        str(code/'build_native_calypso.py'), '--source', str(source), '--output', str(output/'software_beforeGeo'),
        '--jobs', '4', '--source-patch', str(tracked), '--reuse-geometry', str(previous)],
        output=str(output/'jobs/build_beforeGeo'))])
    provenance(root, output)
    print(json.dumps(dict(output=str(output), source=str(source)), indent=2))


def prepare_tasks(root):
    repo = Path(__file__).resolve().parents[1]
    output = root/'common_prior_ckf_chart_v1'
    code = output/'run_code'
    for relative in ('configs/slope_tracking.py', 'scripts/run_slope_task.py'):
        target = code/relative;target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, target)
    folder = output/'tasks';folder.mkdir()
    tasks = []
    for dsid in ('100044', '100116', '100117'):
        source = root/'common_prior_refit_input_v1/scale1_v1/tasks'/('matched_n32_f256_v1_'+dsid+'_beforeGeo_regular_both.json')
        original = json.loads(source.read_text())
        for mode, options in [('native', []), ('update', ['--regular-ckf-update-chart']),
                              ('smoothing', ['--regular-ckf-smoothing-chart']),
                              ('both', ['--regular-ckf-update-chart', '--regular-ckf-smoothing-chart'])]:
            task = dict(original, runner=str(code/'configs/slope_tracking.py'),
                options=[*original['options'], *options], ckf_chart=mode,
                output=str(output/'matched_n32_f256_v1'/dsid/mode), reference_task=str(source))
            path = folder/(dsid+'_'+mode+'.json');dump(path, task)
            tasks.append(dict(kind='analysis', command=['python3',str(code/'scripts/run_slope_task.py'),
                '--task',str(path)], output=task['output']))
    dump(output/'reconstruction_tasks.json', tasks)
    dump(output/'gate_tasks.json', [t for t in tasks if Path(t['output']).name == 'native'])
    dump(output/'control_tasks.json', [t for t in tasks if Path(t['output']).name != 'native'])
    print(json.dumps(dict(tasks=len(tasks), output=str(output)), indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['source', 'tasks'], required=True)
    args = parser.parse_args()
    if args.phase == 'source':
        prepare_source(args.root.resolve())
    else:
        prepare_tasks(args.root.resolve())


if __name__ == '__main__':
    main()
