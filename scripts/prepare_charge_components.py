#!/usr/bin/env python3
"""Prepare native Calypso CKF calls for finite Gaussian prior components."""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_cluster_measurements import dump


HEADER = '''
  struct SeedComponentInput {
    int seed, component;
    FullSeedInput input;
  };
  Gaudi::Property<std::string> m_seedComponentInputFile {this, "SeedComponentInputFile", "",
      "Diagnostic: independent native CKF calls, local event/native seed/component/mean/C/surface"};
  std::map<size_t,std::vector<SeedComponentInput>> m_seedComponentInputs;
'''

INITIALIZE = '''
  if (!m_seedComponentInputFile.value().empty()) {
    if (!m_fullSeedInputFile.value().empty() || m_runOnlyInputSeeds ||
        !m_audit || !m_auditStates || !m_selectedStatesOnly || !m_isolateSmoothing)
      return StatusCode::FAILURE;
    std::ifstream input(m_seedComponentInputFile.value());
    if (!input) return StatusCode::FAILURE;
    size_t entry; int seed, component;
    std::set<std::pair<size_t,int>> seen;
    while (input >> entry >> seed >> component) {
      FullSeedInput value;
      for (int i=0;i<6;++i) input >> value.mean[i];
      for (int i=0;i<6;++i) for (int j=0;j<6;++j) input >> value.covariance(i,j);
      for (int i=0;i<4;++i) for (int j=0;j<4;++j) input >> value.transform.matrix()(i,j);
      if (!input || seed<0 || component<0 || !value.mean.allFinite() ||
          !value.covariance.allFinite() || !value.transform.matrix().allFinite() ||
          !seen.emplace(entry,component).second) return StatusCode::FAILURE;
      m_seedComponentInputs[entry].push_back({seed,component,value});
    }
    if (m_seedComponentInputs.empty()) return StatusCode::FAILURE;
  }
'''

EXECUTE = '''
  // Each component invokes this version's actual finder, calibrator, transport,
  // measurement updates and isolated GainMatrixSmoother. Components never
  // compete in the event-level shared-hit selector. Only CKF posterior states
  // are observed in this opt-in mode; conversion/refit are not invoked.
  if (!m_seedComponentInputFile.value().empty()) {
    const auto inputs=m_seedComponentInputs.find(m_numberOfEvents-1);
    // A native seed finder can legitimately produce no seed for an event.
    // The component diagnostic input then has no rows for that local event;
    // this is an empty CKF call, not a malformed component file.  Preserve
    // the normal successful event flow and record an empty output collection.
    if (inputs==m_seedComponentInputs.end()) {
      m_audit->begin("component_event_end");
      m_audit->out << ",\\\"components\\\":0";
      m_audit->end();
      ATH_CHECK(trackContainer.record(std::move(outputTracks)));
      return StatusCode::SUCCESS;
    }
    for (const auto& component:inputs->second) {
      if (static_cast<size_t>(component.seed)>=initialParameters->size())
        return StatusCode::FAILURE;
      m_audit->inputComponent=component.component;
      const auto& value=component.input;
      Acts::BoundTrackParameters initial(
          Acts::Surface::makeShared<Acts::PlaneSurface>(value.transform),
          value.mean,value.covariance,initialParameters->at(component.seed).particleHypothesis());
      std::vector<const Tracker::FaserSCT_Cluster*> seedHits;
      if (Tool) seedHits=Tool->seeds()->at(component.seed).clusters;
      m_audit->parameters("seed",gctx,initial,component.seed,-1,seedHits);
      tracksTemp.clear();
      selectedTracks.clear();
      auto result=(*m_finder)(initial,options,calibrator,slAccessor,
                             Acts::MeasurementSelector(measurementSelectorCfg),tracksTemp);
      if (!result.ok()) {
        m_audit->begin("component_finder_error",component.seed);
        m_audit->out << ",\\"error\\":" << result.error().value(); m_audit->end();
        continue;
      }
      std::list<TrackInfo> candidates;
      for (const auto& track:result.value()) {
        candidates.push_back(track);
        m_audit->track("found",gctx,track,component.seed,track.index());
      }
      m_numberOfFittedTracks+=candidates.size();
      candidates.sort([](const TrackInfo& left,const TrackInfo& right) {
        if (left.nMeasurements>right.nMeasurements) return true;
        if (left.nMeasurements<right.nMeasurements) return false;
        return left.chi2<right.chi2;
      });
      while (!candidates.empty()) {
        const TrackInfo selected=candidates.front();
        candidates.erase(candidates.begin());
        auto dest=selectedTracks.getTrack(selectedTracks.addTrack());
        dest.copyFrom(tracksTemp.getTrack(selected.index),true);
        const auto status=Acts::GainMatrixSmoother()(
            gctx,selectedTracks.trackStateContainer(),dest.tipIndex(),*m_logger);
        if (!status.ok()) return StatusCode::FAILURE;
        m_audit->track("selected",gctx,dest,component.seed,selected.index);
        ++m_numberOfSelectedTracks;
        candidates.remove_if([&](const TrackInfo& other) {
          return other.nMeasurements<=12 || (other.clusterSet&selected.clusterSet).count()>6;
        });
      }
      m_audit->begin("component_end",component.seed);
      m_audit->out << ",\\"found\\":" << result.value().size()
                   << ",\\"selected\\":" << selectedTracks.size(); m_audit->end();
    }
    m_audit->inputComponent=-1;
    ATH_CHECK(trackContainer.record(std::move(outputTracks)));
    return StatusCode::SUCCESS;
  }

'''


def replace_once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new, 1)


def source(root, output):
    assert json.loads((root/'common_prior_ckf_v1/identity_validation.json').read_text())['all_exact']
    repo = Path(__file__).resolve().parents[1]
    output.mkdir(exist_ok=False)
    code = output/'build_code_v1'
    code.mkdir()
    shutil.copy2(repo/'scripts/build_native_calypso.py', code/'build_native_calypso.py')
    base = Path('/eos/home-x/xcheng/FASER/geoRotation_calypso')
    originals = {
        'beforeGeo':Path('/eos/home-x/xcheng/FASER/geoRotation_calypso_common_prior_before'),
        'afterGeo':Path('/eos/home-x/xcheng/FASER/geoRotation_calypso_charge_fixed_hits')}
    software = {'beforeGeo':root/'common_prior_ckf_v1/software_before',
                'afterGeo':root/'fixed_hit_process_v1/software'}
    jobs, manifests = [], {}
    package = Path('Tracking/Acts/FaserActsKalmanFilter/src')
    for variant, old in originals.items():
        target = base.parent/('geoRotation_calypso_components_'+variant)
        commit = subprocess.check_output(['git','-C',str(old),'rev-parse','HEAD'],text=True).strip()
        if (target/'.git').is_file():
            actual = subprocess.check_output(['git','-C',str(target),'rev-parse','HEAD'],text=True).strip()
            assert actual == commit
        else:
            subprocess.run(['git','-C',str(base),'worktree','add','--detach',str(target),commit],check=True)
        for name in ['FaserActsKalmanFilter','FaserActsGeometry']:
            relative = Path('Tracking/Acts')/name
            shutil.copytree(old/relative,target/relative,dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        changes = {}
        header = (old/package/'CKF2.h').read_text()
        changes['CKF2.h'] = replace_once(header,
            '  std::map<std::pair<size_t,int>,FullSeedInput> m_fullSeedInputs;',
            '  std::map<std::pair<size_t,int>,FullSeedInput> m_fullSeedInputs;'+HEADER)
        implementation = (old/package/'CKF2.cxx').read_text()
        initialize = INITIALIZE
        if variant == 'afterGeo':
            initialize = initialize.replace('    std::ifstream input(m_seedComponentInputFile.value());',
                '    if (!m_seedOverrideFile.value().empty() || m_fixedTruthHits) return StatusCode::FAILURE;\n'
                '    std::ifstream input(m_seedComponentInputFile.value());')
        implementation = replace_once(implementation,
            '  if (!m_fullSeedInputFile.value().empty()) {',
            initialize+'  if (!m_fullSeedInputFile.value().empty()) {')
        implementation = replace_once(implementation,
            '  // Perform the track finding for all initial parameters',
            EXECUTE+'  // Perform the track finding for all initial parameters')
        changes['CKF2.cxx'] = implementation
        audit = (old/package/'TrackingAudit.h').read_text()
        audit = replace_once(audit, '  std::string collection;',
                             '  std::string collection;\n  int inputComponent = -1;')
        changes['TrackingAudit.h'] = replace_once(audit, '        << ",\\"state\\":" << state;',
            '        << ",\\"state\\":" << state;\n'
            '    if (inputComponent>=0) out << ",\\"input_component\\":" << inputComponent;')
        increment = ''
        for name, text in changes.items():
            (target/package/name).write_text(text)
            filename = str(package/name)
            increment += 'diff --git a/'+filename+' b/'+filename+'\n'
            increment += ''.join(difflib.unified_diff((old/package/name).read_text().splitlines(True),
                text.splitlines(True),fromfile='a/'+filename,tofile='b/'+filename))
        patch = repo/'patches'/('0035_native_ckf_components_'+variant+'.patch')
        patch.write_text(increment)
        diff = subprocess.check_output(['git','-C',str(target),'diff','HEAD'],text=True)
        source_patch = output/(variant+'_tracked_source.patch')
        source_patch.write_text(diff)
        command = ['python3',str(code/'build_native_calypso.py'),'--source',str(target),
            '--output',str(output/('software_'+variant)),'--jobs','4',
            '--source-patch',str(source_patch),'--reuse-geometry',str(software[variant])]
        jobs.append(dict(kind='analysis',command=command,output=str(output/'jobs'/('build_'+variant))))
        manifests[variant] = dict(source=str(target),commit=commit,source_base=str(old),
            geometry_reused_from=str(software[variant]),increment=str(patch),tracked_patch=str(source_patch))
    dump(output/'build_tasks.json',jobs)
    dump(output/'manifest.json',dict(variants=manifests,
        scope='Actual native Calypso finder and measurement updates per finite Gaussian component, native selector within each component, private GainMatrixSmoother; CKF posterior only. Default behavior unchanged. No standalone transport or fit.'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--phase',choices=['source'],required=True)
    parser.add_argument('--output-name',default='common_prior_components_v2',
                        help='Output subdirectory below the stage06 result root')
    args = parser.parse_args()
    root = args.root.resolve()
    source(root,root/args.output_name)


if __name__ == '__main__':
    main()
