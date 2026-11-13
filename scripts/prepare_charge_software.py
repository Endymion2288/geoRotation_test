#!/usr/bin/env python3
"""Freeze isolated smoothing and read-only charge observers in Calypso worktrees."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess


def replace(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new)


def instrument(source, isolation):
    subprocess.run(['git', '-C', str(source), 'apply', '-'], input=isolation, text=True, check=True)
    folder = source/'Tracking/Acts/FaserActsKalmanFilter/src'
    path = folder/'CKF2.h'
    text = path.read_text()
    text = replace(text, '#include "xAODEventInfo/EventInfo.h"',
        '#include "xAODEventInfo/EventInfo.h"\n#include "ChargeDiagnostics.h"')
    marker = '  Gaudi::Property<std::string> m_actsLogging {this, "ActsLogging", "VERBOSE"};'
    text = replace(text, marker, marker+'''
  Gaudi::Property<std::string> m_chargeDiagnosticFile {this,"ChargeDiagnosticFile","","Read-only actual candidate lineage CSV"};
  SG::ReadHandleKey<TrackerSimDataCollection> m_chargeSimData {this,"ChargeSimData","SCT_SDO_Map"};
  CalypsoChargeDiagnostics m_chargeObserver;''')
    path.write_text(text)
    path = folder/'CKF2.cxx'
    text = path.read_text()
    text = replace(text, 'StatusCode CKF2::initialize() {', '''StatusCode CKF2::initialize() {
  ATH_CHECK(m_chargeSimData.initialize(!m_chargeDiagnosticFile.value().empty()));
  if (!m_chargeObserver.open(m_chargeDiagnosticFile.value())) return StatusCode::FAILURE;''')
    marker = '  TrackInfo::nClusters = sourceLinks->size();'
    text = replace(text, marker, '''  if (m_chargeObserver.enabled()) {
    SG::ReadHandle<TrackerSimDataCollection> sdo{m_chargeSimData,ctx};
    ATH_CHECK(sdo.isValid());
    m_chargeObserver.begin(m_numberOfEvents-1,ctx.eventID().run_number(),
                           ctx.eventID().event_number(),sdo.cptr());
  }
'''+marker)
    marker = '  for (std::size_t iseed = 0; iseed < (*initialParameters).size(); ++iseed) {'
    text = replace(text, marker, marker+'\n    m_chargeObserver.seed(iseed,(*initialParameters)[iseed]);')
    marker = '      destProxy.copyFrom(track, true);  // make sure we copy track states!'
    text = replace(text, marker, marker+'\n      m_chargeObserver.candidate(destProxy,iseed);')
    marker = '    if(not allTracks.empty()){'
    text = replace(text, marker, '    m_chargeObserver.selected(destProxy,selected.index);\n'+marker)
    marker = '        return ((p.nMeasurements <= 12) || ((p.clusterSet & selected.clusterSet).count() > 6));'
    text = replace(text, marker, '''        const bool rejected=(p.nMeasurements <= 12) || ((p.clusterSet & selected.clusterSet).count() > 6);
        if (rejected && m_chargeObserver.enabled()) m_chargeObserver.rejected(
            tracks.getTrack(p.index),selected.index,
            (p.nMeasurements<=12 ? 1 : 0) | ((p.clusterSet & selected.clusterSet).count()>6 ? 2 : 0));
        return rejected;''')
    marker = '    m_numberOfSelectedTracks++;'
    text = replace(text, marker, '    m_chargeObserver.converted(trk.get(),track.index());\n'+marker)
    marker = '      ATH_MSG_WARNING("Re-Fit either not performed or failed.");\n    }'
    text = replace(text, marker, marker+'''
    m_chargeObserver.final(outputTracks->back(),track.index(),outputTracks->size()-1,
                           refitStatus[track.index()]);''')
    text = replace(text, 'StatusCode CKF2::finalize() {', '''StatusCode CKF2::finalize() {
  if (!m_chargeObserver.close()) return StatusCode::FAILURE;''')
    path.write_text(text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    base = repo.parent
    after = base/'geoRotation_calypso_charge_after'
    before = base/'geoRotation_calypso_charge_before'
    relative = Path('Tracking/Acts/FaserActsKalmanFilter/src')
    header = after/relative/'ChargeDiagnostics.h'
    assert header.is_file()
    isolation = subprocess.check_output(['git','-C',str(base/'geoRotation_calypso_smoothing_fix'),
        'diff','88c457c1','1c06c7b','--',str(relative/'CKF2.cxx'),str(relative/'CKF2.h')],text=True)
    tasks = []
    code = args.root/'build_code'
    (code/'scripts').mkdir(parents=True,exist_ok=False)
    for name in ['build_native_calypso.py','build_slope_diagnostics.py']:
        shutil.copy2(repo/'scripts'/name,code/'scripts'/name)
    for variant, source, reuse in [
        ('beforeGeo',before,base/'geoRotation_test_results/stage05_geometry_comparison/software_before'),
        ('afterGeo',after,base/'geoRotation_test_results/stage04_smoothing_default/20261113_mr/software_clean')]:
        assert not subprocess.check_output(['git','-C',str(source),'diff','HEAD'],text=True)
        if source != after:
            shutil.copy2(header,source/relative/header.name)
        instrument(source,isolation)
        subprocess.run(['git','-C',str(source),'add','-N',str(relative/header.name)],check=True)
        diff = subprocess.check_output(['git','-C',str(source),'diff','HEAD'],text=True)
        patch = args.root/(variant+'_source.patch')
        patch.write_text(diff)
        command = ['python3',str(code/'scripts/build_native_calypso.py'),
            '--source',str(source),'--output',str(args.root/('software_'+variant)),
            '--jobs','4','--source-patch',str(patch),'--reuse-geometry',str(reuse)]
        tasks.append(dict(kind='analysis',command=command,output=str(args.root/'build_jobs'/variant)))
    observer = (base/'geoRotation_calypso'/relative/'CalypsoSlopeDiagnosticsAlg.cxx').read_text()
    marker = '  double m_upRaw[6]{},m_upF[6]{},m_minZ=0,m_boundAngles[2]{},m_angleCov[2][2]{};'
    observer = replace(observer,marker,marker+'\n  double m_qop=NAN,m_qopVariance=NAN;')
    marker = '  m_trackTree->Branch("charge",&m_charge,"charge/D");'
    observer = replace(observer,marker,marker+'''
  m_trackTree->Branch("qop_GeVinv",&m_qop,"qop_GeVinv/D");
  m_trackTree->Branch("qop_variance_GeVinv2",&m_qopVariance,"qop_variance_GeVinv2/D");''')
    marker = '    m_charge=upstream->charge();'
    assert marker in observer
    observer=observer.replace(marker,'''    m_qop=upstream->parameters()[4]*1000.;
    m_qopVariance=upstream->covariance() ? (*upstream->covariance())(4,4)*1e6 : NAN;
'''+marker)
    (after/relative/'CalypsoSlopeDiagnosticsAlg.cxx').write_text(observer)
    tasks.append(dict(kind='analysis',command=['python3',str(code/'scripts/build_slope_diagnostics.py'),
        '--source',str(after),'--output',str(args.root/'observer')],output=str(args.root/'build_jobs/observer')))
    (args.root/'build_tasks.json').write_text(json.dumps(tasks,indent=2)+'\n')
    print(args.root/'build_tasks.json')


if __name__ == '__main__':
    main()
