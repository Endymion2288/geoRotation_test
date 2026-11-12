#!/usr/bin/env python3
"""Read actual GeoModel plate positions, strip pitch, stereo and local bases."""
import argparse
import json
from pathlib import Path


CPP = r'''
#include "TrackerReadoutGeometry/SCT_DetectorManager.h"
#include "TrackerReadoutGeometry/SiDetectorElement.h"
#include "TrackerIdentifier/FaserSCT_ID.h"
#include "ScintReadoutGeometry/VetoNuDetectorManager.h"
#include "ScintReadoutGeometry/ScintDetectorElement.h"
#include "GaudiKernel/Bootstrap.h"
#include "StoreGate/StoreGateSvc.h"
#include <nlohmann/json.hpp>
#include <stdexcept>
std::string calypsoSlopeGeometryJSON() {
  StoreGateSvc* store=nullptr;
  if(Gaudi::svcLocator()->service("DetectorStore",store).isFailure()) throw std::runtime_error("No DetectorStore");
  const TrackerDD::SCT_DetectorManager* tracker=nullptr;
  const ScintDD::VetoNuDetectorManager* veto=nullptr;
  const FaserSCT_ID* ids=nullptr;
  if(store->retrieve(tracker,"SCT").isFailure() || store->retrieve(veto,"VetoNu").isFailure() ||
     store->retrieve(ids,"FaserSCT_ID").isFailure()) throw std::runtime_error("Missing actual detector manager");
  nlohmann::json out;
  out["frame"]="GeoModel/G4 physical FASER frame";
  out["tracker"]=nlohmann::json::array();out["VetoNu"]=nlohmann::json::array();
  for(const auto* e:*tracker->getDetectorElementCollection()) {
    if(!e) continue;
    const auto id=e->identify();const auto& c=e->center();const auto& m=e->transform().matrix();
    nlohmann::json row={{"id",id.get_compact()},{"station",ids->station(id)},{"layer",ids->layer(id)},
      {"eta",ids->eta_module(id)},{"phi",ids->phi_module(id)},{"side",ids->side(id)},
      {"center_F_mm",{c.x(),c.y(),c.z()}},{"pitch_mm",e->phiPitch()},
      {"sin_stereo",e->sinStereo()},{"thickness_mm",e->thickness()}};
    row["local_to_F"]=nlohmann::json::array();
    for(int i=0;i<4;++i) row["local_to_F"].push_back({m(i,0),m(i,1),m(i,2),m(i,3)});
    out["tracker"].push_back(row);
  }
  for(const auto* e:*veto->getDetectorElementCollection()) {
    if(!e) continue;
    const auto& c=e->center();const auto& n=e->normal();
    out["VetoNu"].push_back({{"id",e->identify().get_compact()},
      {"center_F_mm",{c.x(),c.y(),c.z()}},{"normal_F",{n.x(),n.y(),n.z()}},
      {"thickness_mm",e->thickness()}});
  }
  return out.dump();
}
'''


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',required=True)
    parser.add_argument('--output',type=Path,default=Path('slope_geometry.json'))
    args=parser.parse_args()
    from AthenaCommon.Configurable import Configurable
    from CalypsoConfiguration.AllConfigFlags import initConfigFlags
    from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    from FaserGeoModel.FaserGeoModelConfig import FaserGeometryCfg
    import AthenaPython.PyAthena as PyAthena
    import ROOT
    Configurable.configurableRun3Behavior=True
    flags=initConfigFlags()
    flags.Input.Files=[args.input];flags.Input.isMC=True;flags.Input.ProjectName='mc24'
    flags.IOVDb.DatabaseInstance='OFLP200';flags.IOVDb.GlobalTag='OFLCOND-FASER-05'
    flags.GeoModel.FaserVersion='FASERNU-04'
    flags.Concurrency.NumThreads=0;flags.Concurrency.NumConcurrentEvents=1
    flags.lock()
    class Probe(PyAthena.Alg):
        def execute(self):
            if not ROOT.gInterpreter.Declare(CPP):return PyAthena.StatusCode.Failure
            data=json.loads(str(ROOT.calypsoSlopeGeometryJSON()))
            data.update(geometry='FASERNU-04',conditions='OFLCOND-FASER-05',
                        reference_planes_z_F_mm={'st1':0.,'VetoNu':-3112.})
            args.output.write_text(json.dumps(data,indent=2)+'\n')
            return PyAthena.StatusCode.Success
    acc=MainServicesCfg(flags);acc.merge(PoolReadCfg(flags));acc.merge(FaserGeometryCfg(flags))
    acc.addEventAlgo(Probe('ActualSlopeGeometryProbe'))
    return int(acc.run(1).isFailure())


if __name__=='__main__':raise SystemExit(main())
