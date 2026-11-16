#!/usr/bin/env python3
"""Read recorded HepMC continuation vertices via the actual Athena POOL converter."""
import argparse
import json
from pathlib import Path


CPP = r'''
#include "GeneratorObjects/McEventCollection.h"
#include "AtlasHepMC/GenParticle.h"
#include "AtlasHepMC/GenVertex.h"
#include "HepMC3/Units.h"
#include "GaudiKernel/Bootstrap.h"
#include "StoreGate/StoreGateSvc.h"
#include <nlohmann/json.hpp>
std::string latestRecordedVertices() {
  StoreGateSvc* store=nullptr;
  if(Gaudi::svcLocator()->service("StoreGateSvc",store).isFailure())
    throw std::runtime_error("No event store");
  const McEventCollection* events=nullptr;
  if(store->retrieve(events,"TruthEvent").isFailure())
    throw std::runtime_error("No recorded TruthEvent");
  auto out=nlohmann::json::array();
  for(const auto* event:*events) {
    const double ps=event->momentum_unit()==HepMC3::Units::GEV ? 1. : .001;
    const double xs=event->length_unit()==HepMC3::Units::CM ? 10. : 1.;
    for(const auto& p:event->particles()) {
      const int barcode=HepMC::barcode(p);
      if(std::abs(p->pdg_id())!=11 || barcode<=0 || barcode%1000000>=200000 || !p->production_vertex()) continue;
      const auto& m=p->momentum();
      out.push_back({barcode,p->pdg_id(),p->production_vertex()->position().z()*xs,
                     m.rho()*ps,m.pz()*ps});
    }
  }
  return out.dump();
}
'''


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args()
    task=json.loads(args.task.read_text())
    import ROOT
    import numpy as np
    import AthenaPython.PyAthena as PyAthena
    from AthenaCommon.Configurable import Configurable
    from CalypsoConfiguration.AllConfigFlags import initConfigFlags
    from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    class Observer(PyAthena.Alg):
        def initialize(self):
            ROOT.gInterpreter.AddIncludePath('/cvmfs/atlas.cern.ch/repo/sw/software/24.0/sw/lcg/releases/LCG_104d_ATLAS_7/hepmc3/3.2.7/x86_64-el9-gcc13-opt/include')
            ROOT.gInterpreter.Declare('#define HEPMC3 1\n'+CPP)
            self.rows=[];self.entry=0
            return PyAthena.StatusCode.Success
        def execute(self):
            self.rows.extend([[self.entry,*row] for row in json.loads(str(ROOT.latestRecordedVertices()))])
            self.entry+=1
            return PyAthena.StatusCode.Success
        def finalize(self):
            assert self.entry==task['events']
            np.savez_compressed(Path(task['output'])/'vertices.npz',values=np.asarray(self.rows),
                fields=['entry','barcode','pdg','z_F_mm','p_GeV','pz_GeV'])
            return PyAthena.StatusCode.Success
    Configurable.configurableRun3Behavior=True
    flags=initConfigFlags()
    flags.Input.Files=[task['input']];flags.Input.isMC=True;flags.Input.ProjectName='mc24'
    flags.IOVDb.DatabaseInstance='OFLP200'
    flags.GeoModel.FaserVersion='FASERNU-04'
    flags.IOVDb.GlobalTag='OFLCOND-FASER-05'
    flags.lock()
    acc=MainServicesCfg(flags);acc.merge(PoolReadCfg(flags));acc.addEventAlgo(Observer('RecordedVertices'))
    status=acc.run(maxEvents=task['events'])
    return 0 if status.isSuccess() else 1


if __name__=='__main__':
    raise SystemExit(main())
