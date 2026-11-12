"""Read material lengths from the actual initialized Calypso Geant4 geometry."""
import json
import os
from pathlib import Path
import shlex


CPP = r'''
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include <nlohmann/json.hpp>
#include <map>
std::string calypsoActualG4MaterialPaths() {
  auto* world=G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
  if(!world) throw std::runtime_error("Actual Calypso G4 world has not been initialized");
  nlohmann::json out;
  out["geometry"]="FASERNU-04";out["conditions"]="OFLCOND-FASER-05";
  out["definition"]="Straight material rays in the actual initialized Calypso Geant4 geometry, physical FASER frame";
  out["paths"]=nlohmann::json::array();
  for(auto xy:std::vector<std::pair<double,double>>{{0,0},{10,-21},{-30,-30},{-30,30},{30,-30},{30,30}}) {
    G4Navigator navigator;navigator.SetWorldVolume(world);
    G4ThreeVector position(xy.first,xy.second,-3112.),direction(0,0,1);
    std::map<std::string,std::array<double,3>> sums;
    int boundaries=0;
    while(position.z()<0. && boundaries<10000) {
      auto* volume=navigator.LocateGlobalPointAndSetup(position,&direction,false);
      if(!volume) throw std::runtime_error("Ray left actual G4 world");
      auto* material=volume->GetLogicalVolume()->GetMaterial();
      double safety=0.;
      double step=navigator.ComputeStep(position,direction,-position.z(),safety);
      if(!std::isfinite(step) || step>-position.z()) step=-position.z();
      step=std::max(step,1e-6);step=std::min(step,-position.z());
      auto& sum=sums[material->GetName()];sum[0]+=step;sum[1]+=step/material->GetRadlen();
      sum[2]=material->GetDensity()/(CLHEP::g/CLHEP::cm3);
      position+=step*direction;++boundaries;
    }
    if(position.z()<-1e-5) throw std::runtime_error("G4 boundary budget exhausted");
    nlohmann::json materials=nlohmann::json::array();double x0=0.;
    for(const auto& item:sums) {
      materials.push_back({{"material",item.first},{"length_mm",item.second[0]},
        {"X0",item.second[1]},{"density_g_cm3",item.second[2]}});x0+=item.second[1];
    }
    out["paths"].push_back({{"x_F_mm",xy.first},{"y_F_mm",xy.second},{"z_F_mm",{-3112.,0.}},
      {"total_X0",x0},{"boundaries",boundaries},{"materials",materials}});
  }
  return out.dump();
}
'''


def add_material_probe(acc, output):
    import AthenaPython.PyAthena as PyAthena
    import ROOT
    source=Path('/eos/home-x/xcheng/FASER/geoRotation_calypso')
    records=json.loads((source/'build/compile_commands.json').read_text())
    record=next(r for r in records if '/G4FaserAlg/' in r['file'] and r['file'].endswith('.cxx'))
    tokens=shlex.split(record['command'])
    for i,token in enumerate(tokens):
        if token in ['-I','-isystem']:
            ROOT.gInterpreter.AddIncludePath(tokens[i+1])
        elif token.startswith('-I'):
            ROOT.gInterpreter.AddIncludePath(token[2:])
    class Probe(PyAthena.Alg):
        def execute(self):
            if not ROOT.gInterpreter.Declare(CPP):return PyAthena.StatusCode.Failure
            data=json.loads(str(ROOT.calypsoActualG4MaterialPaths()))
            Path(output).write_text(json.dumps(data,indent=2)+'\n')
            return PyAthena.StatusCode.Success
    acc.addEventAlgo(Probe('ActualG4MaterialPathProbe'))
