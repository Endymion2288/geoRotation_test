"""Configure the real SCT alignment conditions chain for nominal controls."""


def configure_alignment(acc, nominal, audit_file=None):
    from AthenaConfiguration.ComponentFactory import CompFactory
    if nominal:
        acc.addCondAlgo(CompFactory.FaserSCT_AlignCondAlg(
            'FaserSCT_AlignCondAlg', UseNominalAlignment=True,
            AlignmentAuditFile=audit_file or ''))
        acc.addCondAlgo(CompFactory.FaserSCT_DetectorElementCondAlg(
            'FaserSCT_DetectorElementCondAlg'))
    elif audit_file:
        acc.getCondAlgo('FaserSCT_AlignCondAlg').AlignmentAuditFile = audit_file


def write_simulation_geometry_audit(filename):
    """Read the actual retained GeoModel used to construct the G4 geometry."""
    import ROOT
    code = r'''
#include "TrackerReadoutGeometry/SCT_DetectorManager.h"
#include "TrackerReadoutGeometry/SiDetectorElement.h"
#include "GaudiKernel/Bootstrap.h"
#include "StoreGate/StoreGateSvc.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>
void calypsoSimulationAlignmentAudit(const std::string& filename) {
  StoreGateSvc* store=nullptr;
  const TrackerDD::SCT_DetectorManager* manager=nullptr;
  if(Gaudi::svcLocator()->service("DetectorStore",store).isFailure() ||
     store->retrieve(manager,"SCT").isFailure()) throw std::runtime_error("No actual SCT GeoModel");
  std::ofstream audit(filename);
  if(!audit) throw std::runtime_error("Cannot write actual simulation geometry audit");
  audit << "wafer_id";
  for(const char* label: {"actual","nominal"})
    for(int i=0;i<4;++i) for(int j=0;j<4;++j) audit << ',' << label << i << j;
  audit << '\n' << std::setprecision(17);
  for(const auto* element:*manager->getDetectorElementCollection()) {
    const Amg::Transform3D actual=element->transform(), nominal=element->defTransform();
    audit << element->identify().get_compact();
    for(const auto* transform: {&actual,&nominal})
      for(int i=0;i<4;++i) for(int j=0;j<4;++j) audit << ',' << transform->matrix()(i,j);
    audit << '\n';
  }
}
'''
    assert ROOT.gInterpreter.Declare(code)
    ROOT.calypsoSimulationAlignmentAudit(str(filename))


def add_simulation_geometry_audit(acc, filename):
    import AthenaPython.PyAthena as PyAthena
    class Probe(PyAthena.Alg):
        def initialize(self):
            self.done = False
            return PyAthena.StatusCode.Success

        def execute(self):
            if not self.done:
                write_simulation_geometry_audit(filename)
                self.done = True
            return PyAthena.StatusCode.Success
    acc.addEventAlgo(Probe('ActualSimulationAlignmentAudit'))
