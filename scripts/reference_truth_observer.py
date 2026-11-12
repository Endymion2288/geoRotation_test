"""Read native Calypso production truth and output order without altering tracking."""
import json
import math
import AthenaPython.PyAthena as PyAthena


# Use the overload exported by the actual loaded Calypso tracking library. The
# declaration avoids importing ACTS template headers into the ROOT interpreter.
FINAL_PURITY_CPP = r'''
#ifndef HEPMC3
#define HEPMC3 1
#endif
#include "EventPrimitives/EventPrimitives.h"
#include "GeoPrimitives/GeoPrimitives.h"
#include "TrkTrack/TrackCollection.h"
#include "TrackerRIO_OnTrack/FaserSCT_ClusterOnTrack.h"
#include "TrackerPrepRawData/FaserSCT_Cluster.h"
#include "TrackerSimData/TrackerSimDataCollection.h"
#include "GaudiKernel/Bootstrap.h"
#include "StoreGate/StoreGateSvc.h"
#include <nlohmann/json.hpp>
struct ParticleHitCount { int particleId; size_t hitCount; };
void identifyContributingParticles(const TrackerSimDataCollection&,
    const std::vector<const Tracker::FaserSCT_Cluster*>, std::vector<ParticleHitCount>&);
std::string referenceFinalPurityJSON() {
  StoreGateSvc* store=nullptr;
  if(Gaudi::svcLocator()->service("StoreGateSvc",store).isFailure())
    throw std::runtime_error("No StoreGateSvc");
  const TrackCollection* tracks=nullptr;
  const TrackerSimDataCollection* sdo=nullptr;
  if(store->retrieve(tracks,"CKFTrackCollectionWithoutIFT").isFailure() ||
     store->retrieve(sdo,"SCT_SDO_Map").isFailure())
    throw std::runtime_error("Missing final tracks or SDO map");
  auto out=nlohmann::json::array();
  for(const auto* track:*tracks) {
    std::vector<const Tracker::FaserSCT_Cluster*> clusters;
    if(track->measurementsOnTrack()) for(const auto* m:*track->measurementsOnTrack()) {
      const auto* cot=dynamic_cast<const Tracker::FaserSCT_ClusterOnTrack*>(m);
      if(cot) clusters.push_back(cot->prepRawData());
    }
    std::vector<ParticleHitCount> counts;
    identifyContributingParticles(*sdo,clusters,counts);
    size_t majority=0; int barcode=0;
    for(const auto& c:counts) if(c.hitCount>majority ||
        (c.hitCount==majority && c.particleId<barcode)) {
      majority=c.hitCount; barcode=c.particleId;
    }
    out.push_back({{"clusters",clusters.size()},{"majority_hits",majority},
        {"majority_barcode",barcode}});
  }
  return out.dump();
}
'''


class ReferenceTruthObserver(PyAthena.Alg):
    def __init__(self, name, skip=0, rotated=True, native_final_purity=False):
        super().__init__(name=name)
        self.skip = skip
        self.rotated = rotated
        self.native_final_purity = native_final_purity

    def initialize(self):
        self.sg = PyAthena.py_svc('StoreGateSvc')
        self.stream = open('production_truth.jsonl', 'w', buffering=1024*1024)
        self.ordinal = 0
        self.purity_declared = False
        return PyAthena.StatusCode.Success

    def execute(self):
        ei = self.sg.retrieve('xAOD::EventInfo', 'EventInfo')
        if self.native_final_purity:
            import ROOT
            if not self.purity_declared:
                # Same HepMC3 headers/ABI as the Athena 24.0.41 compile recipe.
                ROOT.gInterpreter.AddIncludePath('/cvmfs/atlas.cern.ch/repo/sw/software/24.0/sw/lcg/releases/LCG_104d_ATLAS_7/hepmc3/3.2.7/x86_64-el9-gcc13-opt/include')
                if not ROOT.gInterpreter.Declare(FINAL_PURITY_CPP):
                    return PyAthena.StatusCode.Failure
                self.purity_declared = True
            final = json.loads(str(ROOT.referenceFinalPurityJSON()))
            row = dict(source_entry=self.skip+self.ordinal, run=int(ei.runNumber()),
                event=int(ei.eventNumber()), ntracks=len(final), final=final)
            self.stream.write(json.dumps(row, separators=(',', ':'))+'\n')
            self.ordinal += 1
            return PyAthena.StatusCode.Success
        truth = self.sg.retrieve('xAOD::TruthParticleContainer', 'TruthParticles')
        hits = self.sg.retrieve('FaserSiHitCollection', 'SCT_Hits')
        tracks = self.sg.retrieve('TrackCollection', 'CKFTrackCollectionWithoutIFT')
        hit_codes = {int(hit.trackNumber()) for hit in hits}
        particles = {}
        for particle in truth:
            barcode = int(particle.barcode())
            particles[barcode] = dict(p_MeV=math.sqrt(particle.px()**2+particle.py()**2+particle.pz()**2),
                charge=float(particle.charge()), pdg=int(particle.pdgId()), status=int(particle.status()),
                has_sihit=barcode in hit_codes)
        final = []
        for track in tracks:
            params = list(track.trackParameters())
            quality = track.fitQuality()
            def radius(parameter):
                pos = parameter.position()
                return math.hypot(pos.y(), pos.z()) if self.rotated else math.hypot(pos.x(), pos.y())
            final.append(dict(chi2=float(quality.chiSquared()), ndof=float(quality.numberDoF()),
                measurements=int(track.measurementsOnTrack().size()),
                rmax_mm=max((radius(q) for q in params), default=None)))
        row = dict(source_entry=self.skip+self.ordinal, run=int(ei.runNumber()), event=int(ei.eventNumber()),
            ntracks=len(final), truth=particles, final=final)
        self.stream.write(json.dumps(row, allow_nan=False, separators=(',', ':'))+'\n')
        self.ordinal += 1
        return PyAthena.StatusCode.Success

    def finalize(self):
        self.stream.close()
        return PyAthena.StatusCode.Success
