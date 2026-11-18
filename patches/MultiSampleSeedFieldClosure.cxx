#include "GeoPrimitives/GeoPrimitives.h"
#include "FaserActsGeometryInterfaces/IFaserActsExtrapolationTool.h"
#include "FaserActsGeometryInterfaces/IFaserActsTrackingGeometryTool.h"
#include "TrackerReadoutGeometry/SCT_DetectorManager.h"
#include "TrackerIdentifier/FaserSCT_ID.h"
#include "GeneratorObjects/McEventCollection.h"
#include "AtlasHepMC/GenParticle.h"
#include "AtlasHepMC/GenVertex.h"
#include "HepMC3/Units.h"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "CircleFit.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ThreadLocalContext.h"
#include "StoreGate/StoreGateSvc.h"
#include <nlohmann/json.hpp>
#include <cmath>
#include <set>
#include <unordered_map>

struct MultiSampleSeedFieldClosure {
  ToolHandle<IFaserActsExtrapolationTool> tool{
      "FaserActsExtrapolationTool/MultiSampleSeedFieldExtrapolator"};
  Acts::RotationMatrix3 rotation;
  std::array<double, 4> centers{};

  static double pdgCharge(int pdg) {
    const int a = std::abs(pdg);
    if (a == 11 || a == 13 || a == 15) return pdg > 0 ? -1.0 : 1.0;
    if (a == 211 || a == 321 || a == 2212) return pdg > 0 ? 1.0 : -1.0;
    return 0.0;
  }

  MultiSampleSeedFieldClosure() {
    if (tool.retrieve().isFailure()) throw std::runtime_error("No native extrapolation tool");
    rotation << 0, 0, 1, 0, 1, 0, -1, 0, 0;
    StoreGateSvc* detector = nullptr;
    Gaudi::svcLocator()->service("DetectorStore", detector).ignore();
    const TrackerDD::SCT_DetectorManager* manager = nullptr;
    const FaserSCT_ID* ids = nullptr;
    if (detector->retrieve(manager, "SCT").isFailure() ||
        detector->retrieve(ids, "FaserSCT_ID").isFailure()) {
      throw std::runtime_error("No geometry");
    }
    std::array<int, 4> count{};
    for (const auto* element : *manager->getDetectorElementCollection()) {
      int station = ids->station(element->identify());
      centers[station] += element->center().z();
      ++count[station];
    }
    for (int i = 0; i < 4; ++i) centers[i] /= count[i];
  }

  std::string execute() {
    StoreGateSvc* store = nullptr;
    Gaudi::svcLocator()->service("StoreGateSvc", store).ignore();
    const McEventCollection* events = nullptr;
    if (store->retrieve(events, "TruthEvent").isFailure()) throw std::runtime_error("No truth");
    const auto& ctx = Gaudi::Hive::currentContext();
    auto gctx = tool->trackingGeometryTool()->getGeometryContext(ctx).context();
    auto rows = nlohmann::json::array();

    Acts::RotationMatrix3 planeRotation;
    planeRotation.col(0) = Acts::Vector3(0, 1, 0);
    planeRotation.col(1) = Acts::Vector3(0, 0, -1);
    planeRotation.col(2) = Acts::Vector3(-1, 0, 0);
    auto plane = [&](double x) {
      return Acts::Surface::makeShared<Acts::PlaneSurface>(
          Acts::Transform3(Acts::Translation3(x, 0, 0) * planeRotation));
    };

    for (const auto* event : *events) {
      if (!event) continue;
      const double ps = event->momentum_unit() == HepMC3::Units::GEV ? 1. : .001;
      const double xs = event->length_unit() == HepMC3::Units::CM ? 10. : 1.;

      struct StateSnap {
        Acts::Vector3 posF;
        Acts::Vector3 momF;
      };
      std::unordered_map<int, StateSnap> latestPreS1;
      std::unordered_map<int, int> primaryPdg;

      bool hasGenChargedPrimary = false;
      for (const auto& particle : event->particles()) {
        if (!particle) continue;
        const int bc = HepMC::barcode(particle);
        const int a = std::abs(particle->pdg_id());
        if (bc > 0 && bc < 200000 && particle->status() == 1 &&
            (a == 11 || a == 13 || a == 211 || a == 321 || a == 2212)) {
          hasGenChargedPrimary = true;
          break;
        }
      }

      for (const auto& particle : event->particles()) {
        if (!particle || !particle->production_vertex()) continue;
        const int bc = HepMC::barcode(particle);
        const int pdg = particle->pdg_id();
        const int a = std::abs(pdg);
        if (bc <= 0 || (a != 11 && a != 13 && a != 211 && a != 321 && a != 2212)) continue;
        const auto v = particle->production_vertex()->position();
        const auto m = particle->momentum();
        const Acts::Vector3 posF(v.x() * xs, v.y() * xs, v.z() * xs);
        const Acts::Vector3 momF(m.px() * ps, m.py() * ps, m.pz() * ps);
        if (bc < 1000000) {
          bool acceptSignal = false;
          if (hasGenChargedPrimary) {
            acceptSignal = (bc < 200000 && particle->status() == 1);
          } else if (posF.z() < centers[1] - 5.0 && momF.norm() >= 5.0) {
            for (const auto& inP : particle->production_vertex()->particles_in()) {
              if (inP && (HepMC::barcode(inP) % 1000000) < 200000 && pdgCharge(inP->pdg_id()) == 0.0) {
                acceptSignal = true;
                break;
              }
            }
          }
          if (acceptSignal) {
            primaryPdg[bc] = pdg;
            if (posF.z() < centers[1] - 5.0) {
              latestPreS1[bc] = {posF, momF};
            }
          }
        } else {
          const int baseBc = bc % 1000000;
          if (primaryPdg.count(baseBc) && primaryPdg[baseBc] == pdg && posF.z() < centers[1] - 5.0) {
            if (!latestPreS1.count(baseBc) || posF.z() > latestPreS1[baseBc].posF.z()) {
              latestPreS1[baseBc] = {posF, momF};
            }
          }
        }
      }

      for (const auto& kv : latestPreS1) {
        const int bc = kv.first;
        const int pdg = primaryPdg[bc];
        const double charge = pdgCharge(pdg);
        // Start from max(z_F, -50 mm) along straight line if particle originates far upstream in field-free FASERnu
        Acts::Vector3 posF = kv.second.posF;
        Acts::Vector3 momF = kv.second.momF;
        const double p = momF.norm();
        if (p < 5.0 || momF.z() <= 0.0 || charge == 0.0) continue;
        if (posF.z() < -50.0) {
          const double dz = -50.0 - posF.z();
          posF += (dz / momF.z()) * momF;
        }
        const Acts::Vector3 posA = rotation * posF;
        const Acts::Vector3 momA = rotation * momF;
        auto startSurface = plane(posA.x());
        auto local = startSurface->globalToLocal(gctx, posA);
        if (!local.ok()) continue;
        Acts::BoundVector mean = Acts::BoundVector::Zero();
        mean[0] = local.value()[0];
        mean[1] = local.value()[1];
        mean[2] = std::atan2(momA.y(), momA.x());
        mean[3] = std::acos(momA.z() / p);
        mean[4] = charge / p;
        Acts::BoundTrackParameters start(startSurface, mean, std::nullopt, Acts::ParticleHypothesis::muon());

        std::vector<Amg::Vector3D> points;
        std::array<Acts::Vector3, 4> stPosF{};
        std::array<Acts::Vector3, 4> stDirF{};
        bool complete = true;
        for (int st = 1; st <= 3; ++st) {
          auto target = plane(centers[st]);
          auto result = tool->propagate(ctx, start, *target);
          if (!result) {
            complete = false;
            break;
          }
          stPosF[st] = rotation.transpose() * result->position(gctx);
          stDirF[st] = rotation.transpose() * result->direction();
          // Require all 3 stations within physical magnet/SCT aperture r < 100 mm
          if (std::hypot(stPosF[st].x(), stPosF[st].y()) > 100.0) {
            complete = false;
            break;
          }
          points.push_back(stPosF[st]);
          points.push_back(stPosF[st] - 30.0 * stDirF[st]);
          points.push_back(stPosF[st] + 30.0 * stDirF[st]);
        }
        if (!complete) continue;
        const auto circle = CircleFit::circleFit(CircleFit::CircleData(points));
        if (!std::isfinite(circle.r) || circle.r <= 0.0) continue;
        const double bEff = p / (circle.r * 0.3 * 0.001);
        const double p055 = circle.r * 0.3 * 0.55 * 0.001;
        const double p04737 = circle.r * 0.3 * 0.473662 * 0.001;
        const double rS1 = std::hypot(stPosF[1].x(), stPosF[1].y());
        const double rMax = std::max({rS1,
                                      std::hypot(stPosF[2].x(), stPosF[2].y()),
                                      std::hypot(stPosF[3].x(), stPosF[3].y())});
        const double txMrad = 1000.0 * (stDirF[1].x() / stDirF[1].z());
        const double tyMrad = 1000.0 * (stDirF[1].y() / stDirF[1].z());
        rows.push_back({
            {"barcode", bc},
            {"pdg", pdg},
            {"charge", charge},
            {"p_GeV", p},
            {"r_s1_mm", rS1},
            {"r_max_mm", rMax},
            {"tx_s1_mrad", txMrad},
            {"ty_s1_mrad", tyMrad},
            {"slope_s1_mrad", std::hypot(txMrad, tyMrad)},
            {"circle_r_mm", circle.r},
            {"effective_B_T", bEff},
            {"res_B055", p055 / p - 1.0},
            {"res_B0473662", p04737 / p - 1.0},
        });
      }
    }
    return rows.dump();
  }
};

extern "C" void* multi_sample_seed_field_create() {
  return new MultiSampleSeedFieldClosure();
}
extern "C" const char* multi_sample_seed_field_execute(void* observer) {
  static thread_local std::string result;
  result = static_cast<MultiSampleSeedFieldClosure*>(observer)->execute();
  return result.c_str();
}
extern "C" void multi_sample_seed_field_delete(void* observer) {
  delete static_cast<MultiSampleSeedFieldClosure*>(observer);
}
