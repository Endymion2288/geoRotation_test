// Fast binary C++ validation observer for multi-sample Calypso seed, CKF, and refit studies.
#ifndef FASER_LATEST_FAST_VALIDATION_AUDIT_H
#define FASER_LATEST_FAST_VALIDATION_AUDIT_H

#include "CircleFitTrackSeedTool.h"
#include "KalmanFitterTool.h"
#include "TrackingTruthUtils/TrackClassification.h"
#include "GeneratorObjects/McEventCollection.h"
#include "AtlasHepMC/GenParticle.h"
#include "AtlasHepMC/GenVertex.h"
#include "HepMC3/Units.h"
#include "TrackerSimEvent/FaserSiHitCollection.h"
#include "TrackerReadoutGeometry/SCT_DetectorManager.h"
#include "TrackerReadoutGeometry/SiDetectorElement.h"
#include "TrackerIdentifier/FaserSCT_ID.h"
#include "ActsInterop/IndexSourceLink.h"
#include "Acts/EventData/TrackParameters.hpp"
#include "TrkTrack/Track.h"
#include "TrkTrack/TrackCollection.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

class LatestFastValidationAudit {
 public:
  static constexpr size_t kEventCols = 16;
  static constexpr size_t kSignalCols = 38;
  static constexpr size_t kSeedCols = 34;
  static constexpr size_t kTrackCols = 78;

  std::ofstream evOut, sigOut, seedOut, trkOut;
  const FaserSCT_ID* idHelper = nullptr;
  const TrackerDD::SCT_DetectorManager* detMgr = nullptr;

  explicit LatestFastValidationAudit(const std::string& prefix) {
    evOut.open(prefix + "_events.bin", std::ios::binary | std::ios::trunc);
    sigOut.open(prefix + "_signals.bin", std::ios::binary | std::ios::trunc);
    seedOut.open(prefix + "_seeds.bin", std::ios::binary | std::ios::trunc);
    trkOut.open(prefix + "_tracks.bin", std::ios::binary | std::ios::trunc);
    std::ofstream schema(prefix + "_schema.json", std::ios::trunc);
    schema << "{\n"
           << "  \"event_cols\": [\"entry\",\"run\",\"event\",\"n_signal\",\"n_fiducial\","
              "\"raw_s1\",\"raw_s2\",\"raw_s3\",\"gb_s1\",\"gb_s2\",\"gb_s3\","
              "\"raw_comb\",\"raw_nan_chi2\",\"n_seeds\",\"n_seeds_nan_chi2\",\"n_tracks\"],\n"
           << "  \"signal_cols\": [\"entry\",\"barcode\",\"pdg\",\"charge\",\"status\",\"fiducial\","
              "\"hits_s0\",\"hits_s1\",\"hits_s2\",\"hits_s3\",\"r_vtx_mm\",\"r_max_s123_mm\","
              "\"vx_F\",\"vy_F\",\"vz_F\",\"tx_F\",\"ty_F\","
              "\"p_prod_GeV\",\"p_st1_GeV\",\"p_pre_s2_GeV\",\"p_min_trk_GeV\","
              "\"s1_hit_valid\",\"s1_geo_id\",\"s1_loc0\",\"s1_loc1\",\"s1_phi_A\",\"s1_theta_A\","
              "\"n_seed_match\",\"best_seed_purity\",\"best_seed_p_GeV\","
              "\"n_ckf_match\",\"best_ckf_purity\",\"best_ckf_p_GeV\",\"best_ckf_q\",\"best_ckf_chi2\",\"best_ckf_ndf\","
              "\"best_refit_ok\",\"best_refit_p_GeV\"],\n"
           << "  \"seed_cols\": [\"entry\",\"iseed\",\"stations\",\"n_clusters\",\"chi2\",\"red_chi2\","
              "\"chi2_strip\",\"ndf_strip\",\"chi2_is_nan\",\"min_radicand\",\"r_mm\","
              "\"p_seed_GeV\",\"q_seed\",\"loc0\",\"loc1\",\"phi_A\",\"theta_A\","
              "\"max_shared_prev\",\"shared_union_prev\",\"shared_s23_prev\","
              "\"all_same_bc\",\"opp_primary_mix\",\"maj_bc\",\"purity\","
              "\"truth_pdg\",\"truth_q\",\"p_prod_GeV\",\"p_st1_GeV\",\"p_pre_s2_GeV\","
              "\"seed_truth_valid\",\"t_loc0\",\"t_loc1\",\"t_phi_A\",\"t_theta_A\"],\n"
           << "  \"track_cols\": [\"entry\",\"itrk\",\"iseed\",\"n_meas\",\"holes\",\"outliers\","
              "\"chi2_ckf\",\"ndf_ckf\",\"rmax_ckf\",\"refit_ok\",\"chi2_refit\",\"ndf_refit\",\"rmax_refit\","
              "\"seed_maj_bc\",\"trk_maj_bc\",\"trk_purity\",\"truth_pdg\",\"truth_q\","
              "\"p_prod_GeV\",\"p_st1_GeV\",\"p_pre_s2_GeV\",\"p_min_trk_GeV\","
              "\"ckf_tgt_loc0\",\"ckf_tgt_loc1\",\"ckf_tgt_phi\",\"ckf_tgt_theta\",\"ckf_tgt_qop\","
              "\"ckf_tgt_s0\",\"ckf_tgt_s1\",\"ckf_tgt_s2\",\"ckf_tgt_s3\",\"ckf_tgt_s4\","
              "\"ref_tgt_loc0\",\"ref_tgt_loc1\",\"ref_tgt_phi\",\"ref_tgt_theta\",\"ref_tgt_qop\","
              "\"ref_tgt_s0\",\"ref_tgt_s1\",\"ref_tgt_s2\",\"ref_tgt_s3\",\"ref_tgt_s4\","
              "\"ckf_w0_geo\",\"ckf_w0_loc0\",\"ckf_w0_loc1\",\"ckf_w0_phi\",\"ckf_w0_theta\",\"ckf_w0_qop\","
              "\"ckf_w0_s0\",\"ckf_w0_s1\",\"ckf_w0_s2\",\"ckf_w0_s3\",\"ckf_w0_s4\","
              "\"ref_w0_geo\",\"ref_w0_loc0\",\"ref_w0_loc1\",\"ref_w0_phi\",\"ref_w0_theta\",\"ref_w0_qop\","
              "\"ref_w0_s0\",\"ref_w0_s1\",\"ref_w0_s2\",\"ref_w0_s3\",\"ref_w0_s4\","
              "\"hit_w0_valid\",\"hit_w0_loc0\",\"hit_w0_loc1\",\"hit_w0_phi\",\"hit_w0_theta\","
              "\"seed_qop\",\"s1_last_flt_qop\",\"s2_first_prt_qop\",\"s2_first_flt_qop\",\"s2_second_prt_qop\","
              "\"s2_first_smt_qop\",\"s3_first_flt_qop\",\"s3_last_smt_qop\",\"s2_first_pull\"]\n"
           << "}\n";
  }

  bool good() const {
    return evOut.good() && sigOut.good() && seedOut.good() && trkOut.good();
  }

  static double pdgCharge(int pdg) {
    const int a = std::abs(pdg);
    if (a == 11 || a == 13 || a == 15) return pdg > 0 ? -1.0 : 1.0;
    if (a == 211 || a == 321 || a == 2212) return pdg > 0 ? 1.0 : -1.0;
    return 0.0;
  }

  static bool isChargedSignalPdg(int pdg) {
    const int a = std::abs(pdg);
    return a == 11 || a == 13 || a == 211 || a == 321 || a == 2212;
  }

  struct WaferHitTruth {
    Acts::Vector3 posF = Acts::Vector3::Zero();
    Acts::Vector3 dirF = Acts::Vector3::Zero();
    int count = 0;
  };

  struct TruthSummary {
    int primaryBc = 0;
    int pdg = 0;
    int status = 0;
    bool isPrimarySignal = false;
    double charge = 0.0;
    Acts::Vector3 vtxF = Acts::Vector3::Zero();
    Acts::Vector3 momF_GeV = Acts::Vector3::Zero();
    double pProd_GeV = std::numeric_limits<double>::quiet_NaN();
    double pSt1_GeV = std::numeric_limits<double>::quiet_NaN();
    double pPreS2_GeV = std::numeric_limits<double>::quiet_NaN();
    double pMinTrk_GeV = std::numeric_limits<double>::quiet_NaN();
    std::array<int, 4> stHitCounts{0, 0, 0, 0};
    std::array<Acts::Vector3, 4> stPosSum{
        Acts::Vector3::Zero(), Acts::Vector3::Zero(), Acts::Vector3::Zero(), Acts::Vector3::Zero()};
    bool fiducial = false;
    double rVtx = 0.0;
    double rMaxS123 = 0.0;
    std::map<Identifier, WaferHitTruth> waferHits;
  };

  struct RefitRecord {
    bool ok = false;
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    double ndf = std::numeric_limits<double>::quiet_NaN();
    double rmax = std::numeric_limits<double>::quiet_NaN();
    double pTrk_GeV = std::numeric_limits<double>::quiet_NaN();
    double qTrk = std::numeric_limits<double>::quiet_NaN();
    Acts::BoundVector tgtBound = Acts::BoundVector::Constant(std::numeric_limits<double>::quiet_NaN());
    Acts::BoundVector tgtSigma = Acts::BoundVector::Constant(std::numeric_limits<double>::quiet_NaN());
    uint64_t w0Geo = 0;
    Acts::BoundVector w0Bound = Acts::BoundVector::Constant(std::numeric_limits<double>::quiet_NaN());
    Acts::BoundVector w0Sigma = Acts::BoundVector::Constant(std::numeric_limits<double>::quiet_NaN());
  };

  template <size_t N>
  static void writeRow(std::ofstream& s, const std::array<double, N>& row) {
    s.write(reinterpret_cast<const char*>(row.data()), sizeof(double) * N);
  }

  void countSegmentsByStation(const TrackCollection* col, std::array<int, 4>& counts) const {
    counts = {0, 0, 0, 0};
    if (!col || !idHelper) return;
    for (const Trk::Track* trk : *col) {
      if (!trk || !trk->trackStateOnSurfaces()) continue;
      for (const Trk::TrackStateOnSurface* tsos : *(trk->trackStateOnSurfaces())) {
        const auto* cot = dynamic_cast<const Tracker::FaserSCT_ClusterOnTrack*>(tsos->measurementOnTrack());
        if (cot) {
          const int st = idHelper->station(cot->identify());
          if (st >= 0 && st < 4) ++counts[st];
          break;
        }
      }
    }
  }

  void recordEvent(
      uint64_t sourceEntry,
      uint64_t runNo,
      uint64_t eventNo,
      const Acts::GeometryContext& gctx,
      const McEventCollection& mcEvents,
      const FaserSiHitCollection& siHits,
      const TrackerSimDataCollection& simData,
      const std::map<Identifier, Acts::GeometryIdentifier>& idMap,
      const CircleFitTrackSeedTool& seedTool,
      double seedOriginZ,
      bool backward,
      const TrackCollection* rawSegs,
      const TrackCollection* gbSegs,
      const FaserActsTrackContainer& selectedTracks,
      const std::vector<unsigned int>& selectedSeedIndices,
      const std::vector<RefitRecord>& refitRecords) {
    const double qnan = std::numeric_limits<double>::quiet_NaN();
    Acts::RotationMatrix3 rotF2A = Acts::RotationMatrix3::Identity();
    rotF2A.col(0) = Acts::Vector3(0, 0, -1);
    rotF2A.col(1) = Acts::Vector3(0, 1, 0);
    rotF2A.col(2) = Acts::Vector3(1, 0, 0);

    Acts::RotationMatrix3 seedPlaneRot = Acts::RotationMatrix3::Identity();
    seedPlaneRot.col(0) = Acts::Vector3(0, 1, 0);
    seedPlaneRot.col(1) = Acts::Vector3(0, 0, -1);
    seedPlaneRot.col(2) = Acts::Vector3(-1, 0, 0);

    // 1. Build primary + Geant4 neutral-decay / spectrometer-hit truth table
    std::set<int> siHitBaseBcs;
    for (const auto& hit : siHits) {
      if (hit.getStation() >= 1 && hit.getStation() <= 3 && hit.trackNumber() > 0) {
        siHitBaseBcs.insert(hit.trackNumber() % 1000000);
      }
    }

    bool hasGenChargedPrimary = false;
    bool hasNeutralPrimaryDecay = false;
    for (const auto* mc : mcEvents) {
      if (!mc) continue;
      for (const auto& p : mc->particles()) {
        if (!p) continue;
        const int bc = HepMC::barcode(p);
        if (bc <= 0 || !isChargedSignalPdg(p->pdg_id())) continue;
        if (bc < 200000 && p->status() == 1) {
          hasGenChargedPrimary = true;
        } else if (bc >= 200000 && bc < 1000000 && p->production_vertex()) {
          for (const auto& inP : p->production_vertex()->particles_in()) {
            if (inP && HepMC::barcode(inP) > 0 && (HepMC::barcode(inP) % 1000000) < 200000 &&
                pdgCharge(inP->pdg_id()) == 0.0) {
              hasNeutralPrimaryDecay = true;
              break;
            }
          }
        }
      }
    }

    std::unordered_map<int, TruthSummary> primaries;
    std::unordered_map<int, int> anyBcToPrimary;
    struct ContNode {
      double z_mm;
      double p_GeV;
    };
    std::unordered_map<int, std::vector<ContNode>> contChains;

    for (const auto* mc : mcEvents) {
      if (!mc) continue;
      const double pscale = (mc->momentum_unit() == HepMC3::Units::GEV) ? 1.0 : 1e-3;
      const double xscale = (mc->length_unit() == HepMC3::Units::CM) ? 10.0 : 1.0;
      for (const auto& p : mc->particles()) {
        if (!p) continue;
        const int bc = HepMC::barcode(p);
        const int pdg = p->pdg_id();
        if (bc <= 0 || !isChargedSignalPdg(pdg)) continue;
        const int baseBc = bc % 1000000;
        const Acts::Vector3 momGeV(p->momentum().px() * pscale,
                                   p->momentum().py() * pscale,
                                   p->momentum().pz() * pscale);
        const double pMag = momGeV.norm();
        Acts::Vector3 vtxMm = Acts::Vector3::Constant(qnan);
        bool fromNeutralPrimary = false;
        if (p->production_vertex()) {
          const auto v = p->production_vertex()->position();
          vtxMm = Acts::Vector3(v.x() * xscale, v.y() * xscale, v.z() * xscale);
          for (const auto& inP : p->production_vertex()->particles_in()) {
            if (inP && HepMC::barcode(inP) > 0 && (HepMC::barcode(inP) % 1000000) < 200000 &&
                pdgCharge(inP->pdg_id()) == 0.0) {
              fromNeutralPrimary = true;
              break;
            }
          }
        }
        if (bc < 1000000) {
          const bool keepTruth = (bc < 200000) || fromNeutralPrimary || (siHitBaseBcs.count(bc) > 0);
          if (!keepTruth) continue;
          bool isSig = false;
          if (hasGenChargedPrimary) {
            isSig = (bc < 200000 && p->status() == 1);
          } else if (hasNeutralPrimaryDecay) {
            isSig = fromNeutralPrimary && std::isfinite(vtxMm.z()) && (vtxMm.z() < 50.0);
          } else {
            isSig = std::isfinite(vtxMm.z()) && (vtxMm.z() < 0.0) && (pMag >= 5.0) && (siHitBaseBcs.count(bc) > 0);
          }
          TruthSummary& ts = primaries[bc];
          ts.primaryBc = bc;
          ts.pdg = pdg;
          ts.status = p->status();
          ts.isPrimarySignal = isSig;
          ts.charge = pdgCharge(pdg);
          ts.vtxF = vtxMm;
          ts.momF_GeV = momGeV;
          ts.pProd_GeV = pMag;
          anyBcToPrimary[bc] = bc;
          if (std::isfinite(vtxMm.z())) {
            contChains[bc].push_back({vtxMm.z(), pMag});
          }
          // Walk explicit HepMC end_vertex same-PDG continuation daughters if present
          auto cur = p;
          std::set<int> visited{bc};
          while (cur && cur->end_vertex()) {
            const auto ev = cur->end_vertex();
            const double ez = ev->position().z() * xscale;
            decltype(cur) nextP = nullptr;
            double bestP = -1.0;
            for (const auto& outP : ev->particles_out()) {
              if (outP && outP->pdg_id() == pdg) {
                const double op = Acts::Vector3(outP->momentum().px(), outP->momentum().py(), outP->momentum().pz()).norm() * pscale;
                if (op > bestP) {
                  bestP = op;
                  nextP = outP;
                }
              }
            }
            if (!nextP) break;
            const int nbc = HepMC::barcode(nextP);
            if (visited.count(nbc)) break;
            visited.insert(nbc);
            anyBcToPrimary[nbc] = bc;
            contChains[bc].push_back({ez, bestP});
            cur = nextP;
          }
        } else if (baseBc > 0) {
          // Athena Geant4 barcode convention k*1000000 + base_barcode
          auto it = primaries.find(baseBc);
          if (it != primaries.end() && it->second.pdg == pdg) {
            anyBcToPrimary[bc] = baseBc;
            if (std::isfinite(vtxMm.z())) {
              contChains[baseBc].push_back({vtxMm.z(), pMag});
            }
          }
        }
      }
    }

    // Also index any non-signal secondary particles referenced by hits/SDO
    std::unordered_map<int, std::pair<int, double>> allParticlesPdgP;
    for (const auto* mc : mcEvents) {
      if (!mc) continue;
      const double pscale = (mc->momentum_unit() == HepMC3::Units::GEV) ? 1.0 : 1e-3;
      for (const auto& p : mc->particles()) {
        if (!p) continue;
        const int bc = HepMC::barcode(p);
        const double pMag = Acts::Vector3(p->momentum().px(), p->momentum().py(), p->momentum().pz()).norm() * pscale;
        allParticlesPdgP[bc] = {p->pdg_id(), pMag};
      }
    }

    // 2. Accumulate SiHits per primary signal particle and per wafer
    std::map<uint64_t, const TrackerDD::SiDetectorElement*> geoToElem;
    if (detMgr) {
      for (const auto* elem : *detMgr->getDetectorElementCollection()) {
        if (!elem) continue;
        auto it = idMap.find(elem->identify());
        if (it != idMap.end()) {
          geoToElem[it->second.value()] = elem;
        }
      }
    }

    for (const auto& hit : siHits) {
      const int rawBc = hit.trackNumber();
      auto mapIt = anyBcToPrimary.find(rawBc);
      const int pbc = (mapIt != anyBcToPrimary.end()) ? mapIt->second : rawBc;
      auto pit = primaries.find(pbc);
      if (pit == primaries.end()) continue;
      const int st = hit.getStation();
      if (st < 0 || st >= 4 || !detMgr || !idHelper) continue;
      const Identifier waferId = idHelper->wafer_id(
          st, hit.getPlane(), hit.getRow(), hit.getModule(), hit.getSensor());
      const auto* elem = detMgr->getDetectorElement(waferId);
      if (!elem) continue;
      const auto ls = hit.localStartPosition(), le = hit.localEndPosition();
      const Acts::Vector3 pStart(ls.x(), ls.y(), ls.z());
      const Acts::Vector3 pEnd(le.x(), le.y(), le.z());
      const Acts::Vector3 gPosF = elem->transformHit() * (0.5 * (pStart + pEnd));
      const Acts::Vector3 gDirF = (elem->transformHit().rotation() * (pEnd - pStart)).normalized();
      pit->second.stHitCounts[st] += 1;
      pit->second.stPosSum[st] += gPosF;
      auto& wh = pit->second.waferHits[waferId];
      wh.posF += gPosF;
      wh.dirF += gDirF;
      wh.count += 1;
    }

    // Finalize continuation momentum bounds for ALL tracked particles and fiducial flags for primary signals
    size_t nSignal = 0, nFiducial = 0;
    for (auto& kv : primaries) {
      auto& ts = kv.second;
      auto& nodes = contChains[ts.primaryBc];
      std::sort(nodes.begin(), nodes.end(), [](const auto& a, const auto& b) { return a.z_mm < b.z_mm; });
      double pSt1 = ts.pProd_GeV;
      double pPreS2 = ts.pProd_GeV;
      double pMinTrk = ts.pProd_GeV;
      for (const auto& nd : nodes) {
        if (nd.z_mm <= 20.0) pSt1 = std::min(pSt1, nd.p_GeV);
        if (nd.z_mm <= 1200.0) pPreS2 = std::min(pPreS2, nd.p_GeV);
        if (nd.z_mm <= 2500.0) pMinTrk = std::min(pMinTrk, nd.p_GeV);
      }
      ts.pSt1_GeV = pSt1;
      ts.pPreS2_GeV = pPreS2;
      ts.pMinTrk_GeV = pMinTrk;

      ts.rVtx = std::isfinite(ts.vtxF.x()) ? std::hypot(ts.vtxF.x(), ts.vtxF.y()) : 9999.0;
      bool allSt = (ts.stHitCounts[1] > 0 && ts.stHitCounts[2] > 0 && ts.stHitCounts[3] > 0);
      double rMax = 0.0;
      for (int s = 1; s <= 3; ++s) {
        if (ts.stHitCounts[s] > 0) {
          const Acts::Vector3 avg = ts.stPosSum[s] / static_cast<double>(ts.stHitCounts[s]);
          rMax = std::max(rMax, std::hypot(avg.x(), avg.y()));
        } else {
          rMax = 9999.0;
        }
      }
      ts.rMaxS123 = rMax;
      const int hits123 = ts.stHitCounts[1] + ts.stHitCounts[2] + ts.stHitCounts[3];
      ts.fiducial = ts.isPrimarySignal && (ts.rVtx < 100.0 && allSt && rMax < 100.0 && hits123 >= 18);
      if (ts.isPrimarySignal) {
        ++nSignal;
        if (ts.fiducial) ++nFiducial;
      }
    }

    auto resolvePrimaryBc = [&](int rawBc) -> int {
      auto it = anyBcToPrimary.find(rawBc);
      if (it != anyBcToPrimary.end()) return it->second;
      const int base = rawBc % 1000000;
      if (rawBc > 0 && primaries.count(base)) return base;
      return rawBc;
    };

    auto matchClusters = [&](const std::vector<const Tracker::FaserSCT_Cluster*>& clusters,
                             int& majBc, double& purity) {
      majBc = 0;
      purity = 0.0;
      if (clusters.empty()) return;
      std::vector<ParticleHitCount> hc;
      identifyContributingParticles(simData, clusters, hc);
      std::unordered_map<int, size_t> mergedCounts;
      for (const auto& c : hc) {
        mergedCounts[resolvePrimaryBc(c.particleId)] += c.hitCount;
      }
      size_t bestHits = 0;
      for (const auto& kv : mergedCounts) {
        if (kv.second > bestHits || (kv.second == bestHits && (majBc == 0 || kv.first < majBc))) {
          bestHits = kv.second;
          majBc = kv.first;
        }
      }
      purity = static_cast<double>(bestHits) / static_cast<double>(clusters.size());
    };

    auto waferTruthOnSurface = [&](const TruthSummary& ts, const Acts::Surface& surf,
                                   Identifier waferId, double& loc0, double& loc1,
                                   double& phiA, double& thetaA) -> bool {
      auto it = ts.waferHits.find(waferId);
      if (it == ts.waferHits.end() || it->second.count <= 0) return false;
      const Acts::Vector3 posF = it->second.posF / static_cast<double>(it->second.count);
      const Acts::Vector3 dirF = it->second.dirF.normalized();
      const Acts::Vector3 posA = rotF2A * posF;
      const Acts::Vector3 dirA = (rotF2A * dirF).normalized();
      auto lp = surf.globalToLocal(gctx, posA, dirA);
      if (lp.ok()) {
        loc0 = lp.value()[0];
        loc1 = lp.value()[1];
      } else {
        const Acts::Vector3 rel = posA - surf.center(gctx);
        loc0 = rel.dot(surf.transform(gctx).rotation().col(0));
        loc1 = rel.dot(surf.transform(gctx).rotation().col(1));
      }
      phiA = Acts::VectorHelpers::phi(dirA);
      thetaA = Acts::VectorHelpers::theta(dirA);
      return true;
    };

    // Track per-signal matching stats
    struct SigMatchAcc {
      int nSeed = 0;
      double bestSeedPurity = 0.0;
      double bestSeedP = std::numeric_limits<double>::quiet_NaN();
      int nCkf = 0;
      size_t bestMeas = 0;
      double bestChi2 = 1e30;
      double bestCkfPurity = 0.0;
      double bestCkfP = std::numeric_limits<double>::quiet_NaN();
      double bestCkfQ = std::numeric_limits<double>::quiet_NaN();
      double bestCkfChi2 = std::numeric_limits<double>::quiet_NaN();
      double bestCkfNdf = std::numeric_limits<double>::quiet_NaN();
      double bestRefitOk = 0.0;
      double bestRefitP = std::numeric_limits<double>::quiet_NaN();
    };
    std::unordered_map<int, SigMatchAcc> sigMatches;

    // 3. Record selected seeds
    const auto& seeds = *seedTool.seeds();
    const auto& initParams = *seedTool.initialTrackParameters();
    size_t nSeedsNan = 0;
    ClusterSet unionPrev(seeds.empty() ? 0 : seeds.front().clusterSet.size());
    std::vector<ClusterSet> prevS23Sets;

    auto buildS23Mask = [&](const CircleFitTrackSeedTool::Seed& s) {
      ClusterSet mask(s.clusterSet.size());
      for (const auto& seg : s.segs) {
        if (seg.station == 2 || seg.station == 3) {
          mask |= seg.clusterSet;
        }
      }
      return mask;
    };

    for (size_t iseed = 0; iseed < seeds.size(); ++iseed) {
      const auto& s = seeds[iseed];
      if (!std::isfinite(s.chi2)) ++nSeedsNan;
      int maxSharedPrev = 0;
      for (size_t j = 0; j < iseed; ++j) {
        if (s.clusterSet.size() == seeds[j].clusterSet.size()) {
          maxSharedPrev = std::max(maxSharedPrev, static_cast<int>((s.clusterSet & seeds[j].clusterSet).count()));
        }
      }
      int sharedUnionPrev = (unionPrev.size() == s.clusterSet.size())
                                ? static_cast<int>((s.clusterSet & unionPrev).count())
                                : 0;
      ClusterSet s23 = buildS23Mask(s);
      int maxSharedS23Prev = 0;
      for (const auto& ps23 : prevS23Sets) {
        if (ps23.size() == s23.size()) {
          maxSharedS23Prev = std::max(maxSharedS23Prev, static_cast<int>((s23 & ps23).count()));
        }
      }
      if (unionPrev.size() == s.clusterSet.size()) unionPrev |= s.clusterSet;
      prevS23Sets.push_back(std::move(s23));

      double minRad = std::numeric_limits<double>::infinity();
      for (const auto& pos : s.fakePositions) {
        const double dz = pos.z() - s.cx;
        minRad = std::min(minRad, s.r * s.r - dz * dz);
      }

      // Segment barcode composition
      std::set<int> segPrimaryBcs;
      bool hasPlusPrimary = false, hasMinusPrimary = false;
      for (const auto& seg : s.segs) {
        int sbc = 0;
        double spur = 0.0;
        matchClusters(seg.clusters, sbc, spur);
        if (sbc > 0) segPrimaryBcs.insert(sbc);
        auto pit = primaries.find(sbc);
        if (pit != primaries.end() && pit->second.isPrimarySignal) {
          if (pit->second.charge > 0) hasPlusPrimary = true;
          if (pit->second.charge < 0) hasMinusPrimary = true;
        }
      }
      const double allSameBc = (s.segs.size() >= 2 && segPrimaryBcs.size() == 1) ? 1.0 : 0.0;
      const double oppPrimaryMix = (hasPlusPrimary && hasMinusPrimary) ? 1.0 : 0.0;

      int majBc = 0;
      double purity = 0.0;
      matchClusters(s.clusters, majBc, purity);

      double tPdg = qnan, tQ = qnan, pProd = qnan, pSt1 = qnan, pPreS2 = qnan;
      double seedTruthValid = 0.0, tLoc0 = qnan, tLoc1 = qnan, tPhiA = qnan, tThetaA = qnan;
      auto pit = primaries.find(majBc);
      if (pit != primaries.end()) {
        const auto& ts = pit->second;
        tPdg = ts.pdg;
        tQ = ts.charge;
        pProd = ts.pProd_GeV;
        pSt1 = ts.pSt1_GeV;
        pPreS2 = ts.pPreS2_GeV;
        // Local linear truth projection from first Station 1 SiHit to z_F = seedOriginZ (field-free ~10 mm air)
        Acts::Vector3 hitPosF = Acts::Vector3::Zero(), hitDirF = Acts::Vector3::Zero();
        int nS1WaferHits = 0;
        double minHitZ = 1e30;
        for (const auto& wh : ts.waferHits) {
          if (idHelper->station(wh.first) == 1 && wh.second.count > 0) {
            const Acts::Vector3 wp = wh.second.posF / static_cast<double>(wh.second.count);
            if (wp.z() < minHitZ) {
              minHitZ = wp.z();
              hitPosF = wp;
              hitDirF = wh.second.dirF.normalized();
              nS1WaferHits = wh.second.count;
            }
          }
        }
        if (nS1WaferHits > 0 && std::abs(hitDirF.z()) > 1e-6) {
          const Acts::Vector3 extrapF = hitPosF - ((hitPosF.z() - seedOriginZ) / hitDirF.z()) * hitDirF;
          const Acts::Vector3 extrapA = rotF2A * extrapF;
          const Acts::Vector3 dirA = (rotF2A * hitDirF).normalized();
          tLoc0 = extrapA.y();
          tLoc1 = -extrapA.z();
          tPhiA = Acts::VectorHelpers::phi(dirA);
          tThetaA = Acts::VectorHelpers::theta(dirA);
          seedTruthValid = 1.0;
        }
        if (purity >= 0.7) {
          auto& sm = sigMatches[majBc];
          sm.nSeed += 1;
          if (purity > sm.bestSeedPurity) {
            sm.bestSeedPurity = purity;
            sm.bestSeedP = (iseed < initParams.size()) ? initParams[iseed].absoluteMomentum() : s.momentum;
          }
        }
      } else if (allParticlesPdgP.count(majBc)) {
        tPdg = allParticlesPdgP[majBc].first;
        tQ = pdgCharge(static_cast<int>(tPdg));
        pProd = allParticlesPdgP[majBc].second;
      }

      const auto bp = (iseed < initParams.size())
                          ? initParams[iseed].parameters()
                          : s.get_params(seedOriginZ, seedTool.covBound(), backward).parameters();
      const double pSeed = (std::abs(bp[4]) > 0.0) ? (1.0 / std::abs(bp[4])) : s.momentum;
      const double qSeed = (bp[4] >= 0.0) ? 1.0 : -1.0;
      const double redChi2 = s.chi2 / static_cast<double>(
          std::max(s.positions.size() + s.fakePositions.size() - s.constraints, 1UL));

      std::array<double, kSeedCols> row{
          static_cast<double>(sourceEntry),
          static_cast<double>(iseed),
          static_cast<double>(s.stations),
          static_cast<double>(s.clusters.size()),
          s.chi2,
          redChi2,
          s.stripChi2,
          static_cast<double>(s.stripNdf),
          std::isnan(s.chi2) ? 1.0 : 0.0,
          minRad,
          s.r,
          pSeed,
          qSeed,
          bp[0], bp[1], bp[2], bp[3],
          static_cast<double>(maxSharedPrev),
          static_cast<double>(sharedUnionPrev),
          static_cast<double>(maxSharedS23Prev),
          allSameBc,
          oppPrimaryMix,
          static_cast<double>(majBc),
          purity,
          tPdg, tQ, pProd, pSt1, pPreS2,
          seedTruthValid, tLoc0, tLoc1, tPhiA, tThetaA};
      writeRow(seedOut, row);
    }

    // 4. Record selected CKF & post-refit tracks
    size_t itrk = 0;
    for (const auto& trk : selectedTracks) {
      const unsigned int iseed = (itrk < selectedSeedIndices.size()) ? selectedSeedIndices[itrk] : 0U;
      const RefitRecord ref = (itrk < refitRecords.size()) ? refitRecords[itrk] : RefitRecord{};

      // Collect states in forward propagation order (trackStatesReversed is tip->stem = S3->S1)
      struct MeasStateInfo {
        int station = -1;
        uint64_t geoId = 0;
        Identifier waferId{};
        const Acts::Surface* surf = nullptr;
        double prtQop = std::numeric_limits<double>::quiet_NaN();
        double fltQop = std::numeric_limits<double>::quiet_NaN();
        double smtQop = std::numeric_limits<double>::quiet_NaN();
        double pull = std::numeric_limits<double>::quiet_NaN();
        Acts::BoundVector smtBound = Acts::BoundVector::Constant(std::numeric_limits<double>::quiet_NaN());
        Acts::BoundVector smtSigma = Acts::BoundVector::Constant(std::numeric_limits<double>::quiet_NaN());
      };
      std::vector<MeasStateInfo> fwdStates;
      std::vector<const Tracker::FaserSCT_Cluster*> trkClusters;
      double rmaxCkf = 0.0;

      for (const auto& st : trk.trackStatesReversed()) {
        if (!st.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag) || !st.hasUncalibratedSourceLink()) {
          continue;
        }
        const auto sl = st.getUncalibratedSourceLink().template get<IndexSourceLink>();
        const auto* cl = sl.hit();
        trkClusters.push_back(cl);
        if (st.hasSmoothed()) {
          const auto sv = st.smoothed();
          const Acts::Vector3 dirA(std::cos(sv[2]) * std::sin(sv[3]),
                                   std::sin(sv[2]) * std::sin(sv[3]),
                                   std::cos(sv[3]));
          const Acts::Vector3 posA = st.referenceSurface().localToGlobal(gctx, sv.template head<2>(), dirA);
          rmaxCkf = std::max(rmaxCkf, std::hypot(posA.y(), posA.z()));
        }
        MeasStateInfo msi;
        const Identifier cid = cl->identify();
        msi.station = idHelper->station(cid);
        msi.waferId = idHelper->wafer_id(cid);
        msi.geoId = st.referenceSurface().geometryId().value();
        msi.surf = &st.referenceSurface();
        if (st.hasPredicted()) msi.prtQop = st.predicted()[4];
        if (st.hasFiltered()) msi.fltQop = st.filtered()[4];
        if (st.hasSmoothed()) {
          msi.smtQop = st.smoothed()[4];
          msi.smtBound = st.smoothed();
          const auto& scov = st.smoothedCovariance();
          for (int k = 0; k < 6; ++k) {
            msi.smtSigma[k] = (scov(k, k) > 0.0) ? std::sqrt(scov(k, k)) : qnan;
          }
        }
        if (st.hasPredicted() && st.hasCalibrated()) {
          const auto H = st.effectiveProjector();
          const double res = st.effectiveCalibrated()[0] - (H * st.predicted())[0];
          const double var = (H * st.predictedCovariance() * H.transpose())(0, 0) +
                             st.effectiveCalibratedCovariance()(0, 0);
          if (var > 0.0) msi.pull = res / std::sqrt(var);
        }
        fwdStates.push_back(msi);
      }
      std::reverse(fwdStates.begin(), fwdStates.end());

      int trkMajBc = 0;
      double trkPurity = 0.0;
      matchClusters(trkClusters, trkMajBc, trkPurity);

      int seedMajBc = 0;
      double seedPurTmp = 0.0;
      if (iseed < seeds.size()) {
        matchClusters(seeds[iseed].clusters, seedMajBc, seedPurTmp);
      }

      double tPdg = qnan, tQ = qnan, pProd = qnan, pSt1 = qnan, pPreS2 = qnan, pMinTrk = qnan;
      double hitW0Valid = 0.0, hitW0Loc0 = qnan, hitW0Loc1 = qnan, hitW0Phi = qnan, hitW0Theta = qnan;

      auto pit = primaries.find(trkMajBc);
      if (pit != primaries.end()) {
        const auto& ts = pit->second;
        tPdg = ts.pdg;
        tQ = ts.charge;
        pProd = ts.pProd_GeV;
        pSt1 = ts.pSt1_GeV;
        pPreS2 = ts.pPreS2_GeV;
        pMinTrk = ts.pMinTrk_GeV;
        if (!fwdStates.empty() && fwdStates.front().surf) {
          if (waferTruthOnSurface(ts, *fwdStates.front().surf, fwdStates.front().waferId,
                                  hitW0Loc0, hitW0Loc1, hitW0Phi, hitW0Theta)) {
            hitW0Valid = 1.0;
          }
        }
        if (trkPurity >= 0.7) {
          auto& sm = sigMatches[trkMajBc];
          sm.nCkf += 1;
          const size_t nM = trk.nMeasurements();
          const double c2 = trk.chi2();
          if (nM > sm.bestMeas || (nM == sm.bestMeas && c2 < sm.bestChi2)) {
            sm.bestMeas = nM;
            sm.bestChi2 = c2;
            sm.bestCkfPurity = trkPurity;
            const double qopCkf = trk.hasReferenceSurface() ? trk.parameters()[4]
                                  : (!fwdStates.empty() ? fwdStates.front().smtQop : qnan);
            sm.bestCkfP = (std::isfinite(qopCkf) && std::abs(qopCkf) > 0.0) ? (1.0 / std::abs(qopCkf)) : qnan;
            sm.bestCkfQ = std::isfinite(qopCkf) ? (qopCkf >= 0.0 ? 1.0 : -1.0) : qnan;
            sm.bestCkfChi2 = c2;
            sm.bestCkfNdf = static_cast<double>(trk.nDoF());
            sm.bestRefitOk = ref.ok ? 1.0 : 0.0;
            sm.bestRefitP = ref.ok ? ref.pTrk_GeV : sm.bestCkfP;
          }
        }
      } else if (allParticlesPdgP.count(trkMajBc)) {
        tPdg = allParticlesPdgP[trkMajBc].first;
        tQ = pdgCharge(static_cast<int>(tPdg));
        pProd = allParticlesPdgP[trkMajBc].second;
      }

      Acts::BoundVector ckfTgt = Acts::BoundVector::Constant(qnan);
      Acts::BoundVector ckfTgtS = Acts::BoundVector::Constant(qnan);
      if (trk.hasReferenceSurface()) {
        ckfTgt = trk.parameters();
        const auto& cov = trk.covariance();
        for (int k = 0; k < 5; ++k) ckfTgtS[k] = (cov(k, k) > 0.0) ? std::sqrt(cov(k, k)) : qnan;
      }

      uint64_t ckfW0Geo = 0;
      Acts::BoundVector ckfW0 = Acts::BoundVector::Constant(qnan);
      Acts::BoundVector ckfW0S = Acts::BoundVector::Constant(qnan);
      if (!fwdStates.empty()) {
        ckfW0Geo = fwdStates.front().geoId;
        ckfW0 = fwdStates.front().smtBound;
        ckfW0S = fwdStates.front().smtSigma;
      }

      double seedQop = (iseed < initParams.size()) ? initParams[iseed].parameters()[4] : qnan;
      double s1LastFlt = qnan, s2FirstPrt = qnan, s2FirstFlt = qnan, s2SecondPrt = qnan;
      double s2FirstSmt = qnan, s3FirstFlt = qnan, s3LastSmt = qnan, s2FirstPull = qnan;
      int s2Count = 0;
      for (const auto& ms : fwdStates) {
        if (ms.station == 1) {
          s1LastFlt = ms.fltQop;
        } else if (ms.station == 2) {
          if (s2Count == 0) {
            s2FirstPrt = ms.prtQop;
            s2FirstFlt = ms.fltQop;
            s2FirstSmt = ms.smtQop;
            s2FirstPull = ms.pull;
          } else if (s2Count == 1) {
            s2SecondPrt = ms.prtQop;
          }
          ++s2Count;
        } else if (ms.station == 3) {
          if (!std::isfinite(s3FirstFlt)) s3FirstFlt = ms.fltQop;
          s3LastSmt = ms.smtQop;
        }
      }

      std::array<double, kTrackCols> row{
          static_cast<double>(sourceEntry),
          static_cast<double>(itrk),
          static_cast<double>(iseed),
          static_cast<double>(trk.nMeasurements()),
          static_cast<double>(trk.nHoles()),
          static_cast<double>(trk.nOutliers()),
          trk.chi2(),
          static_cast<double>(trk.nDoF()),
          rmaxCkf,
          ref.ok ? 1.0 : 0.0,
          ref.chi2,
          ref.ndf,
          ref.rmax,
          static_cast<double>(seedMajBc),
          static_cast<double>(trkMajBc),
          trkPurity,
          tPdg, tQ, pProd, pSt1, pPreS2, pMinTrk,
          ckfTgt[0], ckfTgt[1], ckfTgt[2], ckfTgt[3], ckfTgt[4],
          ckfTgtS[0], ckfTgtS[1], ckfTgtS[2], ckfTgtS[3], ckfTgtS[4],
          ref.tgtBound[0], ref.tgtBound[1], ref.tgtBound[2], ref.tgtBound[3], ref.tgtBound[4],
          ref.tgtSigma[0], ref.tgtSigma[1], ref.tgtSigma[2], ref.tgtSigma[3], ref.tgtSigma[4],
          static_cast<double>(ckfW0Geo),
          ckfW0[0], ckfW0[1], ckfW0[2], ckfW0[3], ckfW0[4],
          ckfW0S[0], ckfW0S[1], ckfW0S[2], ckfW0S[3], ckfW0S[4],
          static_cast<double>(ref.w0Geo),
          ref.w0Bound[0], ref.w0Bound[1], ref.w0Bound[2], ref.w0Bound[3], ref.w0Bound[4],
          ref.w0Sigma[0], ref.w0Sigma[1], ref.w0Sigma[2], ref.w0Sigma[3], ref.w0Sigma[4],
          hitW0Valid, hitW0Loc0, hitW0Loc1, hitW0Phi, hitW0Theta,
          seedQop, s1LastFlt, s2FirstPrt, s2FirstFlt, s2SecondPrt,
          s2FirstSmt, s3FirstFlt, s3LastSmt, s2FirstPull};
      writeRow(trkOut, row);
      ++itrk;
    }

    // 5. Record primary signal truth denominator rows
    std::vector<int> sortedBcs;
    sortedBcs.reserve(primaries.size());
    for (const auto& kv : primaries) sortedBcs.push_back(kv.first);
    std::sort(sortedBcs.begin(), sortedBcs.end());

    for (int bc : sortedBcs) {
      const auto& ts = primaries[bc];
      if (!ts.isPrimarySignal) continue;
      const auto sm = sigMatches[bc];

      double s1HitValid = 0.0, s1Geo = 0.0, s1L0 = qnan, s1L1 = qnan, s1Phi = qnan, s1Theta = qnan;
      // Find earliest wafer hit in Station 1 for standalone signal truth reference
      double minZ = 1e30;
      for (const auto& wh : ts.waferHits) {
        if (idHelper->station(wh.first) == 1 && wh.second.count > 0) {
          const Acts::Vector3 wp = wh.second.posF / static_cast<double>(wh.second.count);
          if (wp.z() < minZ) {
            minZ = wp.z();
            auto git = idMap.find(wh.first);
            if (git != idMap.end()) s1Geo = static_cast<double>(git->second.value());
            const Acts::Vector3 dirA = (rotF2A * wh.second.dirF.normalized()).normalized();
            const Acts::Vector3 posA = rotF2A * wp;
            s1L0 = posA.y();
            s1L1 = -posA.z();
            s1Phi = Acts::VectorHelpers::phi(dirA);
            s1Theta = Acts::VectorHelpers::theta(dirA);
            s1HitValid = 1.0;
          }
        }
      }

      const double txF = (std::abs(ts.momF_GeV.z()) > 1e-9) ? (ts.momF_GeV.x() / ts.momF_GeV.z()) : qnan;
      const double tyF = (std::abs(ts.momF_GeV.z()) > 1e-9) ? (ts.momF_GeV.y() / ts.momF_GeV.z()) : qnan;
      std::array<double, kSignalCols> srow{
          static_cast<double>(sourceEntry),
          static_cast<double>(bc),
          static_cast<double>(ts.pdg),
          ts.charge,
          static_cast<double>(ts.status),
          ts.fiducial ? 1.0 : 0.0,
          static_cast<double>(ts.stHitCounts[0]),
          static_cast<double>(ts.stHitCounts[1]),
          static_cast<double>(ts.stHitCounts[2]),
          static_cast<double>(ts.stHitCounts[3]),
          ts.rVtx,
          ts.rMaxS123,
          ts.vtxF.x(), ts.vtxF.y(), ts.vtxF.z(),
          txF, tyF,
          ts.pProd_GeV, ts.pSt1_GeV, ts.pPreS2_GeV, ts.pMinTrk_GeV,
          s1HitValid, s1Geo, s1L0, s1L1, s1Phi, s1Theta,
          static_cast<double>(sm.nSeed), sm.bestSeedPurity, sm.bestSeedP,
          static_cast<double>(sm.nCkf), sm.bestCkfPurity, sm.bestCkfP, sm.bestCkfQ, sm.bestCkfChi2, sm.bestCkfNdf,
          sm.bestRefitOk, sm.bestRefitP};
      writeRow(sigOut, srow);
    }

    // 6. Record event summary row
    std::array<int, 4> rawC{0, 0, 0, 0}, gbC{0, 0, 0, 0};
    countSegmentsByStation(rawSegs, rawC);
    countSegmentsByStation(gbSegs, gbC);
    std::array<double, kEventCols> erow{
        static_cast<double>(sourceEntry),
        static_cast<double>(runNo),
        static_cast<double>(eventNo),
        static_cast<double>(nSignal),
        static_cast<double>(nFiducial),
        static_cast<double>(rawC[1]), static_cast<double>(rawC[2]), static_cast<double>(rawC[3]),
        static_cast<double>(gbC[1]), static_cast<double>(gbC[2]), static_cast<double>(gbC[3]),
        static_cast<double>(seedTool.rawSeedCombinationCount()),
        static_cast<double>(seedTool.nanChi2RawSeedCount()),
        static_cast<double>(seeds.size()),
        static_cast<double>(nSeedsNan),
        static_cast<double>(selectedTracks.size())};
    writeRow(evOut, erow);
  }
};

#endif  // FASER_LATEST_FAST_VALIDATION_AUDIT_H
