// Deterministic, opt-in seed and tracking state observer for latest Calypso.
#ifndef FASER_LATEST_SEED_AND_STATE_AUDIT_H
#define FASER_LATEST_SEED_AND_STATE_AUDIT_H

#include "CircleFitTrackSeedTool.h"
#include "CircleFit.h"
#include "TrackingTruthUtils/TrackClassification.h"
#include "GeneratorObjects/McEventCollection.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Units.h"
#include "TrackerSimEvent/FaserSiHitCollection.h"
#include "TrackerReadoutGeometry/SCT_DetectorManager.h"
#include "TrackerReadoutGeometry/SiDetectorElement.h"
#include "TrackerIdentifier/FaserSCT_ID.h"
#include "ActsInterop/IndexSourceLink.h"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "TrkTrack/Track.h"
#include "TrkTrack/TrackCollection.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>

class LatestSeedAndStateAudit {
 public:
  explicit LatestSeedAndStateAudit(const std::string& path, bool states)
      : out(path), fullStates(states) {
    out << std::setprecision(17);
  }
  std::ofstream out;
  bool fullStates;
  bool selectedStatesOnly = false;
  uint64_t sourceEntry = 0;
  uint64_t run = 0, event = 0, lumi = 0;
  std::string collection;
  const TrackerSimDataCollection* simData = nullptr;
  const FaserSCT_ID* idHelper = nullptr;

  void number(double x) {
    if (std::isfinite(x)) out << x;
    else out << "null";
  }
  template <typename M>
  void matrix(const M& m) {
    out << '[';
    for (int i = 0; i < m.rows(); ++i) {
      if (i) out << ',';
      out << '[';
      for (int j = 0; j < m.cols(); ++j) {
        if (j) out << ',';
        number(m(i, j));
      }
      out << ']';
    }
    out << ']';
  }
  template <typename V>
  void vector(const V& v) {
    out << '[';
    for (int i = 0; i < v.size(); ++i) {
      if (i) out << ',';
      number(v[i]);
    }
    out << ']';
  }
  template <typename T>
  void stdvec(const std::vector<T>& v) {
    out << '[';
    for (size_t i = 0; i < v.size(); ++i) {
      if (i) out << ',';
      out << v[i];
    }
    out << ']';
  }
  void begin(const char* stage, int seed = -1, int candidate = -1, int state = -1) {
    out << "{\"schema\":4,\"source_entry\":" << sourceEntry << ",\"run\":" << run
        << ",\"event\":" << event << ",\"lumi\":" << lumi
        << ",\"collection\":\"" << collection << "\",\"stage\":\"" << stage
        << "\",\"seed\":" << seed << ",\"candidate\":" << candidate
        << ",\"state\":" << state;
  }
  void end() { out << "}\n"; }

  void truth(const McEventCollection& events, const FaserSiHitCollection& hits,
             const TrackerDD::SCT_DetectorManager* manager) {
    std::map<int, std::array<Acts::Vector3, 4>> positions, directions;
    std::map<int, std::array<int, 4>> counts;
    for (const auto& hit : hits) {
      const int barcode = hit.trackNumber(), station = hit.getStation();
      if (station < 0 || station >= 4) continue;
      if (!counts.count(barcode)) {
        counts[barcode] = {0, 0, 0, 0};
        for (int s = 0; s < 4; ++s) {
          positions[barcode][s] = Acts::Vector3::Zero();
          directions[barcode][s] = Acts::Vector3::Zero();
        }
      }
      const auto* element = manager->getDetectorElement(idHelper->wafer_id(
          station, hit.getPlane(), hit.getRow(), hit.getModule(), hit.getSensor()));
      const auto start = hit.localStartPosition(), stop = hit.localEndPosition();
      Acts::Vector3 localStart(start.x(), start.y(), start.z());
      Acts::Vector3 localStop(stop.x(), stop.y(), stop.z());
      Acts::Vector3 pos = element->transformHit() * (0.5 * (localStart + localStop));
      Acts::Vector3 dir = (element->transformHit().rotation() * (localStop - localStart)).normalized();
      positions[barcode][station] += pos;
      directions[barcode][station] += dir;
      ++counts[barcode][station];
    }
    for (const auto* mc : events) {
      const double pscale = mc->momentum_unit() == HepMC3::Units::GEV ? 1000. : 1.;
      const double xscale = mc->length_unit() == HepMC3::Units::CM ? 10. : 1.;
      for (const auto& p : mc->particles()) {
        const int barcode = HepMC::barcode(p);
        if (std::abs(p->pdg_id()) != 11 && std::abs(p->pdg_id()) != 13) continue;
        if (barcode % 1000000 >= 200000 && !counts.count(barcode)) continue;
        begin("truth");
        out << ",\"barcode\":" << barcode << ",\"pdg\":" << p->pdg_id()
            << ",\"status\":" << p->status()
            << ",\"signal\":" << (barcode > 0 && barcode < 200000 ? "true" : "false");
        out << ",\"momentum_F_MeV\":";
        vector(Acts::Vector3(p->momentum().px(), p->momentum().py(), p->momentum().pz()) * pscale);
        if (p->production_vertex()) {
          auto v = p->production_vertex()->position();
          out << ",\"vertex_F_mm\":";
          vector(Acts::Vector3(v.x(), v.y(), v.z()) * xscale);
        }
        out << ",\"hit_counts\":[";
        for (int s = 0; s < 4; ++s) {
          if (s) out << ',';
          out << counts[barcode][s];
        }
        out << ']';
        end();
      }
    }
  }

  void countSegmentsByStation(const TrackCollection* col, std::array<int, 4>& counts) {
    counts = {0, 0, 0, 0};
    if (!col || !idHelper) return;
    for (const Trk::Track* trk : *col) {
      if (!trk || !trk->trackStateOnSurfaces()) continue;
      for (const Trk::TrackStateOnSurface* tsos : *(trk->trackStateOnSurfaces())) {
        const auto* cot = dynamic_cast<const Tracker::FaserSCT_ClusterOnTrack*>(tsos->measurementOnTrack());
        if (cot) {
          int st = idHelper->station(cot->identify());
          if (st >= 0 && st < 4) ++counts[st];
          break;
        }
      }
    }
  }

  void seedForensics(const CircleFitTrackSeedTool& tool, double origin,
                     const TrackCollection* rawSegs, const TrackCollection* gbSegs) {
    std::array<int, 4> rawCounts{0, 0, 0, 0}, gbCounts{0, 0, 0, 0};
    countSegmentsByStation(rawSegs, rawCounts);
    countSegmentsByStation(gbSegs, gbCounts);
    const auto seedsPtr = tool.seeds();
    const size_t nSelected = seedsPtr ? seedsPtr->size() : 0;
    begin("seed_event_summary");
    out << ",\"raw_segments_s0123\":[" << rawCounts[0] << ',' << rawCounts[1] << ',' << rawCounts[2] << ',' << rawCounts[3] << ']'
        << ",\"gb_segments_s0123\":[" << gbCounts[0] << ',' << gbCounts[1] << ',' << gbCounts[2] << ',' << gbCounts[3] << ']'
        << ",\"raw_seed_combinations\":" << tool.rawSeedCombinationCount()
        << ",\"raw_nan_chi2_seeds\":" << tool.nanChi2RawSeedCount()
        << ",\"selected_seeds\":" << nSelected;
    end();
    if (!seedsPtr) return;

    Acts::RotationMatrix3 rot = Acts::RotationMatrix3::Identity();
    rot.col(0) = Acts::Vector3(0, 0, -1);
    rot.col(1) = Acts::Vector3(0, 1, 0);
    rot.col(2) = Acts::Vector3(1, 0, 0);
    const double bField = tool.seedBFieldTesla();

    const auto& seeds = *seedsPtr;
    for (size_t iseed = 0; iseed < seeds.size(); ++iseed) {
      const auto& s = seeds[iseed];
      size_t maxSharedPrev = 0;
      for (size_t j = 0; j < iseed; ++j) {
        if (s.clusterSet.size() == seeds[j].clusterSet.size()) {
          maxSharedPrev = std::max(maxSharedPrev, (s.clusterSet & seeds[j].clusterSet).count());
        }
      }
      double minRadicand = std::numeric_limits<double>::infinity();
      for (const auto& pos : s.fakePositions) {
        double rad = -s.cx * s.cx + 2.0 * s.cx * pos.z() + s.r * s.r - pos.z() * pos.z();
        minRadicand = std::min(minRadicand, rad);
      }
      double pCenters3 = std::numeric_limits<double>::quiet_NaN();
      if (s.positions.size() >= 3) {
        CircleFit::CircleData cd(s.positions);
        CircleFit::Circle c3 = CircleFit::circleFit(cd);
        if (std::isfinite(c3.r) && c3.r > 0) {
          pCenters3 = c3.r * 0.3 * bField * 1e-3;
        }
      }
      std::vector<int> segStations, segMajHits, segNClusters;
      std::vector<uint64_t> segBarcodes;
      for (const auto& seg : s.segs) {
        segStations.push_back(seg.station);
        segNClusters.push_back(static_cast<int>(seg.clusters.size()));
        uint64_t topBc = 0;
        int topHits = 0;
        if (simData && !seg.clusters.empty()) {
          std::vector<ParticleHitCount> hc;
          identifyContributingParticles(*simData, seg.clusters, hc);
          std::sort(hc.begin(), hc.end(), [](const auto& a, const auto& b) {
            return a.hitCount != b.hitCount ? a.hitCount > b.hitCount : a.particleId < b.particleId;
          });
          if (!hc.empty()) {
            topBc = hc.front().particleId;
            topHits = hc.front().hitCount;
          }
        }
        segBarcodes.push_back(topBc);
        segMajHits.push_back(topHits);
      }

      const Acts::Vector3 dirChordA = (rot * s.direction).normalized();
      const double phiChord = Acts::VectorHelpers::phi(dirChordA);
      const double thetaChord = Acts::VectorHelpers::theta(dirChordA);

      double phiSeg0 = std::numeric_limits<double>::quiet_NaN();
      double thetaSeg0 = std::numeric_limits<double>::quiet_NaN();
      if (!s.segs.empty() && s.segs[0].momentum.norm() > 0) {
        const Acts::Vector3 dirSeg0A = (rot * s.segs[0].momentum).normalized();
        phiSeg0 = Acts::VectorHelpers::phi(dirSeg0A);
        thetaSeg0 = Acts::VectorHelpers::theta(dirSeg0A);
      }

      double phiCircleTan = std::numeric_limits<double>::quiet_NaN();
      double thetaLine = std::numeric_limits<double>::quiet_NaN();
      if (s.segs.size() > 1 && std::isfinite(s.r) && s.r > 0 && !s.positions.empty()) {
        double rad0 = s.r * s.r - (origin - s.cx) * (origin - s.cx);
        if (rad0 > 0) {
          double sqt0 = std::sqrt(rad0);
          double y1 = s.cy - sqt0, y2 = s.cy + sqt0;
          double y0 = std::abs(y1 - s.positions[0].y()) < std::abs(y2 - s.positions[0].y()) ? y1 : y2;
          if (std::abs(y0 - s.cy) > 1e-9) {
            double dydz = -(origin - s.cx) / (y0 - s.cy);
            double dxdz = s.c1;
            Acts::Vector3 dirTanA = (rot * Acts::Vector3(dxdz, dydz, 1.0)).normalized();
            phiCircleTan = Acts::VectorHelpers::phi(dirTanA);
            thetaLine = Acts::VectorHelpers::theta(dirTanA);
          }
        }
      }

      begin("seed_forensics", static_cast<int>(iseed));
      out << ",\"stations\":" << s.stations << ",\"n_clusters\":" << s.clusters.size()
          << ",\"chi2\":";
      number(s.chi2);
      out << ",\"chi2_is_nan\":" << (std::isnan(s.chi2) ? 1 : 0)
          << ",\"min_radicand\":";
      number(minRadicand);
      out << ",\"circle_r_mm\":";
      number(s.r);
      out << ",\"circle_cx_mm\":";
      number(s.cx);
      out << ",\"circle_cy_mm\":";
      number(s.cy);
      out << ",\"p_fake9_GeV\":";
      number(s.momentum);
      out << ",\"p_centers3_GeV\":";
      number(pCenters3);
      out << ",\"charge_seed\":";
      number(s.charge);
      out << ",\"max_shared_prev\":" << maxSharedPrev
          << ",\"seg_stations\":";
      stdvec(segStations);
      out << ",\"seg_barcodes\":";
      stdvec(segBarcodes);
      out << ",\"seg_maj_hits\":";
      stdvec(segMajHits);
      out << ",\"seg_nclusters\":";
      stdvec(segNClusters);
      out << ",\"phi_chord\":";
      number(phiChord);
      out << ",\"theta_chord\":";
      number(thetaChord);
      out << ",\"phi_seg0\":";
      number(phiSeg0);
      out << ",\"theta_seg0\":";
      number(thetaSeg0);
      out << ",\"phi_circle_tan\":";
      number(phiCircleTan);
      out << ",\"theta_line\":";
      number(thetaLine);
      end();
    }
  }

  void hits(const std::vector<const Tracker::FaserSCT_Cluster*>& clusters) {
    out << ",\"hits\":[";
    std::set<int> layers;
    for (size_t i = 0; i < clusters.size(); ++i) {
      if (i) out << ',';
      out << clusters[i]->identify().get_compact();
      layers.insert(3 * idHelper->station(clusters[i]->identify()) +
                    idHelper->layer(clusters[i]->identify()));
    }
    out << "],\"layers\":" << layers.size() << ",\"matches\":[";
    if (simData) {
      std::vector<ParticleHitCount> counts;
      identifyContributingParticles(*simData, clusters, counts);
      std::sort(counts.begin(), counts.end(), [](const auto& a, const auto& b) {
        return a.hitCount != b.hitCount ? a.hitCount > b.hitCount : a.particleId < b.particleId;
      });
      for (size_t i = 0; i < counts.size(); ++i) {
        if (i) out << ',';
        out << '[' << counts[i].particleId << ',' << counts[i].hitCount << ']';
      }
    }
    out << ']';
  }

  void parameters(const char* stage, const Acts::GeometryContext& gctx,
                  const Acts::BoundTrackParameters& p, int seed, int candidate,
                  const std::vector<const Tracker::FaserSCT_Cluster*>& clusters = {}) {
    begin(stage, seed, candidate);
    out << ",\"bound\":";
    vector(p.parameters());
    if (p.covariance()) {
      out << ",\"variance_qop\":";
      number((*p.covariance())(4, 4));
      if (fullStates) {
        out << ",\"cov\":";
        matrix(*p.covariance());
      }
    }
    out << ",\"position_A_mm\":";
    vector(p.position(gctx));
    out << ",\"momentum_A_GeV\":";
    vector(p.momentum());
    if (fullStates) {
      out << ",\"surface\":";
      matrix(p.referenceSurface().transform(gctx).matrix());
    }
    hits(clusters);
    end();
  }

  template <typename Track>
  void track(const char* stage, const Acts::GeometryContext& gctx,
             const Track& trk, int seed, int candidate) {
    std::vector<const Tracker::FaserSCT_Cluster*> clusters;
    double rmax = 0;
    for (const auto& state : trk.trackStatesReversed()) {
      if (state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag) &&
          state.hasUncalibratedSourceLink()) {
        auto cluster = state.getUncalibratedSourceLink().template get<IndexSourceLink>().hit();
        clusters.push_back(cluster);
        auto pos = cluster->globalPosition();
        rmax = std::max(rmax, std::hypot(pos.x(), pos.y()));
      }
    }
    begin(stage, seed, candidate);
    out << ",\"measurements\":" << trk.nMeasurements() << ",\"holes\":" << trk.nHoles()
        << ",\"outliers\":" << trk.nOutliers() << ",\"ndof\":" << trk.nDoF() << ",\"chi2\":";
    number(trk.chi2());
    out << ",\"rmax_mm\":";
    number(rmax);
    if (trk.hasReferenceSurface()) {
      Acts::BoundTrackParameters p(trk.referenceSurface().getSharedPtr(), trk.parameters(),
                                   trk.covariance(), Acts::ParticleHypothesis::muon());
      out << ",\"bound\":";
      vector(trk.parameters());
      out << ",\"variance_qop\":";
      number(trk.covariance()(4, 4));
      if (fullStates) {
        out << ",\"cov\":";
        matrix(trk.covariance());
      }
      out << ",\"position_A_mm\":";
      vector(p.position(gctx));
      out << ",\"momentum_A_GeV\":";
      vector(p.momentum());
    }
    hits(clusters);
    end();
    if (!fullStates || (selectedStatesOnly && std::string(stage) == "found")) return;
    int index = 0;
    for (const auto& state : trk.trackStatesReversed()) {
      begin((std::string(stage) + "_state").c_str(), seed, candidate, index++);
      out << ",\"geometry_id\":" << state.referenceSurface().geometryId().value()
          << ",\"acts_state_index\":" << state.index()
          << ",\"acts_previous_index\":" << state.previous()
          << ",\"measurement_flag\":" << (state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag) ? "true" : "false")
          << ",\"outlier_flag\":" << (state.typeFlags().test(Acts::TrackStateFlag::OutlierFlag) ? "true" : "false")
          << ",\"hole_flag\":" << (state.typeFlags().test(Acts::TrackStateFlag::HoleFlag) ? "true" : "false");
      out << ",\"surface\":";
      matrix(state.referenceSurface().transform(gctx).matrix());
      out << ",\"path_length_mm\":";
      number(state.pathLength());
      if (state.hasPredicted() && state.referenceSurface().surfaceMaterial()) {
        const auto& mean = state.predicted();
        const Acts::Vector3 direction(std::cos(mean[2]) * std::sin(mean[3]),
                                      std::sin(mean[2]) * std::sin(mean[3]),
                                      std::cos(mean[3]));
        const auto position = state.referenceSurface().localToGlobal(gctx, mean.template head<2>(), direction);
        const auto& slab = state.referenceSurface().surfaceMaterial()->materialSlab(position);
        if (const auto* binned = dynamic_cast<const Acts::BinnedSurfaceMaterial*>(
                state.referenceSurface().surfaceMaterial())) {
          const auto& utility = binned->binUtility();
          out << ",\"material_bin_transform\":";
          matrix(utility.transform().matrix());
          out << ",\"material_bin_local\":";
          vector((utility.transform().inverse() * position).eval());
          out << ",\"material_bin_indices\":";
          vector(utility.binTriple(position));
          out << ",\"material_parameters\":";
          vector(slab.material().parameters());
        }
        out << ",\"material_full_X0\":";
        number(slab.thicknessInX0());
        out << ",\"material_full_thickness_mm\":";
        number(slab.thickness());
        out << ",\"material_path_correction\":";
        number(state.referenceSurface().pathCorrection(gctx, position, direction));
      }
      if (state.hasUncalibratedSourceLink()) {
        const auto sl = state.getUncalibratedSourceLink().template get<IndexSourceLink>();
        const Identifier cid = sl.hit()->identify();
        out << ",\"source_index\":" << sl.index() << ",\"hit\":" << cid.get_compact()
            << ",\"wafer_station\":" << idHelper->station(cid)
            << ",\"wafer_layer\":" << idHelper->layer(cid)
            << ",\"wafer_side\":" << idHelper->side(cid);
      }
      if (state.hasPredicted()) {
        out << ",\"predicted\":";
        vector(state.predicted());
        out << ",\"predicted_cov\":";
        matrix(state.predictedCovariance());
      }
      if (state.hasFiltered()) {
        out << ",\"filtered\":";
        vector(state.filtered());
        out << ",\"filtered_cov\":";
        matrix(state.filteredCovariance());
      }
      if (state.hasSmoothed()) {
        out << ",\"smoothed\":";
        vector(state.smoothed());
        out << ",\"smoothed_cov\":";
        matrix(state.smoothedCovariance());
      }
      if (state.hasJacobian()) {
        out << ",\"jacobian\":";
        matrix(state.jacobian());
      }
      if (state.hasCalibrated()) {
        out << ",\"measurement\":";
        vector(state.effectiveCalibrated());
        out << ",\"measurement_cov\":";
        matrix(state.effectiveCalibratedCovariance());
        out << ",\"projector\":";
        matrix(state.effectiveProjector());
      }
      end();
    }
  }

  void trk(const char* stage, const Trk::Track* track, int seed, int candidate, int refit = -1) {
    if (!track) return;
    begin(stage, seed, candidate);
    out << ",\"refit\":" << refit;
    std::vector<const Tracker::FaserSCT_Cluster*> clusters;
    double rmax = 0;
    if (track->trackParameters() && !track->trackParameters()->empty()) {
      const auto* p = track->trackParameters()->front();
      out << ",\"position_A_mm\":";
      vector(p->position());
      out << ",\"momentum_A_GeV\":";
      vector((p->momentum() * 0.001).eval());
      out << ",\"charge\":";
      number(p->charge());
      out << ",\"bound_trk\":";
      vector(p->parameters());
      if (p->covariance()) {
        out << ",\"variance_qop\":";
        number((*p->covariance())(4, 4) * 1e6);
      }
      for (const auto* q : *track->trackParameters()) {
        rmax = std::max(rmax, std::hypot(q->position().y(), q->position().z()));
      }
    }
    if (track->measurementsOnTrack()) {
      for (const auto* m : *track->measurementsOnTrack()) {
        const auto* cot = dynamic_cast<const Tracker::FaserSCT_ClusterOnTrack*>(m);
        if (cot) clusters.push_back(cot->prepRawData());
      }
    }
    out << ",\"measurements\":" << clusters.size() << ",\"rmax_mm\":";
    number(rmax);
    if (track->fitQuality()) {
      out << ",\"chi2\":";
      number(track->fitQuality()->chiSquared());
      out << ",\"ndof\":";
      number(track->fitQuality()->numberDoF());
    }
    hits(clusters);
    end();
  }
};

#endif
