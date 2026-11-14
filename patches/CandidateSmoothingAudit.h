#ifndef FASER_CANDIDATE_SMOOTHING_AUDIT_H
#define FASER_CANDIDATE_SMOOTHING_AUDIT_H

#include "RefitDirectionChart.h"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "FaserActsKalmanFilter/IndexSourceLink.h"
#include "TrackerPrepRawData/FaserSCT_Cluster.h"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <ostream>

// Read-only instrumentation around the actual Acts::GainMatrixSmoother::calculate.
// Capture its raw covariance through the existing accessor before ACTS validates
// it. No alternative smoothing formula supplies the saved track states.
namespace CandidateSmoothingAudit {
struct Buffer {
  Acts::BoundVector filtered=Acts::BoundVector::Zero();
  Acts::BoundVector predicted=Acts::BoundVector::Zero();
  Acts::BoundVector smoothed=Acts::BoundVector::Zero();
  Acts::BoundSquareMatrix filteredCovariance=Acts::BoundSquareMatrix::Zero();
  Acts::BoundSquareMatrix predictedCovariance=Acts::BoundSquareMatrix::Zero();
  Acts::BoundSquareMatrix smoothedCovariance=Acts::BoundSquareMatrix::Zero();
  Acts::BoundSquareMatrix jacobian=Acts::BoundSquareMatrix::Zero();
  Acts::BoundSquareMatrix rawCovariance=Acts::BoundSquareMatrix::Zero();
  bool capture=false, captured=false;
  unsigned covarianceAccesses=0;
};

template<class Derived>
void eigen(std::ostream& out,const Eigen::MatrixBase<Derived>& value) {
  out << '[';
  for(Eigen::Index i=0;i<value.rows();++i) {
    if(i) out << ',';
    if(value.cols()>1) out << '[';
    for(Eigen::Index j=0;j<value.cols();++j) {
      if(j) out << ',';
      if(std::isfinite(value(i,j))) out << value(i,j); else out << "null";
    }
    if(value.cols()>1) out << ']';
  }
  out << ']';
}

inline Acts::Result<void> run(bool regular,const Acts::GeometryContext&,
    Acts::VectorMultiTrajectory& trajectory,std::size_t tip,const Acts::Logger& logger,
    std::ostream& out,std::size_t event,int seed,int component,std::size_t candidate) {
  Acts::GainMatrixSmoother native;
  Acts::GainMatrixSmoother::GetParameters filtered,predicted,smoothed;
  Acts::GainMatrixSmoother::GetCovariance cf,cp,cs,jacobian;
  filtered.connect([](const void*,void* p) {
    return Eigen::Map<Acts::BoundVector>(static_cast<Buffer*>(p)->filtered.data());
  });
  predicted.connect([](const void*,void* p) {
    return Eigen::Map<Acts::BoundVector>(static_cast<Buffer*>(p)->predicted.data());
  });
  smoothed.connect([](const void*,void* p) {
    return Eigen::Map<Acts::BoundVector>(static_cast<Buffer*>(p)->smoothed.data());
  });
  cf.connect([](const void*,void* p) {
    return Eigen::Map<Acts::BoundSquareMatrix>(static_cast<Buffer*>(p)->filteredCovariance.data());
  });
  cp.connect([](const void*,void* p) {
    return Eigen::Map<Acts::BoundSquareMatrix>(static_cast<Buffer*>(p)->predictedCovariance.data());
  });
  cs.connect([](const void*,void* p) {
    auto* b=static_cast<Buffer*>(p);
    // In the frozen ACTS 32.0.2 calculate(), destination access 1 supplies
    // the assignment target; access 2 reads that raw result for validation.
    // Capture zero matrices as well. Later logging/assignment accesses do
    // not replace the snapshot.
    if(b->capture && ++b->covarianceAccesses==2) {
      b->rawCovariance=b->smoothedCovariance;
      b->captured=true;
    }
    return Eigen::Map<Acts::BoundSquareMatrix>(b->smoothedCovariance.data());
  });
  jacobian.connect([](const void*,void* p) {
    return Eigen::Map<Acts::BoundSquareMatrix>(static_cast<Buffer*>(p)->jacobian.data());
  });
  auto previous=trajectory.getTrackState(tip);
  previous.smoothed()=previous.filtered();
  previous.smoothedCovariance()=previous.filteredCovariance();
  if(!previous.hasPrevious()) return Acts::Result<void>::success();
  std::size_t operation=0;
  std::error_code error;
  trajectory.applyBackwards(previous.previous(),[&](auto state) {
    assert(state.hasFiltered() && state.hasPredicted());
    assert(previous.hasSmoothed() && previous.hasPredicted() && previous.hasJacobian());
    Buffer current,next;
    if(regular) {
      const auto [fm,jf]=refitDirectionChart(state.filtered());
      const auto [pm,jp]=refitDirectionChart(previous.predicted());
      const auto [sm,js]=refitDirectionChart(previous.smoothed());
      current.filtered=fm;
      current.filteredCovariance=jf*state.filteredCovariance()*jf.transpose();
      next.predicted=pm;
      next.predictedCovariance=jp*previous.predictedCovariance()*jp.transpose();
      next.smoothed=sm;
      next.smoothed[Acts::eBoundPhi]=pm[Acts::eBoundPhi]+
          Acts::detail::difference_periodic(sm[Acts::eBoundPhi],pm[Acts::eBoundPhi],2.*std::acos(-1.));
      next.smoothedCovariance=js*previous.smoothedCovariance()*js.transpose();
      next.jacobian=jp*previous.jacobian()*jf.inverse();
    } else {
      current.filtered=state.filtered();
      current.filteredCovariance=state.filteredCovariance();
      next.predicted=previous.predicted();
      next.predictedCovariance=previous.predictedCovariance();
      next.smoothed=previous.smoothed();
      next.smoothedCovariance=previous.smoothedCovariance();
      next.jacobian=previous.jacobian();
    }
    current.capture=true;
    auto status=native.calculate(&current,&next,filtered,cf,smoothed,predicted,cp,cs,jacobian,logger);
    if(!status.ok()) {error=status.error();return false;}
    if(regular) {
      const auto [mean,back]=refitDirectionChart(current.smoothed,true);
      state.smoothed()=mean;
      state.smoothedCovariance()=back*current.smoothedCovariance*back.transpose();
    } else {
      state.smoothed()=current.smoothed;
      state.smoothedCovariance()=current.smoothedCovariance;
    }
    const Acts::BoundSquareMatrix reconstructedGain=current.filteredCovariance*
        next.jacobian.transpose()*next.predictedCovariance.inverse();
    out << std::setprecision(17) << "{\"event\":" << event << ",\"seed\":" << seed
        << ",\"component\":" << component << ",\"candidate\":" << candidate
        << ",\"operation\":" << operation++ << ",\"regular\":" << (regular?"true":"false")
        << ",\"state_index\":" << state.index() << ",\"next_state_index\":" << previous.index()
        << ",\"geometry_id\":" << state.referenceSurface().geometryId().value()
        << ",\"path_length_mm\":" << state.pathLength()
        << ",\"measurement_flag\":" << (state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)?"true":"false")
        << ",\"raw_covariance_captured\":" << (current.captured?"true":"false")
        << ",\"raw_psd\":" << (Acts::detail::covariance_helper<Acts::BoundSquareMatrix>::isSemiPositive(current.rawCovariance)?"true":"false")
        << ",\"validated_psd\":" << (Acts::detail::covariance_helper<Acts::BoundSquareMatrix>::isSemiPositive(current.smoothedCovariance)?"true":"false");
    if(state.hasUncalibratedSourceLink()) {
      const auto sl=state.getUncalibratedSourceLink().template get<IndexSourceLink>();
      out << ",\"hit\":" << sl.hit()->identify().get_compact();
    } else out << ",\"hit\":null";
    auto field=[&](const char* name,const auto& value) {out << ",\"" << name << "\":";eigen(out,value);};
    field("filtered_mean",current.filtered);field("filtered_C",current.filteredCovariance);
    field("next_predicted_mean",next.predicted);field("next_predicted_C",next.predictedCovariance);
    field("next_smoothed_mean",next.smoothed);field("next_smoothed_C",next.smoothedCovariance);
    field("transport_J",next.jacobian);field("reconstructed_gain",reconstructedGain);
    field("calculated_smoothed_mean",current.smoothed);
    field("raw_smoothed_C",current.rawCovariance);field("validated_smoothed_C",current.smoothedCovariance);
    field("native_smoothed_mean",state.smoothed());field("native_smoothed_C",state.smoothedCovariance());
    out << "}\n";
    previous=state;
    return true;
  });
  return error?Acts::Result<void>::failure(error):Acts::Result<void>::success();
}
} // namespace CandidateSmoothingAudit
#endif
