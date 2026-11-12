// T5: actual ACTS EigenStepper + KalmanFitter, known local truth and fixed hits.
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include <fstream>
#include <iomanip>
#include <random>
#include <iostream>

struct ToyHit { const Acts::Surface* surface; double value, sigma; };
struct ToyCalibrator {
  void calibrate(const Acts::GeometryContext&, const Acts::CalibrationContext&,
      const Acts::SourceLink& sl, Acts::VectorMultiTrajectory::TrackStateProxy state) const {
    const auto& hit=sl.get<ToyHit>();
    state.setUncalibratedSourceLink(sl);
    auto m=Acts::makeMeasurement(sl,Acts::ActsVector<1>::Constant(hit.value),
        Acts::ActsSquareMatrix<1>::Constant(hit.sigma*hit.sigma),Acts::eBoundLoc0);
    state.allocateCalibrated(1);state.setCalibrated(m);
  }
  const Acts::Surface* surface(const Acts::SourceLink& sl) const {return sl.get<ToyHit>().surface;}
};

int main(int argc,char** argv) {
  const int replicas=argc>1 ? std::stoi(argv[1]) : 1;
  std::ofstream output(argc>2 ? argv[2] : "toy.csv");
  output << std::setprecision(17) << "p_truth_GeV,charge,replica,mode,success,qop_reco,qop_sigma,p_reco_GeV,chi2,ndof\n";
  Acts::GeometryContext gctx;Acts::MagneticFieldContext mctx;Acts::CalibrationContext cctx;
  auto field=std::make_shared<Acts::ConstantBField>(Acts::Vector3(0,0,-0.55*Acts::UnitConstants::T));
  using Stepper=Acts::EigenStepper<>;
  using Generator=Acts::Propagator<Stepper,Acts::VoidNavigator>;
  Generator generator(Stepper(field),Acts::VoidNavigator{});
  using Propagator=Acts::Propagator<Stepper,Acts::DirectNavigator>;
  using Fitter=Acts::KalmanFitter<Propagator,Acts::VectorMultiTrajectory>;
  auto logger=Acts::getDefaultLogger("ToyKF",Acts::Logging::ERROR);
  Fitter fitter(Propagator(Stepper(field),Acts::DirectNavigator{}),std::move(logger));
  Acts::GainMatrixUpdater updater;Acts::GainMatrixSmoother smoother;ToyCalibrator calibrator;
  Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
  extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&updater);
  extensions.smoother.connect<&Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(&smoother);
  extensions.calibrator.connect<&ToyCalibrator::calibrate>(&calibrator);
  extensions.surfaceAccessor.connect<&ToyCalibrator::surface>(&calibrator);
  auto plane=[](double x,double stereo) {
    Acts::Transform3 t=Acts::Transform3::Identity();
    t.translation()=Acts::Vector3(x,0,0);
    t.linear().col(0)=Acts::Vector3(0,std::cos(stereo),std::sin(stereo));
    t.linear().col(1)=Acts::Vector3(0,std::sin(stereo),-std::cos(stereo));
    t.linear().col(2)=Acts::Vector3(-1,0,0);
    return Acts::Surface::makeShared<Acts::PlaneSurface>(t);
  };
  auto start=plane(0,0);
  std::vector<std::shared_ptr<Acts::PlaneSurface>> owned;
  std::vector<const Acts::Surface*> sequence;
  int index=0;
  for(double station:{0.,1190.,2380.}) for(double layer:{21.415,52.915,84.415}) for(int side:{0,1}) {
    auto s=plane(station+layer+0.885*side,side ? -0.02:0.02);
    s->assignGeometryId(Acts::GeometryIdentifier{}.setVolume(1).setLayer(2*index+2).setSensitive(index+1));
    ++index;sequence.push_back(s.get());owned.push_back(s);
  }
  std::mt19937_64 rng(84001);std::normal_distribution<double> normal(0,1);
  size_t failures=0,total=0;
  for(double p:{10.,50.,100.,250.,1000.}) for(double charge:{-1.,1.}) {
    Acts::BoundVector truth;truth << 0,0,0.002,M_PI/2-0.001,charge/p,0;
    Acts::BoundTrackParameters generated(start,truth,std::nullopt,Acts::ParticleHypothesis::muon());
    std::vector<double> local;
    for(const auto* surface:sequence) {
      Acts::PropagatorOptions<> options(gctx,mctx);options.maxStepSize=10.;
      auto r=generator.propagate(generated,*surface,options);
      if(!r.ok() || !r->endParameters) return 2;
      local.push_back(r->endParameters->parameters()[0]);
    }
    for(int replica=0;replica<replicas;++replica) {
      std::vector<Acts::SourceLink> hits;
      for(size_t i=0;i<sequence.size();++i) hits.emplace_back(ToyHit{sequence[i],local[i]+(replica ? 0.016*normal(rng):0),0.016});
      for(int mode=0;mode<3;++mode) {
        auto mean=truth;
        if(mode<2) {mean[0]+=0.2;mean[1]-=1.;mean[2]+=0.002;mean[3]-=0.001;}
        if(mode==0) mean[4]*=3.;
        Acts::BoundSquareMatrix covariance=Acts::BoundSquareMatrix::Zero();
        covariance.diagonal() << 81.,0.0361,0.1089,0.000841,0.001,1.;
        Acts::BoundTrackParameters initial(start,mean,covariance,Acts::ParticleHypothesis::muon());
        Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> options(gctx,mctx,cctx,extensions,Acts::PropagatorPlainOptions{},start.get());
        options.energyLoss=false;options.multipleScattering=false;
        options.propagatorPlainOptions.maxStepSize=10.;
        options.referenceSurfaceStrategy=Acts::KalmanFitterTargetSurfaceStrategy::first;
        auto tc=std::make_shared<Acts::VectorTrackContainer>();auto sc=std::make_shared<Acts::VectorMultiTrajectory>();
        Acts::TrackContainer<Acts::VectorTrackContainer,Acts::VectorMultiTrajectory,std::shared_ptr> tracks(tc,sc);
        auto result=fitter.fit(hits.begin(),hits.end(),initial,options,sequence,tracks);
        ++total;
        output << p << ',' << charge << ',' << replica << ',' << mode << ',' << result.ok();
        if(result.ok()) {
          const auto& track=result.value();
          double qop=track.parameters()[4],sigma=std::sqrt(track.covariance()(4,4));
          output << ',' << qop << ',' << sigma << ',' << std::abs(1/qop) << ',' << track.chi2() << ',' << track.nDoF();
        } else {++failures;output << ",nan,nan,nan,nan,0";}
        output << '\n';
      }
    }
  }
  std::cout << total << " fits; " << failures << " failures\n";
  return failures ? 1:0;
}
