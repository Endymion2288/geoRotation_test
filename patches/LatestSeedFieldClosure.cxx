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
struct MomentumSeedFieldClosure {
  ToolHandle<IFaserActsExtrapolationTool> tool{
      "FaserActsExtrapolationTool/MomentumSeedFieldExtrapolator"};
  Acts::RotationMatrix3 rotation;
  std::array<double,4> centers{};
  MomentumSeedFieldClosure() {
    if(tool.retrieve().isFailure()) throw std::runtime_error("No native extrapolation tool");
    rotation << 0,0,1, 0,1,0, -1,0,0;
    StoreGateSvc* detector=nullptr;
    Gaudi::svcLocator()->service("DetectorStore",detector).ignore();
    const TrackerDD::SCT_DetectorManager* manager=nullptr;
    const FaserSCT_ID* ids=nullptr;
    if(detector->retrieve(manager,"SCT").isFailure() ||
       detector->retrieve(ids,"FaserSCT_ID").isFailure()) throw std::runtime_error("No geometry");
    std::array<int,4> count{};
    for(const auto* element:*manager->getDetectorElementCollection()) {
      int station=ids->station(element->identify());
      centers[station]+=element->center().z();++count[station];
    }
    for(int i=0;i<4;++i) centers[i]/=count[i];
  }
  std::string execute() {
    StoreGateSvc* store=nullptr;
    Gaudi::svcLocator()->service("StoreGateSvc",store).ignore();
    const McEventCollection* events=nullptr;
    if(store->retrieve(events,"TruthEvent").isFailure()) throw std::runtime_error("No truth");
    const auto& ctx=Gaudi::Hive::currentContext();
    auto gctx=tool->trackingGeometryTool()->getGeometryContext(ctx).context();
    auto rows=nlohmann::json::array();
    for(const auto* event:*events) for(const auto& particle:event->particles()) {
      const int barcode=HepMC::barcode(particle);
      if(barcode<=0 || barcode>=200000 || std::abs(particle->pdg_id())!=11 ||
         !particle->production_vertex()) continue;
      double ps=event->momentum_unit()==HepMC3::Units::GEV ? 1. : .001;
      double xs=event->length_unit()==HepMC3::Units::CM ? 10. : 1.;
      auto v=particle->production_vertex()->position();auto m=particle->momentum();
      Acts::Vector3 position=rotation*Acts::Vector3(v.x()*xs,v.y()*xs,v.z()*xs);
      Acts::Vector3 momentum=rotation*Acts::Vector3(m.px()*ps,m.py()*ps,m.pz()*ps);
      double p=momentum.norm(),charge=particle->pdg_id()>0 ? -1. : 1.;
      if(p<10 || position.x()>=centers[1] || momentum.x()<=0) continue;
      Acts::RotationMatrix3 planeRotation;
      planeRotation.col(0)=Acts::Vector3(0,1,0);
      planeRotation.col(1)=Acts::Vector3(0,0,-1);
      planeRotation.col(2)=Acts::Vector3(-1,0,0);
      auto plane=[&](double x) {
        return Acts::Surface::makeShared<Acts::PlaneSurface>(
          Acts::Transform3(Acts::Translation3(x,0,0)*planeRotation));
      };
      auto startSurface=plane(position.x());
      auto local=startSurface->globalToLocal(gctx,position).value();
      Acts::BoundVector mean=Acts::BoundVector::Zero();
      mean[0]=local[0];mean[1]=local[1];mean[2]=std::atan2(momentum.y(),momentum.x());
      mean[3]=std::acos(momentum.z()/p);mean[4]=charge/p;
      Acts::BoundTrackParameters start(startSurface,mean,std::nullopt,Acts::ParticleHypothesis::muon());
      std::vector<Amg::Vector3D> points;
      nlohmann::json states=nlohmann::json::array();
      bool complete=true;
      for(int station=1;station<=3;++station) {
        auto target=plane(centers[station]);
        auto result=tool->propagate(ctx,start,*target);
        if(!result) {complete=false;break;}
        Acts::Vector3 fpos=rotation.transpose()*result->position(gctx);
        Acts::Vector3 fdir=rotation.transpose()*result->direction();
        points.push_back(fpos);points.push_back(fpos-30*fdir);points.push_back(fpos+30*fdir);
        states.push_back({{"station",station},{"position_F_mm",{fpos.x(),fpos.y(),fpos.z()}},
            {"direction_F",{fdir.x(),fdir.y(),fdir.z()}},{"p_GeV",result->absoluteMomentum()}});
      }
      if(!complete) {rows.push_back({{"barcode",barcode},{"success",false}});continue;}
      const auto circle=CircleFit::circleFit(CircleFit::CircleData(points));
      double seedP=circle.r*.3*.55*.001;
      rows.push_back({{"barcode",barcode},{"success",true},{"p_production_GeV",p},
          {"native_circle_p_GeV",seedP},{"native_circle_residual",seedP/p-1},
          {"circle_radius_mm",circle.r},{"effective_circle_field_T",p/(circle.r*.3*.001)},
          {"states",states}});
    }
    return rows.dump();
  }
};
extern "C" void* latest_seed_field_create() {
  return new MomentumSeedFieldClosure();
}
extern "C" const char* latest_seed_field_execute(void* observer) {
  static thread_local std::string result;
  result=static_cast<MomentumSeedFieldClosure*>(observer)->execute();
  return result.c_str();
}
extern "C" void latest_seed_field_delete(void* observer) {
  delete static_cast<MomentumSeedFieldClosure*>(observer);
}
