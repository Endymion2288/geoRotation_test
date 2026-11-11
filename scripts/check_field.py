#!/usr/bin/env python3
"""Actual field wrapper and map: distinguish grid-boundary and interior derivatives."""
import argparse
import json
from pathlib import Path


CPP = r'''
#include "MagFieldElements/FaserFieldMap.h"
#include "MagFieldConditions/FaserFieldCacheCondObj.h"
#include "FaserActsGeometry/FASERMagneticFieldWrapper.h"
#include "TFile.h"
#include <vector>
#include <cmath>
struct GeoFieldCheck {
  MagField::FaserFieldMap map;
  FaserFieldCacheCondObj conditions;
  FASERMagneticFieldWrapper field;
  GeoFieldCheck(const char* path) {
    TFile input(path,"READ");
    if(!map.initializeMap(&input)) throw std::runtime_error("Field map initialization failed");
    conditions.initialize(1.,&map);
  }
  std::vector<double> point(double x,double y,double z,double h) {
    const FaserFieldCacheCondObj* ptr=&conditions;
    Acts::MagneticFieldContext context(ptr);
    Acts::SquareMatrix3 r; r << 0,0,1, 0,1,0, -1,0,0;
    Acts::Vector3 f(x,y,z), a=r*f;
    auto cache=field.makeCache(context);
    Acts::SquareMatrix3 gradient;
    const Acts::Vector3 b=field.getFieldGradient(a,gradient,cache).value();
    MagField::FaserFieldCache raw(1.,&map);
    double pos[3]={x,y,z}, bf[3], gf[9]; raw.getField(pos,bf,gf);
    Acts::SquareMatrix3 g;
    g << gf[0],gf[1],gf[2],gf[3],gf[4],gf[5],gf[6],gf[7],gf[8];
    Acts::SquareMatrix3 central, forward, backward;
    for(int j=0;j<3;++j) {
      Acts::Vector3 delta=Acts::Vector3::Zero(); delta[j]=h;
      auto cp=field.makeCache(context), cm=field.makeCache(context);
      Acts::Vector3 bp=field.getField(a+delta,cp).value(), bm=field.getField(a-delta,cm).value();
      central.col(j)=(bp-bm)/(2.*h);
      forward.col(j)=(bp-b)/h; backward.col(j)=(b-bm)/h;
    }
    double sideError2=0;
    for(int j=0;j<3;++j) sideError2+=std::min((gradient.col(j)-forward.col(j)).squaredNorm(),(gradient.col(j)-backward.col(j)).squaredNorm());
    std::vector<double> result={(gradient-r*g*r.transpose()*(1000.*Acts::UnitConstants::T)).norm()/Acts::UnitConstants::T,
      (gradient-central).norm()/std::max(gradient.norm(),1e-20),
      (gradient-central).norm()/Acts::UnitConstants::T,
      std::sqrt(sideError2)/std::max(gradient.norm(),1e-20)};
    for(const auto* m: {&gradient,&central,&forward,&backward})
      for(int i=0;i<3;++i) for(int j=0;j<3;++j) result.push_back((*m)(i,j)/Acts::UnitConstants::T);
    // The native mesh also exposes a double-precision interpolation of the
    // same eight field-map corners, independently of its float cache.
    const auto* zone=map.findBFieldZone(x,y,z);
    double bd[3],gd[9]; zone->getB(pos,bd,gd);
    Acts::SquareMatrix3 dg;
    dg << gd[0],gd[1],gd[2],gd[3],gd[4],gd[5],gd[6],gd[7],gd[8];
    Acts::SquareMatrix3 rotatedDouble=r*dg*r.transpose()*(1000.*Acts::UnitConstants::T);
    auto doubleField=[&](const Acts::Vector3& at) {
      Acts::Vector3 ff=r.transpose()*at;
      double pp[3]={ff.x(),ff.y(),ff.z()},bb[3]; zone->getB(pp,bb);
      return Acts::Vector3(r*Acts::Vector3(bb[0],bb[1],bb[2])*(1000.*Acts::UnitConstants::T));
    };
    Acts::SquareMatrix3 dcentral,dplus,dminus;
    const Acts::Vector3 db=doubleField(a);
    for(int j=0;j<3;++j) {
      Acts::Vector3 delta=Acts::Vector3::Zero(); delta[j]=h;
      Acts::Vector3 dp=doubleField(a+delta),dm=doubleField(a-delta);
      dcentral.col(j)=(dp-dm)/(2*h);
      dplus.col(j)=(dp-db)/h;dminus.col(j)=(db-dm)/h;
    }
    double dside2=0;
    for(int j=0;j<3;++j) dside2+=std::min((rotatedDouble.col(j)-dplus.col(j)).squaredNorm(),(rotatedDouble.col(j)-dminus.col(j)).squaredNorm());
    result.push_back((gradient-rotatedDouble).norm()/std::max(rotatedDouble.norm(),1e-20));
    result.push_back((dcentral-rotatedDouble).norm()/std::max(rotatedDouble.norm(),1e-20));
    result.push_back(std::sqrt(dside2)/std::max(rotatedDouble.norm(),1e-20));
    result.push_back((b-db).norm()/Acts::UnitConstants::T);
    result.push_back((dcentral-rotatedDouble).norm()/Acts::UnitConstants::T);
    result.push_back(rotatedDouble.norm()/Acts::UnitConstants::T);
    return result;
  }
};
'''


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--map",default="/cvmfs/faser.cern.ch/repo/sw/software/22.0/faser/offline/ReleaseData/v20/MagneticFieldMaps/FaserFieldTable_v2.root")
    args=parser.parse_args()
    import ROOT
    for library in ["libMagFieldElements","libMagFieldConditions"]:
        assert ROOT.gSystem.Load(library)>=0,library
    assert ROOT.gInterpreter.Declare(CPP)
    checker=ROOT.GeoFieldCheck(args.map)
    rows=[]
    for offset in [False,True]:
        for z in [-1700.,-1000.,0.,500.,1000.,1500.,2000.,2400.]:
            for x in [-50.,0.,50.]:
                f=[x+.371 if offset else x,17.+.293 if offset else 17.,z+.419 if offset else z]
                for h in [.1,.01,.001]:
                    v=list(checker.point(*f,h))
                    rows.append(dict(interior_offset=offset,position_F_mm=f,step_mm=h,
                        rotation_error_T_per_mm=v[0],central_relative_error=v[1],
                        central_error_T_per_mm=v[2],closest_one_sided_relative_error=v[3],
                        gradient_A_T_per_mm=v[4:13],central_A_T_per_mm=v[13:22],
                        forward_A_T_per_mm=v[22:31],backward_A_T_per_mm=v[31:40]))
                    rows[-1].update(gradient_vs_double_relative_error=v[40],double_central_relative_error=v[41],
                        double_one_sided_relative_error=v[42],cached_vs_double_field_error_T=v[43],
                        double_central_error_T_per_mm=v[44],gradient_norm_T_per_mm=v[45])
                    if v[45]<1e-10:
                        for field in ["central_relative_error","closest_one_sided_relative_error",
                                      "gradient_vs_double_relative_error","double_central_relative_error","double_one_sided_relative_error"]:
                            rows[-1][field]=None
    summary={}
    for name,subset in [("grid_points",[r for r in rows if not r["interior_offset"]]),
                        ("interior_points",[r for r in rows if r["interior_offset"]])]:
        summary[name]={"trials":len(subset),"zero_gradient_trials":sum(r["gradient_norm_T_per_mm"]<1e-10 for r in subset),
            **{f"max_{key}":max((r[key] for r in subset if r[key] is not None),default=None) for key in
            ["rotation_error_T_per_mm","central_relative_error","central_error_T_per_mm","closest_one_sided_relative_error",
             "gradient_vs_double_relative_error","double_central_relative_error","double_one_sided_relative_error","cached_vs_double_field_error_T","double_central_error_T_per_mm"]}}
    args.output.write_text(json.dumps(dict(map=args.map,scale=1.,summary=summary,rows=rows),indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":main()
