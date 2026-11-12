#!/usr/bin/env python3
"""Independent circle/curvature closure of the actual recorded seed inputs."""
import argparse
from collections import Counter,defaultdict
import gzip
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def taubin(points):
    """Double-precision transcription of the production CircleFit, with diagnostics."""
    xy=np.asarray(points,dtype=float);origin=xy.mean(axis=0);x,y=(xy-origin).T;z=x*x+y*y
    mxx,myy,mxy,mxz,myz,mzz=[float(np.mean(v)) for v in [x*x,y*y,x*y,x*z,y*z,z*z]]
    mz=mxx+myy;cov=mxx*myy-mxy*mxy;var=mzz-mz*mz
    a3=4*mz;a2=-3*mz*mz-mzz;a1=var*mz+4*cov*mz-mxz*mxz-myz*myz
    a0=mxz*(mxz*myy-myz*mxy)+myz*(myz*mxx-mxz*mxy)-var*cov
    root=0.;value=a0;iterations=0
    for iterations in range(999):
        derivative=a1+root*(2*a2+3*a3*root)
        if not derivative:break
        nxt=root-value/derivative
        if nxt==root or not math.isfinite(nxt):break
        nv=a0+nxt*(a1+nxt*(a2+nxt*a3))
        if abs(nv)>=abs(value):break
        root,value=nxt,nv
    det=root*root-root*mz+cov
    if not math.isfinite(det) or det==0:return dict(valid=False,determinant=det if math.isfinite(det) else None,iterations=iterations)
    cx=(mxz*(myy-root)-myz*mxy)/(2*det);cy=(myz*(mxx-root)-mxz*mxy)/(2*det)
    radius2=cx*cx+cy*cy+mz
    if not math.isfinite(radius2) or radius2<=0:return dict(valid=False,determinant=det,iterations=iterations)
    radius=math.sqrt(radius2)
    condition=float(np.linalg.cond(np.cov(xy.T,bias=True)))
    return dict(valid=True,center=(origin+np.array([cx,cy])).tolist(),radius_mm=radius,
        determinant=det,normalized_determinant=det/max(mz*mz,1e-300),iterations=iterations,
        covariance_condition=condition if math.isfinite(condition) else None)


def quadratic_qop(points,B_T=-.55):
    xy=np.asarray(points);z=xy[:,0]-np.mean(xy[:,0]);scale=np.max(abs(z))
    coefficient=np.linalg.lstsq(np.stack([np.ones(len(z)),z/scale,(z/scale)**2],axis=1),xy[:,1],rcond=None)[0]
    slope=coefficient[1]/scale;curvature=2*coefficient[2]/scale**2/(1+slope*slope)**1.5
    return curvature/(.0003*B_T)


def geometric(points,initial):
    xy=np.asarray(points);origin=xy.mean(axis=0);scale=np.max(np.ptp(xy,axis=0));local=(xy-origin)/scale
    center=(np.asarray(initial['center'])-origin)/scale
    def residual(c):
        distances=np.linalg.norm(local-c,axis=1);return distances-distances.mean()
    fit=least_squares(residual,center,xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=300)
    radius=np.linalg.norm(local-fit.x,axis=1).mean()*scale
    return dict(success=bool(fit.success),center=(origin+fit.x*scale).tolist(),radius_mm=float(radius),
        rms_mm=float(np.sqrt(np.mean(residual(fit.x)**2))*scale))


def ideal_closure():
    results=[]
    for p in [10,50,100,250,1000]:
        for q in [-1,1]:
            radius=p/.000165;z=np.array([50.,1240.,2430.]);cy=-q*radius
            y=cy+q*np.sqrt(radius**2-z**2);slope=-q*z/np.sqrt(radius**2-z**2)
            tangent=np.stack([np.ones(3),slope],axis=1);tangent/=np.linalg.norm(tangent,axis=1)[:,None]
            centers=np.stack([z,y],axis=1)
            fake=np.stack([c+offset*t for c,t in zip(centers,tangent) for offset in [0,-30,30]])
            circle=taubin(centers);fakecircle=taubin(fake)
            results.append(dict(p_GeV=p,charge=q,true_points_p_relative_error=circle['radius_mm']/radius-1,
                fake_points_p_relative_error=fakecircle['radius_mm']/radius-1,
                fake_points_charge=1 if fakecircle['center'][1]<0 else -1,
                qop_quadratic=quadratic_qop(fake),taubin=fakecircle))
    return results


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    rows=[];counts=Counter()
    if args.input:
        for path in sorted(args.input.rglob('stages.jsonl.gz')):
            task=json.loads((path.parent/'task.json').read_text())
            records=[json.loads(line) for line in gzip.open(path,'rt')]
            seeds={r['seed']:r for r in records if r['stage']=='seed'}
            for r in records:
                if r['stage']!='seed_circle' or 'fake_positions_F_mm' not in r:continue
                seed=seeds[r['seed']];counts['seeds']+=1
                row=dict(dsid=task['dsid'],file_index=task['file_index'],source_entry=task['skip'],seed=r['seed'],
                    stations=r['stations'],seed_barcode=seed['matches'][0][0] if seed['matches'] else None,
                    segment_matches=[s['matches'] for s in r['segments']],seed_qop=seed['bound'][4],
                    segment_slopes_y=[s['momentum_F'][1]/s['momentum_F'][2] for s in r['segments']])
                if r['stations']>1:
                    points=np.array(r['fake_positions_F_mm'])[:,[2,1]];fit=taubin(points)
                    row['taubin']=fit
                    if fit['valid'] and r['radius_mm'] is not None:
                        row['production_radius_relative_error']=fit['radius_mm']/r['radius_mm']-1
                        row['geometric']=geometric(points,fit)
                        row['quadratic_qop']=quadratic_qop(points)
                        row['spacepoint_taubin']=taubin(np.array(r['spacepoints_F_mm'])[:,[2,1]]) if len(r['spacepoints_F_mm'])>2 else None
                        row['qop_from_radius']=(1 if fit['center'][1]<0 else -1)/(.000165*fit['radius_mm'])
                    else:counts['invalid_circle']+=1
                else:counts['single_station_default']+=1
                rows.append(row)
    result=dict(counts=dict(counts),ideal_closure=ideal_closure(),seeds=rows)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(counts=dict(counts),ideal_closure=result['ideal_closure']),indent=2))


if __name__=='__main__':main()
