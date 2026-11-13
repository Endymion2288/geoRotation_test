"""Charge-independent truth matching, denominators and source-event rates."""
import numpy as np

from charge_history import rate
from slope_tables import join_truth


STATUS = ['unreconstructed_or_unmatched','correct','wrong','invalid_curvature',
          'low_purity','quality_failed','not_fiducial']


def rank_candidates(truth_index,eligible,tracks,n_truth):
    """Rank on hit evidence and fit quality, then apply all physics selections."""
    rows=np.flatnonzero(eligible)
    chosen=np.full(n_truth,-1,dtype=np.int64)
    duplicates=np.zeros(n_truth,dtype=np.int32)
    if not len(rows):return chosen,duplicates
    matched_hits=np.rint(tracks['purity'][rows]*tracks['n_hits'][rows]).astype(np.int64)
    order=np.lexsort((tracks['ordinal'][rows],tracks['chi2'][rows],
                      -tracks['n_hits'][rows],-matched_hits,truth_index[rows]))
    rows=rows[order]
    values,first,counts=np.unique(truth_index[rows],return_index=True,return_counts=True)
    chosen[values]=rows[first];duplicates[values]=counts
    return chosen,duplicates


def derive(tables,dsid,radius=95):
    t,h,e=tables['tracks'],tables['truth'],tables['events']
    pdg=-13 if dsid=='100116' else 13;charge=-np.sign(pdg)
    primary=(h['pdg']==pdg)&(h['status']==1)&(h['barcode']>0)&(h['barcode']<200000)
    assert len(np.unique(h['sample_entry'][primary]))==primary.sum(), 'More than one primary muon per source event'
    p=np.linalg.norm(h['production_momentum_F_GeV'],axis=1)
    station=h['truth_station_position_F_mm'][:,1:,:2]
    vertex=h['production_position_F_mm'][:,:2]
    fiducial=primary & (np.linalg.norm(vertex,axis=1)<100) & np.all(np.linalg.norm(station,axis=2)<100,axis=1)
    index,matched=join_truth(tables)
    eligible=np.zeros(len(index),dtype=bool)
    eligible[matched]=primary[index[matched]]
    chosen,duplicates=rank_candidates(index,eligible,t,len(h['sample_entry']))
    present=chosen>=0;safe=np.maximum(chosen,0)
    with np.errstate(divide='ignore',invalid='ignore'):
        quality=(t['station_bits']&14==14)&(t['ndof']>0)&(t['chi2']/t['ndof']<25)& \
                (t['n_layers']>=7)&(t['ndof']+5>=12)&(t['rmax_mm']<radius)
    qop=t['qop_GeVinv']
    valid=np.isfinite(qop)&(qop!=0)&np.isin(t['charge'],[-1,1])
    np.testing.assert_array_equal(np.sign(qop[valid]),t['charge'][valid])
    pure=t['purity']>.7
    status=np.zeros(len(chosen),dtype=np.int8)
    if len(t['source_entry']):
        status[present&~quality[safe]]=5
        status[present&quality[safe]&~pure[safe]]=4
        status[present&quality[safe]&pure[safe]&~valid[safe]]=3
        selected=present&quality[safe]&pure[safe]&valid[safe]
        status[selected]=np.where(t['charge'][safe[selected]]*charge<0,2,1)
    status[~fiducial]=6
    selected=(status==1)|(status==2)
    event_track=np.searchsorted(e['sample_entry'],t['sample_entry'])
    event_truth=np.searchsorted(e['sample_entry'],h['sample_entry'])
    np.testing.assert_array_equal(e['sample_entry'][event_track],t['sample_entry'])
    n_good=np.bincount(event_track,weights=quality,minlength=len(e['sample_entry'])).astype(int)
    one=n_good[event_truth]==1
    direction=h['truth_sensor_direction_F']
    with np.errstate(divide='ignore',invalid='ignore'):
        sx=direction[:,0]/direction[:,2];sy=direction[:,1]/direction[:,2]
    features=np.column_stack([np.log(np.maximum(p,1e-20)),h['truth_station_position_F_mm'][:,1,0],
        h['truth_station_position_F_mm'][:,1,1],sx,sy])
    return dict(primary=primary,fiducial=fiducial,truth_p=p,chosen=chosen,safe=safe,
        duplicates=duplicates,track_matched=eligible,track_quality=quality,track_valid=valid,
        status=status,selected=selected,wrong=status==2,one_good=one,features=features,
        truth_H18=np.sum(h['truth_hit_counts'][:,1:],axis=1)>=18,
        charge=charge,missing_track_parameters=int(e['n_tracks'].sum()-len(t['source_entry'])))


def measurement(tables,d,mask=None):
    h,t=tables['truth'],tables['tracks']
    scope=np.ones(len(d['status']),dtype=bool) if mask is None else mask
    fid=d['fiducial']&scope;selected=d['selected']&scope;wrong=d['wrong']&scope
    def efficiency(numerator,denominator):
        result=rate(numerator,denominator)
        result['numerator']=result.pop('wrong');result['denominator']=result.pop('total')
        return result
    result=dict(primary_truth=int((d['primary']&scope).sum()),fiducial_truth=int(fid.sum()),
        conditional_misid=rate(wrong.sum(),selected.sum()),
        truth_correct_efficiency=efficiency(np.count_nonzero(fid&(d['status']==1)),fid.sum()),
        truth_wrong_efficiency=efficiency(np.count_nonzero(fid&wrong),fid.sum()),
        statuses={STATUS[i]:int(np.count_nonzero(scope&(d['status']==i))) for i in range(len(STATUS))},
        fiducial_truth_with_duplicates=int(np.count_nonzero(fid&(d['duplicates']>1))),
        exactly_one_good=rate(np.count_nonzero(wrong&d['one_good']),np.count_nonzero(selected&d['one_good'])),
        H18=rate(np.count_nonzero(wrong&d['truth_H18']),np.count_nonzero(selected&d['truth_H18'])),
        track_parameters_missing=d['missing_track_parameters'])
    return result


def weighted_rate(wrong,selected,weights):
    selected=selected&(weights>0)
    w=weights[selected];y=wrong[selected].astype(float)
    sw=float(w.sum());sw2=float(w@w)
    if not sw:return dict(value=None,sumw=0.,sumw2=0.,effective_n=0.,standard_error=None,
                           interval68=[None,None],interval95=[None,None])
    value=float(w@y/sw)
    variance=float(np.sum(w*w*(y-value)**2)/(sw*sw))
    error=float(np.sqrt(variance));effective=sw*sw/sw2
    # At a boundary retain a nonzero binomial limit based on effective N.
    if value in [0.,1.]:
        from scipy.stats import beta
        interval=[0.,float(beta.ppf(.975,1,effective))] if value==0 else [float(beta.ppf(.025,effective,1)),1.]
        interval68=[0.,float(beta.ppf(.84,1,effective))] if value==0 else [float(beta.ppf(.16,effective,1)),1.]
    else:
        interval=[max(0.,value-1.96*error),min(1.,value+1.96*error)]
        interval68=[max(0.,value-error),min(1.,value+error)]
    return dict(value=value,sumw=sw,sumw2=sw2,effective_n=effective,standard_error=error,
                interval68=interval68,interval95=interval,weighted_wrong=float(w@y),raw_selected=int(selected.sum()),
                raw_wrong=int(np.count_nonzero(wrong&selected)))


def truth_balance(data,derived):
    """Freeze joint truth bins on discovery files; evaluate held-out balance."""
    labels=['100116','100117']
    source={k:derived[k]['features'] for k in labels}
    eligible={k:derived[k]['fiducial']&np.isfinite(source[k]).all(axis=1) for k in labels}
    train={k:eligible[k]&(data[k]['truth']['file_index']<3) for k in labels}
    pooled=np.concatenate([source[k][train[k]] for k in labels])
    assert len(pooled)
    edges=[np.log(np.array([1e-6,50,100,250,500,1000,1e6]))]
    edges.extend([np.unique(np.r_[-np.inf,np.quantile(pooled[:,j],[.25,.5,.75]),np.inf]) for j in range(1,5)])
    shape=tuple(len(e)-1 for e in edges)
    cell={}
    for k in labels:
        positions=[np.clip(np.searchsorted(e,source[k][:,j],side='right')-1,0,len(e)-2) for j,e in enumerate(edges)]
        cell[k]=np.ravel_multi_index(positions,shape)
    counts={k:np.bincount(cell[k][train[k]],minlength=int(np.prod(shape))) for k in labels}
    support=(counts[labels[0]]>=20)&(counts[labels[1]]>=20)
    assert support.any()
    target=.5*(counts[labels[0]]+counts[labels[1]])
    weights={}
    report=dict(variables=['log truth p','st1 truth x','st1 truth y','sensor truth sx','sensor truth sy'],
        discovery_files=[0,1,2],held_out_files=[3,4],minimum_discovery_per_charge_per_cell=20,
        shared_cells=int(support.sum()),cells=int(np.prod(shape)),
        edges=[[None if not np.isfinite(v) else float(v) for v in e] for e in edges],samples={})
    for k in labels:
        factors=np.zeros(len(support));factors[support]=target[support]/counts[k][support]
        weights[k]=factors[cell[k]]*eligible[k]
        report['samples'][k]=dict(eligible_truth=int(eligible[k].sum()),supported_truth=int((weights[k]>0).sum()),
            supported_held_out_truth=int(np.count_nonzero((weights[k]>0)&(data[k]['truth']['file_index']>=3))),
            sumw=float(weights[k].sum()),sumw2=float(weights[k]@weights[k]),
            effective_n=float(weights[k].sum()**2/(weights[k]@weights[k])),
            largest_weight=float(weights[k].max()))
    report['balance']={}
    for split in ['all','held_out']:
        masks={k:eligible[k] & ((data[k]['truth']['file_index']>=3) if split=='held_out' else True) for k in labels}
        before=[];supported=[];after=[]
        for j in range(5):
            a,b=[source[k][masks[k],j] for k in labels]
            scale=np.sqrt(.5*(np.var(a)+np.var(b)))
            before.append(float((np.mean(a)-np.mean(b))/scale) if scale else 0.)
            common_means=[np.mean(source[k][masks[k]&(weights[k]>0),j]) for k in labels]
            supported.append(float((common_means[0]-common_means[1])/scale) if scale else 0.)
            means=[np.average(source[k][masks[k],j],weights=weights[k][masks[k]]) for k in labels]
            after.append(float((means[0]-means[1])/scale) if scale else 0.)
        report['balance'][split]=dict(standardized_mean_difference_before=before,
                                     standardized_mean_difference_common_support=supported,
                                     standardized_mean_difference_after=after)
    return weights,report
