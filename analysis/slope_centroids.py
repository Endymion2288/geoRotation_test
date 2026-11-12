#!/usr/bin/env python3
"""Paired actual Calypso reconstructions with native/truth strip coordinates."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_tables import histogram,join_truth,slope,write_json
from slope_statistics import event_ratio_interval
from slope_controls import archive,save
from slope_origin import harmonic


def collect(paths):
    parts={k:[] for k in ['events','truth','tracks','segments']}
    offset=0
    for path in paths:
        for name in parts:
            d=archive(path/(name+'.npz'))
            begin=archive(path/'events.npz')['source_entry'].min()
            d['sample_entry']=d['source_entry']-begin+offset
            if name=='tracks':
                d['hit_signatures']=np.array([','.join(str(x) for x in np.sort(d['hit_identifiers'][lo:hi]))
                    for lo,hi in zip(d['hit_offsets'][:-1],d['hit_offsets'][1:])])
                del d['hit_offsets'],d['hit_identifiers']
            parts[name].append(d)
        offset+=len(parts['events'][-1]['source_entry'])
    return {name:{k:np.concatenate([p[k] for p in rows]) for k in rows[0]} for name,rows in parts.items()},offset


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['rdo','gun'])
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((args.root/'manifest.json').read_text())
    cases={'production_muplus':[(Path(p['baseline']),Path(p['truth_position'])) for p in manifest['pairs']]} if args.mode=='rdo' else {
        str(pdg):[(Path(s['output'])/'afterGeo',Path(s['output'])/'afterGeo_truth_position')
                 for s in manifest['sources'] if s['pdg']==pdg] for pdg in [-13,13]}
    result={'mode':args.mode,'definition':'Actual full Calypso reconstruction; only measured local strip coordinate replaced with majority SiHit centroid',
            'period':.08/63.,'cases':{}}
    with PdfPages(args.output/'centroid_control.pdf') as report:
        for name,pairs in cases.items():
            tables={}
            for label,col in [('native',0),('truth_position',1)]:
                tables[label],count=collect([p[col] for p in pairs])
            b,a=tables['native'],tables['truth_position']
            truth=b['truth']
            for key in truth:
                if np.issubdtype(truth[key].dtype,np.floating):assert np.allclose(truth[key],a['truth'][key],rtol=1e-13,atol=1e-13,equal_nan=True),key
                else:assert np.array_equal(truth[key],a['truth'][key]),key
            case=dict(source_events=count,truth_identity_and_values_equal=True,configurations={})
            fig,axes=plt.subplots(2,2,figsize=(12,8))
            chosen={}
            for label,t in tables.items():
                r=t['tracks'];index,matched=join_truth(t)
                signal=matched & (np.abs(truth['pdg'][index])==13)
                long=r['station_bits']&14==14
                selected=long & signal & (r['purity']>=.7)
                efficient=np.bincount(index[selected],minlength=len(truth['barcode']))>0
                population=np.abs(truth['pdg'])==13
                c=dict(native_tracks=len(long),long_tracks=int(long.sum()),matched_pure_long=int(selected.sum()),
                    efficiency=event_ratio_interval(truth['sample_entry'],efficient&population,population,(0,count)),
                    wrong_charge=event_ratio_interval(r['sample_entry'],selected & (r['charge']!=np.sign(-truth['pdg'][index])),selected,(0,count)),
                    momentum_tail=event_ratio_interval(r['sample_entry'],selected & (np.linalg.norm(r['upstream_F'][:,3:],axis=1)<.5*np.linalg.norm(truth['production_momentum_F_GeV'][index],axis=1)),selected,(0,count)),
                    planes={})
                values=np.stack([slope(r['projection'][:,i,3:6]) for i in [0,1]],axis=1)
                values=np.where(r['projection_valid'][:,:2]==1,values,np.nan)
                for i,plane in enumerate(['st1','VetoNu']):
                    h=histogram(values[selected,i],np.linspace(-.01,.01,201));c['planes'][plane]=h
                    axes[0,i].stairs(h['counts'],h['edges'],label=f"{label}; N={h['finite']}")
                    axes[0,i].set(title=plane,xlabel='Physical sy',ylabel='Entries / bin')
                si=np.where(r['projection_valid'][:,2]==1,slope(r['projection'][:,2,3:6]),np.nan)
                residual=si-slope(truth['truth_sensor_direction_F'])[index]
                h=histogram(residual[selected],np.linspace(-.0003,.0003,151));c['sensor_residual']=h
                axes[1,0].stairs(h['counts'],h['edges'],label=label)
                axes[1,0].set(xlabel='Same-sensor sy residual',ylabel='Entries / bin')
                segment=t['segments'];ti,linked=join_truth({'truth':truth,'tracks':segment})
                m=linked & (np.abs(truth['pdg'][ti])==13) & (segment['station']==1)
                v=slope(segment['parameters_F'][:,3:])
                h=histogram(v[m],np.linspace(-.008,.008,321));c['segment1']=h
                c['segment1_harmonic']=harmonic(v[m],segment['sample_entry'][m],.08/63.,count)
                axes[1,1].stairs(h['counts'],h['edges'],label=label)
                axes[1,1].set(xlabel='Raw station-1 segment sy',ylabel='Entries / bin')
                best={}
                for row in np.flatnonzero(selected):
                    key=(int(r['sample_entry'][row]),int(r['barcode'][row]))
                    rank=(-r['purity'][row],-r['n_hits'][row],r['chi2'][row]/max(1,r['ndof'][row]))
                    if key not in best or rank<best[key][0]:best[key]=(rank,row)
                chosen[label]={k:v[1] for k,v in best.items()}
                case['configurations'][label]=c
            for ax in axes.flat:ax.legend(fontsize=8);ax.grid(alpha=.2)
            fig.suptitle(f'Native Calypso centroid control | {name} | {count:,} source events\nSame G4 truth and digitized hits; original cluster widths, errors and noise retained')
            save(fig,args.output,name+'_centroid',report)
            keys=sorted(chosen['native'].keys() & chosen['truth_position'].keys())
            rb,ra=b['tracks'],a['tracks']
            bi=np.array([chosen['native'][k] for k in keys]);ai=np.array([chosen['truth_position'][k] for k in keys])
            same=rb['hit_signatures'][bi]==ra['hit_signatures'][ai]
            momentum=np.linalg.norm(rb['upstream_F'][bi,3:],axis=1)
            common=same & (momentum>50)&(momentum<250)
            fig,axes=plt.subplots(1,2,figsize=(11,4.5))
            case['paired']=dict(common_particles=len(keys),same_cluster_identifiers=int(same.sum()),fixed_native_p_window=int(common.sum()),
                matching='One track per majority truth, ranked by purity, hits and chi2/ndof; no charge or momentum matching; same complete cluster identifier set',planes={})
            for i,plane in enumerate(['st1','VetoNu']):
                vb=np.where(rb['projection_valid'][bi,i]==1,slope(rb['projection'][bi,i,3:6]),np.nan)
                va=np.where(ra['projection_valid'][ai,i]==1,slope(ra['projection'][ai,i,3:6]),np.nan)
                finite=common&np.isfinite(vb)&np.isfinite(va)
                vb,va=vb[finite],va[finite]
                ids=rb['sample_entry'][bi[finite]]
                hb,ha=histogram(vb,np.linspace(-.008,.008,321)),histogram(va,np.linspace(-.008,.008,321))
                axes[i].stairs(hb['counts'],hb['edges'],label='Native coordinates')
                axes[i].stairs(ha['counts'],ha['edges'],label='Truth coordinates')
                axes[i].set(title=plane,xlabel='Physical sy',ylabel='Same-hit tracks / bin');axes[i].legend(fontsize=8);axes[i].grid(alpha=.2)
                zb,za=np.exp(2j*np.pi*vb/(.08/63.)),np.exp(2j*np.pi*va/(.08/63.))
                # Paired event bootstrap retains the native/ideal correlation.
                rng=np.random.default_rng(50561);delta=[]
                for _ in range(500):
                    w=rng.multinomial(count,np.full(count,1./count))[ids]
                    if w.sum():delta.append(abs(np.sum(w*za)/w.sum())-abs(np.sum(w*zb)/w.sum()))
                case['paired']['planes'][plane]=dict(native=hb,truth_position=ha,
                    native_harmonic=harmonic(vb,ids,.08/63.,count),truth_position_harmonic=harmonic(va,ids,.08/63.,count),
                    harmonic_amplitude_change=float(abs(za.mean())-abs(zb.mean())),
                    amplitude_change_interval_95=np.quantile(delta,[.025,.975]).tolist(),
                    max_absolute_slope_change=float(np.max(np.abs(va-vb),initial=0)),
                    rms_slope_change=float(np.sqrt(np.mean((va-vb)**2))))
            fig.suptitle(f'Fixed truth and cluster identifiers | {name} | {count:,} source events\nNative total p = 50-250 GeV; candidate membership fixed; only measurement coordinates change')
            save(fig,args.output,name+'_same_hits',report)
            result['cases'][name]=case
    write_json(args.output/'summary.json',result)
    print(args.output/'summary.json')


if __name__=='__main__':main()
