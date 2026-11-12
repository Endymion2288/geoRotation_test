#!/usr/bin/env python3
"""Locate fine slope modulation in real Calypso truth, segments and tracks."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_tables import histogram, join_truth, merge_tables, slope
from slope_analysis import derive
from slope_controls import save


def harmonic(values,entries,period,event_count):
    valid=np.isfinite(values)
    z=np.exp(2j*np.pi*values[valid]/period)
    ids=entries[valid]
    n=len(z)
    if not n:return dict(entries=0)
    mean=z.mean()
    real=np.bincount(ids,weights=z.real-mean.real,minlength=event_count)/n
    imag=np.bincount(ids,weights=z.imag-mean.imag,minlength=event_count)/n
    cov=np.stack([real,imag])@np.stack([real,imag]).T*event_count/max(1,event_count-1)
    return dict(entries=n,period=period,real=float(mean.real),imaginary=float(mean.imag),
                amplitude=float(abs(mean)),phase_rad=float(np.angle(mean)),cluster_covariance=cov.tolist())


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--variant',default='afterGeo')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    tables,coverage=merge_tables(args.root,args.variant)
    d=derive(tables,args.variant,'100044')
    truth,tracks,segments=tables['truth'],tables['tracks'],tables['segments']
    count=coverage['source_events']
    split=json.loads((args.root/'manifest.json').read_text())['held_out_split_sample_entry']
    edges=np.linspace(-.008,.008,321)
    period=.08/63.
    result=dict(source_events=count,variant=args.variant,period_from_geometry=period,
        pitch_mm=.08,three_layer_lever_arm_mm=63.,
        segment_reference='Raw native segment mean at its stored station plane; no fixed-plane propagation; momentum is a fixed 100 GeV convention',
        results={})
    with PdfPages(args.output/'slope_origin.pdf') as report:
        for selection in ['all','truth_p_window']:
            tp=np.linalg.norm(truth['production_momentum_F_GeV'],axis=1)
            population=d['signal_truth'] & (np.sum(truth['truth_hit_counts'][:,1:],axis=1)>=18)
            if selection=='truth_p_window':population &= (tp>50)&(tp<250)
            matched=d['primary'] & d['long'] & d['pure'] & population[d['index']]
            efficient=np.bincount(d['index'][matched],minlength=len(population))>0
            fig,axes=plt.subplots(2,2,figsize=(12,8))
            for i,plane in enumerate(['st1','VetoNu']):
                key=plane+'|'+selection
                series={
                    'All H18 truth':(d['truth_slopes'][:,i],truth['sample_entry'],population),
                    'Reconstructed truth':(d['truth_slopes'][:,i],truth['sample_entry'],population&efficient),
                    'Reconstruction':(d['reconstructed'][:,i],tracks['sample_entry'],matched)}
                result['results'][key]={}
                for name,(values,entries,mask) in series.items():
                    h=histogram(values[mask],edges)
                    result['results'][key][name]=dict(histogram=h,harmonics={})
                    for half,em in [('discovery',entries<split),('validation',entries>=split)]:
                        choose=mask&em
                        result['results'][key][name]['harmonics'][half]=harmonic(values[choose],entries[choose],period,count)
                    axes[0,i].stairs(h['counts'],h['edges'],label=name)
                den=np.asarray(result['results'][key]['All H18 truth']['histogram']['counts'])
                num=np.asarray(result['results'][key]['Reconstructed truth']['histogram']['counts'])
                good=den>0;center=.5*(edges[1:]+edges[:-1]);eff=np.divide(num,den,out=np.zeros_like(num,dtype=float),where=good)
                scale=1+1.96**2/np.maximum(1,den)
                midpoint=(eff+1.96**2/(2*np.maximum(1,den)))/scale
                half=1.96*np.sqrt(eff*(1-eff)/np.maximum(1,den)+1.96**2/(4*np.maximum(1,den)**2))/scale
                lower,upper=midpoint-half,midpoint+half
                axes[1,i].errorbar(center[good],eff[good],yerr=np.maximum(0,np.array([eff[good]-lower[good],upper[good]-eff[good]])),fmt='.',markersize=3,lw=.5)
                axes[0,i].set(title=plane,xlabel='Physical sy',ylabel='Entries / 50 microradian bin')
                axes[1,i].set(xlabel='Same-plane model truth sy',ylabel='Matched pure long-track efficiency',ylim=(0,1.05))
                result['results'][key]['efficiency']=dict(edges=edges.tolist(),numerator=num.tolist(),denominator=den.tolist(),
                    lower_95=lower.tolist(),upper_95=upper.tolist(),
                    interval_method='95% Wilson binomial interval for the one-primary-muon-per-event sample')
                axes[0,i].legend(fontsize=8)
            for ax in axes.flat:ax.grid(alpha=.2)
            fig.suptitle(f'Native Calypso fine slope structure | {args.variant} | {selection} | {count:,} source events\nH18 truth population; native long tracks with purity >= 0.7; fixed 50 microradian bins')
            save(fig,args.output,'truth_efficiency_'+selection,report)

            fig,axes=plt.subplots(1,2,figsize=(11,4.5))
            for i,(plane,ax) in enumerate(zip(['st1','VetoNu'],axes)):
                true_values=d['truth_slopes'][d['index'],i]
                reco_values=d['reconstructed'][:,i]
                keep=matched & np.isfinite(true_values) & np.isfinite(reco_values)
                map_edges=np.linspace(-.01,.01,101)
                counts,xe,ye=np.histogram2d(true_values[keep],reco_values[keep],bins=[map_edges,map_edges])
                image=ax.pcolormesh(xe,ye,counts.T,rasterized=True)
                fig.colorbar(image,ax=ax,label='Entries')
                ax.set(title=plane,xlabel='Same-plane model truth sy',ylabel='Reconstructed sy')
                result['results'][plane+'|'+selection]['truth_reco_map']=dict(edges=map_edges.tolist(),counts=counts.astype(int).tolist(),
                    finite_pairs=int(keep.sum()),out_of_map=int(keep.sum()-counts.sum()),sumw2_equals_unweighted_counts=True)
            fig.suptitle(f'Truth versus reconstruction | {args.variant} | {selection} | {count:,} source events')
            save(fig,args.output,'truth_reco_map_'+selection,report)

            index,linked=join_truth({'truth':truth,'tracks':segments})
            fig,axes=plt.subplots(1,3,figsize=(14,4.5))
            for station,ax in zip([1,2,3],axes):
                mask=linked & population[index] & (segments['station']==station)
                values=slope(segments['parameters_F'][:,3:])
                true_values=slope(truth['truth_station_direction_F'][:,station])
                h=histogram(values[mask],edges);t=histogram(true_values[population],edges)
                key='segment'+str(station)+'|'+selection
                result['results'][key]=dict(histogram=h,truth=t,
                    harmonic=harmonic(values[mask],segments['sample_entry'][mask],period,count),
                    residual=histogram((values-true_values[index])[mask],np.linspace(-.003,.003,241)))
                ax.stairs(h['counts'],edges,label='Raw native segment')
                ax.stairs(t['counts'],edges,label='Station-mean SiHit truth',linestyle='--')
                ax.set(title=f'Station {station}',xlabel='Raw local sy',ylabel='Entries / bin');ax.legend(fontsize=8);ax.grid(alpha=.2)
            fig.suptitle(f'Native segments before seeding | {args.variant} | {selection} | {count:,} source events\nRaw stored station surfaces; multiple segments per particle retained; SiHit station-mean truth reference')
            save(fig,args.output,'segments_'+selection,report)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(args.output/'summary.json')


if __name__=='__main__':main()
