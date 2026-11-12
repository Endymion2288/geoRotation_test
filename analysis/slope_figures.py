#!/usr/bin/env python3
"""English reference-layout figures from frozen S05 counts, with PNG/PDF exports."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


COLORS = dict(beforeGeo='#c62828',afterGeo='#1565c0')


def draw_hist(ax,hist,label,color,linestyle='-'):
    edges = np.asarray(hist['edges'])
    counts = np.asarray(hist['counts'])
    ax.stairs(counts,edges,label=label,color=color,linestyle=linestyle,linewidth=1.3)
    ax.set(xlim=(edges[0],edges[-1]),ylabel='Entries / bin')
    ax.ticklabel_format(axis='x',style='sci',scilimits=(-3,-3),useMathText=True)
    return counts.max(initial=0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--dsid',default='100044')
    parser.add_argument('--compact',action='store_true',help='Reference layouts and truth panels for auxiliary samples')
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    summaries = {v:json.loads((args.input/(v+'_summary.json')).read_text()) for v in COLORS}
    coverage = json.loads((args.input/'coverage.json').read_text())
    source_events = coverage['afterGeo']['source_events']
    pwindow=summaries['afterGeo']['metrics']['recoil_window_GeV']
    p_label=f'{pwindow[0]:g}-{pwindow[1]:g} GeV'
    scope = f'DSID {args.dsid} | {source_events:,} source events per configuration | Native Calypso'
    index = []
    plt.rcParams.update({'font.size':10,'axes.titlesize':11,'legend.fontsize':8,'axes.grid':True,'grid.alpha':.2})
    report = PdfPages(args.output/'slope_reference_report.pdf')

    def save(fig,name,title,reference,selection,definition):
        fig.suptitle(title+'\n'+scope,fontsize=12)
        fig.text(.5,.015,definition,ha='center',fontsize=8)
        fig.tight_layout(rect=(0,.035,1,.92))
        fig.savefig(args.output/(name+'.png'),dpi=180)
        fig.savefig(args.output/(name+'.pdf'))
        report.savefig(fig)
        index.append(dict(name=name,png=str(args.output/(name+'.png')),pdf=str(args.output/(name+'.pdf')),
            reference=reference,selection=selection,definition=definition))
        plt.close(fig)

    for plane,page in [('VetoNu',19),('st1',20)]:
        for mode in ['physical','legacy']:
            fig,axes = plt.subplots(2,2,figsize=(11,8),sharex=True,sharey='row')
            for row,selection in enumerate(['all','reco_p_window' if mode=='physical' else 'legacy_first_pz']):
                maximum=0
                for col,variant in enumerate(COLORS):
                    key = f'{plane}|{selection}|physical' if mode=='physical' else ('legacy_veto' if plane=='VetoNu' else 'raw_st1')+'|'+selection
                    hist=summaries[variant]['histograms'][key]
                    ax=axes[row,col]
                    maximum=max(maximum,draw_hist(ax,hist,variant,COLORS[variant]))
                    flow=hist['underflow']+hist['overflow']
                    ax.set_title(f'{variant} | '+('Full spectrum' if row==0 else ('Physical p = '+p_label if mode=='physical' else 'Legacy first-track pz window')))
                    ax.text(.03,.96,f"N = {hist['raw_entries']:,}\nInvalid = {hist['invalid']:,}; flow = {flow:,}",transform=ax.transAxes,va='top',fontsize=9)
                    ax.set_xlabel(r'$s_y=p_y^F/p_z^F$' if mode=='physical' else 'Native legacy angle [rad]')
                for ax in axes[row]:ax.set_ylim(0,max(1,maximum*1.18))
            definition = ('Physical FASER frame; fixed plane; per-track total upstream p window; raw entries.' if mode=='physical' else
                          'Native atan(raw py / raw pz); st1 uses the stored upstream surface; event-first raw pz cut; failures retained.')
            save(fig,f'{args.dsid}_{plane}_{mode}_reference',f'{plane}: {mode} angle definition',
                f'geometryRotation.pdf p{23 if args.dsid in ["110033","110051"] else page}','All native long tracks; lower row as labelled',definition)

        fig,axes=plt.subplots(1,2,figsize=(11,4.4))
        for ax,selection in zip(axes,['all','reco_p_window']):
            for variant in COLORS:
                hist=summaries[variant]['histograms'][f'{plane}|{selection}|physical']
                draw_hist(ax,hist,f"{variant} (N={hist['finite']:,})",COLORS[variant])
            ax.set(title=selection,xlabel=r'$s_y=p_y^F/p_z^F$')
            ax.legend()
        save(fig,f'{args.dsid}_{plane}_physical_overlay',f'{plane}: before/after overlay',
            f'geometryRotation.pdf p{page}','All native long tracks; physical p slice','Fixed physical plane; common bins; raw entries and efficiency retained.')

    fig,axes=plt.subplots(1,2,figsize=(11,4.4))
    for ax,variant in zip(axes,COLORS):
        hist=summaries[variant]['histograms']['literal_github_macro']
        draw_hist(ax,hist,variant,COLORS[variant])
        ax.set_title(variant+f"; N={hist['raw_entries']:,}, flow={hist['underflow']+hist['overflow']:,}")
        ax.set_xlabel('atan(px / pz) [rad]' if variant=='beforeGeo' else 'Y / sqrt(X^2 + Y^2)')
    save(fig,f'{args.dsid}_literal_github_expressions','Enabled GitHub expressions: different quantities',
        'GitHub root_txt/.../test_tanTheta.C','Before: 100 < raw pz < 250 GeV; after: no p cut',
        'Literal enabled expressions; before 100 bins, after 50 bins. These quantities do not form a physical theta_y comparison.')

    for plane in ['st1','VetoNu']:
        fig,axes=plt.subplots(1,2,figsize=(11,4.4),sharey=True)
        for ax,variant in zip(axes,COLORS):
            hist=summaries[variant]['histograms']
            for selection,color in [('reco_p_window','#1565c0'),('truth_p_window','#2e7d32'),('legacy_first_pz','#c62828')]:
                h=hist[f'{plane}|{selection}|physical']
                draw_hist(ax,h,selection+f"; N={h['raw_entries']:,}",color)
            ax.set(title=variant,xlabel=r'$s_y$');ax.legend()
        save(fig,f'{args.dsid}_{plane}_selection_migration',f'{plane}: momentum selection migration',
            'S05 A2/A5','Native long tracks; truth-p selection additionally requires signal matching',
            'Physical total-p and production-truth p = '+p_label+'; legacy cut uses event-first raw pz = 50-250 GeV. Raw entries.')

    for plane in ['st1','VetoNu','sensor']:
        fig,axes=plt.subplots(2,2,figsize=(11,8),sharex='row',sharey='row')
        for col,variant in enumerate(COLORS):
            hist=summaries[variant]['histograms']
            draw_hist(axes[0,col],hist[f'{plane}|matched_pure|physical'],'Matched pure reco',COLORS[variant])
            draw_hist(axes[0,col],hist[f'{plane}|truth_fiducial'],'Fiducial truth','#222222','--')
            draw_hist(axes[1,col],hist[f'{plane}|matched_pure|residual'],'Same-plane residual',COLORS[variant])
            axes[0,col].set(title=variant,xlabel=r'$s_y$')
            axes[1,col].set(xlabel=r'$s_y^{reco}-s_y^{truth}$')
            axes[0,col].legend();axes[1,col].legend()
        reference='True SiHit chord direction at the same sensor plane' if plane=='sensor' else 'Calypso-transported truth model; not exact local G4 truth'
        save(fig,f'{args.dsid}_{plane}_truth_residual',f'{plane}: truth, reconstruction and residual',
            'S05 A1/A2','Signal, fiducial H18, purity >= 0.7 and final Good95',reference+'; raw counts; truth includes unreconstructed particles.')

        if args.compact:continue

        fig,axes=plt.subplots(1,2,figsize=(11,4.4),sharey=True)
        for ax,variant in zip(axes,COLORS):
            h=summaries[variant]['histograms'][f'{plane}|matched|pull']
            draw_hist(ax,h,variant,COLORS[variant])
            ax.set(title=variant,xlabel=r'$(s_y^{reco}-s_y^{truth})/\sigma(s_y)$')
            ax.text(.03,.95,f"Finite = {h['finite']:,}; flow = {h['underflow']+h['overflow']:,}",transform=ax.transAxes,va='top')
        save(fig,f'{args.dsid}_{plane}_pull',f'{plane}: same-plane slope pull','S05 A2',
            'Native long tracks matched to signal; positive finite projected slope variance',reference+'; covariance from the actual native track; no pull-based selection.')

        fig,axes=plt.subplots(1,2,figsize=(11,4.4),sharey=True)
        for ax,variant in zip(axes,COLORS):
            for width,style in [(1e-4,'-'),(2e-4,'--'),(4e-4,':')]:
                hist=summaries[variant]['histograms'][f'{plane}|bin_{width:g}_offset_0']
                edges=np.asarray(hist['edges']);counts=np.asarray(hist['counts'])/width
                ax.stairs(counts,edges,label=f'Width {width:g}',linestyle=style)
            ax.set(title=variant,xlabel=r'$s_y$',ylabel='Entries / unit slope',xlim=(-.01,.01))
            ax.legend()
        save(fig,f'{args.dsid}_{plane}_binning_scan',f'{plane}: bin-width stability','S05 A2',
            'All native long tracks','Counts divided by bin width for shape comparison; no area normalization.')

        fig,axes=plt.subplots(1,2,figsize=(11,4.4),sharey=True)
        for ax,variant in zip(axes,COLORS):
            for offset,style in [(0.,'-'),(.5,'--')]:
                h=summaries[variant]['histograms'][f'{plane}|bin_0.0002_offset_{offset:g}']
                ax.stairs(np.asarray(h['counts'])/.0002,h['edges'],label=f'Bin offset {offset:g}',linestyle=style)
            ax.set(title=variant,xlabel=r'$s_y$',ylabel='Entries / unit slope',xlim=(-.01,.01));ax.legend()
        save(fig,f'{args.dsid}_{plane}_bin_offset',f'{plane}: half-bin translation','S05 A2',
            'All native long tracks','Same 200 microradian width; original and half-bin shifted edges; raw densities.')

        fig,axes=plt.subplots(1,2,figsize=(11,4.4))
        for variant in COLORS:
            structure=summaries[variant]['structure'][plane+'|all']
            ecdf=structure['ecdf']
            axes[0].plot(ecdf['x'],ecdf['cdf'],label=variant,color=COLORS[variant])
            kde=structure['kde']
            for bandwidth,style in [('0.0001','-'),('0.0002','--'),('0.0004',':')]:
                axes[1].plot(kde['grid'],kde['curves'][bandwidth],label=variant+' h='+bandwidth,color=COLORS[variant],linestyle=style)
        axes[0].set(xlabel=r'$s_y$',ylabel='ECDF',xlim=(-.01,.01))
        axes[1].set(xlabel=r'$s_y$',ylabel='Density',xlim=(-.01,.01))
        for ax in axes:ax.legend()
        save(fig,f'{args.dsid}_{plane}_ecdf_kde',f'{plane}: ECDF and KDE bandwidth scan','S05 A2',
            'All finite native long-track slopes','ECDF retains all finite values, including out-of-plot tails; binned Gaussian KDE uses 1e-5 grid.')

        fig,axes=plt.subplots(1,2,figsize=(11,4.4),sharey=True)
        for ax,variant in zip(axes,COLORS):
            structure=summaries[variant]['structure'][plane+'|reco_p_window']
            kde=structure['kde']
            for bandwidth,style in [('5e-05','-'),('0.0001','--'),('0.0002',':')]:
                ax.plot(kde['grid'],kde['curves'][bandwidth],label='Bandwidth '+bandwidth,linestyle=style)
            for p in summaries[variant]['peaks'].get(plane+'|reco_p_window',[]):
                w=p['window'];ax.axvspan(w['center']-w['half_width'],w['center']+w['half_width'],alpha=.10,color='#666666')
            ax.set(title=variant,xlabel=r'$s_y$',ylabel='Density',xlim=(-.004,.003));ax.legend(fontsize=8)
        save(fig,f'{args.dsid}_{plane}_momentum_slice_zoom',f'{plane}: momentum-slice peak zoom','S05 A2/A5',
            'Native long tracks; physical total p = '+p_label,'KDE bandwidth scan; shaded windows frozen on discovery entries only; no Gaussian-component interpretation.')

        fig,axes=plt.subplots(1,2,figsize=(11,4.4))
        for variant in COLORS:
            peaks=summaries[variant]['peaks'].get(plane+'|reco_p_window',[])
            if not peaks:continue
            centers=np.array([p['window']['center'] for p in peaks])*1000
            area=np.array([p['fraction']['estimate'] for p in peaks])
            ci=np.array([p['fraction']['interval_95'] for p in peaks]).T
            axes[0].errorbar(centers,area,yerr=np.maximum(0,np.array([area-ci[0],ci[1]-area])),fmt='o',label=variant,color=COLORS[variant])
            width=np.array([p['peak_std'] for p in peaks])*1e6
            wci=np.array([p['moment_interval']['std_95'] for p in peaks]).T*1e6
            axes[1].errorbar(centers,width,yerr=np.maximum(0,np.array([width-wci[0],wci[1]-width])),fmt='o',label=variant,color=COLORS[variant])
        axes[0].set(xlabel='Frozen window center [mrad]',ylabel='Validation track fraction')
        axes[1].set(xlabel='Frozen window center [mrad]',ylabel='Conditional std inside window [microrad]')
        for ax in axes:ax.legend()
        save(fig,f'{args.dsid}_{plane}_peak_parameters',f'{plane}: held-out peak-window measurements','S05 A2/A5',
            'Physical p = '+p_label+'; validation entries only','95% source-event intervals; width is an unbinned conditional moment inside each frozen window, not intrinsic peak width.')

        for group in ['production_p_GeV','module_eta','cluster_width','purity']:
            fig,axes=plt.subplots(1,2,figsize=(11,4.4),sharex=True,sharey=True)
            for ax,variant in zip(axes,COLORS):
                boundaries=summaries[variant]['groups'][group]
                for j,(lo,hi) in enumerate(zip(boundaries[:-1],boundaries[1:])):
                    h=summaries[variant]['histograms'][f'{plane}|group_{group}_{j}']
                    if h['finite']:
                        draw_hist(ax,h,f'{lo} to {hi} (N={h["finite"]:,})',plt.cm.tab10(j%10))
                ax.set(title=variant,xlabel=r'$s_y$')
                ax.legend()
            save(fig,f'{args.dsid}_{plane}_group_{group}',f'{plane}: contributions by {group}',
                'S05 A3','Long tracks matched to primary signal','Raw group contributions sum to the displayed selected population; no per-group normalization.')

    for plane in ([] if args.compact else ['st1','VetoNu']):
        for selection in ['all','reco_p_window','truth_p_window']:
            key=plane if selection=='all' else plane+'|'+selection
            available=all(key in summaries[v]['response'] for v in COLORS)
            if not available:continue
            fig,axes=plt.subplots(2,2,figsize=(11,7.5),sharex=True,sharey='row')
            for col,variant in enumerate(COLORS):
                response=summaries[variant]['response'][key]
                edges=np.array(response['edges']);observed=np.array(response['observed_counts']);predicted=np.array(response['predicted_counts'])
                ax=axes[0,col]
                ax.stairs(observed,edges,label='Held-out reconstructed',color=COLORS[variant])
                ax.stairs(predicted,edges,label='Discovery response prediction',color='#222222',linestyle='--')
                lower,upper = np.array(response['predicted_counts_interval_95'])
                ax.fill_between(edges,np.r_[lower,lower[-1]],np.r_[upper,upper[-1]],step='post',alpha=.18,color='#444444',label='Event bootstrap 95%')
                ax.set(title=variant,ylabel='Entries / bin',xlim=(-.01,.01));ax.legend()
                ax.text(.03,.95,f"Uncovered truth = {response['uncovered_validation_truth']:,}\nShape TV = {response['normalized_shape_total_variation']:.3%}",transform=ax.transAxes,va='top',fontsize=8)
                ax=axes[1,col];difference=observed-predicted
                lower,upper=np.array(response['observed_minus_predicted_interval_95'])
                ax.stairs(difference,edges,color=COLORS[variant])
                ax.fill_between(edges,np.r_[lower,lower[-1]],np.r_[upper,upper[-1]],step='post',alpha=.18,color='#444444')
                ax.axhline(0,color='k',lw=.7);ax.set(xlabel=r'$s_y$',ylabel='Observed minus predicted')
            save(fig,f'{args.dsid}_{plane}_{selection}_held_out_response',f'{plane}: independent response prediction ({selection})',
                'S05 A3/A5','Long tracks matched to primary signal; purity >= 0.7; '+selection,
                f"Train <{coverage['afterGeo']['held_out_split_sample_entry']:,}; test >= boundary. Conditions: true p, slope, sensor x/y. Missing support is not an explained population.")
    for plane in ([] if args.compact else ['st1','VetoNu','sensor']):
        for condition in ['production_p_GeV','true_sensor_y_mm','module_eta']:
            fig,axes=plt.subplots(1,2,figsize=(11,4.4),sharex=True,sharey=True)
            vmax=max(np.max(summaries[v]['maps'][plane+'|'+condition]['counts'],initial=0) for v in COLORS)
            for ax,variant in zip(axes,COLORS):
                m=summaries[variant]['maps'][plane+'|'+condition]
                image=ax.pcolormesh(m['slope_edges'],m['condition_edges'],np.asarray(m['counts']).T,rasterized=True,vmin=0,vmax=max(1,vmax))
                fig.colorbar(image,ax=ax,label='Entries')
                ax.set(title=variant,xlabel=r'$s_y$',ylabel=condition)
            save(fig,f'{args.dsid}_{plane}_map_{condition}',f'{plane}: slope versus {condition}',
                'S05 A3/A5','Long tracks matched to primary signal','Common bins; raw entries; invalid and out-of-map entries are counted in the summary.')
    report.close()
    (args.output/'index.json').write_text(json.dumps(dict(dsid=args.dsid,source_events=source_events,
        figure_count=len(index),report=str(args.output/'slope_reference_report.pdf'),figures=index),indent=2)+'\n')
    print(f'Exported {len(index)} PNG/PDF figures')


if __name__=='__main__':
    main()
