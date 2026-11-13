#!/usr/bin/env python3
"""Alma-era charge measurement from all available actual Calypso Ntuples."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta

from slope_tables import histogram, write_json


PBINS=np.array([0,50,100,250,500,1000,10000.])


def rate(wrong,total):
    wrong,total=int(wrong),int(total)
    intervals={}
    for label,alpha in [('68',.32),('95',.05)]:
        intervals[label]=[float(beta.ppf(alpha/2,wrong,total-wrong+1)) if wrong else 0.,
                         float(beta.ppf(1-alpha/2,wrong+1,total-wrong)) if wrong<total else 1.] if total else [None,None]
    return dict(wrong=wrong,total=total,value=wrong/total if total else None,intervals=intervals)


def load_group(root,dsid,tag,inventory):
    sample=next(s for s in inventory['samples'] if s['dsid']==dsid)
    indices={Path(f['path']).name:i for i,f in enumerate(sample['files'])}
    tables={name:[] for name in ['events','tracks']}
    for source in inventory['historical_phys'][dsid+'|'+tag]['files']:
        index=indices[Path(source['path']).name.replace('PHYS.root','RDO.root')]
        folder=root/'history'/dsid/tag/f'f{index:02d}'
        assert (folder/'exit_code.txt').read_text().strip()=='0'
        for name in tables:
            with np.load(folder/(name+'.npz')) as archive:data={k:archive[k] for k in archive.files}
            if name=='events':np.testing.assert_array_equal(data['source_entry'],np.arange(source['entries']))
            data['sample_entry']=100000*index+data['source_entry']
            data['file_index']=np.full(len(data['source_entry']),index,dtype=np.int32)
            tables[name].append(data)
    return {name:{k:np.concatenate([r[k] for r in rows]) for k in rows[0]} for name,rows in tables.items()}


def measure(tables,dsid):
    t,e=tables['tracks'],tables['events']
    event_index=np.searchsorted(e['sample_entry'],t['sample_entry'])
    np.testing.assert_array_equal(e['sample_entry'][event_index],t['sample_entry'])
    pdg=-13 if dsid=='100116' else 13;charge=-np.sign(pdg)
    with np.errstate(divide='ignore',invalid='ignore'):
        base=(t['ndof']>0)&(t['chi2']/t['ndof']<25)&(t['n_layers']>=7)&(t['ndof']+5>=12)
    valid=np.isfinite(t['momentum_GeV'])&(t['momentum_GeV']>0)&np.isin(t['charge'],[-1,1])
    matched=(t['pdg']==pdg)&(t['barcode']>0)&(t['barcode']<200000)
    wrong=t['charge']*charge<0
    results={};selections={};codes={}
    for radius in [100,95]:
        good=base&(t['rmax_mm']<radius)
        n_good=np.bincount(event_index,weights=good,minlength=len(e['sample_entry'])).astype(int)
        one=good&(n_good[event_index]==1)&t['truth_fiducial']
        for label,mask in [('legacy_fiducial_one_good',one),
            ('matched_fiducial_one_good',one&matched),
            ('pure_matched_fiducial_one_good',one&matched&(t['purity_ndof']>.7))]:
            key=f'{label}_r{radius}'
            selected=mask&valid
            row=rate(np.count_nonzero(selected&wrong),selected.sum())
            row.update(invalid_selected=int(np.count_nonzero(mask&~valid)),
                selected_nonprimary_or_wrong_species=int(np.count_nonzero(selected&~matched)),
                exactly_one_good_events=int((n_good==1).sum()),multiple_good_events=int((n_good>1).sum()),
                discovery=rate(np.count_nonzero(selected&wrong&(t['file_index']<3)),
                               np.count_nonzero(selected&(t['file_index']<3))),
                held_out=rate(np.count_nonzero(selected&wrong&(t['file_index']>=3)),
                             np.count_nonzero(selected&(t['file_index']>=3))))
            row['truth_p_bins']=[dict(low_GeV=float(lo),high_GeV=float(hi),
                **rate(np.count_nonzero(selected&wrong&(t['truth_momentum_GeV']>=lo)&(t['truth_momentum_GeV']<hi)),
                np.count_nonzero(selected&(t['truth_momentum_GeV']>=lo)&(t['truth_momentum_GeV']<hi))))
                for lo,hi in zip(PBINS[:-1],PBINS[1:])]
            results[key]=row;selections[key]=selected
            code=np.zeros(len(e['sample_entry']),dtype=np.int8)
            code[event_index[selected]]=np.where(wrong[selected],2,1)
            codes[key]=code
    chosen=selections['legacy_fiducial_one_good_r100']
    with np.errstate(divide='ignore',invalid='ignore'):
        qop=t['charge']/t['momentum_GeV']
    hist=dict(signed_qop=histogram(qop[chosen],np.linspace(-.03,.03,301)),
        reco_p=histogram(t['momentum_GeV'][chosen],np.linspace(0,5000,201)),
        truth_p=histogram(t['truth_momentum_GeV'][chosen],np.linspace(0,5000,201)))
    h,xe,ye=np.histogram2d(t['truth_momentum_GeV'][chosen],t['momentum_GeV'][chosen],
                           bins=[np.linspace(0,2000,101)]*2)
    hist['truth_reco_p']=dict(counts=h.astype(int).tolist(),truth_edges=xe.tolist(),reco_edges=ye.tolist(),
                              selected=int(chosen.sum()),outside=int(chosen.sum()-h.sum()))
    return dict(source_events=len(e['sample_entry']),long_tracks=len(t['source_entry']),
        fiducial_primary_key_unknown=int((~e['primary_fiducial_valid']).sum()),
        source_identity='Production file index + nt entry; repeated run/event retained',
        good_track_definition='chi2/ndof<25, layers>=7, ndof+5>=12; radii shown separately',
        historical_threshold_status='Source-based executable reconstruction; original Alma analysis thresholds not recovered',
        conditional_selections=results,histograms=hist),codes


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    inventory=json.loads((args.root/'input_inventory.json').read_text())
    groups={};payload={};codes={};pairs={}
    for dsid in ['100116','100117']:
        for tag in ['s0012-dev','s0012-r0019']:
            key=dsid+'|'+tag
            groups[key]=load_group(args.root,dsid,tag,inventory)
            payload[key],codes[key]=measure(groups[key],dsid)
        old,new=[dsid+'|'+tag for tag in ['s0012-r0019','s0012-dev']]
        common,oi,ni=np.intersect1d(groups[old]['events']['sample_entry'],groups[new]['events']['sample_entry'],return_indices=True)
        for field in ['run','event']:
            np.testing.assert_array_equal(groups[old]['events'][field][oi],groups[new]['events'][field][ni])
        selections={}
        for selection in codes[old]:
            a,b=codes[old][selection][oi],codes[new][selection][ni]
            matrix=np.bincount(a*3+b,minlength=9).reshape(3,3)
            shared=int(matrix[1:,1:].sum())
            selections[selection]=dict(transition_counts=matrix.tolist(),
                rows_columns=['not_selected','correct','wrong'],
                old_common=rate((a==2).sum(),(a>0).sum()),
                new_common=rate((b==2).sum(),(b>0).sum()),
                old_correct_new_wrong_over_common_selected=rate(matrix[1,2],shared),
                old_correct_new_wrong_over_new_selected=rate(matrix[1,2],(b>0).sum()))
        pairs[dsid]=dict(common_source_events=len(common),identity_checked=True,selections=selections)
    result=dict(groups=payload,paired_history=pairs,
        historical_labels=dict(old='s0012-r0019',new='s0012-dev',source='Alma report p2'),
        published_Alma_percent={'100116':{'old':.101,'new':.197},'100117':{'old':.217,'new':.359}},
        coverage=inventory['historical_phys'],historical_comparison_is_not_geometry_rotation=True)
    write_json(args.output/'summary.json',result)
    plot(args.output,result)
    print(json.dumps({k:v['conditional_selections']['legacy_fiducial_one_good_r100'] for k,v in payload.items()},indent=2))


def plot(output,result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LogNorm
    figures=output/'figures';figures.mkdir()
    index=[]
    with PdfPages(figures/'historical_charge_report.pdf') as pdf:
        def save(fig,name,reference):
            fig.tight_layout(rect=(0,.025,1,.94))
            fig.text(.5,.01,'Calypso PHYS | fiducial + exactly one Good Track | r < 100 mm | source-based 25/7/12 thresholds',ha='center',fontsize=8)
            for suffix in ['png','pdf']:fig.savefig(figures/(name+'.'+suffix),dpi=180)
            pdf.savefig(fig);plt.close(fig)
            index.append(dict(name=name,reference=reference,png=name+'.png',pdf=name+'.pdf'))
        for observable,xlabel,reference in [('signed_qop',r'$q/p$ [GeV$^{-1}$]','Alma pp16,23-25'),
            ('reco_p','Reconstructed momentum [GeV]','Alma pp22,27'),('truth_p','Truth momentum [GeV]','Alma p21')]:
            fig,axes=plt.subplots(1,2,figsize=(11,4.7))
            for dsid,ax in zip(['100116','100117'],axes):
                for tag,color in [('s0012-r0019','#c62828'),('s0012-dev','#1565c0')]:
                    group=result['groups'][dsid+'|'+tag];h=group['histograms'][observable]
                    ax.stairs(h['counts'],h['edges'],label=f"{tag} (N={h['raw_entries']:,})",color=color)
                ax.set(title=r'$\mu^+$' if dsid=='100116' else r'$\mu^-$',xlabel=xlabel,ylabel='Entries / bin',yscale='log')
                ax.legend(fontsize=8);ax.grid(alpha=.2)
            fig.suptitle('Historical Calypso: '+observable)
            save(fig,'historical_'+observable,reference)
        fig,axes=plt.subplots(2,2,figsize=(10,9),sharex=True,sharey=True)
        for i,dsid in enumerate(['100116','100117']):
            for j,tag in enumerate(['s0012-dev','s0012-r0019']):
                h=result['groups'][dsid+'|'+tag]['histograms']['truth_reco_p'];ax=axes[i,j]
                ax.pcolormesh(h['truth_edges'],h['reco_edges'],np.asarray(h['counts']).T,norm=LogNorm(vmin=1),rasterized=True)
                ax.plot([0,2000],[0,2000],color='red',lw=.7)
                ax.set(title=dsid+' | '+tag,xlabel='Truth momentum [GeV]',ylabel='Reconstructed momentum [GeV]')
        fig.suptitle('Historical Calypso: truth versus reconstructed momentum')
        save(fig,'historical_truth_reco_momentum','Alma pp19-20')
        fig,axes=plt.subplots(1,2,figsize=(11,4.7))
        for dsid,ax in zip(['100116','100117'],axes):
            for tag,color in [('s0012-r0019','#c62828'),('s0012-dev','#1565c0')]:
                rows=result['groups'][dsid+'|'+tag]['conditional_selections']['legacy_fiducial_one_good_r100']['truth_p_bins']
                valid=[(i,r) for i,r in enumerate(rows) if r['total']]
                x=np.array([i for i,r in valid]);y=np.array([r['value'] for i,r in valid])*100
                ci=np.array([r['intervals']['68'] for i,r in valid])*100
                ax.errorbar(x,y,yerr=np.stack([y-ci[:,0],ci[:,1]-y]),fmt='o-',color=color,label=tag)
            ax.set(title=dsid,xlabel='Truth momentum bin [GeV]',ylabel='Charge misidentification [%]',xticks=range(6),
                xticklabels=['0-50','50-100','100-250','250-500','500-1000','1000+'])
            ax.tick_params(axis='x',rotation=25);ax.legend(fontsize=8);ax.grid(alpha=.2)
        fig.suptitle('Historical Calypso: charge misidentification (68% intervals)')
        save(fig,'historical_misid_truth_p','Alma pp17-18 plus truth-p conditioning')
    write_json(figures/'figure_index.json',index)
    (figures/'figure_index.md').write_text('\n'.join(f"- [{r['name']}]({r['png']}) · [PDF]({r['pdf']}) · {r['reference']}" for r in index)+'\n')


if __name__=='__main__':
    main()
