#!/usr/bin/env python3
"""Paired checks of actual Calypso writer and candidate smoothing controls."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from slope_controls import archive, save
from slope_tables import histogram
from slope_analysis import derive


def equal_arrays(a,b,exclude=()):
    assert a.keys()==b.keys()
    changed=[]
    for k in a:
        if k not in exclude and not np.array_equal(a[k],b[k],equal_nan=True):
            changed.append(k)
    return changed


def writer(root,output):
    result={'definition':'Actual Calypso NtupleDumper flag false/true on the same RDO entries',
            'pairs':[]}
    with PdfPages(output/'writer_control.pdf') as report:
        for suffix in ['', '_validation']:
            old=root/'afterGeo'/('legacy'+suffix)
            new=root/'afterGeo'/('fixed'+suffix)
            pair={'suffix':suffix,'observer_tables':{}}
            for table in ['events','truth','tracks','segments']:
                a,b=archive(old/(table+'.npz')),archive(new/(table+'.npz'))
                changed=equal_arrays(a,b)
                assert not changed,(table,changed)
                pair['observer_tables'][table]={'all_saved_branches_exactly_equal':True,
                                               'branches':len(a),'rows':len(a['source_entry'])}
            a,b=archive(old/'native_ntuple_branches.npz'),archive(new/'native_ntuple_branches.npz')
            angles=['Track_ThetaX_atVetoNu','Track_ThetaY_atVetoNu']
            assert not equal_arrays(a,b,angles)
            pair['saved_non_angle_native_branches_exactly_equal']=[k for k in a if k not in angles]
            pair['histograms']={}
            fig,axes=plt.subplots(1,2,figsize=(11,4.5))
            for key,ax in zip(angles,axes):
                pair['histograms'][key]={}
                for label,data,color in [('Legacy writer',a,'#c62828'),('Physical writer',b,'#1565c0')]:
                    h=histogram(data[key],np.linspace(-.01,.01,101))
                    pair['histograms'][key][label]=h
                    ax.stairs(h['counts'],h['edges'],label=f"{label}; flow={h['underflow']+h['overflow']}",color=color)
                ax.set(title=key,xlabel='Angle [rad]',ylabel='Entries / bin')
                ax.legend(fontsize=8);ax.grid(alpha=.2)
            fig.suptitle(f"Actual Calypso writer only | {len(a['track_offsets'])-1} source events | {len(a[angles[0]])} native long tracks")
            save(fig,output,'writer'+(suffix or '_smoke'),report)
            result['pairs'].append(pair)
    # Reconstruct only the saved real writer branches as vectors, so that the
    # unchanged original macro executes its own loops, cuts and ROOT Fill.
    import ROOT
    ROOT.gROOT.SetBatch(True)
    names=['110033_ntOut_officialRec_3st_forward_new_PHYS.root',
           '110033_smallerCov_seedParams_3stForward_new_PHYS.root']
    expectation={}
    for variant,name in zip(['beforeGeo','afterGeo'],names):
        data=archive(root/variant/'legacy/native_ntuple_branches.npz')
        f=ROOT.TFile(name,'RECREATE');tree=ROOT.TTree('nt','Saved actual native writer branches')
        keys=[k for k in data if k!='track_offsets']
        unused_double=['t_st1_py','t_st1_pz','t_p','Track_p0']
        unused_int=['Track_charge','t_pdg','t_barcode']
        buffers={k:ROOT.std.vector('double')() for k in keys+unused_double}
        buffers.update({k:ROOT.std.vector('int')() for k in unused_int})
        for key,buffer in buffers.items():tree.Branch(key,buffer)
        for lo,hi in zip(data['track_offsets'][:-1],data['track_offsets'][1:]):
            for buffer in buffers.values():buffer.clear()
            for key in keys:
                for value in data[key][lo:hi]:buffers[key].push_back(float(value))
            tree.Fill()
        tree.Write();f.Close()
        tables={k:archive(root/variant/'legacy'/(k+'.npz')) for k in ['events','truth','tracks','segments']}
        for table in tables.values():table['sample_entry']=table['source_entry']
        d=derive(tables,variant,'100044')
        expectation[variant]=histogram(d['literal'][d['literal_mask']],np.linspace(-.01,.01,101 if variant=='beforeGeo' else 51))
    macro=root.parents[0]/'github_test_tanTheta.C'
    assert macro.exists()
    ROOT.gROOT.ProcessLine('.L '+str(macro))
    ROOT.test_tanTheta()
    result['original_macro']={'macro':str(macro),
        'input':'Small trees made from saved real writer branches; unused truth and p0 branches empty',
        'histograms':{}}
    for variant,name in zip(['beforeGeo','afterGeo'],['h1','h2']):
        h=ROOT.gROOT.FindObject(name)
        assert h
        expected=expectation[variant]
        counts=[int(h.GetBinContent(i)) for i in range(h.GetNbinsX()+2)]
        assert counts==[expected['underflow'],*expected['counts'],expected['overflow']]
        assert int(h.GetEntries())==expected['raw_entries']
        result['original_macro']['histograms'][variant]=dict(counts_with_flows=counts,
            entries=int(h.GetEntries()),python_literal_all_bins_exactly_equal=True)
    return result


def smoothing(root,output):
    baseline=root/'stage_controls_v1/shards'
    revised=root/'smoothing_control_v1/shards'
    result={'definition':'Same native Calypso candidates; only IsolateSmoothing false/true',
            'source_events':2000,'pairs':[]}
    with PdfPages(output/'smoothing_control.pdf') as report:
        for shard in ['s00','s02']:
            a,b=archive(baseline/shard/'states.npz'),archive(revised/shard/'states.npz')
            identities=['entry','seed','candidate','stage','state','geometry','measurement','kind']
            assert all(np.array_equal(a[k],b[k]) for k in identities)
            pair={'shard':shard,'state_rows':len(a['entry']),
                'all_state_identities_exactly_equal':True,
                'stage_mean_changes':{},'angular_identity':{}}
            for code in np.unique(a['stage']):
                mask=a['stage']==code
                delta=np.abs(a['values'][mask]-b['values'][mask])
                pair['stage_mean_changes'][str(code)]={'rows':int(mask.sum()),'changed_rows':int(np.any(delta!=0,axis=1).sum()),'max_abs_delta':float(delta.max(initial=0))}
            r1,r2=archive(baseline/shard/'tracks.npz'),archive(revised/shard/'tracks.npz')
            pair['changed_observer_track_branches']=equal_arrays(r1,r2)
            for key in ['source_entry','barcode','charge','n_hits','hit_offsets','hit_identifiers','upstream_raw','projection','projection_valid']:
                assert np.array_equal(r1[key],r2[key],equal_nan=True),key
            pair['final_mean_charge_hits_and_projections_exactly_equal']=True
            pair['rmax_changed_tracks']=int(np.count_nonzero(r1['rmax_mm']!=r2['rmax_mm']))
            long=r1['station_bits']&14==14
            def good95(r):
                with np.errstate(divide='ignore',invalid='ignore'):
                    return long & (r['n_layers']>=7) & (r['ndof']+5>=12) & (r['rmax_mm']<95) & (r['chi2']/r['ndof']<25)
            pair['Good95_changed_tracks']=int(np.count_nonzero(good95(r1)!=good95(r2)))
            keys=set(zip(a['entry'][a['stage']==11][long],a['seed'][a['stage']==11][long],a['candidate'][a['stage']==11][long]))
            cohort=np.array([(e,s,c) in keys for e,s,c in zip(a['entry'],a['seed'],a['candidate'])])
            mask=cohort & (a['stage']==5)
            fig,axes=plt.subplots(1,2,figsize=(11,4.5))
            for label,table,color in [('Shared slots',a,'#c62828'),('Candidate isolated',b,'#1565c0')]:
                error=np.max(np.abs(table['smoothing_error'][mask,2:4]),axis=1)
                finite=error[np.isfinite(error)]
                pair['angular_identity'][label]=dict(measurement_states=int(mask.sum()),
                    finite_checks=len(finite),max_abs_error_rad=float(finite.max(initial=0)),
                    above_1e_minus_8=int(np.count_nonzero(finite>1e-8)))
                h=histogram(np.tan(table['values'][mask,2]),np.linspace(-.01,.01,201))
                pair.setdefault('histograms',{})[label]=h
                axes[0].stairs(h['counts'],h['edges'],label=label,color=color)
                axes[1].hist(np.log10(np.maximum(finite,1e-18)),bins=np.linspace(-18,0,73),histtype='step',label=label,color=color)
            axes[0].set(xlabel='CKF smoothed physical sy at measurement surfaces',ylabel='Entries / bin')
            axes[1].set(xlabel='log10 max angular smoother identity error [rad]',ylabel='Checks / bin')
            for ax in axes:ax.legend();ax.grid(alpha=.2)
            fig.suptitle(f'Native Calypso smoothing control | {shard} | 1,000 source events\nFinal mean, charge, hits and physical projections remain paired')
            save(fig,output,shard+'_smoothing',report)
            result['pairs'].append(pair)
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['writer','smoothing'])
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    result=globals()[args.mode](args.root,args.output)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(args.output/'summary.json')


if __name__=='__main__':main()
