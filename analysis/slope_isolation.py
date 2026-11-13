#!/usr/bin/env python3
"""Compare full native Calypso angle samples with smoothing OFF and ON."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from slope_tables import histogram, slope, write_json


def tasks(root,variant):
    rows=json.loads((root/(variant+'_tasks.json')).read_text())
    sources=[]
    for row in rows:
        output=Path(row['output'])
        assert (output/'exit_code.txt').read_text().strip()=='0'
        source=json.loads((output/'source_task.json').read_text())
        sources.append((source,output))
    return {(s['input'],s['skip'],s['events']):(s,o) for s,o in sources}


def compare(args):
    inputs={mode:tasks(getattr(args,mode+'_root'),args.variant) for mode in ['off','on']}
    assert inputs['off'].keys()==inputs['on'].keys()
    changed=Counter();maximum={};counts=Counter();histograms={}
    fields_checked=set();examples=[]
    edges=np.linspace(-.01,.01,101)
    for key in sorted(inputs['on']):
        pair={mode:inputs[mode][key][1] for mode in inputs}
        arrays={}
        for mode,output in pair.items():
            run=json.loads((output/'run.json').read_text())
            assert run['success'] and run['isolate_smoothing']==(mode=='on')
        for table in ['events','truth','tracks','segments']:
            arrays[table]={}
            for mode,output in pair.items():
                with np.load(output/(table+'.npz')) as archive:
                    arrays[table][mode]={k:archive[k] for k in archive.files}
            a,b=arrays[table]['off'],arrays[table]['on']
            for name in a.keys() & b.keys():
                tag=table+'/'+name
                assert a[name].shape==b[name].shape, (key,tag)
                same=(a[name]==b[name]) | (np.isnan(a[name]) & np.isnan(b[name]))
                fields_checked.add(tag)
                if not same.all():
                    changed[tag]+=int((~same).sum())
                    finite=np.isfinite(a[name]) & np.isfinite(b[name])
                    maximum[tag]=max(maximum.get(tag,0.),float(np.max(np.abs(
                        a[name][finite].astype(float)-b[name][finite].astype(float)),initial=0)))
                    if len(examples)<20:
                        examples.append(dict(input=key[0],skip=key[1],field=tag,
                            changed_elements=int((~same).sum())))
            for name in ['source_entry','run','event']:
                np.testing.assert_array_equal(a[name],b[name])
            counts[table]+=len(a['source_entry'])
        selection_masks={}
        for mode in ['off','on']:
            t=arrays['tracks'][mode]
            long=t['station_bits'] & 14 == 14
            momentum=np.linalg.norm(t['upstream_F'][:,3:],axis=1)
            with np.errstate(divide='ignore',invalid='ignore'):
                good=long & (t['ndof']>0) & (t['chi2']/t['ndof']<25) & \
                    (t['n_layers']>=7) & (t['ndof']+5>=12) & (t['rmax_mm']<95)
            selection_masks[mode]=dict(all=long,reco_p_window=long&(momentum>50)&(momentum<250),
                                       good95=good)
            for name,mask in selection_masks[mode].items():
                counts[mode+'_'+name]+=int(mask.sum())
                for i,plane in enumerate(['st1','VetoNu']):
                    values=slope(t['projection'][:,i,3:6])
                    values=np.where(t['projection_valid'][:,i]==1,values,np.nan)
                    h=histogram(values[mask],edges)
                    tag='|'.join([mode,plane,name])
                    if tag not in histograms:histograms[tag]=h
                    else:
                        target=histograms[tag]
                        for field in ['counts','sumw','sumw2']:
                            target[field]=(np.asarray(target[field])+h[field]).tolist()
                        for field in ['raw_entries','finite','invalid','underflow','overflow']:
                            target[field]+=h[field]
        for name in selection_masks['off']:
            counts['membership_changes_'+name]+=int(np.count_nonzero(
                selection_masks['off'][name]!=selection_masks['on'][name]))
    result=dict(variant=args.variant,inputs=len(inputs['on']),counts=dict(counts),
        off_root=str(args.off_root),on_root=str(args.on_root),
        fields_checked=sorted(fields_checked),changed_elements=dict(changed),
        maximum_absolute_changes=maximum,examples=examples,histograms=histograms,
        same_identity_and_per_field_comparison=True,
        all_plotted_bins_equal=all(histograms['off|'+p+'|'+s]['counts']==histograms['on|'+p+'|'+s]['counts']
            for p in ['st1','VetoNu'] for s in ['all','reco_p_window','good95']))
    args.output.mkdir(parents=True,exist_ok=False)
    write_json(args.output/'summary.json',result)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for row,selection in enumerate(['all','reco_p_window']):
        for col,plane in enumerate(['st1','VetoNu']):
            ax=axes[row,col]
            for mode,color,style in [('off','#c62828','-'),('on','#1565c0','--')]:
                h=histograms['|'.join([mode,plane,selection])]
                ax.stairs(h['counts'],h['edges'],label=f"Isolated {mode.upper()} (N={h['raw_entries']:,})",
                          color=color,linestyle=style)
            ax.set(title=plane+' | '+('Full spectrum' if row==0 else '50 < total reco p < 250 GeV'),
                   xlabel=r'$s_y=p_y^F/p_z^F$',ylabel='Entries / bin')
            ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle(f"DSID 100044 | {args.variant} | {counts['events']:,} source events | Native Calypso")
    fig.tight_layout(rect=(0,0,1,.95))
    for extension in ['png','pdf']:fig.savefig(args.output/('isolated_vs_shared_angles.'+extension),dpi=180)
    plt.close(fig)
    write_json(args.output/'figure_index.json',[dict(name='isolated_vs_shared_angles',
        reference='geometryRotation.pdf pp19-20; full spectrum and 50-250 GeV',
        bin_counts=str(args.output/'summary.json'),files=['isolated_vs_shared_angles.png','isolated_vs_shared_angles.pdf'])])
    print(json.dumps({k:v for k,v in result.items() if k not in ['histograms','fields_checked','examples']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--off-root',type=Path,required=True)
    parser.add_argument('--on-root',type=Path,required=True)
    parser.add_argument('--variant',required=True)
    parser.add_argument('--output',type=Path,required=True)
    compare(parser.parse_args())
