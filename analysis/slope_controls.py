#!/usr/bin/env python3
"""Assess native Calypso direction scans and full-geometry G4 gun controls."""
import argparse
import json
from pathlib import Path
import shutil
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from PIL import Image

from slope_tables import histogram, join_truth, slope
from slope_statistics import distribution_structure, event_ratio_interval


COLORS = {'beforeGeo':'#c62828','afterGeo':'#1565c0'}


def archive(path):
    with np.load(path,allow_pickle=False) as file:
        return {k:file[k] for k in file.files}


def save(fig,output,name,report):
    fig.tight_layout()
    # Close seek-heavy image/PDF writes locally before copying to EOS.
    with tempfile.TemporaryDirectory(prefix='slope-figure-') as directory:
        local = Path(directory)
        fig.savefig(local/(name+'.png'),dpi=160)
        fig.savefig(local/(name+'.pdf'))
        for suffix in ['.png', '.pdf']:
            shutil.copyfile(local/(name+suffix), output/(name+suffix))
        with Image.open(output/(name+'.png')) as image:
            image.verify()
    report.savefig(fig)
    plt.close(fig)


def direction(root,output,report):
    data = {v:archive(root/v/'transport.npz') for v in COLORS}
    result = {'input_states_per_variant':294,'targets_z_F_mm':[100.,-100.],
              'definition':'Same physical free mean and angular covariance, transported by actual native Calypso',
              'before_pole_covariance':'Undefined phi chart: exact-pole means propagated without covariance'}
    b,a = data['beforeGeo'],data['afterGeo']
    assert np.array_equal(b['candidate'],a['candidate']) and len(b['candidate'])==294
    result['input_momentum_max_difference_GeV'] = float(np.max(np.abs(b['input_momentum_F_GeV']-a['input_momentum_F_GeV'])))
    result['projection_position_max_difference_mm'] = float(np.max(np.abs(b['projection_F'][:,:2,:3]-a['projection_F'][:,:2,:3])))
    result['projection_momentum_max_difference_GeV'] = float(np.max(np.abs(b['projection_F'][:,:2,3:]-a['projection_F'][:,:2,3:])))
    variance = {}
    fig,axes = plt.subplots(2,2,figsize=(11,8))
    for v,d in data.items():
        momentum = d['projection_F'][:,:2,3:]
        phi = np.arctan2(momentum[:,:,1],momentum[:,:,0])
        theta = np.arctan2(np.hypot(momentum[:,:,0],momentum[:,:,1]),momentum[:,:,2])
        g = np.zeros((len(phi),2,2))
        if v=='afterGeo':
            # phi is computed in physical F above; recover native phi from sy.
            native_phi = np.arctan2(momentum[:,:,1],momentum[:,:,2])
            g[:,:,0] = 1/np.cos(native_phi)**2
        else:
            g[:,:,0] = np.cos(phi)*np.tan(theta)
            g[:,:,1] = np.sin(phi)/np.cos(theta)**2
        c = d['projection_covariance_native'][:,:2,2:4,2:4]
        var = np.einsum('...i,...ij,...j->...',g,c,g)
        variance[v] = var
        valid = d['projection_valid'][:,:2]==1
        finite = np.isfinite(var) & (d['input_has_covariance'][:,None]==1)
        result[v] = dict(successful_projections=int(valid.sum()),failed=int((~valid).sum()),
            inputs_with_covariance=int(d['input_has_covariance'].sum()),finite_slope_variances=int(finite.sum()),
            negative_slope_variances=int(np.count_nonzero(finite & (var<0))),
            physical_sy_max_error_between_input_and_projection=float(np.max(np.abs(slope(momentum)-slope(d['input_momentum_F_GeV'])[:,None]))))
        radial = np.hypot(d['input_momentum_F_GeV'][:,0],d['input_momentum_F_GeV'][:,1])/d['input_momentum_F_GeV'][:,2]
        axes[0,0].scatter(radial,d['input_angles_native'][:,0],s=9,label=v,color=COLORS[v])
        axes[0,1].scatter(radial,d['input_angles_native'][:,1],s=9,label=v,color=COLORS[v])
        axes[1,0].scatter(radial,np.abs(slope(momentum[:,0])-slope(d['input_momentum_F_GeV'])),s=9,label=v,color=COLORS[v])
    common = np.isfinite(variance['beforeGeo']) & np.isfinite(variance['afterGeo']) & (variance['afterGeo']>0) & (b['input_has_covariance'][:,None]==1)
    ratio = variance['beforeGeo'][common]/variance['afterGeo'][common]
    result['physical_slope_variance_before_over_after_quantiles'] = np.quantile(ratio,[0,.01,.5,.99,1]).tolist()
    radial = np.hypot(b['input_momentum_F_GeV'][:,0],b['input_momentum_F_GeV'][:,1])/b['input_momentum_F_GeV'][:,2]
    axes[1,1].scatter(np.broadcast_to(radial[:,None],common.shape)[common],ratio,s=9)
    labels = [('Native phi [rad]','Near-pole azimuth chart'),('Native theta [rad]','Near-pole polar chart'),
              ('Absolute change in physical sy','Native field propagation: z = 0 to +100 mm'),('Var(sy) before / after','Equivalent physical input covariance')]
    for ax,(ylabel,title) in zip(axes.flat,labels):
        ax.set(xlabel='Input radial slope',ylabel=ylabel,title=title,xscale='symlog')
        ax.set_xscale('symlog',linthresh=1e-12)
        ax.grid(alpha=.2)
    for ax in axes.flat[:3]:ax.legend()
    fig.suptitle('Native Calypso direction scan | 294 states per geometry | p = 50, 100, 1000 GeV; q = +/-1')
    save(fig,output,'direction_scan',report)
    return result


def gun(root,output,report):
    manifest = json.loads((root/'manifest.json').read_text())
    results = {}
    for pdg in [-13,13]:
        for momentum in [100.,250.]:
            tasks = [t for t in manifest['sources'] if t['pdg']==pdg and t['p_GeV']==momentum]
            name = f'pdg{pdg}_p{int(momentum)}'
            case = {}
            truth_reference = None
            fig,axes = plt.subplots(2,2,figsize=(11,8),sharex='row')
            for variant in COLORS:
                rows = {'tracks':[],'truth':[],'events':[]}
                offset = 0
                for task in tasks:
                    base = Path(task['output'])
                    assert (base/'exit_code.txt').read_text().strip()=='0'
                    for table in rows:
                        d = archive(base/variant/(table+'.npz'))
                        d['sample_entry'] = d['source_entry']+offset
                        rows[table].append(d)
                    offset += task['events']
                # Hit-offset buffers are not used by this analysis.
                tables = {table:{k:np.concatenate([r[k] for r in parts]) for k in parts[0]
                    if k not in ['hit_offsets','hit_identifiers']} for table,parts in rows.items()}
                truth,tracks = tables['truth'],tables['tracks']
                reference = np.column_stack([truth['sample_entry'],truth['barcode'],truth['production_position_F_mm'],truth['production_momentum_F_GeV']])
                if truth_reference is None:truth_reference=reference
                else:assert np.array_equal(reference,truth_reference,equal_nan=True)
                index,matched = join_truth(tables)
                long = tracks['station_bits'] & 14 == 14
                signal = matched & (truth['pdg'][index]==pdg)
                selected = long & signal & (tracks['purity']>=.7)
                population=truth['pdg']==pdg
                accepted=population & (np.sum(truth['truth_hit_counts'][:,1:],axis=1)>=18)
                accepted &= np.all(np.linalg.norm(truth['truth_station_position_F_mm'][:,1:,:2],axis=2)<100,axis=1)
                efficient=np.bincount(index[selected],minlength=len(population))>0
                reco = np.stack([slope(tracks['projection'][:,i,3:6]) for i in [1,2]],axis=1)
                valid = tracks['projection_valid'][:,[1,2]]==1
                reco = np.where(valid,reco,np.nan)
                true = np.column_stack([slope(truth['production_momentum_F_GeV']),slope(truth['truth_sensor_direction_F'])])
                c = dict(source_events=offset,primary_truth=len(truth['barcode']),native_tracks=len(long),long_tracks=int(long.sum()),
                    matched_pure_long_tracks=int(selected.sum()),truth_plane_z_F_mm=-3112.,truth_identical_before_after=True,
                    efficiency_all_truth=event_ratio_interval(truth['sample_entry'],efficient&population,population,(0,offset)),
                    efficiency_accepted_H18_truth=event_ratio_interval(truth['sample_entry'],efficient&accepted,accepted,(0,offset)),
                    upstream_p_quantiles_GeV=np.quantile(np.linalg.norm(tracks['upstream_F'][selected,3:],axis=1),[.16,.5,.84]).tolist(),
                    wrong_charge=event_ratio_interval(tracks['sample_entry'],selected & (tracks['charge']!=np.sign(-pdg)),selected,(0,offset)),
                    momentum_tail=event_ratio_interval(tracks['sample_entry'],selected & (np.linalg.norm(tracks['upstream_F'][:,3:],axis=1)<.5*momentum),selected,(0,offset)))
                for i,plane in enumerate(['VetoNu','sensor']):
                    values = reco[selected,i]
                    residual = values-true[index[selected],i]
                    h=histogram(values,np.linspace(-.01,.01,101))
                    t=histogram(true[:,i],np.linspace(-.01,.01,101))
                    r=histogram(residual,np.linspace(-.01,.01,161) if plane=='VetoNu' else np.linspace(-.0003,.0003,151))
                    finite=residual[np.isfinite(residual)]
                    ids=tracks['sample_entry'][selected][np.isfinite(residual)]
                    rng=np.random.default_rng(50541)
                    widths=[]
                    for _ in range(500):
                        weights=rng.multinomial(offset,np.full(offset,1./offset))
                        draw=np.repeat(finite,weights[ids])
                        q=np.quantile(draw,[.16,.84]);widths.append(float((q[1]-q[0])/2))
                    quantiles=np.quantile(finite,[.01,.16,.5,.84,.99]).tolist() if len(finite) else []
                    c[plane] = dict(reco=h,truth=t,residual=r,
                        residual_quantiles=quantiles,
                        central68_halfwidth=float((quantiles[3]-quantiles[1])/2),
                        central68_halfwidth_interval_95=np.quantile(widths,[.025,.975]).tolist(),
                        central68_bootstrap_widths=widths,bootstrap_unit='source event; 500 replicates',
                        structure=distribution_structure(values),failed=int(np.count_nonzero(selected & ~valid[:,i])))
                    axes[0,i].stairs(h['counts'],h['edges'],label=variant+f"; N={h['finite']}",color=COLORS[variant])
                    if variant=='afterGeo':axes[0,i].stairs(t['counts'],t['edges'],label='Same-plane truth',color='#222222',linestyle='--')
                    axes[1,i].stairs(r['counts'],r['edges'],label=variant,color=COLORS[variant])
                    axes[0,i].set(title=plane,xlabel='Physical sy',ylabel='Entries / bin')
                    axes[1,i].set(xlabel='Reco sy - same-plane truth sy',ylabel='Entries / bin')
                case[variant]=c
            for ax in axes.flat:ax.legend();ax.grid(alpha=.2)
            fig.suptitle(f'Full Calypso G4 + digitization + native tracking | PDG {pdg}, p = {momentum:g} GeV | {offset:,} source events\nUniform launch angles +/-5 mrad and x/y +/-30 mm; exact launch truth at VetoNu')
            save(fig,output,'gun_'+name,report)
            results[name]=case
    results['width_scaling']={}
    for pdg in [-13,13]:
        for variant in COLORS:
            for plane in ['VetoNu','sensor']:
                a=results[f'pdg{pdg}_p100'][variant][plane]
                b=results[f'pdg{pdg}_p250'][variant][plane]
                ratio=np.array(a['central68_bootstrap_widths'])/np.array(b['central68_bootstrap_widths'])
                results['width_scaling'][f'{pdg}|{variant}|{plane}']=dict(
                    sigma100_over_sigma250=a['central68_halfwidth']/b['central68_halfwidth'],
                    interval_95=np.quantile(ratio,[.025,.975]).tolist(),expected_multiple_scattering_ratio=2.5)
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode',choices=['direction','gun'])
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    with PdfPages(args.output/'controls_report.pdf') as report:
        result = globals()[args.mode](args.root,args.output,report)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(args.output/'summary.json')


if __name__ == '__main__':main()
