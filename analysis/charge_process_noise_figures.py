#!/usr/bin/env python3
"""Plot frozen native Calypso charge outcomes and fixed-ID process controls."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

from slope_tables import write_json

COORDINATES = ['native', 'sdo_strips_fixed', 'centroid_fixed']
COORDINATE_LABELS = ['Native X', 'SDO strip X', 'SiHit centroid X']
PROCESSES = ['material_on', 'no_scattering', 'no_energy_loss', 'no_material']
PROCESS_LABELS = ['Material on', 'No scattering', 'No energy loss', 'No material']


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    events = json.loads((source/'events.json').read_text())
    resources = json.loads((source/'resources.json').read_text())
    assert summary['events'] == len(events) == 66 and summary['jobs'] == 792
    pure = []
    for row in events:
        if row['cohort'] != 'pure_wrong':
            continue
        current = dict(dsid=row['dsid'], sample_entry=row['sample_entry'], modes={})
        for mode, state in row['modes'].items():
            current['modes'][mode] = {key:state[key] for key in ['status', 'final_qop_GeVinv',
                'final_qop_variance_GeVinv2', 'same_ordered_hits', 'missing_native_hit_count',
                'native_reference_qop_GeVinv', 'failed_quality_conditions']}
            if row['dsid'] == '100116' and row['sample_entry'] == 51839:
                current['modes'][mode]['qop_trace'] = state.get('qop_trace', [])
        pure.append(current)
    assert len(pure) == 6
    return dict(summary=summary, pure=pure, resources={key:value for key,value in resources.items() if key != 'rows'})


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import ListedColormap
    data = values(source)
    output.mkdir(exist_ok=False)
    index = []
    modes = [c+'|'+p for c in COORDINATES for p in PROCESSES]
    labels = [c+'\n'+p for c in COORDINATE_LABELS for p in PROCESS_LABELS]

    def save(fig, name, title, pdf, top=.93):
        fig.suptitle(title, fontsize=13)
        fig.text(.5, .034, 'Actual geoRotation Calypso | Isolated smoothing ON | Complete original ancestor seed and admitted cluster IDs fixed', ha='center', fontsize=8)
        fig.text(.5, .014, '66 outcome-selected source events: 6 original pure wrong + 60 normal | Only CKF process terms change; refit and observer remain native', ha='center', fontsize=8)
        fig.tight_layout(rect=[0, .085, 1, top])
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
            reference='S06 actual Calypso fixed seed, admitted cluster IDs and CKF material-process control'))

    with PdfPages(output/'fixed_hit_process_report.pdf') as pdf:
        fig, axes = plt.subplots(2, 1, figsize=(15, 8), gridspec_kw={'height_ratios':[3, 1.4]})
        mapping = {'correct':0, 'wrong':1, 'quality_failed':2}
        codes = np.array([[mapping.get(row['modes'][mode]['status'],3) for mode in modes] for row in data['pure']])
        axes[0].imshow(codes, aspect='auto', cmap=ListedColormap(['#3a9257','#ca5555','#e2ae4d','#999999']), vmin=0, vmax=3)
        axes[0].set_yticks(np.arange(6), [r['dsid']+' / '+str(r['sample_entry']) for r in data['pure']], fontsize=8)
        axes[0].set_xticks(np.arange(12), labels, fontsize=7)
        axes[0].set_title('Six original pure wrong: C correct, W wrong, Q quality failure, O other; * same ordered hits')
        for i, row in enumerate(data['pure']):
            for j, mode in enumerate(modes):
                state = row['modes'][mode]
                label = {'correct':'C','wrong':'W','quality_failed':'Q'}.get(state['status'],'O')
                if state['same_ordered_hits']:
                    label += '*'
                axes[0].text(j, i, label, ha='center', va='center', fontsize=8)
        counts = []
        for mode in modes:
            counts.append(sum(data['summary']['groups']['|'.join([dsid,'normal',mode])].get('correct',0) for dsid in ['100116','100117']))
        axes[1].bar(np.arange(12), counts, color='#3a9257')
        axes[1].set_ylim(0,66)
        axes[1].set_xticks(np.arange(12), labels, fontsize=7)
        axes[1].set_ylabel('Normal selected correct / 60')
        axes[1].bar_label(axes[1].containers[0], fontsize=8)
        axes[1].grid(axis='y', alpha=.2)
        save(fig, 'fixed_hit_process_outcomes', 'Charge outcomes with one original ancestor seed and exact admitted cluster IDs', pdf)

        special = next(row for row in data['pure'] if row['dsid']=='100116' and row['sample_entry']==51839)
        fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
        for i, process in enumerate(['material_on','no_scattering']):
            for j, coordinate in enumerate(COORDINATES):
                ax = axes[i,j]
                state = special['modes'][coordinate+'|'+process]
                trace = state.get('qop_trace',[])
                for stage, color in [('predicted','#309342'),('filtered','#cf4545'),('smoothed','#3166b8')]:
                    ax.plot([s['z_F_mm'] for s in trace], [1000*s[stage] for s in trace], '.-',
                            color=color, markersize=3, linewidth=1, label=stage.capitalize())
                ax.axhline(0,color='gray',linewidth=.8)
                ax.set_title(COORDINATE_LABELS[j]+' | '+('Material on' if i==0 else 'CKF scattering off')+'\n'+state['status']+'; missing IDs='+str(state['missing_native_hit_count']),fontsize=9)
                ax.set_xlabel('FASER longitudinal z [mm]')
                ax.set_ylabel(r'$q/p$ [$10^{-3}$ GeV$^{-1}$]')
                ax.grid(alpha=.2)
        axes[0,0].legend(fontsize=8)
        save(fig, 'fixed_hit_51839_states', '100116 / entry 51839: actual predicted, filtered and isolated smoothed q/p', pdf)

        colors = ['#555555','#245fa5','#df8734']
        fig, axes = plt.subplots(3, 2, figsize=(14, 10))
        for ax, row in zip(axes.flat,data['pure']):
            charge = 1 if row['dsid']=='100116' else -1
            for j, coordinate in enumerate(COORDINATES):
                x = np.arange(4)+(j-1)*.14
                states = [row['modes'][coordinate+'|'+p] for p in PROCESSES]
                y = [1000*charge*s['final_qop_GeVinv'] if s['final_qop_GeVinv'] is not None else np.nan for s in states]
                error = [1000*np.sqrt(s['final_qop_variance_GeVinv2']) if s['final_qop_variance_GeVinv2'] is not None and s['final_qop_variance_GeVinv2']>=0 else np.nan for s in states]
                ax.errorbar(x,y,yerr=error,color=colors[j],fmt='o-',linewidth=.9,markersize=3,capsize=2,label=COORDINATE_LABELS[j])
                for position,value,state in zip(x,y,states):
                    if state['status']!='correct' and np.isfinite(value):
                        ax.scatter(position,value,color='black',marker='x',s=35,zorder=5)
            ax.axhline(0,color='gray',linewidth=.8)
            ax.set_xticks(np.arange(4),PROCESS_LABELS,fontsize=7)
            ax.set_ylabel(r'$q_{truth}(q/p)_{final}$ [$10^{-3}$ GeV$^{-1}$]')
            ax.set_title(row['dsid']+' / '+str(row['sample_entry'])+' | reported 1 sigma errors',fontsize=9)
            ax.grid(alpha=.2)
        axes[0,0].legend(fontsize=7)
        save(fig, 'fixed_hit_process_curvatures', 'Six original pure wrong: final signed q/p and reported covariance; crosses mark wrong sign or quality failure', pdf)

        fig, axes = plt.subplots(2, 2, figsize=(14, 8))
        for i, cohort in enumerate(['pure_wrong','normal']):
            for j, statistic in enumerate(['same_ordered_hits','final_correct_sign']):
                ax = axes[i,j]
                rows=[]
                for coordinate in COORDINATES:
                    rows.append([sum(data['summary']['groups']['|'.join([dsid,cohort,coordinate,p])].get(statistic,0) for dsid in ['100116','100117']) for p in PROCESSES])
                matrix=np.asarray(rows)
                ax.imshow(matrix,aspect='auto',cmap='Blues',vmin=0,vmax=6 if cohort=='pure_wrong' else 60)
                ax.set_xticks(np.arange(4),PROCESS_LABELS,fontsize=8)
                ax.set_yticks(np.arange(3),COORDINATE_LABELS,fontsize=8)
                ax.set_title(('Original pure wrong / 6' if cohort=='pure_wrong' else 'Normal / 60')+' | '+statistic.replace('_',' '),fontsize=10)
                for ii in range(3):
                    for jj in range(4):
                        ax.text(jj,ii,str(matrix[ii,jj]),ha='center',va='center',color='black',fontsize=10)
        save(fig, 'fixed_hit_process_paths', 'Ordered-path retention and final correct signs before the quality decision', pdf)
    write_json(output/'plot_values.json',data)
    write_json(output/'figure_index.json',index)


def validate(source,output,pdf_tools):
    from PIL import Image
    sys.path.insert(0,str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics as base_normalized_graphics
    import re

    def normalized_graphics(page):
        """Resolve page-local graphics, marker, image and font aliases."""
        value = base_normalized_graphics(page)
        resources = page['/Resources']
        images = resources.get('/XObject', {})
        def image(match):
            obj = images[match[1]].get_object()
            metadata = {key: value for key, value in obj.items()
                        if key not in ['/Length', '/Filter']}
            return json.dumps(metadata, sort_keys=True, default=str)+' '+obj.get_data().decode('latin1')+' Do'
        return re.sub(r'(/I\d+) Do', image, value)
    assert json.loads((output/'plot_values.json').read_text())==values(source)
    rows=json.loads((output/'figure_index.json').read_text())
    report=PdfReader(output/'fixed_hit_process_report.pdf')
    assert len(rows)==len(report.pages)==4
    for row,page in zip(rows,report.pages):
        with Image.open(output/row['png']) as image:
            image.load()
        single=PdfReader(output/row['pdf']).pages[0]
        assert normalized_graphics(single)==normalized_graphics(page)
        assert list(single.mediabox)==list(page.mediabox) and list(single.cropbox)==list(page.cropbox)
        assert single.rotation==page.rotation
    write_json(output/'validation.json',dict(pages=4,plotted_values_exact=True,PNG_decoded=True,
        normalized_vector_graphics_and_boxes_equal=True,visual_inspection='Pending separate full-size review'))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['plot','validate'],required=True)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--pdf-tools',type=Path)
    args=parser.parse_args()
    if args.phase=='plot':
        plot(args.input,args.output)
    else:
        validate(args.input,args.output,args.pdf_tools)


if __name__=='__main__':
    main()
