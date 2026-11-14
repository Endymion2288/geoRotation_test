#!/usr/bin/env python3
"""Show native backward covariance operations and their exact replay gates."""
import argparse
import json
from pathlib import Path

import numpy as np

from charge_component_refit import audit, metrics
from slope_tables import write_json


SAMPLES=('100044','100116','100117')


def values(source):
    summary=json.loads((source/'summary.json').read_text())
    rows=list(audit(source/'paired_operations.jsonl.gz'))
    groups={}
    for dsid in SAMPLES:
        group=[r for r in rows if r['dsid']==dsid]
        frozen=summary['groups'][dsid]
        assert len(group)==frozen['counts']['paired_backward_operations']
        assert len({(r['entry'],r['component']) for r in group})==frozen['counts']['paired_paths_with_operations']
        assert frozen['counts']['paired_paths_with_operations']+frozen['counts'].get('terminal_only_paths',0)==frozen['counts']['common_selected_components']
        for field,expected in frozen['fields'].items():
            assert metrics([r[field]['direction'] for r in group])==expected
        for variant in ('before','after'):
            assert sum(not r[variant]['psd_correction_exactly_zero'] for r in group)==frozen['counts'][variant+'_corrected_C']
        groups[dsid]=dict(counts=frozen['counts'],fields=frozen['fields'],
            raw=[r['raw_smoothed_C_difference']['direction'] for r in group],
            validated=[r['validated_smoothed_C_difference']['direction'] for r in group])
    extreme=max(rows,key=lambda r:r['validated_smoothed_C_difference']['direction'])
    path=[r for r in rows if (r['dsid'],r['entry'],r['component'])==(
        extreme['dsid'],extreme['entry'],extreme['component'])]
    path.sort(key=lambda r:r['operation'])
    return dict(groups=groups,extreme=extreme,path=path)


def plot(source,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    data=values(source)
    output.mkdir(parents=True,exist_ok=False)
    index=[]

    def save(fig,name,title,pdf):
        fig.suptitle(title,fontsize=14,y=.985)
        fig.text(.5,.044,'Actual Calypso / ACTS 32.0.2 | Isolated smoothing ON | n32, within-C / 256 | Analytic seed coordinate Jacobian',ha='center',fontsize=9)
        fig.text(.5,.022,'60 distinct source events | Correlated components, first-ranked common ordered-hit paths | Counts are not production rates',ha='center',fontsize=9)
        fig.savefig(output/(name+'.png'),dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name,png=name+'.png',pdf=name+'.pdf',reference='S06 actual Calypso backward-operation audit; production reference plots unchanged by exact replay'))

    with PdfPages(output/'smoothing_operations_report.pdf') as pdf:
        fig,ax=plt.subplots(figsize=(15,7))
        ax.axis('off')
        table=[]
        for dsid in SAMPLES:
            g=data['groups'][dsid];c=g['counts']
            table.append([dsid,str(c['common_selected_components']),str(c['paired_backward_operations']),
                f'{c["before_corrected_C"]} / {c["after_corrected_C"]}',
                f'{g["fields"]["raw_smoothed_C_difference"]["max"]:.3e}',
                f'{g["fields"]["validated_smoothed_C_difference"]["max"]:.3e}'])
        cells=ax.table(cellText=table,colLabels=['DSID','Common selected\ncomponents','Paired backward\noperations',
            'PSD-corrected C\nbefore / after','Maximum raw\ndirection-C difference','Maximum validated\ndirection-C difference'],cellLoc='center',loc='center')
        cells.auto_set_font_size(False);cells.set_fontsize(10);cells.scale(1,3)
        for (row,_),cell in cells.get_celld().items():
            if row==0:cell.set_facecolor('#dce6ef')
        ax.text(.5,.82,'Default OFF and observer ON: all six complete native audit streams match the frozen reference exactly',ha='center',transform=ax.transAxes,fontsize=11)
        terminal=' | '.join(f'{dsid}: {data["groups"][dsid]["counts"].get("terminal_only_paths",0)}' for dsid in SAMPLES)
        ax.text(.5,.16,'C difference: maximum elementwise scaled difference in common regular (phi, latitude) coordinates\nTerminal-only selected candidates (zero backward steps): '+terminal+'\nBoth correction counts include all paired backward operations; no outliers removed',ha='center',transform=ax.transAxes,fontsize=10,linespacing=1.6)
        fig.tight_layout(rect=[0,.09,1,.94])
        save(fig,'smoothing_operations_counts','Read-only native smoother observation and paired-path coverage',pdf)

        fig,axes=plt.subplots(1,3,figsize=(16,6))
        for ax,dsid in zip(axes,SAMPLES):
            g=data['groups'][dsid]
            for field,label,color,style in [('raw','Before PSD validation','#326ba4','-'),
                                             ('validated','After PSD validation','#b64d48','--')]:
                all_values=np.asarray(g[field]);ordered=np.sort(all_values[all_values>0])
                ax.step(ordered,(np.arange(len(ordered))+1+(all_values==0).sum())/len(all_values),where='post',label=label,color=color,linestyle=style)
            ax.set_xscale('log');ax.set_ylim(0,1.02)
            ax.set_title(f'{dsid} | N = {len(g["raw"]):,} operations');ax.set_xlabel('Scaled direction-C difference')
            ax.set_ylabel('Fraction of paired operations');ax.grid(alpha=.2);ax.legend(fontsize=9,loc='lower right')
        fig.tight_layout(rect=[0,.10,1,.93])
        save(fig,'smoothing_operations_covariance','Before/after-geometry covariance differences through native PSD validation',pdf)

        path=data['path'];extreme=data['extreme']
        x=[r['operation'] for r in path]
        fig,axes=plt.subplots(2,2,figsize=(15,10))
        for field,label,color,style in [('filtered_C_difference','Current filtered C','#678d45','-'),
            ('next_predicted_C_difference','Next predicted C','#bd943b','-'),
            ('next_smoothed_C_difference','Inherited next smoothed C','#8165a7','-'),
            ('raw_smoothed_C_difference','Calculated raw C','#326ba4','-'),
            ('validated_smoothed_C_difference','Validated C','#b64d48','--')]:
            y=np.array([r[field]['direction'] for r in path])
            axes[0,0].plot(x,np.where(y>0,y,np.nan),label=label,color=color,linestyle=style,marker='.',linewidth=1.3)
        axes[0,0].set_ylabel('Scaled direction-C difference');axes[0,0].legend(fontsize=8)
        for field,label,color in [('filtered_mean_difference','Current filtered','#678d45'),
                                 ('next_smoothed_mean_difference','Inherited next smoothed','#8165a7'),
                                 ('calculated_smoothed_mean_difference','Calculated smoothed','#326ba4')]:
            for ax,name in [(axes[0,1],'direction_rad'),(axes[1,0],'qop_GeVinv')]:
                y=np.array([r[field][name] for r in path])
                ax.plot(x,np.where(y>0,y,np.nan),marker='.',label=label,color=color,linewidth=1.3)
        axes[0,1].set_ylabel('Maximum angular difference [rad]');axes[0,1].legend(fontsize=8)
        axes[1,0].set_ylabel(r'Absolute $q/p$ difference [GeV$^{-1}$]');axes[1,0].legend(fontsize=8)
        for name,label,color in [('before','Before geometry','#326ba4'),('after','After geometry','#b64d48')]:
            axes[1,1].plot(x,[r[name]['predicted_C_condition']['correlation_scaled'] for r in path],marker='.',label=label,color=color,linewidth=1.3)
        axes[1,1].set_ylabel('Condition number of diagonal-scaled predicted C');axes[1,1].legend(fontsize=8)
        for ax in axes.flat:
            ax.set_yscale('log');ax.set_xlabel('Backward operation index (actual execution order)');ax.grid(alpha=.2)
            ax.axvline(extreme['operation'],color='gray',alpha=.4,linestyle=':')
        fig.tight_layout(rect=[0,.09,1,.94])
        save(fig,'smoothing_operations_extreme',f'Largest direction-C path: DSID {extreme["dsid"]}, entry {extreme["entry"]}, component {extreme["component"]}',pdf)
    write_json(output/'figure_index.json',index)
    write_json(output/'plot_values.json',data)


def validate(source,output,pdf_tools):
    import re
    import sys
    from PIL import Image
    sys.path.insert(0,str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics
    assert json.loads((output/'plot_values.json').read_text())==values(source)
    def graphics(page):
        fonts=page['/Resources']['/Font']
        return re.sub(r'/F\d+(?=\s|-)',lambda m:re.sub(r'^/[A-Z]{6}\+','/',
            str(fonts[m[0]].get_object()['/BaseFont'])),normalized_graphics(page))
    index=json.loads((output/'figure_index.json').read_text())
    pages=PdfReader(output/'smoothing_operations_report.pdf').pages
    assert len(index)==len(pages)==3
    for row,page in zip(index,pages):
        with Image.open(output/row['png']) as png:png.load()
        original=PdfReader(output/row['pdf']).pages[0]
        assert graphics(original)==graphics(page)
        assert (list(original.mediabox),list(original.cropbox),original.rotation)==(list(page.mediabox),list(page.cropbox),page.rotation)
    write_json(output/'validation.json',dict(pages=3,operation_metrics_and_counts_recomputed=True,
        plotted_values_exact=True,PNG_decoded=True,normalized_PDF_content_and_boxes_equal=True,
        visual_inspection='Pending original-resolution inspection'))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--phase',choices=['plot','validate'],required=True)
    parser.add_argument('--pdf-tools',type=Path)
    args=parser.parse_args()
    if args.phase=='plot':plot(args.input,args.output)
    else:validate(args.input,args.output,args.pdf_tools)


if __name__=='__main__':
    main()
