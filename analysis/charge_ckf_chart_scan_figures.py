#!/usr/bin/env python3
"""Plot native finite-prior chart checks and the first residual growth stage."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path

import numpy as np

from charge_ckf_chart_scan import FIELDS, STAGES
from charge_component_refit import metrics, sign
from slope_tables import write_json


SAMPLES = ('100044', '100116', '100117')
SHORT = ('32 / 256', '32 / 1024', '64 / 1024', '64 / 4096')


def compressed_rows(path):
    with gzip.open(path, 'rt') as stream:
        return [json.loads(line) for line in stream]


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    components = json.loads((source/'paired_components.json').read_text())
    states = compressed_rows(source/'paired_state_metrics.jsonl.gz')
    growth = compressed_rows(source/'first_growth_native_states.jsonl.gz')
    counts = Counter((r['series'], r['dsid'], r['field'], r['first_stage']) for r in growth)
    tables, graph = [], {}
    for series in summary['series']:
        for dsid in SAMPLES:
            group = summary['groups'][series+'|'+dsid]
            rows = [r for r in components if r['series']==series and r['dsid']==dsid]
            c = group['counts']
            assert len(rows) == c['components']
            for side in ('before', 'after'):
                assert sum(r[side+'_selected'] for r in rows) == c[side+'_selected']
                assert sum(r[side+'_refit_success'] is True for r in rows) == c[side+'_refit_success']
            common = [r for r in rows if r['before_selected'] and r['after_selected']]
            for stage in ('selected', 'converted', 'final'):
                assert sum(sign(r['before_'+stage+'_qop_GeVinv']) != sign(r['after_'+stage+'_qop_GeVinv'])
                           for r in common) == c[stage+'_sign_migrations']
            available = [r for r in states if r['series']==series and r['dsid']==dsid and r['same_ordered_CKF_hits']]
            frozen = group['physical_state_metrics']
            assert len(available)==frozen['n_states']
            assert len({(r['entry'],r['component']) for r in available})==frozen['components']
            for stage in STAGES:
                for field in FIELDS:
                    assert metrics([r[stage][field] for r in available])==frozen['metrics'][stage][field]
            tables.append([SHORT[summary['series'].index(series)],dsid,str(c['components']),
                f'{c["before_selected"]}/{c["after_selected"]}',str(c['same_selected_hit_sequence']),
                f'{c["selected_sign_migrations"]}/{c["final_sign_migrations"]}',
                f'{c["before_refit_success"]}/{c["after_refit_success"]}'])
            graph[series+'|'+dsid] = dict(metrics=frozen['metrics'], components=frozen['components'],
                growth={field:[counts[(series,dsid,field,s)] for s in STAGES]
                        for field in summary['diagnostic_thresholds']})
    return dict(series=summary['series'], tables=tables, graph=graph,
                diagnostic_thresholds=summary['diagnostic_thresholds'])


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    data = values(source)
    output.mkdir(parents=True, exist_ok=False)
    index = []

    def save(fig, name, title, caption, pdf):
        fig.suptitle(title,fontsize=14)
        fig.text(.5,.04,'Actual Calypso | Before regular CKF+refit vs after native | Isolated smoothing ON | Refit diagonal scale 1',ha='center',fontsize=8)
        fig.text(.5,.018,caption,ha='center',fontsize=8)
        fig.tight_layout(rect=[0,.09,1,.94])
        fig.savefig(output/(name+'.png'),dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name,png=name+'.png',pdf=name+'.pdf',reference='S06 actual native finite-prior chart width scan'))

    with PdfPages(output/'ckf_chart_scan_report.pdf') as pdf:
        fig,ax=plt.subplots(figsize=(14,8))
        ax.axis('off')
        table=ax.table(cellText=data['tables'],colLabels=['n / C divisor','DSID','Inputs','Selected B/A',
            'Same ordered hits','Sign differences CKF/final','Refit success B/A'],cellLoc='center',loc='center')
        table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,2.1)
        for (row,column),cell in table.get_celld().items():
            if row==0:cell.set_facecolor('#dce6ef')
            elif (row-1)//3%2==0:cell.set_facecolor('#f3f5f7')
        save(fig,'ckf_chart_scan_coverage','Finite-prior width scan: exact coverage, paths, charge and refit counts',
            '60 distinct source events reused in four settings | Components are correlated diagnostics; counts are not production misidentification rates',pdf)

        fig,axes=plt.subplots(2,3,figsize=(15,9))
        for j,dsid in enumerate(SAMPLES):
            for i,(field,label) in enumerate([('direction_covariance_scaled_max','Smoothed direction-C discrepancy'),
                    ('qop_difference_GeVinv',r'Smoothed absolute $q/p$ difference [GeV$^{-1}$]')]):
                ax=axes[i,j]
                measured=[data['graph'][s+'|'+dsid]['metrics']['smoothed'][field] for s in data['series']]
                for key,marker,color in [('median','o','#326ba4'),('p95','^','#d29a31'),('max','x','#b64d48')]:
                    y=np.array([r[key] for r in measured],dtype=float)
                    ax.plot(np.arange(4),np.where(y>0,y,np.nan),marker+'-',color=color,label=key)
                ax.set_yscale('log');ax.set_xticks(np.arange(4),SHORT,fontsize=8)
                ax.set_xlabel('Number of components / within-C divisor')
                ax.set_ylabel(label);ax.set_title(dsid);ax.grid(alpha=.2)
                if j==0:ax.legend(fontsize=8)
        save(fig,'ckf_chart_scan_residuals','Physical state differences across four finite-prior settings',
            'Same ordered CKF paths within each before/after pair | Each setting retains its own complete cohort | Direction-C metric is dimensionless in FASER frame',pdf)

        fig,axes=plt.subplots(3,3,figsize=(15,11))
        labels={'direction_difference_norm':'Direction norm','qop_difference_GeVinv':'q/p [GeV^-1]',
                'direction_covariance_scaled_max':'Direction-C metric'}
        for i,(field,threshold) in enumerate(data['diagnostic_thresholds'].items()):
            for j,dsid in enumerate(SAMPLES):
                ax=axes[i,j];bottom=np.zeros(4)
                for k,(stage,color) in enumerate(zip(STAGES,['#326ba4','#d29a31','#b64d48'])):
                    rates=np.array([100*data['graph'][s+'|'+dsid]['growth'][field][k]/data['graph'][s+'|'+dsid]['components'] for s in data['series']])
                    ax.bar(np.arange(4),rates,bottom=bottom,color=color,label=stage)
                    bottom+=rates
                ax.set_xticks(np.arange(4),SHORT,fontsize=8)
                ax.set_ylim(0,max(1.,1.35*max(bottom)))
                ax.set_ylabel('Same-hit diagnostic components [%]')
                ax.set_title(f'{dsid} | {labels[field]} > {threshold:g}',fontsize=9)
                ax.grid(axis='y',alpha=.2)
                if i==0 and j==0:ax.legend(fontsize=7)
        save(fig,'ckf_chart_scan_first_growth','First measurement location along the path with a residual above a descriptive threshold',
            'P/F/S are tested at that location; this is not causal execution order | Thresholds locate residuals; all states enter the numerical summaries',pdf)
    write_json(output/'plot_values.json',data)
    write_json(output/'figure_index.json',index)


def validate(source, output, pdf_tools):
    import re
    import sys
    from PIL import Image
    sys.path.insert(0,str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics
    assert json.loads((output/'plot_values.json').read_text())==values(source)
    def normalize(page):
        fonts=page['/Resources']['/Font']
        return re.sub(r'/F\d+(?=\s|-)',lambda m:re.sub(r'^/[A-Z]{6}\+','/',
            str(fonts[m[0]].get_object()['/BaseFont'])),normalized_graphics(page))
    index=json.loads((output/'figure_index.json').read_text())
    pages=PdfReader(output/'ckf_chart_scan_report.pdf').pages
    assert len(index)==len(pages)==3
    for row,page in zip(index,pages):
        with Image.open(output/row['png']) as png:png.load()
        original=PdfReader(output/row['pdf']).pages[0]
        assert normalize(original)==normalize(page)
        assert list(original.mediabox)==list(page.mediabox) and list(original.cropbox)==list(page.cropbox)
        assert original.rotation==page.rotation
    write_json(output/'validation.json',dict(pages=3,native_counts_and_all_state_metrics_recomputed=True,
        plotted_values_exact=True,PNG_decoded=True,normalized_PDF_content_and_boxes_equal=True,
        visual_inspection='Pending original-resolution inspection'))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['plot','validate'],required=True)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--pdf-tools',type=Path)
    args=parser.parse_args()
    (plot(args.input,args.output) if args.phase=='plot' else validate(args.input,args.output,args.pdf_tools))


if __name__=='__main__':
    main()
