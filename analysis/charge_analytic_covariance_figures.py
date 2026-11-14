#!/usr/bin/env python3
"""Plot the isolated covariance-map precision test with auditable paired inputs."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from charge_analytic_covariance import key, read_rows
from charge_component_refit import metrics, sign
from slope_tables import write_json


SAMPLES = ('100044', '100116', '100117')
SHORT = ('32 / 256', '32 / 1024', '64 / 1024', '64 / 4096')
MODES = ('finite_difference', 'analytic')


def values(source):
    data = json.loads((source/'summary.json').read_text())
    validation = json.loads((source/'validation.json').read_text())
    seed_rows = read_rows(source/'seed_input_metrics.jsonl.gz')
    controls = json.loads((source/'after_control_paired_components.json').read_text())
    old = {key(r): r for r in read_rows(Path(validation['reference_analysis'])/'paired_state_metrics.jsonl.gz')}
    new = {key(r): r for r in read_rows(Path(validation['analytic_analysis'])/'paired_state_metrics.jsonl.gz')}
    paired = defaultdict(list)
    for k in old.keys() & new.keys():
        if old[k]['same_ordered_CKF_hits'] and new[k]['same_ordered_CKF_hits']:
            paired[(k[0], k[1])].append((old[k], new[k]))
    groups, table = {}, []
    series_list = list(dict.fromkeys(k.split('|')[0] for k in data['groups']))
    for series in series_list:
        for dsid in SAMPLES:
            name = series+'|'+dsid
            group = data['groups'][name]
            seeds = [r for r in seed_rows if (r['series'],r['dsid']) == (series,dsid)]
            assert len(seeds) == group['input_components']
            for mode in MODES:
                for field, frozen in group['seed_input_metrics'][mode].items():
                    assert metrics([r[mode][field] for r in seeds]) == frozen
            cohort = group['common_state_cohort']
            states = paired[(series,dsid)]
            assert len(states) == cohort['states']
            assert len({(r[0]['entry'],r[0]['component']) for r in states}) == cohort['components']
            for stage, fields in cohort['metrics'].items():
                for field, modes in fields.items():
                    for i, mode in enumerate(MODES):
                        assert metrics([r[i][stage][field] for r in states]) == modes[mode]
            rows = [r for r in controls if (r['series'],r['dsid']) == (series,dsid)]
            c = group['after_only_control']['counts']
            assert len(rows) == c['components']
            common = [r for r in rows if r['before_selected'] and r['after_selected']]
            for stage in ['selected','converted','final']:
                assert sum(sign(r['before_'+stage+'_qop_GeVinv']) != sign(r['after_'+stage+'_qop_GeVinv'])
                           for r in common) == c[stage+'_sign_migrations']
            assert sum(r['same_selected_hit_sequence'] for r in common) == c['same_selected_hit_sequence']
            vc = group['version_counts']
            table.append([SHORT[series_list.index(series)],dsid,str(c['components']),
                f'{c["before_selected"]}/{c["after_selected"]}',
                f'{len(common)-c["same_selected_hit_sequence"]}/{c["final_sign_migrations"]}',
                f'{vc["finite_difference"]["final_sign_migrations"]}/{vc["analytic"]["final_sign_migrations"]}',
                f'{vc["analytic"]["before_refit_success"]}/{vc["analytic"]["after_refit_success"]}'])
            groups[name] = dict(seed=group['seed_input_metrics'],state=cohort['metrics'],
                                 states=cohort['states'],components=cohort['components'])
    return dict(series=series_list,groups=groups,table=table)


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D
    data = values(source)
    output.mkdir(parents=True, exist_ok=False)
    index = []

    def save(fig, name, title, caption, pdf):
        fig.suptitle(title,fontsize=14,y=.985)
        fig.text(.5,.044,'Actual Calypso | Isolated smoothing ON | Before regular CKF + refit; after native | Refit diagonal scale 1',ha='center',fontsize=8)
        fig.text(.5,.022,caption,ha='center',fontsize=8)
        fig.savefig(output/(name+'.png'),dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name,png=name+'.png',pdf=name+'.pdf',reference='S06 actual Calypso seed covariance Jacobian control'))

    with PdfPages(output/'analytic_covariance_report.pdf') as pdf:
        fig,ax=plt.subplots(figsize=(15,8))
        ax.axis('off')
        table=ax.table(cellText=data['table'],colLabels=['n / C divisor','DSID','Inputs',
            'After selected\nFD / analytic','After changes\nhits / final sign',
            'Version final sign\nFD / analytic','Analytic refit\nsuccess B / A'],cellLoc='center',loc='center')
        table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,2.2)
        for (row,col),cell in table.get_celld().items():
            if row==0:cell.set_facecolor('#dce6ef');cell.set_height(cell.get_height()*1.4)
            elif (row-1)//3%2==0:cell.set_facecolor('#f2f5f8')
        fig.tight_layout(rect=[0,.09,1,.94])
        save(fig,'analytic_covariance_outcomes','Changing only the after-geometry seed covariance differential',
            '60 distinct source events | Four correlated finite-prior settings | FD: historical central difference | Counts are not production rates',pdf)

        for family, name, title, rows in [
            ('seed','analytic_covariance_seed','Physical covariance agreement at the unchanged seed centers',
             [('physical_covariance_scaled_max','Seed full physical-C discrepancy'),
              ('direction_covariance_scaled_max','Seed direction-C discrepancy')]),
            ('state','analytic_covariance_states','Native smoothed-state residuals on the common ordered-hit cohort',
             [('direction_covariance_scaled_max','Smoothed direction-C discrepancy'),
              ('qop_difference_GeVinv',r'Smoothed absolute $q/p$ difference [GeV$^{-1}$]')])]:
            fig,axes=plt.subplots(2,3,figsize=(16,10))
            for i,(field,label) in enumerate(rows):
                for j,dsid in enumerate(SAMPLES):
                    ax=axes[i,j]
                    for mode,style in [('finite_difference','--'),('analytic','-')]:
                        for stat,marker,color in [('median','o','#326ba4'),('p95','^','#c68d27'),('max','x','#b64d48')]:
                            measured=[]
                            for series in data['series']:
                                g=data['groups'][series+'|'+dsid][family]
                                m=(g[mode][field] if family=='seed' else g['smoothed'][field][mode])
                                measured.append(m[stat])
                            y=np.asarray(measured)
                            ax.plot(np.arange(4),np.where(y>0,y,np.nan),linestyle=style,marker=marker,color=color,linewidth=1.2,markersize=5)
                    ax.set_yscale('log');ax.set_xticks(np.arange(4),SHORT,fontsize=8)
                    ax.set_xlabel('Number of components / within-C divisor')
                    ax.set_title(dsid);ax.set_ylabel(label);ax.grid(alpha=.2)
            handles=[Line2D([0],[0],color=c,marker=m,label=s) for s,m,c in
                     [('median','o','#326ba4'),('p95','^','#c68d27'),('maximum','x','#b64d48')]]
            handles += [Line2D([0],[0],color='black',linestyle=style,label=label)
                        for style,label in [('--','Finite difference'),('-','Analytic')]]
            fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.95),ncol=5,fontsize=9)
            fig.tight_layout(rect=[0,.09,1,.89])
            save(fig,name,title,'Common FASER frame | All maxima retained | Zero differences omitted on log axes | Finite Gaussian posterior convergence is not tested',pdf)
    write_json(output/'plot_values.json',data)
    write_json(output/'figure_index.json',index)


def validate(source, output, pdf_tools):
    import re
    import sys
    from PIL import Image
    sys.path.insert(0,str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics
    assert json.loads((output/'plot_values.json').read_text()) == values(source)
    def normalize(page):
        fonts=page['/Resources']['/Font']
        return re.sub(r'/F\d+(?=\s|-)',lambda m:re.sub(r'^/[A-Z]{6}\+','/',
            str(fonts[m[0]].get_object()['/BaseFont'])),normalized_graphics(page))
    index=json.loads((output/'figure_index.json').read_text())
    pages=PdfReader(output/'analytic_covariance_report.pdf').pages
    assert len(index)==len(pages)==3
    for row,page in zip(index,pages):
        with Image.open(output/row['png']) as png:png.load()
        original=PdfReader(output/row['pdf']).pages[0]
        assert normalize(original)==normalize(page)
        assert (list(original.mediabox),list(original.cropbox),original.rotation)==(list(page.mediabox),list(page.cropbox),page.rotation)
    write_json(output/'validation.json',dict(pages=3,all_seed_and_state_metrics_recomputed=True,
        after_control_counts_recomputed=True,plotted_values_exact=True,PNG_decoded=True,
        normalized_PDF_content_and_boxes_equal=True,visual_inspection='Pending original-resolution inspection'))


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
