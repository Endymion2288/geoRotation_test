#!/usr/bin/env python3
"""Make compact English figures for the finite common-prior CKF control."""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    rows = json.loads((args.input/'components.json').read_text())
    pairs = json.loads((args.input/'paired.json').read_text())
    manifest = json.loads((args.input.parent/'manifest.json').read_text())
    components = manifest['n_components']
    fraction = manifest['within_component_covariance_fraction']
    args.output.mkdir(parents=True, exist_ok=False)
    index=[]
    def save(fig,name):
        fig.savefig(args.output/(name+'.png'),dpi=160)
        fig.savefig(args.output/(name+'.pdf'))
        index.append(dict(name=name,png=name+'.png',pdf=name+'.pdf',
            reference='Actual Calypso finite common-prior CKF; isolated smoothing ON'))
    with PdfPages(args.output/'common_prior_report.pdf') as pdf:
        fig,axes=plt.subplots(1,2,figsize=(13,5))
        labels=['100044 μ+','100116 μ+','100117 μ−']
        keys=['100044','100116','100117']
        for ax,metric,title in zip(axes,['first_selected_sign_correct','selected'],['First selected sign','Selected candidate exists']):
            x=np.arange(3); width=.35
            for j,var in enumerate(['beforeGeo','afterGeo']):
                vals=[]
                for ds in keys:
                    a=[r for r in rows if r['dsid']==ds and r['variant']==var]
                    if metric=='first_selected_sign_correct': v=sum(r[metric] is True for r in a)/len(a)
                    else: v=sum(r[metric]>0 for r in a)/len(a)
                    vals.append(100*v)
                ax.bar(x+(j-.5)*width,vals,width,label=var)
                for i,v in enumerate(vals): ax.text(i+(j-.5)*width,v+1,f'{v:.1f}',ha='center',fontsize=8)
            ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_ylim(0,105); ax.set_ylabel('Components [%]'); ax.set_title(title); ax.grid(axis='y',alpha=.2); ax.legend()
        fig.suptitle('Finite common physical-prior CKF component comparison')
        fig.text(.5,.01,f'{components} moment-matched components per recorded seed | within C fraction {fraction:g} | actual Calypso CKF only',ha='center',fontsize=8)
        fig.tight_layout(rect=[0,.04,1,.93]); pdf.savefig(fig);save(fig,'component_outcomes');plt.close(fig)
        fig,ax=plt.subplots(figsize=(10,5))
        data=[]; labels=[]
        for ds in keys:
            for var,color in [('beforeGeo','#c54b4b'),('afterGeo','#2f67b1')]:
                v=[r['first_selected_qop_GeVinv']*1000 for r in rows if r['dsid']==ds and r['variant']==var and np.isfinite(r['first_selected_qop_GeVinv'])]
                data.append(v);labels.append(ds+' '+var)
        ax.boxplot(data,labels=labels,showfliers=False);ax.axhline(0,color='gray',lw=.8);ax.set_ylabel('First selected q/p [10⁻³ GeV⁻¹]');ax.set_title('Component posterior q/p');ax.tick_params(axis='x',rotation=25);ax.grid(axis='y',alpha=.2)
        fig.tight_layout();pdf.savefig(fig);save(fig,'component_qop');plt.close(fig)
        fig,ax=plt.subplots(figsize=(9,5))
        changes=[(int(r['after_correct'])-int(r['before_correct'])) if r['before_correct'] is not None and r['after_correct'] is not None else np.nan for r in pairs]
        values=[sum(x==v for x in changes) for v in [-1,0,1]]
        ax.bar(['correct→wrong','unchanged','wrong→correct'],values,color=['#c54b4b','#888888','#3a9257']);ax.set_ylabel('Paired components');ax.set_title('BeforeGeo → afterGeo first-candidate sign transitions');ax.grid(axis='y',alpha=.2)
        for i,v in enumerate(values):ax.text(i,v+max(values)*.02,str(v),ha='center')
        fig.tight_layout();pdf.savefig(fig);save(fig,'component_transitions');plt.close(fig)
    (args.output/'figure_index.json').write_text(json.dumps(index,indent=2)+'\n')
    print(json.dumps(dict(pages=len(index),output=str(args.output)),indent=2))


if __name__ == '__main__':
    main()
