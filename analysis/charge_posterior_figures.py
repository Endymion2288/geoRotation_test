#!/usr/bin/env python3
"""Visualize complete coverage and native fixed-observation likelihoods."""
import argparse
import json
from pathlib import Path
import re

import numpy as np

from charge_component_refit import audit
from slope_tables import write_json


SAMPLES = ('100044', '100116', '100117')
VARIANTS = ('beforeGeo', 'afterGeo')
COLORS = dict(beforeGeo='#b44e3d', afterGeo='#236ca4')


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    components = list(audit(source/'components.jsonl.gz'))
    events = list(audit(source/'events.jsonl.gz'))
    comparisons = list(audit(source/'comparisons.jsonl.gz'))
    groups, pairs, examples = summary['groups'], {}, {}
    for dsid in SAMPLES:
        records = {(r['variant'], r['entry'], r['component']):r for r in components if r['dsid'] == dsid}
        paired = []
        for (variant, entry, component), row in records.items():
            if variant != 'beforeGeo':
                continue
            other = records['afterGeo', entry, component]
            if row['complete'] and other['complete']:
                paired.append(dict(entry=entry, component=component,
                    before_logL=row['log_likelihood'], after_logL=other['log_likelihood'],
                    difference=other['log_likelihood']-row['log_likelihood']))
        pairs[dsid] = paired
        complete_entries = [r['entry'] for r in comparisons if r['dsid'] == dsid and r['complete_both']]
        if complete_entries:
            candidates = [r for r in events if r['dsid'] == dsid and r['variant'] == 'afterGeo'
                          and r['entry'] in complete_entries]
            chosen = min(candidates, key=lambda r:(r['posterior']['effective_components'], r['entry']))
            entry = chosen['entry']
            selected = {v:[r for r in components if r['dsid'] == dsid and r['entry'] == entry
                           and r['variant'] == v] for v in VARIANTS}
            low = min(r['smoothed_mean'][4]-5*np.sqrt(r['smoothed_cov'][4][4])
                      for rows in selected.values() for r in rows)
            high = max(r['smoothed_mean'][4]+5*np.sqrt(r['smoothed_cov'][4][4])
                       for rows in selected.values() for r in rows)
            x = np.linspace(low, high, 1001)
            density = {}
            for variant, rows in selected.items():
                y = np.zeros_like(x)
                for row in rows:
                    mean, var = row['smoothed_mean'][4], row['smoothed_cov'][4][4]
                    y += row['posterior_weight']*np.exp(-.5*(x-mean)**2/var)/np.sqrt(2*np.pi*var)
                density[variant] = y.tolist()
            examples[dsid] = dict(entry=entry, x=x.tolist(), density=density,
                event_results={v:next(r for r in events if r['dsid'] == dsid and r['entry'] == entry
                                     and r['variant'] == v) for v in VARIANTS})
    return dict(groups=groups, paired_components=pairs, event_comparisons=comparisons,
        events=events, examples=examples)


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    data = values(source)
    output.mkdir(parents=True, exist_ok=False)
    index = []

    def save(fig, name, title, pdf):
        fig.suptitle(title, fontsize=15, y=.985)
        fig.text(.5, .04, 'Actual Calypso | Isolated smoothing ON | CKF step limit 10 mm | One frozen observation sequence per event',
            ha='center', fontsize=9)
        fig.text(.5, .018, '60 source events including 5 without a native seed | n32, within-C / 256 | Diagnostic likelihood; not a production misID rate',
            ha='center', fontsize=9)
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
            reference='S06 Q5 common finite prior and fixed-observation native likelihood'))

    with PdfPages(output/'posterior_report.pdf') as pdf:
        fig, ax = plt.subplots(figsize=(16, 8))
        ax.axis('off')
        cells = []
        for dsid in SAMPLES:
            for variant in VARIANTS:
                c = data['groups'][dsid+'|'+variant]
                cells.append([dsid, variant, str(c['input_components']), str(c['selected_candidates']),
                    str(c.get('complete', 0)), str(c['input_components']-c.get('complete', 0)),
                    str(c.get('terminal_only', 0)), str(c.get('refit_success', 0)),
                    f"{c.get('complete_event_posteriors', 0)} / {20-c.get('events_without_prior_components', 0)}"])
        table = ax.table(cellText=cells, colLabels=['DSID', 'Geometry', 'All prior\ncomponents', 'Selected\ncandidates',
            'Complete\nobservations', 'Incomplete\ncomponents', 'Terminal\nonly', 'Refit\nsuccess', 'Complete posterior /\nevents with a prior'],
            cellLoc='center', bbox=[0., .30, 1., .48])
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        for (row, _), cell in table.get_celld().items():
            if row == 0:
                cell.set_facecolor('#dce6ef')
        ax.text(.5, .86, 'Every prior component is retained. A missing observation path has an unmeasured likelihood.',
            ha='center', fontsize=11, transform=ax.transAxes)
        ax.text(.5, .14, 'The native CKF selector uses chi2Max = 1e30 and nMax = 1 only for this observation control.\n'
            'A full finite-prior posterior is normalized only when every component has the complete same observations.\n'
            'Refit fallback is retained; the refit reuses measurements and is not multiplied in as an independent likelihood.',
            ha='center', fontsize=10, linespacing=1.8, transform=ax.transAxes)
        fig.tight_layout(rect=[0, .085, 1, .94])
        save(fig, 'posterior_coverage', 'Fixed-observation native CKF: complete input and likelihood coverage', pdf)

        fig, axes = plt.subplots(2, 3, figsize=(16, 9))
        for column, dsid in enumerate(SAMPLES):
            rows = data['paired_components'][dsid]
            a, b = (np.array([r[k] for r in rows]) for k in ('before_logL', 'after_logL'))
            ax = axes[0, column]
            ax.scatter(a, b-a, s=8, alpha=.55)
            ax.axhline(0., color='black', lw=.7)
            ax.set_xlabel('Before-geometry log likelihood')
            ax.set_ylabel('After minus before log likelihood')
            ax.set_title(dsid+f' | {len(rows)} complete component pairs')
            ax.grid(alpha=.2)
            ax = axes[1, column]
            if len(a):
                x = np.sort(np.maximum(np.abs(b-a), 1e-16))
                ax.step(x, (np.arange(len(x))+1)/len(x), where='post')
            ax.set_xscale('log')
            ax.set_ylim(0., 1.025)
            ax.set_xlabel('Absolute log-likelihood difference (zero at 1e-16)')
            ax.set_ylabel('Fraction of complete component pairs')
            ax.grid(alpha=.2)
        fig.tight_layout(rect=[0, .09, 1, .94])
        save(fig, 'posterior_likelihood', 'Native Gaussian innovations: common observations, both geometry versions', pdf)

        fig, axes = plt.subplots(2, 3, figsize=(16, 9))
        for column, dsid in enumerate(SAMPLES):
            rows = [r for r in data['event_comparisons'] if r['dsid'] == dsid and r['complete_both']]
            for ax, field, label in [(axes[0, column], 'weight_total_variation', 'Posterior-weight total variation'),
                (axes[1, column], 'qop_mean_difference_over_posterior_sigma', 'q/p mean difference / posterior sigma')]:
                if rows:
                    ax.plot([r['entry'] for r in rows], np.maximum([r[field] for r in rows], 1e-18), 'o', ms=5)
                else:
                    ax.text(.5, .5, 'No complete event posterior', ha='center', transform=ax.transAxes)
                ax.set_yscale('log')
                ax.set_xlabel('Local source entry')
                ax.set_ylabel(label+' (floor: 1e-18)')
                ax.set_title(dsid+f' | {len(rows)} complete event pairs')
                ax.grid(alpha=.2)
        fig.tight_layout(rect=[0, .09, 1, .94])
        save(fig, 'posterior_comparison', 'Full finite-prior posterior comparison: no missing mass removed', pdf)

        fig, axes = plt.subplots(1, 3, figsize=(16, 7))
        from matplotlib.ticker import MaxNLocator
        for ax, dsid in zip(axes, SAMPLES):
            example = data['examples'].get(dsid)
            if example:
                for variant in VARIANTS:
                    event = example['event_results'][variant]
                    ax.plot(example['x'], example['density'][variant], color=COLORS[variant],
                        ls='-' if variant == 'beforeGeo' else '--', lw=1.7,
                        label=variant+f" | ESS={event['posterior']['effective_components']:.2f}")
                ax.set_title(dsid+f" | entry {example['entry']}")
                ax.legend(fontsize=9)
            else:
                ax.text(.5, .5, 'No complete event posterior', ha='center', transform=ax.transAxes)
                ax.set_title(dsid)
            ax.set_xlabel('Signed q/p at the first common sensor [GeV$^{-1}$]')
            ax.set_ylabel('Finite-mixture posterior density [GeV]')
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
            ax.axvline(0., color='black', lw=.7)
            ax.grid(alpha=.2)
        fig.text(.5, .10, 'Example per DSID: the complete event with the lowest after-geometry effective component count.',
            ha='center', fontsize=10)
        fig.tight_layout(rect=[0, .13, 1, .93])
        save(fig, 'posterior_qop', 'Signed q/p: complete native-component Gaussian mixtures', pdf)
    write_json(output/'figure_index.json', index)
    write_json(output/'plot_values.json', data)


def validate(source, output, pdf_tools):
    import sys
    from PIL import Image
    sys.path.insert(0, str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics
    saved = json.loads((output/'plot_values.json').read_text())
    recomputed = values(source)
    density_checks = {}
    # SIMD exp implementations can differ in the final bits across worker CPUs.
    # All recorded inputs, event results, x coordinates and likelihoods remain
    # exact; only the newly evaluated analytic density curve uses a tolerance.
    for dsid, example in recomputed['examples'].items():
        for variant, density in example['density'].items():
            a = np.asarray(saved['examples'][dsid]['density'][variant])
            b = np.asarray(density)
            np.testing.assert_allclose(a, b, rtol=5e-13, atol=1e-300)
            mask = np.maximum(np.abs(a), np.abs(b)) > 1e-300
            density_checks[dsid+'|'+variant] = dict(max_absolute_difference=float(np.abs(a-b).max()),
                max_relative_difference=float(np.max(np.abs((a[mask]-b[mask])/b[mask]))) if mask.any() else 0.)
            saved['examples'][dsid]['density'][variant] = density
    assert saved == recomputed

    def graphics(page):
        fonts = page['/Resources']['/Font']
        return re.sub(r'/F\d+(?=\s|-)', lambda match: re.sub(r'^/[A-Z]{6}\+', '/',
            str(fonts[match[0]].get_object()['/BaseFont'])), normalized_graphics(page))
    index = json.loads((output/'figure_index.json').read_text())
    pages = PdfReader(output/'posterior_report.pdf').pages
    assert len(index) == len(pages) == 4
    for row, page in zip(index, pages):
        with Image.open(output/row['png']) as png:
            png.load()
        original = PdfReader(output/row['pdf']).pages[0]
        assert graphics(original) == graphics(page)
        assert (list(original.mediabox), list(original.cropbox), original.rotation) == (
            list(page.mediabox), list(page.cropbox), page.rotation)
    write_json(output/'validation.json', dict(pages=4, recorded_plot_inputs_exact=True,
        analytic_density_rtol=5e-13, analytic_density_atol=1e-300, density_checks=density_checks,
        PNG_decoded=True, normalized_PDF_content_and_boxes_equal=True,
        visual_inspection='Pending original-resolution inspection'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--phase', choices=['plot', 'validate'], required=True)
    parser.add_argument('--pdf-tools', type=Path)
    args = parser.parse_args()
    if args.phase == 'plot':
        plot(args.input, args.output)
    else:
        validate(args.input, args.output, args.pdf_tools)


if __name__ == '__main__':
    main()
