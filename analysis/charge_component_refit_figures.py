#!/usr/bin/env python3
"""Plot paired finite-prior diagnostics from the actual Calypso full chain."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys

import numpy as np

from slope_tables import write_json


SAMPLES = ('100044', '100116', '100117')
LABELS = ('100044 mu+', '100116 mu+', '100117 mu-')
COLORS = ('#aa4b49', '#d49a37', '#326ba4')


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    pairs = json.loads((source/'paired_components.json').read_text())
    # Validate the stage populations and differences independently from the
    # plotted summary; components are correlated diagnostics, not events.
    for series, info in summary['series'].items():
        for dsid, groups in info['groups'].items():
            for comparison, group in groups.items():
                rows = [r for r in pairs if r['series'] == series and r['dsid'] == dsid
                        and r['comparison'] == comparison]
                counts = group['counts']
                assert len(rows) == counts['components']
                for side in ('before', 'after'):
                    assert sum(r[side+'_selected'] for r in rows) == counts[side+'_selected']
                    for stage in ('selected', 'converted', 'final'):
                        for correct in (True, False):
                            key = side+'_'+stage+('_correct' if correct else '_wrong')
                            assert sum(r[side+'_'+stage+'_correct'] is correct for r in rows) == counts.get(key, 0)
                for stage in ('selected', 'converted', 'final'):
                    differences = [abs(r['after_'+stage+'_qop_GeVinv']-r['before_'+stage+'_qop_GeVinv'])
                                   for r in rows if r['before_'+stage+'_qop_GeVinv'] is not None
                                   and r['after_'+stage+'_qop_GeVinv'] is not None]
                    measured = group['absolute_qop_differences_GeVinv'].get(stage)
                    if differences:
                        assert measured['n'] == len(differences)
                        for key, number in [('max', max(differences)), ('median', np.median(differences)),
                                            ('p95', np.quantile(differences, .95))]:
                            assert measured[key] == float(number)
    return dict(summary=summary, pairs=pairs)


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator, NullLocator
    from matplotlib.backends.backend_pdf import PdfPages

    data = values(source)
    series_list = list(data['summary']['series'])
    output.mkdir(parents=True, exist_ok=False)
    index = []

    def save(fig, name, title, pdf):
        fig.suptitle(title, fontsize=14)
        fig.text(.5, .035, 'Actual Calypso CKF / conversion / refit | Isolated smoothing ON | 20 source events per sample',
                 ha='center', fontsize=8)
        fig.text(.5, .016, 'Moment-matched finite prior components; first-order within-component covariance mapping | Counts are diagnostic ensembles, not production rates',
                 ha='center', fontsize=7)
        fig.tight_layout(rect=[0, .075, 1, .935])
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
                          reference='S06 common physical-prior component control in actual Calypso'))

    def series_label(series):
        info = data['summary']['series'][series]
        fraction = Fraction(info['within_fraction']).limit_denominator()
        return f"n={info['n_components']}; within C fraction={fraction}"

    def curvature_ticks(ax, limits):
        ticks = sorted([0.] + [sign*10.**power for sign in (-1, 1)
                              for power in range(-5, 5)
                              if limits[0] <= sign*10.**power <= limits[1]])

        def label(value, position):
            if value == 0:
                return '$0$'
            sign = '-' if value < 0 else ''
            return f'${sign}10^{{{int(round(np.log10(abs(value))))}}}$'

        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(FixedLocator(ticks))
            axis.set_major_formatter(FuncFormatter(label))
            axis.set_minor_locator(NullLocator())
        ax.tick_params(axis='both', labelsize=7)

    with PdfPages(output/'component_refit_report.pdf') as pdf:
        fig, axes = plt.subplots(len(series_list), 2, figsize=(14, 4.1*len(series_list)), squeeze=False)
        for i, series in enumerate(series_list):
            for j, stage in enumerate(('selected', 'final')):
                ax = axes[i, j]
                for mode, comparison, side, label, color, offset in [
                        ('before_native', 'full_native', 'before', 'Before native', COLORS[0], -.25),
                        ('before_regular', 'full_regular_before', 'before', 'Before regular refit', COLORS[1], 0.),
                        ('after', 'full_native', 'after', 'After native', COLORS[2], .25)]:
                    counts = [data['summary']['series'][series]['groups'][ds][comparison]['counts'] for ds in SAMPLES]
                    wrong = [c.get(side+'_'+stage+'_wrong', 0) for c in counts]
                    total = [c[side+'_selected'] for c in counts]
                    fractions = [100*w/n if n else 0. for w, n in zip(wrong, total)]
                    bars = ax.bar(np.arange(3)+offset, fractions, .24, label=label, color=color)
                    ax.bar_label(bars, labels=[f'{w}/{n}' for w, n in zip(wrong, total)], fontsize=7, padding=3)
                ax.set_xticks(np.arange(3), LABELS)
                ax.set_ylabel('Wrong-sign components / selected [%]')
                ax.set_ylim(0, 105)
                ax.set_title(('CKF selected' if stage == 'selected' else 'Final after native refit')+' | '+series_label(series), fontsize=10)
                ax.grid(axis='y', alpha=.2)
                if i == 0 and j == 0:
                    ax.legend(fontsize=8)
        save(fig, 'component_full_chain_outcomes', 'Paired finite prior: CKF and full-chain charge outcomes', pdf)

        fig, axes = plt.subplots(len(series_list), 3, figsize=(15, 4.1*len(series_list)), squeeze=False)
        for i, series in enumerate(series_list):
            for j, dsid in enumerate(SAMPLES):
                ax = axes[i, j]
                for comparison, label, color in [('full_native', 'Before native', COLORS[0]),
                                                  ('full_regular_before', 'Before regular refit', COLORS[1])]:
                    rows = [r for r in data['pairs'] if r['series'] == series and r['dsid'] == dsid
                            and r['comparison'] == comparison and r['before_selected'] and r['after_selected']]
                    x = [r['after_final_qop_GeVinv'] for r in rows]
                    y = [r['before_final_qop_GeVinv'] for r in rows]
                    ax.scatter(x, y, s=7, alpha=.55, color=color, label=label)
                ax.axhline(0, color='gray', linewidth=.7)
                ax.axvline(0, color='gray', linewidth=.7)
                ax.set_xscale('symlog', linthresh=1e-5)
                ax.set_yscale('symlog', linthresh=1e-5)
                limits = [min(*ax.get_xlim(), *ax.get_ylim()), max(*ax.get_xlim(), *ax.get_ylim())]
                ax.plot(limits, limits, ':', color='gray', linewidth=.8)
                ax.set_xlim(limits)
                ax.set_ylim(limits)
                curvature_ticks(ax, limits)
                ax.set_xlabel('After native final q/p [GeV^-1]')
                ax.set_ylabel('Before final q/p [GeV^-1]')
                ax.set_title(LABELS[j]+' | '+series_label(series), fontsize=9)
                ax.grid(alpha=.2)
                if i == 0 and j == 0:
                    ax.legend(fontsize=7)
        save(fig, 'component_final_qop_pairs', 'Same component: native and regular-chart refit curvature', pdf)

        fig, axes = plt.subplots(len(series_list), 3, figsize=(15, 4.1*len(series_list)), squeeze=False)
        for i, series in enumerate(series_list):
            for j, dsid in enumerate(SAMPLES):
                ax = axes[i, j]
                for comparison, stage, label, color, refit_only in [
                        ('CKF', 'selected', 'CKF selected', '#448557', False),
                        ('full_native', 'final', 'Final native', COLORS[0], False),
                        ('full_regular_before', 'final', 'Final regular before', COLORS[1], False),
                        ('full_regular_before', 'final', 'Regular; both refit OK', COLORS[2], True)]:
                    rows = [r for r in data['pairs'] if r['series'] == series and r['dsid'] == dsid
                            and r['comparison'] == comparison and r['before_selected'] and r['after_selected']]
                    if refit_only:
                        rows = [r for r in rows if r['before_refit_success'] and r['after_refit_success']]
                    if not rows:
                        continue
                    delta = np.sort([abs(r['after_'+stage+'_qop_GeVinv']-r['before_'+stage+'_qop_GeVinv']) for r in rows])
                    positive = delta > 0
                    rank = np.arange(1, len(delta)+1)/len(delta)
                    ax.step(delta[positive], rank[positive], where='post', label=label+f'; exact zero={sum(~positive)}', color=color)
                ax.set_xscale('log')
                ax.set_ylim(0, 1.02)
                ax.set_ylabel('Cumulative paired components')
                ax.set_xlabel('|After - before q/p| [GeV^-1]')
                ax.set_title(LABELS[j]+' | '+series_label(series), fontsize=9)
                ax.grid(alpha=.2)
                ax.legend(fontsize=7, loc='lower right')
        save(fig, 'component_qop_difference_cdf', 'Residual coordinate differences before and after the refit chart control', pdf)

        fig, axes = plt.subplots(len(series_list), 3, figsize=(15, 4.1*len(series_list)), squeeze=False)
        for i, series in enumerate(series_list):
            for j, dsid in enumerate(SAMPLES):
                ax = axes[i, j]
                categories = ['CKF q/p', 'Converted q/p', 'Final q/p']
                for comparison, label, color, offset in [
                        ('full_native', 'Before native -> after', COLORS[0], -.22),
                        ('before_refit_chart', 'Before native -> regular', COLORS[1], 0.),
                        ('full_regular_before', 'Before regular -> after', COLORS[2], .22)]:
                    c = data['summary']['series'][series]['groups'][dsid][comparison]['counts']
                    counts = [c.get(stage+'_sign_migrations', 0) for stage in ('selected', 'converted', 'final')]
                    bars = ax.bar(np.arange(3)+offset, counts, .21, color=color, label=label)
                    ax.bar_label(bars, fontsize=7, padding=2)
                ax.set_xticks(np.arange(3), categories, fontsize=8)
                ax.set_ylabel('Sign migrations [components]')
                ax.set_title(LABELS[j]+' | '+series_label(series), fontsize=9)
                ax.set_ylim(bottom=0, top=max(1., ax.get_ylim()[1]*1.2))
                ax.yaxis.set_major_locator(MaxNLocator(integer=True))
                ax.grid(axis='y', alpha=.2)
                if i == 0 and j == 0:
                    ax.legend(fontsize=7)
        save(fig, 'component_stage_sign_migrations', 'Locate sign differences along the actual CKF / conversion / refit chain', pdf)

        fig, axes = plt.subplots(len(series_list), 3, figsize=(15, 4.1*len(series_list)), squeeze=False)
        for i, series in enumerate(series_list):
            for j, dsid in enumerate(SAMPLES):
                ax = axes[i, j]
                text = [LABELS[j], series_label(series)]
                for comparison, label in [('full_native', 'Native before / native after'),
                                          ('full_regular_before', 'Regular before / native after')]:
                    c = data['summary']['series'][series]['groups'][dsid][comparison]['counts']
                    text += ['', label, f"Input components: {c['components']}",
                             f"Selected both: {c['selected_in_both']}",
                             f"Before only / after only: {c['before_only_selected']} / {c['after_only_selected']}",
                             f"Same selected hit path: {c['same_selected_hit_sequence']}",
                             f"Refit success both: {c.get('refit_success_both', 0)}",
                             f"Refit success before / after: {c.get('before_refit_success', 0)} / {c.get('after_refit_success', 0)}",
                             f"Refit fallback before / after: {c.get('before_refit_fallback', 0)} / {c.get('after_refit_fallback', 0)}"]
                ax.axis('off')
                ax.text(.02, .98, '\n'.join(text), va='top', fontsize=9, family='monospace')
        save(fig, 'component_coverage_and_refit', 'Component coverage, selected paths and actual refit outcomes', pdf)

        fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))
        labels = [f"n={data['summary']['series'][s]['n_components']}\nC/{Fraction(data['summary']['series'][s]['within_fraction']).limit_denominator().denominator}" for s in series_list]
        for j, dsid in enumerate(SAMPLES):
            ax = axes[j]
            selected, wrong, refitted = [], [], []
            for series in series_list:
                c = data['summary']['series'][series]['groups'][dsid]['full_native']['counts']
                selected.append(100*c['after_selected']/c['components'])
                wrong.append(100*c.get('after_final_wrong', 0)/c['after_selected'])
                refitted.append(100*c.get('after_refit_success', 0)/c['after_selected'])
            for numbers, label, color in [(selected, 'Selected / input components', '#448557'),
                                           (wrong, 'Final wrong / selected', COLORS[0]),
                                           (refitted, 'Refit success / selected', COLORS[2])]:
                ax.plot(np.arange(len(series_list)), numbers, '.-', label=label, color=color)
            ax.set_xticks(np.arange(len(series_list)), labels, fontsize=7)
            ax.set_ylim(0, 105)
            ax.set_ylabel('Diagnostic component fraction [%]')
            ax.set_title(LABELS[j])
            ax.grid(alpha=.2)
            if j == 0:
                ax.legend(fontsize=7)
        save(fig, 'component_sampling_sensitivity', 'Component count and covariance split: approximation sensitivity, afterGeo', pdf)
    write_json(output/'plot_values.json', data)
    write_json(output/'figure_index.json', index)


def validate(source, output, pdf_tools):
    from PIL import Image
    sys.path.insert(0, str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics as base_normalized_graphics
    import re

    def fonts(page):
        return {str(name): obj.get_object() for name, obj in page['/Resources']['/Font'].items()}

    def font_name(font):
        return re.sub(r'^/[A-Z]{6}\+', '/', str(font['/BaseFont']))

    def normalized_graphics(page):
        aliases = fonts(page)
        return re.sub(r'/F\d+(?=\s|-)', lambda match: font_name(aliases[match[0]]),
                      base_normalized_graphics(page))

    def check_font_glyphs(single, page):
        # PdfPages may assign /F1 to Sans and /F2 to Mono, while the single
        # coverage page first uses Mono. Resolve aliases and check the actual
        # Type3 glyph outlines and character codes used by each font subset.
        report_fonts = {font_name(font): font for font in fonts(page).values()}
        for font in fonts(single).values():
            other = report_fonts[font_name(font)]
            assert font['/Subtype'] == other['/Subtype'] == '/Type3'
            glyphs = font['/CharProcs']
            for name, glyph in glyphs.items():
                assert glyph.get_object().get_data() == other['/CharProcs'][name].get_object().get_data()

            def codes(current):
                result, code = {}, 0
                for item in current['/Encoding']['/Differences']:
                    if isinstance(item, int):
                        code = item
                    else:
                        result[code] = str(item)
                        code += 1
                return result

            a, b = codes(font), codes(other)
            for code, glyph in a.items():
                if glyph in glyphs:
                    assert b[code] == glyph
                    assert font['/Widths'][code-int(font['/FirstChar'])] == other['/Widths'][code-int(other['/FirstChar'])]

    assert json.loads((output/'plot_values.json').read_text()) == values(source)
    index = json.loads((output/'figure_index.json').read_text())
    report = PdfReader(output/'component_refit_report.pdf')
    assert len(index) == len(report.pages) == 6
    for row, page in zip(index, report.pages):
        with Image.open(output/row['png']) as image:
            image.load()
        single = PdfReader(output/row['pdf']).pages[0]
        check_font_glyphs(single, page)
        assert normalized_graphics(single) == normalized_graphics(page), row['name']
        assert list(single.mediabox) == list(page.mediabox)
        assert list(single.cropbox) == list(page.cropbox) and single.rotation == page.rotation
    write_json(output/'validation.json', dict(pages=len(index), plotted_values_exact=True,
        stage_counts_and_qop_metrics_recomputed_from_pairs=True, PNG_decoded=True,
        normalized_vector_graphics_and_boxes_equal=True, font_glyph_outlines_codes_and_widths_equal=True,
        visual_inspection='Pending full-size review'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', choices=['plot', 'validate'], required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--pdf-tools', type=Path)
    args = parser.parse_args()
    if args.phase == 'plot':
        plot(args.input, args.output)
    else:
        validate(args.input, args.output, args.pdf_tools)


if __name__ == '__main__':
    main()
