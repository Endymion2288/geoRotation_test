#!/usr/bin/env python3
"""Plot the recorded native refit priors and diagonal-inflation control."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

from slope_tables import write_json


SAMPLES = ('100044', '100116', '100117')
LABELS = ('100044 mu+', '100116 mu+', '100117 mu-')
COLORS = ('#aa4b49', '#d49a37', '#326ba4')
SERIES = 'matched_n32_f256_v1'


def values(source):
    summary = json.loads((source/'summary.json').read_text())
    rows = json.loads((source/'paired_refit_inputs.json').read_text())
    groups, paired = {}, {}
    for dsid in SAMPLES:
        groups[dsid] = {}
        for mode in ('native', 'regular_both'):
            for control in (False, True):
                name = mode+('_scale1' if control else '')
                selected = [r for r in rows if r['series'] == SERIES and r['dsid'] == dsid
                            and r['mode'] == name]
                measured = summary['scale1_groups' if control else 'groups'][SERIES+'|'+dsid+'|'+mode]
                counts = measured['counts']
                assert len(selected) == counts['selected_both']
                for side in ('before', 'after'):
                    assert sum(r[side+'_has_refit_input'] for r in selected) == counts['refit_input_'+side]
                assert sum(r['before_has_refit_input'] and r['after_has_refit_input'] for r in selected) == counts['refit_input_both']
                migration = [r for r in selected if r['final_sign_difference']]
                assert len(migration) == len(measured['final_sign_difference_examples'])
                for cohort, predicate in [
                        ('all_input_pairs', lambda r: True),
                        ('same_CKF_hits', lambda r: r['same_CKF_hits']),
                        ('both_refit_success', lambda r: r['before_refit_success'] and r['after_refit_success']),
                        ('same_hits_both_refit_success', lambda r: r['same_CKF_hits'] and r['before_refit_success'] and r['after_refit_success'])]:
                    available = [r for r in selected if 'physical_covariance_scaled_max' in r and predicate(r)]
                    frozen = measured['cohorts'][cohort]
                    assert len(available) == frozen['n']
                    for field, metric in frozen['metrics'].items():
                        numbers = [r[field] for r in available]
                        assert metric['n'] == len(numbers)
                        if numbers:
                            assert metric['max'] == float(max(numbers))
                            assert metric['median'] == float(np.median(numbers))
                            assert metric['p95'] == float(np.quantile(numbers, .95))
                groups[dsid][name] = dict(selected=len(selected),
                    before_success=sum(r['before_refit_success'] for r in selected),
                    after_success=sum(r['after_refit_success'] for r in selected),
                    sign_differences=len(migration), same_hits_both_success_sign_differences=sum(
                        r['same_CKF_hits'] and r['before_refit_success'] and r['after_refit_success'] for r in migration))
        default = {(r['entry'], r['component']): r for r in rows if r['dsid'] == dsid and r['mode'] == 'native'}
        control = {(r['entry'], r['component']): r for r in rows if r['dsid'] == dsid and r['mode'] == 'native_scale1'}
        assert default.keys() == control.keys()
        paired[dsid] = []
        for key, a in default.items():
            b = control[key]
            if not all(r['same_CKF_hits'] and r['before_refit_success'] and r['after_refit_success'] for r in (a, b)):
                continue
            item = dict(entry=key[0], component=key[1], default_sign_difference=a['final_sign_difference'])
            for field in ('physical_covariance_scaled_max', 'direction_covariance_scaled_max'):
                item[field] = [a[field], b[field]]
            for field in ('position_difference_mm', 'direction_difference_norm', 'qop_difference_GeVinv'):
                assert a[field] == b[field]
            paired[dsid].append(item)
    cases = []
    for dsid, entry, component in [('100044', 1, 23), ('100044', 3, 19),
                                  ('100044', 12, 2), ('100117', 0, 25)]:
        case = dict(dsid=dsid, entry=entry, component=component, modes={})
        for mode in ('native', 'regular_both', 'native_scale1', 'regular_both_scale1'):
            matched = [r for r in rows if r['dsid'] == dsid and r['entry'] == entry
                       and r['component'] == component and r['mode'] == mode]
            assert len(matched) == 1
            case['modes'][mode] = {k: matched[0][k] for k in ['before_final_qop_GeVinv',
                'after_final_qop_GeVinv', 'same_CKF_hits', 'before_refit_success', 'after_refit_success']}
        cases.append(case)
    return dict(groups=groups, paired_covariances=paired, cases=cases,
                source_summary=str(source/'summary.json'), source_pairs=str(source/'paired_refit_inputs.json'))


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.ticker import MaxNLocator

    data = values(source)
    output.mkdir(parents=True, exist_ok=False)
    index = []

    def save(fig, name, title, caption, pdf):
        fig.suptitle(title, fontsize=14)
        fig.text(.5, .037, 'Actual native Calypso refit | Isolated smoothing ON | n=32, within C/256 | 20 source events per sample', ha='center', fontsize=8)
        fig.text(.5, .017, caption, ha='center', fontsize=8)
        fig.tight_layout(rect=[0, .085, 1, .94])
        fig.savefig(output/(name+'.png'), dpi=150)
        fig.savefig(output/(name+'.pdf'))
        pdf.savefig(fig)
        plt.close(fig)
        index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
                          reference='S06 recorded native refit inputs and diagonal covariance scaling control'))

    with PdfPages(output/'refit_input_report.pdf') as pdf:
        fig, axes = plt.subplots(2, 3, figsize=(14, 9))
        for j, dsid in enumerate(SAMPLES):
            rows = data['paired_covariances'][dsid]
            for i, (field, label) in enumerate([
                    ('physical_covariance_scaled_max', 'Full physical covariance'),
                    ('direction_covariance_scaled_max', 'Direction covariance')]):
                ax = axes[i, j]
                x, y = np.array([r[field] for r in rows]).T
                limit = 1.06*max(max(x), max(y))
                ax.plot([0, limit], [0, limit], '--', color='gray', linewidth=.8)
                ax.scatter(x, y, color=COLORS[j], s=19, alpha=.65)
                special = [r[field] for r in rows if r['default_sign_difference']]
                if special:
                    ax.scatter(*np.array(special).T, edgecolors='black', facecolors='none', s=90,
                               label='Native final sign differs at scale 10')
                    ax.legend(fontsize=7, loc='upper left')
                ax.set(xlim=(0, limit), ylim=(0, limit), xlabel='Before/after discrepancy: diagonal scale 10',
                       ylabel='Before/after discrepancy: diagonal scale 1')
                ax.set_title(LABELS[j]+' | '+label, fontsize=10)
                ax.text(.97, .06, f'N={len(rows)}; paired medians\n{np.median(x):.3f} -> {np.median(y):.3f}',
                        transform=ax.transAxes, ha='right', fontsize=9)
                ax.grid(alpha=.2)
        save(fig, 'refit_input_covariance', 'Diagonal inflation amplifies the physical refit-prior discrepancy',
             'Same CKF hits and both refits successful in both settings; max |Cb-Ca| / sqrt(mean Cii * mean Cjj), common FASER frame', pdf)

        fig, axes = plt.subplots(2, 2, figsize=(13, 9))
        for ax, case in zip(axes.flat, data['cases']):
            charge = -1 if case['dsid'] == '100117' else 1
            for side, base, label, color, marker in [
                    ('before', 'native', 'Before native', COLORS[0], 'o'),
                    ('before', 'regular_both', 'Before regular refit', COLORS[1], 's'),
                    ('after', 'native', 'After native', COLORS[2], '^')]:
                states = [case['modes'][base+suffix] for suffix in ('', '_scale1')]
                y = [1000*charge*s[side+'_final_qop_GeVinv'] for s in states]
                ax.plot([0, 1], y, marker+'-', color=color, label=label, linewidth=1.2)
                for x, value, state in zip((0, 1), y, states):
                    if not state[side+'_refit_success']:
                        ax.scatter([x], [value], marker='x', color='black', s=75, zorder=5)
            native = case['modes']['native']
            status = 'Same CKF hits; both refits successful' if native['same_CKF_hits'] else 'Different CKF hits; after refit fallback'
            ax.set_title(f"{case['dsid']} / entry {case['entry']} / component {case['component']}\n"+status, fontsize=10)
            ax.axhline(0, color='gray', linewidth=.8)
            ax.set_xticks([0, 1], ['Diagonal scale 10', 'Diagonal scale 1'])
            ax.set_ylabel(r'$q_{truth}(q/p)_{final}$ [$10^{-3}$ GeV$^{-1}$]')
            ax.grid(alpha=.2)
        axes[0, 0].legend(fontsize=8)
        save(fig, 'refit_input_sign_cases', 'Four original native sign differences: identical signs can still be wrong',
             'Positive = correct sign; negative = wrong sign; crosses = refit fallback | Only extra diagonal inflation changes; complete C x10 remains', pdf)

        fig, axes = plt.subplots(2, 3, figsize=(15, 9))
        for j, dsid in enumerate(SAMPLES):
            groups = data['groups'][dsid]
            ax = axes[0, j]
            modes = ['native', 'native_scale1', 'regular_both', 'regular_both_scale1']
            same = [groups[m]['same_hits_both_success_sign_differences'] for m in modes]
            other = [groups[m]['sign_differences']-n for m, n in zip(modes, same)]
            ax.bar(np.arange(4), same, color='#aa4b49', label='Same hits, both refits successful')
            ax.bar(np.arange(4), other, bottom=same, color='#326ba4', label='Different hits or refit outcome')
            for k, m in enumerate(modes):
                ax.text(k, groups[m]['sign_differences']+.08, str(groups[m]['sign_differences']), ha='center', fontsize=9)
            ax.set_xticks(np.arange(4), ['Native\n10', 'Native\n1', 'Regular\n10', 'Regular\n1'], fontsize=8)
            ax.set_ylim(0, 4.4)
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.set_title(LABELS[j]+' | before/after sign differences')
            ax.set_ylabel('Paired diagnostic components')
            ax.grid(axis='y', alpha=.2)
            if j == 0:
                ax.legend(fontsize=7)
            ax = axes[1, j]
            successful, totals = [], []
            for suffix in ('', '_scale1'):
                for side, mode in [('before', 'native'), ('before', 'regular_both'), ('after', 'native')]:
                    group = groups[mode+suffix]
                    successful.append(group[side+'_success'])
                    totals.append(group['selected'])
            fallback = np.array(totals)-successful
            bars = ax.bar(np.arange(6), successful, color='#448557', label='Refit successful')
            ax.bar(np.arange(6), fallback, bottom=successful, color='#d7d7d7', label='Converted CKF fallback')
            ax.bar_label(bars, fontsize=8, padding=2)
            ax.set_xticks(np.arange(6), ['Before\n10', 'Regular\n10', 'After\n10', 'Before\n1', 'Regular\n1', 'After\n1'], fontsize=8)
            ax.set_ylim(0, max(totals)*1.07)
            ax.set_ylabel('Selected diagnostic components')
            ax.set_title('Actual refit success and fallback; labels show success counts', fontsize=9)
            if j == 0:
                ax.legend(fontsize=8, loc='upper left')
        save(fig, 'refit_input_migrations_coverage', 'Residual sign differences and actual refit coverage',
             'Native/regular describes the beforeGeo refit chart; afterGeo refit remains native | Components are correlated controls, not production rates', pdf)
    write_json(output/'plot_values.json', data)
    write_json(output/'figure_index.json', index)


def validate(source, output, pdf_tools):
    import re
    from PIL import Image
    sys.path.insert(0, str(pdf_tools))
    from pypdf import PdfReader
    from validate_charge_gun_figures import normalized_graphics as base_graphics

    def fonts(page):
        return {str(k): v.get_object() for k, v in page['/Resources']['/Font'].items()}

    def font_name(font):
        return re.sub(r'^/[A-Z]{6}\+', '/', str(font['/BaseFont']))

    def graphics(page):
        aliases = fonts(page)
        return re.sub(r'/F\d+(?=\s|-)', lambda m: font_name(aliases[m[0]]), base_graphics(page))

    assert json.loads((output/'plot_values.json').read_text()) == values(source)
    index = json.loads((output/'figure_index.json').read_text())
    report = PdfReader(output/'refit_input_report.pdf')
    assert len(index) == len(report.pages) == 3
    for row, page in zip(index, report.pages):
        with Image.open(output/row['png']) as png:
            png.load()
        single = PdfReader(output/row['pdf']).pages[0]
        assert graphics(single) == graphics(page), row['name']
        other_fonts = {font_name(f): f for f in fonts(page).values()}
        for font in fonts(single).values():
            other = other_fonts[font_name(font)]
            assert font['/Subtype'] == other['/Subtype'] == '/Type3'
            for name, glyph in font['/CharProcs'].items():
                assert glyph.get_object().get_data() == other['/CharProcs'][name].get_object().get_data()
        assert list(single.mediabox) == list(page.mediabox)
        assert list(single.cropbox) == list(page.cropbox) and single.rotation == page.rotation
    write_json(output/'validation.json', dict(pages=3, plotted_values_exact=True,
        counts_and_covariance_metrics_recomputed=True, PNG_decoded=True,
        normalized_vector_graphics_boxes_and_glyphs_equal=True,
        visual_inspection='Pending original-resolution review'))


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
