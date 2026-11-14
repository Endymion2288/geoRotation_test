#!/usr/bin/env python3
"""Independently compare joint/transport plot payloads to measured summaries."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from PIL import Image

from slope_tables import write_json
from validate_charge_gun_figures import normalized_graphics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kind', choices=['joint', 'transport'], required=True)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--figures', type=Path, required=True)
    parser.add_argument('--pdf-tools', type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.pdf_tools))
    from pypdf import PdfReader
    summary = json.loads(args.summary.read_text())
    counts = json.loads((args.figures/'bin_counts.json').read_text())
    index = json.loads((args.figures/'figure_index.json').read_text())
    comparisons, histogram_count = 0, 0
    if args.kind == 'joint':
        for dsid in ['100116', '100117', '110033', '110051']:
            modes = ['3forward', '3backward', '3merge', '4forward', '4backward', '4merge']
            for control in ['native', 'both']:
                groups = [summary['groups'][f'{dsid}|{m}|{control}'] for m in modes]
                for metric in ['conditional_misid', 'correct_efficiency',
                        'momentum_below_half_production', 'angle_above_1mrad']:
                    assert counts[f'{dsid}|{metric}|{control}'] == [g[metric] for g in groups]
                    comparisons += 1
                assert counts[f'{dsid}|angle68|{control}'] == [g['angle']['absolute_residual_68']*1e6 for g in groups]
                assert counts[f'{dsid}|CPU_ms|{control}'] == [summary['resources'][f'{dsid}|{m}|{control}']['process_cpu_ms_per_source_event'] for m in modes]
                fallback = []
                for mode in modes:
                    r = summary['refit_counters'][f'{dsid}|{mode}|{control}']
                    row = r['merged']['selected_merged_output'] if mode.endswith('merge') else r['selected']
                    fallback.append(100*row.get('fallback', 0)/row['total'] if row['total'] else None)
                np.testing.assert_allclose(counts[f'{dsid}|selected_fallback|{control}'], fallback)
                comparisons += 3
                for mode, group in zip(modes, groups):
                    for field in ['relative_p', 'signed_qop', 'angle']:
                        expected = group['angle']['residual_fine'] if field == 'angle' else group['histograms'][field]
                        assert counts[f'{dsid}|{mode}|{control}|{field}'] == expected
                        histogram_count += 1
                        comparisons += 1
            for mode in modes:
                row = summary['comparisons'][dsid+'|'+mode]
                assert counts[f'{dsid}|{mode}|paired_changes'] == dict(
                    correct_efficiency=row['correct_efficiency_difference'],
                    conditional_misid=row['conditional_misid_difference'])
                comparisons += 1
        report, pages = 'joint_regression_report.pdf', 18
    else:
        for dsid in ['100044', '100116', '100117']:
            for control in [r['key'] for r in summary.get('plot_controls',
                    [dict(key='material_on'), dict(key='material_off')])]:
                for target, group in enumerate(summary['groups'][dsid+'|'+control]):
                    for field in ['position_mm', 'direction_rad', 'relative_p']:
                        assert counts[f'{dsid}|{control}|{target}|{field}'] == group['histograms'][field]
                        comparisons += 1
                        histogram_count += 1
        report, pages = 'physical_prior_transport_report.pdf', 3
    for row in counts.values():
        if isinstance(row, dict) and 'counts' in row and 'edges' in row:
            assert sum(row['counts'])+row['underflow']+row['overflow'] == row['finite']
            assert sum(row['sumw2']) == sum(row['counts'])
    pdf = PdfReader(args.figures/report)
    assert len(pdf.pages) == len(index) == pages
    page_checks = []
    for page, row in zip(pdf.pages, index):
        with Image.open(args.figures/row['png']) as image:
            image.load()
        single = PdfReader(args.figures/row['pdf'])
        assert len(single.pages) == 1
        other = single.pages[0]
        assert normalized_graphics(page) == normalized_graphics(other)
        assert (list(page.mediabox), list(page.cropbox), page.rotation) == \
               (list(other.mediabox), list(other.cropbox), other.rotation)
        page_checks.append(dict(figure=row['name'], PNG_fully_decoded=True,
            vector_graphics_and_page_boxes_equal=True))
    write_json(args.figures/'validation.json', dict(pages=pages, all_plotted_data_exact=True,
        summary_comparisons=comparisons, histogram_count=histogram_count,
        page_checks=page_checks, visual_inspection='Pending separate image review'))
    print('Verified', comparisons, 'summary comparisons,', histogram_count, 'histograms,', pages, 'PDF pages')


if __name__ == '__main__':
    main()
