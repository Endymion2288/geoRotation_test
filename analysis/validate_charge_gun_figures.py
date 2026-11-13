#!/usr/bin/env python3
"""Compare every gun plot input and resolve equivalent PDF resource names."""
import argparse
import json
from pathlib import Path
import re
import sys

from PIL import Image

from slope_tables import write_json


def normalized_graphics(page):
    """Resolve Matplotlib's page-local graphics-state and marker aliases."""
    resources = page['/Resources']
    def state(match):
        return json.dumps(dict(resources['/ExtGState'][match[1]].get_object()), sort_keys=True)+' gs'
    def marker(match):
        obj = resources['/XObject'][match[1]].get_object()
        metadata = {key: value for key, value in obj.items() if key not in ['/Length', '/Filter']}
        return json.dumps(metadata, sort_keys=True, default=str)+' '+obj.get_data().decode()+' Do'
    content = page.get_contents().get_data().decode()
    content = re.sub(r'(/A\d+) gs', state, content)
    return re.sub(r'(/M\d+) Do', marker, content)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--figures', type=Path, required=True)
    parser.add_argument('--pdf-tools', type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.pdf_tools))
    from pypdf import PdfReader
    folder = args.figures
    index = json.loads((folder/'figure_index.json').read_text())
    counts = json.loads((folder/'bin_counts.json').read_text())
    data = json.loads((args.root/'gun_analysis_v2/summary.json').read_text())
    checks = 0
    modes = ['beforeGeo', 'beforeGeo_regular_both', 'afterGeo']
    for dsid in ['100116', '100117']:
        for momentum in [100, 1000]:
            groups = [data['groups'][f'{dsid}|{momentum}|{mode}'] for mode in modes]
            for metric in ['conditional_misid', 'truth_correct_efficiency']:
                assert counts[f'{dsid}|{momentum}|{metric}'] == [group[metric] for group in groups]
                checks += 1
            assert counts[f'{dsid}|{momentum}|angular_summary'] == [group['angular_histograms']['residual_summary'] for group in groups]
            checks += 1
            assert counts[f'{dsid}|{momentum}|refit_counts'] == [group['refit_counts'] for group in groups]
            checks += 1
            assert counts[f'{dsid}|{momentum}|observed_wall_ms'] == [1000*group['reconstruction_wall_seconds_per_source_event'] for group in groups]
            checks += 1
            for mode, group in zip(modes, groups):
                for observable in ['signed_qop', 'relative_p', 'qop_pull', 'st1_sy_minus_launch_fine', 'st1_sy_pull']:
                    family = 'angular_histograms' if observable.startswith('st1_') else 'histograms'
                    assert counts[f'{dsid}|{momentum}|{mode}|{observable}'] == group[family][observable]
                    checks += 1
                regions = ['xminus_yminus', 'xminus_yplus', 'xplus_yminus', 'xplus_yplus']
                assert counts[f'{dsid}|{momentum}|{mode}|production_regions'] == [group['production_regions'][region]['conditional_misid'] for region in regions]
                checks += 1
    assert counts['charge_differences'] == data['charge_differences']
    assert counts['migrations'] == data['migrations']
    assert counts['nominal_launch_precision'] == data['nominal_launch_precision']
    pdf = PdfReader(folder/'charge_gun_report.pdf')
    assert len(pdf.pages) == len(index) == 11
    pages = []
    for page, row in zip(pdf.pages, index):
        with Image.open(folder/row['png']) as image:
            image.verify()
        single = PdfReader(folder/row['pdf'])
        assert len(single.pages) == 1
        original = single.pages[0]
        assert normalized_graphics(page) == normalized_graphics(original)
        assert (list(page.mediabox), list(page.cropbox), page.rotation) == (list(original.mediabox), list(original.cropbox), original.rotation)
        pages.append(dict(figure=row['name'], vector_operations_equal_with_resource_aliases_resolved=True,
                          PNG_integrity=True, page_boxes_equal=True))
    write_json(folder/'validation.json', dict(figure_count=11, summary_comparisons=checks,
        all_plotted_data_exact=True, pdf_pages=pages,
        visual='Pending independent inspection of exported images; numeric/PDF checks do not replace it'))
    print('Verified', checks, 'plot-summary comparisons and 11 vector PDF pages')


if __name__ == '__main__':
    main()
