#!/usr/bin/env python3
"""Check refit-control figure counts, exports, and unchanged measured values."""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess

import numpy as np
from PIL import Image


def numerical_closure(previous, current, path=''):
    if isinstance(previous, dict):
        return sum(numerical_closure(value, current[key], path + '/' + key)
                   for key, value in previous.items())
    if isinstance(previous, list):
        assert len(previous) == len(current), path
        return sum(numerical_closure(a, b, path + '/' + str(i))
                   for i, (a, b) in enumerate(zip(previous, current)))
    if isinstance(previous, (int, float)):
        assert previous == current or (isinstance(previous, float)
                                      and math.isnan(previous) and math.isnan(current)), path
        return 1
    return 0


def histograms(item, path=''):
    if isinstance(item, dict):
        if {'edges', 'counts', 'sumw', 'sumw2', 'raw_entries'} <= item.keys():
            yield path, item
        else:
            for key, value in item.items():
                yield from histograms(value, path + '/' + key)
    elif isinstance(item, list):
        for i, value in enumerate(item):
            yield from histograms(value, path + '/' + str(i))


def pages(path):
    text = subprocess.check_output(['pdfinfo', str(path)], text=True)
    return int(re.search(r'Pages:\s+(\d+)', text).group(1))


def page_text(path, page=1):
    return subprocess.check_output(
        ['pdftotext', '-f', str(page), '-l', str(page), str(path), '-'], text=True)


def resolve_panel(case, item):
    if 'reference_path' not in item:
        return item
    value = case
    for part in item['reference_path'].split('/'):
        value = value[part]
    return value


def validate(path, previous=None):
    data = json.loads(path.read_text())
    control = data['control']
    folder = path.parent
    report = folder / (control + '_refit_report.pdf')
    expected = 3 * len(data['cases'])
    assert pages(report) == expected, report
    checked = 0
    for name, h in histograms(data):
        edges, counts = np.array(h['edges']), np.array(h['counts'])
        assert len(edges) == len(counts) + 1 and np.all(np.diff(edges) > 0), name
        assert (counts >= 0).all(), name
        np.testing.assert_array_equal(counts, h['sumw'], err_msg=name)
        np.testing.assert_array_equal(counts, h['sumw2'], err_msg=name)
        assert int(counts.sum()) + h['underflow'] + h['overflow'] == h['finite'], name
        assert h['finite'] + h['invalid'] == h['raw_entries'], name
        checked += 1
    assert checked == (54 if control.startswith('regular_') else 48) * len(data['cases']), checked
    page = 0
    for case in data['cases'].values():
        for plane in ['st1', 'VetoNu', 'sensor']:
            panels = case['figure_histograms'][plane]['panels']
            for mode in ['native', control]:
                full = resolve_panel(case, panels['full'][mode])
                for view in ['zoom', 'residual']:
                    assert panels[view][mode]['raw_entries'] == full['raw_entries']
            if control.startswith('regular_'):
                for view in ['full', 'zoom']:
                    assert panels[view]['reference']['raw_entries'] == full['raw_entries']
            stem = case['label'] + '_' + plane + '_' + control + '_refit'
            png, pdf = (folder / (stem + suffix) for suffix in ['.png', '.pdf'])
            with Image.open(png) as bitmap:
                bitmap.verify()
            with Image.open(png) as bitmap:
                bitmap.load()
            assert pages(pdf) == 1, pdf
            page += 1
            assert page_text(pdf) == page_text(report, page), (report, page)
    assert len(list(folder.glob('*.png'))) == expected
    result = dict(summary=str(path), cases=len(data['cases']), histograms=checked,
                  pngs=expected, report_pages=expected, png_integrity=True,
                  counts_flows_invalid_sumw2_consistent=True,
                  all_figure_panels_have_counts=True, single_report_page_text_equal=True)
    if previous is not None:
        result['previous_summary'] = str(previous)
        result['original_numerical_fields_exactly_equal'] = numerical_closure(
            json.loads(previous.read_text()), data)
    return result


def main():
    parser = argparse.ArgumentParser()
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--summary', type=Path)
    inputs.add_argument('--comparisons', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.comparisons:
        results = []
        for row in json.loads(args.comparisons.read_text()):
            assert (Path(row['job_output']) / 'exit_code.txt').read_text().strip() == '0'
            results.append(validate(Path(row['new_summary']), Path(row['source_summary'])))
    else:
        results = [validate(args.summary)]
    args.output.write_text(json.dumps(dict(reports=results, histograms=sum(
        r['histograms'] for r in results), figure_count=sum(r['pngs'] for r in results),
        original_numerical_fields_exactly_equal=sum(
            r.get('original_numerical_fields_exactly_equal', 0) for r in results)), indent=2) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
