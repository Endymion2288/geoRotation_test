#!/usr/bin/env python3
"""Assemble the measured S05 reports and a checked figure/source index."""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from PIL import Image

from slope_tables import write_json


def pages(path):
    info = subprocess.check_output(['pdfinfo', str(path)], text=True)
    return int(re.search(r'Pages:\s+(\d+)', info).group(1))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--catalog', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    catalog = json.loads(args.catalog.read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    reports, figures = [], []
    page = 1
    for item in catalog['reports']:
        report = args.root/item['path']
        n = pages(report)
        pngs = sorted(report.parent.glob('*.png'))
        assert len(pngs) == n, (report, len(pngs), n)
        rows = []
        for png in pngs:
            with Image.open(png) as image:
                image.verify()
            pdf = png.with_suffix('.pdf')
            assert pages(pdf) == 1, pdf
            row = dict(name=png.stem, png=str(png), pdf=str(pdf), report=str(report),
                       source=item['source'], reference=item.get('reference'), png_integrity=True)
            figures.append(row)
            rows.append(row)
        if item.get('index'):
            original = json.loads((args.root/item['index']).read_text())
            assert original['figure_count'] == n
            assert {r['png'] for r in original['figures']} == {r['png'] for r in rows}
            for row in rows:
                match = next(r for r in original['figures'] if r['png'] == row['png'])
                row['reference'] = match['reference']
                row['selection'] = match['selection']
                row['definition'] = match['definition']
        reports.append(dict(**item, pdf=str(report), pages=n, combined_first_page=page,
                            combined_last_page=page+n-1))
        page += n
    evidence = []
    for item in catalog['evidence']:
        source = args.root/item
        with source.open() as stream:
            data = json.load(stream)
        evidence.append(dict(path=str(source), valid_json=True, top_level_keys=list(data)))
    destination = args.output/'stage05_analysis_report.pdf'
    with tempfile.TemporaryDirectory(prefix='s05_report_') as scratch:
        local = Path(scratch)/destination.name
        subprocess.run(['pdfunite', *[r['pdf'] for r in reports], str(local)], check=True)
        assert pages(local) == page-1
        shutil.copyfile(local, destination)
    assert pages(destination) == page-1
    write_json(args.output/'index.json', dict(report=str(destination), pages=page-1,
        figure_count=len(figures), reports=reports, figures=figures, evidence=evidence,
        validation='All required evidence JSON parses; every indexed PNG opens; every single PDF has one page; each report and combined page count closes. Physical conclusions remain in the workbook.'))
    lines = ['# S05 图集与证据索引', '', f'[完整图集]({destination})，{page-1}页，{len(figures)}张图。', '',
             '| 图组 | 汇编页码 | 输入/统计来源 |', '| --- | --- | --- |']
    for row in reports:
        lines.append(f"| [{row['name']}]({row['pdf']}) | {row['combined_first_page']}–{row['combined_last_page']} | {row['source']} |")
    lines += ['', '逐图PNG/PDF、参考页、选择和完整来源保存在同目录index.json；物理结论与未解释余量见仓库workbook。', '']
    (args.output/'index.md').write_text('\n'.join(lines))
    print(json.dumps(dict(report=str(destination), pages=page-1, figures=len(figures))))


if __name__ == '__main__':
    main()
