#!/usr/bin/env python3
"""Assemble independently indexed vector charge figures and validate exports."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import tempfile

from PIL import Image

from slope_tables import write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--figures', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--pdf-tools', type=Path, help='Recorded existing pypdf installation')
    args = parser.parse_args()
    if args.pdf_tools:
        sys.path.insert(0, str(args.pdf_tools))
    from pypdf import PdfReader, PdfWriter
    args.output.mkdir(parents=True, exist_ok=False)
    writer, index, checks, contents = PdfWriter(), [], [], []
    for folder in args.figures:
        rows = json.loads((folder/'figure_index.json').read_text())
        assert isinstance(rows, list)
        for row in rows:
            png, pdf = folder/row['png'], folder/row['pdf']
            with Image.open(png) as image:
                image.verify()
            reader = PdfReader(pdf)
            assert len(reader.pages) == 1
            page = reader.pages[0]
            contents.append((page.get_contents().get_data(), list(page.mediabox), list(page.cropbox), page.rotation))
            writer.add_page(page)
            index.append(dict(row, png=str(png), pdf=str(pdf), report_page=len(index)+1))
            checks.append(dict(figure=row['name'], png_valid=True, vector_pdf_pages=1))
    target = args.output/'charge_stage06_report.pdf'
    with tempfile.TemporaryDirectory(prefix='charge-report-') as directory:
        local = Path(directory)/target.name
        writer.write(local)
        merged = PdfReader(local)
        assert len(merged.pages) == len(contents)
        for page, original in zip(merged.pages, contents):
            assert (page.get_contents().get_data(), list(page.mediabox), list(page.cropbox), page.rotation) == original
        shutil.copy2(local, target)
    write_json(args.output/'figure_index.json', index)
    write_json(args.output/'validation.json', dict(figure_count=len(index), pdf_pages=len(contents),
        merged_vector_contents_and_boxes_equal=True, PNG_integrity=checks,
        visual_inspection='Pending separate inspection of main figures and the contact sheet'))
    (args.output/'figure_index.md').write_text('\n'.join(
        f"- Page {r['report_page']}: [{r['name']}]({r['png']}) · [PDF]({r['pdf']}) · {r['reference']}" for r in index)+'\n')
    print(json.dumps(dict(report=str(target), pages=len(index)), indent=2))


if __name__ == '__main__':
    main()
