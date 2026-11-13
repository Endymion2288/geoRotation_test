#!/usr/bin/env python3
"""Compare the assembled S05 atlas to every original report and single PDF."""
import argparse
import json
from pathlib import Path
import sys

from PIL import Image


def page_values(page):
    return (page.get_contents().get_data(), page.extract_text(),
            tuple(page.mediabox), tuple(page.cropbox), page.rotation)


def resolved_commands(page):
    """Resolve per-document resource names before comparing single/group pages."""
    from pypdf.generic import ArrayObject, DictionaryObject, NameObject, StreamObject

    def value(obj):
        obj = obj.get_object() if hasattr(obj, 'get_object') else obj
        if isinstance(obj, StreamObject):
            return ('stream', obj.get_data(), tuple(sorted(
                (str(k), value(v)) for k, v in obj.items()
                if k not in ['/Length', '/Filter', '/DecodeParms'])))
        if isinstance(obj, DictionaryObject):
            return tuple(sorted((str(k), value(v)) for k, v in obj.items()))
        if isinstance(obj, (ArrayObject, list)):
            return tuple(value(v) for v in obj)
        return str(obj)

    resources = page['/Resources']
    names = {}
    for category in ['/Font', '/ExtGState', '/XObject', '/Pattern', '/Shading']:
        for name, resource in resources.get(category, {}).items():
            resource = resource.get_object()
            if category == '/Font':
                # PdfPages embeds a superset of the same font's glyphs.  Its
                # subset size and /F ordinal are not drawing differences.
                descriptor = (str(resource.get('/BaseFont')), str(resource.get('/Subtype')))
            else:
                descriptor = value(resource)
            names[str(name)] = (category, descriptor)

    def operand(obj):
        if isinstance(obj, NameObject) and str(obj) in names:
            return names[str(obj)]
        if isinstance(obj, (ArrayObject, list)):
            return tuple(operand(v) for v in obj)
        return value(obj)

    return tuple((tuple(operand(v) for v in operands), operator)
                 for operands, operator in page.get_contents().operations)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', type=Path, required=True)
    parser.add_argument('--catalog', type=Path, required=True)
    parser.add_argument('--pdf-tools', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.pdf_tools))
    from pypdf import PdfReader

    index = json.loads(args.index.read_text())
    catalog = json.loads(args.catalog.read_text())
    combined = PdfReader(index['report'])
    assert len(combined.pages) == index['pages'] == index['figure_count']
    assert len(index['reports']) == len(catalog['reports'])
    assert len(index['evidence']) == len(catalog['evidence'])
    bookmarks = [item for item in combined.outline if not isinstance(item, list)]
    assert len(bookmarks) == len(index['reports'])

    reports = []
    next_page = 1
    for item, source_item, bookmark in zip(index['reports'], catalog['reports'], bookmarks):
        for key, value in source_item.items():
            assert item[key] == value, (item['name'], key)
        assert item['combined_first_page'] == next_page
        source = PdfReader(item['pdf'])
        assert len(source.pages) == item['pages']
        assert item['combined_last_page'] == next_page + len(source.pages) - 1
        assert bookmark.title == item['name']
        assert combined.get_destination_page_number(bookmark) == next_page - 1
        original = [page_values(page) for page in source.pages]
        normalized = [resolved_commands(page) for page in source.pages]
        for i, values in enumerate(original):
            assert values == page_values(combined.pages[next_page - 1 + i]), (item['name'], i)

        figures = [row for row in index['figures'] if row['report'] == item['pdf']]
        assert len(figures) == len(source.pages)
        remaining = list(range(len(original)))
        for figure in figures:
            with Image.open(figure['png']) as png:
                png.verify()
            with Image.open(figure['png']) as png:
                png.load()
            single = PdfReader(figure['pdf'])
            assert len(single.pages) == 1
            values = page_values(single.pages[0])
            matches = [i for i in remaining if values == original[i]]
            if not matches:
                commands = resolved_commands(single.pages[0])
                matches = [i for i in remaining if values[1:] == original[i][1:]
                           and commands == normalized[i]]
            assert matches, (figure['pdf'], 'single PDF differs from report pages')
            matched = matches[0]
            remaining.remove(matched)
            figure['combined_page'] = next_page + matched
        assert not remaining
        reports.append(dict(name=item['name'], pages=len(original),
                            first_page=next_page, last_page=item['combined_last_page'],
                            decoded_stream_text_boxes_rotation_equal=True,
                            single_pdf_resource_definitions_and_commands_equal=True,
                            single_pdfs_and_pngs_checked=len(figures), bookmark_equal=True))
        next_page += len(original)

    assert next_page - 1 == len(combined.pages)
    assert len({row['combined_page'] for row in index['figures']}) == index['pages']
    for item in index['evidence']:
        data = json.loads(Path(item['path']).read_text())
        assert list(data) == item['top_level_keys']
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(
        index=str(args.index), catalog=str(args.catalog), report=index['report'],
        reports=reports, pages=index['pages'], figure_count=index['figure_count'],
        evidence_files=len(index['evidence']), all_pages_exactly_equal_to_sources=True,
        pngs_fully_decoded=True, bookmarks_and_index_ranges_equal=True,
        figures=[dict(name=row['name'], combined_page=row['combined_page'],
                      png=row['png'], pdf=row['pdf']) for row in index['figures']],
        physics_validation='Counts, selections, flows and uncertainties are validated by the catalog evidence and the A1-A5 workbook audit.'
    ), indent=2) + '\n')
    print(json.dumps(dict(reports=len(reports), pages=index['pages'],
                          figures=index['figure_count'], evidence=len(index['evidence']))))


if __name__ == '__main__':
    main()
