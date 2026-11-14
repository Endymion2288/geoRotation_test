import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from validate_charge_gun_figures import normalized_graphics


class Object(dict):
    def get_object(self):
        return self


class Page(dict):
    def __init__(self, alias, whitespace):
        super().__init__({'/Resources': {'/ExtGState': {alias: Object({'/ca': .4, '/CA': .4})}}})
        self.data = (alias+whitespace+'gs').encode()

    def get_contents(self):
        return self

    def get_data(self):
        return self.data


class PDFResourceAliases(unittest.TestCase):
    def test_equivalent_alpha_resource_aliases_across_line_breaks(self):
        self.assertEqual(normalized_graphics(Page('/A2', '\n')),
                         normalized_graphics(Page('/A5', ' ')))


if __name__ == '__main__':
    unittest.main()
