"""Boundary regression for the explicitly requested strict final purity gate."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from native_purity import passes_purity


class FinalPurityTest(unittest.TestCase):
    def test_strict_boundary(self):
        self.assertFalse(passes_purity(7, 10))
        self.assertTrue(passes_purity(13, 18))
        self.assertFalse(passes_purity(12, 18))

    def test_missing_associations_stay_in_denominator(self):
        self.assertFalse(passes_purity(0, 0))
        self.assertFalse(passes_purity(0, 18))
        self.assertFalse(passes_purity(7, 18))


if __name__ == '__main__':
    unittest.main()
