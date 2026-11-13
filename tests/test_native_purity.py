"""Boundary regression for the explicitly requested strict final purity gate."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from native_purity import passes_purity, seed_cuts


class FinalPurityTest(unittest.TestCase):
    def test_strict_boundary(self):
        self.assertFalse(passes_purity(7, 10))
        self.assertTrue(passes_purity(13, 18))
        self.assertFalse(passes_purity(12, 18))

    def test_missing_associations_stay_in_denominator(self):
        self.assertFalse(passes_purity(0, 0))
        self.assertFalse(passes_purity(0, 18))
        self.assertFalse(passes_purity(7, 18))


class SeedCutTest(unittest.TestCase):
    def test_p_magnitude_and_charge_are_separate(self):
        self.assertEqual(seed_cuts(.01, .01, .2, .2), (True, True, True))
        self.assertEqual(seed_cuts(-.01, .01, .2, .2), (True, False, True))
        self.assertFalse(seed_cuts(.01, .012, .2, .2)[0])

    def test_theta_units_and_invalid_seed(self):
        self.assertTrue(seed_cuts(.01, .01, 0., .000099)[2])
        self.assertFalse(seed_cuts(.01, .01, 0., .0001)[2])
        self.assertEqual(seed_cuts(0., .01, .2, .2), (False, False, False))
        self.assertEqual(seed_cuts(float('nan'), .01, .2, .2), (False, False, False))


if __name__ == '__main__':
    unittest.main()
