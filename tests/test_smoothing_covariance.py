"""Regression checks for physical left-tail thresholds and paired event errors."""
from collections import Counter
import math
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from reference_truth import accumulate
from smoothing_covariance import paired_rate


class SmoothingCovarianceTests(unittest.TestCase):
    def test_residual_units_tail_thresholds_and_invalid_reconstruction(self):
        counter = Counter()
        # truth p=1000 GeV -> truth q/p=1e-6 MeV^-1.
        for reco in [.001, .002, .01, 0., math.nan]:
            accumulate(counter, 'smt_native', reco, 1e-6)
        self.assertEqual(counter['smt_native_finite'], 3)
        self.assertEqual(counter['smt_native_below_minus_half'], 1)
        self.assertEqual(counter['smt_native_below_minus_point8'], 1)
        self.assertEqual(counter['smt_native_invalid_reco'], 2)

    def test_identical_paired_events_have_zero_difference_and_error(self):
        first = np.array([[4, 1, 0], [2, 2, 1], [3, 0, 0]])
        result = paired_rate(first, first.copy(), 1)
        self.assertEqual(result['delta_pp'], 0)
        self.assertEqual(result['lower95_pp'], 0)
        self.assertEqual(result['upper95_pp'], 0)

    def test_pairing_retains_event_covariance(self):
        first = np.array([[10, 0], [10, 10], [10, 0], [10, 10]])
        same = paired_rate(first, first, 1)
        shuffled = paired_rate(first, first[::-1], 1)
        self.assertEqual(same['first_rate'], .5)
        self.assertEqual(shuffled['delta_pp'], 0)
        self.assertGreater(shuffled['upper95_pp'], same['upper95_pp'])


if __name__ == '__main__':
    unittest.main()
