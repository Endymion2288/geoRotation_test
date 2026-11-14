import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from charge_joint_regression import cluster_rate, paired_rate


class ChargeJointRegressionTests(unittest.TestCase):
    def test_replicating_correlated_daughters_does_not_improve_precision(self):
        event = np.arange(4)
        wrong = np.array([True, False, False, False])
        single = cluster_rate(event, wrong, np.ones(4, dtype=bool), 4)
        paired = cluster_rate(np.repeat(event, 2), np.repeat(wrong, 2), np.ones(8, dtype=bool), 4)
        self.assertEqual(single['value'], paired['value'])
        self.assertEqual(single['standard_error'], paired['standard_error'])
        self.assertAlmostEqual(paired['standard_error'], .25)

    def test_no_wrong_tracks_has_a_nonzero_upper_bound(self):
        result = cluster_rate(np.repeat(np.arange(20), 2), np.zeros(40, dtype=bool),
                              np.ones(40, dtype=bool), 20)
        self.assertEqual(result['value'], 0.)
        self.assertGreater(result['interval95'][1], 0.)

    def test_exactly_paired_changes_cancel_event_fluctuations(self):
        event = np.repeat(np.arange(4), 2)
        wrong = np.array([True, True, False, False, False, False, False, False])
        selected = np.ones(8, dtype=bool)
        identical = paired_rate(event, wrong, selected, wrong, selected, 4)
        self.assertEqual(identical['value'], 0.)
        self.assertEqual(identical['standard_error'], 0.)
        self.assertLess(identical['interval95'][0], 0.)
        self.assertGreater(identical['interval95'][1], 0.)
        moved = np.roll(wrong, 2)
        shifted = paired_rate(event, wrong, selected, moved, selected, 4)
        self.assertEqual(shifted['value'], 0.)
        self.assertGreater(shifted['standard_error'], 0.)


if __name__ == '__main__':
    unittest.main()
