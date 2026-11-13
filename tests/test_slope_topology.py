import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from slope_topology import pair_topology


class TopologyTest(unittest.TestCase):
    def truth(self):
        return dict(sample_entry=np.array([0, 1, 0, 1]), barcode=np.array([2, 2, 3, 3]),
            pdg=np.array([11, 11, -11, -11]), status=np.ones(4, dtype=int),
            production_position_F_mm=np.array([[0, 0, -100], [0, 0, -500], [0, 0, -100], [0, 0, -500]]),
            production_momentum_F_GeV=np.array([[0, 0, 100], [0, 0, 50], [0, 0, 50], [50, 0, 0]]),
            truth_station_position_F_mm=np.zeros((4, 4, 3)), truth_hit_counts=np.ones((4, 4), dtype=int))

    def test_source_identity_and_opening_units(self):
        truth = self.truth()
        truth['truth_station_position_F_mm'][2, 1, :2] = [3, 4]
        fields, counts = pair_topology(truth)
        self.assertEqual(counts['valid_pairs'], 2)
        np.testing.assert_allclose(fields['opening_mrad'][[0, 2]], 0)
        np.testing.assert_allclose(fields['opening_mrad'][[1, 3]], np.pi*500)
        np.testing.assert_allclose(fields['st1_centroid_separation_mm'][[0, 2]], 5)

    def test_missing_hits_do_not_become_zero_separation(self):
        truth = self.truth()
        truth['truth_hit_counts'][2, 1] = 0
        fields, counts = pair_topology(truth)
        self.assertTrue(np.isnan(fields['st1_centroid_separation_mm'][[0, 2]]).all())
        self.assertEqual(counts['missing_st1_centroids'], 1)
        self.assertEqual(counts['valid_pairs'], 2)

    def test_different_vertices_are_not_a_decay_pair(self):
        truth = self.truth()
        truth['production_position_F_mm'][2, 2] = -101
        fields, counts = pair_topology(truth)
        self.assertEqual(counts['different_vertices'], 1)
        self.assertTrue(np.isnan(fields['opening_mrad'][[0, 2]]).all())


if __name__ == '__main__':
    unittest.main()
