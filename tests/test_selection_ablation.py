"""Legacy matching and quality ordering used for smoothing cut ablations."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from selection_ablation import FULL, final_selection


class SelectionAblationTests(unittest.TestCase):
    def setUp(self):
        self.truth = {1: dict(signal=True, status=1, pdg=11, vertex_F_mm=[0., 0., 0.],
            station_position_F_mm=[None, [0., 0., 0.], [0., 0., 1000.], [0., 0., 2000.]],
            hit_counts=[0, 6, 6, 6])}

    def track(self, candidate=1, **changes):
        row = dict(seed=0, candidate=candidate, matches=[[1, 18]], hits=list(range(18)),
            momentum_A_GeV=[100., 0., 0.], chi2=13., ndof=13, layers=9, rmax_mm=94.)
        row.update(changes)
        return row

    def test_final_purity_uses_actual_final_hit_denominator(self):
        low = self.track(matches=[[1, 12]])
        high = self.track(2, matches=[[1, 13]])
        self.assertEqual(final_selection([low, high], self.truth, {'purity'}), {(0, 2)})

    def test_unique_ranking_occurs_before_good95(self):
        best = self.track(rmax_mm=100.)
        lower = self.track(2, matches=[[1, 17]], hits=list(range(17)))
        self.assertEqual(final_selection([best, lower], self.truth, FULL), set())
        self.assertEqual(final_selection([best, lower], self.truth, FULL-{'unique'}), {(0, 2)})

    def test_removing_purity_recomputes_the_best_candidate(self):
        clean = self.track(matches=[[1, 12]], hits=list(range(15)))
        mixed = self.track(2, matches=[[1, 13]], hits=list(range(20)))
        self.assertEqual(final_selection([clean, mixed], self.truth, FULL), {(0, 1)})
        self.assertEqual(final_selection([clean, mixed], self.truth, FULL-{'purity'}), {(0, 2)})

    def test_quality_boundaries_are_original_good95(self):
        radius_boundary = self.track(rmax_mm=95.)
        chi2_boundary = self.track(2, chi2=25*13)
        ndof_boundary = self.track(3, ndof=6)
        self.assertEqual(final_selection([radius_boundary], self.truth, FULL), set())
        self.assertEqual(final_selection([chi2_boundary], self.truth, FULL), set())
        self.assertEqual(final_selection([ndof_boundary], self.truth, FULL), set())
        self.assertEqual(final_selection([self.track(ndof=7)], self.truth, FULL), {(0, 1)})

    def test_signal_and_hit_coverage_are_separate_filters(self):
        self.truth[1]['hit_counts'] = [0, 6, 6, 5]
        self.assertEqual(final_selection([self.track()], self.truth, FULL), set())
        self.assertEqual(final_selection([self.track()], self.truth, FULL-{'H18'}), {(0, 1)})
        self.truth[1]['signal'] = False
        self.assertEqual(final_selection([self.track()], self.truth, FULL-{'H18'}), set())


if __name__ == '__main__':
    unittest.main()
