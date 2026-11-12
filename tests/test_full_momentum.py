import copy
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"analysis"))
from momentum_full_sample import add, empty, finalize, histogram, merge_values, metrics, step_metrics, verify_coverage


class FullMomentumChecks(unittest.TestCase):
    def test_full_coverage_rejects_duplicate_input_paths_gaps_and_overlaps(self):
        manifest = dict(samples=[dict(dsid="110051", files=[dict(path="first.root", entries=3), dict(path="second.root", entries=2)])])
        tasks = [dict(dsid="110051", file_index=0, input="first.root", skip=0, events=2),
                 dict(dsid="110051", file_index=0, input="first.root", skip=2, events=1),
                 dict(dsid="110051", file_index=1, input="second.root", skip=0, events=2)]
        verify_coverage(manifest, tasks)
        for field, value, index in [("input", "first.root", 2), ("skip", 1, 1), ("skip", 3, 1)]:
            broken = copy.deepcopy(tasks)
            broken[index][field] = value
            with self.subTest(field=field, value=value):
                with self.assertRaises(AssertionError):
                    verify_coverage(manifest, broken)

    def test_binning_at_edges_preserves_original_counts_and_overflow(self):
        edges = np.linspace(-1.2, 1.2, 101)
        values = np.concatenate([edges, np.nextafter(edges, -np.inf), np.nextafter(edges, np.inf), [-5., 10.]])
        h = histogram(edges)
        for x in values:
            add(h, x)
        expected = np.histogram(values[(values >= edges[0]) & (values < edges[-1])], edges)[0]
        self.assertEqual(h["counts"], expected.tolist())
        self.assertEqual(h["n"], sum(h["counts"])+h["underflow"]+h["overflow"])

    def test_signed_qop_tail_folding_and_recovery_use_same_valid_states(self):
        m = metrics()
        def bound(q):
            return [0., 0., 0., 1.57, q, 0.]
        for q in [.03, -.03]:
            step = dict(predicted=bound(q), filtered=bound(q), smoothed=bound(.01))
            step_metrics(m, step, 100., 1.)
        step_metrics(m, dict(predicted=bound(None), filtered=bound(.03), smoothed=bound(.01)), 100., 1.)
        self.assertEqual(m["counts"]["all_three_finite"], 2)
        self.assertEqual(m["counts"]["F_S_repaired"], 2)
        self.assertEqual(m["counts"]["missing_or_invalid_parameter"], 1)
        self.assertEqual(sum(m["masks"]), 2)

    def test_shard_merge_keeps_integer_counts_and_maximum_closure_error(self):
        left, right = empty("100116"), empty("100116")
        for v in ["current", "strict"]:
            for result, values in [(left, [-.8, 0., .2]), (right, [-.7, .1, 8.])]:
                h = result["variants"][v]["histograms"]["seed_candidates"]
                for value in values:
                    add(h, value)
            left["variants"][v]["metrics"]["good95"][0][0]["max_update_error_over_truth"] = 1e-9
            right["variants"][v]["metrics"]["good95"][0][0]["max_update_error_over_truth"] = 3e-9
        merged = finalize(merge_values(copy.deepcopy(left), right))
        finalize(merged)
        h = merged["variants"]["strict"]["histograms"]["seed_candidates"]
        self.assertEqual((h["n"], h["tail_half"], h["overflow"]), (6, 2, 1))
        self.assertEqual(merged["variants"]["strict"]["metrics"]["good95"][0][0]["max_update_error_over_truth"], 3e-9)
        self.assertNotIn("tail_half", merged["variants"]["strict"]["metrics"]["good95"][0][0]["sigma"]["filtered"])


if __name__ == "__main__":
    unittest.main()
