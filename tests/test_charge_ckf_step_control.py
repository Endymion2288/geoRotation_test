import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from charge_ckf_step_control import compare_path


def component(selected=True, final_qop=None):
    row = dict(candidate=3, measurements=18, chi2=2., hits=[11, 12],
               bound=[1., 2., .3, .003, .02, 0.])
    state = dict(acts_state_index=0, acts_previous_index=4294967295,
        path_length_mm=10., geometry_id=17, hit=11, measurement_flag=True,
        outlier_flag=False, hole_flag=False)
    for stage in ('predicted', 'filtered', 'smoothed'):
        state[stage] = row['bound']
        state[stage+'_cov'] = np.eye(6)
    result = dict(selected={3: row} if selected else {}, states={3:[state]},
        truth_charge=1, final={})
    if final_qop is not None:
        result['final'][3] = dict(bound_trk=[0., 0., .3, .003, final_qop/1000.], refit=0)
    return result


class InternalStepComparison(unittest.TestCase):
    def test_one_sided_selection_remains_in_input_and_failure_counts(self):
        rows = []
        counts = compare_path(component(final_qop=-.02), component(False),
            'beforeGeo', 'beforeGeo', (1, 2), {}, rows)
        self.assertEqual(counts['inputs'], 1)
        self.assertEqual(counts['one_sided_selected'], 1)
        self.assertEqual(counts['before_final_wrong'], 1)
        self.assertEqual(counts['before_refit_success'], 0)
        self.assertEqual(rows, [])

    def test_same_hits_with_changed_passive_topology_are_counted_without_false_state_pairing(self):
        a, b, rows = component(), component(), []
        b['states'][3][0]['geometry_id'] = 18
        counts = compare_path(a, b, 'afterGeo', 'afterGeo', (1, 2), {}, rows)
        self.assertEqual(counts['same_ordered_hits'], 1)
        self.assertEqual(counts['different_full_state_signatures'], 1)
        self.assertEqual(counts['paired_states'], 0)
        self.assertEqual(rows, [])

    def test_terminal_only_path_is_retained_and_identical_state_difference_is_zero(self):
        a, b, rows = component(final_qop=.02), component(final_qop=-.02), []
        counts = compare_path(a, b, 'beforeGeo', 'beforeGeo', (1, 2), {}, rows)
        self.assertEqual(counts['terminal_only'], 1)
        self.assertEqual(counts['paired_states'], 1)
        self.assertEqual(counts['final_sign_difference'], 1)
        self.assertEqual(rows[0]['smoothed']['direction_scaled_max'], 0.)

    def test_zero_sign_is_invalid_and_not_counted_as_wrong(self):
        counts = compare_path(component(final_qop=0.), component(False),
            'beforeGeo', 'beforeGeo', (1, 2), {}, [])
        self.assertEqual(counts['before_final_invalid'], 1)
        self.assertEqual(counts['before_final_wrong'], 0)


if __name__ == '__main__':
    unittest.main()
