import math
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from charge_posterior import component_path, event_posterior, innovation, mixture, normalized_weights


class NativeLikelihoodAnalysis(unittest.TestCase):
    def test_likelihood_keeps_determinant_and_normalization(self):
        row = dict(measurement=[2.], projector=[[1., 0.]], predicted=[1., 8.],
            predicted_cov=[[3., 1.], [1., 2.]], measurement_cov=[[1.]])
        value = innovation(row)
        self.assertAlmostEqual(value['chi2'], .25)
        self.assertAlmostEqual(value['log_likelihood'], -.5*(.25+math.log(4.)+math.log(2*math.pi)))
        row['measurement_cov'] = [[13.]]
        self.assertLess(innovation(row)['log_likelihood'], value['log_likelihood'])

    def test_log_weights_remain_normalized_with_large_negative_log_likelihoods(self):
        weights, evidence = normalized_weights([-10000., -10001.])
        self.assertAlmostEqual(float(weights.sum()), 1.)
        self.assertAlmostEqual(weights[0], 1/(1+math.exp(-1)))
        self.assertAlmostEqual(evidence, -10000+math.log1p(math.exp(-1))-math.log(2))

    def test_mixture_includes_between_component_covariance(self):
        value = mixture([[-1.], [3.]], [[[2.]], [[4.]]], [.25, .75])
        self.assertAlmostEqual(value['mean'][0], 2.)
        self.assertAlmostEqual(value['within_cov'][0][0], 3.5)
        self.assertAlmostEqual(value['between_cov'][0][0], 3.)
        self.assertAlmostEqual(value['cov'][0][0], 6.5)

    def test_missing_component_has_no_conditional_posterior(self):
        value = event_posterior([dict(complete=True), dict(complete=False)], 1)
        self.assertEqual(value['input_components'], 2)
        self.assertEqual(value['missing_prior_mass'], .5)
        self.assertFalse(value['complete'])
        self.assertIsNone(value['posterior'])

    def test_wrong_sign_probability_and_refit_fallback_are_both_kept(self):
        row = dict(complete=True, log_likelihood=-2., smoothed_mean=[0., 0., 0., 0., -1., 0.],
            smoothed_cov=np.eye(6).tolist(), final_qop_GeVinv=-2.,
            final_qop_variance_GeVinv2=2., refit_success=False)
        value = event_posterior([row], 1)
        self.assertAlmostEqual(value['posterior']['wrong_charge_probability'], .5*math.erfc(-1/math.sqrt(2)))
        self.assertEqual(value['native_final_output_mixture']['refit_success_CKF_posterior_weight'], 0.)
        self.assertEqual(row['posterior_weight'], 1.)

    def test_multiple_paths_are_not_replaced_by_a_first_candidate(self):
        value, measurements = component_path(dict(selected={1:{}, 2:{}}), [11], 'afterGeo')
        self.assertFalse(value['complete'])
        self.assertEqual(value['reason'], 'multiple_selected')
        self.assertIsNone(measurements)

    def test_same_hit_set_in_another_order_does_not_define_the_same_observations(self):
        states = [dict(acts_state_index=22, acts_previous_index=11, hit=1),
                  dict(acts_state_index=11, acts_previous_index=4294967295, hit=2)]
        for row in states:
            row.update(measurement_flag=True, outlier_flag=False, hole_flag=False)
        item = dict(selected={7:{}}, states={7:states}, final={})
        value, measurements = component_path(item, [1, 2], 'afterGeo')
        self.assertEqual(value['actual_hits'], [2, 1])
        self.assertEqual(value['reason'], 'different_observation_path')
        self.assertFalse(value['complete'])
        self.assertIsNone(measurements)

    def test_invalid_innovation_is_not_repaired_by_a_covariance_floor(self):
        row = dict(measurement=[1.], projector=[[1.]], predicted=[0.],
            predicted_cov=[[-2.]], measurement_cov=[[1.]])
        with self.assertRaises(np.linalg.LinAlgError):
            innovation(row)


if __name__ == '__main__':
    unittest.main()
