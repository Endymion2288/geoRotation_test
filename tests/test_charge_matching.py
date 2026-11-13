import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
from charge_tables import rank_candidates, weighted_rate
from charge_analysis import paired_difference
from charge_lineage_analysis import classify_path


class ChargeMatchingTest(unittest.TestCase):
    def test_candidate_choice_does_not_favour_correct_charge_or_momentum(self):
        tracks=dict(purity=np.array([1.,1.,1.]),n_hits=np.array([18,17,18]),
            chi2=np.array([5.,1.,6.]),ordinal=np.array([0,1,2]),
            charge=np.array([-1,1,1]),momentum_GeV=np.array([900.,100.,100.]))
        chosen,duplicates=rank_candidates(np.array([0,0,0]),np.array([True]*3),tracks,1)
        self.assertEqual(chosen.tolist(),[0]);self.assertEqual(duplicates.tolist(),[3])
        tracks['charge']*=-1;tracks['momentum_GeV']*=20
        self.assertEqual(rank_candidates(np.array([0,0,0]),np.array([True]*3),tracks,1)[0].tolist(),[0])

    def test_zero_wrong_charge_retains_a_nonzero_upper_bound(self):
        result=weighted_rate(np.zeros(50,dtype=bool),np.ones(50,dtype=bool),np.ones(50))
        self.assertEqual(result['effective_n'],50)
        self.assertEqual(result['value'],0)
        self.assertGreater(result['interval95'][1],0)

    def test_paired_rotation_difference_uses_shared_event_covariance(self):
        before=dict(selected=np.ones(6,dtype=bool),wrong=np.array([1,0,0,0,0,0],dtype=bool))
        same=paired_difference(before,before,np.ones(6,dtype=bool))
        self.assertEqual(same['value'],0.)
        self.assertEqual(same['standard_error'],0.)
        after=dict(selected=np.ones(6,dtype=bool),wrong=np.array([0,1,0,0,0,0],dtype=bool))
        exchanged=paired_difference(before,after,np.ones(6,dtype=bool))
        self.assertEqual(exchanged['value'],0.)
        self.assertGreater(exchanged['standard_error'],0.)

    def test_seed_recovery_and_later_refit_error_are_distinguished(self):
        first,final=classify_path([1,2,3,4,5,6],np.array([-.01,.01,.01,.01,.01,-.01]),1)
        self.assertEqual(first,1)
        self.assertEqual(final,6)

    def test_last_forward_error_can_recover_before_native_reference(self):
        first,final=classify_path([1,2,3,4,5,6],np.array([.01,.01,-.01,.01,.01,.01]),1)
        self.assertEqual(first,3)
        self.assertEqual(final,0)


if __name__=='__main__':unittest.main()
