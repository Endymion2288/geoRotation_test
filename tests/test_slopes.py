import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from slope_tables import angle, event_first_mask, histogram, join_truth, physical, slope
from slope_statistics import conditional_response, event_ratio_interval


class SlopeDefinitions(unittest.TestCase):
    def test_frame_near_beam_and_sign(self):
        f = np.array([[0.,0.,100.], [.01,.2,100.], [-.01,-.2,100.]])
        a = np.stack([f[:,2], f[:,1], -f[:,0]], axis=1)
        np.testing.assert_array_equal(physical(a, True), f)
        np.testing.assert_array_equal(slope(f), a[:,1]/a[:,0])
        np.testing.assert_allclose(angle(f), np.arctan(slope(f)), atol=1e-16)
        self.assertTrue(np.all(np.isfinite(slope(physical(a,True)))))

    def test_pdf_first_long_track_cut_is_event_cut(self):
        entries = np.array([0,0,0,1,1,2])
        pz = np.array([5.,100.,500.,300.,100.,70.])
        long = np.array([False,True,True,True,True,True])
        np.testing.assert_array_equal(event_first_mask(entries,pz,long),
                                      [False,True,True,False,False,True])

    def test_invalid_and_root_upper_flow(self):
        h = histogram(np.array([-2.,-1.,0.,1.,np.nan,np.inf]),np.array([-1.,0.,1.]))
        self.assertEqual(h['counts'],[1,1])
        self.assertEqual((h['underflow'],h['overflow'],h['invalid']),(1,1,2))
        self.assertEqual(h['raw_entries'],sum(h['counts'])+h['underflow']+h['overflow']+h['invalid'])

    def test_truth_identity_includes_input_entry(self):
        tables = dict(truth=dict(sample_entry=np.array([1,0]),barcode=np.array([4,4])),
                      tracks=dict(sample_entry=np.array([0,1,2]),barcode=np.array([4,4,4])))
        index,valid = join_truth(tables)
        np.testing.assert_array_equal(valid,[True,True,False])
        np.testing.assert_array_equal(index[valid],[1,0])

    def test_event_bootstrap_keeps_correlated_tracks_together(self):
        result = event_ratio_interval(np.array([0,0,1,1]),np.array([1,1,0,0]),np.ones(4),(0,2),repeats=2000)
        self.assertEqual(result['estimate'],.5)
        self.assertEqual(result['interval_95'],[0.,1.])

    def test_validation_reco_cannot_change_learned_response(self):
        n=1200
        entries=np.arange(n)
        true=np.full(n,.001)
        reco=true+.0003+2e-5*np.sin(entries)
        args=[true,np.full(n,100.),np.zeros(n),np.zeros(n),entries,reco,
              entries,np.ones(n,dtype=bool),entries,600,np.linspace(-.01,.01,101)]
        windows=[dict(center=.0013,half_width=.0001)]
        first=conditional_response(*args,windows=windows)
        args[5]=reco.copy()
        args[5][600:]+=.001
        second=conditional_response(*args,windows=windows)
        self.assertEqual(first['predicted_counts'],second['predicted_counts'])
        self.assertEqual(first['predicted_counts_interval_95'],second['predicted_counts_interval_95'])
        self.assertGreater(max(first['observed_counts_interval_95'][0]),500)
        self.assertGreater(first['peak_comparison'][0]['predicted_std'],1e-5)
        self.assertGreater(first['peak_comparison'][0]['observed_std'],1e-5)
        self.assertNotEqual(first['observed_counts'],second['observed_counts'])

    def test_response_compares_only_the_same_covered_truth_cells(self):
        n=1200;entries=np.arange(n);true=np.full(n,.001)
        momentum=np.full(n,100.);momentum[900:]=10.
        result=conditional_response(true,momentum,np.zeros(n),np.zeros(n),entries,
            true,entries,np.ones(n,dtype=bool),entries,600,np.linspace(-.01,.01,101),
            bootstrap_repeats=10,windows=[dict(center=.001,half_width=.0001)])
        self.assertEqual(result['validation_reco'],600)
        self.assertEqual(result['covered_validation_reco'],300)
        self.assertEqual(result['uncovered_validation_reco'],300)
        np.testing.assert_array_equal(np.array(result['covered_observed_counts'])+
            result['uncovered_validation_reco_counts'],result['observed_counts'])
        self.assertEqual(result['covered_peak_comparison'][0]['observed_entries'],300)
        self.assertEqual(result['covered_peak_comparison'][0]['uncovered_observed_entries'],300)

    def test_common_validation_support_does_not_restrict_training(self):
        entries = np.arange(1200)
        true = np.full(1200, .001)
        common = entries % 2 == 0
        result = conditional_response(true, np.full(1200, 100.), np.zeros(1200), np.zeros(1200),
            entries, true+.0003, entries, np.ones(1200, bool), entries, 600,
            np.linspace(-.01, .01, 101), phase_period=.08/63., conditioning={'phase_bins': 4},
            validation_truth_mask=common, bootstrap_repeats=10)
        self.assertEqual(result['discovery_reco'], 600)
        self.assertEqual(result['validation_reco'], 300)
        self.assertEqual(result['covered_validation_reco'], 300)
        self.assertEqual(result['cells'][0]['training_truth'], 600)
        self.assertEqual(result['cells'][0]['training_reco'], 600)
        self.assertEqual(sum(result['predicted_counts']), 300)


if __name__ == '__main__':
    unittest.main()
