"""Checks of matching independence, frame interpretation, and truth denominators."""
import math
from copy import deepcopy
import unittest

from analysis.pilot import candidate_rank, event_groups, good_track, match, pnorm, process_event, summarize


class PhysicsDefinitions(unittest.TestCase):
    def test_nonfinite_candidate_does_not_hide_valid_tie(self):
        base=dict(matches=[[10001,18]],hits=list(range(18)),candidate=0,
                  ndof=13,layers=9,rmax_mm=10)
        invalid=dict(base,chi2=None)
        valid=dict(base,chi2=25,candidate=1)
        self.assertIs(max([invalid,valid],key=candidate_rank),valid)
        self.assertFalse(good_track(invalid))
        self.assertTrue(good_track(valid))
        for value in [None,float('inf'),float('nan')]:
            self.assertIsNone(pnorm(dict(momentum_A_GeV=[1,value,0])))
        self.assertIsNone(pnorm(dict(momentum_A_GeV=[0.,0.,0.])))

    def test_reused_production_ids_remain_separate_source_events(self):
        records=[dict(schema=2,stage="event",event=0),dict(schema=2,stage="truth",event=0),
                 dict(schema=2,stage="event",event=0),dict(schema=2,stage="truth",event=0)]
        events=list(event_groups(records))
        self.assertEqual([len(e) for e in events],[2,2])

    def test_wrong_charge_and_large_tail_remain_matched(self):
        self.assertEqual(match({"matches": [[10001, 17]], "hits": list(range(18)),
                                "charge": -1, "momentum_A_GeV": [1,0,0]}),10001)
        self.assertIsNone(match({"matches": [[10001, 6]], "hits": list(range(18))}))

    def test_rotated_slope_and_missing_track_denominator(self):
        truth=dict(stage="truth",barcode=10001,signal=True,pdg=-13,status=1,
                   momentum_F_MeV=[0,300,150000],vertex_F_mm=[0,0,-100],
                   hit_counts=[0,6,6,6],station_position_F_mm=[None,[0,0,10],[0,0,1000],[0,0,2000]])
        track=dict(stage="final",matches=[[10001,18]],hits=list(range(18)),
                   layers=9,ndof=13,chi2=13,rmax_mm=10,momentum_A_GeV=[150,.3,0],
                   charge=1,seed=0,candidate=0,refit=1,measurements=18)
        found,counts=process_event((100116,0,1),[truth,track],"100116")
        missing,_=process_event((100116,0,2),[truth],"100116")
        self.assertAlmostEqual(found[0]["final"]["slope_y"],.002)
        self.assertIsNone(found[0]["final"]["legacy_angle"])
        result=summarize(found+missing)
        self.assertEqual(result["good95_efficiency"]["denominator"],2)
        self.assertEqual(result["good95_efficiency"]["numerator"],1)
        self.assertEqual(counts["final"],1)
        decayed=deepcopy(truth)
        decayed.update(barcode=10002,status=2)
        stable,_=process_event((100116,0,1),[truth,decayed,track],"100116")
        self.assertEqual([r["barcode"] for r in stable],[10001])

    def test_first_change_uses_final_candidate_ancestry(self):
        truth=dict(stage="truth",barcode=10001,signal=True,pdg=-13,status=1,
                   momentum_F_MeV=[0,0,150000],vertex_F_mm=[0,0,-100],
                   hit_counts=[0,6,6,6],station_position_F_mm=[None,[0,0,10],[0,0,1000],[0,0,2000]])
        final=dict(stage="final",matches=[[10001,18]],hits=list(range(18)),layers=9,
                   ndof=13,chi2=13,rmax_mm=10,momentum_A_GeV=[30,0,0],charge=1,
                   seed=1,candidate=1,refit=1,measurements=18)
        found=dict(final,stage="found")
        seed=dict(final,stage="seed",candidate=-1,momentum_A_GeV=[150,0,0])
        unrelated=dict(seed,seed=0,momentum_A_GeV=[1,0,0],hits=list(range(19)),matches=[[10001,19]])
        rows,_=process_event((100116,0,1),[truth,unrelated,seed,found,final],"100116")
        self.assertLess(rows[0]["seed"]["p_ratio"],.5)
        self.assertEqual(rows[0]["final_lineage"]["seed"]["p_ratio"],1.)
        self.assertEqual(summarize(rows)["good95_final_lineage"]["first_low_p_stage_for_final_tail"],{"found":1})


if __name__=="__main__":
    unittest.main()
