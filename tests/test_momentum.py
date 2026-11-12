import math
import gzip
from pathlib import Path
import tempfile
import unittest

import numpy as np

from analysis.seed_geometry import taubin
from analysis.momentum_simulation import sensor_bounds,t1_selection
from scripts.audit_reduce import reduce_event,trace_record


def state(index,qpred,qfiltered,qsmoothed):
    mean=lambda q:[0.,0.,0.,math.pi/2,q,0.]
    identity=np.eye(6).tolist()
    return dict(schema=2,run=1,event=1,lumi=0,collection='test',stage='selected_state',seed=0,candidate=0,
        state=index,geometry_id=index,acts_state_index=index,acts_previous_index=index-1,
        measurement_flag=False,outlier_flag=False,hole_flag=False,surface=np.eye(4).tolist(),
        predicted=mean(qpred),filtered=mean(qfiltered),smoothed=mean(qsmoothed),
        predicted_cov=identity,filtered_cov=identity,smoothed_cov=identity,jacobian=identity)


class MomentumChecks(unittest.TestCase):
    def test_sensor_truth_does_not_merge_repeated_g4_event_ids(self):
        header='event,source_entry,g4_event,track,parent,station,layer,row,module,sensor,z_mm,px_MeV,py_MeV,pz_MeV,post_px_MeV,post_py_MeV,post_pz_MeV\n'
        data='1,0,1,1,0,1,0,1,1,0,20,0,0,100000,0,0,90000\n2,1,1,1,0,1,0,1,1,0,20,0,0,1000000,0,0,800000\n'
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'truth.csv.gz'
            with gzip.open(path,'wt') as stream:stream.write(header+data)
            bounds=sensor_bounds(path)
            self.assertEqual(set(bounds),{1,2})
            self.assertEqual(bounds[1][0]['lower_GeV'],90.)
            self.assertEqual(bounds[2][0]['upper_GeV'],1000.)

    def test_good_seed_selection_keeps_failed_and_pruned_truth_in_denominator(self):
        final=lambda seed,p:dict(seed=seed,p_ratio=p,good95=True,wrong=False)
        rows=[dict(source_entry=i,barcode=1,fiducial=True,truth_H18=True,final=f)
              for i,f in enumerate([final(0,1.),None,final(9,.1),final(0,.1),final(0,1.)])]
        injections={(i,0):dict(event_ordinal=i,barcode=1,seed=0,eligible=True,good_seed_T1=i<4) for i in range(5)}
        result=t1_selection(rows,injections)
        self.assertEqual(result['mother_events'],4)
        self.assertEqual(result['good95'],2)
        self.assertEqual(result['tails']['0.5'],1)
        self.assertEqual(result['lineage_counts'],dict(final_from_good_seed=2,final_from_other_seed=1,no_matched_final=1))

    def test_collinear_seed_is_classified_without_infinite_radius(self):
        fit=taubin([[0.,1.],[1000.,1.],[2000.,1.]])
        self.assertFalse(fit['valid'])

    def test_known_circle_large_radius_and_sign(self):
        z=np.array([20.,1250.,2440.]);radius=1000/.000165
        for charge in [-1,1]:
            cy=-charge*radius;y=cy+charge*np.sqrt(radius**2-z**2)
            fit=taubin(np.stack([z,y],axis=1))
            self.assertLess(abs(fit['radius_mm']/radius-1),1e-7)
            self.assertEqual(1 if fit['center'][1]<0 else -1,charge)

    def test_compact_trace_detects_shared_smoother_overwrite(self):
        a=state(0,.001,.001,.002);b=state(1,.001,.002,.002)
        trace=trace_record([b,a]);self.assertAlmostEqual(trace['steps'][0]['smoothing_identity_error'],0.)
        a['smoothed'][4]=.003
        trace=trace_record([b,a]);self.assertAlmostEqual(trace['steps'][0]['smoothing_identity_error'],.001)
        refit=dict(a,stage='refit_state')
        reduced=list(reduce_event([a,refit]))
        self.assertEqual(len({(r['stage'],r['seed'],r['candidate'],r['state']) for r in reduced}),2)

    def test_compact_trace_closes_gain_times_innovation(self):
        row=state(0,.001,.011,.011);cov=np.eye(6);cov[0,4]=cov[4,0]=.2
        row.update(measurement_flag=True,predicted_cov=cov.tolist(),measurement=[.1],
                   projector=[[1,0,0,0,0,0]],measurement_cov=[[1.]])
        step=trace_record([row])['steps'][0]
        self.assertAlmostEqual(step['update_gain_qop'][0],.1)
        self.assertAlmostEqual(step['update_identity_error'],0.,places=14)


if __name__=='__main__':unittest.main()
