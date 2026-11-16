"""Regression checks for attribution and source identity, not a tracking model."""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
from latest_momentum import STATE_FIELDS, covariance_health, seed_good
from latest_momentum_study import transitions, compare_paths, native_formula_residual, model_gate_validation


class LatestMomentumTests(unittest.TestCase):
    def data(self):
        names=['file',*STATE_FIELDS]
        f={name:i for i,name in enumerate(names)}
        rows=[]
        # Numerical accounting fixtures only: two independent files share entry.
        for file in [0,1]:
            for station,measurement,predicted,filtered,previous in [
                    (1,1,.002,.003,.001),(1,0,.0031,.0031,.003),
                    (2,1,.0032,.004,.0031)]:
                row=np.zeros(len(names))
                values=dict(file=file,entry=4,candidate=0,seed=2,station=station,layer=0,
                    measurement=measurement,github=1,purity=1,good_seed=1,same_seed_truth=1,
                    p_production_GeV=1000,truth_charge=1,track_barcode=7,
                    seed_qop=.001,q_predicted=predicted,q_filtered=filtered,q_smoothed=.001,
                    dq_transport=predicted-previous,dq_update=filtered-predicted)
                for k,v in values.items():row[f[k]]=v
                rows.append(row)
        return dict(values=np.asarray(rows),fields=names,hits=np.asarray([1,0,2,1,0,2]),
                    coverage=[dict(file_index=i,events=10) for i in [0,1]])

    def test_transport_update_identity_with_passive_state(self):
        rows,fields,paths=transitions(self.data())
        c={name:i for i,name in enumerate(fields)}
        self.assertEqual(len(rows),2)
        np.testing.assert_allclose(rows[:,c['sum_transport']],.0012,atol=1e-15)
        np.testing.assert_allclose(rows[:,c['sum_updates']],.001,atol=1e-15)
        self.assertEqual(paths[0],paths[1])

    def test_file_identity_is_part_of_matching(self):
        a=self.data();b=self.data()
        b['values']=b['values'][3:];b['hits']=b['hits'][3:]
        result=compare_paths(a,b)
        self.assertEqual(result['same_seed_hits_truth'],1)

    def test_second_sensor_prediction_is_distinct_from_station_entry(self):
        data=self.data();f={name:i for i,name in enumerate(data['fields'])}
        rows=[];hits=[]
        for begin in [0,3]:
            rows.extend(data['values'][begin:begin+3]);hits.extend(data['hits'][begin:begin+3])
            second=data['values'][begin+2].copy()
            second[f['q_predicted']]=.004;second[f['q_filtered']]=.005
            second[f['dq_transport']]=0;second[f['dq_update']]=.001
            rows.append(second);hits.append(3)
        data.update(values=np.asarray(rows),hits=np.asarray(hits))
        tracks,columns,_=transitions(data);c={name:i for i,name in enumerate(columns)}
        np.testing.assert_allclose(tracks[:,c['s2_predicted']],.0032)
        np.testing.assert_allclose(tracks[:,c['s2_first_filtered']],.004)
        np.testing.assert_allclose(tracks[:,c['s2_second_predicted']],.004)

    def test_covariance_diagnostic_detects_indefinite_input(self):
        c=np.eye(6);c[0,4]=c[4,0]=1.2
        self.assertLess(covariance_health(c),0)
        self.assertEqual(covariance_health(np.eye(6)),1)

    def test_original_good_seed_gate_does_not_cut_phi(self):
        values=[0]*14
        values[4]=.2;values[10]=0
        values[5]=values[11]=np.pi/2
        values[6]=values[12]=.001
        self.assertTrue(seed_good(values))

    def test_reference_zero_division_is_not_algebraically_cancelled(self):
        values=native_formula_residual(np.asarray([.001,0.,.001]),np.asarray([0.,.001,np.nan]))
        np.testing.assert_equal(values,np.asarray([0.,-1.,np.nan]))

    def test_seed_model_gate_evaluation_excludes_calibration_file(self):
        data=self.data()
        seeds=[]
        for file in [0,1]:
            seed=np.zeros(16)
            seed[:3]=[file,4,2]
            seed[8]=.001/1.16;seed[14]=.001
            seeds.append(seed)
        data['seeds']=np.asarray(seeds)
        closure=dict(residual=dict(quantiles=dict(p50=.16)),source_events=100)
        result,_,_,masks=model_gate_validation(data,closure)
        self.assertEqual(result['evaluation_source_events'],10)
        self.assertEqual(result['selections']['original']['tracks'],0)
        self.assertEqual(result['selections']['native_model_calibrated']['tracks'],1)
        self.assertEqual(masks['native_model_calibrated'].tolist(),[False,True])


if __name__=='__main__':
    unittest.main()
