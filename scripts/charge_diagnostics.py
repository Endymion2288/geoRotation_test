"""Reduce read-only native Calypso charge lineage on the worker."""
from pathlib import Path
import shutil

import numpy as np


STAGES = {0:'seed',1:'all_candidate_reference',2:'selected_predicted',
          3:'selected_filtered',4:'selected_isolated_smoothed',5:'selected_reference',
          6:'converted',7:'final',8:'all_candidate_measurement',10:'rejected_candidate'}
DTYPE = [(k,t) for k,t in [
    ('local_event_index','i8'),('run','u8'),('event','u8'),('stage','i1'),('seed','i4'),
    ('candidate','i4'),('selected','i4'),('ordinal','i4'),('state_order','i4'),
    ('geometry_id','u8'),('hit_id','u8'),('qop_GeVinv','f8'),('c44_GeVinv2','f8'),
    ('chi2','f8'),('n_measurements','i4'),('status','i4'),('phi','f8'),('theta','f8'),
    ('loc0','f8'),('loc1','f8'),('truth_barcodes','U256')]]


def extract(output, task, arrays):
    records = np.atleast_1d(np.loadtxt('charge_lineage.csv',delimiter=',',skiprows=1,
        dtype=DTYPE,encoding='utf-8'))
    fields = {name:records[name] for name in records.dtype.names}
    fields['source_entry'] = fields['local_event_index'] + task['skip']
    assert ((fields['source_entry']>=task['skip']) &
            (fields['source_entry']<task['skip']+task['events'])).all()
    assert np.isin(fields['stage'],list(STAGES)).all()
    final = fields['stage']==7
    # Invalid final parameter collections are retained as status=-1, and have
    # no row in the pre-existing slope observer. Check the valid rows only.
    valid_final = final & (fields['status']>=0)
    for name in ['source_entry','run','event','ordinal']:
        np.testing.assert_array_equal(fields[name][valid_final],arrays['tracks'][name])
    np.testing.assert_array_equal(fields['qop_GeVinv'][valid_final],arrays['tracks']['qop_GeVinv'])
    np.testing.assert_array_equal(fields['c44_GeVinv2'][valid_final],arrays['tracks']['qop_variance_GeVinv2'])
    np.savez_compressed('charge_lineage.npz',**fields)
    shutil.copy2('charge_lineage.npz',Path(output)/'charge_lineage.npz')
    return dict(rows=len(records),stages={STAGES[int(k)]:int((fields['stage']==k).sum())
        for k in np.unique(fields['stage'])},final_identity_qop_covariance_exact=True,
        invalid_final_parameters=int((final & (fields['status']<0)).sum()),
        refit_fallback=int((final & (fields['status']==0)).sum()),
        refit_success=int((final & (fields['status']==1)).sum()),
        seed_identity='Actual iseed recorded independently; native trackGroup is unchanged',
        state_order='Reverse traversal, last measurement first; filtering order is descending state_order',
        truth_matching='Unique SDO barcode contributions per measured cluster; no reco sign/momentum cut')


def refit_counters(output, task):
    """Keep only actual final refit outcomes; input CSV lives on scratch."""
    result = {}
    for path in sorted(Path.cwd().glob('mode_refit_*.csv')):
        rows = np.atleast_1d(np.loadtxt(path, delimiter=',', skiprows=1,
                                      dtype=DTYPE, encoding='utf-8'))
        rows = rows[rows['stage'] == 7]
        fields = {k: rows[k] for k in ['local_event_index', 'run', 'event', 'candidate',
                                      'ordinal', 'qop_GeVinv', 'status']}
        fields['source_entry'] = fields['local_event_index']+task['skip']
        fields['sample_entry'] = fields['local_event_index']+task['sample_entry']
        np.savez_compressed(Path(output)/(path.stem+'.npz'), **fields)
        result[path.stem] = dict(final=len(rows), success=int((rows['status'] == 1).sum()),
            fallback=int((rows['status'] == 0).sum()), invalid=int((rows['status'] < 0).sum()))
    assert len(result) == (2 if '--merge-directions' in task.get('options', []) else 1)
    return result
