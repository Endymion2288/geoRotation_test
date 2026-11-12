#!/usr/bin/env python3
"""Check the native Ntuple, observer invariance, and covariance mean closure."""
import argparse
import json
from pathlib import Path

import numpy as np


def load(path):
    with np.load(path, allow_pickle=False) as file:
        return {key:file[key] for key in file.files}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    result = dict(source_events_per_check=20,variants={})
    for variant in ['beforeGeo','afterGeo']:
        folders = {mode:args.root/('smoke_'+mode+'_v1')/'shards'/variant/'s0000'
                   for mode in ['identity','covariance','native']}
        for folder in folders.values():
            assert (folder/'exit_code.txt').read_text().strip() == '0', folder
        native = {mode:load(folder/'native_ntuple_branches.npz') for mode,folder in folders.items()}
        for mode in ['identity','covariance']:
            assert native[mode].keys()==native['native'].keys()
            for key in native[mode]:
                np.testing.assert_array_equal(native[mode][key],native['native'][key],err_msg=variant+' '+mode+' '+key)
        identity,covariance = (load(folders[mode]/'tracks.npz') for mode in ['identity','covariance'])
        for key in ['source_entry','ordinal','barcode','charge','station_bits','hit_identifiers',
                    'upstream_raw','upstream_F','projection_valid']:
            np.testing.assert_array_equal(identity[key],covariance[key])
        columns = [0,1,2,3,4,5,7]
        a,b = (x['projection'][:,:,columns] for x in [identity,covariance])
        np.testing.assert_allclose(a,b,rtol=2e-11,atol=2e-10,equal_nan=True)
        diff=np.abs(a-b)
        validations = {mode:json.loads((folders[mode]/'ntuple_validation.json').read_text()) for mode in ['identity','covariance']}
        variance=covariance['projection'][:,:,6]
        valid=covariance['projection_valid']==1
        result['variants'][variant] = dict(native_Ntuple_all_saved_branches_exact=True,
            tracks_hits_charge_and_covariance_only_means_equal=True,
            maximum_projection_mean_difference=float(np.nanmax(diff,initial=0)),
            propagated_variance_negative=int(np.count_nonzero(valid & (variance<0))),
            propagated_variance_nonfinite=int(np.count_nonzero(valid & ~np.isfinite(variance))),
            ntuple_validation=validations)
    (args.root/'smoke_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
