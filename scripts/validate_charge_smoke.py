#!/usr/bin/env python3
"""Check read-only diagnostics against unchanged actual Calypso reconstruction."""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args();result={}
    for variant in ['beforeGeo','afterGeo']:
        for dsid in ['100044','100116','100117']:
            base=args.root/'charge_smoke'/variant/dsid
            checked=[]
            for mode in ['on','off']:
                folder=base/mode
                assert (folder/'exit_code.txt').read_text().strip()=='0'
                run=json.loads((folder/'run.json').read_text())
                expected={str((args.root/('software_'+variant)/'lib'/name).resolve()) for name in [
                    'libFaserActsGeometry.so','libFaserActsGeometryLib.so','libFaserActsKalmanFilter.so']}
                expected.add(str((args.root/'observer/lib/libCalypsoSlopeDiagnostics.so').resolve()))
                assert set(run['loaded_acts_libraries'])==expected
                assert run['isolate_smoothing'] and run['charge_diagnostics']==(mode=='on')
            for name in ['events','truth','tracks','segments']:
                with np.load(base/'on'/(name+'.npz')) as a,np.load(base/'off'/(name+'.npz')) as b:
                    assert a.files==b.files
                    for field in a.files:
                        np.testing.assert_array_equal(a[field],b[field],err_msg=variant+'/'+dsid+'/'+name+'/'+field)
                        checked.append(name+'/'+field)
                    if variant=='afterGeo':
                        with np.load(args.root/'initial_smoke'/dsid/'default'/(name+'.npz')) as clean:
                            for field in clean.files:np.testing.assert_array_equal(a[field],clean[field])
            with np.load(base/'on/charge_lineage.npz') as a:
                seed_keys=set(zip(a['local_event_index'][a['stage']==0],a['seed'][a['stage']==0]))
                candidates=a['stage']==1
                assert all(k in seed_keys for k in zip(a['local_event_index'][candidates],a['seed'][candidates]))
                ref=np.flatnonzero(a['stage']==5);converted=np.flatnonzero(a['stage']==6)
                for field in ['local_event_index','selected','candidate','seed']:
                    np.testing.assert_array_equal(a[field][ref],a[field][converted])
                valid=np.isfinite(a['qop_GeVinv'][ref])&(a['qop_GeVinv'][ref]!=0)
                np.testing.assert_allclose(a['qop_GeVinv'][ref][valid],a['qop_GeVinv'][converted][valid],rtol=1e-13,atol=0)
                np.testing.assert_array_equal(np.sign(a['qop_GeVinv'][ref][valid]),np.sign(a['qop_GeVinv'][converted][valid]))
                with np.load(base/'on/tracks.npz') as tracks:
                    q=tracks['qop_GeVinv'];finite=np.isfinite(q)&(q!=0)
                    np.testing.assert_array_equal(np.sign(q[finite]),tracks['charge'][finite])
                    np.testing.assert_allclose(1/np.abs(q[finite]),np.linalg.norm(tracks['upstream_F'][finite,3:],axis=1),rtol=1e-13)
                result[variant+'|'+dsid]=dict(equal_diagnostic_on_off_arrays=len(checked),
                    actual_libraries_match=True,all_candidate_seed_indices_resolved=True,
                    conversion_signed_qop_preserved=int(valid.sum()),
                    extra_observer_qop_unit_sign_checks=True,
                    clean_tracking_agreement=variant=='afterGeo')
    (args.root/'charge_smoke/validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
