#!/usr/bin/env python3
"""Extract historical Calypso Ntuple branches without charge-dependent matching."""
import argparse
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


TRACKS = dict(chi2='Track_Chi2',ndof='Track_nDoF',n_layers='Track_nLayers',
    rmax_mm='Track_r_atMaxRadius',charge='Track_charge',momentum_GeV='Track_p0',
    barcode='t_barcode',pdg='t_pdg',purity_ndof='t_truthHitRatio',
    truth_fiducial='isFiducial',truth_momentum_GeV='t_p',
    x_F_mm='Track_x0',y_F_mm='Track_y0',z_F_mm='Track_z0',
    px_F_GeV='Track_px0',py_F_GeV='Track_py0',pz_F_GeV='Track_pz0')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args()
    task=json.loads(args.task.read_text())
    output=Path(task['output']);output.mkdir(parents=True,exist_ok=True)
    events,tracks=[],[]
    source_entry=0
    with uproot.open(task['input']) as file:
        tree=file['nt']
        branches=['run','eventID','longTracks','truthParticleBarcode','truthParticleIsFiducial',
                  *TRACKS.values()]
        for data in tree.iterate(branches,step_size=20000,library='ak'):
            size=len(data)
            entries=np.arange(source_entry,source_entry+size,dtype=np.int64)
            primary=(data['truthParticleBarcode']>0)&(data['truthParticleBarcode']<200000)
            assert ak.all(ak.sum(primary,axis=1)<=1)
            event=dict(source_entry=entries,run=ak.to_numpy(data['run']),event=ak.to_numpy(data['eventID']),
                primary_barcode=ak.to_numpy(ak.fill_none(ak.firsts(data['truthParticleBarcode'][primary]),0)),
                primary_fiducial_valid=ak.to_numpy(ak.sum(primary,axis=1)==1),
                primary_fiducial=ak.to_numpy(ak.any(data['truthParticleIsFiducial']&primary,axis=1)),
                long_tracks=ak.to_numpy(data['longTracks']))
            lengths=ak.to_numpy(ak.num(data['Track_charge']))
            assert np.array_equal(lengths,event['long_tracks'])
            row={name:ak.to_numpy(ak.flatten(data[branch])) for name,branch in TRACKS.items()}
            for name in list(row):
                if name.endswith('_GeV'): row[name]=row[name]*.001
            row['source_entry']=np.repeat(entries,lengths)
            row['ordinal']=ak.to_numpy(ak.flatten(ak.local_index(data['Track_charge'])))
            row['run']=np.repeat(event['run'],lengths);row['event']=np.repeat(event['event'],lengths)
            events.append(event);tracks.append(row);source_entry+=size
    assert source_entry==task['events']
    for name,rows in [('events',events),('tracks',tracks)]:
        np.savez_compressed(output/(name+'.npz'),**{k:np.concatenate([r[k] for r in rows]) for k in rows[0]})
    (output/'source_task.json').write_text(json.dumps(task,indent=2)+'\n')
    (output/'summary.json').write_text(json.dumps(dict(source_events=source_entry,
        tracks=sum(len(r['source_entry']) for r in tracks),input_tree='nt',
        truth_purity_denominator='Historical t_truthHitRatio uses ndof+5, not measured cluster count',
        truth_coverage='truthParticleBarcode is subject to the historical MinMomentum gate; absent primary keys are unknown fiducial status, not nonfiducial truth',
        matching='Stored SDO majority barcode and PDG; no charge or reconstructed momentum selection',
        good_track_reproduction='25/7/12 thresholds with r<100 mm and r<95 mm kept separately; original Alma analysis thresholds not yet recovered'),indent=2)+'\n')


if __name__=='__main__':
    main()
