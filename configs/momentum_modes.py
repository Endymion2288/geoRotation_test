#!/usr/bin/env python3
"""Run the four production CKF configurations and both real merge algorithms."""
import argparse
import gzip
import json
import math
import os
from itertools import zip_longest
from pathlib import Path
import sys
import time

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'scripts'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tasks', type=Path, required=True)
    parser.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    task = json.loads(args.tasks.read_text())[args.index]
    destination = Path(task['output'])
    started = time.monotonic()

    from AthenaCommon.Configurable import Configurable
    from AthenaConfiguration.ComponentFactory import CompFactory
    from CalypsoConfiguration.AllConfigFlags import initConfigFlags
    from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    from FaserGeoModel.FaserGeoModelConfig import FaserGeometryCfg
    from TrackerPrepRawDataFormation.TrackerPrepRawDataFormationConfig import FaserSCT_ClusterizationCfg
    from TrackerSpacePointFormation.TrackerSpacePointFormationConfig import TrackerSpacePointFinderCfg
    from TrackerSegmentFit.TrackerSegmentFitConfig import SegmentFitAlgCfg
    from FaserActsKalmanFilter.GhostBustersConfig import GhostBustersCfg
    from FaserActsKalmanFilter.CKF2Config import CKF2Cfg
    from Reconstruction.xAODTruthCnvAlgConfig import xAODTruthCnvAlgCfg
    import AthenaPython.PyAthena as PyAthena

    Configurable.configurableRun3Behavior = True
    flags = initConfigFlags()
    flags.Input.isMC = True
    flags.Input.ProjectName = 'mc24'
    flags.Input.Files = [task['input']]
    flags.IOVDb.DatabaseInstance = 'OFLP200'
    flags.GeoModel.FaserVersion = 'FASERNU-04'
    flags.IOVDb.GlobalTag = 'OFLCOND-FASER-05'
    flags.TrackingGeometry.MaterialSource = task['material']
    flags.Exec.SkipEvents = task['skip']
    flags.Concurrency.NumThreads = 0
    flags.Concurrency.NumConcurrentEvents = 1
    flags.Output.doWriteESD = False
    flags.lock()
    acc = MainServicesCfg(flags)
    acc.merge(PoolReadCfg(flags))
    acc.merge(FaserGeometryCfg(flags))
    acc.merge(FaserSCT_ClusterizationCfg(flags, DataObjectName='SCT_RDOs', checkBadChannels=True))
    acc.merge(TrackerSpacePointFinderCfg(flags))
    acc.merge(SegmentFitAlgCfg(flags, ReducedChi2Cut=25., SharedHitFraction=.61,
                               MinClustersPerFit=5, TanThetaXZCut=.083))
    acc.merge(GhostBustersCfg(flags))
    acc.merge(xAODTruthCnvAlgCfg(flags))
    modes = [('forward3', False, False, 'CKFTrackCollectionWithoutIFT'),
             ('backward3', False, True, 'CKFTrackCollectionBackwardWithoutIFT'),
             ('forward4', True, False, 'CKFTrackCollection'),
             ('backward4', True, True, 'CKFTrackCollectionBackward')]
    for mode, ift, backward, collection in modes:
        name = 'CKF_' + mode
        cfg = CKF2Cfg(flags, actsOutputTag=mode, name=name,
                      maskedLayers=[] if ift else [0, 1, 2, 3, 4, 5],
                      OutputCollection=collection, BackwardPropagation=backward,
                      AuditFile=mode + '.jsonl', AuditStates=False, noDiagnostics=True)
        cfg.getEventAlgo(name).TrackSeed.ResetEventMaps = True
        acc.merge(cfg)
    merges = {'merge3': ('CKFTrackCollectionWithoutIFT', 'CKFTrackCollectionBackwardWithoutIFT'),
              'merge4': ('CKFTrackCollection', 'CKFTrackCollectionBackward')}
    for mode, (forward, backward) in merges.items():
        acc.addEventAlgo(CompFactory.MergeTracksAlg('Merge_' + mode,
            TrackCollectionForward=forward, TrackCollectionBackward=backward,
            MergeTrackCollection=mode, acceptTrackSeenByOne=True))

    class MergedWriter(PyAthena.Alg):
        def initialize(self):
            self.sg = PyAthena.py_svc('StoreGateSvc')
            self.stream = open('merged.jsonl', 'w')
            return PyAthena.StatusCode.Success

        def execute(self):
            result = {}
            for mode in merges:
                rows = []
                for track in self.sg.retrieve('TrackCollection', mode):
                    p = list(track.trackParameters())[0]
                    hits = [int(m.identify().get_compact()) for m in track.measurementsOnTrack()]
                    momentum = p.momentum()
                    rows.append(dict(hits=sorted(hits), momentum_A_GeV=[float(momentum[i])*.001 for i in range(3)],
                                     charge=float(p.charge()), fitter=int(track.info().trackFitter())))
                result[mode] = rows
            self.stream.write(json.dumps(result, allow_nan=False) + '\n')
            return PyAthena.StatusCode.Success

        def finalize(self):
            self.stream.close()
            return PyAthena.StatusCode.Success

    acc.addEventAlgo(MergedWriter('MergedWriter'))
    acc.getService('MessageSvc').defaultLimit = 20
    status = acc.run(maxEvents=task['events'])
    if status.isFailure():
        return 1

    def events(path):
        records = []
        with path.open() as stream:
            for line in stream:
                row = json.loads(line)
                if row['stage'] == 'event' and records:
                    yield records
                    records = []
                records.append(row)
        if records:
            yield records

    writers = {}
    for mode in [m[0] for m in modes] + list(merges):
        output = destination / mode
        output.mkdir()
        writers[mode] = gzip.open(output / 'stages.jsonl.gz', 'wt')
        record = dict(task, shard=task['shard'], output=str(output), variant=mode)
        (output / 'task.json').write_text(json.dumps(record, indent=2) + '\n')
        manifest = dict(success=True, file_index=task['file_index'], input=task['input'],
                        wall_seconds=time.monotonic()-started, args=dict(events=task['events'],skip=task['skip']),
                        material_source=task['material'], library_override=os.environ.get('GEOROTATION_LIBRARY_DIR'))
        (output / 'run.json').write_text(json.dumps(manifest, indent=2) + '\n')
        (output / 'resource.txt').write_text('Shared six-mode process; see parent resource.txt\n')
    validation = dict(events=0, merged_tracks=0, unmatched_merge_copies=0, ambiguous_origin=0)
    streams = [events(Path(mode + '.jsonl')) for mode, *_ in modes]
    with open('merged.jsonl') as merged:
        for group in zip_longest(*streams, merged):
            assert all(item is not None for item in group), 'Mode event coverage differs'
            bundles, line = group[:-1], group[-1]
            validation['events'] += 1
            merged_event = json.loads(line)
            for (mode, *_), records in zip(modes, bundles):
                for row in records:
                    writers[mode].write(json.dumps(row, separators=(',', ':'),allow_nan=False)+'\n')
            for mode, indices in [('merge3', (0, 1)), ('merge4', (2, 3))]:
                rows = [dict(x, collection=mode) for x in bundles[indices[0]]
                        if x['stage'] in ['event', 'truth', 'truth_hit']]
                originals = [x for i in indices for x in bundles[i] if x['stage'] == 'final']
                for j, copy in enumerate(merged_event[mode]):
                    matches = [x for x in originals if sorted(x['hits']) == copy['hits']
                               and x['charge'] == copy['charge']
                               and all(math.isclose(a, b, rel_tol=2e-14, abs_tol=1e-12)
                                       for a, b in zip(x['momentum_A_GeV'], copy['momentum_A_GeV']))]
                    validation['merged_tracks'] += 1
                    if not matches:
                        validation['unmatched_merge_copies'] += 1
                        raise AssertionError(('Merge changed parameters or hit IDs', mode, copy))
                    validation['ambiguous_origin'] += len(matches) > 1
                    row = dict(matches[0], collection=mode, candidate=j,
                               origin_collections=sorted({m['collection'] for m in matches}),
                               merge_fitter=copy['fitter'])
                    rows.append(row)
                for row in rows:
                    writers[mode].write(json.dumps(row, separators=(',', ':'),allow_nan=False)+'\n')
    for stream in writers.values():
        stream.close()
    assert validation['events'] == task['events']
    (destination/'merge_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    for mode in writers:
        (destination/mode/'exit_code.txt').write_text('0\n')
    print(json.dumps(validation))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
