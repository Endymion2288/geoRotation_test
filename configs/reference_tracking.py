#!/usr/bin/env python3
"""Run the native Calypso CKF and its original ROOT states writer."""
import argparse
import json
import os
import socket
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--events', type=int, required=True)
    parser.add_argument('--skip', type=int, default=0)
    parser.add_argument('--material', required=True)
    parser.add_argument('--variant', required=True)
    parser.add_argument('--mask', default='0,1,2,3,4,5')
    parser.add_argument('--isolate-smoothing', action='store_true')
    parser.add_argument('--native-final-purity', action='store_true',
        help='Observe final-track SDO purity without reading alternative truth references')
    parser.add_argument('--audit-states', action='store_true',
        help='Diagnostic only: record pre/post-isolation states with an audit-enabled library')
    parser.add_argument('--initial-variance-inflation', type=float, nargs=5,
        help='Override only the five seed variance inflation factors')
    args = parser.parse_args()
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
    from reference_truth_observer import ReferenceTruthObserver
    Configurable.configurableRun3Behavior = True
    flags = initConfigFlags()
    flags.Input.isMC = True
    flags.Input.ProjectName = 'mc24'
    flags.Input.Files = [args.input]
    flags.IOVDb.DatabaseInstance = 'OFLP200'
    flags.GeoModel.FaserVersion = 'FASERNU-04'
    flags.IOVDb.GlobalTag = 'OFLCOND-FASER-05'
    flags.TrackingGeometry.MaterialSource = args.material
    flags.Exec.SkipEvents = args.skip
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
    masked = [int(x) for x in args.mask.split(',') if x]
    cfg = CKF2Cfg(flags, actsOutputTag='reference_3station_forward', name='CKF_Reference',
        maskedLayers=masked, OutputCollection='CKFTrackCollectionWithoutIFT',
        BackwardPropagation=False, noDiagnostics=False, seedWriter=False,
        StatesWriter=True, SummaryWriter=False)
    ckf = cfg.getEventAlgo('CKF_Reference')
    if args.isolate_smoothing:
        ckf.IsolateCandidateSmoothing = True
    if args.audit_states:
        ckf.AuditFile = 'stages.jsonl'
        ckf.AuditStates = True
        ckf.AuditSelectedStatesOnly = True
    if args.initial_variance_inflation is not None:
        sigmas = [.9, .019, .033, .00029, .0001]
        for prop, factor, sigma in zip(
                ['covLoc0', 'covLoc1', 'covPhi', 'covTheta', 'covQOverP'],
                args.initial_variance_inflation, sigmas):
            setattr(ckf.TrackSeed, prop, factor*sigma*sigma)
    writer = ckf.RootTrajectoryStatesWriterTool.clone('ReferenceCKFStates')
    writer.FilePath = 'track_states.root'
    writer.noDiagnostics = False
    ckf.RootTrajectoryStatesWriterTool = writer
    fitter = ckf.KalmanFitterTool1
    fitter.noDiagnostics = True
    fitter.StatesWriter = False
    fitter.SummaryWriter = False
    acc.merge(cfg)
    acc.addEventAlgo(ReferenceTruthObserver('ReferenceTruthObserver', skip=args.skip,
        rotated=args.variant != 'beforeGeo', native_final_purity=args.native_final_purity))
    acc.addService(CompFactory.AthenaEventLoopMgr(EventPrintoutInterval=1000))
    acc.getService('MessageSvc').defaultLimit = 20
    covariance = {p: getattr(ckf.TrackSeed, p) for p in ['covLoc0','covLoc1','covPhi','covTheta','covQOverP']}
    manifest = dict(args=vars(args), calypso_dir=os.environ.get('Calypso_DIR'),
        library=os.environ.get('GEOROTATION_LIBRARY_DIR'), python_path=os.environ.get('PYTHONPATH'),
        geometry='FASERNU-04', conditions='OFLCOND-FASER-05', seed_covariance=covariance,
        audit=args.audit_states, reset_event_maps=False, isolate_smoothing=args.isolate_smoothing,
        hostname=socket.gethostname(), athena=os.environ.get('AtlasVersion'))
    Path('run.json').write_text(json.dumps(manifest, indent=2)+'\n')
    begin = time.monotonic()
    status = acc.run(maxEvents=args.events)
    loaded = sorted({line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines()
        if '/libFaserActs' in line and '.so' in line})
    manifest.update(success=status.isSuccess(), wall_s=time.monotonic()-begin,
        loaded_acts_libraries=loaded)
    Path('run.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0 if status.isSuccess() else 1


if __name__ == '__main__':
    raise SystemExit(main())
