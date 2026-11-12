#!/usr/bin/env python3
"""Native Calypso reconstruction with a separate, opt-in slope observer."""
import argparse
import json
import os
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
    parser.add_argument('--ntuple-check', action='store_true')
    parser.add_argument('--actual-covariance', action='store_true')
    parser.add_argument('--no-seed-writer', action='store_true')
    parser.add_argument('--no-observer', action='store_true', help='Native Ntuple baseline for the read-only observer check')
    parser.add_argument('--stage-audit', action='store_true', help='Read-only detailed states; requires the recorded Calypso audit library')
    parser.add_argument('--native-states', action='store_true', help='Enable existing Calypso CKF/refit ROOT state writers for both native branches')
    parser.add_argument('--isolate-smoothing', action='store_true', help='Controlled candidate-specific smoothing; changes saved states')
    parser.add_argument('--physical-ntuple-angles', action='store_true', help='Opt-in physical VetoNu angles in the actual rotated NtupleDumper')
    parser.add_argument('--truth-cluster-position', action='store_true', help='Calypso diagnostic: replace only local measured strip coordinates by majority SiHit centroids')
    parser.add_argument('--refit-covariance-scale', type=float, help='Controlled change to native KalmanFitterTool1.SeedCovarianceScale; CKF remains native')
    args = parser.parse_args()
    from AthenaCommon.Configurable import Configurable
    from AthenaConfiguration.ComponentFactory import CompFactory
    from CalypsoConfiguration.AllConfigFlags import initConfigFlags
    from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    from FaserGeoModel.FaserGeoModelConfig import FaserGeometryCfg
    from FaserActsGeometry.ActsGeometryConfig import ActsTrackingGeometryToolCfg
    from TrackerPrepRawDataFormation.TrackerPrepRawDataFormationConfig import FaserSCT_ClusterizationCfg
    from TrackerSpacePointFormation.TrackerSpacePointFormationConfig import TrackerSpacePointFinderCfg
    from TrackerSegmentFit.TrackerSegmentFitConfig import SegmentFitAlgCfg
    from FaserActsKalmanFilter.GhostBustersConfig import GhostBustersCfg
    from FaserActsKalmanFilter.CKF2Config import CKF2Cfg
    from Reconstruction.xAODTruthCnvAlgConfig import xAODTruthCnvAlgCfg
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
    cluster=FaserSCT_ClusterizationCfg(flags, DataObjectName='SCT_RDOs', checkBadChannels=True)
    if args.truth_cluster_position:
        cluster.getEventAlgo('Tracker::FaserSCT_Clusterization').SCT_ClusteringTool.globalPosAlg.UseTruthClusterPosition=True
    acc.merge(cluster)
    acc.merge(TrackerSpacePointFinderCfg(flags))
    acc.merge(SegmentFitAlgCfg(flags, ReducedChi2Cut=25., SharedHitFraction=.61,
                             MinClustersPerFit=5, TanThetaXZCut=.083))
    acc.merge(GhostBustersCfg(flags))
    acc.merge(xAODTruthCnvAlgCfg(flags))
    cfg = CKF2Cfg(flags, actsOutputTag='slope_3station_forward', name='CKF_Slope',
        maskedLayers=[0, 1, 2, 3, 4, 5], OutputCollection='CKFTrackCollectionWithoutIFT',
        BackwardPropagation=False, noDiagnostics=False, seedWriter=not args.no_seed_writer,
        StatesWriter=args.native_states, SummaryWriter=False)
    ckf = cfg.getEventAlgo('CKF_Slope')
    if args.stage_audit:
        ckf.AuditFile = 'stages.jsonl'
        ckf.AuditStates = True
        ckf.AuditSelectedStatesOnly = True
    if args.isolate_smoothing:
        ckf.IsolateCandidateSmoothing = True
    ckf.KalmanFitterTool1.noDiagnostics = not args.native_states
    ckf.KalmanFitterTool1.StatesWriter = args.native_states
    ckf.KalmanFitterTool1.RootTrajectoryStatesWriterTool.noDiagnostics = not args.native_states
    ckf.KalmanFitterTool1.SummaryWriter = False
    if args.refit_covariance_scale is not None:
        ckf.KalmanFitterTool1.SeedCovarianceScale=args.refit_covariance_scale
    acc.merge(cfg)
    result, geometry = ActsTrackingGeometryToolCfg(flags)
    acc.merge(result)
    propagator = CompFactory.FaserActsExtrapolationTool('SlopeNativeExtrapolator',
        TrackingGeometryTool=geometry, MaxSteps=10000)
    if not args.no_observer:
        acc.addEventAlgo(CompFactory.CalypsoSlopeDiagnosticsAlg('CalypsoSlopeDiagnostics',
            ExtrapolationTool=propagator, Rotated=args.variant != 'beforeGeo', Skip=args.skip,
            ActualCovariance=args.actual_covariance))
    if args.ntuple_check:
        from ntuple_tracking_check import EmptyWaveforms
        from NtupleDumper.NtupleDumperConfig import NtupleDumperAlgCfg
        acc.addEventAlgo(EmptyWaveforms('SlopeEmptyWaveforms'))
        ntuple_options = {'RotatedPhysicalAngles':True} if args.physical_ntuple_angles else {}
        acc.merge(NtupleDumperAlgCfg(flags, 'native_ntuple.root', McEventCollection='TruthEvent',
            TrackCollection='CKFTrackCollectionWithoutIFT', TrackCollectionWithoutIFT='CKFTrackCollectionWithoutIFT',
            RunVertexing=False, DoBlinding=False, DoTrackFilter=False, StableOnly=False, DoCalib=False, **ntuple_options))
    acc.addService(CompFactory.AthenaEventLoopMgr(EventPrintoutInterval=1000))
    acc.getService('MessageSvc').defaultLimit = 20
    manifest = dict(args=vars(args), geometry='FASERNU-04', conditions='OFLCOND-FASER-05',
        calypso_dir=os.environ.get('Calypso_DIR'), library=os.environ.get('GEOROTATION_LIBRARY_DIR'),
        seed_covariance={p: getattr(ckf.TrackSeed, p) for p in ['covLoc0','covLoc1','covPhi','covTheta','covQOverP']},
        refit_covariance_scale=ckf.KalmanFitterTool1.SeedCovarianceScale,
        audit=args.stage_audit, native_states=args.native_states, reset_event_maps=False, isolate_smoothing=args.isolate_smoothing, observer=not args.no_observer,
        propagation_material_energy_loss=propagator.InteractionEloss,
        propagation_material_scattering=propagator.InteractionMultiScatering)
    Path('run.json').write_text(json.dumps(manifest, indent=2)+'\n')
    begin = time.monotonic()
    status = acc.run(maxEvents=args.events)
    manifest.update(success=status.isSuccess(), wall_s=time.monotonic()-begin)
    Path('run.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0 if status.isSuccess() else 1


if __name__ == '__main__':
    raise SystemExit(main())
