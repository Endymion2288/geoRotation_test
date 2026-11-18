#!/usr/bin/env python3
"""Run latest Calypso CKF with opt-in seed forensics and single-factor covariance/seed-B controls."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import socket
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--events', type=int, required=True)
    parser.add_argument('--skip', type=int, default=0)
    parser.add_argument('--material', required=True)
    parser.add_argument('--variant', default='afterGeo')
    parser.add_argument('--mask', default='0,1,2,3,4,5')
    smoothing = parser.add_mutually_exclusive_group()
    smoothing.add_argument('--isolate-smoothing', action='store_true')
    smoothing.add_argument('--shared-smoothing', action='store_true')
    parser.add_argument('--native-final-purity', action='store_true')
    parser.add_argument('--seed-unbiased', action='store_true')
    parser.add_argument('--initial-variance-inflation', type=float, nargs=5)
    parser.add_argument('--swap-cov-loc', action='store_true',
                        help='Assign fine y_F sigma (0.019 mm) to covLoc0 and stereo x_F sigma (0.9 mm) to covLoc1')
    parser.add_argument('--seed-b-tesla', type=float,
                        help='Override CircleFitTrackSeedTool SeedBFieldTesla (default 0.55 T)')
    parser.add_argument('--momentum-audit', action='store_true')
    parser.add_argument('--ckf-scattering', type=int, choices=[0, 1], default=1)
    parser.add_argument('--ckf-energy-loss', type=int, choices=[0, 1], default=1)
    parser.add_argument('--ckf-max-step-mm', type=float, default=0.)
    parser.add_argument('--seed-override-file')
    parser.add_argument('--seed-override-mode', choices=['qop', 'full'], default='qop')
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
    if importlib.util.find_spec('TrackerSegmentFit.GhostBustersConfig'):
        from TrackerSegmentFit.GhostBustersConfig import GhostBustersCfg
    else:
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
    if args.seed_unbiased:
        ckf.seedWriter = True
    if args.seed_b_tesla is not None:
        ckf.TrackSeed.SeedBFieldTesla = float(args.seed_b_tesla)
    if args.momentum_audit:
        ckf.MomentumAuditFile = 'momentum_audit.jsonl'
        ckf.MomentumScattering = bool(args.ckf_scattering)
        ckf.MomentumEnergyLoss = bool(args.ckf_energy_loss)
        ckf.MomentumMaxStepSize = args.ckf_max_step_mm
        if args.seed_override_file:
            ckf.MomentumSeedOverride = args.seed_override_file
            ckf.MomentumSeedOverrideMode = args.seed_override_mode
    if args.isolate_smoothing:
        ckf.IsolateCandidateSmoothing = True
    elif args.shared_smoothing:
        ckf.IsolateCandidateSmoothing = False

    if args.initial_variance_inflation is not None or args.swap_cov_loc:
        factors = args.initial_variance_inflation or [100., 100., 100., 10000., 100000.]
        sigmas = [.019, .9, .033, .00029, .0001] if args.swap_cov_loc else [.9, .019, .033, .00029, .0001]
        for prop, factor, sigma in zip(
                ['covLoc0', 'covLoc1', 'covPhi', 'covTheta', 'covQOverP'],
                factors, sigmas):
            setattr(ckf.TrackSeed, prop, factor * sigma * sigma)

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
                                            rotated=args.variant != 'beforeGeo',
                                            native_final_purity=args.native_final_purity))
    acc.addService(CompFactory.AthenaEventLoopMgr(EventPrintoutInterval=1000))
    acc.getService('MessageSvc').defaultLimit = 20

    covariance = {p: getattr(ckf.TrackSeed, p) for p in ['covLoc0', 'covLoc1', 'covPhi', 'covTheta', 'covQOverP']}
    manifest = dict(
        args=vars(args),
        calypso_dir=os.environ.get('Calypso_DIR'),
        library=os.environ.get('GEOROTATION_LIBRARY_DIR'),
        python_path=os.environ.get('PYTHONPATH'),
        geometry='FASERNU-04',
        conditions='OFLCOND-FASER-05',
        seed_covariance=covariance,
        swap_cov_loc=bool(args.swap_cov_loc),
        seed_b_tesla=float(args.seed_b_tesla) if args.seed_b_tesla is not None else 0.55,
        isolate_smoothing=bool(getattr(ckf, 'IsolateCandidateSmoothing', False)),
        hostname=socket.gethostname(),
        athena=os.environ.get('AtlasVersion'),
    )
    Path('run.json').write_text(json.dumps(manifest, indent=2) + '\n')
    begin = time.monotonic()
    status = acc.run(maxEvents=args.events)
    loaded = sorted({
        line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines()
        if any('/lib' + name in line for name in
               ['FaserActs', 'RootWriterTools', 'TrackingTruthUtils', 'TrackerSegmentFit']) and '.so' in line
    })
    manifest.update(success=status.isSuccess(), wall_s=time.monotonic() - begin,
                    loaded_acts_libraries=loaded)
    Path('run.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return 0 if status.isSuccess() else 1


if __name__ == '__main__':
    raise SystemExit(main())
