#!/usr/bin/env python3
"""High-throughput Calypso reconstruction runner for latest seed fixes, particle-hypothesis check, and seeding checks."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import socket
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--events', type=int, required=True)
    parser.add_argument('--skip', type=int, default=0)
    parser.add_argument('--material', required=True)
    parser.add_argument('--mask', default='0,1,2,3,4,5')
    parser.add_argument('--audit-prefix', default='fast_audit')
    parser.add_argument('--write-states-root', action='store_true',
                        help='Write track_states.root when running ckf_residual_perLayer.C')

    # Step-by-step ablation & worktree fix flags (defaults match fixed worktree C5)
    parser.add_argument('--legacy-swapped-cov-loc', action='store_true')
    parser.add_argument('--initial-variance-inflation', type=float, nargs=5,
                        default=[100., 100., 100., 100., 1000.])
    parser.add_argument('--seed-b-tesla', type=float, default=0.473662)
    parser.add_argument('--fix-gety-nan', type=int, choices=[0, 1], default=1)
    parser.add_argument('--use-local-tangent', type=int, choices=[0, 1], default=1)

    # Exploratory Check 1: Particle hypothesis (13=muon, 11=electron, 211=pion) with identical material effects
    parser.add_argument('--abs-pdg-hypothesis', type=int, choices=[11, 13, 211], default=13)
    parser.add_argument('--ckf-scattering', type=int, choices=[0, 1], default=1)
    parser.add_argument('--ckf-energy-loss', type=int, choices=[0, 1], default=1)

    # Exploratory Check 2: Seed combinatorics prefilter, deduplication, and 18-cluster field-profile strip fit
    parser.add_argument('--seed-prefilter', type=int, choices=[0, 1], default=0)
    parser.add_argument('--seed-improved-dedup', type=int, choices=[0, 1], default=0)
    parser.add_argument('--seed-strip-field-fit', type=int, choices=[0, 1], default=0)
    parser.add_argument('--seed-max-reduced-chi2', type=float, default=35.0)
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

    masked = [int(x) for x in args.mask.split(',') if x]
    cfg = CKF2Cfg(
        flags,
        actsOutputTag='validation_3station_forward',
        name='CKF_Validation',
        maskedLayers=masked,
        OutputCollection='CKFTrackCollectionWithoutIFT',
        BackwardPropagation=False,
        noDiagnostics=not args.write_states_root,
        seedWriter=False,
        StatesWriter=bool(args.write_states_root),
        SummaryWriter=False,
        legacySwappedCovLoc=bool(args.legacy_swapped_cov_loc),
        initialVarianceInflation=list(args.initial_variance_inflation),
        absPdgHypothesis=int(args.abs_pdg_hypothesis),
    )
    ckf = cfg.getEventAlgo('CKF_Validation')
    ckf.IsolateCandidateSmoothing = True
    ckf.CkfMultipleScattering = bool(args.ckf_scattering)
    ckf.CkfEnergyLoss = bool(args.ckf_energy_loss)
    if args.audit_prefix:
        ckf.FastBinaryAuditPrefix = args.audit_prefix

    ts = ckf.TrackSeed
    ts.SeedBFieldTesla = float(args.seed_b_tesla)
    ts.FixGetYNaN = bool(args.fix_gety_nan)
    ts.UseLocalTangent = bool(args.use_local_tangent)
    ts.AbsPdgHypothesis = int(args.abs_pdg_hypothesis)
    ts.SeedPreFilter = bool(args.seed_prefilter)
    ts.SeedImprovedDedup = bool(args.seed_improved_dedup)
    ts.SeedStripFieldFit = bool(args.seed_strip_field_fit)
    ts.SeedMaxReducedChi2 = float(args.seed_max_reduced_chi2)

    if args.write_states_root:
        writer = ckf.RootTrajectoryStatesWriterTool.clone('ValidationCKFStates')
        writer.FilePath = 'track_states.root'
        writer.noDiagnostics = False
        ckf.RootTrajectoryStatesWriterTool = writer

    fitter = ckf.KalmanFitterTool1
    fitter.noDiagnostics = True
    fitter.StatesWriter = False
    fitter.SummaryWriter = False

    acc.merge(cfg)
    acc.addService(CompFactory.AthenaEventLoopMgr(EventPrintoutInterval=1000))
    acc.getService('MessageSvc').defaultLimit = 20

    covariance = {p: getattr(ts, p) for p in ['covLoc0', 'covLoc1', 'covPhi', 'covTheta', 'covQOverP']}
    manifest = dict(
        args=vars(args),
        calypso_dir=os.environ.get('Calypso_DIR'),
        library=os.environ.get('GEOROTATION_LIBRARY_DIR'),
        python_path=os.environ.get('PYTHONPATH'),
        geometry='FASERNU-04',
        conditions='OFLCOND-FASER-05',
        seed_covariance=covariance,
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
