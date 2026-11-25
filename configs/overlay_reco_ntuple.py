#!/usr/bin/env python3
"""Single-step inline Athena overlay reconstruction + TwoTrackOverlayNTuplerAlg.

Reconstructs an existing OverlaidRDO POOL file using unmodified Calypso
components (`FaserSCT_Clusterization`, `TrackerSpacePointFinder`,
`SegmentFitAlg`, `GhostBusters`, 3-station forward `CKF2`) and feeds
StoreGate directly into `TwoTrackOverlayNTuplerAlg` to write `/HIST/nt`
without serializing intermediate xAOD files to disk.

Coordinate-frame note:
  - When running on rotated GeoRotation Calypso (`--geometry-rotation true`),
    `CKFTrackCollectionWithoutIFT` is produced in the rotated ACTS frame
    v_A = (z_F, y_F, -x_F), while `Orig_CKFTrackCollectionWithoutIFT` (from
    unrotated EOS OverlaidRDOs) and `SegmentFit`/`Segments` are in the physical
    FASER frame v_F = (x_F, y_F, z_F).
  - The raw output `nt` tree must therefore be normalized by
    `overlay_acts_to_faser` (`actsToFaser` on the 12 `Track_{x,y,z,px,py,pz}{0,1}`
    branches) before downstream analysis.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time


def parse_bool(value: str) -> bool:
    norm = value.strip().lower()
    if norm in ("1", "true", "yes", "on"):
        return True
    if norm in ("0", "false", "no", "off"):
        return False
    raise argparse.ArgumentTypeError(f"Expected boolean string, got {value!r}")


def resolve_geometry(runtype: str) -> str:
    mapping = {
        "TI12Data03": "FASERNU-03",
        "TI12MC03": "FASERNU-03",
        "TI12Data04": "FASERNU-04",
        "TI12MC04": "FASERNU-04",
        "TI12Data05": "FASERNU-05",
        "TI12MC05": "FASERNU-05",
        "TI12DataCaloNu": "FASERNU-CALONU",
        "TI12MCCaloNu": "FASERNU-CALONU",
    }
    if runtype not in mapping:
        raise ValueError(f"Unsupported geometry runtype {runtype!r}; expected one of {sorted(mapping)}")
    return mapping[runtype]


def main() -> int:
    t_start = time.time()

    parser = argparse.ArgumentParser(
        description="Single-step FASER overlay tracking reconstruction and ntupling (zero xAOD)."
    )
    parser.add_argument("file_path", help="Input OverlaidRDO POOL file path")
    parser.add_argument("-o", "--output", default="nt_raw.root", help="Output raw ROOT ntuple path")
    parser.add_argument("-g", "--geom", default="TI12MC04", help="Geometry preset (default: TI12MC04 -> FASERNU-04)")
    parser.add_argument(
        "-c",
        "--cond",
        default="OFLCOND-FASER-05",
        help="Global conditions tag (default: OFLCOND-FASER-05 for Pawan baseline parity)",
    )
    parser.add_argument("--isMC", action="store_true", help="Input is MC overlay (uses OFLP200)")
    parser.add_argument("-n", "--nevents", type=int, default=-1, help="Max events to process (default: -1 all)")
    parser.add_argument("--skip", type=int, default=0, help="Events to skip (default: 0)")
    parser.add_argument(
        "--materialMap",
        default="",
        help="Override ACTS material map JSON/.json.gz path (default: RotatedInput when rotated, Input when unrotated)",
    )
    parser.add_argument(
        "--geometry-rotation",
        type=parse_bool,
        default=True,
        help="Whether the active Calypso build uses rotated ACTS geometry (true/false, default: true)",
    )
    parser.add_argument(
        "--keep-backward",
        action="store_true",
        help="Also run CKF_Back_woIFT (off by default since TwoTrackOverlayNTuplerAlg reads forward CKF only)",
    )
    parser.add_argument(
        "--seed-cov-phi",
        type=float,
        default=None,
        help="Optional override for CircleFitTrackSeedTool.covPhi [rad^2] (default: None -> unmodified CKF2Cfg default)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable VERBOSE Athena output")
    args = parser.parse_args()

    filepath = Path(args.file_path)
    faser_version = resolve_geometry(args.geom)
    filestem = filepath.stem
    if filestem.endswith("-RDO"):
        filestem = filestem[:-4]

    print(
        f"Starting inline overlay reco+ntuple of {filepath.name} "
        f"(geom={args.geom} -> {faser_version}, cond={args.cond}, isMC={args.isMC}, "
        f"geometryRotation={args.geometry_rotation}, keepBackward={args.keep_backward}, "
        f"seedCovPhi={args.seed_cov_phi})"
    )
    if args.nevents > 0:
        print(f"Reconstructing {args.nevents} events by command-line option")
    if args.skip > 0:
        print(f"Skipping {args.skip} events by command-line option")
    if not args.geometry_rotation:
        print(
            "NOTE: --geometry-rotation false selected: using unrotated material source and "
            "legacySwappedCovLoc=True; must be executed with a matching unrotated Calypso build."
        )

    from AthenaCommon.Configurable import Configurable
    from AthenaCommon.Constants import INFO, VERBOSE
    from AthenaCommon.Logging import log
    from AthenaConfiguration.ComponentFactory import CompFactory
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    from CalypsoConfiguration.AllConfigFlags import initConfigFlags
    from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
    from FaserActsKalmanFilter.CKF2Config import CKF2Cfg
    from FaserGeoModel.FaserGeoModelConfig import FaserGeometryCfg
    from TrackerPrepRawDataFormation.TrackerPrepRawDataFormationConfig import (
        FaserSCT_ClusterizationCfg,
    )
    try:
        from TrackerSegmentFit.GhostBustersConfig import GhostBustersCfg
    except ImportError:
        from FaserActsKalmanFilter.GhostBustersConfig import GhostBustersCfg
    from TrackerSegmentFit.TrackerSegmentFitConfig import SegmentFitAlgCfg
    from TrackerSpacePointFormation.TrackerSpacePointFormationConfig import (
        TrackerSpacePointFinderCfg,
    )

    Configurable.configurableRun3Behavior = True

    configFlags = initConfigFlags()
    configFlags.Input.isMC = args.isMC
    configFlags.IOVDb.DatabaseInstance = "OFLP200" if args.isMC else "CONDBR3"
    configFlags.Input.ProjectName = "data20"
    configFlags.Exec.SkipEvents = args.skip
    configFlags.GeoModel.FaserVersion = faser_version
    configFlags.IOVDb.GlobalTag = args.cond

    if args.materialMap:
        material_source = args.materialMap
    elif args.geometry_rotation:
        material_source = "RotatedInput"
    else:
        material_source = "Input"
    configFlags.TrackingGeometry.MaterialSource = material_source
    print(f">>> Using ACTS material source: {material_source} <<<")

    # Keep original string (pathlib collapses double slashes in XRootD URLs).
    configFlags.Input.Files = [args.file_path]

    if not args.isMC:
        # Match faser_reco.py metadata handling for real-data overlay RDOs.
        configFlags.Input.TypedCollections = []
        configFlags.Input.MetadataItems = []

    configFlags.lock()
    if args.verbose:
        configFlags.dump()

    acc = MainServicesCfg(configFlags)
    acc.merge(PoolReadCfg(configFlags))
    acc.merge(FaserGeometryCfg(configFlags))

    # Tracker reconstruction (exact parity with faser_reco.py --isOverlay)
    acc.merge(
        FaserSCT_ClusterizationCfg(
            configFlags,
            DataObjectName="SCT_RDOs",
            checkBadChannels=True,
        )
    )
    acc.merge(TrackerSpacePointFinderCfg(configFlags))
    acc.merge(
        SegmentFitAlgCfg(
            configFlags,
            ReducedChi2Cut=25.0,
            SharedHitFraction=0.61,
            MinClustersPerFit=5,
            TanThetaXZCut=0.083,
        )
    )
    acc.merge(GhostBustersCfg(configFlags))

    ckf_extra = {}
    if not args.geometry_rotation:
        ckf_extra["legacySwappedCovLoc"] = True
    if args.seed_cov_phi is not None:
        ckf_extra["seedToolKwargs"] = {"covPhi": float(args.seed_cov_phi)}

    acc.merge(
        CKF2Cfg(
            configFlags,
            maskedLayers=[0, 1, 2, 3, 4, 5],
            name="CKF_woIFT",
            actsOutputTag=f"{filestem}_3station_forward",
            OutputCollection="CKFTrackCollectionWithoutIFT",
            BackwardPropagation=False,
            noDiagnostics=True,
            **ckf_extra,
        )
    )
    ckf_alg = acc.getEventAlgo("CKF_woIFT")
    print(
        f">>> CKF_woIFT TrackSeed covariance: covLoc0={ckf_alg.TrackSeed.covLoc0}, "
        f"covLoc1={ckf_alg.TrackSeed.covLoc1}, covPhi={ckf_alg.TrackSeed.covPhi}, "
        f"covTheta={ckf_alg.TrackSeed.covTheta}, covQOverP={ckf_alg.TrackSeed.covQOverP} <<<"
    )
    if args.keep_backward:
        acc.merge(
            CKF2Cfg(
                configFlags,
                maskedLayers=[0, 1, 2, 3, 4, 5],
                name="CKF_Back_woIFT",
                actsOutputTag=f"{filestem}_3station_backward",
                OutputCollection="CKFTrackCollectionBackwardWithoutIFT",
                BackwardPropagation=True,
                noDiagnostics=True,
                **ckf_extra,
            )
        )

    # Inline TwoTrackOverlayNTuplerAlg directly on StoreGate (zero xAOD written)
    acc.addEventAlgo(CompFactory.TwoTrackOverlayNTuplerAlg("TwoTrackOverlayNTuplerAlg"))
    thist_svc = CompFactory.THistSvc()
    thist_svc.Output += [f"HIST DATAFILE='{args.output}' OPT='RECREATE'"]
    acc.addService(thist_svc)

    # Exact Data DBReplicaSvc settings matching faser_reco.py:430-435
    if not args.isMC:
        replica_svc = acc.getService("DBReplicaSvc")
        replica_svc.COOLSQLiteVetoPattern = ""
        replica_svc.UseCOOLSQLite = True
        replica_svc.UseCOOLFrontier = False
        replica_svc.UseGeomSQLite = True

    if args.verbose:
        acc.foreach_component("*").OutputLevel = VERBOSE
    else:
        acc.foreach_component("*").OutputLevel = INFO
        event_loop = CompFactory.AthenaEventLoopMgr()
        event_loop.EventPrintoutInterval = 1000
        acc.addService(event_loop)

    acc.foreach_component("*ClassID*").OutputLevel = INFO
    acc.getService("MessageSvc").Format = "% F%40W%S%7W%R%T %0W%M"

    sc = acc.run(maxEvents=args.nevents)

    elapsed = time.time() - t_start
    log.info(f"Finish execution in {elapsed} seconds")
    if sc.isSuccess():
        log.info("Execution succeeded")
        return 0
    log.error("Execution failed, return 1")
    return 1


if __name__ == "__main__":
    sys.exit(main())
