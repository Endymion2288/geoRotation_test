#!/usr/bin/env python3
"""Tracking-only reconstruction with explicit conditions and compact event output."""
import argparse
import json
import os
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dsid", required=True)
    parser.add_argument("--file-index", type=int, default=0)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--events", type=int, default=10)
    parser.add_argument("--geometry", default="FASERNU-04")
    parser.add_argument("--conditions", default="OFLCOND-FASER-05")
    parser.add_argument("--material", default="/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json")
    parser.add_argument("--diagnostics", action="store_true")
    parser.add_argument("--audit", action="store_true")
    parser.add_argument("--states", action="store_true")
    parser.add_argument("--closure", action="store_true")
    parser.add_argument("--correct-interfaces", action="store_true")
    parser.add_argument("--reset-event-maps", action="store_true")
    parser.add_argument("--compact-only", action="store_true", help="Use the C++ audit without the independent Python observer")
    args = parser.parse_args()
    started = time.monotonic()
    samples = json.loads((REPO / "configs/samples.json").read_text())
    sample = next(s for s in samples["samples"] if s["dsid"] == args.dsid)
    input_file = sample["files"][args.file_index]["path"]

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

    Configurable.configurableRun3Behavior = True
    flags = initConfigFlags()
    flags.Input.isMC = True
    flags.Input.ProjectName = "mc24"
    flags.Input.Files = [input_file]
    flags.IOVDb.DatabaseInstance = "OFLP200"
    flags.GeoModel.FaserVersion = args.geometry
    flags.IOVDb.GlobalTag = args.conditions
    flags.TrackingGeometry.MaterialSource = args.material
    flags.Exec.SkipEvents = args.skip
    flags.Concurrency.NumThreads = 0
    flags.Concurrency.NumConcurrentEvents = 1
    flags.Output.doWriteESD = False
    flags.lock()
    acc = MainServicesCfg(flags)
    acc.merge(PoolReadCfg(flags))
    acc.merge(FaserGeometryCfg(flags))
    acc.merge(FaserSCT_ClusterizationCfg(flags, DataObjectName="SCT_RDOs", checkBadChannels=True))
    acc.merge(TrackerSpacePointFinderCfg(flags))
    acc.merge(SegmentFitAlgCfg(flags, ReducedChi2Cut=25., SharedHitFraction=0.61,
                               MinClustersPerFit=5, TanThetaXZCut=0.083))
    acc.merge(GhostBustersCfg(flags))
    acc.merge(xAODTruthCnvAlgCfg(flags))
    kwargs = dict(name="CKF_Audit", maskedLayers=[0, 1, 2, 3, 4, 5],
                  OutputCollection="CKFTrackCollectionWithoutIFT", BackwardPropagation=False,
                  noDiagnostics=not args.diagnostics, seedWriter=args.diagnostics,
                  SummaryWriter=args.diagnostics, StatesWriter=args.states)
    if args.audit:
        kwargs.update(AuditFile="stages.jsonl", AuditStates=args.states, AuditClosure=args.closure)
    cfg = CKF2Cfg(flags, actsOutputTag="forward3", **kwargs)
    ckf = cfg.getEventAlgo("CKF_Audit")
    if args.reset_event_maps:
        ckf.TrackSeed.ResetEventMaps = True
    fitter = ckf.KalmanFitterTool1
    # Distinct names prevent the two private writers from being merged.
    for prop, name, path in [
        ("RootTrajectoryStatesWriterTool", "RefitStates", "refit_states.root"),
        ("RootTrajectorySummaryWriterTool", "RefitSummary", "refit_summary.root"),
    ]:
        writer = getattr(fitter, prop)
        writer = writer.clone(name)
        writer.FilePath = path
        writer.noDiagnostics = not args.diagnostics
        setattr(fitter, prop, writer)
    fitter.noDiagnostics = not args.diagnostics
    fitter.StatesWriter = args.states
    fitter.SummaryWriter = args.diagnostics
    if args.correct_interfaces:
        ckf.CreateTrkTrackTool = CompFactory.CreateTrkTrackTool("CorrectedConverter", CorrectCovariance=True)
        fitter.CreateTrkTrackTool = CompFactory.CreateTrkTrackTool("CorrectedRefitConverter", CorrectCovariance=True)
        fitter.CorrectInterfaces = True
    acc.merge(cfg)
    if not args.compact_only:
        from event_audit import CompactEventWriter
        acc.addEventAlgo(CompactEventWriter(name="CompactEventWriter", output="events.jsonl"))
    loop = CompFactory.AthenaEventLoopMgr(EventPrintoutInterval=1000)
    acc.addService(loop)
    acc.getService("MessageSvc").defaultLimit = 20
    acc.getService("MessageSvc").Format = "% F%40W%S%7W%R%T %0W%M"
    manifest = {"schema_version": 1, "dsid": args.dsid, "input": input_file,
                "file_index": args.file_index, "args": vars(args), "masked_layers": [0, 1, 2, 3, 4, 5],
                "calypso_dir": os.environ.get("Calypso_DIR"), "athena": os.environ.get("AtlasVersion"),
                "threads": 1, "collection": "CKFTrackCollectionWithoutIFT"}
    manifest["material_source"] = os.environ.get("GEOROTATION_MATERIAL_SOURCE", args.material)
    Path("run.json").write_text(json.dumps(manifest, indent=2) + "\n")
    sc = acc.run(maxEvents=args.events)
    manifest.update(success=sc.isSuccess(), wall_seconds=time.monotonic() - started)
    Path("run.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return 0 if sc.isSuccess() else 1


if __name__ == "__main__":
    sys.exit(main())
