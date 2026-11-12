#!/usr/bin/env python3
"""Tracking-only reconstruction with explicit conditions and compact event output."""
import argparse
import json
import math
import os
import socket
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dsid", required=True)
    parser.add_argument("--file-index", type=int, default=0)
    parser.add_argument("--input", help="Explicit frozen input, e.g. a reserved validation file")
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
    parser.add_argument("--correct-covariance", action="store_true")
    parser.add_argument("--correct-refit-surface", action="store_true")
    parser.add_argument("--selected-states-only", action="store_true")
    parser.add_argument("--isolate-smoothing", action="store_true")
    parser.add_argument("--no-energy-loss", action="store_true")
    parser.add_argument("--no-scattering", action="store_true")
    parser.add_argument("--bypass-refit", action="store_true")
    parser.add_argument("--fixed-truth-hits", action="store_true")
    parser.add_argument("--electron-hypothesis", action="store_true")
    parser.add_argument("--seed-override", type=str)
    parser.add_argument("--seed-override-mode", choices=["qop", "mean"], default="mean")
    parser.add_argument("--max-step", type=float)
    parser.add_argument("--measurement-chi2", type=float)
    parser.add_argument("--max-branches", type=int)
    parser.add_argument("--seed-covariance-scale", type=float)
    parser.add_argument("--refit-diagonal-scale", type=float)
    parser.add_argument("--refit-covariance-scale", type=float)
    parser.add_argument("--backward", action="store_true")
    parser.add_argument("--include-ift", action="store_true")
    parser.add_argument("--reset-event-maps", action="store_true")
    parser.add_argument("--compact-only", action="store_true", help="Use the C++ audit without the independent Python observer")
    parser.add_argument("--ntuple-check", action="store_true", help="Run the actual NtupleDumper and compare its momentum branches with audit final tracks")
    args = parser.parse_args()
    started = time.monotonic()
    samples = json.loads((REPO / "configs/samples.json").read_text())
    sample = next(s for s in samples["samples"] if s["dsid"] == args.dsid)
    input_file = args.input or sample["files"][args.file_index]["path"]

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
    masked_layers = [] if args.include_ift else [0, 1, 2, 3, 4, 5]
    kwargs = dict(name="CKF_Audit", maskedLayers=masked_layers,
                  OutputCollection="CKFTrackCollectionWithoutIFT", BackwardPropagation=args.backward,
                  noDiagnostics=not args.diagnostics, seedWriter=args.diagnostics,
                  SummaryWriter=args.diagnostics, StatesWriter=args.states)
    if args.audit:
        kwargs.update(AuditFile="stages.jsonl", AuditStates=args.states, AuditClosure=args.closure)
    cfg = CKF2Cfg(flags, actsOutputTag="forward3", **kwargs)
    ckf = cfg.getEventAlgo("CKF_Audit")
    if args.reset_event_maps:
        ckf.TrackSeed.ResetEventMaps = True
    fitter = ckf.KalmanFitterTool1
    for enabled, prop in [(args.selected_states_only, "AuditSelectedStatesOnly"),
                          (args.isolate_smoothing, "IsolateCandidateSmoothing"),
                          (args.bypass_refit, "BypassRefit"), (args.fixed_truth_hits, "FixedTruthHits")]:
        if enabled:setattr(ckf, prop, True)
    for enabled, prop, value in [(args.no_energy_loss,"MaterialEnergyLoss",False),
                                (args.no_scattering,"MaterialScattering",False),
                                (args.electron_hypothesis,"ElectronHypothesis",True)]:
        if enabled:
            setattr(ckf, prop, value);setattr(fitter, prop, value)
    if args.max_step is not None:
        ckf.MaxStepSize=args.max_step;fitter.MaxStepSize=args.max_step
    if args.measurement_chi2 is not None:ckf.chi2Max=args.measurement_chi2
    if args.max_branches is not None:ckf.nMax=args.max_branches
    if args.seed_override:
        ckf.SeedOverrideFile=args.seed_override;ckf.SeedOverrideMode=args.seed_override_mode
    if args.seed_covariance_scale is not None:
        for prop in ["covLoc0","covLoc1","covPhi","covTheta","covQOverP"]:
            setattr(ckf.TrackSeed,prop,getattr(ckf.TrackSeed,prop)*args.seed_covariance_scale)
    if args.refit_diagonal_scale is not None:fitter.SeedCovarianceScale=args.refit_diagonal_scale
    if args.refit_covariance_scale is not None:fitter.RefitCovarianceMultiplier=args.refit_covariance_scale
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
    if args.correct_covariance:
        ckf.CreateTrkTrackTool = CompFactory.CreateTrkTrackTool("ChartConverter", CorrectCovariance=True)
        fitter.CreateTrkTrackTool = CompFactory.CreateTrkTrackTool("RefitChartConverter", CorrectCovariance=True)
        fitter.CorrectRefitCovariance = True
    if args.correct_refit_surface:fitter.CorrectRefitSurface = True
    acc.merge(cfg)
    if args.ntuple_check:
        from ntuple_tracking_check import EmptyWaveforms
        acc.addEventAlgo(EmptyWaveforms(name='TrackingOnlyWaveforms'))
        from NtupleDumper.NtupleDumperConfig import NtupleDumperAlgCfg
        acc.merge(NtupleDumperAlgCfg(flags,'output_check.root',McEventCollection='TruthEvent',
            TrackCollection='CKFTrackCollectionWithoutIFT',TrackCollectionWithoutIFT='CKFTrackCollectionWithoutIFT',
            RunVertexing=False,DoBlinding=False,DoTrackFilter=False,StableOnly=False,DoCalib=False))
    if not args.compact_only:
        from event_audit import CompactEventWriter
        acc.addEventAlgo(CompactEventWriter(name="CompactEventWriter", output="events.jsonl"))
    loop = CompFactory.AthenaEventLoopMgr(EventPrintoutInterval=1000)
    acc.addService(loop)
    acc.getService("MessageSvc").defaultLimit = 20
    acc.getService("MessageSvc").Format = "% F%40W%S%7W%R%T %0W%M"
    manifest = {"schema_version": 1, "dsid": args.dsid, "input": input_file,
                "file_index": args.file_index, "args": vars(args), "masked_layers": masked_layers,
                "calypso_dir": os.environ.get("Calypso_DIR"), "athena": os.environ.get("AtlasVersion"),
                "threads": 1, "collection": "CKFTrackCollectionWithoutIFT"}
    manifest["material_source"] = os.environ.get("GEOROTATION_MATERIAL_SOURCE", args.material)
    manifest["hostname"] = socket.gethostname()
    manifest["software_manifest"] = os.environ.get("GEOROTATION_SOFTWARE_MANIFEST")
    manifest["library_override"] = os.environ.get("GEOROTATION_LIBRARY_DIR")
    Path("run.json").write_text(json.dumps(manifest, indent=2) + "\n")
    sc = acc.run(maxEvents=args.events)
    manifest.update(success=sc.isSuccess(), wall_seconds=time.monotonic() - started)
    Path("run.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if sc.isSuccess() and args.ntuple_check:
        import ROOT
        source=ROOT.TFile.Open('output_check.root')
        tree=next((source.Get(key.GetName()) for key in source.GetListOfKeys()
                   if source.Get(key.GetName()).InheritsFrom('TTree')),None)
        assert tree is not None,'NtupleDumper produced no TTree'
        audits=[];event=[]
        with Path('stages.jsonl').open() as stream:
            for line in stream:
                row=json.loads(line)
                if row['stage']=='event' and event:audits.append(event);event=[]
                event.append(row)
        if event:audits.append(event)
        assert tree.GetEntries()==len(audits),(tree.GetEntries(),len(audits))
        check=dict(events=len(audits),ntuple_tracks=0,good95_final_tracks=0,max_norm_relative_error=0.,max_final_relative_error=0.)
        for ordinal,records in enumerate(audits):
            tree.GetEntry(ordinal)
            finals=[r for r in records if r['stage']=='final']
            copies=[]
            for i,p in enumerate(tree.Track_p0):
                vector=[float(getattr(tree,'Track_'+c+'0')[i]) for c in ['px','py','pz']]
                norm=math.sqrt(sum(x*x for x in vector))
                check['max_norm_relative_error']=max(check['max_norm_relative_error'],abs(float(p)/norm-1))
                matches=[r for r in finals if all(math.isclose(a,b*.001,rel_tol=2e-6,abs_tol=1e-8)
                                               for a,b in zip(r['momentum_A_GeV'],vector))]
                assert matches,('Ntuple momentum has no final-track source',ordinal,vector)
                copies.extend(matches)
                expected=math.sqrt(sum(x*x for x in matches[0]['momentum_A_GeV']))
                check['max_final_relative_error']=max(check['max_final_relative_error'],abs(float(p)*.001/expected-1))
                check['ntuple_tracks']+=1
            for r in finals:
                if r['ndof']>0 and r['chi2'] is not None and r['chi2']/r['ndof']<25 and r['layers']>=7 and r['ndof']+5>=12 and r['rmax_mm']<95:
                    check['good95_final_tracks']+=1
                    assert r in copies,('Good95 final missing from Ntuple',ordinal,r['candidate'])
        assert check['max_norm_relative_error']<2e-6 and check['max_final_relative_error']<2e-6
        Path('ntuple_check.json').write_text(json.dumps(check,indent=2)+'\n')
        source.Close()
    return 0 if sc.isSuccess() else 1


if __name__ == "__main__":
    sys.exit(main())
