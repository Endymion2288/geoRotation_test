#!/usr/bin/env python3
"""Native Calypso reconstruction with a separate, opt-in slope observer."""
import argparse
import json
import os
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))


def configure_refit(ckf, args, factory, suffix=''):
    """Apply the requested refit control to each actual CKF direction."""
    fitter = ckf.KalmanFitterTool1
    if args.refit_covariance_scale is not None:
        fitter.SeedCovarianceScale = args.refit_covariance_scale
    for option, property_name in [
            ('wrap_refit_smoothing_phi', 'WrapRefitSmoothingPhi'),
            ('canonical_refit_smoothing_angles', 'CanonicalRefitSmoothingAngles'),
            ('canonical_refit_update_angles', 'CanonicalRefitUpdateAngles'),
            ('regular_refit_smoothing_chart', 'RegularRefitSmoothingChart'),
            ('regular_refit_update_chart', 'RegularRefitUpdateChart'),
            ('correct_refit_basis', 'CorrectCovarianceBasis'),
            ('refit_input_reference_plane', 'UseInputReferencePlane')]:
        if getattr(args, option):
            setattr(fitter, property_name, True)
    if args.correct_conversion_basis:
        ckf.CreateTrkTrackTool = factory.CreateTrkTrackTool('SlopeCKFConverter'+suffix,
                                                         CorrectCovarianceBasis=True)
        fitter.CreateTrkTrackTool = factory.CreateTrkTrackTool('SlopeRefitConverter'+suffix,
                                                            CorrectCovarianceBasis=True)
    return {name: getattr(fitter, name, False) for name in [
        'SeedCovarianceScale', 'WrapRefitSmoothingPhi', 'CanonicalRefitSmoothingAngles',
        'CanonicalRefitUpdateAngles', 'RegularRefitSmoothingChart',
        'RegularRefitUpdateChart', 'CorrectCovarianceBasis', 'UseInputReferencePlane']}


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
    parser.add_argument('--charge-audit', action='store_true', help='Retain r/S/K and full C for the frozen small charge replay sample')
    parser.add_argument('--prior-component-audit', action='store_true',
                        help='Retain actual CKF seed and state matrices, including events with no selected component')
    parser.add_argument('--native-states', action='store_true', help='Enable existing Calypso CKF/refit ROOT state writers for both native branches')
    smoothing = parser.add_mutually_exclusive_group()
    smoothing.add_argument('--isolate-smoothing', action='store_true', help='Recompute each selected candidate on its private state path')
    smoothing.add_argument('--shared-smoothing', action='store_true', help='Explicit historical shared-state control')
    parser.add_argument('--physical-ntuple-angles', action='store_true', help='Opt-in physical VetoNu angles in the actual rotated NtupleDumper')
    position = parser.add_mutually_exclusive_group()
    position.add_argument('--truth-cluster-position', action='store_true', help='Calypso diagnostic: replace only local measured strip coordinates by majority SiHit centroids')
    position.add_argument('--sdo-supported-strip-position', action='store_true',
                          help='Calypso diagnostic: average native strip centres supported by majority SDO deposits; retain strip quantization, width and error')
    parser.add_argument('--cluster-coordinate-audit', action='store_true',
                        help='Read-only native cluster channels, actual strip geometry and same-wafer SiHit endpoints')
    parser.add_argument('--measurement-variance-scale', type=float,
                        help='Calypso diagnostic: scale only cluster local-X variance; requires the recorded cluster control library')
    parser.add_argument('--ckf-measurement-variance-scale', type=float,
                        help='Calypso diagnostic: scale actual CKF measurement variance, retaining cluster errors and seed covariance')
    parser.add_argument('--ckf-cluster-measurement-variance', action='store_true',
                        help='Calypso diagnostic: use cluster local-X variance in CKF instead of the native fixed resolution')
    parser.add_argument('--nominal-alignment', action='store_true', help='Use default GeoModel sensor transforms with a real nominal SCT conditions store')
    parser.add_argument('--alignment-audit', action='store_true', help='Write actual and default sensor transforms; requires the recorded conditions control library')
    parser.add_argument('--refit-covariance-scale', type=float, help='Controlled change to native KalmanFitterTool1.SeedCovarianceScale; CKF remains native')
    parser.add_argument('--wrap-refit-smoothing-phi', action='store_true', help='Calypso diagnostic: use a periodic phi subtraction in the actual ACTS refit smoother')
    parser.add_argument('--canonical-refit-smoothing-angles', action='store_true', help='Calypso diagnostic: normalize smoothed phi/theta and transform the covariance on pole reflection')
    parser.add_argument('--canonical-refit-update-angles', action='store_true', help='Calypso diagnostic: normalize the native refit update output and transform theta covariance before transport')
    parser.add_argument('--regular-refit-smoothing-chart', action='store_true', help='Calypso diagnostic: use a beam-regular angular chart for the actual native smoother, transforming all covariance and transport Jacobian blocks')
    parser.add_argument('--regular-refit-update-chart', action='store_true', help='Calypso diagnostic: perform the native measurement update in beam-regular angles and transform filtered mean/covariance back before actual native transport')
    parser.add_argument('--regular-ckf-update-chart', action='store_true',
                        help='Calypso diagnostic: apply the existing beam-regular direction adapter to actual CKF updates')
    parser.add_argument('--regular-ckf-smoothing-chart', action='store_true',
                        help='Calypso diagnostic: apply the existing beam-regular direction adapter to CKF and isolated smoothing')
    parser.add_argument('--correct-conversion-basis', action='store_true', help='Calypso diagnostic: transform ACTS surface covariance to the actual Trk curvilinear local basis')
    parser.add_argument('--correct-refit-basis', action='store_true', help='Calypso diagnostic: transform input curvilinear covariance to the actual refit plane basis')
    parser.add_argument('--refit-input-reference-plane', action='store_true', help='Calypso diagnostic: keep the refit plane at the input point, independently of covariance changes')
    parser.add_argument('--refit-outcomes', action='store_true', help='Read-only actual Calypso selected-track refit success/fallback CSV')
    parser.add_argument('--charge-diagnostics', action='store_true', help='Read-only seed, candidate, selected state and conversion/refit lineage')
    parser.add_argument('--mode-refit-counters', action='store_true',
                        help='Read-only compact refit counters for every native direction, including merge')
    parser.add_argument('--backward', action='store_true', help='Actual Calypso backward CKF path')
    parser.add_argument('--four-station', action='store_true', help='Use the native four-station mask [0,1,2]')
    parser.add_argument('--merge-directions', action='store_true', help='Run native forward/backward CKF and MergeTracksAlg')
    parser.add_argument('--match-faser-seed-prior', action='store_true', help='Use the beforeGeo seed covariance transformed with the actual surface and angular Jacobians')
    parser.add_argument('--seed-prior-audit', action='store_true', help='Read-only full matched seed covariance and Jacobian CSV')
    parser.add_argument('--truth-source-link-barcode', type=int, default=0,
                        help='Actual CKF control: retain source links supported by this SDO barcode')
    parser.add_argument('--reverse-charge-and-field', action='store_true',
                        help='Reverse seed q/p and the CKF/refit/extrapolation field together')
    parser.add_argument('--seed-override-file',
                        help='Actual CKF seed overrides: local event index, seed index, six bound parameters per row; preserve covariance')
    parser.add_argument('--seed-override-mode', choices=['qop', 'mean'], default='qop',
                        help='Replace q/p only or the complete mean; requires the recorded Calypso hit-control library')
    parser.add_argument('--full-seed-input-file',
                        help='Actual CKF diagnostic: complete seed mean, covariance and reference-plane transform')
    parser.add_argument('--seed-component-input-file',
                        help='Actual CKF diagnostic: independent native finder/update/smoother calls for finite prior components')
    parser.add_argument('--component-refit', action='store_true',
                        help='Extend each finite-prior component through the native converter and refit')
    parser.add_argument('--run-only-input-seeds', action='store_true',
                        help='Run the native CKF only for seed indices listed in the full seed input')
    parser.add_argument('--source-link-selection-file',
                        help='Actual CKF diagnostic: admit only the listed original cluster identifiers')
    parser.add_argument('--ckf-no-scattering', action='store_true',
                        help='Disable only the native CKF material scattering term; refit and observer retain their native settings')
    parser.add_argument('--ckf-no-energy-loss', action='store_true',
                        help='Disable only the native CKF material energy-loss term; refit and observer retain their native settings')
    args = parser.parse_args()
    if args.charge_audit and not args.stage_audit:
        parser.error('--charge-audit requires --stage-audit')
    if args.prior_component_audit and args.stage_audit:
        parser.error('prior-component and stage audit use separate output reductions')
    if args.mode_refit_counters and args.charge_diagnostics:
        parser.error('--mode-refit-counters and --charge-diagnostics use different compact outputs')
    if args.full_seed_input_file and args.seed_override_file:
        parser.error('full seed input and mean-only overrides are separate controls')
    if args.seed_component_input_file and (
            not args.prior_component_audit or not args.isolate_smoothing or
            args.full_seed_input_file or args.seed_override_file or args.run_only_input_seeds):
        parser.error('seed components require prior-component audit and isolated smoothing, without other seed overrides')
    if args.component_refit and not args.seed_component_input_file:
        parser.error('--component-refit requires --seed-component-input-file')
    if (args.regular_ckf_update_chart or args.regular_ckf_smoothing_chart) and args.variant != 'beforeGeo':
        parser.error('the CKF direction-chart adapter is defined for the native beforeGeo frame')
    if args.run_only_input_seeds and not args.full_seed_input_file:
        parser.error('--run-only-input-seeds requires --full-seed-input-file')
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
    if args.nominal_alignment:
        flags.GeoModel.Align.Disable = True
    flags.TrackingGeometry.MaterialSource = args.material
    flags.Exec.SkipEvents = args.skip
    flags.Concurrency.NumThreads = 0
    flags.Concurrency.NumConcurrentEvents = 1
    flags.Output.doWriteESD = False
    flags.lock()
    acc = MainServicesCfg(flags)
    acc.merge(PoolReadCfg(flags))
    acc.merge(FaserGeometryCfg(flags))
    if args.nominal_alignment or args.alignment_audit:
        from calypso_alignment import configure_alignment
        configure_alignment(acc, args.nominal_alignment, 'alignment.csv' if args.alignment_audit else None)
    cluster=FaserSCT_ClusterizationCfg(flags, DataObjectName='SCT_RDOs', checkBadChannels=True)
    if args.truth_cluster_position:
        cluster.getEventAlgo('Tracker::FaserSCT_Clusterization').SCT_ClusteringTool.globalPosAlg.UseTruthClusterPosition=True
    if args.sdo_supported_strip_position:
        cluster.getEventAlgo('Tracker::FaserSCT_Clusterization').SCT_ClusteringTool.globalPosAlg.UseSDOSupportedStripPosition=True
    if args.cluster_coordinate_audit:
        cluster.getEventAlgo('Tracker::FaserSCT_Clusterization').SCT_ClusteringTool.globalPosAlg.CoordinateAuditFile='cluster_coordinates.jsonl'
    if args.measurement_variance_scale is not None:
        cluster.getEventAlgo('Tracker::FaserSCT_Clusterization').SCT_ClusteringTool.globalPosAlg.MeasurementVarianceScale=args.measurement_variance_scale
    acc.merge(cluster)
    acc.merge(TrackerSpacePointFinderCfg(flags))
    acc.merge(SegmentFitAlgCfg(flags, ReducedChi2Cut=25., SharedHitFraction=.61,
                             MinClustersPerFit=5, TanThetaXZCut=.083))
    acc.merge(GhostBustersCfg(flags))
    acc.merge(xAODTruthCnvAlgCfg(flags))
    mask = [0, 1, 2] if args.four_station else [0, 1, 2, 3, 4, 5]
    if args.merge_directions and (args.backward or args.charge_diagnostics or args.stage_audit or args.prior_component_audit or args.native_states
                                  or args.truth_source_link_barcode or args.reverse_charge_and_field
                                  or args.seed_override_file or args.full_seed_input_file
                                  or args.source_link_selection_file or args.ckf_no_scattering
                                  or args.ckf_no_energy_loss):
        parser.error('--merge-directions uses two native CKF paths; single-path audit flags do not apply')
    primary_tag = f'slope_{4 if args.four_station else 3}station_'+('backward' if args.backward else 'forward')
    cfg = CKF2Cfg(flags, actsOutputTag=primary_tag, name='CKF_Slope',
        maskedLayers=mask, OutputCollection='CKFTrackCollectionWithoutIFT',
        BackwardPropagation=args.backward, noDiagnostics=False, seedWriter=not args.no_seed_writer,
        StatesWriter=args.native_states, SummaryWriter=False)
    ckf = cfg.getEventAlgo('CKF_Slope')
    if args.ckf_measurement_variance_scale is not None:
        ckf.TrackSeed.MeasurementVarianceScale=args.ckf_measurement_variance_scale
    if args.ckf_cluster_measurement_variance:
        ckf.TrackSeed.UseClusterMeasurementVariance=True
    if args.seed_override_file:
        ckf.SeedOverrideFile = args.seed_override_file
        ckf.SeedOverrideMode = args.seed_override_mode
    if args.full_seed_input_file:
        ckf.FullSeedInputFile = args.full_seed_input_file
        ckf.RunOnlyInputSeeds = args.run_only_input_seeds
    if args.seed_component_input_file:
        ckf.SeedComponentInputFile = args.seed_component_input_file
    if args.component_refit:
        ckf.ComponentRefit = True
    if args.regular_ckf_update_chart:
        ckf.RegularCkfUpdateChart = True
    if args.regular_ckf_smoothing_chart:
        ckf.RegularCkfSmoothingChart = True
    if args.source_link_selection_file:
        ckf.SourceLinkSelectionFile = args.source_link_selection_file
    if args.ckf_no_scattering:
        ckf.MaterialScattering = False
    if args.ckf_no_energy_loss:
        ckf.MaterialEnergyLoss = False
    if args.truth_source_link_barcode:
        ckf.TruthSourceLinkBarcode = args.truth_source_link_barcode
    if args.reverse_charge_and_field:
        ckf.ReverseChargeHypothesis = True
        ckf.MagneticFieldScale = -1.
        ckf.KalmanFitterTool1.MagneticFieldScale = -1.
    if args.match_faser_seed_prior:
        if args.variant == 'beforeGeo':
            parser.error('--match-faser-seed-prior requires the rotated Calypso control library')
        # Exact beforeGeo CKF2Config.py expressions, including its local swap.
        for name, inflation, sigma in zip(
                ['covLoc0', 'covLoc1', 'covPhi', 'covTheta', 'covQOverP'],
                [100, 100, 100, 100, 1000], [.9, .019, .033, .00029, .0001]):
            setattr(ckf.TrackSeed, name, inflation*sigma*sigma)
        ckf.TrackSeed.MatchFaserSeedCovariance = True
    if args.seed_prior_audit:
        if not args.match_faser_seed_prior:
            parser.error('--seed-prior-audit requires --match-faser-seed-prior')
        ckf.TrackSeed.SeedPriorAuditFile = 'seed_covariance.csv'
    if args.charge_diagnostics:
        ckf.ChargeDiagnosticFile = 'charge_lineage.csv'
    if args.mode_refit_counters:
        ckf.ChargeDiagnosticFile = 'mode_refit_forward.csv'
    if args.stage_audit or args.prior_component_audit:
        ckf.AuditFile = 'stages.jsonl'
        ckf.AuditStates = True
        ckf.AuditSelectedStatesOnly = True
    if args.isolate_smoothing:
        ckf.IsolateCandidateSmoothing = True
    if args.shared_smoothing:
        ckf.IsolateCandidateSmoothing = False
    ckf.KalmanFitterTool1.noDiagnostics = not args.native_states
    ckf.KalmanFitterTool1.StatesWriter = args.native_states
    ckf.KalmanFitterTool1.RootTrajectoryStatesWriterTool.noDiagnostics = not args.native_states
    ckf.KalmanFitterTool1.SummaryWriter = False
    refit_controls = {'backward' if args.backward else 'forward':
                      configure_refit(ckf, args, CompFactory)}
    if args.refit_outcomes:
        ckf.RefitOutcomeFile='refit_outcomes.csv'
    acc.merge(cfg)
    observer_collection = 'CKFTrackCollectionWithoutIFT'
    if args.merge_directions:
        backward_cfg = CKF2Cfg(flags, actsOutputTag='slope_native_backward', name='CKF_Slope_Backward',
            maskedLayers=mask, OutputCollection='CKFTrackCollectionBackwardWithoutIFT',
            BackwardPropagation=True, noDiagnostics=True, seedWriter=False,
            StatesWriter=False, SummaryWriter=False)
        backward_ckf = backward_cfg.getEventAlgo('CKF_Slope_Backward')
        if args.regular_ckf_update_chart:
            backward_ckf.RegularCkfUpdateChart = True
        if args.regular_ckf_smoothing_chart:
            backward_ckf.RegularCkfSmoothingChart = True
        if args.isolate_smoothing:
            backward_ckf.IsolateCandidateSmoothing=True
        if args.shared_smoothing:
            backward_ckf.IsolateCandidateSmoothing=False
        if args.mode_refit_counters:
            backward_ckf.ChargeDiagnosticFile = 'mode_refit_backward.csv'
        refit_controls['backward'] = configure_refit(backward_ckf, args, CompFactory, 'Backward')
        acc.merge(backward_cfg)
        from FaserActsKalmanFilter.MergeTracksConfig import MergeTrackCfg
        acc.merge(MergeTrackCfg(flags, TrackCollectionForward='CKFTrackCollectionWithoutIFT',
            TrackCollectionBackward='CKFTrackCollectionBackwardWithoutIFT',
            MergeTrackCollection='CKFTrackCollectionMergeWithoutIFT', acceptTrackSeenByOne=True))
        observer_collection = 'CKFTrackCollectionMergeWithoutIFT'
    result, geometry = ActsTrackingGeometryToolCfg(flags)
    acc.merge(result)
    propagator = CompFactory.FaserActsExtrapolationTool('SlopeNativeExtrapolator',
        TrackingGeometryTool=geometry, MaxSteps=10000)
    if args.reverse_charge_and_field:
        propagator.MagneticFieldScale = -1.
    if not args.no_observer:
        observer_options = {'TruthChargeScale': -1.} if args.reverse_charge_and_field else {}
        acc.addEventAlgo(CompFactory.CalypsoSlopeDiagnosticsAlg('CalypsoSlopeDiagnostics',
            ExtrapolationTool=propagator, Rotated=args.variant != 'beforeGeo', Skip=args.skip,
            ActualCovariance=args.actual_covariance, TrackCollection=observer_collection, **observer_options))
    if args.ntuple_check:
        from ntuple_tracking_check import EmptyWaveforms
        from NtupleDumper.NtupleDumperConfig import NtupleDumperAlgCfg
        acc.addEventAlgo(EmptyWaveforms('SlopeEmptyWaveforms'))
        ntuple_options = {'RotatedPhysicalAngles':True} if args.physical_ntuple_angles else {}
        acc.merge(NtupleDumperAlgCfg(flags, 'native_ntuple.root', McEventCollection='TruthEvent',
            TrackCollection=observer_collection, TrackCollectionWithoutIFT=observer_collection,
            RunVertexing=False, DoBlinding=False, DoTrackFilter=False, StableOnly=False, DoCalib=False, **ntuple_options))
    acc.addService(CompFactory.AthenaEventLoopMgr(EventPrintoutInterval=1000))
    acc.getService('MessageSvc').defaultLimit = 20
    manifest = dict(args=vars(args), geometry='FASERNU-04', conditions='OFLCOND-FASER-05',
        calypso_dir=os.environ.get('Calypso_DIR'), library=os.environ.get('GEOROTATION_LIBRARY_DIR'),
        seed_covariance={p: getattr(ckf.TrackSeed, p) for p in ['covLoc0','covLoc1','covPhi','covTheta','covQOverP']},
        refit_covariance_scale=ckf.KalmanFitterTool1.SeedCovarianceScale,
        refit_smoothing_phi_periodic=args.wrap_refit_smoothing_phi,
        refit_smoothing_angles_canonical=args.canonical_refit_smoothing_angles,
        refit_update_angles_canonical=args.canonical_refit_update_angles,
        refit_smoothing_chart_regular=args.regular_refit_smoothing_chart,
        refit_update_chart_regular=args.regular_refit_update_chart,
        ckf_update_chart_regular=args.regular_ckf_update_chart,
        ckf_smoothing_chart_regular=args.regular_ckf_smoothing_chart,
        actual_refit_controls_by_direction=refit_controls,
        conversion_covariance_basis_corrected=args.correct_conversion_basis,
        refit_covariance_basis_corrected=args.correct_refit_basis,
        refit_uses_input_reference_plane=args.refit_input_reference_plane,
        refit_outcome_observer=args.refit_outcomes,
        charge_diagnostics=args.charge_diagnostics,
        mode_refit_counters=args.mode_refit_counters,
        matched_faser_seed_prior=args.match_faser_seed_prior,
        full_seed_prior_audit=args.seed_prior_audit,
        truth_source_link_barcode=args.truth_source_link_barcode,
        magnetic_field_scale=-1. if args.reverse_charge_and_field else 1.,
        seed_charge_reversed=args.reverse_charge_and_field,
        seed_override_file=args.seed_override_file, seed_override_mode=args.seed_override_mode,
        full_seed_input_file=args.full_seed_input_file, run_only_input_seeds=args.run_only_input_seeds,
        seed_component_input_file=args.seed_component_input_file,
        component_refit=args.component_refit,
        source_link_selection_file=args.source_link_selection_file,
        ckf_material_scattering=getattr(ckf, 'MaterialScattering', True),
        ckf_material_energy_loss=getattr(ckf, 'MaterialEnergyLoss', True),
        backward=args.backward, masked_layers=mask, merge_directions=args.merge_directions,
        output_collection=observer_collection,
        audit=args.stage_audit, native_states=args.native_states, reset_event_maps=False,
        prior_component_audit=args.prior_component_audit,
        charge_matrix_audit=args.charge_audit,
        isolate_smoothing=bool(getattr(ckf, 'IsolateCandidateSmoothing', False)), observer=not args.no_observer,
        propagation_material_energy_loss=propagator.InteractionEloss,
        propagation_material_scattering=propagator.InteractionMultiScatering)
    manifest['alignment'] = dict(nominal_alignment=args.nominal_alignment,
        conditions_displacements_disabled=bool(flags.GeoModel.Align.Disable),
        tracker_alignment_folders=[folder for folder in acc.getService('IOVDbSvc').Folders
                                   if '/Tracker/Align' in folder])
    Path('run.json').write_text(json.dumps(manifest, indent=2)+'\n')
    begin = time.monotonic()
    status = acc.run(maxEvents=args.events)
    mappings = Path('/proc/self/maps').read_text().splitlines()
    loaded = sorted({line.split()[-1] for line in mappings
                     if any(name in line for name in ['libFaserActsGeometry',
                         'libFaserActsKalmanFilter', 'libCalypsoSlopeDiagnostics'])})
    manifest.update(success=status.isSuccess(), wall_s=time.monotonic()-begin,
                    loaded_acts_libraries=loaded,
                    loaded_cluster_libraries=sorted({line.split()[-1] for line in mappings
                        if 'libFaserSiClusterizationTool' in line and '.so' in line}),
                    loaded_conditions_libraries=sorted({line.split()[-1] for line in mappings
                        if 'libFaserSCT_ConditionsAlgorithms' in line and '.so' in line}))
    Path('run.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0 if status.isSuccess() else 1


if __name__ == '__main__':
    raise SystemExit(main())
