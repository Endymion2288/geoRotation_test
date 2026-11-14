#!/usr/bin/env python3
"""Transport identical finite physical prior samples with native Calypso tools."""
import argparse
import json
import os
from pathlib import Path
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--material', required=True)
    parser.add_argument('--variant', choices=['beforeGeo', 'afterGeo'], required=True)
    parser.add_argument('--states', required=True)
    parser.add_argument('--skip', type=int, default=0)
    parser.add_argument('--events', type=int, default=20)
    parser.add_argument('--disable-material', action='store_true')
    parser.add_argument('--trace-material', action='store_true')
    parser.add_argument('--max-step-size-m', type=float,
                        help='Explicit native extrapolator step limit in metres')
    args = parser.parse_args()
    from AthenaCommon.Configurable import Configurable
    from AthenaConfiguration.ComponentFactory import CompFactory
    from CalypsoConfiguration.AllConfigFlags import initConfigFlags
    from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    from FaserGeoModel.FaserGeoModelConfig import FaserGeometryCfg
    from FaserActsGeometry.ActsGeometryConfig import ActsTrackingGeometryToolCfg
    from MagFieldServices.MagFieldServicesConfig import MagneticFieldSvcCfg
    Configurable.configurableRun3Behavior = True
    flags = initConfigFlags()
    flags.Input.Files = [args.input]
    flags.Input.isMC = True
    flags.Input.ProjectName = 'mc24'
    flags.IOVDb.DatabaseInstance = 'OFLP200'
    flags.IOVDb.GlobalTag = 'OFLCOND-FASER-05'
    flags.GeoModel.FaserVersion = 'FASERNU-04'
    flags.TrackingGeometry.MaterialSource = args.material
    flags.Exec.SkipEvents = args.skip
    flags.Concurrency.NumThreads = 0
    flags.Concurrency.NumConcurrentEvents = 1
    flags.lock()
    acc = MainServicesCfg(flags)
    acc.merge(PoolReadCfg(flags))
    acc.merge(FaserGeometryCfg(flags))
    acc.merge(MagneticFieldSvcCfg(flags))
    cfg, geometry = ActsTrackingGeometryToolCfg(flags)
    acc.merge(cfg)
    propagator = CompFactory.FaserActsExtrapolationTool('ChargePriorExtrapolator',
        TrackingGeometryTool=geometry, MaxSteps=10000,
        InteractionEloss=not args.disable_material,
        InteractionMultiScatering=not args.disable_material,
        InteractionRecord=args.trace_material)
    if args.max_step_size_m is not None:
        propagator.MaxStepSize = args.max_step_size_m
    acc.addEventAlgo(CompFactory.CalypsoSlopeTransportProbeAlg('ChargePriorTransport',
        ExtrapolationTool=propagator, Rotated=args.variant != 'beforeGeo',
        DirectionScan=False, InputStates=args.states, Skip=args.skip,
        MaterialTraceOutput='material_trace.csv' if args.trace_material else ''))
    acc.getService('MessageSvc').defaultLimit = 20
    manifest = dict(args=vars(args), geometry='FASERNU-04', conditions='OFLCOND-FASER-05',
        library=os.environ.get('GEOROTATION_LIBRARY_DIR'),
        interaction_energy_loss=propagator.InteractionEloss,
        interaction_scattering=propagator.InteractionMultiScatering,
        material_trace=args.trace_material,
        max_step_size_m=propagator.MaxStepSize,
        prior='Native bound inputs: legacy kinds carry no covariance; kind 2 supplies complete '
              '6x6 bound covariance and the actual target geometry ID. Input files define the mode.',
        scope='Actual Calypso extrapolation only. A tangent covariance coordinate map does not '
              'establish finite-distribution or CKF posterior equivalence.')
    Path('run.json').write_text(json.dumps(manifest, indent=2)+'\n')
    begin = time.monotonic()
    status = acc.run(args.events)
    mappings = Path('/proc/self/maps').read_text().splitlines()
    manifest.update(success=status.isSuccess(), wall_s=time.monotonic()-begin,
        loaded_acts_libraries=sorted({line.split()[-1] for line in mappings if any(
            name in line for name in ['libFaserActsGeometry', 'libFaserActsKalmanFilter', 'libCalypsoSlopeDiagnostics'])}))
    Path('run.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return int(status.isFailure())


if __name__ == '__main__':
    raise SystemExit(main())
