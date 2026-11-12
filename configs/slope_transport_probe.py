#!/usr/bin/env python3
"""Propagate archived Calypso states with the native field and tracking geometry."""
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--material', required=True)
    parser.add_argument('--variant', choices=['beforeGeo', 'afterGeo'], required=True)
    parser.add_argument('--states', default='')
    parser.add_argument('--direction-scan', action='store_true')
    parser.add_argument('--skip', type=int, default=0)
    parser.add_argument('--events', type=int, default=1)
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
    propagator = CompFactory.FaserActsExtrapolationTool(
        'SlopeTransportExtrapolator', TrackingGeometryTool=geometry, MaxSteps=10000)
    acc.addEventAlgo(CompFactory.CalypsoSlopeTransportProbeAlg(
        'SlopeTransportProbe', ExtrapolationTool=propagator,
        Rotated=args.variant != 'beforeGeo', DirectionScan=args.direction_scan,
        InputStates=args.states, Skip=args.skip))
    acc.getService('MessageSvc').defaultLimit = 20
    return int(acc.run(args.events).isFailure())


if __name__ == '__main__':
    raise SystemExit(main())
