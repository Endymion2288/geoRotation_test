#!/usr/bin/env python3
"""Compare native Calypso circle seeding with its full-field extrapolation tool."""
import argparse
import json
from pathlib import Path




def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();task=json.loads(args.task.read_text())
    import ctypes
    import AthenaPython.PyAthena as PyAthena
    from AthenaCommon.Configurable import Configurable
    from AthenaConfiguration.ComponentFactory import CompFactory
    from CalypsoConfiguration.AllConfigFlags import initConfigFlags
    from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    from FaserSCT_GeoModel.FaserSCT_GeoModelConfig import FaserSCT_GeometryCfg
    from FaserGeoModel.FaserGeoModelConfig import FaserGeometryCfg
    from MagFieldServices.MagFieldServicesConfig import MagneticFieldSvcCfg
    from FaserActsGeometry.ActsGeometryConfig import ActsTrackingGeometrySvcCfg,ActsTrackingGeometryToolCfg
    class Observer(PyAthena.Alg):
        def initialize(self):
            self.library=ctypes.CDLL(task['helper_library'])
            self.library.latest_seed_field_create.restype=ctypes.c_void_p
            self.library.latest_seed_field_execute.argtypes=[ctypes.c_void_p]
            self.library.latest_seed_field_execute.restype=ctypes.c_char_p
            self.library.latest_seed_field_delete.argtypes=[ctypes.c_void_p]
            self.closure=self.library.latest_seed_field_create()
            self.entry=0;self.rows=[]
            return PyAthena.StatusCode.Success
        def execute(self):
            self.rows.append(dict(entry=self.entry,particles=json.loads(self.library.latest_seed_field_execute(self.closure))))
            self.entry+=1
            return PyAthena.StatusCode.Success
        def finalize(self):
            self.library.latest_seed_field_delete(self.closure)
            assert self.entry==task['events']
            output=dict(source_task=task,events=self.entry,rows=self.rows,
                loaded_libraries=sorted({line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines() if any(name in line for name in ['libFaserActs','libLatestSeedField','libActsCore']) and '.so' in line}),
                definition='Native Calypso extrapolation in actual geometry/full field, then native CircleFit with original +/-30 mm points and 0.55 T formula; no detector measurements or stochastic G4 propagation')
            (Path(task['output'])/'closure.json').write_text(json.dumps(output,indent=2)+'\n')
            return PyAthena.StatusCode.Success
    Configurable.configurableRun3Behavior=True
    flags=initConfigFlags();flags.Input.Files=[task['input']]
    flags.Input.isMC=True;flags.Input.ProjectName='mc24';flags.IOVDb.DatabaseInstance='OFLP200'
    flags.GeoModel.FaserVersion='FASERNU-04';flags.IOVDb.GlobalTag='OFLCOND-FASER-05'
    flags.TrackingGeometry.MaterialSource=task['material'];flags.lock()
    acc=MainServicesCfg(flags);acc.merge(PoolReadCfg(flags));acc.merge(FaserGeometryCfg(flags))
    acc.merge(FaserSCT_GeometryCfg(flags))
    acc.merge(MagneticFieldSvcCfg(flags));acc.merge(ActsTrackingGeometrySvcCfg(flags))
    geo,tool=ActsTrackingGeometryToolCfg(flags);acc.merge(geo)
    extrapolator=CompFactory.FaserActsExtrapolationTool('MomentumSeedFieldExtrapolator',
        TrackingGeometryTool=tool,MaxSteps=1000,MaxStepSize=.01,
        InteractionMultiScatering=False,InteractionEloss=False)
    acc.addPublicTool(extrapolator);acc.addEventAlgo(Observer('SeedFieldClosure'))
    return 0 if acc.run(maxEvents=task['events']).isSuccess() else 1


if __name__=='__main__':
    raise SystemExit(main())
