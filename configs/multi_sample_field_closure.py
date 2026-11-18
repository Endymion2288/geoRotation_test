#!/usr/bin/env python3
"""Run multi-sample deterministic full-field extrapolation + native CircleFit closure."""
import argparse
import ctypes
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    task = json.loads(args.task.read_text())

    import AthenaPython.PyAthena as PyAthena
    from AthenaCommon.Configurable import Configurable
    from AthenaConfiguration.ComponentFactory import CompFactory
    from CalypsoConfiguration.AllConfigFlags import initConfigFlags
    from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    from FaserSCT_GeoModel.FaserSCT_GeoModelConfig import FaserSCT_GeometryCfg
    from FaserGeoModel.FaserGeoModelConfig import FaserGeometryCfg
    from MagFieldServices.MagFieldServicesConfig import MagneticFieldSvcCfg
    from FaserActsGeometry.ActsGeometryConfig import ActsTrackingGeometrySvcCfg, ActsTrackingGeometryToolCfg

    class Observer(PyAthena.Alg):
        def initialize(self):
            self.library = ctypes.CDLL(task['helper_library'])
            self.library.multi_sample_seed_field_create.restype = ctypes.c_void_p
            self.library.multi_sample_seed_field_execute.argtypes = [ctypes.c_void_p]
            self.library.multi_sample_seed_field_execute.restype = ctypes.c_char_p
            self.library.multi_sample_seed_field_delete.argtypes = [ctypes.c_void_p]
            self.closure = self.library.multi_sample_seed_field_create()
            self.entry = 0
            self.particles = []
            return PyAthena.StatusCode.Success

        def execute(self):
            recs = json.loads(self.library.multi_sample_seed_field_execute(self.closure))
            for r in recs:
                r['entry'] = self.entry + int(task.get('skip', 0))
                self.particles.append(r)
            self.entry += 1
            return PyAthena.StatusCode.Success

        def finalize(self):
            self.library.multi_sample_seed_field_delete(self.closure)
            assert self.entry == task['events']
            out_dir = Path(task['output'])
            out_dir.mkdir(parents=True, exist_ok=True)
            payload = dict(
                sample=task['sample'],
                events=self.entry,
                particle_count=len(self.particles),
                particles=self.particles,
            )
            (out_dir / 'field_closure.json').write_text(json.dumps(payload, indent=2) + '\n')
            return PyAthena.StatusCode.Success

    Configurable.configurableRun3Behavior = True
    flags = initConfigFlags()
    flags.Input.Files = [task['input']]
    flags.Input.isMC = True
    flags.Input.ProjectName = 'mc24'
    flags.IOVDb.DatabaseInstance = 'OFLP200'
    flags.GeoModel.FaserVersion = 'FASERNU-04'
    flags.IOVDb.GlobalTag = 'OFLCOND-FASER-05'
    flags.TrackingGeometry.MaterialSource = task['material']
    flags.Exec.SkipEvents = int(task.get('skip', 0))
    flags.lock()

    acc = MainServicesCfg(flags)
    acc.merge(PoolReadCfg(flags))
    acc.merge(FaserGeometryCfg(flags))
    acc.merge(FaserSCT_GeometryCfg(flags))
    acc.merge(MagneticFieldSvcCfg(flags))
    acc.merge(ActsTrackingGeometrySvcCfg(flags))
    geo, tool = ActsTrackingGeometryToolCfg(flags)
    acc.merge(geo)
    extrapolator = CompFactory.FaserActsExtrapolationTool(
        'MultiSampleSeedFieldExtrapolator',
        TrackingGeometryTool=tool,
        MaxSteps=1000,
        MaxStepSize=.01,
        InteractionMultiScatering=False,
        InteractionEloss=False,
    )
    acc.addPublicTool(extrapolator)
    acc.addEventAlgo(Observer('MultiSampleFieldClosureAlg'))
    return 0 if acc.run(maxEvents=task['events']).isSuccess() else 1


if __name__ == '__main__':
    raise SystemExit(main())
