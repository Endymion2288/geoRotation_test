#!/usr/bin/env python3
"""Tracker-only G4 or digitization with an exact truth vertex on the seed plane."""
import argparse
import json
import math
import os
from pathlib import Path
import time


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['simulate','digitize'])
    parser.add_argument('--pdg',type=int,default=-13)
    parser.add_argument('--p',type=float,default=100.,help='Total momentum [GeV]')
    parser.add_argument('--events',type=int,default=1)
    parser.add_argument('--seed',type=int,default=80401)
    parser.add_argument('--input',default='diagnostic.HITS.pool.root')
    parser.add_argument('--output')
    parser.add_argument('--no-sidecar',action='store_true')
    args=parser.parse_args();started=time.monotonic()
    from AthenaCommon.Configurable import Configurable
    from AthenaConfiguration.Enums import BeamType,ProductionStep
    from CalypsoConfiguration.AllConfigFlags import initConfigFlags
    from CalypsoConfiguration.DetectorConfigFlags import setupDetectorsFromList
    from CalypsoConfiguration.MainServicesConfig import MainServicesCfg
    Configurable.configurableRun3Behavior=True
    flags=initConfigFlags()
    flags.Input.isMC=True;flags.Input.ProjectName='mc24'
    flags.IOVDb.DatabaseInstance='OFLP200';flags.IOVDb.GlobalTag='OFLCOND-FASER-05'
    flags.GeoModel.FaserVersion='FASERNU-04'
    flags.Concurrency.NumThreads=0;flags.Concurrency.NumConcurrentEvents=1
    setupDetectorsFromList(flags,['FaserSCT','Dipole'],toggle_geometry=True)
    flags.Exec.MaxEvents=args.events
    if args.mode=='simulate':
        from Campaigns.Utils import Campaign
        if not args.no_sidecar:os.environ['FASER_SCT_LOCAL_TRUTH_FILE']=str(Path('sensor_truth.csv').resolve())
        flags.Common.ProductionStep=ProductionStep.Simulation
        flags.Input.RunNumbers=[1000001];flags.Input.LumiBlockNumbers=[1];flags.Input.OverrideRunNumber=True
        flags.Input.MCCampaign=Campaign.Unknown;flags.Input.TypedCollections=[];flags.Input.MetadataItems=[]
        flags.Beam.Type=BeamType.Collisions
        flags.Sim.PhysicsList='FTFP_BERT';flags.Sim.ReleaseGeoModel=False;flags.Sim.IncludeParentsInG4Event=True
        mass=105.6583755 if abs(args.pdg)==13 else .51099895
        energy=math.hypot(args.p*1000,mass)
        gun=dict(Generator='SingleParticle',pid=args.pdg,energy=energy,mass=mass,x=20.,y=-20.,z=11.415,
                 theta=.002,phi=0.,t=0.,axialTiming=False,randomSeed=args.seed)
        flags.addFlag('Sim.Gun',gun)
        flags.addFlag('Input.InitialTimeStamp',0)
        flags.Output.HITSFileName=args.output or 'diagnostic.HITS.pool.root'
        flags.lock();acc=MainServicesCfg(flags)
        from FaserParticleGun.FaserParticleGunConfig import FaserParticleGunCfg
        from McEventSelector.McEventSelectorConfig import McEventSelectorCfg
        from xAODEventInfoCnv.xAODEventInfoCnvConfig import EventInfoCnvAlgCfg
        from G4FaserAlg.G4FaserAlgConfig import G4FaserAlgCfg,getG4FaserAlgItemList
        from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
        from xAODMetaDataCnv.InfileMetaDataConfig import SetupMetaDataForStreamCfg
        acc.merge(FaserParticleGunCfg(flags),sequenceName='AthBeginSeq');acc.merge(McEventSelectorCfg(flags))
        acc.merge(EventInfoCnvAlgCfg(flags,disableBeamSpot=True));acc.merge(G4FaserAlgCfg(flags))
        acc.merge(OutputStreamCfg(flags,'HITS',ItemList=getG4FaserAlgItemList(flags),disableEventTag=True,AcceptAlgs=['G4FaserAlg']))
        acc.merge(SetupMetaDataForStreamCfg(flags,'HITS',AcceptAlgs=['G4FaserAlg']))
    else:
        flags.Input.Files=[str(Path(args.input).resolve())]
        flags.Output.RDOFileName=args.output or 'diagnostic.RDO.pool.root'
        flags.Digitization.EnableTruth=True;flags.Beam.NumberOfCollisions=0.
        flags.lock();acc=MainServicesCfg(flags)
        from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
        from FaserGeoModel.FaserGeoModelConfig import FaserGeometryCfg
        from FaserSCT_Digitization.FaserSCT_DigitizationConfigNew import FaserSCT_DigitizationCfg
        from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
        acc.merge(PoolReadCfg(flags));acc.merge(FaserGeometryCfg(flags));acc.merge(FaserSCT_DigitizationCfg(flags))
        acc.merge(OutputStreamCfg(flags,'RDO',ItemList=['McEventCollection#*','FaserSiHitCollection#SCT_Hits']))
    acc.getService('MessageSvc').defaultLimit=20
    sc=acc.run(maxEvents=args.events)
    Path(args.mode+'_run.json').write_text(json.dumps(dict(args=vars(args),success=sc.isSuccess(),
        wall_seconds=time.monotonic()-started,geometry='FASERNU-04',conditions='OFLCOND-FASER-05',
        simulated_detectors=['FaserSCT','Dipole'],truth_plane_z_F_mm=11.415),indent=2)+'\n')
    return int(sc.isFailure())


if __name__=='__main__':raise SystemExit(main())
