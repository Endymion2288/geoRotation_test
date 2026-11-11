#!/usr/bin/env python3
"""One-event G4 validation of the optional SCT local-momentum sidecar."""
import os
from pathlib import Path
import sys

if "--no-sidecar" in sys.argv:
    os.environ.pop("FASER_SCT_LOCAL_TRUTH_FILE", None)
else:
    os.environ["FASER_SCT_LOCAL_TRUTH_FILE"] = str(Path("sensor_truth.csv").resolve())
script = Path("/eos/home-x/xcheng/FASER/geoRotation_calypso/Simulation/G4Faser/G4FaserAlg/test/G4FaserAlgConfig_Test.py")
gun = {"Generator": "SingleParticle", "pid": -13, "energy": 100000., "x": 0., "y": 0.,
       "z": -10., "theta": 0.002, "phi": 0., "randomSeed": 80217}
sys.argv = [str(script), "--evtMax=1", "Output.HITSFileName=local_truth.HITS.pool.root", f"Sim.Gun={gun!r}"]
os.execv(sys.executable, [sys.executable, *sys.argv])
