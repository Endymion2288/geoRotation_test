#!/usr/bin/env bash
# Generate Gaudi metadata from the actual newly built plugin, then install it.
set -euo pipefail
georotation_source=/eos/home-x/xcheng/FASER/geoRotation_calypso
georotation_build="$georotation_source/build"
georotation_platform=x86_64-el9-gcc13-opt
georotation_gen="$georotation_build/Tracking/Acts/FaserActsKalmanFilter/genConf"
georotation_external=/cvmfs/atlas.cern.ch/repo/sw/software/24.0/AthenaExternals/24.0.41/InstallArea/$georotation_platform
env LD_LIBRARY_PATH="$georotation_build/$georotation_platform/lib:$LD_LIBRARY_PATH" \
    "$georotation_external/bin/genconf" -o "$georotation_gen" -p FaserActsKalmanFilter --no-init \
    -i "$georotation_build/$georotation_platform/lib/libFaserActsKalmanFilter.so"
cp "$georotation_gen/FaserActsKalmanFilterConf.py" "$georotation_build/$georotation_platform/python/FaserActsKalmanFilter/"
cp "$georotation_gen/FaserActsKalmanFilter.confdb" "$georotation_build/$georotation_platform/lib/libFaserActsKalmanFilter.confdb"
touch "$georotation_gen/FaserActsKalmanFilter.confdb2_part"
python3 /cvmfs/atlas.cern.ch/repo/sw/software/24.0/Athena/24.0.41/InstallArea/$georotation_platform/cmake/modules/scripts/mergeConfdb2.py \
    "$georotation_build/$georotation_platform/lib/Calypso.confdb2" "$georotation_build/CMakeFiles/CalypsoConfdb2MergeFiles.txt"
cmake -DCMAKE_INSTALL_PREFIX="$georotation_source/run" -P "$georotation_build/Tracking/Acts/FaserActsKalmanFilter/cmake_install.cmake"
cp -p "$georotation_build/$georotation_platform/lib/Calypso.confdb2" "$georotation_source/run/lib/"
python3 - <<'PY'
from AthenaConfiguration.ComponentFactory import CompFactory
ckf=CompFactory.CKF2('MetadataCheck')
for name in ['AuditSelectedStatesOnly','MaterialEnergyLoss','FixedTruthHits','SeedOverrideFile','IsolateCandidateSmoothing']:
    print(name,getattr(ckf,name))
PY
