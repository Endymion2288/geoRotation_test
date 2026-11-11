#!/usr/bin/env bash
set -euo pipefail
georotation_repo=$1
georotation_output=$2
georotation_dsid=$3
georotation_file=$4
georotation_skip=$5
georotation_events=$6
georotation_mode=$7
georotation_shard=$8
mkdir -p "$georotation_output/$georotation_shard"
georotation_destination="$georotation_output/$georotation_shard"
trap 'georotation_exit=$?; printf "%s\n" "$georotation_exit" > "$georotation_destination/exit_code.txt"' EXIT
georotation_scratch=${_CONDOR_SCRATCH_DIR:-$PWD}
cd "$georotation_scratch"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
set +eu
source /eos/home-x/xcheng/FASER/geoRotation_calypso/setup_environment.sh > "$georotation_destination/environment.log" 2>&1
georotation_setup_status=$?
set -eu
if [[ $georotation_setup_status -ne 0 ]]; then
  exit "$georotation_setup_status"
fi
if [[ -f "$georotation_repo/../software.json" ]]; then
  export GEOROTATION_SOFTWARE_MANIFEST="$georotation_repo/../software.json"
  export GEOROTATION_LIBRARY_DIR="$georotation_repo/../lib"
  export LD_LIBRARY_PATH="$GEOROTATION_LIBRARY_DIR:$LD_LIBRARY_PATH"
fi
case "$georotation_mode" in
  current) georotation_material=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage02_closure/material-current-source.json ;;
  historical) georotation_material=/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json ;;
  consistent) georotation_material=/eos/home-x/xcheng/FASER/geoRotation_test_results/stage02_closure/material-current-rotated.json ;;
  *) exit 2 ;;
esac
cp "$georotation_material" material.json
export GEOROTATION_MATERIAL_SOURCE="$georotation_material"
set +e
/usr/bin/time -v -o resource.txt python3 "$georotation_repo/configs/tracking_runner.py" \
  --dsid "$georotation_dsid" --file-index "$georotation_file" --skip "$georotation_skip" \
  --events "$georotation_events" --material "$PWD/material.json" --audit --compact-only \
  > reconstruction.log 2>&1
georotation_status=$?
set -e
for georotation_artifact in stages.jsonl run.json resource.txt reconstruction.log; do
  if [[ -f "$georotation_artifact" ]]; then
    case "$georotation_artifact" in
      *.jsonl|*.log) gzip -c "$georotation_artifact" > "$georotation_destination/$georotation_artifact.gz.partial"
        mv "$georotation_destination/$georotation_artifact.gz.partial" "$georotation_destination/$georotation_artifact.gz" ;;
      *) cp "$georotation_artifact" "$georotation_destination/$georotation_artifact" ;;
    esac
  fi
done
exit "$georotation_status"
