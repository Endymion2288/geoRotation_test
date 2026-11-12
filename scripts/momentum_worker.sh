#!/usr/bin/env bash
set -euo pipefail
georotation_tasks=$1
georotation_index=$2
georotation_library=$3
georotation_task_script=$4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
set +eu
source /eos/home-x/xcheng/FASER/geoRotation_calypso/setup_environment.sh > environment.log 2>&1
georotation_status=$?
set -eu
if [[ $georotation_status -ne 0 ]]; then exit "$georotation_status"; fi
export GEOROTATION_LIBRARY_DIR="$georotation_library"
export LD_LIBRARY_PATH="$GEOROTATION_LIBRARY_DIR:$LD_LIBRARY_PATH"
georotation_python="${GEOROTATION_LIBRARY_DIR%/lib}/python"
if [[ -d "$georotation_python" ]]; then export PYTHONPATH="$georotation_python:$PYTHONPATH"; fi
python3 "$georotation_task_script" --tasks "$georotation_tasks" --index "$georotation_index"
