#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="/eos/home-x/xcheng/FASER/geoRotation_test"
OUT_BIN_DIR="${1:-/eos/home-x/xcheng/FASER/geoRotation_test_results/stage08_overlay_2track/run_20261125_v1/bin}"

mkdir -p "$OUT_BIN_DIR"

CXX="${CXX:-g++}"
ROOT_CFLAGS="$(root-config --cflags)"
ROOT_LIBS="$(root-config --libs)"

echo "[build_overlay_analysis] Compiler: $($CXX --version | head -n 1)"
echo "[build_overlay_analysis] ROOT: $(root-config --version)"
echo "[build_overlay_analysis] Output dir: $OUT_BIN_DIR"

set -x
$CXX -O3 -std=c++17 -Wall -Wextra -Wno-unused-parameter -fopenmp \
    "$REPO_ROOT/analysis/overlay_signal_110038.cpp" \
    $ROOT_CFLAGS $ROOT_LIBS -lm \
    -o "$OUT_BIN_DIR/overlay_signal_110038"

$CXX -O3 -std=c++17 -Wall -Wextra -Wno-unused-parameter -fopenmp \
    "$REPO_ROOT/analysis/overlay_analyze.cpp" \
    $ROOT_CFLAGS $ROOT_LIBS -lm \
    -o "$OUT_BIN_DIR/overlay_analyze"

$CXX -O3 -std=c++17 -Wall -Wextra -Wno-unused-parameter -fopenmp \
    "$REPO_ROOT/analysis/overlay_tier2_analyze.cpp" \
    $ROOT_CFLAGS $ROOT_LIBS -lm \
    -o "$OUT_BIN_DIR/overlay_tier2_analyze"
set +x

ls -lh "$OUT_BIN_DIR/overlay_signal_110038" "$OUT_BIN_DIR/overlay_analyze" "$OUT_BIN_DIR/overlay_tier2_analyze"
echo "[build_overlay_analysis] SUCCESS"
