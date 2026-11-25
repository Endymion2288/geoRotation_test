#!/usr/bin/env bash
set -euo pipefail

TASKS_JSON=$1
TASK_INDEX=$2

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1

WORK_ROOT="${_CONDOR_SCRATCH_DIR:-$PWD}/work_overlay_${TASK_INDEX}_$$"
mkdir -p "$WORK_ROOT/env_setup" "$WORK_ROOT/reco" "$WORK_ROOT/conv"
cleanup() {
    rm -rf "$WORK_ROOT"
}
trap cleanup EXIT

# Parse frozen task metadata using system python3 before Athena setup
read -r TASK_ID DATASET POS_BIN NEG_BIN IS_MC INPUT_RDO OUT_DIR COND_TAG GEOM GEOM_ROT \
        MATERIAL_MAP CALYPSO_SRC CALYPSO_RUN CALYPSO_COMMIT RECO_SCRIPT CONV_BIN NEVENTS SEED_COV_PHI < <(
python3 - "$TASKS_JSON" "$TASK_INDEX" <<'PY'
import json, sys
tasks = json.loads(open(sys.argv[1]).read())
t = tasks[int(sys.argv[2])]
scp = t.get("seed_cov_phi")
fields = [
    t["task_id"],
    t["dataset"],
    str(t["pos_bin"]),
    str(t["neg_bin"]),
    "1" if t["is_mc"] else "0",
    t["input_rdo"],
    t["output_dir"],
    t["cond_tag"],
    t["geom"],
    "1" if t.get("geometry_rotation", True) else "0",
    t["material_map"],
    t["calypso_source"],
    t["calypso_run"],
    t["calypso_commit"],
    t["overlay_reco_script"],
    t["acts_to_faser_bin"],
    str(t.get("nevents", -1)),
    "NONE" if scp is None else repr(float(scp)),
]
print("\t".join(fields))
PY
)

mkdir -p "$OUT_DIR"
FINAL_NTP="${OUT_DIR}/${TASK_ID}_nt.root"
SUMMARY_JSON="${OUT_DIR}/summary.json"

# Idempotent fast-exit if shard already completed cleanly
if [[ -s "$FINAL_NTP" && -s "$SUMMARY_JSON" ]]; then
    ALREADY_OK=$(python3 - "$SUMMARY_JSON" <<'PY'
import json, sys
try:
    d = json.loads(open(sys.argv[1]).read())
    print("yes" if d.get("status") == "ok" else "no")
except Exception:
    print("no")
PY
)
    if [[ "$ALREADY_OK" == "yes" ]]; then
        echo "[$(date -u +%FT%TZ)] Shard ${TASK_ID} already completed (status=ok); skipping."
        exit 0
    fi
fi

echo "[$(date -u +%FT%TZ)] Starting ${TASK_ID} on $(hostname -f): rdo=${INPUT_RDO} cond=${COND_TAG} isMC=${IS_MC} rot=${GEOM_ROT}"

# Stage input OverlaidRDO locally to avoid mid-run EOS stalls
LOCAL_RDO="${WORK_ROOT}/input_rdo.pool.root"
if ! cp -f "$INPUT_RDO" "$LOCAL_RDO"; then
    echo "Direct cp failed for ${INPUT_RDO}; retrying via xrdcp..." >&2
    xrdcp -f "root://eospublic.cern.ch/${INPUT_RDO}" "$LOCAL_RDO"
fi
if [[ ! -s "$LOCAL_RDO" ]]; then
    echo "ERROR: Local staged RDO is empty: $LOCAL_RDO" >&2
    exit 2
fi

# Source Athena 24.0.41 + Calypso run/ inside private scratch cwd (isolates .asetup.save)
set +eu
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
(
    true
)
cd "$WORK_ROOT/env_setup"
source "${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" --quiet > "$WORK_ROOT/env.log" 2>&1
ENV_RC=$?
if [[ $ENV_RC -eq 0 ]]; then
    source "${AtlasSetup}/scripts/asetup.sh" --input="${CALYPSO_SRC}/asetup.faser" Athena,24.0.41 >> "$WORK_ROOT/env.log" 2>&1
    ENV_RC=$?
fi
if [[ $ENV_RC -eq 0 ]]; then
    source "${CALYPSO_RUN}/setup.sh" >> "$WORK_ROOT/env.log" 2>&1
    ENV_RC=$?
fi
set -eu

if [[ $ENV_RC -ne 0 ]]; then
    echo "ERROR: Environment setup failed rc=${ENV_RC}" >&2
    cp -f "$WORK_ROOT/env.log" "${OUT_DIR}/failed_env.log" 2>/dev/null || true
    exit "$ENV_RC"
fi

if [[ "$GEOM_ROT" == "1" ]]; then
    export FASER_ACTS_MATERIAL_MAP="$MATERIAL_MAP"
    ROT_ARG="true"
else
    ROT_ARG="false"
fi

RAW_NTP="${WORK_ROOT}/reco/nt_raw.root"
FASER_NTP="${WORK_ROOT}/conv/nt_faser.root"
CONV_STATS="${WORK_ROOT}/conv/conv_stats.json"

RECO_CMD=(
    python3 "$RECO_SCRIPT"
    "$LOCAL_RDO"
    --output "$RAW_NTP"
    --geom "$GEOM"
    --cond "$COND_TAG"
    --materialMap "$MATERIAL_MAP"
    --geometry-rotation "$ROT_ARG"
    --nevents "$NEVENTS"
)
if [[ "$IS_MC" == "1" ]]; then
    RECO_CMD+=(--isMC)
fi
if [[ "$SEED_COV_PHI" != "NONE" ]]; then
    RECO_CMD+=(--seed-cov-phi "$SEED_COV_PHI")
fi

T0=$(date +%s)
set +e
(
    cd "$WORK_ROOT/reco"
    /usr/bin/time -v -o "$WORK_ROOT/reco.time" "${RECO_CMD[@]}" > "$WORK_ROOT/reco.log" 2>&1
)
RC_RECO=$?
T1=$(date +%s)
RECO_WALL=$((T1 - T0))

RC_CONV=99
CONV_WALL=0
if [[ $RC_RECO -eq 0 && -s "$RAW_NTP" ]]; then
    T2=$(date +%s)
    CONV_CMD=(
        "$CONV_BIN"
        --input "$RAW_NTP"
        --output "$FASER_NTP"
        --rdo "$LOCAL_RDO"
        --rotate "$GEOM_ROT"
        --stats-json "$CONV_STATS"
        --meta-frame "FASER"
        --meta-calypso-commit "$CALYPSO_COMMIT"
        --meta-cond-tag "$COND_TAG"
        --meta-material-map "$MATERIAL_MAP"
        --meta-source-rdo "$INPUT_RDO"
        --meta-dataset "$DATASET"
        --meta-pos-bin "$POS_BIN"
        --meta-neg-bin "$NEG_BIN"
    )
    if [[ "$NEVENTS" != "-1" && "$NEVENTS" -gt 0 ]]; then
        CONV_CMD+=(--expected-entries "$NEVENTS")
    fi
    (
        cd "$WORK_ROOT/conv"
        /usr/bin/time -v -o "$WORK_ROOT/conv.time" "${CONV_CMD[@]}" > "$WORK_ROOT/conv.log" 2>&1
    )
    RC_CONV=$?
    T3=$(date +%s)
    CONV_WALL=$((T3 - T2))
fi
set -e

if [[ $RC_RECO -ne 0 || $RC_CONV -ne 0 || ! -s "$FASER_NTP" ]]; then
    echo "FAILED ${TASK_ID}: rc_reco=${RC_RECO} (${RECO_WALL}s), rc_conv=${RC_CONV} (${CONV_WALL}s)" >&2
    cp -f "$WORK_ROOT/reco.log" "${OUT_DIR}/failed_reco.log" 2>/dev/null || true
    cp -f "$WORK_ROOT/reco.time" "${OUT_DIR}/failed_reco.time" 2>/dev/null || true
    cp -f "$WORK_ROOT/conv.log" "${OUT_DIR}/failed_conv.log" 2>/dev/null || true
    cp -f "$WORK_ROOT/conv.time" "${OUT_DIR}/failed_conv.time" 2>/dev/null || true
    cp -f "$CONV_STATS" "${OUT_DIR}/failed_conv_stats.json" 2>/dev/null || true
    exit 1
fi

# Atomic staging of compact FASER-frame ntuple to EOS (zero xAOD written)
NTP_BYTES=$(stat -c%s "$FASER_NTP")
cp -f "$FASER_NTP" "${FINAL_NTP}.tmp"
mv -f "${FINAL_NTP}.tmp" "$FINAL_NTP"

# Write compact summary.json atomically
python3 - "$SUMMARY_JSON" "$TASK_ID" "$DATASET" "$POS_BIN" "$NEG_BIN" "$IS_MC" "$INPUT_RDO" \
          "$FINAL_NTP" "$COND_TAG" "$GEOM" "$GEOM_ROT" "$MATERIAL_MAP" "$CALYPSO_COMMIT" \
          "$RECO_WALL" "$CONV_WALL" "$NTP_BYTES" "$WORK_ROOT/reco.log" "$WORK_ROOT/reco.time" \
          "$CONV_STATS" "$SEED_COV_PHI" <<'PY'
import datetime, json, os, re, socket, sys
from pathlib import Path

(
    summary_path, task_id, dataset, pos_bin, neg_bin, is_mc, input_rdo,
    final_ntp, cond_tag, geom, geom_rot, material_map, calypso_commit,
    reco_wall, conv_wall, ntp_bytes, reco_log_p, reco_time_p, conv_stats_p,
    seed_cov_phi_s,
) = sys.argv[1:21]

reco_txt = Path(reco_log_p).read_text(errors="replace")
time_txt = Path(reco_time_p).read_text(errors="replace")
conv_stats = json.loads(Path(conv_stats_p).read_text())

m_al = re.search(r"resolved to (\S+) for folder /Tracker/Align", reco_txt)
m_mat = re.search(r"Configured to use material input:\s*(\S+)", reco_txt)
m_gt = re.search(r"Global tag:\s*(\S+)", reco_txt)
m_ev = re.search(r"CKF_woIFT\s+INFO\s+(\d+)\s+events processed", reco_txt)
m_rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", time_txt)
m_athena = re.search(r"Finish execution in ([0-9.]+) seconds", reco_txt)
m_cov = re.search(r"CKF_woIFT TrackSeed covariance:.*", reco_txt)

summary = {
    "status": "ok",
    "task_id": task_id,
    "dataset": dataset,
    "pos_bin": int(pos_bin),
    "neg_bin": int(neg_bin),
    "is_mc": is_mc == "1",
    "input_rdo": input_rdo,
    "output_ntuple": final_ntp,
    "frame": "FASER",
    "geometry_rotation": geom_rot == "1",
    "seed_cov_phi": None if seed_cov_phi_s == "NONE" else float(seed_cov_phi_s),
    "seed_cov_log_line": m_cov.group(0).strip() if m_cov else None,
    "geom": geom,
    "cond_tag_requested": cond_tag,
    "cond_tag_resolved": m_gt.group(1) if m_gt else None,
    "align_tag": m_al.group(1) if m_al else None,
    "material_map_requested": material_map,
    "material_map_loaded": m_mat.group(1) if m_mat else None,
    "calypso_commit": calypso_commit,
    "ckf_events_processed": int(m_ev.group(1)) if m_ev else None,
    "reco_wall_s": int(reco_wall),
    "athena_wall_s": float(m_athena.group(1)) if m_athena else None,
    "conv_wall_s": int(conv_wall),
    "peak_rss_kb": int(m_rss.group(1)) if m_rss else None,
    "ntuple_bytes": int(ntp_bytes),
    "conv_stats": conv_stats,
    "hostname": socket.getfqdn(),
    "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
}
tmp_path = Path(summary_path + ".tmp")
tmp_path.write_text(json.dumps(summary, indent=2) + "\n")
os.replace(tmp_path, summary_path)
print(json.dumps({"status": "ok", "task_id": task_id, "entries": conv_stats["nt_entries"], "reco_wall_s": int(reco_wall), "peak_rss_kb": summary["peak_rss_kb"]}))
PY

echo "[$(date -u +%FT%TZ)] Finished ${TASK_ID} cleanly."
