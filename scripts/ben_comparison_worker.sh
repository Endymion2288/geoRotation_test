#!/usr/bin/env bash
set -euo pipefail

TASKS_JSON=$1
TASK_INDEX=$2

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1

CALYPSO_LATEST="/eos/home-x/xcheng/FASER/geoRotation_calypso_latest"
CALYPSO_FIX_V1="/eos/home-x/xcheng/FASER/geoRotation_calypso_latest_fix_v1"
FIX_TRK="/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_latest_validation/20261118_fix_v1/software_fix_v1"
FIX_NTP="/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_waveform_reduced_charge/software_ntuple_fix_v1_clean"

WORK_ROOT="${_CONDOR_SCRATCH_DIR:-$PWD}/work_${TASK_INDEX}_$$"
mkdir -p "$WORK_ROOT/reco" "$WORK_ROOT/ntp"

set +eu
source "${CALYPSO_LATEST}/setup_environment.sh" > "$WORK_ROOT/environment.log" 2>&1
ENV_RC=$?
set -eu
if [[ $ENV_RC -ne 0 ]]; then
    echo "Environment setup failed rc=$ENV_RC" >&2
    cat "$WORK_ROOT/environment.log" >&2
    exit "$ENV_RC"
fi

export LD_LIBRARY_PATH="${FIX_TRK}/lib:${FIX_NTP}/lib:${LD_LIBRARY_PATH}"
export PYTHONPATH="${FIX_TRK}/python:${FIX_NTP}/python:${PYTHONPATH}"
export JOBOPTSEARCHPATH="${FIX_TRK}/lib:${FIX_NTP}/lib:${JOBOPTSEARCHPATH:-}"

# Parse task fields from frozen JSON
read -r TASK_ID PERIOD RUN SEG RAW_URL COND PRESET OUT_DIR < <(
python3 - "$TASKS_JSON" "$TASK_INDEX" <<'PY'
import json, sys
tasks = json.loads(open(sys.argv[1]).read())
t = tasks[int(sys.argv[2])]
print(t["task_id"], t["period"], t["run"], t["segment"], t["raw_url"], t["cond_tag"], t["timing_preset"], t["output_dir"])
PY
)

mkdir -p "$OUT_DIR"
FINAL_PHYS="${OUT_DIR}/${TASK_ID}_phys.root"
SUMMARY_JSON="${OUT_DIR}/summary.json"
if [[ -s "$FINAL_PHYS" && -s "$SUMMARY_JSON" ]]; then
    echo "Shard ${TASK_ID} already completed; skipping."
    exit 0
fi

RAW_BASENAME="$(basename "$RAW_URL" .raw)"
XAOD_LOCAL="${WORK_ROOT}/reco/${RAW_BASENAME}-xAOD.root"
PHYS_LOCAL="${WORK_ROOT}/ntp/${TASK_ID}_phys.root"

echo "[$(date -u +%FT%TZ)] Starting ${TASK_ID} on $(hostname): raw=${RAW_URL} cond=${COND} preset=${PRESET}"

T0=$(date +%s)
set +e
(
    cd "$WORK_ROOT/reco"
    python3 "${CALYPSO_LATEST}/Control/CalypsoExample/Reconstruction/scripts/faser_reco.py" \
        "$RAW_URL" \
        --nevents -1 \
        --geom TI12Data04 \
        --noIFT \
        --noBackward \
        --cond "$COND" \
        > "$WORK_ROOT/reco.log" 2>&1
)
RC_RECO=$?
T1=$(date +%s)
RECO_WALL=$((T1 - T0))

RC_NTP=99
NTP_WALL=0
XAOD_BYTES=0
if [[ $RC_RECO -eq 0 && -f "$XAOD_LOCAL" ]]; then
    XAOD_BYTES=$(stat -c%s "$XAOD_LOCAL")
    T2=$(date +%s)
    (
        cd "$WORK_ROOT/ntp"
        python3 "${CALYPSO_FIX_V1}/PhysicsAnalysis/NtupleDumper/scripts/faser_ntuple_maker.py" \
            "$XAOD_LOCAL" \
            --outfile "$PHYS_LOCAL" \
            --geom TI12Data04 \
            --cond "$COND" \
            --timingPreset "$PRESET" \
            --no_stable \
            > "$WORK_ROOT/ntp.log" 2>&1
    )
    RC_NTP=$?
    T3=$(date +%s)
    NTP_WALL=$((T3 - T2))
fi
set -e

# Always clean up large intermediate xAOD immediately
rm -f "$XAOD_LOCAL"

if [[ $RC_RECO -ne 0 || $RC_NTP -ne 0 || ! -s "$PHYS_LOCAL" ]]; then
    echo "FAILED ${TASK_ID}: rc_reco=${RC_RECO} (${RECO_WALL}s), rc_ntp=${RC_NTP} (${NTP_WALL}s)" >&2
    cp -f "$WORK_ROOT/reco.log" "$OUT_DIR/reco_failed.log" 2>/dev/null || true
    cp -f "$WORK_ROOT/ntp.log" "$OUT_DIR/ntp_failed.log" 2>/dev/null || true
    rm -rf "$WORK_ROOT"
    exit 1
fi

# Verify produced PHYS tree & capture compact metadata before staging to EOS
PHYS_STATS=$(python3 - "$PHYS_LOCAL" <<'PY'
import json, sys, ROOT
f = ROOT.TFile.Open(sys.argv[1])
t = f.Get("nt")
n = int(t.GetEntries()) if t else -1
e0 = e1 = -1
if n > 0:
    t.GetEntry(0); e0 = int(t.eventID)
    t.GetEntry(n - 1); e1 = int(t.eventID)
f.Close()
print(json.dumps({"entries": n, "eventID_min": e0, "eventID_max": e1}))
PY
)

PHYS_BYTES=$(stat -c%s "$PHYS_LOCAL")
cp -f "$PHYS_LOCAL" "${FINAL_PHYS}.tmp"
mv -f "${FINAL_PHYS}.tmp" "$FINAL_PHYS"

# Write compact summary manifest
python3 - "$SUMMARY_JSON" "$TASK_ID" "$PERIOD" "$RUN" "$SEG" "$RAW_URL" "$COND" "$PRESET" \
          "$RECO_WALL" "$NTP_WALL" "$XAOD_BYTES" "$PHYS_BYTES" "$PHYS_STATS" "$WORK_ROOT/reco.log" "$WORK_ROOT/ntp.log" <<'PY'
import datetime, json, re, socket, sys
from pathlib import Path
summary_path, task_id, period, run, seg, raw_url, cond, preset = sys.argv[1:9]
reco_wall, ntp_wall, xaod_bytes, phys_bytes = map(int, sys.argv[9:13])
pstats = json.loads(sys.argv[13])
reco_txt = Path(sys.argv[14]).read_text(errors="replace")
ntp_txt = Path(sys.argv[15]).read_text(errors="replace")
m_ev = re.search(r"CKF_woIFT\s+INFO\s+(\d+)\s+events processed", reco_txt)
m_al = re.search(r"resolved to (\S+) for folder /Tracker/Align", reco_txt)
m_rc = re.search(r"Reduced VetoNu charge configuration.*", ntp_txt)
summary = {
    "task_id": task_id,
    "period": period,
    "run": int(run),
    "segment": seg,
    "raw_url": raw_url,
    "cond_tag": cond,
    "align_tag": m_al.group(1) if m_al else None,
    "timing_preset": preset,
    "reduced_charge_config_line": m_rc.group(0).strip() if m_rc else None,
    "raw_events_processed": int(m_ev.group(1)) if m_ev else None,
    "phys_entries": pstats["entries"],
    "eventID_min": pstats["eventID_min"],
    "eventID_max": pstats["eventID_max"],
    "reco_wall_s": reco_wall,
    "ntp_wall_s": ntp_wall,
    "scratch_xaod_bytes": xaod_bytes,
    "phys_bytes": phys_bytes,
    "hostname": socket.gethostname(),
    "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "calypso_latest_commit": "7fd40df8cd3b60e81b4095a4f84310c4a54700a0",
    "calypso_fix_v1_tracking_commit": "e47b85fd",
    "calypso_fix_v1_ntuple_commit": "1c0f3f45",
}
Path(summary_path).write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary))
PY

rm -rf "$WORK_ROOT"
echo "[$(date -u +%FT%TZ)] Finished ${TASK_ID} cleanly."
