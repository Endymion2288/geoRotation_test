#!/usr/bin/env bash
set -euo pipefail

TASKS_JSON=$1
TASK_INDEX=$2

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1

WORK_ROOT="${_CONDOR_SCRATCH_DIR:-$PWD}/work_sig110038_${TASK_INDEX}_$$"
mkdir -p "$WORK_ROOT/env_setup" "$WORK_ROOT/reco" "$WORK_ROOT/ntp"
cleanup() {
    rm -rf "$WORK_ROOT"
}
trap cleanup EXIT

read -r TASK_ID INPUT_RDO OUT_DIR COND_TAG GEOM MATERIAL_MAP CALYPSO_SRC CALYPSO_RUN CALYPSO_COMMIT NEVENTS < <(
python3 - "$TASKS_JSON" "$TASK_INDEX" <<'PY'
import json, sys
tasks = json.loads(open(sys.argv[1]).read())
t = tasks[int(sys.argv[2])]
fields = [
    t["task_id"],
    t["input_rdo"],
    t["output_dir"],
    t["cond_tag"],
    t["geom"],
    t["material_map"],
    t["calypso_source"],
    t["calypso_run"],
    t["calypso_commit"],
    str(t.get("nevents", -1)),
]
print("\t".join(fields))
PY
)

mkdir -p "$OUT_DIR"
FINAL_PHYS="${OUT_DIR}/${TASK_ID}_phys.root"
SUMMARY_JSON="${OUT_DIR}/summary.json"

if [[ -s "$FINAL_PHYS" && -s "$SUMMARY_JSON" ]]; then
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
        echo "[$(date -u +%FT%TZ)] Signal shard ${TASK_ID} already completed (status=ok); skipping."
        exit 0
    fi
fi

echo "[$(date -u +%FT%TZ)] Starting signal shard ${TASK_ID} on $(hostname -f): rdo=${INPUT_RDO} cond=${COND_TAG}"

# Stage 1.8 GB RDO locally with xrdcp fallback
RDO_BASENAME="$(basename "$INPUT_RDO")"
LOCAL_RDO="${WORK_ROOT}/${RDO_BASENAME}"
if ! cp -f "$INPUT_RDO" "$LOCAL_RDO"; then
    echo "Direct cp failed for ${INPUT_RDO}; retrying via xrdcp..." >&2
    xrdcp -f "root://eospublic.cern.ch/${INPUT_RDO}" "$LOCAL_RDO"
fi
if [[ ! -s "$LOCAL_RDO" ]]; then
    echo "ERROR: Staged RDO is empty: $LOCAL_RDO" >&2
    exit 2
fi

# Clean Athena 24.0.41 + Calypso run/ setup in private scratch cwd
set +eu
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
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

export FASER_ACTS_MATERIAL_MAP="$MATERIAL_MAP"

RDO_STEM="${RDO_BASENAME%.root}"
if [[ "$RDO_STEM" == *-RDO ]]; then
    RDO_STEM="${RDO_STEM%-RDO}"
fi
XAOD_LOCAL="${WORK_ROOT}/reco/${RDO_STEM}-xAOD.root"
PHYS_LOCAL="${WORK_ROOT}/ntp/${TASK_ID}_phys.root"

T0=$(date +%s)
set +e
(
    cd "$WORK_ROOT/reco"
    /usr/bin/time -v -o "$WORK_ROOT/reco.time" python3 "${CALYPSO_SRC}/Control/CalypsoExample/Reconstruction/scripts/faser_reco.py" \
        "$LOCAL_RDO" \
        --isMC \
        --geom "$GEOM" \
        --cond "$COND_TAG" \
        --noIFT \
        --noBackward \
        --materialMap "$MATERIAL_MAP" \
        --nevents "$NEVENTS" \
        > "$WORK_ROOT/reco.log" 2>&1
)
RC_RECO=$?
T1=$(date +%s)
RECO_WALL=$((T1 - T0))

RC_NTP=99
NTP_WALL=0
XAOD_BYTES=0
if [[ $RC_RECO -eq 0 && -s "$XAOD_LOCAL" ]]; then
    XAOD_BYTES=$(stat -c%s "$XAOD_LOCAL")
    T2=$(date +%s)
    (
        cd "$WORK_ROOT/ntp"
        /usr/bin/time -v -o "$WORK_ROOT/ntp.time" python3 "${CALYPSO_SRC}/PhysicsAnalysis/NtupleDumper/scripts/faser_ntuple_maker.py" \
            "$XAOD_LOCAL" \
            --outfile "$PHYS_LOCAL" \
            --isMC \
            --geom "$GEOM" \
            --cond "$COND_TAG" \
            --materialMap "$MATERIAL_MAP" \
            --NoTrackFilt \
            --nevents "$NEVENTS" \
            > "$WORK_ROOT/ntp.log" 2>&1
    )
    RC_NTP=$?
    T3=$(date +%s)
    NTP_WALL=$((T3 - T2))
fi
set -e

# Always remove bulky scratch xAOD immediately
rm -f "$XAOD_LOCAL"

if [[ $RC_RECO -ne 0 || $RC_NTP -ne 0 || ! -s "$PHYS_LOCAL" ]]; then
    echo "FAILED ${TASK_ID}: rc_reco=${RC_RECO} (${RECO_WALL}s), rc_ntp=${RC_NTP} (${NTP_WALL}s)" >&2
    cp -f "$WORK_ROOT/reco.log" "${OUT_DIR}/failed_reco.log" 2>/dev/null || true
    cp -f "$WORK_ROOT/ntp.log" "${OUT_DIR}/failed_ntp.log" 2>/dev/null || true
    exit 1
fi

# Verify FASER-frame nt tree entries & sanity before staging to EOS
PHYS_STATS=$(python3 - "$PHYS_LOCAL" "$LOCAL_RDO" "$NEVENTS" <<'PY'
import json, math, statistics, sys, ROOT
phys_path, rdo_path, nev_str = sys.argv[1:4]
nev = int(nev_str)

f_rdo = ROOT.TFile.Open(rdo_path, "READ")
col = f_rdo.Get("CollectionTree")
rdo_entries = int(col.GetEntries()) if col else -1
f_rdo.Close()
expected = nev if nev > 0 else rdo_entries

f = ROOT.TFile.Open(phys_path, "READ")
t = f.Get("nt")
n = int(t.GetEntries()) if t else -1
z0s, z1s, pz0s = [], [], []
n_lt2 = 0
for i in range(max(0, n)):
    t.GetEntry(i)
    lt = int(t.longTracks)
    if lt == 2:
        n_lt2 += 1
    for z0 in t.Track_z0:
        z0s.append(float(z0))
    for z1 in t.Track_z1:
        z1s.append(float(z1))
    for pz0 in t.Track_pz0:
        pz0s.append(float(pz0))
f.Close()

med_z0 = statistics.median(z0s) if z0s else None
med_z1 = statistics.median(z1s) if z1s else None
pos_pz0_frac = (sum(1 for v in pz0s if v > 0.0) / len(pz0s)) if pz0s else 1.0
ok = (n == expected) and (med_z0 is None or (-20.0 <= med_z0 <= 80.0)) and (med_z1 is None or (2300.0 <= med_z1 <= 2550.0)) and (pos_pz0_frac >= 0.99)
stats = {
    "passed": ok,
    "nt_entries": n,
    "expected_entries": expected,
    "n_longTracks_eq_2": n_lt2,
    "total_tracks": len(z0s),
    "median_Track_z0_mm": med_z0,
    "median_Track_z1_mm": med_z1,
    "positive_pz0_fraction": pos_pz0_frac,
}
print(json.dumps(stats))
sys.exit(0 if ok else 10)
PY
)

PHYS_BYTES=$(stat -c%s "$PHYS_LOCAL")
cp -f "$PHYS_LOCAL" "${FINAL_PHYS}.tmp"
mv -f "${FINAL_PHYS}.tmp" "$FINAL_PHYS"

python3 - "$SUMMARY_JSON" "$TASK_ID" "$INPUT_RDO" "$FINAL_PHYS" "$COND_TAG" "$GEOM" \
          "$MATERIAL_MAP" "$CALYPSO_COMMIT" "$RECO_WALL" "$NTP_WALL" "$XAOD_BYTES" \
          "$PHYS_BYTES" "$PHYS_STATS" "$WORK_ROOT/reco.log" "$WORK_ROOT/reco.time" "$WORK_ROOT/ntp.time" <<'PY'
import datetime, json, os, re, socket, sys
from pathlib import Path

(
    summary_path, task_id, input_rdo, final_phys, cond_tag, geom,
    material_map, calypso_commit, reco_wall, ntp_wall, xaod_bytes,
    phys_bytes, phys_stats_s, reco_log_p, reco_time_p, ntp_time_p,
) = sys.argv[1:17]

reco_txt = Path(reco_log_p).read_text(errors="replace")
reco_time_txt = Path(reco_time_p).read_text(errors="replace")
ntp_time_txt = Path(ntp_time_p).read_text(errors="replace")
phys_stats = json.loads(phys_stats_s)

m_mat = re.search(r"Configured to use material input:\s*(\S+)", reco_txt)
m_gt = re.search(r"Global tag:\s*(\S+)", reco_txt)
m_rss_r = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", reco_time_txt)
m_rss_n = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", ntp_time_txt)

summary = {
    "status": "ok",
    "task_id": task_id,
    "dataset": "signal_110038",
    "is_mc": True,
    "input_rdo": input_rdo,
    "output_phys": final_phys,
    "frame": "FASER",
    "geom": geom,
    "cond_tag_requested": cond_tag,
    "cond_tag_resolved": m_gt.group(1) if m_gt else None,
    "material_map_requested": material_map,
    "material_map_loaded": m_mat.group(1) if m_mat else None,
    "calypso_commit": calypso_commit,
    "reco_wall_s": int(reco_wall),
    "ntp_wall_s": int(ntp_wall),
    "scratch_xaod_bytes": int(xaod_bytes),
    "phys_bytes": int(phys_bytes),
    "reco_peak_rss_kb": int(m_rss_r.group(1)) if m_rss_r else None,
    "ntp_peak_rss_kb": int(m_rss_n.group(1)) if m_rss_n else None,
    "phys_stats": phys_stats,
    "hostname": socket.getfqdn(),
    "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
}
tmp_path = Path(summary_path + ".tmp")
tmp_path.write_text(json.dumps(summary, indent=2) + "\n")
os.replace(tmp_path, summary_path)
print(json.dumps({"status": "ok", "task_id": task_id, "entries": phys_stats["nt_entries"], "reco_wall_s": int(reco_wall), "ntp_wall_s": int(ntp_wall)}))
PY

echo "[$(date -u +%FT%TZ)] Finished signal shard ${TASK_ID} cleanly."
