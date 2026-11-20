#!/usr/bin/env python3
"""Vectorized multi-process comparison of geoRotation_calypso_latest_fix_v1 full
reconstruction (RAW -> xAOD -> PHYS) against Ben's baseline (official r0021/r0022
PHYS + El9CalypsoForWaveForm reduced charge) across all 5 periods:
  2022 (007961), 2023 (011523), 2024_pre (015387), 2024_caloNu (015849), 2024_post (016963).
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from time import perf_counter
from typing import Any

import awkward as ak
import numpy as np
import uproot

DEFAULT_RUN_ROOT = Path(
    "/eos/home-x/xcheng/FASER/geoRotation_test_results/stage04_ben_comparison/run_20261120_v1"
)

# Histogram specifications matching FASER_DQ_Analyis (Fig. 5 & Fig. 6) + nu_muAnalysis cuts.yaml
HIST_SPECS: dict[str, tuple[int, float, float]] = {
    # Reduced charge & VetoNu (DQ + nu_mu signal-region zoom)
    "VetoNu0_reduced_charge": (60, 0.0, 300.0),
    "VetoNu1_reduced_charge": (60, 0.0, 300.0),
    "VetoNu0_reduced_charge_zoom": (60, 0.0, 60.0),
    "VetoNu1_reduced_charge_zoom": (60, 0.0, 60.0),
    "VetoNu0_reduced_charge_numu": (60, 0.0, 300.0),
    "VetoNu1_reduced_charge_numu": (60, 0.0, 300.0),
    "VetoNu0_reduced_charge_numu_zoom": (60, 0.0, 60.0),
    "VetoNu1_reduced_charge_numu_zoom": (60, 0.0, 60.0),
    "VetoNu0_charge": (50, 0.01, 300.0),
    "VetoNu1_charge": (50, 0.01, 300.0),
    "VetoNu0_raw_charge": (50, 0.01, 300.0),
    "VetoNu1_raw_charge": (50, 0.01, 300.0),
    # Scintillator & Calorimeter DQ (Fig. 5)
    "Veto20_charge": (50, 0.01, 300.0),
    "Veto21_charge": (50, 0.01, 300.0),
    "Timing0_charge": (50, 1.0, 80.0),
    "Timing1_charge": (50, 1.0, 80.0),
    "Timing2_charge": (50, 1.0, 80.0),
    "Timing3_charge": (50, 1.0, 80.0),
    "Timing_charge_total": (60, 0.0, 300.0),
    "Preshower0_charge": (50, 0.01, 15.0),
    "Preshower1_charge": (50, 0.01, 15.0),
    "Calo_total_E_EM": (60, 0.0, 150.0),  # GeV
    # Track DQ (Fig. 6 — in common physical FASER frame)
    "Track_Chi2": (50, 0.0, 50.0),
    "Track_nDoF": (20, 0.0, 20.0),
    "Track_chi2_per_dof": (50, 0.0, 25.0),
    "Track_nLayers": (7, 4.5, 11.5),
    "Track_pz_charge0": (100, -500.0, 500.0),  # GeV
    "Track_x0": (100, -100.0, 100.0),  # mm
    "Track_y0": (100, -100.0, 100.0),  # mm
    "Track_theta_x0": (50, -10.0, 10.0),  # mrad
    "Track_theta_y0": (50, -5.0, 5.0),  # mrad
    "Track_theta_x0_pos": (50, -10.0, 10.0),
    "Track_theta_x0_neg": (50, -10.0, 10.0),
    "Track_theta_y0_pos": (50, -5.0, 5.0),
    "Track_theta_y0_neg": (50, -5.0, 5.0),
    # nu_mu selection observables (LeadTrack, in common physical FASER frame)
    "LeadTrack_pz0": (60, 0.0, 1500.0),  # GeV
    "LeadTrack_Chi2_NDF": (50, 0.0, 25.0),
    "LeadTrack_nLayers": (7, 4.5, 11.5),
    "LeadTrack_nDoF": (20, 0.0, 20.0),
    "LeadTrack_Theta": (60, 0.0, 30.0),  # mrad
    "LeadTrack_r_atMaxRadius": (60, 0.0, 150.0),  # mm
    "LeadTrack_rIFT": (60, 0.0, 150.0),  # mm
    "LeadTrack_rVetoNu": (60, 0.0, 150.0),  # mm
    "longTracks": (7, -0.5, 6.5),
    # Event-by-event diffs & residuals
    "delta_redQ0": (80, -20.0, 20.0),
    "delta_redQ1": (80, -20.0, 20.0),
    "delta_redQ0_st0": (80, -20.0, 20.0),
    "delta_redQ0_st512": (80, -20.0, 20.0),
    "delta_qoverp_rel": (80, -0.5, 0.5),
    "delta_chi2": (80, -20.0, 20.0),
    "delta_x0": (80, -2.0, 2.0),
    "delta_y0": (80, -2.0, 2.0),
    "delta_thetax": (80, -1.0, 1.0),
    "delta_thetay": (80, -1.0, 1.0),
    # Gained vs lost tracks (only_fixv1 vs only_off)
    "gain_chi2_ndf": (50, 0.0, 25.0),
    "gain_pz0": (60, 0.0, 1000.0),
    "gain_nLayers": (7, 4.5, 11.5),
    "loss_chi2_ndf": (50, 0.0, 25.0),
    "loss_pz0": (60, 0.0, 1000.0),
    "loss_nLayers": (7, 4.5, 11.5),
    # Raw (un-rotated) NtupleDumper frame audit histograms
    "raw_Track_pz0": (100, -50.0, 50.0),  # GeV (shows rotated -px_F!)
    "raw_Track_px0": (100, -500.0, 500.0),  # GeV (shows rotated +pz_F!)
    "raw_Track_X_atVetoNu": (80, -150.0, 150.0),  # mm (shows -x_F sign flip!)
}


def make_hist(vals: np.ndarray, key: str) -> dict[str, Any]:
    nb, lo, hi = HIST_SPECS[key]
    v = np.asarray(vals, dtype=np.float64)
    v = v[np.isfinite(v)]
    counts, edges = np.histogram(v, bins=nb, range=(lo, hi))
    under = int(np.sum(v < lo))
    over = int(np.sum(v >= hi))
    return {
        "counts": counts.astype(int).tolist(),
        "edges": edges.tolist(),
        "n_total": int(v.size),
        "n_in_range": int(counts.sum()),
        "underflow": under,
        "overflow": over,
        "mean": float(np.mean(v)) if v.size else None,
        "std": float(np.std(v)) if v.size else None,
        "median": float(np.median(v)) if v.size else None,
    }


def load_nt_branches(paths: list[str], is_fixv1: bool, lo_eid: int | None = None, hi_eid: int | None = None) -> dict[str, Any]:
    """Read required branches from one or more PHYS ROOT files into contiguous arrays."""
    scalar_base = [
        "eventID", "eventTime", "distanceToCollidingBCID", "TAP", "inputBits", "longTracks",
        "VetoNu0_charge", "VetoNu1_charge", "VetoNu0_raw_charge", "VetoNu1_raw_charge",
        "VetoNu0_status", "VetoNu1_status",
        "Veto20_charge", "Veto21_charge", "Veto20_raw_charge", "Veto21_raw_charge",
        "Veto20_status", "Veto21_status",
        "Timing0_charge", "Timing1_charge", "Timing2_charge", "Timing3_charge",
        "Timing0_status", "Timing1_status", "Timing2_status", "Timing3_status",
        "Preshower0_charge", "Preshower1_charge", "Preshower0_raw_charge", "Preshower1_raw_charge",
        "Preshower0_status", "Preshower1_status",
    ]
    vec_base = [
        "Track_x0", "Track_y0", "Track_z0",
        "Track_px0", "Track_py0", "Track_pz0", "Track_p0",
        "Track_charge", "Track_Chi2", "Track_nDoF", "Track_nLayers",
        "Track_r_atMaxRadius",
        "Track_X_atVetoNu", "Track_Y_atVetoNu", "Track_ThetaX_atVetoNu",
        "Track_X_atVetoStation2", "Track_Y_atVetoStation2",
        "Track_X_atTrig", "Track_Y_atTrig",
    ]
    fixv1_extra = [
        "VetoNu0_reduced_charge", "VetoNu1_reduced_charge",
        "VetoNu0_reduced_charge_valid", "VetoNu1_reduced_charge_valid",
        "VetoNu_reduced_charge_fallback", "VetoNu_expected_localtime",
    ]

    with uproot.open(paths[0]) as f0:
        avail = set(f0["nt"].keys())
    calo_branch = "CaloLo_total_E_EM" if "CaloLo_total_E_EM" in avail else "Calo_total_E_EM"
    want = scalar_base + [calo_branch] + vec_base + (fixv1_extra if is_fixv1 else [])
    want = [b for b in want if b in avail]

    chunks = []
    for p in paths:
        with uproot.open(p) as f:
            t = f["nt"]
            if t.num_entries == 0:
                continue
            if lo_eid is not None and hi_eid is not None:
                eids = t["eventID"].array(library="np")
                mask = (eids >= lo_eid) & (eids <= hi_eid)
                if not np.any(mask):
                    continue
                arr = t.arrays(want, library="ak")[mask]
            else:
                arr = t.arrays(want, library="ak")
            chunks.append(arr)
    if not chunks:
        raise RuntimeError(f"No events loaded from {paths}")
    full = ak.concatenate(chunks, axis=0)

    out: dict[str, Any] = {}
    for b in want:
        if b in vec_base:
            out[b] = full[b]
        else:
            out[b] = ak.to_numpy(full[b])
    out["Calo_total_E_EM_MeV"] = out.pop(calo_branch)
    return out


def load_ben_wf(wf_paths: list[str]) -> dict[str, np.ndarray]:
    """Load Ben's El9CalypsoForWaveForm tree into sorted numpy arrays indexed by event."""
    want = [
        "event",
        "myVetoNu0_rawCharge", "myVetoNu0_rawCharge35",
        "myVetoNu1_rawCharge", "myVetoNu1_rawCharge35",
    ]
    opt = ["fallback_mode", "myVetoNu0_rawCharge35_valid", "myVetoNu1_rawCharge35_valid"]
    with uproot.open(wf_paths[0]) as f0:
        avail = set(f0["tree"].keys())
    cols = want + [c for c in opt if c in avail]
    parts = []
    for p in wf_paths:
        with uproot.open(p) as f:
            parts.append(f["tree"].arrays(cols, library="np"))
    merged = {c: np.concatenate([pt[c] for pt in parts]) for c in cols}
    order = np.argsort(merged["event"])
    return {c: merged[c][order] for c in cols}


def lookup_sorted(sorted_keys: np.ndarray, values: np.ndarray, query: np.ndarray, fill: float = np.nan):
    idx = np.searchsorted(sorted_keys, query)
    in_b = idx < sorted_keys.size
    safe_idx = np.where(in_b, idx, 0)
    hit = in_b & (sorted_keys[safe_idx] == query)
    res = np.full(query.shape, fill, dtype=np.float64)
    res[hit] = values[safe_idx[hit]]
    return res, hit


def grl_mask(event_time: np.ndarray, grl_json: str, run: int) -> tuple[np.ndarray, np.ndarray]:
    entry = json.loads(Path(grl_json).read_text())[str(run)]
    good = np.zeros(event_time.shape, dtype=bool)
    for iv in entry["stable_list"]:
        good |= (event_time >= iv["start_utime"]) & (event_time <= iv["stop_utime"])
    excl = np.zeros(event_time.shape, dtype=bool)
    for iv in entry.get("excluded_list", []):
        excl |= (event_time >= iv["start_utime"]) & (event_time <= iv["stop_utime"])
    return good, excl


def extract_physical_frame(data: dict[str, Any], is_fixv1: bool) -> dict[str, Any]:
    """Return track vectors in the common physical FASER frame.
    For geoRotation_calypso_latest_fix_v1, NtupleDumperAlg writes:
      Track_{x,y,z}0, Track_{px,py,pz}0 in rotated ACTS frame (x_rot=z_F, y_rot=y_F, z_rot=-x_F)
      Track_X_at<Station> as +z_rot = -x_F.
    Applying R^T maps fix_v1 strictly back to the physical FASER frame.
    """
    if is_fixv1:
        x0 = -data["Track_z0"]
        y0 = data["Track_y0"]
        z0 = data["Track_x0"]
        px0 = -data["Track_pz0"]
        py0 = data["Track_py0"]
        pz0 = data["Track_px0"]
        x_vn = -data["Track_X_atVetoNu"]
        y_vn = data["Track_Y_atVetoNu"]
        x_v2 = -data["Track_X_atVetoStation2"]
        y_v2 = data["Track_Y_atVetoStation2"]
        x_tr = -data["Track_X_atTrig"]
        y_tr = data["Track_Y_atTrig"]
    else:
        x0 = data["Track_x0"]
        y0 = data["Track_y0"]
        z0 = data["Track_z0"]
        px0 = data["Track_px0"]
        py0 = data["Track_py0"]
        pz0 = data["Track_pz0"]
        x_vn = data["Track_X_atVetoNu"]
        y_vn = data["Track_Y_atVetoNu"]
        x_v2 = data["Track_X_atVetoStation2"]
        y_v2 = data["Track_Y_atVetoStation2"]
        x_tr = data["Track_X_atTrig"]
        y_tr = data["Track_Y_atTrig"]

    p0 = data["Track_p0"]
    q = data["Track_charge"]
    chi2 = data["Track_Chi2"]
    ndof = data["Track_nDoF"]
    nlay = data["Track_nLayers"]
    rmax = data["Track_r_atMaxRadius"]
    r_vn = np.hypot(x_vn, y_vn)
    r_ift = np.hypot(x_v2, y_v2)
    pt0 = np.hypot(px0, py0)
    denom = np.maximum(p0, 1e-9)
    theta0 = np.arcsin(np.minimum(np.maximum(pt0 / denom, -1.0), 1.0))
    theta_x0 = np.arcsin(np.minimum(np.maximum(px0 / denom, -1.0), 1.0))
    theta_y0 = np.arcsin(np.minimum(np.maximum(py0 / denom, -1.0), 1.0))

    return {
        "x0": x0, "y0": y0, "z0": z0,
        "px0": px0, "py0": py0, "pz0": pz0, "p0": p0,
        "q": q, "chi2": chi2, "ndof": ndof, "nlay": nlay,
        "rmax": rmax, "r_vn": r_vn, "r_ift": r_ift,
        "x_vn": x_vn, "y_vn": y_vn, "x_tr": x_tr, "y_tr": y_tr,
        "theta0": theta0, "theta_x0": theta_x0, "theta_y0": theta_y0,
    }


def dedup_track_mask(p0_ak: ak.Array) -> ak.Array:
    """Replicate FASER_DQ RemoveDuplicates(Track_p0) (drop exact duplicate momenta per event)."""
    counts = ak.to_numpy(ak.num(p0_ak))
    flat = ak.to_numpy(ak.flatten(p0_ak))
    keep = np.ones(flat.shape, dtype=bool)
    # Fast vectorized check: multi-track events are rare (~5%), only loop over those with n >= 2
    offsets = np.zeros(counts.size + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    multi = np.flatnonzero(counts >= 2)
    for idx in multi:
        a, b = offsets[idx], offsets[idx + 1]
        sub = flat[a:b]
        for j in range(1, sub.size):
            if np.any(np.abs(sub[:j] - sub[j]) <= 1e-15):
                keep[a + j] = False
    return ak.unflatten(keep, counts)


def compute_lead_arrays(tf: dict[str, Any], n_ev: int) -> dict[str, np.ndarray]:
    """Extract leading-track (highest pz0_F) scalar arrays per event."""
    has_trk = ak.to_numpy(ak.num(tf["pz0"])) > 0
    lead_idx = ak.argmax(tf["pz0"], axis=1, keepdims=True)

    def first_or_nan(arr):
        sel = ak.firsts(arr[lead_idx])
        return ak.to_numpy(ak.fill_none(sel, np.nan)).astype(np.float64)

    return {
        "has_trk": has_trk,
        "pz0_GeV": first_or_nan(tf["pz0"]) * 1e-3,
        "p0_GeV": first_or_nan(tf["p0"]) * 1e-3,
        "px0_GeV": first_or_nan(tf["px0"]) * 1e-3,
        "py0_GeV": first_or_nan(tf["py0"]) * 1e-3,
        "q": first_or_nan(tf["q"]),
        "x0": first_or_nan(tf["x0"]),
        "y0": first_or_nan(tf["y0"]),
        "chi2": first_or_nan(tf["chi2"]),
        "ndof": first_or_nan(tf["ndof"]),
        "chi2_ndf": first_or_nan(tf["chi2"]) / np.maximum(first_or_nan(tf["ndof"]), 1e-9),
        "nlay": first_or_nan(tf["nlay"]),
        "theta_mrad": first_or_nan(tf["theta0"]) * 1e3,
        "theta_x_mrad": first_or_nan(tf["theta_x0"]) * 1e3,
        "theta_y_mrad": first_or_nan(tf["theta_y0"]) * 1e3,
        "rmax": first_or_nan(tf["rmax"]),
        "r_ift": first_or_nan(tf["r_ift"]),
        "r_vn": first_or_nan(tf["r_vn"]),
        "y_tr": first_or_nan(tf["y_tr"]),
        "x_vn": first_or_nan(tf["x_vn"]),
    }


def evaluate_side(
    data: dict[str, Any],
    tf: dict[str, Any],
    lead: dict[str, np.ndarray],
    red0: np.ndarray,
    red1: np.ndarray,
    aux_ok: np.ndarray,
    grl_json: str,
    run: int,
    is_calonu: bool,
) -> tuple[dict[str, Any], dict[str, int], dict[str, Any]]:
    """Compute exact FASER DQ histograms + exact nu_muAnalysis cutflow + nu_mu pre-VetoNu histograms."""
    good_t, excl_t = grl_mask(data["eventTime"], grl_json, run)
    colliding = data["distanceToCollidingBCID"] == 0
    no_t_sat = (
        ((data["Timing0_status"] & 4) == 0)
        & ((data["Timing1_status"] & 4) == 0)
        & ((data["Timing2_status"] & 4) == 0)
        & ((data["Timing3_status"] & 4) == 0)
    )
    tap4 = (data["TAP"] & 4) != 0
    dq_ev = good_t & (~excl_t) & colliding & tap4 & no_t_sat

    # GoodTracks (FASER_DQ_RDF.py definition in physical FASER frame)
    dedup = dedup_track_mask(tf["p0"])
    chi2_ndf_v = tf["chi2"] / ak.where(tf["ndof"] > 0, tf["ndof"], 1.0)
    good_trk_mask = (
        (tf["nlay"] >= 7)
        & (chi2_ndf_v < 25.0)
        & ((tf["ndof"] + 5) >= 12)
        & (tf["pz0"] > 20000.0)
        & dedup
    )
    dq_trk = good_trk_mask & dq_ev
    n_good_trk = int(ak.sum(dq_trk))
    n_good_pos = int(ak.sum(dq_trk & (tf["q"] > 0)))
    n_good_neg = int(ak.sum(dq_trk & (tf["q"] < 0)))

    # Flatten DQ track quantities (note: FASER DQ Fig. 6 fills ALL tracks in DQ events, matching analysis.py!)
    dq_all_trk_ev = dq_ev
    def flat_dq(arr):
        return ak.to_numpy(ak.flatten(arr[dq_all_trk_ev]))

    t_chi2 = flat_dq(tf["chi2"])
    t_ndof = flat_dq(tf["ndof"])
    t_chi2_ndf = flat_dq(chi2_ndf_v)
    t_nlay = flat_dq(tf["nlay"])
    t_pzq = flat_dq(tf["pz0"] * tf["q"]) * 1e-3
    t_x0 = flat_dq(tf["x0"])
    t_y0 = flat_dq(tf["y0"])
    t_thx = flat_dq(tf["theta_x0"]) * 1e3
    t_thy = flat_dq(tf["theta_y0"]) * 1e3
    t_q = flat_dq(tf["q"])

    # Exact nu_muAnalysis cutflow (cuts.yaml)
    trig_numu = (
        ((data["inputBits"] & 0x8) == 0x8)
        | ((data["inputBits"] & 0x10) == 0x10)
        | ((data["inputBits"] & 0x20) == 0x20)
        | ((data["inputBits"] & 0x40) == 0x40)
    )
    bad_v2 = (data["Veto20_status"] == 528) | (data["Veto21_status"] == 528)
    good_vn_st = (
        (((data["VetoNu0_status"] == 0) | (data["VetoNu0_status"] == 1))
         & ((data["VetoNu1_status"] == 0) | (data["VetoNu1_status"] == 1)))
        | ((data["VetoNu0_status"] == 0) & (data["VetoNu1_status"] == 16))
    )
    calonu_clean = np.ones_like(colliding) if not is_calonu else ((~bad_v2) & good_vn_st)
    v2_cut = (data["Veto20_charge"] > 40.0) & (data["Veto21_charge"] > 40.0)
    t_bot = data["Timing0_charge"] + data["Timing1_charge"]
    t_top = data["Timing2_charge"] + data["Timing3_charge"]
    t_tot = t_top + t_bot
    hits_timing = (
        ((lead["y_tr"] > 20.0) & (t_top > 20.0))
        | ((lead["y_tr"] < -20.0) & (t_bot > 20.0))
        | ((np.abs(lead["y_tr"]) <= 20.0) & (t_tot > 20.0))
    )
    ps_cut = (data["Preshower0_charge"] > 2.5) & (data["Preshower1_charge"] > 2.5)

    steps = [
        ("All_PHYS_events", np.ones_like(colliding)),
        ("Colliding", colliding),
        ("Trigger", trig_numu),
        ("GoodTimes_notExcluded", good_t & (~excl_t)),
        ("longTracks_gt0", data["longTracks"] > 0),
        ("CaloNu_status_clean", calonu_clean),
        ("Veto2_gt40pC", v2_cut),
        ("Timing_gt20pC", hits_timing),
        ("Preshower_gt2.5pC", ps_cut),
        ("LeadTrack_pz_gt100GeV", lead["pz0_GeV"] > 100.0),
        ("LeadTrack_quality", (lead["nlay"] >= 7) & (lead["ndof"] >= 9) & (lead["chi2_ndf"] < 15.0)),
        ("LeadTrack_fiducial_angle", (lead["rmax"] < 95.0) & (lead["r_ift"] < 95.0) & (lead["r_vn"] < 120.0) & (lead["theta_mrad"] < 25.0)),
        ("RedQ_available", aux_ok & np.isfinite(red0) & np.isfinite(red1)),
        ("VetoNu0_redQ_lt30", red0 < 30.0),
        ("VetoNu1_redQ_lt30", red1 < 30.0),
    ]
    cutflow: dict[str, int] = {}
    cur = np.ones_like(colliding)
    pre_vetonu_mask = None
    pre_kin_mask = None
    for sname, cond in steps:
        cur = cur & cond
        cutflow[sname] = int(cur.sum())
        if sname == "Preshower_gt2.5pC":
            pre_kin_mask = cur.copy()
        if sname == "RedQ_available":
            pre_vetonu_mask = cur.copy()

    # Book all histograms
    hists: dict[str, Any] = {}
    # DQ event-level scintillator & calo
    hists["VetoNu0_reduced_charge"] = make_hist(red0[dq_ev & aux_ok], "VetoNu0_reduced_charge")
    hists["VetoNu1_reduced_charge"] = make_hist(red1[dq_ev & aux_ok], "VetoNu1_reduced_charge")
    hists["VetoNu0_reduced_charge_zoom"] = make_hist(red0[dq_ev & aux_ok], "VetoNu0_reduced_charge_zoom")
    hists["VetoNu1_reduced_charge_zoom"] = make_hist(red1[dq_ev & aux_ok], "VetoNu1_reduced_charge_zoom")
    hists["VetoNu0_reduced_charge_numu"] = make_hist(red0[pre_vetonu_mask], "VetoNu0_reduced_charge_numu")
    hists["VetoNu1_reduced_charge_numu"] = make_hist(red1[pre_vetonu_mask], "VetoNu1_reduced_charge_numu")
    hists["VetoNu0_reduced_charge_numu_zoom"] = make_hist(red0[pre_vetonu_mask], "VetoNu0_reduced_charge_numu_zoom")
    hists["VetoNu1_reduced_charge_numu_zoom"] = make_hist(red1[pre_vetonu_mask], "VetoNu1_reduced_charge_numu_zoom")
    hists["VetoNu0_charge"] = make_hist(data["VetoNu0_charge"][dq_ev], "VetoNu0_charge")
    hists["VetoNu1_charge"] = make_hist(data["VetoNu1_charge"][dq_ev], "VetoNu1_charge")
    hists["VetoNu0_raw_charge"] = make_hist(data["VetoNu0_raw_charge"][dq_ev], "VetoNu0_raw_charge")
    hists["VetoNu1_raw_charge"] = make_hist(data["VetoNu1_raw_charge"][dq_ev], "VetoNu1_raw_charge")
    hists["Veto20_charge"] = make_hist(data["Veto20_charge"][dq_ev], "Veto20_charge")
    hists["Veto21_charge"] = make_hist(data["Veto21_charge"][dq_ev], "Veto21_charge")
    for k in ("Timing0_charge", "Timing1_charge", "Timing2_charge", "Timing3_charge"):
        hists[k] = make_hist(data[k][dq_ev], k)
    hists["Timing_charge_total"] = make_hist(t_tot[dq_ev], "Timing_charge_total")
    hists["Preshower0_charge"] = make_hist(data["Preshower0_charge"][dq_ev], "Preshower0_charge")
    hists["Preshower1_charge"] = make_hist(data["Preshower1_charge"][dq_ev], "Preshower1_charge")
    hists["Calo_total_E_EM"] = make_hist(data["Calo_total_E_EM_MeV"][dq_ev] * 1e-3, "Calo_total_E_EM")

    # DQ track-level (common physical FASER frame)
    hists["Track_Chi2"] = make_hist(t_chi2, "Track_Chi2")
    hists["Track_nDoF"] = make_hist(t_ndof, "Track_nDoF")
    hists["Track_chi2_per_dof"] = make_hist(t_chi2_ndf, "Track_chi2_per_dof")
    hists["Track_nLayers"] = make_hist(t_nlay, "Track_nLayers")
    hists["Track_pz_charge0"] = make_hist(t_pzq, "Track_pz_charge0")
    hists["Track_x0"] = make_hist(t_x0, "Track_x0")
    hists["Track_y0"] = make_hist(t_y0, "Track_y0")
    hists["Track_theta_x0"] = make_hist(t_thx, "Track_theta_x0")
    hists["Track_theta_y0"] = make_hist(t_thy, "Track_theta_y0")
    hists["Track_theta_x0_pos"] = make_hist(t_thx[t_q > 0], "Track_theta_x0_pos")
    hists["Track_theta_x0_neg"] = make_hist(t_thx[t_q < 0], "Track_theta_x0_neg")
    hists["Track_theta_y0_pos"] = make_hist(t_thy[t_q > 0], "Track_theta_y0_pos")
    hists["Track_theta_y0_neg"] = make_hist(t_thy[t_q < 0], "Track_theta_y0_neg")

    # nu_mu LeadTrack observables (after scintillator preselection, before kinematic cuts)
    hists["LeadTrack_pz0"] = make_hist(lead["pz0_GeV"][pre_kin_mask], "LeadTrack_pz0")
    hists["LeadTrack_Chi2_NDF"] = make_hist(lead["chi2_ndf"][pre_kin_mask], "LeadTrack_Chi2_NDF")
    hists["LeadTrack_nLayers"] = make_hist(lead["nlay"][pre_kin_mask], "LeadTrack_nLayers")
    hists["LeadTrack_nDoF"] = make_hist(lead["ndof"][pre_kin_mask], "LeadTrack_nDoF")
    hists["LeadTrack_Theta"] = make_hist(lead["theta_mrad"][pre_kin_mask], "LeadTrack_Theta")
    hists["LeadTrack_r_atMaxRadius"] = make_hist(lead["rmax"][pre_kin_mask], "LeadTrack_r_atMaxRadius")
    hists["LeadTrack_rIFT"] = make_hist(lead["r_ift"][pre_kin_mask], "LeadTrack_rIFT")
    hists["LeadTrack_rVetoNu"] = make_hist(lead["r_vn"][pre_kin_mask], "LeadTrack_rVetoNu")
    hists["longTracks"] = make_hist(data["longTracks"][dq_ev], "longTracks")

    # Raw (un-rotated) branches for explicit frame audit
    hists["raw_Track_pz0"] = make_hist(ak.to_numpy(ak.flatten(data["Track_pz0"][dq_ev])) * 1e-3, "raw_Track_pz0")
    hists["raw_Track_px0"] = make_hist(ak.to_numpy(ak.flatten(data["Track_px0"][dq_ev])) * 1e-3, "raw_Track_px0")
    hists["raw_Track_X_atVetoNu"] = make_hist(ak.to_numpy(ak.flatten(data["Track_X_atVetoNu"][dq_ev])), "raw_Track_X_atVetoNu")

    dq_counts = {
        "n_phys_total": int(data["eventID"].size),
        "n_dq_events": int(dq_ev.sum()),
        "n_dq_all_tracks": int(t_chi2.size),
        "n_dq_good_tracks": n_good_trk,
        "n_dq_good_pos": n_good_pos,
        "n_dq_good_neg": n_good_neg,
        "dq_ev_mask": dq_ev,
        "pre_vetonu_mask": pre_vetonu_mask,
    }
    return hists, cutflow, dq_counts


def analyze_period(item: tuple[str, dict[str, Any], list[dict[str, Any]], bool]) -> tuple[str, dict[str, Any]]:
    period, pspec, shards, allow_partial = item
    t0 = perf_counter()
    run = int(pspec["run"])
    is_calonu = (period == "2024_caloNu")

    # Collect completed shards
    phys_paths = []
    shard_summaries = []
    missing_shards = []
    for sh in shards:
        p = Path(sh["output_dir"]) / f"{sh['task_id']}_phys.root"
        s = Path(sh["output_dir"]) / "summary.json"
        if not p.is_file() or not s.is_file():
            if allow_partial:
                missing_shards.append(sh["task_id"])
                continue
            raise FileNotFoundError(f"Incomplete shard {sh['task_id']}: {p}")
        phys_paths.append(str(p))
        shard_summaries.append(json.loads(s.read_text()))

    # Load Ben's exact 2 chunks first to establish exact eventID window [lo, hi] of the 10 RAW segments
    ben_wf = load_ben_wf(pspec["ben_wf_files"])
    lo_eid = int(ben_wf["event"][0])
    hi_eid = int(ben_wf["event"][-1])

    # Load our fix_v1 PHYS shards and check duplicate eventIDs + range closure
    my_data = load_nt_branches(phys_paths, is_fixv1=True)
    my_eids = my_data["eventID"]
    n_unique_my = int(np.unique(my_eids).size)
    assert n_unique_my == my_eids.size, f"{period}: duplicate eventIDs in fix_v1 shards!"
    assert int(my_eids.min()) >= lo_eid and int(my_eids.max()) <= hi_eid, (
        f"{period}: fix_v1 eventID [{my_eids.min()}..{my_eids.max()}] outside Ben [{lo_eid}..{hi_eid}]"
    )
    # Sort fix_v1 by eventID
    ord_my = np.argsort(my_eids)
    for k, v in list(my_data.items()):
        my_data[k] = v[ord_my]
    my_eids = my_data["eventID"]

    # Load official PHYS slice filtered to exact same [lo_eid, hi_eid] (and covered shard windows if partial)
    off_data = load_nt_branches([pspec["off_phys"]], is_fixv1=False, lo_eid=lo_eid, hi_eid=hi_eid)
    if missing_shards:
        cov = np.zeros(off_data["eventID"].shape, dtype=bool)
        for s in shard_summaries:
            cov |= (off_data["eventID"] >= s["eventID_min"]) & (off_data["eventID"] <= s["eventID_max"])
        for k, v in list(off_data.items()):
            off_data[k] = v[cov]
    ord_off = np.argsort(off_data["eventID"])
    for k, v in list(off_data.items()):
        off_data[k] = v[ord_off]
    off_eids = off_data["eventID"]

    # Join Ben's reduced charge onto official PHYS (Ben's exact analysis pipeline!)
    ben_red0_on_off, ben_hit_off = lookup_sorted(ben_wf["event"], ben_wf["myVetoNu0_rawCharge35"], off_eids)
    ben_red1_on_off, _ = lookup_sorted(ben_wf["event"], ben_wf["myVetoNu1_rawCharge35"], off_eids)
    ben_raw0_on_off, _ = lookup_sorted(ben_wf["event"], ben_wf["myVetoNu0_rawCharge"], off_eids)
    ben_raw1_on_off, _ = lookup_sorted(ben_wf["event"], ben_wf["myVetoNu1_rawCharge"], off_eids)
    # Apply Ben's fallback rule from cuts.yaml/Analyis.cxx (if myVetoNu*_rawCharge35 < 0 -> raw charge)
    ben_red0_eff = np.where(ben_red0_on_off >= 0.0, ben_red0_on_off, ben_raw0_on_off)
    ben_red1_eff = np.where(ben_red1_on_off >= 0.0, ben_red1_on_off, ben_raw1_on_off)

    # Physical FASER-frame kinematics
    my_tf = extract_physical_frame(my_data, is_fixv1=True)
    off_tf = extract_physical_frame(off_data, is_fixv1=False)
    my_lead = compute_lead_arrays(my_tf, my_eids.size)
    off_lead = compute_lead_arrays(off_tf, off_eids.size)

    my_red0 = my_data["VetoNu0_reduced_charge"].astype(np.float64)
    my_red1 = my_data["VetoNu1_reduced_charge"].astype(np.float64)
    my_red_ok = np.isfinite(my_red0) & np.isfinite(my_red1) & (my_red0 >= -0.5) & (my_red1 >= -0.5)

    h_my, cf_my, dq_my = evaluate_side(
        my_data, my_tf, my_lead, my_red0, my_red1, my_red_ok, pspec["grl_json"], run, is_calonu
    )
    h_ben, cf_ben, dq_ben = evaluate_side(
        off_data, off_tf, off_lead, ben_red0_eff, ben_red1_eff, ben_hit_off, pspec["grl_json"], run, is_calonu
    )

    # ---- EVENT-BY-EVENT MATCHED COMPARISON ----
    # 1) Reduced charge comparison against Ben WF on ALL fix_v1 DQ events!
    dq_m = dq_my["dq_ev_mask"]
    ben_r0_on_my, hit_b_my = lookup_sorted(ben_wf["event"], ben_wf["myVetoNu0_rawCharge35"], my_eids)
    ben_r1_on_my, _ = lookup_sorted(ben_wf["event"], ben_wf["myVetoNu1_rawCharge35"], my_eids)
    ben_raw0_on_my, _ = lookup_sorted(ben_wf["event"], ben_wf["myVetoNu0_rawCharge"], my_eids)
    ben_raw1_on_my, _ = lookup_sorted(ben_wf["event"], ben_wf["myVetoNu1_rawCharge"], my_eids)
    ben_r0_my_eff = np.where(ben_r0_on_my >= 0.0, ben_r0_on_my, ben_raw0_on_my)
    ben_r1_my_eff = np.where(ben_r1_on_my >= 0.0, ben_r1_on_my, ben_raw1_on_my)

    m_rc = dq_m & hit_b_my & np.isfinite(my_red0) & np.isfinite(ben_r0_my_eff)
    d_r0 = my_red0[m_rc] - ben_r0_my_eff[m_rc]
    d_r1 = my_red1[m_rc] - ben_r1_my_eff[m_rc]
    fb_arr = my_data["VetoNu_reduced_charge_fallback"][m_rc]
    v20_st = my_data["Veto20_status"][m_rc]
    v21_st = my_data["Veto21_status"][m_rc]
    is_st512 = ((v20_st & 512) != 0) | ((v21_st & 512) != 0)
    is_st0 = (~is_st512) & (v20_st == 0) & (v21_st == 0)

    def summarize_diff(arr: np.ndarray) -> dict[str, Any]:
        if arr.size == 0:
            return {"n": 0}
        a = np.abs(arr)
        return {
            "n": int(arr.size),
            "median": float(np.median(arr)),
            "mean": float(np.mean(arr)),
            "std": float(np.std(arr)),
            "q01": float(np.quantile(arr, 0.01)),
            "q05": float(np.quantile(arr, 0.05)),
            "q95": float(np.quantile(arr, 0.95)),
            "q99": float(np.quantile(arr, 0.99)),
            "frac_lt_0p01": float(np.mean(a < 0.01)),
            "frac_lt_0p1": float(np.mean(a < 0.1)),
            "frac_lt_1p0": float(np.mean(a < 1.0)),
            "frac_gt_5p0": float(np.mean(a > 5.0)),
        }

    # Check fallback breakdown & outlier diagnostics for |d_r0| > 1.0 pC
    big_idx = np.flatnonzero(np.abs(d_r0) > 1.0)
    big_reasons = {
        "n_gt1pC": int(big_idx.size),
        "fb0_raw_fallback": int(np.sum(fb_arr[big_idx] == 0)),
        "fb1_veto2": int(np.sum(fb_arr[big_idx] == 1)),
        "fb2_ps_dig1": int(np.sum(fb_arr[big_idx] == 2)),
        "fb3_ps_dig2": int(np.sum(fb_arr[big_idx] == 3)),
        "ben_was_neg1_fallback": int(np.sum(ben_r0_on_my[m_rc][big_idx] < 0.0)),
        "v2_has_st512": int(np.sum(is_st512[big_idx])),
        "v2_has_overflow4": int(np.sum(((v20_st[big_idx] & 4) != 0) | ((v21_st[big_idx] & 4) != 0))),
    }

    # Check reduced-charge < 30 pC agreement on nu_mu pre-VetoNu selection!
    pre_m = dq_my["pre_vetonu_mask"] & hit_b_my
    my_pass30 = (my_red0[pre_m] < 30.0) & (my_red1[pre_m] < 30.0)
    ben_pass30 = (ben_r0_my_eff[pre_m] < 30.0) & (ben_r1_my_eff[pre_m] < 30.0)

    redq_ebye = {
        "all_vn0": summarize_diff(d_r0),
        "all_vn1": summarize_diff(d_r1),
        "fb_counts": {str(k): int(np.sum(fb_arr == k)) for k in (0, 1, 2, 3)},
        "st0_vn0": summarize_diff(d_r0[is_st0]),
        "st512_vn0": summarize_diff(d_r0[is_st512]),
        "fb1_vn0": summarize_diff(d_r0[fb_arr == 1]),
        "fb2_vn0": summarize_diff(d_r0[fb_arr == 2]),
        "outliers_gt1pC": big_reasons,
        "numu_pre_n": int(pre_m.sum()),
        "numu_both_pass30": int(np.sum(my_pass30 & ben_pass30)),
        "numu_only_my_pass30": int(np.sum(my_pass30 & (~ben_pass30))),
        "numu_only_ben_pass30": int(np.sum((~my_pass30) & ben_pass30)),
    }

    # 2) Event-by-event tracking comparison between fix_v1 and official PHYS on DQ events
    dq_my_eids = my_eids[dq_my["dq_ev_mask"]]
    dq_off_eids = off_eids[dq_ben["dq_ev_mask"]]
    common_eids, i_my_sub, i_off_sub = np.intersect1d(dq_my_eids, dq_off_eids, return_indices=True)
    idx_my_dq = np.flatnonzero(dq_my["dq_ev_mask"])
    idx_off_dq = np.flatnonzero(dq_ben["dq_ev_mask"])
    im = idx_my_dq[i_my_sub]
    io = idx_off_dq[i_off_sub]

    only_my_mask = ~np.isin(dq_my_eids, dq_off_eids)
    only_off_mask = ~np.isin(dq_off_eids, dq_my_eids)
    idx_only_my = idx_my_dq[only_my_mask]
    idx_only_off = idx_off_dq[only_off_mask]

    both_lead = my_lead["has_trk"][im] & off_lead["has_trk"][io]
    im_b = im[both_lead]
    io_b = io[both_lead]

    qp_my = my_lead["q"][im_b] / np.maximum(my_lead["p0_GeV"][im_b], 1e-6)
    qp_off = off_lead["q"][io_b] / np.maximum(off_lead["p0_GeV"][io_b], 1e-6)
    # Relative signed curvature shift: (q/p_my - q/p_off) / |1/p_off|
    d_qp_rel = (qp_my - qp_off) * off_lead["p0_GeV"][io_b]
    d_chi2 = my_lead["chi2"][im_b] - off_lead["chi2"][io_b]
    valid_ndof = (my_lead["ndof"][im_b] >= 0) & (off_lead["ndof"][io_b] >= 0)
    d_ndof = my_lead["ndof"][im_b][valid_ndof] - off_lead["ndof"][io_b][valid_ndof]
    d_x0 = my_lead["x0"][im_b] - off_lead["x0"][io_b]
    d_y0 = my_lead["y0"][im_b] - off_lead["y0"][io_b]
    d_thx = my_lead["theta_x_mrad"][im_b] - off_lead["theta_x_mrad"][io_b]
    d_thy = my_lead["theta_y_mrad"][im_b] - off_lead["theta_y_mrad"][io_b]
    q_same = (my_lead["q"][im_b] == off_lead["q"][io_b])

    # Charge agreement & curvature shift binned by official momentum
    p_bins = [(0, 50), (50, 150), (150, 400), (400, 1000), (1000, 1e9)]
    p_binned = []
    p_off_b = off_lead["p0_GeV"][io_b]
    for plo, phi in p_bins:
        pm = (p_off_b >= plo) & (p_off_b < phi)
        if np.any(pm):
            p_binned.append({
                "p_range_GeV": f"[{plo},{phi if phi < 1e8 else 'inf'})",
                "n": int(pm.sum()),
                "q_agree_frac": float(np.mean(q_same[pm])),
                "d_qp_rel_med": float(np.median(d_qp_rel[pm])),
                "d_qp_rel_iqr": float(np.quantile(d_qp_rel[pm], 0.75) - np.quantile(d_qp_rel[pm], 0.25)),
                "p_ratio_med": float(np.median(my_lead["p0_GeV"][im_b][pm] / p_off_b[pm])),
            })

    # Raw frame audit on matched tracks (demonstrating NtupleDumperAlg.cxx frame swap)
    raw_my_px0 = ak.to_numpy(ak.fill_none(ak.firsts(my_data["Track_px0"][im_b]), np.nan)) * 1e-3
    raw_my_pz0 = ak.to_numpy(ak.fill_none(ak.firsts(my_data["Track_pz0"][im_b]), np.nan)) * 1e-3
    raw_my_xvn = ak.to_numpy(ak.fill_none(ak.firsts(my_data["Track_X_atVetoNu"][im_b]), np.nan))
    raw_off_xvn = ak.to_numpy(ak.fill_none(ak.firsts(off_data["Track_X_atVetoNu"][io_b]), np.nan))
    raw_my_thxvn = ak.to_numpy(ak.fill_none(ak.firsts(my_data["Track_ThetaX_atVetoNu"][im_b]), np.nan))
    frame_audit = {
        "n_matched": int(im_b.size),
        "raw_px0_med_GeV": float(np.nanmedian(raw_my_px0)),
        "raw_pz0_med_GeV": float(np.nanmedian(raw_my_pz0)),
        "off_pz0_med_GeV": float(np.nanmedian(off_lead["pz0_GeV"][io_b])),
        "raw_X_atVetoNu_plus_off_med_mm": float(np.nanmedian(raw_my_xvn + raw_off_xvn)),
        "raw_X_atVetoNu_minus_off_std_mm": float(np.nanstd(raw_my_xvn - raw_off_xvn)),
        "raw_ThetaX_atVetoNu_abs_med_rad": float(np.nanmedian(np.abs(raw_my_thxvn))),
    }

    # Extra event-by-event histograms
    h_ebye = {
        "delta_redQ0": make_hist(d_r0, "delta_redQ0"),
        "delta_redQ1": make_hist(d_r1, "delta_redQ1"),
        "delta_redQ0_st0": make_hist(d_r0[is_st0], "delta_redQ0_st0"),
        "delta_redQ0_st512": make_hist(d_r0[is_st512], "delta_redQ0_st512"),
        "delta_qoverp_rel": make_hist(d_qp_rel, "delta_qoverp_rel"),
        "delta_chi2": make_hist(d_chi2, "delta_chi2"),
        "delta_x0": make_hist(d_x0, "delta_x0"),
        "delta_y0": make_hist(d_y0, "delta_y0"),
        "delta_thetax": make_hist(d_thx, "delta_thetax"),
        "delta_thetay": make_hist(d_thy, "delta_thetay"),
        "gain_chi2_ndf": make_hist(my_lead["chi2_ndf"][idx_only_my], "gain_chi2_ndf"),
        "gain_pz0": make_hist(my_lead["pz0_GeV"][idx_only_my], "gain_pz0"),
        "gain_nLayers": make_hist(my_lead["nlay"][idx_only_my], "gain_nLayers"),
        "loss_chi2_ndf": make_hist(off_lead["chi2_ndf"][idx_only_off], "loss_chi2_ndf"),
        "loss_pz0": make_hist(off_lead["pz0_GeV"][idx_only_off], "loss_pz0"),
        "loss_nLayers": make_hist(off_lead["nlay"][idx_only_off], "loss_nLayers"),
    }

    trk_ebye = {
        "dq_events_fixv1": int(dq_my_eids.size),
        "dq_events_off": int(dq_off_eids.size),
        "dq_events_common": int(common_eids.size),
        "dq_events_only_fixv1": int(idx_only_my.size),
        "dq_events_only_off": int(idx_only_off.size),
        "net_dq_event_gain_pct": float(100.0 * (dq_my_eids.size - dq_off_eids.size) / max(dq_off_eids.size, 1)),
        "net_good_trk_gain_pct": float(
            100.0 * (dq_my["n_dq_good_tracks"] - dq_ben["n_dq_good_tracks"]) / max(dq_ben["n_dq_good_tracks"], 1)
        ),
        "q_agree_all": float(np.mean(q_same)) if q_same.size else None,
        "d_chi2_med": float(np.median(d_chi2)) if d_chi2.size else None,
        "d_ndof_mean": float(np.mean(d_ndof)) if d_ndof.size else None,
        "d_x0_med_mm": float(np.median(d_x0)) if d_x0.size else None,
        "d_x0_std_mm": float(np.std(d_x0[np.abs(d_x0) < 5.0])) if d_x0.size else None,
        "d_y0_med_mm": float(np.median(d_y0)) if d_y0.size else None,
        "d_y0_std_mm": float(np.std(d_y0[np.abs(d_y0) < 5.0])) if d_y0.size else None,
        "d_thx_med_mrad": float(np.median(d_thx)) if d_thx.size else None,
        "d_thy_med_mrad": float(np.median(d_thy)) if d_thy.size else None,
        "d_qp_rel_med": float(np.median(d_qp_rel)) if d_qp_rel.size else None,
        "by_p_bin": p_binned,
        "only_fixv1_good_chi2_frac": float(np.mean(my_lead["chi2_ndf"][idx_only_my] < 15.0)) if idx_only_my.size else None,
        "only_off_good_chi2_frac": float(np.mean(off_lead["chi2_ndf"][idx_only_off] < 15.0)) if idx_only_off.size else None,
    }

    result = {
        "period": period,
        "run": run,
        "segments": pspec["segments"],
        "ben_parts": pspec["ben_parts"],
        "cond_tag": pspec["cond_tag"],
        "timing_preset": pspec["timing_preset"],
        "eventID_window": [lo_eid, hi_eid],
        "n_shards": len(shard_summaries),
        "raw_events_total": int(sum(s["raw_events_processed"] or 0 for s in shard_summaries)),
        "mean_reco_wall_s": float(np.mean([s["reco_wall_s"] for s in shard_summaries])),
        "mean_ntp_wall_s": float(np.mean([s["ntp_wall_s"] for s in shard_summaries])),
        "dq_counts_fixv1": {k: v for k, v in dq_my.items() if not isinstance(v, np.ndarray)},
        "dq_counts_ben": {k: v for k, v in dq_ben.items() if not isinstance(v, np.ndarray)},
        "cutflow_fixv1": cf_my,
        "cutflow_ben": cf_ben,
        "redq_ebye": redq_ebye,
        "trk_ebye": trk_ebye,
        "frame_audit": frame_audit,
        "hists_fixv1": h_my,
        "hists_ben": h_ben,
        "hists_ebye": h_ebye,
        "analysis_wall_s": round(perf_counter() - t0, 2),
    }
    return period, result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()

    t_start = perf_counter()
    periods = json.loads((args.run_root / "periods.json").read_text())
    tasks = json.loads((args.run_root / "tasks.json").read_text())
    by_period = {}
    for t in tasks:
        by_period.setdefault(t["period"], []).append(t)

    order = ["2022", "2023", "2024_pre", "2024_caloNu", "2024_post"]
    work = [(p, periods[p], by_period[p], args.allow_partial) for p in order if p in periods]

    results = {}
    with ProcessPoolExecutor(max_workers=min(args.workers, len(work))) as pool:
        for p, res in pool.map(analyze_period, work):
            results[p] = res
            print(
                f"[{p:12s}] run={res['run']:06d} | DQ_ev: fix_v1={res['dq_counts_fixv1']['n_dq_events']} vs Ben={res['dq_counts_ben']['n_dq_events']} "
                f"({res['trk_ebye']['net_dq_event_gain_pct']:+.2f}%) | GoodTrk: {res['dq_counts_fixv1']['n_dq_good_tracks']} vs {res['dq_counts_ben']['n_dq_good_tracks']} "
                f"({res['trk_ebye']['net_good_trk_gain_pct']:+.2f}%) | RedQ0 |d|<0.01: {100*res['redq_ebye']['all_vn0']['frac_lt_0p01']:.2f}% "
                f"| wall={res['analysis_wall_s']}s",
                flush=True,
            )

    out_dir = args.run_root / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    full_json = out_dir / "ben_vs_fixv1_full.json"
    full_json.write_text(json.dumps(results, indent=2) + "\n")

    # Compact summary without histogram bin arrays for quick reading
    summary = {
        p: {k: v for k, v in r.items() if not k.startswith("hists_")}
        for p, r in results.items()
    }
    (out_dir / "ben_vs_fixv1_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Saved full results to {full_json} in {perf_counter() - t_start:.1f}s")


if __name__ == "__main__":
    main()
