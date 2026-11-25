#!/usr/bin/env python3
"""Shared visual styling, color palettes, axes builders, and figure helpers for Stage 08 overlay plots.

Design goals:
- Echo Pawan's FASER Tracking Meeting report layout (2-pad efficiency + ratio, dual-axis
  green log10 A' signal step histogram, red MC / blue Data markers with +/- half-bin
  x-error bars and Clopper-Pearson / McNemar y-error bars, crisp FASER Preliminary header).
- All text strictly in English.
- Zero label clipping or legend overlap; compatible with Matplotlib 3.4+.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.ticker import AutoMinorLocator, LogLocator, MultipleLocator, NullFormatter
import numpy as np


# ---------------------------------------------------------------------------
# Color palette & marker definitions
# ---------------------------------------------------------------------------
COLORS: Dict[str, str] = {
    "mc": "#d62728",           # Crisp red (Muon MC Overlay)
    "data": "#1557d6",         # Crisp blue (Data Overlay)
    "mc_base": "#b82020",      # Baseline MC (open triangle)
    "data_base": "#1245a8",    # Baseline Data (open triangle)
    "may30_mc": "#c44e52",     # Baseline-May30 MC reference
    "may30_data": "#4c72b0",   # Baseline-May30 Data reference
    "signal": "#1b8a2a",       # Forest green (A' MC signal distribution)
    "signal_latest": "#0d5e18",  # Darker green dashed (A' Latest GeoRotation check)
    "ratio": "#111111",        # Black (MC / Data ratio)
    "ref_line": "#777777",     # Grey dashed reference line
    "grid": "#d9d9d9",
    "footer": "#333333",
    "slate": "#555555",
}

# Colorblind-safe Okabe-Ito / Tol qualitative palette for diagnostics (Suite C)
DIAG_PALETTE: List[str] = [
    "#0072B2",  # Blue
    "#D55E00",  # Vermilion
    "#009E73",  # Bluish green
    "#CC79A7",  # Reddish purple
    "#E69F00",  # Orange
    "#56B4E9",  # Sky blue
    "#333333",  # Charcoal
    "#882255",  # Wine
]

SIGNAL_LABEL = r"A' MC ($m = 100\ \mathrm{MeV},\ \varepsilon = 10^{-5}$)"
SIGNAL_LATEST_LABEL = r"A' MC Latest GeoRot ($m = 100\ \mathrm{MeV},\ \varepsilon = 10^{-5}$)"
HEADER_PREFIX = r"$\mathbf{FASER}\ \mathit{Preliminary}$"


def apply_style() -> None:
    """Apply consistent publication-grade Matplotlib rcParams (Matplotlib 3.4+ safe)."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
            "mathtext.default": "regular",
            "font.size": 10.0,
            "axes.titlesize": 10.5,
            "axes.labelsize": 10.5,
            "xtick.labelsize": 9.2,
            "ytick.labelsize": 9.2,
            "legend.fontsize": 8.4,
            "axes.linewidth": 0.9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 5.0,
            "xtick.minor.size": 2.6,
            "ytick.major.size": 5.0,
            "ytick.minor.size": 2.6,
            "xtick.top": True,
            "ytick.right": False,
            "axes.grid": False,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "pdf.fonttype": 42,  # TrueType fonts embedded in vector PDF
            "ps.fonttype": 42,
        }
    )


def format_faser_header(
    ax: Axes,
    subtitle: str = "Overlay Efficiency",
    right_text: str = "",
    fontsize: float = 9.8,
    right_fontsize: float = 8.2,
) -> None:
    """Place the 'FASER Preliminary - <subtitle>' header cleanly above the top spine."""
    title_str = f"{HEADER_PREFIX} — {subtitle}" if subtitle else HEADER_PREFIX
    ax.text(
        0.0,
        1.018,
        title_str,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=fontsize,
        color="#111111",
    )
    if right_text:
        ax.text(
            1.0,
            1.018,
            right_text,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=right_fontsize,
            color="#333333",
        )


def add_figure_banners(
    fig: Figure,
    title: str,
    scope_line: str,
    footer_line: str,
    title_y: float = 0.975,
    scope_y: float = 0.938,
    footer_y: float = 0.012,
    footer_fontsize: float = 7.4,
) -> None:
    """Add standardized top title, dataset/metric scope subtitle, and bottom provenance footer."""
    fig.text(
        0.5,
        title_y,
        title,
        ha="center",
        va="top",
        fontsize=12.5,
        fontweight="bold",
        color="#111111",
    )
    if scope_line:
        fig.text(
            0.5,
            scope_y,
            scope_line,
            ha="center",
            va="top",
            fontsize=8.8,
            color="#333333",
        )
    if footer_line:
        fig.text(
            0.5,
            footer_y,
            footer_line,
            ha="center",
            va="bottom",
            fontsize=footer_fontsize,
            color=COLORS["footer"],
        )


# ---------------------------------------------------------------------------
# Multi-panel grid builders
# ---------------------------------------------------------------------------
def make_two_column_ratio_figure(
    figsize: Tuple[float, float] = (14.2, 6.1),
    height_ratios: Sequence[float] = (3.0, 1.15),
    left: float = 0.065,
    right: float = 0.935,
    top: float = 0.855,
    bottom: float = 0.115,
    wspace: float = 0.34,
    hspace: float = 0.055,
) -> Tuple[Figure, List[Tuple[Axes, ...]]]:
    """Build a 2-column figure where each column has stacked (top, bottom[, third]) panels.

    Returns:
        fig: Matplotlib Figure.
        columns: List of 2 tuples, each containing the stacked Axes from top to bottom
                 sharing the x-axis.
    """
    fig = plt.figure(figsize=figsize)
    outer = GridSpec(1, 2, figure=fig, left=left, right=right, top=top, bottom=bottom, wspace=wspace)
    columns: List[Tuple[Axes, ...]] = []
    n_rows = len(height_ratios)
    for col_idx in range(2):
        inner = GridSpecFromSubplotSpec(
            n_rows, 1, subplot_spec=outer[0, col_idx], height_ratios=height_ratios, hspace=hspace
        )
        col_axes: List[Axes] = []
        for row_idx in range(n_rows):
            sharex = col_axes[0] if row_idx > 0 else None
            ax = fig.add_subplot(inner[row_idx, 0], sharex=sharex)
            if row_idx < n_rows - 1:
                plt.setp(ax.get_xticklabels(), visible=False)
            col_axes.append(ax)
        columns.append(tuple(col_axes))
    return fig, columns


def make_single_ratio_figure(
    figsize: Tuple[float, float] = (8.6, 6.2),
    height_ratios: Sequence[float] = (2.8, 1.15),
    left: float = 0.11,
    right: float = 0.94,
    top: float = 0.855,
    bottom: float = 0.12,
    hspace: float = 0.06,
) -> Tuple[Figure, Tuple[Axes, ...]]:
    """Build a single-column stacked (top, bottom) ratio figure."""
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(
        len(height_ratios),
        1,
        figure=fig,
        height_ratios=height_ratios,
        left=left,
        right=right,
        top=top,
        bottom=bottom,
        hspace=hspace,
    )
    axes: List[Axes] = []
    for row_idx in range(len(height_ratios)):
        sharex = axes[0] if row_idx > 0 else None
        ax = fig.add_subplot(gs[row_idx, 0], sharex=sharex)
        if row_idx < len(height_ratios) - 1:
            plt.setp(ax.get_xticklabels(), visible=False)
        axes.append(ax)
    return fig, tuple(axes)


# ---------------------------------------------------------------------------
# Axis styling & plotting primitives
# ---------------------------------------------------------------------------
def configure_efficiency_axis(
    ax: Axes,
    xlim: Tuple[float, float],
    ylim: Tuple[float, float] = (0.80, 1.00),
    ylabel: str = "Overlay Efficiency",
    major_ytick: float = 0.05,
    minor_ytick: float = 0.01,
    has_right_twin: bool = True,
) -> None:
    """Style the primary efficiency top axis matching Pawan's tick conventions."""
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel)
    ax.yaxis.set_major_locator(MultipleLocator(major_ytick))
    ax.yaxis.set_minor_locator(MultipleLocator(minor_ytick))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=(not has_right_twin))
    ax.grid(axis="y", which="major", color=COLORS["grid"], linestyle=":", linewidth=0.6, alpha=0.65)


def configure_ratio_axis(
    ax: Axes,
    xlim: Tuple[float, float],
    xlabel: str,
    ylabel: str = "MC / Data",
    ylim: Tuple[float, float] = (0.95, 1.05),
    ref_value: float = 1.0,
    major_ytick: Optional[float] = 0.02,
    minor_ytick: Optional[float] = 0.005,
) -> None:
    """Style the bottom ratio or paired-delta axis."""
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if major_ytick is not None and major_ytick > 0:
        ax.yaxis.set_major_locator(MultipleLocator(major_ytick))
    if minor_ytick is not None and minor_ytick > 0:
        ax.yaxis.set_minor_locator(MultipleLocator(minor_ytick))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax.axhline(ref_value, color=COLORS["ref_line"], linestyle="--", linewidth=0.9, zorder=1)
    ax.grid(axis="y", which="major", color=COLORS["grid"], linestyle=":", linewidth=0.55, alpha=0.55)


def draw_signal_overlay(
    ax: Axes,
    edges: np.ndarray,
    density: np.ndarray,
    label: str = SIGNAL_LABEL,
    color: str = COLORS["signal"],
    linestyle: str = "-",
    ylim: Tuple[float, float] = (8e-5, 0.7),
    existing_twin: Optional[Axes] = None,
    linewidth: float = 1.05,
    alpha: float = 0.90,
) -> Tuple[Axes, Any]:
    """Draw normalized log10 A' signal step histogram on a right-hand twin axis.

    Returns (ax_right, line_handle) for unified legend assembly.
    """
    ax_r = existing_twin if existing_twin is not None else ax.twinx()
    edges = np.asarray(edges, dtype=float)
    vals = np.asarray(density, dtype=float)
    # Mask zero/negative bins slightly below ylim floor so log step curve drops cleanly to bottom
    floor_val = ylim[0] * 0.5
    safe_vals = np.where(np.isfinite(vals) & (vals > 0.0), vals, floor_val)

    # Build step vertices (stairs-compatible across matplotlib 3.4+)
    step_y = np.concatenate([safe_vals, safe_vals[-1:]])
    (handle,) = ax_r.step(
        edges,
        step_y,
        where="post",
        color=color,
        linestyle=linestyle,
        linewidth=linewidth,
        alpha=alpha,
        label=label,
        zorder=2,
    )

    if existing_twin is None:
        ax_r.set_yscale("log")
        ax_r.set_ylim(*ylim)
        ax_r.set_ylabel("A' Signal Distribution [Normalized]", color=color, fontsize=9.8)
        ax_r.tick_params(axis="y", which="both", direction="in", colors=color)
        ax_r.spines["right"].set_color(color)
        ax_r.yaxis.set_major_locator(LogLocator(base=10.0, numticks=8))
        ax_r.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=12))
        ax_r.yaxis.set_minor_formatter(NullFormatter())
        # Keep efficiency markers visually on top of the green A' step curve
        ax.set_zorder(ax_r.get_zorder() + 1)
        ax.patch.set_visible(False)

    return ax_r, handle


def safe_asymmetric_yerr(
    val: np.ndarray,
    err_lo: Optional[np.ndarray] = None,
    err_hi: Optional[np.ndarray] = None,
    ci_lo: Optional[np.ndarray] = None,
    ci_hi: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Return non-negative (2, N) error array for Matplotlib errorbar.

    Accepts either direct half-widths (err_lo, err_hi) or interval bounds (ci_lo, ci_hi).
    """
    v = np.asarray(val, dtype=float)
    if err_lo is not None and err_hi is not None:
        lo = np.asarray(err_lo, dtype=float)
        hi = np.asarray(err_hi, dtype=float)
    elif ci_lo is not None and ci_hi is not None:
        lo = v - np.asarray(ci_lo, dtype=float)
        hi = np.asarray(ci_hi, dtype=float) - v
    else:
        lo = np.zeros_like(v)
        hi = np.zeros_like(v)
    lo = np.where(np.isfinite(lo), np.maximum(lo, 0.0), 0.0)
    hi = np.where(np.isfinite(hi), np.maximum(hi, 0.0), 0.0)
    return np.vstack([lo, hi])


def plot_binned_errorbar(
    ax: Axes,
    edges: np.ndarray,
    values: np.ndarray,
    yerr: Any,
    mask: Optional[np.ndarray] = None,
    label: str = "",
    color: str = COLORS["mc"],
    marker: str = "o",
    markersize: float = 3.2,
    markerfacecolor: Optional[str] = None,
    elinewidth: float = 0.85,
    alpha: float = 0.95,
    zorder: int = 5,
    show_xerr: bool = True,
    x_offset: float = 0.0,
):
    """Plot binned points with +/- half-bin horizontal error bars and vertical error bars."""
    edges = np.asarray(edges, dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:]) + x_offset
    half_w = 0.5 * (edges[1:] - edges[:-1])
    vals = np.asarray(values, dtype=float)

    valid = np.isfinite(vals)
    if mask is not None:
        valid = valid & np.asarray(mask, dtype=bool)
    if not np.any(valid):
        return None

    x = centers[valid]
    y = vals[valid]
    xe = half_w[valid] if show_xerr else None

    if yerr is None:
        ye = None
    else:
        ye_arr = np.asarray(yerr, dtype=float)
        if ye_arr.ndim == 2:
            ye = np.maximum(np.where(np.isfinite(ye_arr[:, valid]), ye_arr[:, valid], 0.0), 0.0)
        else:
            ye = np.maximum(np.where(np.isfinite(ye_arr[valid]), ye_arr[valid], 0.0), 0.0)

    mfc = color if markerfacecolor is None else markerfacecolor
    return ax.errorbar(
        x,
        y,
        yerr=ye,
        xerr=xe,
        fmt=marker,
        color=color,
        markerfacecolor=mfc,
        markeredgecolor=color,
        markeredgewidth=0.9,
        markersize=markersize,
        elinewidth=elinewidth,
        capsize=0.0,
        alpha=alpha,
        label=label,
        zorder=zorder,
    )


def choose_eff_ylim(
    series_list: Sequence[Tuple[np.ndarray, Optional[np.ndarray]]],
    default: Tuple[float, float] = (0.80, 1.00),
    min_floor: float = 0.60,
) -> Tuple[float, float]:
    """Keep Pawan's default [0.80, 1.00] unless valid points dip below it; step down in 0.05 increments."""
    lo_default, hi_default = default
    min_val = lo_default
    for vals, mask in series_list:
        arr = np.asarray(vals, dtype=float)
        valid = np.isfinite(arr)
        if mask is not None:
            valid = valid & np.asarray(mask, dtype=bool)
        if np.any(valid):
            min_val = min(min_val, float(np.nanmin(arr[valid])))
    if min_val >= lo_default - 0.005:
        return (lo_default, hi_default)
    # Round down to next 0.05 multiple with clean breathing room
    floored = max(min_floor, float(np.floor((min_val - 0.025) * 20.0) / 20.0))
    return (floored, hi_default)


def choose_ratio_ylim(
    ratio_arrays: Sequence[Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]],
    default: Tuple[float, float] = (0.95, 1.05),
    max_half_span: float = 0.14,
) -> Tuple[float, float, float, float]:
    """Keep Pawan's [0.95, 1.05] unless well-measured points dip outside; return (lo, hi, major, minor)."""
    max_dev = 0.0
    for vals, yerr, mask in ratio_arrays:
        arr = np.asarray(vals, dtype=float)
        valid = np.isfinite(arr) & (arr > 0.0)
        if mask is not None:
            valid = valid & np.asarray(mask, dtype=bool)
        if yerr is not None:
            ye = np.asarray(yerr, dtype=float)
            err_mag = np.max(ye, axis=0) if ye.ndim == 2 else ye
            # Ignore ultra-noisy low-stat bins (err > 0.025) when deciding axis expansion
            valid = valid & np.isfinite(err_mag) & (err_mag <= 0.025)
        if np.any(valid):
            dev = float(np.nanmax(np.abs(arr[valid] - 1.0)))
            max_dev = max(max_dev, dev)

    if max_dev <= 0.044:
        return (default[0], default[1], 0.02, 0.005)
    target = min(max_half_span, max_dev * 1.18)
    for half_span, major, minor in (
        (0.08, 0.04, 0.01),
        (0.10, 0.05, 0.01),
        (0.12, 0.04, 0.01),
        (0.14, 0.05, 0.01),
    ):
        if target <= half_span:
            return (1.0 - half_span, 1.0 + half_span, major, minor)
    return (1.0 - max_half_span, 1.0 + max_half_span, 0.05, 0.01)


def choose_delta_ylim(
    delta_arrays: Sequence[Tuple[np.ndarray, Optional[np.ndarray]]],
    min_half_span: float = 1.0,
    max_half_span: float = 35.0,
    err_arrays: Optional[Sequence[np.ndarray]] = None,
) -> Tuple[float, float, float, float]:
    """Choose clean y-limits (in percentage points) and (major, minor) ticks for paired delta plots.

    Preserves narrow deep features (using true min/max of well-measured bins rather than percentiles)
    and supports clean asymmetric ranges when shifts are strongly one-sided.
    """
    vmin = 0.0
    vmax = 0.0
    for idx, (vals, mask) in enumerate(delta_arrays):
        arr = np.asarray(vals, dtype=float)
        valid = np.isfinite(arr)
        if mask is not None:
            valid = valid & np.asarray(mask, dtype=bool)
        if err_arrays is not None and idx < len(err_arrays):
            e_arr = np.asarray(err_arrays[idx], dtype=float)
            valid = valid & np.isfinite(e_arr) & (e_arr <= 2.5)
        if np.any(valid):
            vmin = min(vmin, float(np.nanmin(arr[valid])))
            vmax = max(vmax, float(np.nanmax(arr[valid])))

    peak = max(abs(vmin), abs(vmax))
    # Strongly one-sided deep feature (e.g., -25 pp crater with +2 pp plateau): use asymmetric range
    if peak > 8.0 and (abs(vmin) > 2.2 * max(1.0, abs(vmax)) or abs(vmax) > 2.2 * max(1.0, abs(vmin))):
        lo_raw = min(-2.0, vmin * 1.14 - 0.8)
        hi_raw = max(2.0, vmax * 1.25 + 1.5)
        span = hi_raw - lo_raw
        major = 5.0 if span <= 26.0 else 10.0
        minor = 1.0 if major == 5.0 else 2.0
        lo = max(-max_half_span, float(np.floor(lo_raw / 2.0) * 2.0))
        hi = min(max_half_span, float(np.ceil(hi_raw / 2.0) * 2.0))
        # Ensure room for top-right legend above zero line when feature is negative
        if abs(vmin) > abs(vmax):
            hi = max(hi, 6.0 if span <= 26.0 else 8.0)
        return (lo, hi, major, minor)

    target = max(min_half_span, min(max_half_span, peak * 1.22))
    for span, major, minor in (
        (0.5, 0.25, 0.05),
        (1.0, 0.5, 0.1),
        (1.5, 0.5, 0.1),
        (2.0, 1.0, 0.2),
        (3.0, 1.0, 0.25),
        (5.0, 2.0, 0.5),
        (8.0, 2.0, 0.5),
        (10.0, 5.0, 1.0),
        (15.0, 5.0, 1.0),
        (20.0, 5.0, 1.0),
        (30.0, 10.0, 2.0),
        (35.0, 10.0, 2.0),
    ):
        if target <= span:
            return (-span, span, major, minor)
    return (-max_half_span, max_half_span, 10.0, 2.0)


def combine_legend(
    ax: Axes,
    extra_handles: Optional[Sequence[Any]] = None,
    loc: str = "lower right",
    ncol: int = 1,
    fontsize: float = 8.2,
    bbox_to_anchor: Optional[Tuple[float, float]] = None,
) -> None:
    """Create a legible, non-overlapping legend merging primary and twin-axis handles."""
    handles, labels = ax.get_legend_handles_labels()
    if extra_handles:
        for h in extra_handles:
            if h is not None:
                lbl = h.get_label()
                if lbl and not lbl.startswith("_"):
                    handles.append(h)
                    labels.append(lbl)
    if not handles:
        return
    kwargs: Dict[str, Any] = {
        "loc": loc,
        "ncol": ncol,
        "fontsize": fontsize,
        "frameon": True,
        "facecolor": "white",
        "edgecolor": "#dddddd",
        "framealpha": 0.90,
        "borderpad": 0.4,
        "handletextpad": 0.5,
        "labelspacing": 0.32,
    }
    if bbox_to_anchor is not None:
        kwargs["bbox_to_anchor"] = bbox_to_anchor
    leg = ax.legend(handles, labels, **kwargs)
    leg.set_zorder(20)


def format_count_compact(n: int) -> str:
    """Format large event counts cleanly (e.g., 3,487,558 -> '3.49M' or exact comma-formatted)."""
    return f"{int(n):,}"


def save_figure_outputs(
    fig: Figure,
    out_dir: Path,
    stem: str,
    dpi: int = 300,
) -> Tuple[Path, Path]:
    """Save figure as both 300 DPI PNG and vector PDF, closing the figure afterwards."""
    png_dir = out_dir / "png"
    pdf_dir = out_dir / "pdf"
    png_dir.mkdir(parents=True, exist_ok=True)
    pdf_dir.mkdir(parents=True, exist_ok=True)
    png_path = png_dir / f"{stem}.png"
    pdf_path = pdf_dir / f"{stem}.pdf"
    fig.savefig(png_path, dpi=dpi)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path, pdf_path
