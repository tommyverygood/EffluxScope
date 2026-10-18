#!/usr/bin/env python3
"""Render the frozen ABCB1 structural-shift/OOD audit as a four-panel figure.

The plotting code consumes only the locked structural assignments and the
completed metric/bootstrap outputs.  It does not fit, recalibrate, or modify
any model or analysis result.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"
OUT_DIR = PROJECT / "outputs" / "20260807_abcb1_ood_structural_shift_v1" / "figures"

ANALYSIS_ROWS = PROCESSED / "abcb1_ood_analysis_rows_v1.csv"
STRATUM_METRICS = PROCESSED / "abcb1_ood_stratum_metrics_v1.csv"
BOOTSTRAP_SUMMARY = PROCESSED / "abcb1_ood_bootstrap_summary_v1.csv"
STRUCTURE_LOCK = PROCESSED / "abcb1_ood_structure_assignment_lock_v1.json"

PNG_OUT = OUT_DIR / "Figure2_abcb1_ood_structural_shift.png"
PDF_OUT = OUT_DIR / "Figure2_abcb1_ood_structural_shift.pdf"

PRIMARY_BOOTSTRAP_SCHEME = "murcko_cluster"

# Color-blind-friendly, restrained journal palette.
NAVY = "#2B5C8A"
ORANGE = "#D9852F"
DARK_GREY = "#4D4D4D"
MID_GREY = "#8A8A8A"
LIGHT_GREY = "#D9D9D9"
VERY_LIGHT_GREY = "#F2F2F2"
BLACK = "#1A1A1A"

QUINTILES = [
    ("external_Q1_nearest", "Q1"),
    ("external_Q2", "Q2"),
    ("external_Q3", "Q3"),
    ("external_Q4", "Q4"),
    ("external_Q5_farthest", "Q5"),
]

PRESSURE_SETS = [
    ("reference_extreme_OOD", "Extreme OOD"),
    ("strict_LoHi", "Strict LoHi"),
    ("dual_OOD", "Dual OOD"),
    ("kmeans_OOD20", "KMeans OOD20"),
    ("external_Q5_farthest", "Q5 farthest"),
]


def configure_style() -> None:
    """Set a compact, editable publication style."""
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.0,
            "axes.labelsize": 8.0,
            "axes.titlesize": 8.5,
            "axes.titleweight": "normal",
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 7.0,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "lines.linewidth": 1.15,
            "lines.markersize": 4.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def require_inputs() -> None:
    missing = [path for path in [ANALYSIS_ROWS, STRATUM_METRICS, BOOTSTRAP_SUMMARY, STRUCTURE_LOCK] if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Required frozen OOD outputs are missing: " + ", ".join(str(path) for path in missing)
        )


def metric_row(metrics: pd.DataFrame, subset: str, layer: str) -> pd.Series:
    match = metrics.loc[metrics["subset_id"].eq(subset) & metrics["layer"].eq(layer)]
    if len(match) != 1:
        raise RuntimeError(f"Expected one metric row for {subset}/{layer}; found {len(match)}")
    return match.iloc[0]


def bootstrap_row(bootstrap: pd.DataFrame, subset: str, layer: str, metric: str) -> pd.Series:
    match = bootstrap.loc[
        bootstrap["subset_id"].eq(subset)
        & bootstrap["layer"].eq(layer)
        & bootstrap["metric"].eq(metric)
        & bootstrap["bootstrap_scheme"].eq(PRIMARY_BOOTSTRAP_SCHEME)
    ]
    if len(match) != 1:
        raise RuntimeError(
            f"Expected one {PRIMARY_BOOTSTRAP_SCHEME} bootstrap row for "
            f"{subset}/{layer}/{metric}; found {len(match)}"
        )
    return match.iloc[0]


def finite_ci(row: pd.Series) -> tuple[float, float, float]:
    point = float(row["point_estimate"])
    low = float(row["ci95_low"])
    high = float(row["ci95_high"])
    if not all(math.isfinite(value) for value in [point, low, high]):
        return point, math.nan, math.nan
    return point, low, high


def asymmetric_errors(points: np.ndarray, lows: np.ndarray, highs: np.ndarray) -> np.ndarray:
    """Return Matplotlib-compatible errors while preserving missing intervals."""
    lower = np.where(np.isfinite(lows), np.maximum(points - lows, 0.0), 0.0)
    upper = np.where(np.isfinite(highs), np.maximum(highs - points, 0.0), 0.0)
    return np.vstack([lower, upper])


def clean_axis(ax: plt.Axes, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(BLACK)
    ax.spines["bottom"].set_color(BLACK)
    ax.tick_params(colors=BLACK)
    if grid_axis is not None:
        ax.grid(axis=grid_axis, color=LIGHT_GREY, linewidth=0.45, alpha=0.65, zorder=0)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=10.5,
        fontweight="bold",
        ha="left",
        va="top",
        color=BLACK,
    )


def plot_distance_distribution(ax: plt.Axes, rows: pd.DataFrame, lock: dict) -> None:
    distance = rows["joint_model_d5"].to_numpy(dtype=float)
    thresholds = lock["reference_d5_thresholds"]["joint_model_connectivity_unique"]
    q50, q80, q95 = (float(thresholds[key]) for key in ["q50", "q80", "q95"])

    left = min(0.52, float(np.nanmin(distance)) - 0.01)
    right = max(0.88, float(np.nanmax(distance)) + 0.01)
    bounds = [left, q50, q80, q95, right]
    fills = [VERY_LIGHT_GREY, "#EAF1F7", "#F8EFE4", "#F3E3E3"]
    for index in range(4):
        ax.axvspan(bounds[index], bounds[index + 1], color=fills[index], ec="none", zorder=0)

    bins = np.linspace(left, right, 27)
    counts, _, _ = ax.hist(
        distance,
        bins=bins,
        color=NAVY,
        alpha=0.78,
        edgecolor="white",
        linewidth=0.35,
        zorder=2,
    )
    for value, label in [(q50, "Q50"), (q80, "Q80"), (q95, "Q95")]:
        ax.axvline(value, color=DARK_GREY, linestyle=(0, (2.2, 2.0)), linewidth=0.8, zorder=3)
        ax.text(
            value,
            max(counts) * 1.025,
            label,
            rotation=90,
            ha="right",
            va="bottom",
            fontsize=6.2,
            color=DARK_GREY,
        )

    observed = rows["reference_adaptive_stratum"].value_counts().to_dict()
    ax.text(
        0.035,
        0.78,
        "Strata counts\n"
        f"ID-like {int(observed.get('ID_like', 0)):,};  transition {int(observed.get('transition', 0)):,}\n"
        f"OOD {int(observed.get('OOD', 0)):,};  extreme OOD {int(observed.get('extreme_OOD', 0)):,}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.5,
        linespacing=1.25,
        color=BLACK,
    )

    ax.set_xlim(left, right)
    ax.set_ylim(0, max(counts) * 1.18)
    ax.set_xlabel(r"Structural distance, $D_5=1-$mean top-5 Tanimoto similarity")
    ax.set_ylabel("External compounds")
    ax.set_title("External distance distribution; reference-LOO strata", loc="left", pad=6)
    clean_axis(ax)


def plot_quintile_bacc(ax: plt.Axes, metrics: pd.DataFrame, bootstrap: pd.DataFrame) -> None:
    x = np.arange(1, len(QUINTILES) + 1, dtype=float)
    values: dict[str, list[float]] = {"full": [], "full_low": [], "full_high": [], "accepted": [], "accepted_low": [], "accepted_high": []}
    tick_labels: list[str] = []
    for subset, short in QUINTILES:
        full = bootstrap_row(bootstrap, subset, "full_forced_binary", "balanced_accuracy")
        accepted = bootstrap_row(bootstrap, subset, "selective_accepted", "balanced_accuracy")
        full_point, full_low, full_high = finite_ci(full)
        acc_point, acc_low, acc_high = finite_ci(accepted)
        values["full"].append(full_point)
        values["full_low"].append(full_low)
        values["full_high"].append(full_high)
        values["accepted"].append(acc_point)
        values["accepted_low"].append(acc_low)
        values["accepted_high"].append(acc_high)
        row = metric_row(metrics, subset, "selective_accepted")
        tick_labels.append(f"{short}\n{int(row['n_total'])}/{int(row['accepted_n'])}")

    full_values = np.asarray(values["full"])
    accepted_values = np.asarray(values["accepted"])
    full_low = np.asarray(values["full_low"])
    full_high = np.asarray(values["full_high"])
    accepted_low = np.asarray(values["accepted_low"])
    accepted_high = np.asarray(values["accepted_high"])

    ax.axhline(0.5, color=MID_GREY, linestyle=(0, (3, 2.5)), linewidth=0.75, zorder=1)
    ax.errorbar(
        x - 0.055,
        full_values,
        yerr=asymmetric_errors(full_values, full_low, full_high),
        color=DARK_GREY,
        marker="o",
        markerfacecolor="white",
        markeredgecolor=DARK_GREY,
        capsize=2.0,
        elinewidth=0.8,
        label="Forced binary (100%)",
        zorder=3,
    )
    ax.errorbar(
        x + 0.055,
        accepted_values,
        yerr=asymmetric_errors(accepted_values, accepted_low, accepted_high),
        color=NAVY,
        marker="o",
        markerfacecolor=NAVY,
        markeredgecolor=NAVY,
        capsize=2.0,
        elinewidth=0.8,
        label="v0.4 accepted",
        zorder=4,
    )
    ax.set_xticks(x, tick_labels)
    ax.set_xlim(0.65, 5.35)
    all_bounds = np.r_[full_low, accepted_low, full_high, accepted_high]
    finite_bounds = all_bounds[np.isfinite(all_bounds)]
    lower = min(0.46, float(np.min(finite_bounds)) - 0.025) if finite_bounds.size else 0.46
    upper = max(0.84, float(np.max(finite_bounds)) + 0.025) if finite_bounds.size else 0.84
    ax.set_ylim(max(0.35, lower), min(1.0, upper))
    ax.set_xlabel("External $D_5$ quintile (n total/accepted)")
    ax.set_ylabel("Balanced accuracy")
    ax.set_title("Performance across external distance quintiles", loc="left", pad=6)
    ax.legend(loc="lower left", frameon=False, handlelength=2.0, borderaxespad=0.2)
    clean_axis(ax, grid_axis="y")


def plot_quintile_coverage(ax: plt.Axes, metrics: pd.DataFrame, bootstrap: pd.DataFrame) -> None:
    x = np.arange(1, len(QUINTILES) + 1, dtype=float)
    points, lows, highs = [], [], []
    accepted_counts: list[int] = []
    for subset, _ in QUINTILES:
        row = bootstrap_row(bootstrap, subset, "selective_framework", "coverage")
        point, low, high = finite_ci(row)
        points.append(point * 100.0)
        lows.append(low * 100.0)
        highs.append(high * 100.0)
        accepted_counts.append(int(metric_row(metrics, subset, "selective_accepted")["accepted_n"]))

    points_array = np.asarray(points)
    lows_array = np.asarray(lows)
    highs_array = np.asarray(highs)
    ax.errorbar(
        x,
        points_array,
        yerr=asymmetric_errors(points_array, lows_array, highs_array),
        color=NAVY,
        marker="o",
        markerfacecolor=NAVY,
        markeredgecolor=NAVY,
        capsize=2.0,
        elinewidth=0.8,
        zorder=3,
    )
    for xpos, point, accepted in zip(x, points_array, accepted_counts):
        ax.text(
            xpos,
            point + 1.25,
            f"{point:.1f}%\n(n={accepted})",
            ha="center",
            va="bottom",
            fontsize=6.4,
            color=BLACK,
        )
    ax.set_xticks(x, [short for _, short in QUINTILES])
    ax.set_xlim(0.65, 5.35)
    finite_highs = highs_array[np.isfinite(highs_array)]
    top = max(27.0, float(np.max(finite_highs)) + 5.0) if finite_highs.size else 27.0
    ax.set_ylim(0.0, top)
    ax.set_xlabel("External $D_5$ quintile (nearest → farthest)")
    ax.set_ylabel("v0.4 coverage (%)")
    ax.set_title("Selective coverage across chemical-space shift", loc="left", pad=6)
    clean_axis(ax, grid_axis="y")


def plot_pressure_sets(ax: plt.Axes, metrics: pd.DataFrame, bootstrap: pd.DataFrame) -> None:
    base_y = np.arange(len(PRESSURE_SETS) - 1, -1, -1, dtype=float)
    offset = 0.12
    accepted_points, accepted_lows, accepted_highs = [], [], []
    rejected_points, rejected_lows, rejected_highs = [], [], []
    tick_labels: list[str] = []
    accepted_ns: list[int] = []

    for subset, label in PRESSURE_SETS:
        accepted = bootstrap_row(bootstrap, subset, "selective_accepted", "balanced_accuracy")
        rejected = bootstrap_row(bootstrap, subset, "forced_rejected", "balanced_accuracy")
        acc_point, acc_low, acc_high = finite_ci(accepted)
        rej_point, rej_low, rej_high = finite_ci(rejected)
        accepted_points.append(acc_point)
        accepted_lows.append(acc_low)
        accepted_highs.append(acc_high)
        rejected_points.append(rej_point)
        rejected_lows.append(rej_low)
        rejected_highs.append(rej_high)
        metric = metric_row(metrics, subset, "selective_accepted")
        accepted_n = int(metric["accepted_n"])
        rejected_n = int(metric["rejected_n"])
        accepted_ns.append(accepted_n)
        tick_labels.append(f"{label}\n{accepted_n}/{rejected_n}")

    accepted_points_array = np.asarray(accepted_points)
    accepted_lows_array = np.asarray(accepted_lows)
    accepted_highs_array = np.asarray(accepted_highs)
    rejected_points_array = np.asarray(rejected_points)
    rejected_lows_array = np.asarray(rejected_lows)
    rejected_highs_array = np.asarray(rejected_highs)

    ax.axvline(0.5, color=MID_GREY, linestyle=(0, (3, 2.5)), linewidth=0.75, zorder=1)
    for ypos, left, right in zip(base_y, rejected_points_array, accepted_points_array):
        ax.plot([left, right], [ypos, ypos], color=LIGHT_GREY, linewidth=0.8, zorder=1)

    ax.errorbar(
        rejected_points_array,
        base_y - offset,
        xerr=asymmetric_errors(rejected_points_array, rejected_lows_array, rejected_highs_array),
        fmt="o",
        color=DARK_GREY,
        markerfacecolor="white",
        markeredgecolor=DARK_GREY,
        capsize=1.8,
        elinewidth=0.8,
        label="Forced rejected",
        zorder=3,
    )
    ax.errorbar(
        accepted_points_array,
        base_y + offset,
        xerr=asymmetric_errors(accepted_points_array, accepted_lows_array, accepted_highs_array),
        fmt="o",
        color=NAVY,
        markerfacecolor=NAVY,
        markeredgecolor=NAVY,
        capsize=1.8,
        elinewidth=0.8,
        label="v0.4 accepted",
        zorder=4,
    )

    # The extreme-OOD accepted subset has only five observations.  Keep it
    # visible as a descriptive estimate, but explicitly flag the sample size.
    if accepted_ns[0] < 15:
        ax.scatter(
            [accepted_points_array[0]],
            [base_y[0] + offset],
            s=25,
            facecolors="white",
            edgecolors=NAVY,
            linewidths=0.9,
            zorder=5,
        )
        ax.text(
            min(0.98, accepted_points_array[0] + 0.018),
            base_y[0] + offset,
            "descriptive",
            va="center",
            ha="left",
            fontsize=6.0,
            color=NAVY,
        )

    ax.set_yticks(base_y, tick_labels)
    ax.set_xlabel("Balanced accuracy")
    ax.set_ylabel("OOD pressure set\n(n accepted/rejected)")
    all_bounds = np.r_[accepted_lows_array, accepted_highs_array, rejected_lows_array, rejected_highs_array]
    finite_bounds = all_bounds[np.isfinite(all_bounds)]
    lower = min(0.45, float(np.min(finite_bounds)) - 0.03) if finite_bounds.size else 0.45
    upper = max(0.84, float(np.max(finite_bounds)) + 0.08) if finite_bounds.size else 0.84
    ax.set_xlim(max(0.30, lower), min(1.02, upper))
    # Reserve a narrow empty band below the last pressure set so the legend
    # does not conceal the Q5 estimates or confidence intervals.
    ax.set_ylim(-1.18, len(PRESSURE_SETS) - 0.45)
    ax.set_title("Accepted versus forced-rejected OOD compounds", loc="left", pad=6)
    ax.legend(loc="lower right", frameon=False, handlelength=1.8, borderaxespad=0.2)
    clean_axis(ax, grid_axis="x")


def main() -> None:
    require_inputs()
    configure_style()
    rows = pd.read_csv(ANALYSIS_ROWS, low_memory=False)
    metrics = pd.read_csv(STRATUM_METRICS, low_memory=False)
    bootstrap = pd.read_csv(BOOTSTRAP_SUMMARY, low_memory=False)
    lock = json.loads(STRUCTURE_LOCK.read_text(encoding="utf-8"))

    if len(rows) != 2624:
        raise RuntimeError(f"Expected 2,624 identity-clean external rows; found {len(rows)}")
    if lock.get("status") != "structure_only_assignments_frozen_before_new_metric_generation":
        raise RuntimeError("Structure-only assignment lock is not valid")

    # 180 mm width approximates a two-column journal figure.
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 6.28), constrained_layout=False)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.965, bottom=0.10, wspace=0.37, hspace=0.44)

    plot_distance_distribution(axes[0, 0], rows, lock)
    plot_quintile_bacc(axes[0, 1], metrics, bootstrap)
    plot_quintile_coverage(axes[1, 0], metrics, bootstrap)
    plot_pressure_sets(axes[1, 1], metrics, bootstrap)
    for ax, label in zip(axes.flat, ["a", "b", "c", "d"]):
        panel_label(ax, label)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(PDF_OUT, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)

    print(
        json.dumps(
            {
                "status": "complete",
                "bootstrap_scheme": PRIMARY_BOOTSTRAP_SCHEME,
                "rows_n": int(len(rows)),
                "png": str(PNG_OUT),
                "pdf": str(PDF_OUT),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
