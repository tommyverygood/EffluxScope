#!/usr/bin/env python3
"""Create supplementary figures for deep ABCB1 calibration and error analysis."""

from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"
FIGURE_DIR = PROJECT / "outputs" / "20260805_abcb1_external_calibration_error_v2" / "figures"
MANIFEST = PROCESSED / "abcb1_external_calibration_error_figure_manifest_v2.csv"

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#CC79A7"
GRAY = "#6B7280"
LIGHT = "#D1D5DB"
DARK = "#1F2937"
RED = "#B91C1C"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color=LIGHT, linewidth=0.6, alpha=0.65)
    ax.set_axisbelow(True)


def panel(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.05, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def save(fig: plt.Figure, stem: str, title: str, sources: str, rows: list[dict[str, object]]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "svg", "pdf"]:
        path = FIGURE_DIR / f"{stem}.{extension}"
        fig.savefig(path, dpi=300 if extension == "png" else None)
        rows.append({
            "figure": stem.split("_")[0].replace("FigureS", "Figure S"),
            "title": title,
            "format": extension,
            "path": str(path.relative_to(PROJECT)),
            "sha256": sha256(path),
            "size_bytes": path.stat().st_size,
            "source_tables": sources,
        })
    plt.close(fig)


def figure_s5(rows: list[dict[str, object]]) -> None:
    bins = pd.read_csv(PROCESSED / "abcb1_external_calibration_bins_v2.csv")
    bins = bins[(bins["population"] == "identity_clean") & (bins["strategy"] == "equal_width")]
    metrics = pd.read_csv(PROCESSED / "abcb1_external_calibration_metrics_v2.csv")
    metrics = metrics[metrics["population"] == "identity_clean"].set_index("method")
    counts = pd.read_csv(PROCESSED / "abcb1_external_error_counts_v2.csv")
    counts = counts[counts["population"] == "identity_clean"]
    strata = pd.read_csv(PROCESSED / "abcb1_external_error_strata_v2.csv")
    bands = strata[(strata["population"] == "identity_clean") & (strata["stratum_type"] == "probability_band")].copy()

    fig, axes_grid = plt.subplots(2, 2, figsize=(11.8, 8.2))
    axes = axes_grid.ravel()
    method_specs = [
        ("raw_frozen_score", BLUE, "Raw frozen score"),
        ("cross_fitted_platt", GREEN, "Cross-fitted Platt"),
        ("cross_fitted_isotonic", ORANGE, "Cross-fitted isotonic"),
    ]
    axes[0].plot([0, 1], [0, 1], color=GRAY, linestyle="--", linewidth=1, label="Ideal")
    for method, color, label in method_specs:
        frame = bins[bins["method"] == method].sort_values("mean_predicted_probability")
        axes[0].errorbar(
            frame["mean_predicted_probability"], frame["observed_positive_fraction"],
            yerr=np.vstack([
                frame["observed_positive_fraction"] - frame["observed_positive_wilson95_low"],
                frame["observed_positive_wilson95_high"] - frame["observed_positive_fraction"],
            ]),
            color=color, marker="o", markersize=4, linewidth=1.5, capsize=2, label=label,
        )
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)
    axes[0].set_xlabel("Mean predicted probability")
    axes[0].set_ylabel("Observed ER ≥ 2 fraction")
    axes[0].set_title("Identity-clean calibration curves")
    axes[0].legend(frameon=False, loc="upper left")
    clean_axis(axes[0])
    panel(axes[0], "A")

    metric_order = ["brier_score", "ece_equal_width_10", "log_loss"]
    metric_names = ["Brier", "ECE", "Log loss"]
    x = np.arange(len(metric_order))
    width = 0.24
    for method_index, (method, color, label) in enumerate(method_specs):
        values = metrics.loc[method, metric_order].to_numpy(dtype=float)
        positions = x + (method_index - 1) * width
        axes[1].bar(positions, values, width, color=color, alpha=0.82, label=label)
        for position, value in zip(positions, values):
            axes[1].text(position, value + 0.012, f"{value:.3f}", ha="center", va="bottom", fontsize=7.5, rotation=90 if value < 0.05 else 0)
    axes[1].set_xticks(x, metric_names)
    axes[1].set_ylim(0, 0.78)
    axes[1].set_ylabel("Metric value (lower is better)")
    axes[1].set_title("Same-source cross-fitted calibration diagnostic")
    axes[1].legend(frameon=False, loc="upper left")
    clean_axis(axes[1])
    panel(axes[1], "B")

    outcome_order = ["true_positive", "true_negative", "false_positive", "false_negative"]
    outcome_names = ["True positive", "True negative", "False positive", "False negative"]
    outcome_colors = [GREEN, BLUE, RED, ORANGE]
    frame = counts.set_index("error_outcome").loc[outcome_order]
    axes[2].barh(np.arange(4), frame["n"], color=outcome_colors, alpha=0.85)
    axes[2].set_yticks(np.arange(4), outcome_names)
    axes[2].invert_yaxis()
    for index, value in enumerate(frame["n"]):
        axes[2].text(value + 15, index, f"{int(value):,}", va="center", fontsize=8)
    axes[2].set_xlim(0, frame["n"].max() * 1.18)
    axes[2].set_xlabel("Structures")
    axes[2].set_title("Forced-binary error taxonomy")
    clean_axis(axes[2])
    panel(axes[2], "C")

    bands["order"] = bands["stratum_level"].str.extract(r"\[(\d+\.\d+)")[0].astype(float)
    bands = bands.sort_values("order")
    x = np.arange(len(bands))
    y = bands["error_rate"].to_numpy()
    lower = y - bands["error_rate_wilson95_low"].to_numpy()
    upper = bands["error_rate_wilson95_high"].to_numpy() - y
    axes[3].errorbar(x, y, yerr=np.vstack([lower, upper]), color=PURPLE, marker="o", linewidth=1.8, capsize=3)
    axes[3].axhline(0.5, color=GRAY, linestyle="--", linewidth=1)
    for index, row in bands.reset_index(drop=True).iterrows():
        axes[3].text(index, row["error_rate"] + upper[index] + 0.025, f"{row['error_rate']:.1%}\nn={int(row['n'])}", ha="center", fontsize=7.5)
    axes[3].set_xticks(x, bands["stratum_level"], rotation=20)
    axes[3].set_ylim(0, 0.70)
    axes[3].set_xlabel("Frozen score band")
    axes[3].set_ylabel("Forced-decision error rate")
    axes[3].set_title("Error is highest just above the 0.5 threshold")
    clean_axis(axes[3])
    panel(axes[3], "D")

    fig.suptitle("Deep external calibration and forced-decision error analysis", fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    save(
        fig,
        "FigureS5_external_calibration_error_v2",
        "Deep external calibration and forced-decision error analysis",
        "abcb1_external_calibration_bins_v2.csv|abcb1_external_calibration_metrics_v2.csv|abcb1_external_error_counts_v2.csv|abcb1_external_error_strata_v2.csv",
        rows,
    )


def effect_panel(ax: plt.Axes, features: pd.DataFrame, comparison: str, title: str, panel_label: str) -> None:
    frame = features[(features["population"] == "identity_clean") & (features["comparison"] == comparison) & (features["feature"] != "log10_efflux_ratio")].copy()
    frame["absolute_delta"] = frame["cliffs_delta"].abs()
    frame = frame.sort_values("absolute_delta", ascending=False).head(10).sort_values("cliffs_delta")
    display = {
        "raw_probability": "Raw score",
        "bootstrap_p10": "Bootstrap p10",
        "bootstrap_p90": "Bootstrap p90",
        "bootstrap_interval_width": "Bootstrap width",
        "direct_probability": "Direct score",
        "primary_max_tanimoto": "Primary similarity",
        "direct_max_tanimoto": "Direct similarity",
        "primary_local_positive_fraction": "Primary local positive",
        "direct_local_positive_fraction": "Direct local positive",
        "primary_local_conflict_index": "Primary local conflict",
        "direct_local_conflict_index": "Direct local conflict",
        "physchem_max_abs_robust_z": "Physchem max |z|",
        "molecular_weight": "Molecular weight",
        "xlogp": "XlogP",
        "tpsa": "TPSA",
        "hbond_donors": "H-bond donors",
        "hbond_acceptors": "H-bond acceptors",
        "rotatable_bonds": "Rotatable bonds",
        "formal_charge": "Formal charge",
        "heavy_atoms": "Heavy atoms",
    }
    colors = np.where(frame["cliffs_delta"] >= 0, ORANGE, BLUE)
    ax.barh(np.arange(len(frame)), frame["cliffs_delta"], color=colors, alpha=0.82)
    ax.set_yticks(np.arange(len(frame)), frame["feature"].map(display).fillna(frame["feature"]))
    ax.axvline(0, color=DARK, linewidth=0.8)
    ax.set_xlim(-0.58, 0.58)
    ax.set_xlabel("Cliff delta (error group minus correct group)")
    ax.set_title(title)
    clean_axis(ax)
    panel(ax, panel_label)


def figure_s6(rows: list[dict[str, object]]) -> None:
    features = pd.read_csv(PROCESSED / "abcb1_external_error_feature_comparisons_v2.csv")
    strata = pd.read_csv(PROCESSED / "abcb1_external_error_strata_v2.csv")
    hotspots = pd.read_csv(PROCESSED / "abcb1_external_error_scaffold_hotspots_v2.csv")

    fig, axes_grid = plt.subplots(2, 2, figsize=(12.2, 8.5))
    axes = axes_grid.ravel()
    effect_panel(axes[0], features, "false_positive_vs_true_positive", "False positives versus true positives", "A")
    effect_panel(axes[1], features, "false_negative_vs_true_negative", "False negatives versus true negatives", "B")

    selected_levels = [
        ("selection_layer", "v0_4_accepted", "v0.4 accepted"),
        ("selection_layer", "v0_4_rejected", "v0.4 rejected"),
        ("support_combination", "both_supported", "Both supported"),
        ("support_combination", "direct_only", "Direct only"),
        ("support_combination", "primary_only", "Primary only"),
        ("support_combination", "neither_supported", "Neither supported"),
    ]
    selected_rows = []
    identity = strata[strata["population"] == "identity_clean"]
    for stratum_type, level, label in selected_levels:
        row = identity[(identity["stratum_type"] == stratum_type) & (identity["stratum_level"] == level)].iloc[0].copy()
        row["display"] = label
        selected_rows.append(row)
    frame = pd.DataFrame(selected_rows)
    y = np.arange(len(frame))
    values = frame["error_rate"].to_numpy(dtype=float)
    lower = values - frame["error_rate_wilson95_low"].to_numpy(dtype=float)
    upper = frame["error_rate_wilson95_high"].to_numpy(dtype=float) - values
    colors = [GREEN, ORANGE, PURPLE, BLUE, ORANGE, RED]
    axes[2].barh(y, values, color=colors, alpha=0.82, xerr=np.vstack([lower, upper]), error_kw={"elinewidth": 0.8, "capsize": 2})
    axes[2].set_yticks(y, frame["display"])
    axes[2].invert_yaxis()
    for position, (_, row) in enumerate(frame.iterrows()):
        axes[2].text(row["error_rate"] + upper[position] + 0.012, position, f"{row['error_rate']:.1%} | n={int(row['n'])}", va="center", fontsize=7.5)
    axes[2].set_xlim(0, 0.62)
    axes[2].set_xlabel("Forced-decision error rate (Wilson 95% CI)")
    axes[2].set_title("Error localization by selection and support")
    clean_axis(axes[2])
    panel(axes[2], "C")

    frame = hotspots[(hotspots["population"] == "identity_clean") & (hotspots["n"] >= 4)].sort_values("shrunken_error_rate", ascending=False).head(10).copy()
    frame = frame.sort_values("shrunken_error_rate")
    labels = [text if len(text) <= 32 else text[:29] + "…" for text in frame["bootstrap_cluster"].astype(str)]
    y = np.arange(len(frame))
    axes[3].barh(y, frame["shrunken_error_rate"], color=GRAY, alpha=0.75, label="Shrunken error rate")
    axes[3].scatter(frame["raw_error_rate"], y, color=RED, s=24, zorder=3, label="Raw error rate")
    axes[3].set_yticks(y, labels)
    for index, (_, row) in enumerate(frame.iterrows()):
        axes[3].text(max(row["raw_error_rate"], row["shrunken_error_rate"]) + 0.015, index, f"n={int(row['n'])}", va="center", fontsize=7.2)
    axes[3].set_xlim(0, 0.90)
    axes[3].set_xlabel("Error rate")
    axes[3].set_title("Murcko-scaffold error hotspots (descriptive)")
    axes[3].legend(frameon=False, loc="lower right")
    clean_axis(axes[3])
    panel(axes[3], "D")

    fig.suptitle("External error phenotypes and concentration hotspots", fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    save(
        fig,
        "FigureS6_external_error_phenotypes_v2",
        "External error phenotypes and concentration hotspots",
        "abcb1_external_error_feature_comparisons_v2.csv|abcb1_external_error_strata_v2.csv|abcb1_external_error_scaffold_hotspots_v2.csv",
        rows,
    )


def main() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    rows: list[dict[str, object]] = []
    figure_s5(rows)
    figure_s6(rows)
    pd.DataFrame(rows).to_csv(MANIFEST, index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
