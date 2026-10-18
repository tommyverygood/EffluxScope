#!/usr/bin/env python3
"""Plot public ABCB1 models across frozen project-reference distance strata.

This is a descriptive projection onto this project's ABCB1 reference-relative
D5 quintiles. It does not define OOD relative to any public model's training
chemistry and does not refit, recalibrate, or retune a model.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"
INPUT = PROCESSED / "abcb1_public_predictor_ood_stratum_metrics_v1.csv"
OUT_DIR = PROJECT / "outputs" / "20260807_abcb1_ood_structural_shift_v1" / "figures"
PNG_OUT = OUT_DIR / "Figure3_public_models_across_reference_relative_distance.png"
PDF_OUT = OUT_DIR / "Figure3_public_models_across_reference_relative_distance.pdf"

MODEL_STYLE = {
    "astrazeneca_gnn_mtl_mdck_er_2025": {
        "label": "AZ MDCK ER",
        "color": "#2B5C8A",
        "marker": "o",
    },
    "astrazeneca_gnn_mtl_nih_mdck_er_2025": {
        "label": "AZ NIH-MDCK ER",
        "color": "#D9852F",
        "marker": "s",
    },
    "deeppk_pgp_substrate_2026_08_05": {
        "label": "Deep-PK substrate",
        "color": "#3C8D76",
        "marker": "^",
    },
}

QUINTILES = ["Q1_nearest", "Q2", "Q3", "Q4", "Q5_farthest"]
TICK_LABELS = ["Q1", "Q2", "Q3", "Q4", "Q5"]
METRICS = [
    ("balanced_accuracy", "Balanced accuracy", (0.48, 0.80), 0.5),
    ("roc_auc", "ROC AUC", (0.64, 0.87), 0.5),
    ("spearman_rho", "Spearman $\\rho$", (0.30, 0.69), 0.0),
]


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7.4,
            "axes.labelsize": 7.6,
            "axes.titlesize": 8.2,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 7.0,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "lines.linewidth": 1.25,
            "lines.markersize": 4.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.45, alpha=0.70, zorder=0)
    ax.set_axisbelow(True)


def load_rows() -> pd.DataFrame:
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)
    data = pd.read_csv(INPUT)
    rows = data.loc[
        data["comparator_id"].isin(MODEL_STYLE)
        & data["evaluation_layer"].eq("full_coverage")
        & data["stratum_type"].eq("external_d5_quintile")
        & data["stratum_value"].isin(QUINTILES)
        & data["substrate_ranking_eligible"].eq(True)
    ].copy()
    expected = len(MODEL_STYLE) * len(QUINTILES)
    if len(rows) != expected:
        raise RuntimeError(f"Expected {expected} model/quintile rows; found {len(rows)}")
    if rows.duplicated(["comparator_id", "stratum_value"]).any():
        raise RuntimeError("Duplicate model/quintile metric rows")
    return rows


def main() -> None:
    configure_style()
    rows = load_rows()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    x = np.arange(1, 6)

    figure, axes = plt.subplots(1, 3, figsize=(7.10, 2.85), constrained_layout=False)
    for panel_index, (axis, (metric, ylabel, ylim, reference)) in enumerate(zip(axes, METRICS)):
        for model_id, style in MODEL_STYLE.items():
            model_rows = rows.loc[rows["comparator_id"].eq(model_id)].set_index("stratum_value")
            values = model_rows.loc[QUINTILES, metric].to_numpy(dtype=float)
            axis.plot(
                x,
                values,
                color=style["color"],
                marker=style["marker"],
                markerfacecolor="white",
                markeredgewidth=0.9,
                label=style["label"],
                zorder=3,
            )
        if reference >= ylim[0]:
            axis.axhline(reference, color="#9A9A9A", linestyle=(0, (3, 2.5)), linewidth=0.75, zorder=1)
        axis.set_xlim(0.72, 5.28)
        axis.set_ylim(*ylim)
        axis.set_xticks(x, TICK_LABELS)
        axis.set_xlabel(r"External $D_5$ quintile (nearest $\rightarrow$ farthest)")
        axis.set_ylabel(ylabel)
        axis.set_title(
            ["Frozen-threshold classification", "Continuous discrimination", "Continuous ER association"][panel_index],
            loc="left",
            pad=5,
        )
        clean_axis(axis)

    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.50, 0.985),
        frameon=False,
        ncol=3,
        columnspacing=1.25,
        handlelength=1.8,
        borderaxespad=0.1,
    )
    figure.text(
        0.01,
        0.015,
        "Descriptive distance strata are relative to this project's frozen ABCB1 reference chemistry, not AZ or Deep-PK training space. "
        "ADMET-AI inhibition was endpoint-mismatched; ADMETlab had no available output.",
        ha="left",
        va="bottom",
        fontsize=6.1,
        color="#4D4D4D",
    )
    figure.subplots_adjust(left=0.075, right=0.992, bottom=0.225, top=0.755, wspace=0.38)
    for panel_index, axis in enumerate(axes):
        position = axis.get_position()
        figure.text(
            position.x0 - 0.032,
            position.y1 + 0.018,
            "abc"[panel_index],
            fontsize=10.3,
            fontweight="bold",
            ha="right",
            va="bottom",
        )
    figure.savefig(PNG_OUT, dpi=600, bbox_inches="tight")
    figure.savefig(PDF_OUT, bbox_inches="tight")
    plt.close(figure)
    print(PNG_OUT)
    print(PDF_OUT)


if __name__ == "__main__":
    main()
