#!/usr/bin/env python3
"""Plot public ABCB1 model performance across reference-relative distance.

This script is presentation-only.  It consumes the completed, frozen public
predictor stratum metrics and does not modify model outputs or analysis code.
"""

from __future__ import annotations

import json
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

BLACK = "#1A1A1A"
DARK_GREY = "#555555"
LIGHT_GREY = "#D9D9D9"
NAVY = "#2B5C8A"
ORANGE = "#D9822B"
TEAL = "#2B7A78"

QUINTILE_ORDER = ["Q1_nearest", "Q2", "Q3", "Q4", "Q5_farthest"]
QUINTILE_LABELS = ["Q1", "Q2", "Q3", "Q4", "Q5"]

MODEL_STYLE = {
    "astrazeneca_gnn_mtl_mdck_er_2025": {
        "label": "AZ MDCK ER",
        "color": NAVY,
        "marker": "o",
        "linestyle": "-",
    },
    "astrazeneca_gnn_mtl_nih_mdck_er_2025": {
        "label": "AZ NIH-MDCK ER",
        "color": ORANGE,
        "marker": "s",
        "linestyle": (0, (4.0, 2.0)),
    },
    "deeppk_pgp_substrate_2026_08_05": {
        "label": "Deep-PK substrate",
        "color": TEAL,
        "marker": "^",
        "linestyle": (0, (1.4, 1.7)),
    },
}

METRICS = [
    ("balanced_accuracy", "Balanced accuracy", (0.49, 0.80), [0.50, 0.60, 0.70, 0.80]),
    ("roc_auc", "ROC–AUC", (0.65, 0.88), [0.65, 0.70, 0.75, 0.80, 0.85]),
    ("spearman_rho", "Spearman $\\rho$ with LogER", (0.33, 0.70), [0.35, 0.45, 0.55, 0.65]),
]


def configure_style() -> None:
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
            "legend.fontsize": 7.2,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "lines.linewidth": 1.25,
            "lines.markersize": 4.4,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(BLACK)
    ax.spines["bottom"].set_color(BLACK)
    ax.tick_params(colors=BLACK)
    ax.grid(axis="y", color=LIGHT_GREY, linewidth=0.45, alpha=0.70, zorder=0)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.20,
        1.10,
        label,
        transform=ax.transAxes,
        fontsize=10.5,
        fontweight="bold",
        ha="left",
        va="top",
        color=BLACK,
    )


def load_rankable_quintiles() -> pd.DataFrame:
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)
    frame = pd.read_csv(INPUT, low_memory=False)
    required_columns = {
        "comparator_id",
        "substrate_ranking_eligible",
        "stratum_type",
        "stratum_value",
        "evaluation_layer",
        "n",
        "balanced_accuracy",
        "roc_auc",
        "spearman_rho",
    }
    missing = sorted(required_columns.difference(frame.columns))
    if missing:
        raise RuntimeError(f"Missing required columns: {missing}")

    selected = frame.loc[
        frame["substrate_ranking_eligible"].eq(True)
        & frame["stratum_type"].eq("external_d5_quintile")
        & frame["evaluation_layer"].eq("full_coverage")
    ].copy()
    selected = selected.loc[selected["comparator_id"].isin(MODEL_STYLE)].copy()
    selected["quintile_index"] = selected["stratum_value"].map(
        {value: index for index, value in enumerate(QUINTILE_ORDER, start=1)}
    )
    if selected["quintile_index"].isna().any():
        raise RuntimeError("Unexpected distance quintile label in input")
    selected["quintile_index"] = selected["quintile_index"].astype(int)

    expected_pairs = {(model, quintile) for model in MODEL_STYLE for quintile in QUINTILE_ORDER}
    observed_pairs = set(zip(selected["comparator_id"], selected["stratum_value"]))
    if observed_pairs != expected_pairs or len(selected) != 15:
        missing_pairs = sorted(expected_pairs.difference(observed_pairs))
        extra_pairs = sorted(observed_pairs.difference(expected_pairs))
        raise RuntimeError(
            f"Expected one row for each of 15 model/quintile pairs; "
            f"missing={missing_pairs}, extra={extra_pairs}, rows={len(selected)}"
        )
    if not selected["n"].between(524, 525).all():
        raise RuntimeError("Unexpected full-coverage quintile sample counts")
    return selected.sort_values(["comparator_id", "quintile_index"]).reset_index(drop=True)


def main() -> None:
    configure_style()
    data = load_rankable_quintiles()

    fig, axes = plt.subplots(1, 3, figsize=(7.09, 3.42), constrained_layout=False)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.75, bottom=0.31, wspace=0.34)
    x = np.arange(1, 6, dtype=float)
    legend_handles = []

    for panel_index, (ax, (metric, ylabel, ylim, yticks)) in enumerate(zip(axes, METRICS)):
        for model_id, style in MODEL_STYLE.items():
            model = data.loc[data["comparator_id"].eq(model_id)].sort_values("quintile_index")
            values = model[metric].to_numpy(dtype=float)
            if len(values) != 5 or not np.isfinite(values).all():
                raise RuntimeError(f"Non-finite or incomplete {metric} series for {model_id}")
            line, = ax.plot(
                x,
                values,
                color=style["color"],
                marker=style["marker"],
                markerfacecolor="white" if model_id.endswith("nih_mdck_er_2025") else style["color"],
                markeredgecolor=style["color"],
                markeredgewidth=0.9,
                linestyle=style["linestyle"],
                label=style["label"],
                zorder=3,
            )
            if panel_index == 0:
                legend_handles.append(line)

        if metric == "balanced_accuracy":
            ax.axhline(0.5, color="#A6A6A6", linestyle=(0, (3.0, 2.5)), linewidth=0.75, zorder=1)
            ax.text(5.03, 0.502, "chance", fontsize=6.0, color=DARK_GREY, ha="right", va="bottom")
        ax.set_xticks(x, QUINTILE_LABELS)
        ax.set_xlim(0.78, 5.22)
        ax.set_ylim(*ylim)
        ax.set_yticks(yticks)
        ax.set_ylabel(ylabel)
        ax.set_title(f"Full coverage: {ylabel}", loc="left", pad=6)
        clean_axis(ax)
        panel_label(ax, chr(ord("a") + panel_index))

    fig.legend(
        handles=legend_handles,
        labels=[MODEL_STYLE[key]["label"] for key in MODEL_STYLE],
        loc="upper center",
        bbox_to_anchor=(0.53, 0.965),
        ncol=3,
        frameon=False,
        handlelength=2.5,
        columnspacing=1.8,
    )
    fig.text(
        0.535,
        0.205,
        r"Reference-relative $D_5$ quintile (nearest $\rightarrow$ farthest; n=524–525 per quintile)",
        ha="center",
        va="center",
        fontsize=8.0,
        color=BLACK,
    )
    fig.text(
        0.085,
        0.105,
        "Distance strata are relative only to this study's frozen ABCB1 reference space; they do not establish OOD status "
        "relative to AZ or Deep-PK training chemistry.",
        ha="left",
        va="center",
        fontsize=6.2,
        color=DARK_GREY,
    )
    fig.text(
        0.085,
        0.052,
        "Excluded from substrate/ER ranking: ADMET-AI (P-gp inhibition endpoint mismatch); "
        "ADMETlab 3.0 (no usable benchmark output).",
        ha="left",
        va="center",
        fontsize=6.2,
        color=DARK_GREY,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(PDF_OUT, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)

    print(
        json.dumps(
            {
                "status": "complete",
                "input_rows": int(len(data)),
                "models": [MODEL_STYLE[key]["label"] for key in MODEL_STYLE],
                "evaluation_layer": "full_coverage",
                "png": str(PNG_OUT),
                "pdf": str(PDF_OUT),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
