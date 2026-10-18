#!/usr/bin/env python3
"""Plot the frozen ABCB1 evidence-channel ablation results."""

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
FIGURE_DIR = PROJECT / "outputs" / "20260805_abcb1_external_evidence_channel_ablation_v1" / "figures"
MANIFEST = PROCESSED / "abcb1_external_evidence_channel_ablation_figure_manifest_v1.csv"

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#CC79A7"
GRAY = "#6B7280"
LIGHT = "#D1D5DB"
DARK = "#1F2937"


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
    metrics = pd.read_csv(PROCESSED / "abcb1_external_evidence_channel_ablation_metrics_v1.csv")
    identity = metrics[metrics["population"].eq("identity_clean")].set_index("variant_id")
    scaffold = metrics[metrics["population"].eq("scaffold_clean")].set_index("variant_id")

    fig, axes_grid = plt.subplots(2, 2, figsize=(11.8, 8.1))
    axes = axes_grid.ravel()

    group_color = {
        "V00_forced_binary_all": GRAY,
        "V01_probability_margin_only": GRAY,
        "V02_prior_stability_only": GRAY,
        "V03_primary_support_no_local": BLUE,
        "V04_primary_channel_complete": BLUE,
        "V05_direct_support_no_local": ORANGE,
        "V06_direct_channel_complete": ORANGE,
        "V07_dual_support_only": PURPLE,
        "V08_dual_cross_model_no_local": GREEN,
        "V09_dual_cross_primary_local": GREEN,
        "V10_dual_cross_direct_local": GREEN,
        "V11_full_v0_4": DARK,
    }
    for variant_id, row in identity.iterrows():
        color = group_color[variant_id]
        axes[0].scatter(row["coverage"], row["balanced_accuracy"], color=color, s=48, zorder=3)
        scaffold_row = scaffold.loc[variant_id]
        axes[0].scatter(scaffold_row["coverage"], scaffold_row["balanced_accuracy"], facecolors="white", edgecolors=color, s=40, linewidth=1.2, zorder=2)
        axes[0].plot([row["coverage"], scaffold_row["coverage"]], [row["balanced_accuracy"], scaffold_row["balanced_accuracy"]], color=color, linewidth=0.6, alpha=0.45)
    annotation_specs = {
        "V00_forced_binary_all": ("Forced all", (-58, 7)),
        "V02_prior_stability_only": ("Prior + stability", (6, -15)),
        "V04_primary_channel_complete": ("Primary channel", (6, -14)),
        "V05_direct_support_no_local": ("Direct, no local", (6, 5)),
        "V06_direct_channel_complete": ("Direct complete", (7, 7)),
        "V08_dual_cross_model_no_local": ("Dual + cross, no local", (6, -15)),
        "V11_full_v0_4": ("Full v0.4", (7, -18)),
    }
    for variant_id, (label, offset) in annotation_specs.items():
        row = identity.loc[variant_id]
        axes[0].annotate(label, (row["coverage"], row["balanced_accuracy"]), xytext=offset, textcoords="offset points", fontsize=7.5, color=group_color[variant_id])
    axes[0].set_xlim(0.12, 1.03)
    axes[0].set_ylim(0.55, 0.80)
    axes[0].set_xlabel("Coverage")
    axes[0].set_ylabel("Balanced accuracy")
    axes[0].set_title("Evidence-channel variants: coverage versus reliability")
    axes[0].scatter([], [], color=DARK, label="Identity-clean")
    axes[0].scatter([], [], facecolors="white", edgecolors=DARK, label="Scaffold-clean")
    axes[0].legend(frameon=False, loc="lower left")
    clean_axis(axes[0])
    panel(axes[0], "A")

    selected = [
        "V02_prior_stability_only",
        "V04_primary_channel_complete",
        "V06_direct_channel_complete",
        "V07_dual_support_only",
        "V08_dual_cross_model_no_local",
        "V11_full_v0_4",
    ]
    names = ["Prior +\nstability", "Primary\nchannel", "Direct\nchannel", "Dual\nsupport", "Dual + cross\nno local", "Full\nv0.4"]
    colors = [GRAY, BLUE, ORANGE, PURPLE, GREEN, DARK]
    frame = identity.loc[selected]
    x = np.arange(len(frame))
    values = frame["balanced_accuracy"].to_numpy()
    lower = values - frame["balanced_accuracy_ci95_low"].to_numpy()
    upper = frame["balanced_accuracy_ci95_high"].to_numpy() - values
    axes[1].bar(x, values, color=colors, alpha=0.85, yerr=np.vstack([lower, upper]), error_kw={"elinewidth": 0.8, "capsize": 2})
    for index, (_, row) in enumerate(frame.iterrows()):
        axes[1].text(index, row["balanced_accuracy"] + upper[index] + 0.012, f"{row['balanced_accuracy']:.3f}\n{row['coverage']:.1%}", ha="center", va="bottom", fontsize=7.5)
    axes[1].set_xticks(x, names)
    axes[1].set_ylim(0.50, 0.86)
    axes[1].set_ylabel("Balanced accuracy (95% scaffold CI)")
    axes[1].set_title("Frozen channel-level ablation")
    clean_axis(axes[1])
    panel(axes[1], "B")

    performance_variants = [
        "V04_primary_channel_complete",
        "V06_direct_channel_complete",
        "V08_dual_cross_model_no_local",
        "V11_full_v0_4",
    ]
    performance_names = ["Primary", "Direct", "Dual + cross\nno local", "Full v0.4"]
    performance = identity.loc[performance_variants]
    x = np.arange(len(performance))
    width = 0.34
    axes[2].bar(x - width / 2, performance["sensitivity"], width, color=BLUE, alpha=0.80, label="Sensitivity")
    axes[2].bar(x + width / 2, performance["specificity"], width, color=ORANGE, alpha=0.80, label="Specificity")
    for index, (_, row) in enumerate(performance.iterrows()):
        axes[2].text(index - width / 2, row["sensitivity"] + 0.015, f"{row['sensitivity']:.2f}", ha="center", fontsize=7.5)
        axes[2].text(index + width / 2, row["specificity"] + 0.015, f"{row['specificity']:.2f}", ha="center", fontsize=7.5)
    axes[2].set_xticks(x, performance_names)
    axes[2].set_ylim(0, 1.0)
    axes[2].set_ylabel("Metric value")
    axes[2].set_title("Channel effects on sensitivity and specificity")
    axes[2].legend(frameon=False, loc="lower right")
    clean_axis(axes[2])
    panel(axes[2], "C")

    local_variants = [
        "V08_dual_cross_model_no_local",
        "V09_dual_cross_primary_local",
        "V10_dual_cross_direct_local",
        "V11_full_v0_4",
    ]
    local_names = ["No local\nguard", "+ Primary\nlocal", "+ Direct\nlocal", "+ Both\n(full)"]
    local = identity.loc[local_variants]
    x = np.arange(len(local))
    local_colors = [GREEN, BLUE, ORANGE, DARK]
    axes[3].bar(x, local["balanced_accuracy"], color=local_colors, alpha=0.85)
    axes[3].axhline(local.iloc[0]["balanced_accuracy"], color=GREEN, linestyle="--", linewidth=1)
    for index, (_, row) in enumerate(local.iterrows()):
        delta = row["balanced_accuracy"] - local.iloc[0]["balanced_accuracy"]
        axes[3].text(index, row["balanced_accuracy"] + 0.008, f"{row['balanced_accuracy']:.3f}\nΔ {delta:+.3f}", ha="center", va="bottom", fontsize=7.5)
    axes[3].set_xticks(x, local_names)
    axes[3].set_ylim(0.67, 0.80)
    axes[3].set_ylabel("Balanced accuracy")
    axes[3].set_title("Primary versus direct local-neighborhood guards")
    clean_axis(axes[3])
    panel(axes[3], "D")

    fig.suptitle("Post-hoc external ablation of ABCB1 evidence channels", fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for extension in ["png", "svg", "pdf"]:
        path = FIGURE_DIR / f"FigureS4_external_evidence_channel_ablation.{extension}"
        fig.savefig(path, dpi=300 if extension == "png" else None)
        rows.append({
            "figure": "Figure S4",
            "title": "Post-hoc external evidence-channel ablation",
            "format": extension,
            "path": str(path.relative_to(PROJECT)),
            "sha256": sha256(path),
            "size_bytes": path.stat().st_size,
            "source_tables": "abcb1_external_evidence_channel_ablation_metrics_v1.csv",
        })
    plt.close(fig)
    pd.DataFrame(rows).to_csv(MANIFEST, index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
