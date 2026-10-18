#!/usr/bin/env python3
"""Create the public-predictor rejection benchmark contrasts and figure."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "数据" / "processed"
OUTPUT_DIR = ROOT / "outputs" / "20260805_abcb1_public_predictor_benchmark_v1"
FIGURE_DIR = OUTPUT_DIR / "figures"
METRICS = PROCESSED / "abcb1_public_predictor_metrics_v1.csv"
COVERAGE = PROCESSED / "abcb1_public_predictor_coverage_curve_v1.csv"
BOOTSTRAP = PROCESSED / "abcb1_public_predictor_cluster_bootstrap_v1.csv"
CONTRASTS = PROCESSED / "abcb1_public_predictor_rejection_contrasts_v1.csv"
MANIFEST = PROCESSED / "abcb1_public_predictor_figure_manifest_v1.csv"

DISPLAY = {
    "internal_v0_4_primary_score": "Internal weak-label",
    "astrazeneca_gnn_mtl_mdck_er_2025": "AZ MDCK ER",
    "astrazeneca_gnn_mtl_nih_mdck_er_2025": "AZ NIH-MDCK ER",
    "deeppk_pgp_substrate_2026_08_05": "Deep-PK substrate",
    "admet_ai_v2_pgp_broccatelli_inhibition": "ADMET-AI inhibition*",
}
ELIGIBLE = [
    "internal_v0_4_primary_score",
    "astrazeneca_gnn_mtl_mdck_er_2025",
    "astrazeneca_gnn_mtl_nih_mdck_er_2025",
    "deeppk_pgp_substrate_2026_08_05",
]
COLORS = {
    "internal_v0_4_primary_score": "#355C7D",
    "astrazeneca_gnn_mtl_mdck_er_2025": "#2A9D8F",
    "astrazeneca_gnn_mtl_nih_mdck_er_2025": "#66A182",
    "deeppk_pgp_substrate_2026_08_05": "#E76F51",
    "admet_ai_v2_pgp_broccatelli_inhibition": "#9B9B9B",
}
POP_LABEL = {
    "all_identity_clean": "All (100%)",
    "v0_4_accepted_mask": "Frozen v0.4 mask (17.3%)",
    "native_margin_matched_453": "Native margin (17.3%)",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def build_contrasts(bootstrap: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for comparator_id, group in bootstrap.groupby("comparator_id"):
        pivot = group.pivot(index="replicate", columns="population", values="balanced_accuracy")
        definitions = {
            "v0_4_mask_minus_all": pivot["v0_4_accepted_mask"] - pivot["all_identity_clean"],
            "v0_4_mask_minus_native_margin": pivot["v0_4_accepted_mask"]
            - pivot["native_margin_matched"],
        }
        for contrast, values in definitions.items():
            rows.append(
                {
                    "comparator_id": comparator_id,
                    "display_name": DISPLAY[comparator_id],
                    "contrast": contrast,
                    "mean_difference": float(values.mean()),
                    "median_difference": float(values.median()),
                    "ci95_low": float(values.quantile(0.025)),
                    "ci95_high": float(values.quantile(0.975)),
                    "probability_difference_gt_zero": float((values > 0).mean()),
                    "n_bootstrap": int(len(values)),
                }
            )
    contrasts = pd.DataFrame(rows)
    contrasts.to_csv(CONTRASTS, index=False)
    return contrasts


def metric_ci(bootstrap: pd.DataFrame, comparator_id: str, population: str, metric: str) -> tuple[float, float]:
    values = bootstrap.loc[
        bootstrap["comparator_id"].eq(comparator_id) & bootstrap["population"].eq(population), metric
    ].dropna()
    return float(values.quantile(0.025)), float(values.quantile(0.975))


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    metrics = pd.read_csv(METRICS)
    coverage = pd.read_csv(COVERAGE)
    bootstrap = pd.read_csv(BOOTSTRAP)
    contrasts = build_contrasts(bootstrap)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "figure.dpi": 160,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2))
    ax_a, ax_b, ax_c, ax_d = axes.flat

    # A: full-coverage balanced accuracy with scaffold-cluster bootstrap intervals.
    order_a = list(reversed(list(DISPLAY)))
    for y, comparator_id in enumerate(order_a):
        row = metrics.loc[
            metrics["comparator_id"].eq(comparator_id)
            & metrics["population"].eq("all_identity_clean")
        ].iloc[0]
        low, high = metric_ci(bootstrap, comparator_id, "all_identity_clean", "balanced_accuracy")
        ax_a.errorbar(
            row["balanced_accuracy"],
            y,
            xerr=[[row["balanced_accuracy"] - low], [high - row["balanced_accuracy"]]],
            fmt="o",
            color=COLORS[comparator_id],
            capsize=3,
            markersize=7,
        )
        ax_a.text(high + 0.008, y, f"{row['balanced_accuracy']:.3f}", va="center", fontsize=8)
    ax_a.axvline(0.5, color="#777777", linestyle="--", linewidth=1)
    ax_a.set_yticks(range(len(order_a)), [DISPLAY[x] for x in order_a])
    ax_a.set_xlim(0.47, 0.80)
    ax_a.set_xlabel("Balanced accuracy (ER ≥ 2)")
    ax_a.set_title("A  Full-coverage weak binary performance", loc="left", fontweight="bold")
    ax_a.grid(axis="x", alpha=0.2)

    # B: same 17.3% coverage comparison.
    populations = ["all_identity_clean", "v0_4_accepted_mask", "native_margin_matched_453"]
    x = np.arange(len(ELIGIBLE))
    width = 0.23
    pop_colors = ["#D9D9D9", "#355C7D", "#F4A261"]
    for offset, (population, color) in enumerate(zip(populations, pop_colors)):
        values = []
        for comparator_id in ELIGIBLE:
            row = metrics.loc[
                metrics["comparator_id"].eq(comparator_id) & metrics["population"].eq(population)
            ].iloc[0]
            values.append(row["balanced_accuracy"])
        bars = ax_b.bar(x + (offset - 1) * width, values, width, color=color, label=POP_LABEL[population])
        for bar, value in zip(bars, values):
            ax_b.text(bar.get_x() + bar.get_width() / 2, value + 0.009, f"{value:.2f}", ha="center", fontsize=7)
    ax_b.axhline(0.5, color="#777777", linestyle="--", linewidth=1)
    ax_b.set_ylim(0.48, 0.96)
    ax_b.set_xticks(x, [DISPLAY[i].replace(" ", "\n", 1) for i in ELIGIBLE])
    ax_b.set_ylabel("Balanced accuracy")
    ax_b.set_title("B  Frozen mask versus predictor-native margin", loc="left", fontweight="bold")
    ax_b.legend(frameon=False, fontsize=8, loc="upper left")
    ax_b.grid(axis="y", alpha=0.2)

    # C: paired bootstrap rejection contrasts.
    contrast_order = list(reversed(ELIGIBLE))
    y_positions = np.arange(len(contrast_order))
    offsets = {"v0_4_mask_minus_all": -0.10, "v0_4_mask_minus_native_margin": 0.10}
    labels = {
        "v0_4_mask_minus_all": "Frozen mask − all",
        "v0_4_mask_minus_native_margin": "Frozen mask − native margin",
    }
    contrast_colors = {"v0_4_mask_minus_all": "#355C7D", "v0_4_mask_minus_native_margin": "#F4A261"}
    for contrast in offsets:
        rows = contrasts.set_index(["comparator_id", "contrast"])
        means, lows, highs = [], [], []
        for comparator_id in contrast_order:
            row = rows.loc[(comparator_id, contrast)]
            means.append(row["mean_difference"])
            lows.append(row["ci95_low"])
            highs.append(row["ci95_high"])
        means = np.array(means)
        ax_c.errorbar(
            means,
            y_positions + offsets[contrast],
            xerr=[means - np.array(lows), np.array(highs) - means],
            fmt="o",
            color=contrast_colors[contrast],
            capsize=3,
            label=labels[contrast],
        )
    ax_c.axvline(0, color="#555555", linewidth=1)
    ax_c.set_yticks(y_positions, [DISPLAY[i] for i in contrast_order])
    ax_c.set_xlim(-0.18, 0.22)
    ax_c.set_xlabel("Paired BAcc difference")
    ax_c.set_title("C  Scaffold-cluster bootstrap contrasts", loc="left", fontweight="bold")
    ax_c.legend(frameon=False, fontsize=8, loc="upper left")
    ax_c.grid(axis="x", alpha=0.2)

    # D: native-margin coverage curves and frozen-mask points.
    curve_models = [
        "internal_v0_4_primary_score",
        "astrazeneca_gnn_mtl_mdck_er_2025",
        "deeppk_pgp_substrate_2026_08_05",
    ]
    for comparator_id in curve_models:
        data = coverage.loc[coverage["comparator_id"].eq(comparator_id)].sort_values("coverage")
        ax_d.plot(
            data["coverage"] * 100,
            data["balanced_accuracy"],
            marker="o",
            color=COLORS[comparator_id],
            label=f"{DISPLAY[comparator_id]} native margin",
        )
        accepted = metrics.loc[
            metrics["comparator_id"].eq(comparator_id)
            & metrics["population"].eq("v0_4_accepted_mask")
        ].iloc[0]
        ax_d.scatter(
            [accepted["coverage_of_2624"] * 100],
            [accepted["balanced_accuracy"]],
            color=COLORS[comparator_id],
            marker="*",
            s=110,
            edgecolor="white",
            linewidth=0.7,
            zorder=5,
        )
    ax_d.axhline(0.5, color="#777777", linestyle="--", linewidth=1)
    ax_d.set_xlim(0, 103)
    ax_d.set_ylim(0.48, 0.98)
    ax_d.set_xlabel("Coverage retained by native score margin (%)")
    ax_d.set_ylabel("Balanced accuracy")
    ax_d.set_title("D  Coverage–reliability curves (* = frozen v0.4 mask)", loc="left", fontweight="bold")
    ax_d.legend(frameon=False, fontsize=7, loc="upper right")
    ax_d.grid(alpha=0.2)

    fig.text(
        0.01,
        0.01,
        "Post-hoc comparison on 2,624 identity-clean Biogen molecules. AZ training identities are proprietary and not auditable. "
        "*ADMET-AI predicts P-gp inhibition, not substrate status. ADMETlab 3.0 was unavailable because of an official server-side API failure.",
        fontsize=7.5,
        color="#444444",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])

    outputs = []
    stem = FIGURE_DIR / "Figure6_public_predictor_rejection_benchmark_v1"
    for suffix in ["png", "svg", "pdf"]:
        path = stem.with_suffix(f".{suffix}")
        fig.savefig(path, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)

    manifest_rows = [
        {
            "figure_id": "Figure 6",
            "title": "Public-predictor rejection benchmark",
            "file": str(path.relative_to(ROOT)),
            "format": path.suffix.lstrip("."),
            "sha256": sha256(path),
            "status": "figure_ready",
            "claim_boundary": "post-hoc comparator analysis; ADMET-AI endpoint mismatch; ADMETlab pending",
        }
        for path in outputs
    ]
    pd.DataFrame(manifest_rows).to_csv(MANIFEST, index=False)
    print(
        json.dumps(
            {
                "contrasts": str(CONTRASTS.relative_to(ROOT)),
                "figure_outputs": [str(p.relative_to(ROOT)) for p in outputs],
                "manifest": str(MANIFEST.relative_to(ROOT)),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
