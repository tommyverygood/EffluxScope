#!/usr/bin/env python3
"""Create publication-ready figures for ABCB1 external risk experiments."""

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
FIGURE_DIR = PROJECT / "outputs" / "20260805_abcb1_external_risk_experiments_v1" / "figures"
MANIFEST = PROCESSED / "abcb1_external_risk_figure_manifest_v1.csv"

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#CC79A7"
GRAY = "#6B7280"
LIGHT_GRAY = "#D1D5DB"
DARK = "#1F2937"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.edgecolor": DARK,
        "axes.linewidth": 0.8,
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


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color=LIGHT_GRAY, linewidth=0.6, alpha=0.65)
    ax.set_axisbelow(True)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.14, 1.06, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def save_figure(fig: plt.Figure, stem: str, title: str, sources: str, manifest: list[dict[str, object]]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "svg", "pdf"]:
        path = FIGURE_DIR / f"{stem}.{extension}"
        kwargs = {"dpi": 300} if extension == "png" else {}
        fig.savefig(path, **kwargs)
        manifest.append({
            "figure": stem,
            "title": title,
            "format": extension,
            "path": str(path.relative_to(PROJECT)),
            "sha256": sha256(path),
            "size_bytes": path.stat().st_size,
            "source_tables": sources,
        })
    plt.close(fig)


def figure2(manifest: list[dict[str, object]]) -> None:
    curve = pd.read_csv(PROCESSED / "abcb1_external_risk_coverage_curve_v1.csv")
    operating = pd.read_csv(PROCESSED / "abcb1_external_risk_operating_points_v1.csv")
    operating = operating[operating["population"].eq("identity_clean")].set_index("label")

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.1))
    for ax, metric, ylabel in [
        (axes[0], "balanced_accuracy", "Balanced accuracy"),
        (axes[1], "mcc", "Matthews correlation coefficient"),
    ]:
        for population, color, line_style, label in [
            ("identity_clean", BLUE, "-", "Identity-clean"),
            ("scaffold_clean", ORANGE, "--", "Scaffold-clean"),
        ]:
            frame = curve[curve["population"].eq(population)].sort_values("coverage")
            ax.plot(frame["coverage"], frame[metric], color=color, linestyle=line_style, linewidth=1.8, marker="o", markersize=3.2, label=label)
            low_col = f"{metric}_ci95_low"
            high_col = f"{metric}_ci95_high"
            ax.fill_between(frame["coverage"], frame[low_col], frame[high_col], color=color, alpha=0.10, linewidth=0)
        full = operating.loc["forced_binary_all"]
        selective = operating.loc["frozen_full_v0_4"]
        ax.scatter([full["coverage"]], [full[metric]], marker="s", s=45, color=DARK, zorder=5, label="Forced binary" if metric == "balanced_accuracy" else None)
        ax.scatter([selective["coverage"]], [selective[metric]], marker="*", s=125, color=GREEN, edgecolor="white", linewidth=0.6, zorder=6, label="Full v0.4 selective" if metric == "balanced_accuracy" else None)
        ax.annotate(f"v0.4\n{selective[metric]:.3f}", (selective["coverage"], selective[metric]), xytext=(9, -23), textcoords="offset points", color=GREEN, fontsize=8)
        ax.set_xlim(0, 1.03)
        ax.set_xlabel("Coverage")
        ax.set_ylabel(ylabel)
        ax.set_title("Confidence-only retention curve")
        clean_axis(ax)
    axes[0].legend(frameon=False, loc="lower right")
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")

    labels = [
        "forced_binary_all",
        "frozen_probability_margin_only",
        "frozen_full_v0_4",
        "forced_binary_on_v0_4_rejected_layer",
    ]
    short = ["Forced\nall", "Margin\nonly", "Full v0.4\naccepted", "Rejected\nforced"]
    colors = [GRAY, BLUE, GREEN, ORANGE]
    x = np.arange(len(labels))
    width = 0.36
    bacc = operating.loc[labels, "balanced_accuracy"].to_numpy()
    mcc = operating.loc[labels, "mcc"].to_numpy()
    bacc_low = bacc - operating.loc[labels, "balanced_accuracy_ci95_low"].to_numpy()
    bacc_high = operating.loc[labels, "balanced_accuracy_ci95_high"].to_numpy() - bacc
    mcc_low = mcc - operating.loc[labels, "mcc_ci95_low"].to_numpy()
    mcc_high = operating.loc[labels, "mcc_ci95_high"].to_numpy() - mcc
    axes[2].bar(x - width / 2, bacc, width, color=colors, alpha=0.85, label="Balanced accuracy", yerr=np.vstack([bacc_low, bacc_high]), error_kw={"elinewidth": 0.8, "capsize": 2})
    axes[2].bar(x + width / 2, mcc, width, color=colors, alpha=0.42, hatch="//", label="MCC", yerr=np.vstack([mcc_low, mcc_high]), error_kw={"elinewidth": 0.8, "capsize": 2})
    for index, label in enumerate(labels):
        coverage = operating.loc[label, "coverage"]
        axes[2].text(index, max(bacc[index] + bacc_high[index], mcc[index] + mcc_high[index]) + 0.035, f"cov. {coverage:.1%}", ha="center", va="bottom", fontsize=7, color=DARK)
    axes[2].set_xticks(x, short)
    axes[2].set_ylim(0, 0.9)
    axes[2].set_ylabel("Metric value")
    axes[2].set_title("Frozen operating-point comparison")
    axes[2].legend(frameon=False, loc="upper right")
    clean_axis(axes[2])
    panel_label(axes[2], "C")

    fig.suptitle("External coverage–reliability trade-off on the identity-clean Biogen benchmark", fontsize=12, fontweight="bold", y=1.025)
    fig.tight_layout()
    save_figure(
        fig,
        "Figure2_external_coverage_reliability",
        "External coverage-reliability trade-off",
        "abcb1_external_risk_coverage_curve_v1.csv|abcb1_external_risk_operating_points_v1.csv",
        manifest,
    )


def figure3(manifest: list[dict[str, object]]) -> None:
    cumulative = pd.read_csv(PROCESSED / "abcb1_external_rejector_cumulative_ablation_v1.csv")
    cumulative = cumulative[cumulative["population"].eq("identity_clean")].sort_values("stage_index")
    loo = pd.read_csv(PROCESSED / "abcb1_external_rejector_leave_one_out_v1.csv")
    loo = loo[loo["population"].eq("identity_clean")]
    burden = pd.read_csv(PROCESSED / "abcb1_external_rejector_rule_burden_v1.csv")
    burden = burden[burden["population"].eq("identity_clean")]

    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.3))
    axes[0].plot(cumulative["coverage"], cumulative["balanced_accuracy"], color=BLUE, marker="o", linewidth=1.8)
    short_stage = ["All", "Margin", "Bootstrap", "AD", "Model", "Local", "Stereo"]
    for (_, row), text in zip(cumulative.iterrows(), short_stage):
        offset = (4, 7) if text not in {"All", "Local"} else (5, -14)
        axes[0].annotate(text, (row["coverage"], row["balanced_accuracy"]), xytext=offset, textcoords="offset points", fontsize=7.5)
    axes[0].scatter(cumulative.iloc[-1]["coverage"], cumulative.iloc[-1]["balanced_accuracy"], marker="*", s=125, color=GREEN, edgecolor="white", linewidth=0.6, zorder=5)
    axes[0].set_xlim(0.14, 1.03)
    axes[0].set_ylim(0.54, 0.80)
    axes[0].set_xlabel("Coverage after cumulative rules")
    axes[0].set_ylabel("Balanced accuracy")
    axes[0].set_title("Cumulative rejection-rule path")
    clean_axis(axes[0])
    panel_label(axes[0], "A")

    order = [
        "full_rule_set",
        "without_bootstrap_direction_stability",
        "without_applicability_domain",
        "without_supported_model_direction_agreement",
        "without_local_neighborhood_agreement",
        "without_stereo_collision_guard",
    ]
    names = ["Full rules", "– Bootstrap", "– AD", "– Model agreement", "– Local agreement", "– Stereo guard"]
    frame = loo.set_index("label").loc[order].reset_index()
    y = np.arange(len(frame))
    bar_colors = [GREEN, BLUE, BLUE, ORANGE, PURPLE, BLUE]
    axes[1].barh(y, frame["balanced_accuracy"], color=bar_colors, alpha=0.82)
    axes[1].axvline(float(frame.iloc[0]["balanced_accuracy"]), color=DARK, linewidth=0.9, linestyle="--")
    for index, row in frame.iterrows():
        axes[1].text(row["balanced_accuracy"] + 0.006, index, f"{row['balanced_accuracy']:.3f} | cov. {row['coverage']:.1%}", va="center", fontsize=7.5)
    axes[1].set_yticks(y, names)
    axes[1].invert_yaxis()
    axes[1].set_xlim(0.55, 0.81)
    axes[1].set_xlabel("Balanced accuracy")
    axes[1].set_title("Leave-one-rule-out analysis")
    clean_axis(axes[1])
    panel_label(axes[1], "B")

    burden_names = {
        "probability_margin": "Probability gray zone",
        "bootstrap_direction_stability": "Bootstrap instability",
        "applicability_domain": "Outside AD",
        "supported_model_direction_agreement": "Model disagreement",
        "local_neighborhood_agreement": "Local disagreement",
        "stereo_collision_guard": "Stereo collision",
    }
    burden = burden.assign(display=burden["rule_family"].map(burden_names)).sort_values("triggered_n")
    axes[2].barh(np.arange(len(burden)), burden["triggered_n"], color=GRAY, alpha=0.82)
    axes[2].set_yticks(np.arange(len(burden)), burden["display"])
    for index, row in burden.reset_index(drop=True).iterrows():
        axes[2].text(row["triggered_n"] + 20, index, f"{int(row['triggered_n']):,}", va="center", fontsize=8)
    axes[2].set_xlim(0, max(burden["triggered_n"].max() * 1.16, 50))
    axes[2].set_xlabel("Triggered structures (overlapping)")
    axes[2].set_title("Frozen rule burden")
    clean_axis(axes[2])
    panel_label(axes[2], "C")

    fig.suptitle("External decomposition of the v0.4 rejector", fontsize=12, fontweight="bold", y=1.025)
    fig.tight_layout()
    save_figure(
        fig,
        "Figure3_external_rejector_ablation",
        "External rejector ablation",
        "abcb1_external_rejector_cumulative_ablation_v1.csv|abcb1_external_rejector_leave_one_out_v1.csv|abcb1_external_rejector_rule_burden_v1.csv",
        manifest,
    )


def figure4(manifest: list[dict[str, object]]) -> None:
    raw = pd.read_csv(PROCESSED / "abcb1_quantity_biogen_analysis_rows_v1.csv", low_memory=False)
    raw = raw[raw["identity_clean"].eq("yes")].copy()
    raw["primary_calibrated_probability"] = pd.to_numeric(raw["primary_calibrated_probability"], errors="coerce")
    raw["log10_mdr1_mdck_er"] = pd.to_numeric(raw["log10_mdr1_mdck_er"], errors="coerce")
    raw = raw.dropna(subset=["primary_calibrated_probability", "log10_mdr1_mdck_er"])
    raw["selection_group"] = np.where(
        raw["risk_call"].eq("high_ABCB1_efflux_risk"), "Accepted high",
        np.where(raw["risk_call"].eq("low_ABCB1_efflux_risk"), "Accepted low", "Rejected"),
    )
    bins = pd.read_csv(PROCESSED / "abcb1_external_calibration_bins_v1.csv")
    bins = bins[bins["strategy"].eq("equal_width")]

    fig, axes = plt.subplots(1, 3, figsize=(13.3, 4.2))
    calibration_specs = [
        ("identity_clean_all", BLUE, "All"),
        ("identity_clean_v0_4_accepted", GREEN, "v0.4 accepted"),
        ("identity_clean_v0_4_rejected", ORANGE, "v0.4 rejected"),
    ]
    axes[0].plot([0, 1], [0, 1], color=GRAY, linestyle="--", linewidth=1, label="Ideal")
    for population, color, label in calibration_specs:
        frame = bins[bins["population"].eq(population)].sort_values("mean_predicted_probability")
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
    axes[0].set_title("External calibration diagnostic")
    axes[0].legend(frameon=False, loc="upper left")
    clean_axis(axes[0])
    panel_label(axes[0], "A")

    hb = axes[1].hexbin(
        raw["primary_calibrated_probability"], raw["log10_mdr1_mdck_er"],
        gridsize=32, mincnt=1, cmap="Blues", linewidths=0.15,
    )
    quantile_bins = pd.qcut(raw["primary_calibrated_probability"], 10, duplicates="drop")
    median = raw.groupby(quantile_bins, observed=True).agg(
        probability=("primary_calibrated_probability", "mean"),
        endpoint=("log10_mdr1_mdck_er", "median"),
    )
    axes[1].plot(median["probability"], median["endpoint"], color=ORANGE, marker="o", markersize=3.5, linewidth=1.5, label="Decile median")
    axes[1].axhline(np.log10(2.0), color=DARK, linestyle="--", linewidth=1, label="ER = 2")
    axes[1].axvline(0.5, color=GRAY, linestyle=":", linewidth=1)
    axes[1].set_xlabel("Primary model probability")
    axes[1].set_ylabel("log10 MDR1-MDCK efflux ratio")
    axes[1].set_title("Probability versus continuous endpoint")
    axes[1].legend(frameon=False, loc="upper left")
    colorbar = fig.colorbar(hb, ax=axes[1], fraction=0.046, pad=0.02)
    colorbar.set_label("Structures per hexagon", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)
    clean_axis(axes[1])
    panel_label(axes[1], "B")

    order = ["Accepted low", "Accepted high", "Rejected"]
    data = [raw.loc[raw["selection_group"].eq(group), "log10_mdr1_mdck_er"].to_numpy() for group in order]
    box = axes[2].boxplot(data, labels=order, patch_artist=True, showfliers=False, widths=0.6, medianprops={"color": DARK, "linewidth": 1.4})
    for patch, color in zip(box["boxes"], [BLUE, GREEN, ORANGE]):
        patch.set_facecolor(color)
        patch.set_alpha(0.58)
        patch.set_edgecolor(color)
    for index, values in enumerate(data, start=1):
        axes[2].text(index, axes[2].get_ylim()[1], f"n={len(values):,}", ha="center", va="top", fontsize=7.5)
    axes[2].axhline(np.log10(2.0), color=DARK, linestyle="--", linewidth=1)
    axes[2].set_ylabel("log10 MDR1-MDCK efflux ratio")
    axes[2].set_title("Endpoint enrichment by selective layer")
    axes[2].tick_params(axis="x", rotation=15)
    clean_axis(axes[2])
    panel_label(axes[2], "C")

    fig.suptitle("External calibration and error-layer diagnostics", fontsize=12, fontweight="bold", y=1.025)
    fig.tight_layout()
    save_figure(
        fig,
        "Figure4_external_calibration_error",
        "External calibration and error-layer diagnostics",
        "abcb1_external_calibration_bins_v1.csv|abcb1_quantity_biogen_analysis_rows_v1.csv",
        manifest,
    )


def main() -> None:
    style()
    manifest: list[dict[str, object]] = []
    figure2(manifest)
    figure3(manifest)
    figure4(manifest)
    pd.DataFrame(manifest).to_csv(MANIFEST, index=False)
    print(pd.DataFrame(manifest).to_string(index=False))


if __name__ == "__main__":
    main()
