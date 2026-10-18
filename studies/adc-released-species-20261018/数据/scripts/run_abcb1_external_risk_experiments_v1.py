#!/usr/bin/env python3
"""Run post-hoc external risk experiments for the frozen ABCB1 v0.4 screener.

This script does not refit a model or select an operating threshold. It uses the
locked Biogen external-analysis rows to quantify coverage/reliability trade-offs,
reconstruct and ablate the frozen rejection rules, audit calibration, and
stratify forced-decision errors. All confidence intervals use the inherited
Murcko-scaffold cluster bootstrap.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, roc_auc_score


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"
INPUT = PROCESSED / "abcb1_quantity_biogen_analysis_rows_v1.csv"
MODEL = PROJECT / "模型" / "ABCB1_selective_screener_v0_4_dev" / "ABCB1_selective_uncertainty_v0_4_dev.joblib"
UPSTREAM_PROTOCOL = PROCESSED / "abcb1_quantity_external_validation_protocol_lock_v1.json"
PROTOCOL = PROCESSED / "abcb1_external_risk_experiment_protocol_v1.json"

CURVE = PROCESSED / "abcb1_external_risk_coverage_curve_v1.csv"
OPERATING_POINTS = PROCESSED / "abcb1_external_risk_operating_points_v1.csv"
CUMULATIVE = PROCESSED / "abcb1_external_rejector_cumulative_ablation_v1.csv"
LEAVE_ONE_OUT = PROCESSED / "abcb1_external_rejector_leave_one_out_v1.csv"
RULE_BURDEN = PROCESSED / "abcb1_external_rejector_rule_burden_v1.csv"
CALIBRATION_BINS = PROCESSED / "abcb1_external_calibration_bins_v1.csv"
CALIBRATION_SUMMARY = PROCESSED / "abcb1_external_calibration_summary_v1.csv"
ERROR_STRATA = PROCESSED / "abcb1_external_error_strata_v1.csv"
ENRICHMENT = PROCESSED / "abcb1_external_selective_enrichment_bootstrap_v1.csv"
SUMMARY = PROCESSED / "abcb1_external_risk_experiment_summary_v1.json"

HIGH = "high_ABCB1_efflux_risk"
LOW = "low_ABCB1_efflux_risk"
ABSTAIN = "abstain"
P_COL = "primary_calibrated_probability"
Y_COL = "frozen_er2_binary_label"
CLUSTER_COL = "bootstrap_cluster"
SEED = 20260805
N_BOOT = 2000

RULE_FAMILIES = {
    "probability_margin": {"probability_gray_zone"},
    "bootstrap_direction_stability": {"bootstrap_direction_unstable"},
    "applicability_domain": {"outside_both_reference_panels", "outside_physicochemical_envelope"},
    "supported_model_direction_agreement": {"supported_model_direction_disagreement"},
    "local_neighborhood_agreement": {
        "primary_local_neighborhood_strongly_opposite",
        "direct_local_neighborhood_strongly_opposite",
    },
    "stereo_collision_guard": {"stereo_blind_fingerprint_collision"},
}
RULE_ORDER = list(RULE_FAMILIES)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean_frame(df: pd.DataFrame, population: str) -> pd.DataFrame:
    if population == "identity_clean":
        work = df[df["identity_clean"].eq("yes")].copy()
    elif population == "scaffold_clean":
        work = df[df["scaffold_clean"].eq("yes")].copy()
    else:
        raise ValueError(population)
    work[P_COL] = pd.to_numeric(work[P_COL], errors="coerce")
    work[Y_COL] = pd.to_numeric(work[Y_COL], errors="coerce")
    work = work.dropna(subset=[P_COL, Y_COL, CLUSTER_COL]).copy()
    work[Y_COL] = work[Y_COL].astype(int)
    work["forced_prediction"] = (work[P_COL] >= 0.5).astype(int)
    work["confidence_margin"] = (work[P_COL] - 0.5).abs()
    work["reason_set"] = work["abstention_reasons_json"].fillna("[]").map(parse_reasons)
    for family, reasons in RULE_FAMILIES.items():
        work[f"pass_{family}"] = work["reason_set"].map(lambda value: not bool(value & reasons))
    return work.reset_index(drop=True)


def parse_reasons(value: object) -> set[str]:
    try:
        parsed = json.loads(str(value))
        return {str(item) for item in parsed}
    except (json.JSONDecodeError, TypeError):
        return {"unparseable_reason_json"}


def divide(num: np.ndarray | float, den: np.ndarray | float) -> np.ndarray | float:
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.asarray(den) > 0, np.asarray(num) / np.asarray(den), np.nan)


def counts_from_arrays(y: np.ndarray, pred: np.ndarray) -> np.ndarray:
    return np.array([
        np.sum((y == 1) & (pred == 1)),
        np.sum((y == 0) & (pred == 0)),
        np.sum((y == 0) & (pred == 1)),
        np.sum((y == 1) & (pred == 0)),
    ], dtype=float)


def metrics_from_counts(counts: np.ndarray) -> dict[str, float]:
    tp, tn, fp, fn = [float(x) for x in counts]
    n = tp + tn + fp + fn
    sensitivity = float(divide(tp, tp + fn))
    specificity = float(divide(tn, tn + fp))
    accuracy = float(divide(tp + tn, n))
    bacc = float(np.nanmean([sensitivity, specificity])) if np.isfinite(sensitivity) and np.isfinite(specificity) else np.nan
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / denominator if denominator > 0 else np.nan
    return {
        "n": int(n),
        "n_positive": int(tp + fn),
        "n_negative": int(tn + fp),
        "tp": int(tp),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "prevalence": float(divide(tp + fn, n)),
        "accuracy": accuracy,
        "error_rate": 1.0 - accuracy if np.isfinite(accuracy) else np.nan,
        "balanced_accuracy": bacc,
        "mcc": mcc,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "high_risk_precision": float(divide(tp, tp + fp)),
        "low_risk_npv": float(divide(tn, tn + fn)),
    }


def metrics_for_frame(df: pd.DataFrame, prediction_col: str = "forced_prediction") -> dict[str, float]:
    y = df[Y_COL].astype(int).to_numpy()
    pred = df[prediction_col].astype(int).to_numpy()
    return metrics_from_counts(counts_from_arrays(y, pred))


def grouped_confusion(df: pd.DataFrame, prediction_col: str = "forced_prediction") -> np.ndarray:
    rows = []
    for _, group in df.groupby(CLUSTER_COL, sort=False):
        rows.append(counts_from_arrays(group[Y_COL].astype(int).to_numpy(), group[prediction_col].astype(int).to_numpy()))
    return np.asarray(rows, dtype=float)


def vector_metrics(counts: np.ndarray) -> dict[str, np.ndarray]:
    tp, tn, fp, fn = [counts[:, index] for index in range(4)]
    n = tp + tn + fp + fn
    sensitivity = divide(tp, tp + fn)
    specificity = divide(tn, tn + fp)
    accuracy = divide(tp + tn, n)
    bacc = (sensitivity + specificity) / 2.0
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = divide(tp * tn - fp * fn, denominator)
    return {"accuracy": accuracy, "balanced_accuracy": bacc, "mcc": mcc}


def bootstrap_metric_intervals(df: pd.DataFrame, seed_offset: int) -> dict[str, float | int]:
    grouped = grouped_confusion(df)
    if len(grouped) < 2:
        return {
            f"{metric}_{suffix}": np.nan
            for metric in ["accuracy", "balanced_accuracy", "mcc"]
            for suffix in ["ci95_low", "ci95_high"]
        } | {"bootstrap_valid_replicates": 0, "n_bootstrap_clusters": int(len(grouped))}
    rng = np.random.default_rng(SEED + seed_offset)
    values = {"accuracy": [], "balanced_accuracy": [], "mcc": []}
    remaining = N_BOOT
    while remaining:
        batch = min(100, remaining)
        sampled = rng.integers(0, len(grouped), size=(batch, len(grouped)))
        totals = grouped[sampled].sum(axis=1)
        batch_metrics = vector_metrics(totals)
        for metric in values:
            values[metric].extend(batch_metrics[metric][np.isfinite(batch_metrics[metric])].tolist())
        remaining -= batch
    output: dict[str, float | int] = {
        "bootstrap_valid_replicates": min(len(values[metric]) for metric in values),
        "n_bootstrap_clusters": int(len(grouped)),
    }
    for metric, metric_values in values.items():
        if metric_values:
            output[f"{metric}_ci95_low"] = float(np.quantile(metric_values, 0.025))
            output[f"{metric}_ci95_high"] = float(np.quantile(metric_values, 0.975))
        else:
            output[f"{metric}_ci95_low"] = np.nan
            output[f"{metric}_ci95_high"] = np.nan
    return output


def result_row(
    selected: pd.DataFrame,
    denominator_n: int,
    *,
    population: str,
    experiment: str,
    label: str,
    seed_offset: int,
) -> dict[str, object]:
    output: dict[str, object] = {
        "population": population,
        "experiment": experiment,
        "label": label,
        "denominator_n": int(denominator_n),
        "accepted_n": int(len(selected)),
        "coverage": float(len(selected) / denominator_n) if denominator_n else np.nan,
    }
    output.update(metrics_for_frame(selected))
    output.update(bootstrap_metric_intervals(selected, seed_offset))
    return output


def coverage_curve(work: pd.DataFrame, population: str, target_coverages: Iterable[float]) -> list[dict[str, object]]:
    ordered = work.sort_values(["confidence_margin", "external_row_id"], ascending=[False, True], kind="stable")
    rows = []
    for index, target in enumerate(target_coverages):
        take = len(work) if float(target) >= 1 else max(2, int(round(float(target) * len(work))))
        selected = ordered.iloc[:take].copy()
        row = result_row(
            selected,
            len(work),
            population=population,
            experiment="confidence_ranked_probability_margin_curve",
            label=f"target_coverage_{float(target):.3f}",
            seed_offset=1000 + 100 * (population == "scaffold_clean") + index,
        )
        row["target_coverage"] = float(target)
        row["minimum_retained_absolute_probability_margin"] = float(selected["confidence_margin"].min())
        rows.append(row)
    return rows


def operating_points(work: pd.DataFrame, population: str) -> list[dict[str, object]]:
    margin = work[(work[P_COL] <= 0.35) | (work[P_COL] >= 0.7)].copy()
    accepted = work[work["risk_call"].isin([HIGH, LOW])].copy()
    accepted["frozen_selective_prediction"] = accepted["risk_call"].map({HIGH: 1, LOW: 0}).astype(int)
    rejected = work[work["risk_call"].eq(ABSTAIN)].copy()
    specs = [
        (work, "forced_binary_all", "forced_prediction"),
        (margin, "frozen_probability_margin_only", "forced_prediction"),
        (accepted, "frozen_full_v0_4", "frozen_selective_prediction"),
        (rejected, "forced_binary_on_v0_4_rejected_layer", "forced_prediction"),
    ]
    rows = []
    for index, (selected, label, pred_col) in enumerate(specs):
        frame = selected.copy()
        frame["analysis_prediction"] = frame[pred_col]
        rows.append(result_row(
            frame,
            len(work),
            population=population,
            experiment="frozen_operating_point",
            label=label,
            seed_offset=2000 + 100 * (population == "scaffold_clean") + index,
        ))
    return rows


def cumulative_ablation(work: pd.DataFrame, population: str) -> list[dict[str, object]]:
    rows = []
    current = np.ones(len(work), dtype=bool)
    previous_n = len(work)
    base = result_row(
        work.copy(), len(work), population=population, experiment="cumulative_rejector_ablation",
        label="00_forced_binary_all", seed_offset=3000 + 100 * (population == "scaffold_clean"),
    )
    base.update({"stage_index": 0, "rule_added": "none", "newly_rejected_n": 0})
    rows.append(base)
    for index, family in enumerate(RULE_ORDER, start=1):
        current &= work[f"pass_{family}"].to_numpy(dtype=bool)
        selected = work[current].copy()
        row = result_row(
            selected, len(work), population=population, experiment="cumulative_rejector_ablation",
            label=f"{index:02d}_plus_{family}", seed_offset=3000 + 100 * (population == "scaffold_clean") + index,
        )
        row.update({"stage_index": index, "rule_added": family, "newly_rejected_n": int(previous_n - len(selected))})
        rows.append(row)
        previous_n = len(selected)
    return rows


def leave_one_out(work: pd.DataFrame, population: str) -> list[dict[str, object]]:
    rows = []
    margin_pass = work["pass_probability_margin"].to_numpy(dtype=bool)
    downstream = RULE_ORDER[1:]
    specifications: list[tuple[str, list[str]]] = [
        ("probability_margin_only", []),
        ("full_rule_set", downstream),
    ]
    specifications.extend((f"without_{family}", [rule for rule in downstream if rule != family]) for family in downstream)
    full_mask = margin_pass.copy()
    for family in downstream:
        full_mask &= work[f"pass_{family}"].to_numpy(dtype=bool)
    full_n = int(full_mask.sum())
    for index, (label, included_rules) in enumerate(specifications):
        mask = margin_pass.copy()
        for family in included_rules:
            mask &= work[f"pass_{family}"].to_numpy(dtype=bool)
        selected = work[mask].copy()
        row = result_row(
            selected, len(work), population=population, experiment="rejector_leave_one_rule_out",
            label=label, seed_offset=4000 + 100 * (population == "scaffold_clean") + index,
        )
        row["included_downstream_rules"] = "|".join(included_rules)
        row["rows_recovered_vs_full"] = int(len(selected) - full_n)
        rows.append(row)
    return rows


def rule_burden(work: pd.DataFrame, population: str) -> list[dict[str, object]]:
    rows = []
    for family in RULE_ORDER:
        triggered = work[~work[f"pass_{family}"]].copy()
        row = {
            "population": population,
            "rule_family": family,
            "triggered_n": int(len(triggered)),
            "triggered_fraction": float(len(triggered) / len(work)),
            "note": "overlapping rule triggers; rows are not mutually exclusive",
        }
        row.update(metrics_for_frame(triggered) if len(triggered) else metrics_from_counts(np.zeros(4)))
        rows.append(row)
    return rows


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return np.nan, np.nan
    p = successes / total
    denominator = 1.0 + z * z / total
    center = (p + z * z / (2.0 * total)) / denominator
    half = z * math.sqrt(p * (1.0 - p) / total + z * z / (4.0 * total * total)) / denominator
    return center - half, center + half


def calibration_rows(work: pd.DataFrame, population_label: str, strategy: str, n_bins: int = 10) -> list[dict[str, object]]:
    frame = work.copy()
    if strategy == "equal_width":
        edges = np.linspace(0.0, 1.0, n_bins + 1)
        frame["bin"] = pd.cut(frame[P_COL], bins=edges, labels=False, include_lowest=True, right=True)
    elif strategy == "adaptive":
        frame["bin"] = pd.qcut(frame[P_COL], q=n_bins, labels=False, duplicates="drop")
    else:
        raise ValueError(strategy)
    rows = []
    for bin_id, group in frame.dropna(subset=["bin"]).groupby("bin", sort=True):
        n = len(group)
        observed = float(group[Y_COL].mean())
        lower, upper = wilson_interval(int(group[Y_COL].sum()), n)
        rows.append({
            "population": population_label,
            "strategy": strategy,
            "bin_id": int(bin_id),
            "n": int(n),
            "probability_min": float(group[P_COL].min()),
            "probability_max": float(group[P_COL].max()),
            "mean_predicted_probability": float(group[P_COL].mean()),
            "observed_positive_fraction": observed,
            "observed_positive_wilson95_low": lower,
            "observed_positive_wilson95_high": upper,
            "absolute_calibration_gap": abs(observed - float(group[P_COL].mean())),
        })
    return rows


def ece_from_bins(rows: list[dict[str, object]]) -> float:
    total = sum(int(row["n"]) for row in rows)
    if total == 0:
        return np.nan
    return float(sum(int(row["n"]) * float(row["absolute_calibration_gap"]) for row in rows) / total)


def calibration_metrics(work: pd.DataFrame, population_label: str) -> tuple[dict[str, object], list[dict[str, object]]]:
    y = work[Y_COL].astype(int).to_numpy()
    p = work[P_COL].astype(float).clip(1e-6, 1 - 1e-6).to_numpy()
    equal_rows = calibration_rows(work, population_label, "equal_width")
    adaptive_rows = calibration_rows(work, population_label, "adaptive")
    intercept = np.nan
    slope = np.nan
    if len(np.unique(y)) == 2:
        logit = np.log(p / (1.0 - p)).reshape(-1, 1)
        model = LogisticRegression(C=np.inf, l1_ratio=0.0, solver="lbfgs", max_iter=5000)
        model.fit(logit, y)
        intercept = float(model.intercept_[0])
        slope = float(model.coef_[0, 0])
    summary = {
        "population": population_label,
        "n": int(len(work)),
        "n_positive": int(y.sum()),
        "n_negative": int(len(y) - y.sum()),
        "observed_prevalence": float(y.mean()),
        "mean_predicted_probability": float(p.mean()),
        "brier_score": float(brier_score_loss(y, p)),
        "log_loss": float(log_loss(y, p, labels=[0, 1])),
        "ece_equal_width_10": ece_from_bins(equal_rows),
        "ece_adaptive_10": ece_from_bins(adaptive_rows),
        "calibration_intercept": intercept,
        "calibration_slope": slope,
        "roc_auc": float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else np.nan,
        "average_precision": float(average_precision_score(y, p)) if len(np.unique(y)) == 2 else np.nan,
    }
    return summary, equal_rows + adaptive_rows


def error_strata(work: pd.DataFrame, population: str) -> list[dict[str, object]]:
    rows = []

    def add(stratum_type: str, level: str, frame: pd.DataFrame, overlap_note: str = "") -> None:
        if frame.empty:
            return
        row = {
            "population": population,
            "stratum_type": stratum_type,
            "stratum_level": level,
            "fraction_of_population": float(len(frame) / len(work)),
            "overlap_note": overlap_note,
        }
        row.update(metrics_for_frame(frame))
        rows.append(row)

    accepted_mask = work["risk_call"].isin([HIGH, LOW])
    add("selection_layer", "v0_4_accepted", work[accepted_mask])
    add("selection_layer", "v0_4_rejected_forced_binary", work[~accepted_mask])
    for value, frame in work.groupby("risk_call", dropna=False):
        add("risk_call", str(value), frame)
    for column in ["applicability_status", "uncertainty_status"]:
        for value, frame in work.groupby(column, dropna=False):
            add(column, str(value), frame)
    all_reasons = sorted({reason for values in work["reason_set"] for reason in values})
    for reason in all_reasons:
        frame = work[work["reason_set"].map(lambda value: reason in value)]
        add("abstention_reason", reason, frame, "overlapping reason sets; do not sum rows")
    return rows


def joint_selective_enrichment(work: pd.DataFrame, population: str) -> list[dict[str, object]]:
    accepted_mask = work["risk_call"].isin([HIGH, LOW]).to_numpy()
    accepted = work[accepted_mask]
    rejected = work[~accepted_mask]
    accepted_point = metrics_for_frame(accepted)
    rejected_point = metrics_for_frame(rejected)

    grouped_rows = []
    for _, group in work.groupby(CLUSTER_COL, sort=False):
        group_mask = group["risk_call"].isin([HIGH, LOW]).to_numpy()
        a = group[group_mask]
        r = group[~group_mask]
        grouped_rows.append(np.concatenate([
            counts_from_arrays(a[Y_COL].astype(int).to_numpy(), a["forced_prediction"].astype(int).to_numpy()),
            counts_from_arrays(r[Y_COL].astype(int).to_numpy(), r["forced_prediction"].astype(int).to_numpy()),
        ]))
    grouped = np.asarray(grouped_rows, dtype=float)
    rng = np.random.default_rng(SEED + 9000 + 100 * (population == "scaffold_clean"))
    boot_values = {"accuracy": [], "balanced_accuracy": [], "mcc": [], "error_rate": []}
    remaining = N_BOOT
    while remaining:
        batch = min(100, remaining)
        sampled = rng.integers(0, len(grouped), size=(batch, len(grouped)))
        totals = grouped[sampled].sum(axis=1)
        accepted_metrics = vector_metrics(totals[:, :4])
        rejected_metrics = vector_metrics(totals[:, 4:])
        for metric in ["accuracy", "balanced_accuracy", "mcc"]:
            delta = accepted_metrics[metric] - rejected_metrics[metric]
            boot_values[metric].extend(delta[np.isfinite(delta)].tolist())
        error_delta = (1.0 - rejected_metrics["accuracy"]) - (1.0 - accepted_metrics["accuracy"])
        boot_values["error_rate"].extend(error_delta[np.isfinite(error_delta)].tolist())
        remaining -= batch

    rows = []
    for metric in ["accuracy", "balanced_accuracy", "mcc", "error_rate"]:
        values = np.asarray(boot_values[metric], dtype=float)
        if metric == "error_rate":
            point = rejected_point["error_rate"] - accepted_point["error_rate"]
            contrast = "rejected_minus_accepted"
        else:
            point = accepted_point[metric] - rejected_point[metric]
            contrast = "accepted_minus_rejected"
        rows.append({
            "population": population,
            "metric": metric,
            "contrast": contrast,
            "accepted_value": accepted_point[metric],
            "rejected_forced_value": rejected_point[metric],
            "favorable_difference": point,
            "ci95_low": float(np.quantile(values, 0.025)),
            "ci95_high": float(np.quantile(values, 0.975)),
            "bootstrap_probability_favorable_gt_0": float(np.mean(values > 0)),
            "valid_bootstrap_replicates": int(len(values)),
            "n_bootstrap_clusters": int(len(grouped)),
        })
    return rows


def round_frame(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    for column in result.select_dtypes(include=["number"]).columns:
        result[column] = result[column].round(8)
    return result


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text())
    assert protocol["status"].startswith("locked_posthoc_secondary_analysis")
    for entry in protocol["frozen_inputs"].values():
        path = PROJECT / entry["path"]
        assert path.exists(), path
        assert sha256(path) == entry["sha256"], f"Hash mismatch: {path}"
    assert protocol["experiments"]["cluster_bootstrap"]["replicates"] == N_BOOT
    assert protocol["experiments"]["cluster_bootstrap"]["seed"] == SEED

    raw = pd.read_csv(INPUT, low_memory=False)
    populations = {name: clean_frame(raw, name) for name in ["identity_clean", "scaffold_clean"]}
    target_coverages = protocol["experiments"]["probability_margin_coverage_reliability"]["target_coverages"]

    curve_rows: list[dict[str, object]] = []
    operating_rows: list[dict[str, object]] = []
    cumulative_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []
    burden_rows: list[dict[str, object]] = []
    calibration_summary_rows: list[dict[str, object]] = []
    calibration_bin_rows: list[dict[str, object]] = []
    error_rows: list[dict[str, object]] = []
    enrichment_rows: list[dict[str, object]] = []

    for population, work in populations.items():
        curve_rows.extend(coverage_curve(work, population, target_coverages))
        operating_rows.extend(operating_points(work, population))
        cumulative_rows.extend(cumulative_ablation(work, population))
        loo_rows.extend(leave_one_out(work, population))
        burden_rows.extend(rule_burden(work, population))
        error_rows.extend(error_strata(work, population))
        enrichment_rows.extend(joint_selective_enrichment(work, population))

    identity = populations["identity_clean"]
    scaffold = populations["scaffold_clean"]
    calibration_populations = {
        "identity_clean_all": identity,
        "identity_clean_v0_4_accepted": identity[identity["risk_call"].isin([HIGH, LOW])],
        "identity_clean_v0_4_rejected": identity[identity["risk_call"].eq(ABSTAIN)],
        "scaffold_clean_all": scaffold,
    }
    for label, frame in calibration_populations.items():
        calibration_summary, calibration_bins = calibration_metrics(frame, label)
        calibration_summary_rows.append(calibration_summary)
        calibration_bin_rows.extend(calibration_bins)

    outputs = {
        CURVE: round_frame(pd.DataFrame(curve_rows)),
        OPERATING_POINTS: round_frame(pd.DataFrame(operating_rows)),
        CUMULATIVE: round_frame(pd.DataFrame(cumulative_rows)),
        LEAVE_ONE_OUT: round_frame(pd.DataFrame(loo_rows)),
        RULE_BURDEN: round_frame(pd.DataFrame(burden_rows)),
        CALIBRATION_BINS: round_frame(pd.DataFrame(calibration_bin_rows)),
        CALIBRATION_SUMMARY: round_frame(pd.DataFrame(calibration_summary_rows)),
        ERROR_STRATA: round_frame(pd.DataFrame(error_rows)),
        ENRICHMENT: round_frame(pd.DataFrame(enrichment_rows)),
    }
    for path, frame in outputs.items():
        frame.to_csv(path, index=False)

    identity_full = outputs[OPERATING_POINTS]
    identity_full = identity_full[(identity_full["population"] == "identity_clean")]
    point_lookup = identity_full.set_index("label").to_dict(orient="index")
    identity_final = outputs[CUMULATIVE]
    identity_final = identity_final[(identity_final["population"] == "identity_clean")].sort_values("stage_index").iloc[-1]
    frozen_accepted_ids = set(identity.loc[identity["risk_call"].isin([HIGH, LOW]), "external_row_id"].astype(str))
    reconstructed_mask = np.ones(len(identity), dtype=bool)
    for family in RULE_ORDER:
        reconstructed_mask &= identity[f"pass_{family}"].to_numpy(dtype=bool)
    reconstructed_ids = set(identity.loc[reconstructed_mask, "external_row_id"].astype(str))
    exact_reproduction = frozen_accepted_ids == reconstructed_ids

    enrichment_identity = outputs[ENRICHMENT]
    enrichment_identity = enrichment_identity[enrichment_identity["population"] == "identity_clean"].set_index("metric")
    calibration_identity = outputs[CALIBRATION_SUMMARY]
    calibration_identity = calibration_identity[calibration_identity["population"] == "identity_clean_all"].iloc[0]
    summary = {
        "analysis_version": "ABCB1_EXTERNAL_RISK_EXPERIMENT_V1",
        "analysis_date": "2026-08-05",
        "status": "completed_frozen_model_posthoc_secondary_external_risk_experiments",
        "protocol_path": str(PROTOCOL.relative_to(PROJECT)),
        "protocol_sha256": sha256(PROTOCOL),
        "analysis_populations": {
            "identity_clean_n": int(len(identity)),
            "identity_clean_positive_n": int(identity[Y_COL].sum()),
            "identity_clean_negative_n": int(len(identity) - identity[Y_COL].sum()),
            "scaffold_clean_n": int(len(scaffold)),
        },
        "frozen_rule_reconstruction": {
            "exactly_reproduced": bool(exact_reproduction),
            "frozen_accepted_n": int(len(frozen_accepted_ids)),
            "reconstructed_accepted_n": int(len(reconstructed_ids)),
            "symmetric_difference_n": int(len(frozen_accepted_ids ^ reconstructed_ids)),
            "final_cumulative_accepted_n": int(identity_final["accepted_n"]),
        },
        "identity_clean_operating_points": {
            label: {
                key: value.item() if isinstance(value, np.generic) else value
                for key, value in row.items()
                if key in [
                    "accepted_n", "coverage", "accuracy", "error_rate", "balanced_accuracy", "mcc",
                    "sensitivity", "specificity", "high_risk_precision", "low_risk_npv",
                    "balanced_accuracy_ci95_low", "balanced_accuracy_ci95_high", "mcc_ci95_low", "mcc_ci95_high",
                ]
            }
            for label, row in point_lookup.items()
        },
        "selective_enrichment_identity_clean": {
            metric: {
                key: value.item() if isinstance(value, np.generic) else value
                for key, value in enrichment_identity.loc[metric].to_dict().items()
                if key in ["contrast", "accepted_value", "rejected_forced_value", "favorable_difference", "ci95_low", "ci95_high", "bootstrap_probability_favorable_gt_0"]
            }
            for metric in enrichment_identity.index
        },
        "external_calibration_identity_clean": {
            key: value.item() if isinstance(value, np.generic) else value
            for key, value in calibration_identity.to_dict().items()
            if key in ["n", "observed_prevalence", "mean_predicted_probability", "brier_score", "log_loss", "ece_equal_width_10", "ece_adaptive_10", "calibration_intercept", "calibration_slope", "roc_auc", "average_precision"]
        },
        "interpretation": {
            "allowed": protocol["interpretation_rules"]["allowed"],
            "not_allowed": protocol["interpretation_rules"]["not_allowed"],
            "analysis_status_note": "These are post-hoc secondary analyses performed after the primary external outcome was observed; they are not a new preregistered validation.",
        },
        "output_hashes": {str(path.relative_to(PROJECT)): sha256(path) for path in outputs},
    }
    assert exact_reproduction, summary["frozen_rule_reconstruction"]
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")

    print(json.dumps({
        "analysis_populations": summary["analysis_populations"],
        "frozen_rule_reconstruction": summary["frozen_rule_reconstruction"],
        "identity_clean_operating_points": summary["identity_clean_operating_points"],
        "selective_enrichment_identity_clean": summary["selective_enrichment_identity_clean"],
        "external_calibration_identity_clean": summary["external_calibration_identity_clean"],
        "outputs": [str(path.relative_to(PROJECT)) for path in [*outputs, SUMMARY]],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
