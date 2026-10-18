#!/usr/bin/env python3
"""Run the frozen retrospective ABCB1 structural-shift/OOD audit."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"

PROTOCOL = PROCESSED / "abcb1_ood_structural_shift_protocol_v1.json"
STRUCTURE_LOCK = PROCESSED / "abcb1_ood_structure_assignment_lock_v1.json"
STRUCTURE_ROWS = PROCESSED / "abcb1_ood_structural_assignments_v1.csv"
BIOGEN_ANALYSIS = PROCESSED / "abcb1_quantity_biogen_analysis_rows_v1.csv"

ANALYSIS_ROWS = PROCESSED / "abcb1_ood_analysis_rows_v1.csv"
STRATUM_METRICS = PROCESSED / "abcb1_ood_stratum_metrics_v1.csv"
BOOTSTRAP_SUMMARY = PROCESSED / "abcb1_ood_bootstrap_summary_v1.csv"
CONTRASTS = PROCESSED / "abcb1_ood_contrasts_v1.csv"
ABSTENTION_REASONS = PROCESSED / "abcb1_ood_abstention_reasons_v1.csv"
SUMMARY = PROCESSED / "abcb1_ood_structural_shift_summary_v1.json"

SEED = 20260807
N_BOOTSTRAP = 2000
MIN_SUBSET_TOTAL = 50
MIN_EACH_CLASS = 20
MIN_ACCEPTED_EACH_CLASS = 15
MIN_OOD_COVERAGE = 0.05

HIGH_CALL = "high_ABCB1_efflux_risk"
LOW_CALL = "low_ABCB1_efflux_risk"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_float(function: Callable[[], float]) -> float:
    try:
        value = float(function())
        return value if math.isfinite(value) else math.nan
    except Exception:
        return math.nan


def json_safe(value: Any) -> Any:
    """Recursively replace non-finite/numpy scalars before strict JSON output."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    return value


def confusion_metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    if len(y) == 0:
        return {key: math.nan for key in [
            "balanced_accuracy", "mcc", "sensitivity", "specificity", "ppv", "npv", "accuracy"
        ]}
    y = y.astype(int)
    pred = pred.astype(int)
    tp = int(np.sum((y == 1) & (pred == 1)))
    tn = int(np.sum((y == 0) & (pred == 0)))
    fp = int(np.sum((y == 0) & (pred == 1)))
    fn = int(np.sum((y == 1) & (pred == 0)))
    sensitivity = tp / (tp + fn) if tp + fn else math.nan
    specificity = tn / (tn + fp) if tn + fp else math.nan
    bacc = np.nanmean([sensitivity, specificity]) if not (math.isnan(sensitivity) and math.isnan(specificity)) else math.nan
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / denominator) if denominator else math.nan
    return {
        "balanced_accuracy": float(bacc),
        "mcc": float(mcc),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "ppv": float(tp / (tp + fp)) if tp + fp else math.nan,
        "npv": float(tn / (tn + fn)) if tn + fn else math.nan,
        "accuracy": float((tp + tn) / len(y)),
    }


def ece_equal_width(y: np.ndarray, score: np.ndarray, n_bins: int = 10) -> float:
    if len(y) == 0:
        return math.nan
    clipped = np.clip(score.astype(float), 0.0, 1.0)
    bin_index = np.minimum((clipped * n_bins).astype(int), n_bins - 1)
    ece = 0.0
    for index in range(n_bins):
        mask = bin_index == index
        if not np.any(mask):
            continue
        ece += float(np.mean(mask)) * abs(float(np.mean(clipped[mask])) - float(np.mean(y[mask])))
    return float(ece)


def layer_arrays(frame: pd.DataFrame, layer: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    y_all = frame["frozen_er2_binary_label"].to_numpy(dtype=int)
    score_all = frame["primary_calibrated_probability"].to_numpy(dtype=float)
    forced_all = (score_all >= 0.5).astype(int)
    accepted_mask = frame["risk_call"].isin([HIGH_CALL, LOW_CALL]).to_numpy()
    if layer == "full_forced_binary":
        return y_all, forced_all, score_all
    if layer == "selective_accepted":
        pred = np.where(frame.loc[accepted_mask, "risk_call"].eq(HIGH_CALL), 1, 0).astype(int)
        return y_all[accepted_mask], pred, score_all[accepted_mask]
    if layer == "forced_rejected":
        rejected_mask = ~accepted_mask
        return y_all[rejected_mask], forced_all[rejected_mask], score_all[rejected_mask]
    raise ValueError(layer)


def metrics_for_layer(frame: pd.DataFrame, layer: str) -> dict[str, Any]:
    y, pred, score = layer_arrays(frame, layer)
    confusion = confusion_metrics(y, pred)
    continuous = frame.loc[
        np.ones(len(frame), dtype=bool)
        if layer == "full_forced_binary"
        else (frame["risk_call"].isin([HIGH_CALL, LOW_CALL]).to_numpy() if layer == "selective_accepted" else ~frame["risk_call"].isin([HIGH_CALL, LOW_CALL]).to_numpy()),
        "log10_mdr1_mdck_er",
    ].to_numpy(dtype=float)
    return {
        "n_evaluated": int(len(y)),
        "positive_n": int(np.sum(y == 1)),
        "negative_n": int(np.sum(y == 0)),
        **confusion,
        "roc_auc": safe_float(lambda: roc_auc_score(y, score)) if len(np.unique(y)) == 2 else math.nan,
        "spearman_probability_log10_er": safe_float(lambda: spearmanr(score, continuous).statistic),
        "brier_score": float(np.mean((score - y) ** 2)) if len(y) else math.nan,
        "ece_10_equal_width": ece_equal_width(y, score),
    }


def subset_definitions(frame: pd.DataFrame) -> list[tuple[str, pd.Series, str]]:
    definitions: list[tuple[str, pd.Series, str]] = [
        ("overall_identity_clean", pd.Series(True, index=frame.index), "primary external structural-shift population"),
        ("scaffold_clean", frame["scaffold_clean"].eq("yes"), "Murcko-clean sensitivity only"),
    ]
    for value in ["ID_like", "transition", "OOD", "extreme_OOD"]:
        definitions.append((f"reference_{value}", frame["reference_adaptive_stratum"].eq(value), "reference-LOO D5 stratum"))
    for value in ["Q1_nearest", "Q2", "Q3", "Q4", "Q5_farthest"]:
        definitions.append((f"external_{value}", frame["external_d5_quintile"].eq(value), "outcome-free external D5 quintile"))
    definitions.extend(
        [
            ("strict_LoHi", frame["strict_lohi_ood"].eq("yes"), "S1 joint-model reference <=0.40"),
            ("dual_OOD", frame["dual_ood"].eq("yes"), "beyond primary and direct-model Q80 D5"),
            ("kmeans_OOD20", frame["kmeans_ood20"].eq("yes"), "label-free farthest KMeans clusters >=20%"),
            ("kmeans_non_OOD80", frame["kmeans_ood20"].eq("no"), "remaining KMeans clusters"),
        ]
    )
    return definitions


def adequate_subset(frame: pd.DataFrame) -> bool:
    y = frame["frozen_er2_binary_label"].to_numpy(dtype=int)
    return len(frame) >= MIN_SUBSET_TOTAL and np.sum(y == 0) >= MIN_EACH_CLASS and np.sum(y == 1) >= MIN_EACH_CLASS


def adequate_accepted(frame: pd.DataFrame) -> bool:
    accepted = frame[frame["risk_call"].isin([HIGH_CALL, LOW_CALL])]
    y = accepted["frozen_er2_binary_label"].to_numpy(dtype=int)
    return np.sum(y == 0) >= MIN_ACCEPTED_EACH_CLASS and np.sum(y == 1) >= MIN_ACCEPTED_EACH_CLASS


def make_groups(frame: pd.DataFrame, cluster_column: str) -> list[np.ndarray]:
    return [np.asarray(indices, dtype=int) for indices in frame.groupby(cluster_column, sort=False).indices.values()]


def percentile_interval(values: list[float]) -> tuple[float, float, int]:
    clean = np.asarray([value for value in values if math.isfinite(value)], dtype=float)
    if clean.size == 0:
        return math.nan, math.nan, 0
    return float(np.quantile(clean, 0.025)), float(np.quantile(clean, 0.975)), int(clean.size)


def bootstrap_subset(
    frame: pd.DataFrame,
    subset_id: str,
    cluster_column: str,
    scheme: str,
    seed_offset: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    point_layers = {layer: metrics_for_layer(frame, layer) for layer in ["full_forced_binary", "selective_accepted", "forced_rejected"]}
    coverage_point = float(frame["risk_call"].isin([HIGH_CALL, LOW_CALL]).mean()) if len(frame) else math.nan
    rows: list[dict[str, Any]] = []
    contrast_rows: list[dict[str, Any]] = []
    if not adequate_subset(frame):
        for layer, metrics in point_layers.items():
            for metric in ["balanced_accuracy", "mcc", "spearman_probability_log10_er", "brier_score"]:
                rows.append(
                    {
                        "subset_id": subset_id,
                        "bootstrap_scheme": scheme,
                        "layer": layer,
                        "metric": metric,
                        "point_estimate": metrics[metric],
                        "ci95_low": math.nan,
                        "ci95_high": math.nan,
                        "valid_replicates": 0,
                        "inference_status": "insufficient_subset_or_class_count",
                    }
                )
        rows.append(
            {
                "subset_id": subset_id,
                "bootstrap_scheme": scheme,
                "layer": "selective_framework",
                "metric": "coverage",
                "point_estimate": coverage_point,
                "ci95_low": math.nan,
                "ci95_high": math.nan,
                "valid_replicates": 0,
                "inference_status": "insufficient_subset_or_class_count",
            }
        )
        return rows, contrast_rows

    groups = make_groups(frame, cluster_column)
    rng = np.random.default_rng(SEED + seed_offset)
    values: dict[tuple[str, str], list[float]] = {}
    for layer in ["full_forced_binary", "selective_accepted", "forced_rejected"]:
        for metric in ["balanced_accuracy", "mcc", "spearman_probability_log10_er", "brier_score"]:
            values[(layer, metric)] = []
    coverage_values: list[float] = []
    bacc_contrast: list[float] = []
    mcc_contrast: list[float] = []
    for _ in range(N_BOOTSTRAP):
        selected = rng.integers(0, len(groups), size=len(groups))
        indices = np.concatenate([groups[index] for index in selected])
        sample = frame.iloc[indices]
        sampled_metrics = {layer: metrics_for_layer(sample, layer) for layer in ["full_forced_binary", "selective_accepted", "forced_rejected"]}
        for key in values:
            values[key].append(float(sampled_metrics[key[0]][key[1]]))
        coverage_values.append(float(sample["risk_call"].isin([HIGH_CALL, LOW_CALL]).mean()))
        bacc_contrast.append(
            float(sampled_metrics["selective_accepted"]["balanced_accuracy"])
            - float(sampled_metrics["forced_rejected"]["balanced_accuracy"])
        )
        mcc_contrast.append(
            float(sampled_metrics["selective_accepted"]["mcc"])
            - float(sampled_metrics["forced_rejected"]["mcc"])
        )

    for (layer, metric), metric_values in values.items():
        low, high, valid = percentile_interval(metric_values)
        rows.append(
            {
                "subset_id": subset_id,
                "bootstrap_scheme": scheme,
                "layer": layer,
                "metric": metric,
                "point_estimate": point_layers[layer][metric],
                "ci95_low": low,
                "ci95_high": high,
                "valid_replicates": valid,
                "inference_status": "adequate" if valid >= int(0.9 * N_BOOTSTRAP) else "unstable_replicates",
            }
        )
    low, high, valid = percentile_interval(coverage_values)
    rows.append(
        {
            "subset_id": subset_id,
            "bootstrap_scheme": scheme,
            "layer": "selective_framework",
            "metric": "coverage",
            "point_estimate": coverage_point,
            "ci95_low": low,
            "ci95_high": high,
            "valid_replicates": valid,
            "inference_status": "adequate",
        }
    )
    for metric, point, metric_values in [
        (
            "accepted_minus_forced_rejected_balanced_accuracy",
            point_layers["selective_accepted"]["balanced_accuracy"] - point_layers["forced_rejected"]["balanced_accuracy"],
            bacc_contrast,
        ),
        (
            "accepted_minus_forced_rejected_mcc",
            point_layers["selective_accepted"]["mcc"] - point_layers["forced_rejected"]["mcc"],
            mcc_contrast,
        ),
    ]:
        low, high, valid = percentile_interval(metric_values)
        contrast_rows.append(
            {
                "contrast_id": f"{subset_id}::{metric}",
                "subset_id": subset_id,
                "bootstrap_scheme": scheme,
                "metric": metric,
                "point_estimate": point,
                "ci95_low": low,
                "ci95_high": high,
                "valid_replicates": valid,
                "inference_status": "adequate" if adequate_accepted(frame) else "accepted_class_count_below_confirmatory_minimum",
            }
        )
    return rows, contrast_rows


def paired_between_subset_contrast(
    frame: pd.DataFrame,
    subset_a: pd.Series,
    subset_b: pd.Series,
    contrast_id: str,
    cluster_column: str,
    scheme: str,
    seed_offset: int,
) -> list[dict[str, Any]]:
    a = frame.loc[subset_a]
    b = frame.loc[subset_b]
    metrics_a = metrics_for_layer(a, "full_forced_binary")
    metrics_b = metrics_for_layer(b, "full_forced_binary")
    coverage_a = float(a["risk_call"].isin([HIGH_CALL, LOW_CALL]).mean()) if len(a) else math.nan
    coverage_b = float(b["risk_call"].isin([HIGH_CALL, LOW_CALL]).mean()) if len(b) else math.nan
    adequate = adequate_subset(a) and adequate_subset(b)
    point = {
        "full_balanced_accuracy_a_minus_b": metrics_a["balanced_accuracy"] - metrics_b["balanced_accuracy"],
        "selective_coverage_a_minus_b": coverage_a - coverage_b,
    }
    if not adequate:
        return [
            {
                "contrast_id": f"{contrast_id}::{metric}",
                "subset_id": contrast_id,
                "bootstrap_scheme": scheme,
                "metric": metric,
                "point_estimate": value,
                "ci95_low": math.nan,
                "ci95_high": math.nan,
                "valid_replicates": 0,
                "inference_status": "insufficient_subset_or_class_count",
            }
            for metric, value in point.items()
        ]

    groups = make_groups(frame, cluster_column)
    rng = np.random.default_rng(SEED + seed_offset)
    values = {key: [] for key in point}
    a_col = subset_a.to_numpy()
    b_col = subset_b.to_numpy()
    for _ in range(N_BOOTSTRAP):
        selected = rng.integers(0, len(groups), size=len(groups))
        indices = np.concatenate([groups[index] for index in selected])
        sample = frame.iloc[indices]
        sample_a = sample.loc[a_col[indices]]
        sample_b = sample.loc[b_col[indices]]
        ma = metrics_for_layer(sample_a, "full_forced_binary")
        mb = metrics_for_layer(sample_b, "full_forced_binary")
        values["full_balanced_accuracy_a_minus_b"].append(ma["balanced_accuracy"] - mb["balanced_accuracy"])
        values["selective_coverage_a_minus_b"].append(
            float(sample_a["risk_call"].isin([HIGH_CALL, LOW_CALL]).mean())
            - float(sample_b["risk_call"].isin([HIGH_CALL, LOW_CALL]).mean())
        )
    output: list[dict[str, Any]] = []
    for metric, metric_values in values.items():
        low, high, valid = percentile_interval(metric_values)
        output.append(
            {
                "contrast_id": f"{contrast_id}::{metric}",
                "subset_id": contrast_id,
                "bootstrap_scheme": scheme,
                "metric": metric,
                "point_estimate": point[metric],
                "ci95_low": low,
                "ci95_high": high,
                "valid_replicates": valid,
                "inference_status": "adequate",
            }
        )
    return output


def distance_error_bootstrap(frame: pd.DataFrame, cluster_column: str, scheme: str, seed_offset: int) -> dict[str, Any]:
    d5 = frame["joint_model_d5"].to_numpy(dtype=float)
    y = frame["frozen_er2_binary_label"].to_numpy(dtype=float)
    score = frame["primary_calibrated_probability"].to_numpy(dtype=float)
    loss = (score - y) ** 2
    point = safe_float(lambda: spearmanr(d5, loss).statistic)
    groups = make_groups(frame, cluster_column)
    rng = np.random.default_rng(SEED + seed_offset)
    values: list[float] = []
    for _ in range(N_BOOTSTRAP):
        selected = rng.integers(0, len(groups), size=len(groups))
        indices = np.concatenate([groups[index] for index in selected])
        values.append(safe_float(lambda: spearmanr(d5[indices], loss[indices]).statistic))
    low, high, valid = percentile_interval(values)
    return {
        "contrast_id": "distance_error_spearman",
        "subset_id": "overall_identity_clean",
        "bootstrap_scheme": scheme,
        "metric": "spearman_joint_model_d5_vs_row_brier_loss",
        "point_estimate": point,
        "ci95_low": low,
        "ci95_high": high,
        "valid_replicates": valid,
        "inference_status": "adequate",
    }


def find_bootstrap_row(rows: pd.DataFrame, subset: str, layer: str, metric: str, scheme: str = "murcko_cluster") -> pd.Series | None:
    match = rows[
        rows["subset_id"].eq(subset)
        & rows["layer"].eq(layer)
        & rows["metric"].eq(metric)
        & rows["bootstrap_scheme"].eq(scheme)
    ]
    return None if match.empty else match.iloc[0]


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    lock = json.loads(STRUCTURE_LOCK.read_text(encoding="utf-8"))
    if lock.get("status") != "structure_only_assignments_frozen_before_new_metric_generation":
        raise RuntimeError("structure-only assignment lock is not valid")
    if lock.get("protocol_sha256") != sha256(PROTOCOL):
        raise RuntimeError("protocol hash differs from the structure lock")
    expected_structure_hash = lock["outputs"][str(STRUCTURE_ROWS.relative_to(PROJECT))]
    if sha256(STRUCTURE_ROWS) != expected_structure_hash:
        raise RuntimeError("structural assignment file hash differs from the lock")
    expected_biogen_hash = protocol["input_hashes_sha256"]["数据/processed/abcb1_quantity_biogen_analysis_rows_v1.csv"]
    if sha256(BIOGEN_ANALYSIS) != expected_biogen_hash:
        raise RuntimeError("Biogen analysis input hash differs from the protocol")

    structural = pd.read_csv(STRUCTURE_ROWS, low_memory=False)
    outcomes = pd.read_csv(BIOGEN_ANALYSIS, low_memory=False)
    outcome_columns = [
        "external_row_id",
        "log10_mdr1_mdck_er",
        "mdr1_mdck_er",
        "frozen_er2_binary_label",
        "risk_call",
        "primary_calibrated_probability",
        "abstention_reasons_json",
        "applicability_status",
        "uncertainty_status",
    ]
    outcomes = outcomes.loc[outcomes["identity_clean"].eq("yes"), outcome_columns].copy()
    frame = structural.merge(outcomes, on="external_row_id", how="left", validate="one_to_one")
    if len(frame) != 2624 or frame["frozen_er2_binary_label"].isna().any():
        raise RuntimeError("failed to merge all 2624 identity-clean outcomes")
    frame["frozen_er2_binary_label"] = frame["frozen_er2_binary_label"].astype(int)
    frame["forced_binary_prediction"] = (frame["primary_calibrated_probability"] >= 0.5).astype(int)
    frame["accepted_by_v0_4"] = np.where(frame["risk_call"].isin([HIGH_CALL, LOW_CALL]), "yes", "no")
    frame["row_brier_loss"] = (
        frame["primary_calibrated_probability"] - frame["frozen_er2_binary_label"]
    ) ** 2
    frame.to_csv(ANALYSIS_ROWS, index=False)

    metric_rows: list[dict[str, Any]] = []
    definitions = subset_definitions(frame)
    subset_frames: dict[str, pd.DataFrame] = {}
    for subset_id, mask, role in definitions:
        subset = frame.loc[mask].copy().reset_index(drop=True)
        subset_frames[subset_id] = subset
        accepted_mask = subset["risk_call"].isin([HIGH_CALL, LOW_CALL])
        total_positive = int(subset["frozen_er2_binary_label"].sum())
        total_negative = int(len(subset) - total_positive)
        accepted_positive = int(subset.loc[accepted_mask, "frozen_er2_binary_label"].sum())
        accepted_negative = int(accepted_mask.sum() - accepted_positive)
        coverage = float(accepted_mask.mean()) if len(subset) else math.nan
        positive_coverage = accepted_positive / total_positive if total_positive else math.nan
        negative_coverage = accepted_negative / total_negative if total_negative else math.nan
        for layer in ["full_forced_binary", "selective_accepted", "forced_rejected"]:
            metrics = metrics_for_layer(subset, layer)
            metric_rows.append(
                {
                    "subset_id": subset_id,
                    "subset_role": role,
                    "layer": layer,
                    "n_total": int(len(subset)),
                    "total_positive_n": total_positive,
                    "total_negative_n": total_negative,
                    "accepted_n": int(accepted_mask.sum()),
                    "rejected_n": int((~accepted_mask).sum()),
                    "coverage": coverage,
                    "positive_coverage": positive_coverage,
                    "negative_coverage": negative_coverage,
                    "inference_count_status": "adequate" if adequate_subset(subset) else "insufficient_subset_or_class_count",
                    **metrics,
                }
            )
    metric_frame = pd.DataFrame(metric_rows)
    metric_frame.to_csv(STRATUM_METRICS, index=False)

    bootstrap_rows: list[dict[str, Any]] = []
    contrast_rows: list[dict[str, Any]] = []
    for subset_index, (subset_id, _, _) in enumerate(definitions):
        subset = subset_frames[subset_id]
        for scheme_index, (cluster_column, scheme) in enumerate(
            [("bootstrap_cluster", "murcko_cluster"), ("butina_cluster_id", "butina_ecfp4_cluster")]
        ):
            rows, within = bootstrap_subset(
                subset,
                subset_id,
                cluster_column,
                scheme,
                seed_offset=1000 * scheme_index + 31 * subset_index,
            )
            bootstrap_rows.extend(rows)
            contrast_rows.extend(within)

    for scheme_index, (cluster_column, scheme) in enumerate(
        [("bootstrap_cluster", "murcko_cluster"), ("butina_cluster_id", "butina_ecfp4_cluster")]
    ):
        contrast_rows.append(distance_error_bootstrap(frame, cluster_column, scheme, 8000 + scheme_index))
        contrast_rows.extend(
            paired_between_subset_contrast(
                frame,
                frame["reference_adaptive_stratum"].eq("extreme_OOD"),
                frame["reference_adaptive_stratum"].eq("ID_like"),
                "reference_extreme_minus_ID_like",
                cluster_column,
                scheme,
                9000 + scheme_index,
            )
        )
        contrast_rows.extend(
            paired_between_subset_contrast(
                frame,
                frame["external_d5_quintile"].eq("Q5_farthest"),
                frame["external_d5_quintile"].eq("Q1_nearest"),
                "external_Q5_minus_Q1",
                cluster_column,
                scheme,
                10000 + scheme_index,
            )
        )

    bootstrap_frame = pd.DataFrame(bootstrap_rows)
    bootstrap_frame.to_csv(BOOTSTRAP_SUMMARY, index=False)
    contrast_frame = pd.DataFrame(contrast_rows)
    contrast_frame.to_csv(CONTRASTS, index=False)

    reason_rows: list[dict[str, Any]] = []
    for _, row in frame.loc[frame["accepted_by_v0_4"].eq("no")].iterrows():
        try:
            reasons = json.loads(row.get("abstention_reasons_json", "[]"))
        except Exception:
            reasons = ["unparseable_reason_json"]
        if not reasons:
            reasons = ["unspecified"]
        for reason in reasons:
            reason_rows.append(
                {
                    "external_row_id": row["external_row_id"],
                    "reference_adaptive_stratum": row["reference_adaptive_stratum"],
                    "external_d5_quintile": row["external_d5_quintile"],
                    "strict_lohi_ood": row["strict_lohi_ood"],
                    "dual_ood": row["dual_ood"],
                    "kmeans_ood20": row["kmeans_ood20"],
                    "abstention_reason": reason,
                }
            )
    reason_frame = pd.DataFrame(reason_rows)
    reason_summary = (
        reason_frame.groupby(["reference_adaptive_stratum", "external_d5_quintile", "abstention_reason"], dropna=False)
        .size()
        .reset_index(name="reason_count")
    )
    reason_summary.to_csv(ABSTENTION_REASONS, index=False)

    primary_bootstrap = bootstrap_frame[bootstrap_frame["bootstrap_scheme"].eq("murcko_cluster")]
    primary_contrasts = contrast_frame[contrast_frame["bootstrap_scheme"].eq("murcko_cluster")]
    distance_row = primary_contrasts[primary_contrasts["metric"].eq("spearman_joint_model_d5_vs_row_brier_loss")].iloc[0]
    reference_contrasts = primary_contrasts[primary_contrasts["subset_id"].eq("reference_extreme_minus_ID_like")]
    reference_bacc = reference_contrasts[reference_contrasts["metric"].eq("full_balanced_accuracy_a_minus_b")].iloc[0]
    reference_coverage = reference_contrasts[reference_contrasts["metric"].eq("selective_coverage_a_minus_b")].iloc[0]

    if reference_bacc["inference_status"] != "adequate":
        gate1 = "not_evaluable_ID_like_count_below_frozen_minimum"
    else:
        gate1 = "pass" if distance_row["ci95_low"] > 0 and reference_bacc["ci95_high"] < 0 else "fail"
    if reference_coverage["inference_status"] != "adequate":
        gate2 = "not_evaluable_ID_like_count_below_frozen_minimum"
    else:
        gate2 = "pass" if reference_coverage["ci95_high"] < 0 else "fail"

    pressure_sets = ["reference_extreme_OOD", "strict_LoHi", "dual_OOD", "kmeans_OOD20", "external_Q5_farthest"]
    gate3_candidates: list[dict[str, Any]] = []
    for subset_id in pressure_sets:
        subset = subset_frames[subset_id]
        accepted_bacc = find_bootstrap_row(primary_bootstrap, subset_id, "selective_accepted", "balanced_accuracy")
        accepted_mcc = find_bootstrap_row(primary_bootstrap, subset_id, "selective_accepted", "mcc")
        contrast = primary_contrasts[
            primary_contrasts["subset_id"].eq(subset_id)
            & primary_contrasts["metric"].eq("accepted_minus_forced_rejected_balanced_accuracy")
        ]
        coverage = float(subset["risk_call"].isin([HIGH_CALL, LOW_CALL]).mean()) if len(subset) else 0.0
        adequate = adequate_subset(subset) and adequate_accepted(subset)
        passed = bool(
            adequate
            and coverage >= MIN_OOD_COVERAGE
            and accepted_bacc is not None
            and accepted_mcc is not None
            and not contrast.empty
            and float(accepted_bacc["ci95_low"]) > 0.5
            and float(accepted_mcc["ci95_low"]) > 0.0
            and float(contrast.iloc[0]["ci95_low"]) > 0.0
        )
        gate3_candidates.append(
            {
                "subset_id": subset_id,
                "n": int(len(subset)),
                "coverage": coverage,
                "adequate_counts": bool(adequate),
                "pass": passed,
            }
        )
    gate3 = "pass" if any(row["pass"] for row in gate3_candidates) else "fail_or_not_evaluable"
    if gate1 == "pass" and gate2 == "pass" and gate3 == "pass":
        overall_gate = "go_for_retrospective_OOD_selective_reliability_claim"
    elif gate3 == "pass":
        overall_gate = "partial_support_OOD_accepted_reliability_only_main_distance_gradient_not_evaluable"
    else:
        overall_gate = "no_go_for_OOD_reliability_claim_descriptive_structural_shift_only"

    def metric_anchor(subset_id: str, layer: str) -> dict[str, Any]:
        row = metric_frame[metric_frame["subset_id"].eq(subset_id) & metric_frame["layer"].eq(layer)].iloc[0]
        return {
            "n_total": int(row["n_total"]),
            "n_evaluated": int(row["n_evaluated"]),
            "coverage": float(row["coverage"]),
            "balanced_accuracy": float(row["balanced_accuracy"]),
            "mcc": float(row["mcc"]),
            "spearman_probability_log10_er": float(row["spearman_probability_log10_er"]),
        }

    summary = {
        "analysis_id": "ABCB1_OOD_STRUCTURAL_SHIFT_AUDIT_V1",
        "status": "complete_retrospective_post_outcome_frozen_OOD_audit",
        "protocol_sha256": sha256(PROTOCOL),
        "structure_lock_sha256": sha256(STRUCTURE_LOCK),
        "timing": protocol["analysis_timing"],
        "structural_shift": {
            "identity_clean_n": int(len(frame)),
            "scaffold_clean_n": int(frame["scaffold_clean"].eq("yes").sum()),
            "reference_strata": frame["reference_adaptive_stratum"].value_counts().sort_index().to_dict(),
            "strict_lohi_n": int(frame["strict_lohi_ood"].eq("yes").sum()),
            "dual_ood_n": int(frame["dual_ood"].eq("yes").sum()),
            "kmeans_ood20_n": int(frame["kmeans_ood20"].eq("yes").sum()),
            "joint_model_s1_median": float(frame["joint_model_s1"].median()),
            "joint_model_s5_median": float(frame["joint_model_s5"].median()),
            "joint_model_d5_median": float(frame["joint_model_d5"].median()),
        },
        "metric_anchors": {
            "overall_full": metric_anchor("overall_identity_clean", "full_forced_binary"),
            "overall_accepted": metric_anchor("overall_identity_clean", "selective_accepted"),
            "Q1_full": metric_anchor("external_Q1_nearest", "full_forced_binary"),
            "Q5_full": metric_anchor("external_Q5_farthest", "full_forced_binary"),
            "Q1_accepted": metric_anchor("external_Q1_nearest", "selective_accepted"),
            "Q5_accepted": metric_anchor("external_Q5_farthest", "selective_accepted"),
            "dual_OOD_accepted": metric_anchor("dual_OOD", "selective_accepted"),
            "kmeans_OOD20_accepted": metric_anchor("kmeans_OOD20", "selective_accepted"),
        },
        "gates": {
            "gate1_distance_marks_failure": gate1,
            "gate2_rejector_responds_to_shift": gate2,
            "gate3_OOD_accepted_reliability": gate3,
            "gate3_candidates": gate3_candidates,
            "overall": overall_gate,
        },
        "key_contrasts": {
            "distance_error_spearman": distance_row.to_dict(),
            "reference_extreme_minus_ID_like_full_BAcc": reference_bacc.to_dict(),
            "reference_extreme_minus_ID_like_coverage": reference_coverage.to_dict(),
            "external_Q5_minus_Q1": primary_contrasts[
                primary_contrasts["subset_id"].eq("external_Q5_minus_Q1")
            ].to_dict(orient="records"),
        },
        "outputs": {},
        "interpretation_boundaries": protocol["interpretation_boundaries"],
    }
    for path in [ANALYSIS_ROWS, STRATUM_METRICS, BOOTSTRAP_SUMMARY, CONTRASTS, ABSTENTION_REASONS]:
        summary["outputs"][str(path.relative_to(PROJECT))] = sha256(path)
    safe_summary = json_safe(summary)
    SUMMARY.write_text(json.dumps(safe_summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(safe_summary, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
