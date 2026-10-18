#!/usr/bin/env python3
"""Independently validate the frozen ABCB1 structural-shift/OOD audit.

This script intentionally does not import the analysis runner.  It re-derives
structural counts, subset membership, point metrics, bootstrap completeness,
output hashes, and gate decisions from the frozen artifacts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from datetime import datetime, timezone
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
ANALYSIS_ROWS = PROCESSED / "abcb1_ood_analysis_rows_v1.csv"
STRATUM_METRICS = PROCESSED / "abcb1_ood_stratum_metrics_v1.csv"
BOOTSTRAP_SUMMARY = PROCESSED / "abcb1_ood_bootstrap_summary_v1.csv"
CONTRASTS = PROCESSED / "abcb1_ood_contrasts_v1.csv"
ABSTENTION_REASONS = PROCESSED / "abcb1_ood_abstention_reasons_v1.csv"
SUMMARY = PROCESSED / "abcb1_ood_structural_shift_summary_v1.json"
VALIDATION = PROCESSED / "abcb1_ood_structural_shift_validation_v1.json"

RESULT_FILES = [
    ANALYSIS_ROWS,
    STRATUM_METRICS,
    BOOTSTRAP_SUMMARY,
    CONTRASTS,
    ABSTENTION_REASONS,
    SUMMARY,
]

HIGH_CALL = "high_ABCB1_efflux_risk"
LOW_CALL = "low_ABCB1_efflux_risk"
LAYERS = ["full_forced_binary", "selective_accepted", "forced_rejected"]
BOOTSTRAP_SCHEMES = ["murcko_cluster", "butina_ecfp4_cluster"]
PRESSURE_SETS = [
    "reference_extreme_OOD",
    "strict_LoHi",
    "dual_OOD",
    "kmeans_OOD20",
    "external_Q5_farthest",
]

EXPECTED_STRUCTURE_COUNTS = {
    "identity_clean": 2624,
    "scaffold_clean": 2562,
    "reference_adaptive_stratum": {
        "ID_like": 7,
        "transition": 1173,
        "OOD": 1249,
        "extreme_OOD": 195,
    },
    "strict_lohi_ood": 2581,
    "dual_ood": 609,
    "kmeans_ood20": 544,
    "kmeans_ood20_clusters": [3, 15],
    "external_d5_quintile": {
        "Q1_nearest": 525,
        "Q2": 525,
        "Q3": 524,
        "Q4": 525,
        "Q5_farthest": 525,
    },
}

REQUIRED_COLUMNS = {
    "structure_rows": {
        "external_row_id",
        "standardized_inchikey",
        "connectivity_block",
        "murcko_scaffold_smiles",
        "bootstrap_cluster",
        "identity_clean",
        "scaffold_clean",
        "primary443_d5",
        "direct63_d5",
        "joint_model_s1",
        "joint_model_d5",
        "reference_adaptive_stratum",
        "strict_lohi_ood",
        "dual_ood",
        "external_d5_quintile",
        "kmeans_cluster_id",
        "kmeans_ood20",
        "butina_cluster_id",
    },
    "analysis_rows": {
        "external_row_id",
        "frozen_er2_binary_label",
        "log10_mdr1_mdck_er",
        "risk_call",
        "primary_calibrated_probability",
        "forced_binary_prediction",
        "accepted_by_v0_4",
        "row_brier_loss",
        "scaffold_clean",
        "reference_adaptive_stratum",
        "strict_lohi_ood",
        "dual_ood",
        "external_d5_quintile",
        "kmeans_ood20",
        "joint_model_d5",
    },
    "stratum_metrics": {
        "subset_id",
        "subset_role",
        "layer",
        "n_total",
        "total_positive_n",
        "total_negative_n",
        "accepted_n",
        "rejected_n",
        "coverage",
        "n_evaluated",
        "positive_n",
        "negative_n",
        "balanced_accuracy",
        "mcc",
        "roc_auc",
        "spearman_probability_log10_er",
        "brier_score",
        "ece_10_equal_width",
    },
    "bootstrap_summary": {
        "subset_id",
        "bootstrap_scheme",
        "layer",
        "metric",
        "point_estimate",
        "ci95_low",
        "ci95_high",
        "valid_replicates",
        "inference_status",
    },
    "contrasts": {
        "contrast_id",
        "subset_id",
        "bootstrap_scheme",
        "metric",
        "point_estimate",
        "ci95_low",
        "ci95_high",
        "valid_replicates",
        "inference_status",
    },
    "abstention_reasons": {
        "reference_adaptive_stratum",
        "external_d5_quintile",
        "abstention_reason",
        "reason_count",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_or_none(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): finite_or_none(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_or_none(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def same_number(left: Any, right: Any, tolerance: float = 1e-10) -> bool:
    try:
        a = float(left)
        b = float(right)
    except (TypeError, ValueError):
        return left == right
    if math.isnan(a) and math.isnan(b):
        return True
    return math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def safe_stat(function: Callable[[], float]) -> float:
    try:
        value = float(function())
        return value if math.isfinite(value) else math.nan
    except Exception:
        return math.nan


def confusion_metrics(y: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    y = np.asarray(y, dtype=int)
    prediction = np.asarray(prediction, dtype=int)
    if y.size == 0:
        return {name: math.nan for name in [
            "balanced_accuracy", "mcc", "sensitivity", "specificity", "ppv", "npv", "accuracy"
        ]}
    tp = int(np.sum((y == 1) & (prediction == 1)))
    tn = int(np.sum((y == 0) & (prediction == 0)))
    fp = int(np.sum((y == 0) & (prediction == 1)))
    fn = int(np.sum((y == 1) & (prediction == 0)))
    sensitivity = tp / (tp + fn) if tp + fn else math.nan
    specificity = tn / (tn + fp) if tn + fp else math.nan
    if math.isnan(sensitivity) and math.isnan(specificity):
        balanced_accuracy = math.nan
    else:
        balanced_accuracy = float(np.nanmean([sensitivity, specificity]))
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {
        "balanced_accuracy": balanced_accuracy,
        "mcc": (tp * tn - fp * fn) / denominator if denominator else math.nan,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "ppv": tp / (tp + fp) if tp + fp else math.nan,
        "npv": tn / (tn + fn) if tn + fn else math.nan,
        "accuracy": (tp + tn) / len(y),
    }


def ece_equal_width(y: np.ndarray, score: np.ndarray, bins: int = 10) -> float:
    y = np.asarray(y, dtype=float)
    score = np.clip(np.asarray(score, dtype=float), 0.0, 1.0)
    if y.size == 0:
        return math.nan
    bin_index = np.minimum((score * bins).astype(int), bins - 1)
    total = 0.0
    for index in range(bins):
        mask = bin_index == index
        if np.any(mask):
            total += float(np.mean(mask)) * abs(float(np.mean(score[mask])) - float(np.mean(y[mask])))
    return total


def accepted_mask(frame: pd.DataFrame) -> np.ndarray:
    return frame["risk_call"].isin([HIGH_CALL, LOW_CALL]).to_numpy()


def independent_layer_metrics(frame: pd.DataFrame, layer: str) -> dict[str, Any]:
    y_all = frame["frozen_er2_binary_label"].to_numpy(dtype=int)
    score_all = frame["primary_calibrated_probability"].to_numpy(dtype=float)
    forced_all = (score_all >= 0.5).astype(int)
    accepted = accepted_mask(frame)
    if layer == "full_forced_binary":
        mask = np.ones(len(frame), dtype=bool)
        prediction = forced_all
    elif layer == "selective_accepted":
        mask = accepted
        prediction = np.where(frame.loc[mask, "risk_call"].eq(HIGH_CALL), 1, 0).astype(int)
    elif layer == "forced_rejected":
        mask = ~accepted
        prediction = forced_all[mask]
    else:
        raise ValueError(f"unknown layer: {layer}")
    y = y_all[mask]
    score = score_all[mask]
    continuous = frame.loc[mask, "log10_mdr1_mdck_er"].to_numpy(dtype=float)
    confusion = confusion_metrics(y, prediction)
    return {
        "n_evaluated": int(len(y)),
        "positive_n": int(np.sum(y == 1)),
        "negative_n": int(np.sum(y == 0)),
        **confusion,
        "roc_auc": safe_stat(lambda: roc_auc_score(y, score)) if len(np.unique(y)) == 2 else math.nan,
        "spearman_probability_log10_er": safe_stat(lambda: spearmanr(score, continuous).statistic),
        "brier_score": float(np.mean((score - y) ** 2)) if len(y) else math.nan,
        "ece_10_equal_width": ece_equal_width(y, score),
    }


def subset_masks(frame: pd.DataFrame) -> dict[str, pd.Series]:
    masks: dict[str, pd.Series] = {
        "overall_identity_clean": pd.Series(True, index=frame.index),
        "scaffold_clean": frame["scaffold_clean"].eq("yes"),
    }
    for value in ["ID_like", "transition", "OOD", "extreme_OOD"]:
        masks[f"reference_{value}"] = frame["reference_adaptive_stratum"].eq(value)
    for value in ["Q1_nearest", "Q2", "Q3", "Q4", "Q5_farthest"]:
        masks[f"external_{value}"] = frame["external_d5_quintile"].eq(value)
    masks.update(
        {
            "strict_LoHi": frame["strict_lohi_ood"].eq("yes"),
            "dual_OOD": frame["dual_ood"].eq("yes"),
            "kmeans_OOD20": frame["kmeans_ood20"].eq("yes"),
            "kmeans_non_OOD80": frame["kmeans_ood20"].eq("no"),
        }
    )
    return masks


def wait_for_results(wait_seconds: int, poll_seconds: int) -> list[str]:
    deadline = time.monotonic() + max(0, wait_seconds)
    while True:
        missing = [str(path.relative_to(PROJECT)) for path in RESULT_FILES if not path.exists()]
        if not missing or time.monotonic() >= deadline:
            return missing
        time.sleep(max(1, poll_seconds))


class Recorder:
    def __init__(self) -> None:
        self.checks: list[dict[str, Any]] = []

    def add(self, check_id: str, passed: bool, observed: Any = None, expected: Any = None, detail: str = "") -> None:
        self.checks.append(
            {
                "check_id": check_id,
                "status": "pass" if passed else "fail",
                "observed": finite_or_none(observed),
                "expected": finite_or_none(expected),
                "detail": detail,
            }
        )

    def warn(self, check_id: str, observed: Any, detail: str) -> None:
        self.checks.append(
            {
                "check_id": check_id,
                "status": "warning",
                "observed": finite_or_none(observed),
                "expected": None,
                "detail": detail,
            }
        )


def adequate_subset(frame: pd.DataFrame, protocol: dict[str, Any]) -> bool:
    minimums = protocol["minimum_inference_counts"]
    y = frame["frozen_er2_binary_label"].to_numpy(dtype=int)
    return bool(
        len(frame) >= int(minimums["subset_total"])
        and np.sum(y == 0) >= int(minimums["each_true_class"])
        and np.sum(y == 1) >= int(minimums["each_true_class"])
    )


def adequate_accepted(frame: pd.DataFrame, protocol: dict[str, Any]) -> bool:
    minimum = int(protocol["minimum_inference_counts"]["accepted_each_true_class"])
    selected = frame.loc[accepted_mask(frame)]
    y = selected["frozen_er2_binary_label"].to_numpy(dtype=int)
    return bool(np.sum(y == 0) >= minimum and np.sum(y == 1) >= minimum)


def selected_row(frame: pd.DataFrame, **conditions: str) -> pd.Series | None:
    mask = pd.Series(True, index=frame.index)
    for column, value in conditions.items():
        mask &= frame[column].eq(value)
    match = frame.loc[mask]
    return None if match.empty else match.iloc[0]


def expected_gates(
    analysis: pd.DataFrame,
    bootstrap: pd.DataFrame,
    contrasts: pd.DataFrame,
    protocol: dict[str, Any],
) -> dict[str, Any]:
    primary_bootstrap = bootstrap.loc[bootstrap["bootstrap_scheme"].eq("murcko_cluster")]
    primary_contrasts = contrasts.loc[contrasts["bootstrap_scheme"].eq("murcko_cluster")]
    distance = selected_row(
        primary_contrasts,
        metric="spearman_joint_model_d5_vs_row_brier_loss",
    )
    reference_bacc = selected_row(
        primary_contrasts,
        subset_id="reference_extreme_minus_ID_like",
        metric="full_balanced_accuracy_a_minus_b",
    )
    reference_coverage = selected_row(
        primary_contrasts,
        subset_id="reference_extreme_minus_ID_like",
        metric="selective_coverage_a_minus_b",
    )
    if distance is None or reference_bacc is None or reference_bacc["inference_status"] != "adequate":
        gate1 = "not_evaluable_ID_like_count_below_frozen_minimum"
    else:
        gate1 = "pass" if float(distance["ci95_low"]) > 0 and float(reference_bacc["ci95_high"]) < 0 else "fail"
    if reference_coverage is None or reference_coverage["inference_status"] != "adequate":
        gate2 = "not_evaluable_ID_like_count_below_frozen_minimum"
    else:
        gate2 = "pass" if float(reference_coverage["ci95_high"]) < 0 else "fail"

    masks = subset_masks(analysis)
    minimum_coverage = float(protocol["minimum_inference_counts"]["minimum_selective_coverage_for_OOD_claim"])
    candidates: list[dict[str, Any]] = []
    for subset_id in PRESSURE_SETS:
        subset = analysis.loc[masks[subset_id]]
        coverage = float(np.mean(accepted_mask(subset))) if len(subset) else 0.0
        counts_ok = adequate_subset(subset, protocol) and adequate_accepted(subset, protocol)
        accepted_bacc = selected_row(
            primary_bootstrap,
            subset_id=subset_id,
            layer="selective_accepted",
            metric="balanced_accuracy",
        )
        accepted_mcc = selected_row(
            primary_bootstrap,
            subset_id=subset_id,
            layer="selective_accepted",
            metric="mcc",
        )
        contrast = selected_row(
            primary_contrasts,
            subset_id=subset_id,
            metric="accepted_minus_forced_rejected_balanced_accuracy",
        )
        passed = bool(
            counts_ok
            and coverage >= minimum_coverage
            and accepted_bacc is not None
            and accepted_mcc is not None
            and contrast is not None
            and float(accepted_bacc["ci95_low"]) > 0.5
            and float(accepted_mcc["ci95_low"]) > 0.0
            and float(contrast["ci95_low"]) > 0.0
        )
        candidates.append(
            {
                "subset_id": subset_id,
                "n": int(len(subset)),
                "coverage": coverage,
                "adequate_counts": bool(counts_ok),
                "pass": passed,
            }
        )
    gate3 = "pass" if any(row["pass"] for row in candidates) else "fail_or_not_evaluable"
    if gate1 == "pass" and gate2 == "pass" and gate3 == "pass":
        overall = "go_for_retrospective_OOD_selective_reliability_claim"
    elif gate3 == "pass":
        overall = "partial_support_OOD_accepted_reliability_only_main_distance_gradient_not_evaluable"
    else:
        overall = "no_go_for_OOD_reliability_claim_descriptive_structural_shift_only"
    return {
        "gate1_distance_marks_failure": gate1,
        "gate2_rejector_responds_to_shift": gate2,
        "gate3_OOD_accepted_reliability": gate3,
        "gate3_candidates": candidates,
        "overall": overall,
    }


def validate() -> dict[str, Any]:
    recorder = Recorder()
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    lock = json.loads(STRUCTURE_LOCK.read_text(encoding="utf-8"))
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))

    protocol_hash = sha256(PROTOCOL)
    structure_lock_hash = sha256(STRUCTURE_LOCK)
    recorder.add(
        "protocol_hash_matches_structure_lock",
        lock.get("protocol_sha256") == protocol_hash,
        lock.get("protocol_sha256"),
        protocol_hash,
    )
    recorder.add(
        "protocol_hash_matches_summary",
        summary.get("protocol_sha256") == protocol_hash,
        summary.get("protocol_sha256"),
        protocol_hash,
    )
    recorder.add(
        "structure_lock_hash_matches_summary",
        summary.get("structure_lock_sha256") == structure_lock_hash,
        summary.get("structure_lock_sha256"),
        structure_lock_hash,
    )
    recorder.add(
        "structure_lock_status",
        lock.get("status") == "structure_only_assignments_frozen_before_new_metric_generation",
        lock.get("status"),
        "structure_only_assignments_frozen_before_new_metric_generation",
    )
    structure_relative = str(STRUCTURE_ROWS.relative_to(PROJECT))
    recorder.add(
        "structure_rows_hash_matches_lock",
        lock.get("outputs", {}).get(structure_relative) == sha256(STRUCTURE_ROWS),
        sha256(STRUCTURE_ROWS),
        lock.get("outputs", {}).get(structure_relative),
    )

    structure = pd.read_csv(STRUCTURE_ROWS, low_memory=False)
    analysis = pd.read_csv(ANALYSIS_ROWS, low_memory=False)
    metrics = pd.read_csv(STRATUM_METRICS, low_memory=False)
    bootstrap = pd.read_csv(BOOTSTRAP_SUMMARY, low_memory=False)
    contrasts = pd.read_csv(CONTRASTS, low_memory=False)
    reasons = pd.read_csv(ABSTENTION_REASONS, low_memory=False)
    tables = {
        "structure_rows": structure,
        "analysis_rows": analysis,
        "stratum_metrics": metrics,
        "bootstrap_summary": bootstrap,
        "contrasts": contrasts,
        "abstention_reasons": reasons,
    }
    for name, frame in tables.items():
        missing = sorted(REQUIRED_COLUMNS[name] - set(frame.columns))
        recorder.add(f"{name}_required_columns", not missing, missing, [])

    recorder.add("structure_row_count", len(structure) == 2624, len(structure), 2624)
    recorder.add("analysis_row_count", len(analysis) == 2624, len(analysis), 2624)
    recorder.add(
        "structure_external_row_id_unique",
        structure["external_row_id"].is_unique and structure["external_row_id"].notna().all(),
        int(structure["external_row_id"].nunique()),
        2624,
    )
    recorder.add(
        "analysis_external_row_id_unique",
        analysis["external_row_id"].is_unique and analysis["external_row_id"].notna().all(),
        int(analysis["external_row_id"].nunique()),
        2624,
    )
    recorder.add(
        "analysis_ids_match_structure_ids",
        set(analysis["external_row_id"]) == set(structure["external_row_id"]),
        len(set(analysis["external_row_id"]) ^ set(structure["external_row_id"])),
        0,
    )

    observed_structure_counts = {
        "identity_clean": int(structure["identity_clean"].eq("yes").sum()),
        "scaffold_clean": int(structure["scaffold_clean"].eq("yes").sum()),
        "reference_adaptive_stratum": {
            key: int(value)
            for key, value in structure["reference_adaptive_stratum"].value_counts().to_dict().items()
        },
        "strict_lohi_ood": int(structure["strict_lohi_ood"].eq("yes").sum()),
        "dual_ood": int(structure["dual_ood"].eq("yes").sum()),
        "kmeans_ood20": int(structure["kmeans_ood20"].eq("yes").sum()),
        "kmeans_ood20_clusters": sorted(
            int(value) for value in structure.loc[structure["kmeans_ood20"].eq("yes"), "kmeans_cluster_id"].unique()
        ),
        "external_d5_quintile": {
            key: int(value)
            for key, value in structure["external_d5_quintile"].value_counts().to_dict().items()
        },
    }
    recorder.add(
        "frozen_structural_counts",
        observed_structure_counts == EXPECTED_STRUCTURE_COUNTS,
        observed_structure_counts,
        EXPECTED_STRUCTURE_COUNTS,
    )
    recorder.add(
        "structural_counts_match_lock",
        lock.get("external_counts", {}).get("identity_clean") == observed_structure_counts["identity_clean"]
        and lock.get("external_counts", {}).get("scaffold_clean") == observed_structure_counts["scaffold_clean"]
        and lock.get("external_counts", {}).get("reference_adaptive_strata")
        == observed_structure_counts["reference_adaptive_stratum"]
        and lock.get("external_counts", {}).get("strict_lohi_ood") == observed_structure_counts["strict_lohi_ood"]
        and lock.get("external_counts", {}).get("dual_ood") == observed_structure_counts["dual_ood"]
        and lock.get("external_counts", {}).get("kmeans_ood20") == observed_structure_counts["kmeans_ood20"],
        observed_structure_counts,
        lock.get("external_counts"),
    )

    structural_compare_columns = [
        "scaffold_clean",
        "reference_adaptive_stratum",
        "strict_lohi_ood",
        "dual_ood",
        "external_d5_quintile",
        "kmeans_ood20",
        "joint_model_d5",
    ]
    joined = analysis[["external_row_id", *structural_compare_columns]].merge(
        structure[["external_row_id", *structural_compare_columns]],
        on="external_row_id",
        suffixes=("_analysis", "_structure"),
        validate="one_to_one",
    )
    mismatches = 0
    for column in structural_compare_columns:
        left = joined[f"{column}_analysis"]
        right = joined[f"{column}_structure"]
        if pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right):
            mismatches += int((~np.isclose(
                left.to_numpy(float),
                right.to_numpy(float),
                rtol=0,
                atol=1e-12,
                equal_nan=True,
            )).sum())
        else:
            mismatches += int((left.fillna("<NA>") != right.fillna("<NA>")).sum())
    recorder.add("analysis_preserves_frozen_structure_assignments", mismatches == 0, mismatches, 0)

    forced_expected = (analysis["primary_calibrated_probability"].to_numpy(float) >= 0.5).astype(int)
    accepted_expected = np.where(accepted_mask(analysis), "yes", "no")
    brier_expected = (
        analysis["primary_calibrated_probability"].to_numpy(float)
        - analysis["frozen_er2_binary_label"].to_numpy(int)
    ) ** 2
    recorder.add(
        "derived_analysis_columns",
        np.array_equal(forced_expected, analysis["forced_binary_prediction"].to_numpy(int))
        and np.array_equal(accepted_expected, analysis["accepted_by_v0_4"].to_numpy(str))
        and np.allclose(brier_expected, analysis["row_brier_loss"].to_numpy(float), atol=1e-12, rtol=0),
        {
            "forced_prediction_mismatches": int(np.sum(forced_expected != analysis["forced_binary_prediction"].to_numpy(int))),
            "accepted_flag_mismatches": int(np.sum(accepted_expected != analysis["accepted_by_v0_4"].to_numpy(str))),
            "max_brier_difference": float(np.max(np.abs(brier_expected - analysis["row_brier_loss"].to_numpy(float)))),
        },
        {"all_mismatches": 0},
    )

    masks = subset_masks(analysis)
    expected_subset_ids = set(masks)
    recorder.add(
        "metrics_composite_key_unique",
        not metrics.duplicated(["subset_id", "layer"]).any(),
        int(metrics.duplicated(["subset_id", "layer"]).sum()),
        0,
    )
    recorder.add(
        "metrics_expected_rows",
        len(metrics) == len(expected_subset_ids) * len(LAYERS)
        and set(metrics["subset_id"]) == expected_subset_ids
        and set(metrics["layer"]) == set(LAYERS),
        {"rows": len(metrics), "subset_ids": sorted(set(metrics["subset_id"])), "layers": sorted(set(metrics["layer"]))},
        {"rows": 45, "subset_ids": sorted(expected_subset_ids), "layers": sorted(LAYERS)},
    )

    metric_fields = [
        "n_evaluated",
        "positive_n",
        "negative_n",
        "balanced_accuracy",
        "mcc",
        "sensitivity",
        "specificity",
        "ppv",
        "npv",
        "accuracy",
        "roc_auc",
        "spearman_probability_log10_er",
        "brier_score",
        "ece_10_equal_width",
    ]
    metric_mismatches: list[dict[str, Any]] = []
    subset_count_mismatches: list[dict[str, Any]] = []
    recomputed_anchors: dict[str, Any] = {}
    for subset_id, mask in masks.items():
        subset = analysis.loc[mask]
        accepted = accepted_mask(subset)
        total_positive = int(subset["frozen_er2_binary_label"].sum())
        total_negative = int(len(subset) - total_positive)
        expected_counts = {
            "n_total": int(len(subset)),
            "total_positive_n": total_positive,
            "total_negative_n": total_negative,
            "accepted_n": int(np.sum(accepted)),
            "rejected_n": int(np.sum(~accepted)),
            "coverage": float(np.mean(accepted)) if len(subset) else math.nan,
        }
        subset_rows = metrics.loc[metrics["subset_id"].eq(subset_id)]
        for field, expected in expected_counts.items():
            values = subset_rows[field].tolist()
            if len(values) != len(LAYERS) or not all(same_number(value, expected) for value in values):
                subset_count_mismatches.append(
                    {"subset_id": subset_id, "field": field, "observed": values, "expected": expected}
                )
        recomputed_anchors[subset_id] = {}
        for layer in LAYERS:
            row = selected_row(metrics, subset_id=subset_id, layer=layer)
            independent = independent_layer_metrics(subset, layer)
            recomputed_anchors[subset_id][layer] = independent
            if row is None:
                metric_mismatches.append({"subset_id": subset_id, "layer": layer, "field": "row", "issue": "missing"})
                continue
            for field in metric_fields:
                if not same_number(row[field], independent[field]):
                    metric_mismatches.append(
                        {
                            "subset_id": subset_id,
                            "layer": layer,
                            "field": field,
                            "observed": row[field],
                            "expected": independent[field],
                        }
                    )
    recorder.add("subset_counts_consistent", not subset_count_mismatches, subset_count_mismatches[:20], [])
    recorder.add("point_metrics_independently_recomputed", not metric_mismatches, metric_mismatches[:20], [])

    expected_bootstrap_rows = len(expected_subset_ids) * len(BOOTSTRAP_SCHEMES) * (len(LAYERS) * 4 + 1)
    recorder.add(
        "bootstrap_composite_key_unique",
        not bootstrap.duplicated(["subset_id", "bootstrap_scheme", "layer", "metric"]).any(),
        int(bootstrap.duplicated(["subset_id", "bootstrap_scheme", "layer", "metric"]).sum()),
        0,
    )
    recorder.add(
        "bootstrap_expected_design_rows",
        len(bootstrap) == expected_bootstrap_rows
        and set(bootstrap["subset_id"]) == expected_subset_ids
        and set(bootstrap["bootstrap_scheme"]) == set(BOOTSTRAP_SCHEMES),
        {"rows": len(bootstrap), "subsets": len(set(bootstrap["subset_id"])), "schemes": sorted(set(bootstrap["bootstrap_scheme"]))},
        {"rows": expected_bootstrap_rows, "subsets": len(expected_subset_ids), "schemes": BOOTSTRAP_SCHEMES},
    )
    recorder.add(
        "contrast_composite_key_unique",
        not contrasts.duplicated(["contrast_id", "bootstrap_scheme"]).any(),
        int(contrasts.duplicated(["contrast_id", "bootstrap_scheme"]).sum()),
        0,
    )
    recorder.add(
        "abstention_reason_key_unique",
        not reasons.duplicated(["reference_adaptive_stratum", "external_d5_quintile", "abstention_reason"]).any(),
        int(reasons.duplicated(["reference_adaptive_stratum", "external_d5_quintile", "abstention_reason"]).sum()),
        0,
    )

    point_mismatch_bootstrap: list[dict[str, Any]] = []
    for _, row in bootstrap.iterrows():
        subset_id = str(row["subset_id"])
        layer = str(row["layer"])
        metric = str(row["metric"])
        if layer == "selective_framework" and metric == "coverage":
            expected = float(np.mean(accepted_mask(analysis.loc[masks[subset_id]])))
        else:
            expected = recomputed_anchors[subset_id][layer][metric]
        if not same_number(row["point_estimate"], expected):
            point_mismatch_bootstrap.append(
                {
                    "subset_id": subset_id,
                    "scheme": row["bootstrap_scheme"],
                    "layer": layer,
                    "metric": metric,
                    "observed": row["point_estimate"],
                    "expected": expected,
                }
            )
    recorder.add("bootstrap_point_estimates_match_analysis", not point_mismatch_bootstrap, point_mismatch_bootstrap[:20], [])

    adequate_bootstrap = bootstrap.loc[bootstrap["inference_status"].eq("adequate")]
    unstable_bootstrap = bootstrap.loc[bootstrap["inference_status"].eq("unstable_replicates")]
    adequate_contrasts = contrasts.loc[contrasts["inference_status"].eq("adequate")]
    bootstrap_valid_min = int(adequate_bootstrap["valid_replicates"].min()) if len(adequate_bootstrap) else 0
    contrast_valid_min = int(adequate_contrasts["valid_replicates"].min()) if len(adequate_contrasts) else 0
    recorder.add(
        "adequate_bootstrap_rows_have_at_least_1800_valid_replicates",
        len(adequate_bootstrap) > 0 and bootstrap_valid_min >= 1800,
        {
            "adequate_rows": int(len(adequate_bootstrap)),
            "minimum_valid_replicates": bootstrap_valid_min,
        },
        {"minimum_valid_replicates": 1800},
    )
    if not unstable_bootstrap.empty:
        recorder.warn(
            "bootstrap_rows_explicitly_marked_unstable",
            unstable_bootstrap[[
                "subset_id",
                "bootstrap_scheme",
                "layer",
                "metric",
                "valid_replicates",
            ]].to_dict(orient="records"),
            (
                "These rows are not treated as adequate inference. In this run they arise from the "
                "five accepted extreme-OOD rows and must not support a selective-reliability claim."
            ),
        )
    recorder.add(
        "adequate_contrast_rows_have_at_least_1800_valid_replicates",
        len(adequate_contrasts) > 0 and contrast_valid_min >= 1800,
        {"adequate_rows": int(len(adequate_contrasts)), "minimum_valid_replicates": contrast_valid_min},
        {"minimum_valid_replicates": 1800},
    )
    invalid_replicate_counts = int(((bootstrap["valid_replicates"] < 0) | (bootstrap["valid_replicates"] > 2000)).sum())
    invalid_replicate_counts += int(((contrasts["valid_replicates"] < 0) | (contrasts["valid_replicates"] > 2000)).sum())
    recorder.add("replicate_counts_within_frozen_range", invalid_replicate_counts == 0, invalid_replicate_counts, 0)

    summary_hash_mismatches: list[dict[str, str | None]] = []
    for path in [ANALYSIS_ROWS, STRATUM_METRICS, BOOTSTRAP_SUMMARY, CONTRASTS, ABSTENTION_REASONS]:
        relative = str(path.relative_to(PROJECT))
        observed = sha256(path)
        expected = summary.get("outputs", {}).get(relative)
        if observed != expected:
            summary_hash_mismatches.append({"path": relative, "observed": observed, "expected": expected})
    recorder.add("summary_output_hashes", not summary_hash_mismatches, summary_hash_mismatches, [])

    summary_shift = summary.get("structural_shift", {})
    summary_shift_ok = (
        summary_shift.get("identity_clean_n") == EXPECTED_STRUCTURE_COUNTS["identity_clean"]
        and summary_shift.get("scaffold_clean_n") == EXPECTED_STRUCTURE_COUNTS["scaffold_clean"]
        and summary_shift.get("reference_strata") == EXPECTED_STRUCTURE_COUNTS["reference_adaptive_stratum"]
        and summary_shift.get("strict_lohi_n") == EXPECTED_STRUCTURE_COUNTS["strict_lohi_ood"]
        and summary_shift.get("dual_ood_n") == EXPECTED_STRUCTURE_COUNTS["dual_ood"]
        and summary_shift.get("kmeans_ood20_n") == EXPECTED_STRUCTURE_COUNTS["kmeans_ood20"]
    )
    recorder.add("summary_structural_counts", summary_shift_ok, summary_shift, EXPECTED_STRUCTURE_COUNTS)

    anchor_map = {
        "overall_full": ("overall_identity_clean", "full_forced_binary"),
        "overall_accepted": ("overall_identity_clean", "selective_accepted"),
        "Q1_full": ("external_Q1_nearest", "full_forced_binary"),
        "Q5_full": ("external_Q5_farthest", "full_forced_binary"),
        "Q1_accepted": ("external_Q1_nearest", "selective_accepted"),
        "Q5_accepted": ("external_Q5_farthest", "selective_accepted"),
        "dual_OOD_accepted": ("dual_OOD", "selective_accepted"),
        "kmeans_OOD20_accepted": ("kmeans_OOD20", "selective_accepted"),
    }
    anchor_mismatches: list[dict[str, Any]] = []
    for anchor_name, (subset_id, layer) in anchor_map.items():
        observed = summary.get("metric_anchors", {}).get(anchor_name, {})
        metric_row = selected_row(metrics, subset_id=subset_id, layer=layer)
        if metric_row is None:
            anchor_mismatches.append({"anchor": anchor_name, "issue": "metrics row missing"})
            continue
        expected = {
            "n_total": int(metric_row["n_total"]),
            "n_evaluated": int(metric_row["n_evaluated"]),
            "coverage": float(metric_row["coverage"]),
            "balanced_accuracy": float(metric_row["balanced_accuracy"]),
            "mcc": float(metric_row["mcc"]),
            "spearman_probability_log10_er": float(metric_row["spearman_probability_log10_er"]),
        }
        if set(observed) != set(expected) or any(not same_number(observed.get(key), value) for key, value in expected.items()):
            anchor_mismatches.append({"anchor": anchor_name, "observed": observed, "expected": expected})
    recorder.add("summary_metric_anchors", not anchor_mismatches, anchor_mismatches, [])

    gates_expected = expected_gates(analysis, bootstrap, contrasts, protocol)
    gates_observed = summary.get("gates", {})
    gate_scalar_fields = [
        "gate1_distance_marks_failure",
        "gate2_rejector_responds_to_shift",
        "gate3_OOD_accepted_reliability",
        "overall",
    ]
    scalar_ok = all(gates_observed.get(field) == gates_expected[field] for field in gate_scalar_fields)
    candidates_observed = gates_observed.get("gate3_candidates", [])
    candidates_ok = len(candidates_observed) == len(gates_expected["gate3_candidates"])
    if candidates_ok:
        for observed, expected in zip(candidates_observed, gates_expected["gate3_candidates"]):
            for field in ["subset_id", "n", "adequate_counts", "pass"]:
                if observed.get(field) != expected[field]:
                    candidates_ok = False
            if not same_number(observed.get("coverage"), expected["coverage"]):
                candidates_ok = False
    recorder.add(
        "summary_gates_independently_rederived",
        scalar_ok and candidates_ok,
        gates_observed,
        gates_expected,
    )

    failures = [check for check in recorder.checks if check["status"] == "fail"]
    warnings = [check for check in recorder.checks if check["status"] == "warning"]
    result = {
        "validation_id": "ABCB1_OOD_STRUCTURAL_SHIFT_VALIDATION_V1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "pass" if not failures else "fail",
        "independence": "does_not_import_or_execute_the_primary_analysis_runner",
        "validated_analysis_id": summary.get("analysis_id"),
        "checks_total": len(recorder.checks),
        "checks_passed": sum(check["status"] == "pass" for check in recorder.checks),
        "checks_failed": len(failures),
        "warnings": len(warnings),
        "key_observations": {
            "identity_clean_n": len(analysis),
            "scaffold_clean_n": int(analysis["scaffold_clean"].eq("yes").sum()),
            "reference_ID_like_n": int(analysis["reference_adaptive_stratum"].eq("ID_like").sum()),
            "adequate_bootstrap_rows": int(len(adequate_bootstrap)),
            "minimum_valid_bootstrap_replicates": bootstrap_valid_min,
            "adequate_contrast_rows": int(len(adequate_contrasts)),
            "minimum_valid_contrast_replicates": contrast_valid_min,
            "independently_rederived_gates": gates_expected,
        },
        "checks": recorder.checks,
        "validated_file_hashes": {
            str(path.relative_to(PROJECT)): sha256(path)
            for path in [PROTOCOL, STRUCTURE_LOCK, STRUCTURE_ROWS, *RESULT_FILES]
        },
        "interpretation_boundary": (
            "Validation confirms computational integrity and frozen gate logic only; "
            "it does not convert the retrospective post-outcome audit into prospective OOD or ADC validation."
        ),
    }
    return finite_or_none(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wait-seconds", type=int, default=0, help="wait for primary result files before validating")
    parser.add_argument("--poll-seconds", type=int, default=10, help="poll interval while waiting")
    args = parser.parse_args()

    prerequisites = [PROTOCOL, STRUCTURE_LOCK, STRUCTURE_ROWS]
    missing_prerequisites = [str(path.relative_to(PROJECT)) for path in prerequisites if not path.exists()]
    if missing_prerequisites:
        raise SystemExit(f"Missing frozen prerequisites: {missing_prerequisites}")
    missing_results = wait_for_results(args.wait_seconds, args.poll_seconds)
    if missing_results:
        incomplete = {
            "validation_id": "ABCB1_OOD_STRUCTURAL_SHIFT_VALIDATION_V1",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "status": "incomplete_primary_outputs_missing",
            "missing_results": missing_results,
        }
        VALIDATION.write_text(json.dumps(incomplete, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(incomplete, ensure_ascii=False, indent=2))
        raise SystemExit(2)

    result = validate()
    VALIDATION.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "checks_total": result["checks_total"],
        "checks_passed": result["checks_passed"],
        "checks_failed": result["checks_failed"],
        "warnings": result["warnings"],
        "output": str(VALIDATION.relative_to(PROJECT)),
        "failed_check_ids": [check["check_id"] for check in result["checks"] if check["status"] == "fail"],
    }, ensure_ascii=False, indent=2))
    raise SystemExit(0 if result["status"] == "pass" else 1)


if __name__ == "__main__":
    main()
