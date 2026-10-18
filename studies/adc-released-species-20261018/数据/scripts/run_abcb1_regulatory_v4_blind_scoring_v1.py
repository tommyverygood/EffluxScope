#!/usr/bin/env python3
"""Run the one-time, pre-locked strict-pool v4 ABCB1 evaluation.

The model, threshold, applicability rules, external rows, labels, uncertainty
units, and claim gates are fixed in the protocol lock.  This script refuses to
overwrite an existing prediction output and never refits or recalibrates.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import joblib
import numpy as np
from rdkit import Chem, DataStructs
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    log_loss,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT = SCRIPT_DIR.parents[1]
PROCESSED = PROJECT / "数据/processed"
MODEL_PATH = PROJECT / "模型/ABCB1_screener_v0_1/ABCB1_screener_v0_1.joblib"
POOL_PATH = PROCESSED / "abcb1_regulatory_strict_external_pool_v4.csv"
LOCK_PATH = PROCESSED / "abcb1_regulatory_v4_blind_scoring_protocol_lock_v1.json"
sys.path.insert(0, str(SCRIPT_DIR))

from abcb1_screener_core_v0_1 import FP_GENERATOR, predict_smiles  # noqa: E402


PREDICTIONS_PATH = PROCESSED / "abcb1_regulatory_v4_blind_predictions_v1.csv"
METRICS_PATH = PROCESSED / "abcb1_regulatory_v4_blind_metric_summary_v1.csv"
UNCERTAINTY_PATH = PROCESSED / "abcb1_regulatory_v4_blind_uncertainty_v1.csv"
PAIRED_PATH = PROCESSED / "abcb1_regulatory_v4_blind_paired_comparison_v1.csv"
CLUSTER_PATH = PROCESSED / "abcb1_regulatory_v4_blind_cluster_sensitivity_v1.csv"
ERRORS_PATH = PROCESSED / "abcb1_regulatory_v4_blind_error_audit_v1.csv"
SUMMARY_PATH = PROCESSED / "abcb1_regulatory_v4_blind_summary_v1.json"

OUTPUT_PATHS = [
    PREDICTIONS_PATH,
    METRICS_PATH,
    UNCERTAINTY_PATH,
    PAIRED_PATH,
    CLUSTER_PATH,
    ERRORS_PATH,
    SUMMARY_PATH,
]

METRIC_NAMES = [
    "mcc",
    "balanced_accuracy",
    "sensitivity",
    "specificity",
    "precision",
    "negative_predictive_value",
    "f1",
    "roc_auc",
    "average_precision",
    "brier",
    "log_loss",
    "ece_10bin",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite(value: float) -> float | str:
    return round(float(value), 8) if math.isfinite(float(value)) else ""


def expected_calibration_error(y: np.ndarray, probability: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0.0, 1.0, bins + 1)
    total = len(y)
    error = 0.0
    for idx in range(bins):
        if idx == bins - 1:
            mask = (probability >= edges[idx]) & (probability <= edges[idx + 1])
        else:
            mask = (probability >= edges[idx]) & (probability < edges[idx + 1])
        count = int(mask.sum())
        if not count:
            continue
        error += (count / total) * abs(float(y[mask].mean()) - float(probability[mask].mean()))
    return error


def metric_dict(y: np.ndarray, probability: np.ndarray, threshold: float = 0.5) -> dict:
    prediction = (probability >= threshold).astype(int)
    both = len(set(y.tolist())) == 2
    tn = int(((y == 0) & (prediction == 0)).sum())
    fp = int(((y == 0) & (prediction == 1)).sum())
    fn = int(((y == 1) & (prediction == 0)).sum())
    tp = int(((y == 1) & (prediction == 1)).sum())
    specificity = tn / (tn + fp) if tn + fp else math.nan
    npv = tn / (tn + fn) if tn + fn else math.nan
    clipped = np.clip(probability, 1e-8, 1 - 1e-8)
    return {
        "n": len(y),
        "n_negative": int((y == 0).sum()),
        "n_positive": int((y == 1).sum()),
        "threshold": threshold,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "tp": tp,
        "mcc": float(matthews_corrcoef(y, prediction)) if both else math.nan,
        "balanced_accuracy": float(balanced_accuracy_score(y, prediction)) if both else math.nan,
        "sensitivity": float(recall_score(y, prediction, zero_division=0)),
        "specificity": specificity,
        "precision": float(precision_score(y, prediction, zero_division=0)),
        "negative_predictive_value": npv,
        "f1": float(f1_score(y, prediction, zero_division=0)),
        "roc_auc": float(roc_auc_score(y, probability)) if both else math.nan,
        "average_precision": float(average_precision_score(y, probability)) if both else math.nan,
        "brier": float(brier_score_loss(y, probability)),
        "log_loss": float(log_loss(y, clipped, labels=[0, 1])),
        "ece_10bin": expected_calibration_error(y, probability),
    }


def training_baseline(package: dict) -> tuple[list[int], list]:
    audit = read_csv(PROCESSED / "abcb1_screener_v0_1_training_audit.csv")
    if len(audit) != package["training_n"]:
        raise RuntimeError("training audit length mismatch")
    keys = [row.get("derived_standardized_inchikey", "") for row in audit]
    if keys != package["training_structure_keys"]:
        raise RuntimeError("training audit order/key mismatch")
    return [int(row["target_label"]) for row in audit], package["training_fingerprints"]


def similarity_baselines(mol: Chem.Mol, training_labels: list[int], training_fps: list) -> dict:
    fp = FP_GENERATOR.GetFingerprint(mol)
    similarities = np.asarray(DataStructs.BulkTanimotoSimilarity(fp, training_fps), dtype=float)
    order = np.argsort(-similarities)
    nearest = int(order[0])
    top = order[:5]
    weights = np.maximum(similarities[top], 1e-8)
    labels = np.asarray(training_labels, dtype=float)
    return {
        "nearest_training_index": nearest,
        "nearest_training_tanimoto": round(float(similarities[nearest]), 6),
        "one_nn_probability": float(labels[nearest]),
        "knn5_weighted_probability": float(np.average(labels[top], weights=weights)),
    }


def entity_bootstrap(
    y: np.ndarray,
    probabilities: dict[str, np.ndarray],
    iterations: int,
    seed: int,
) -> tuple[dict[str, dict[str, tuple[float, float]]], int]:
    negative = np.flatnonzero(y == 0)
    positive = np.flatnonzero(y == 1)
    rng = np.random.default_rng(seed)
    storage: dict[str, dict[str, list[float]]] = {
        model: defaultdict(list) for model in probabilities
    }
    for _ in range(iterations):
        sampled = np.concatenate(
            [
                rng.choice(negative, size=len(negative), replace=True),
                rng.choice(positive, size=len(positive), replace=True),
            ]
        )
        for model, probability in probabilities.items():
            metrics = metric_dict(y[sampled], probability[sampled])
            for name in METRIC_NAMES:
                if math.isfinite(metrics[name]):
                    storage[model][name].append(metrics[name])
    intervals = {
        model: {
            name: (float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)))
            for name, values in metrics.items()
            if values
        }
        for model, metrics in storage.items()
    }
    return intervals, iterations


def cluster_bootstrap(
    y: np.ndarray,
    clusters: np.ndarray,
    probabilities: dict[str, np.ndarray],
    iterations: int,
    seed: int,
) -> tuple[dict[str, dict[str, tuple[float, float]]], int]:
    unique_clusters = sorted(set(clusters.tolist()))
    cluster_indices = {cluster: np.flatnonzero(clusters == cluster) for cluster in unique_clusters}
    rng = np.random.default_rng(seed)
    storage: dict[str, dict[str, list[float]]] = {
        model: defaultdict(list) for model in probabilities
    }
    usable = 0
    for _ in range(iterations):
        chosen = rng.choice(unique_clusters, size=len(unique_clusters), replace=True)
        sampled = np.concatenate([cluster_indices[cluster] for cluster in chosen])
        if len(set(y[sampled].tolist())) < 2:
            continue
        usable += 1
        for model, probability in probabilities.items():
            metrics = metric_dict(y[sampled], probability[sampled])
            for name in METRIC_NAMES:
                if math.isfinite(metrics[name]):
                    storage[model][name].append(metrics[name])
    intervals = {
        model: {
            name: (float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)))
            for name, values in metrics.items()
            if values
        }
        for model, metrics in storage.items()
    }
    return intervals, usable


def paired_comparison(
    y: np.ndarray,
    clusters: np.ndarray,
    probabilities: dict[str, np.ndarray],
    entity_iterations: int,
    cluster_iterations: int,
    entity_seed: int,
    cluster_seed: int,
) -> list[dict]:
    frozen = probabilities["frozen_ABCB1_v0_1"]
    point_frozen = metric_dict(y, frozen)
    output = []
    negative = np.flatnonzero(y == 0)
    positive = np.flatnonzero(y == 1)
    unique_clusters = sorted(set(clusters.tolist()))
    cluster_indices = {cluster: np.flatnonzero(clusters == cluster) for cluster in unique_clusters}
    metric_names = ["mcc", "balanced_accuracy", "roc_auc", "brier"]
    for comparator in ["one_nearest_neighbor", "five_nearest_neighbor_weighted", "training_class_prior"]:
        point_comparator = metric_dict(y, probabilities[comparator])
        for unit, iterations, seed in [
            ("entity_stratified", entity_iterations, entity_seed + 17),
            ("upstream_experiment_cluster", cluster_iterations, cluster_seed + 17),
        ]:
            rng = np.random.default_rng(seed)
            storage = defaultdict(list)
            usable = 0
            for _ in range(iterations):
                if unit == "entity_stratified":
                    sampled = np.concatenate(
                        [
                            rng.choice(negative, size=len(negative), replace=True),
                            rng.choice(positive, size=len(positive), replace=True),
                        ]
                    )
                else:
                    chosen = rng.choice(unique_clusters, size=len(unique_clusters), replace=True)
                    sampled = np.concatenate([cluster_indices[cluster] for cluster in chosen])
                    if len(set(y[sampled].tolist())) < 2:
                        continue
                usable += 1
                frozen_metrics = metric_dict(y[sampled], frozen[sampled])
                comparator_metrics = metric_dict(y[sampled], probabilities[comparator][sampled])
                for name in metric_names:
                    difference = frozen_metrics[name] - comparator_metrics[name]
                    if math.isfinite(difference):
                        storage[name].append(difference)
            for name in metric_names:
                values = storage[name]
                output.append(
                    {
                        "comparison": f"frozen_ABCB1_v0_1_minus_{comparator}",
                        "resampling_unit": unit,
                        "metric": name,
                        "frozen_value": finite(point_frozen[name]),
                        "comparator_value": finite(point_comparator[name]),
                        "paired_difference": finite(point_frozen[name] - point_comparator[name]),
                        "ci95_lower": finite(float(np.quantile(values, 0.025))),
                        "ci95_upper": finite(float(np.quantile(values, 0.975))),
                        "requested_iterations": iterations,
                        "usable_iterations": usable,
                        "direction_note": "positive favors frozen" if name != "brier" else "negative favors frozen because lower Brier is better",
                        "claim_boundary": "pre-locked paired diagnostic; no model selection or tuning",
                    }
                )
    return output


def main() -> None:
    existing = [str(path) for path in OUTPUT_PATHS if path.exists()]
    if existing:
        raise RuntimeError(f"one-time output already exists; refusing overwrite: {existing}")

    protocol = json.loads(LOCK_PATH.read_text(encoding="utf-8"))
    if protocol["lock_status"] != "locked_after_label_and_identity_freeze_before_any_v4_model_scoring":
        raise RuntimeError("protocol is not in the expected pre-score locked state")
    if sha256_file(MODEL_PATH) != protocol["frozen_model"]["sha256"]:
        raise RuntimeError("frozen model hash mismatch")
    if sha256_file(POOL_PATH) != protocol["frozen_external_pool"]["sha256"]:
        raise RuntimeError("strict pool v4 hash mismatch")
    for item in protocol["supporting_locks"]:
        if sha256_file(PROJECT / item["path"]) != item["sha256"]:
            raise RuntimeError(f"supporting lock hash mismatch: {item['path']}")
    for item in [protocol["blindness_audit"]["audit_csv"], protocol["blindness_audit"]["audit_summary"]]:
        if sha256_file(PROJECT / item["path"]) != item["sha256"]:
            raise RuntimeError(f"blindness audit hash mismatch: {item['path']}")

    package = joblib.load(MODEL_PATH)
    if float(package.get("binary_threshold", 0.5)) != 0.5:
        raise RuntimeError("frozen package threshold is not 0.5")
    training_labels, training_fps = training_baseline(package)
    prior_probability = package["training_n_positive"] / package["training_n"]
    pool = read_csv(POOL_PATH)
    if len(pool) != protocol["frozen_external_pool"]["n"]:
        raise RuntimeError("strict pool row count differs from lock")

    predictions = []
    for row in pool:
        prediction = predict_smiles(row["standardized_isomeric_smiles"], package)
        if prediction["standardization_status"] not in {"ok", "matched_frozen_training_structure"}:
            raise RuntimeError(f"prediction standardization failed: {row['entity_name']}")
        if prediction["standardized_inchikey"].split("-")[0] != row["standardized_inchikey"].split("-")[0]:
            raise RuntimeError(f"prediction identity drift: {row['entity_name']}")
        mol = Chem.MolFromSmiles(prediction["standardized_isomeric_smiles"])
        if mol is None:
            raise RuntimeError(f"standardized molecule unresolved: {row['entity_name']}")
        baselines = similarity_baselines(mol, training_labels, training_fps)
        label = int(row["directional_label"])
        binary = int(prediction["binary_prediction"])
        if label == 0 and binary == 1:
            error_type = "false_positive_substrate_call"
        elif label == 1 and binary == 0:
            error_type = "false_negative_non_substrate_call"
        else:
            error_type = "correct"
        predictions.append(
            {
                **row,
                "frozen_model_sha256": protocol["frozen_model"]["sha256"],
                "frozen_probability": prediction["calibrated_probability"],
                "frozen_raw_probability": prediction["raw_probability"],
                "frozen_binary_prediction": binary,
                "fixed_threshold": 0.5,
                "prediction_correct": "yes" if binary == label else "no",
                "error_type": error_type,
                "max_training_tanimoto": prediction["max_training_morgan_tanimoto"],
                "applicability_band": prediction["applicability_band"],
                "screening_call": prediction["screening_call"],
                "training_class_prior_probability": prior_probability,
                **baselines,
                "blind_scoring_role": "primary_all_33_prelocked",
                "resampling_cluster": row["upstream_experiment_cluster"] or row["source_id"],
                "result_use_boundary": "ordinary-drug direct-assay external evaluation; not ADC release-species validation",
            }
        )

    y = np.asarray([int(row["directional_label"]) for row in predictions], dtype=int)
    clusters = np.asarray([row["resampling_cluster"] for row in predictions])
    probabilities = {
        "frozen_ABCB1_v0_1": np.asarray([float(row["frozen_probability"]) for row in predictions]),
        "training_class_prior": np.asarray([float(row["training_class_prior_probability"]) for row in predictions]),
        "one_nearest_neighbor": np.asarray([float(row["one_nn_probability"]) for row in predictions]),
        "five_nearest_neighbor_weighted": np.asarray([float(row["knn5_weighted_probability"]) for row in predictions]),
    }

    entity_config = protocol["uncertainty"]["entity_bootstrap"]
    cluster_config = protocol["uncertainty"]["cluster_bootstrap"]
    entity_intervals, entity_usable = entity_bootstrap(
        y,
        probabilities,
        int(entity_config["iterations"]),
        int(entity_config["seed"]),
    )
    cluster_intervals, cluster_usable = cluster_bootstrap(
        y,
        clusters,
        probabilities,
        int(cluster_config["iterations"]),
        int(cluster_config["seed"]),
    )

    uncertainty_rows = []
    for model, probability in probabilities.items():
        point = metric_dict(y, probability)
        for name in METRIC_NAMES:
            entity_lower, entity_upper = entity_intervals.get(model, {}).get(name, (math.nan, math.nan))
            cluster_lower, cluster_upper = cluster_intervals.get(model, {}).get(name, (math.nan, math.nan))
            uncertainty_rows.append(
                {
                    "benchmark": "strict_regulatory_primary_study_pool_v4",
                    "subset": "all_33_prelocked",
                    "model": model,
                    "metric": name,
                    "value": finite(point[name]),
                    "entity_ci95_lower": finite(entity_lower),
                    "entity_ci95_upper": finite(entity_upper),
                    "entity_requested_iterations": entity_config["iterations"],
                    "entity_usable_iterations": entity_usable,
                    "cluster_ci95_lower": finite(cluster_lower),
                    "cluster_ci95_upper": finite(cluster_upper),
                    "cluster_requested_iterations": cluster_config["iterations"],
                    "cluster_usable_iterations": cluster_usable,
                    "n": point["n"],
                    "n_negative": point["n_negative"],
                    "n_positive": point["n_positive"],
                    "tn": point["tn"],
                    "fp": point["fp"],
                    "fn": point["fn"],
                    "tp": point["tp"],
                    "threshold": 0.5,
                    "claim_boundary": "prelocked blind external metric; cluster interval is primary dependence sensitivity",
                }
            )

    subset_definitions = {
        "all_33_prelocked": predictions,
        "within_or_caution_nonabstain": [row for row in predictions if row["applicability_band"] != "outside_empirical_support"],
        "within_empirical_support": [row for row in predictions if row["applicability_band"] == "within_empirical_support"],
        "low_empirical_support": [row for row in predictions if row["applicability_band"] == "low_empirical_support"],
        "outside_empirical_support": [row for row in predictions if row["applicability_band"] == "outside_empirical_support"],
        "max_training_tanimoto_below_0_5": [row for row in predictions if float(row["max_training_tanimoto"]) < 0.5],
        "max_training_tanimoto_at_least_0_5": [row for row in predictions if float(row["max_training_tanimoto"]) >= 0.5],
    }
    metric_rows = []
    for subset_name, subset in subset_definitions.items():
        if not subset:
            continue
        subset_y = np.asarray([int(row["directional_label"]) for row in subset], dtype=int)
        for model, field in [
            ("frozen_ABCB1_v0_1", "frozen_probability"),
            ("training_class_prior", "training_class_prior_probability"),
            ("one_nearest_neighbor", "one_nn_probability"),
            ("five_nearest_neighbor_weighted", "knn5_weighted_probability"),
        ]:
            probability = np.asarray([float(row[field]) for row in subset])
            metrics = metric_dict(subset_y, probability)
            for name in METRIC_NAMES:
                metric_rows.append(
                    {
                        "benchmark": "strict_regulatory_primary_study_pool_v4",
                        "subset": subset_name,
                        "model": model,
                        "metric": name,
                        "value": finite(metrics[name]),
                        "n": metrics["n"],
                        "n_negative": metrics["n_negative"],
                        "n_positive": metrics["n_positive"],
                        "tn": metrics["tn"],
                        "fp": metrics["fp"],
                        "fn": metrics["fn"],
                        "tp": metrics["tp"],
                        "threshold": 0.5,
                        "uncertainty_location": "see blind_uncertainty_v1 for all_33; secondary subsets are point diagnostics",
                        "claim_boundary": "prelocked external result; applicability subsets are secondary",
                    }
                )

    cluster_rows = []
    unique_clusters = sorted(set(clusters.tolist()))
    for cluster in unique_clusters:
        for mode, subset in [
            ("cluster_only", [row for row in predictions if row["resampling_cluster"] == cluster]),
            ("leave_one_cluster_out", [row for row in predictions if row["resampling_cluster"] != cluster]),
        ]:
            subset_y = np.asarray([int(row["directional_label"]) for row in subset], dtype=int)
            probability = np.asarray([float(row["frozen_probability"]) for row in subset])
            metrics = metric_dict(subset_y, probability)
            cluster_rows.append(
                {
                    "mode": mode,
                    "upstream_experiment_cluster": cluster,
                    "cluster_size_in_full_pool": int((clusters == cluster).sum()),
                    **{key: finite(value) if isinstance(value, float) else value for key, value in metrics.items()},
                    "claim_boundary": "cluster sensitivity; single-class cluster-only metrics are descriptive",
                }
            )

    paired_rows = paired_comparison(
        y,
        clusters,
        probabilities,
        int(entity_config["iterations"]),
        int(cluster_config["iterations"]),
        int(entity_config["seed"]),
        int(cluster_config["seed"]),
    )

    errors = [row for row in predictions if row["prediction_correct"] == "no"]
    errors.sort(key=lambda row: (row["error_type"], -abs(float(row["frozen_probability"]) - 0.5)))

    frozen_point = metric_dict(y, probabilities["frozen_ABCB1_v0_1"])
    frozen_uncertainty = {
        row["metric"]: row
        for row in uncertainty_rows
        if row["model"] == "frozen_ABCB1_v0_1"
    }
    leave_out_bacc = [
        float(row["balanced_accuracy"])
        for row in cluster_rows
        if row["mode"] == "leave_one_cluster_out" and row["balanced_accuracy"] != ""
    ]
    nonabstain_n = len(subset_definitions["within_or_caution_nonabstain"])
    nonabstain_coverage = nonabstain_n / len(predictions)
    gate_rows = [
        {
            "gate": "entity_MCC_lower_95CI_gt_0",
            "value": frozen_uncertainty["mcc"]["entity_ci95_lower"],
            "rule": ">0",
            "status": "pass" if float(frozen_uncertainty["mcc"]["entity_ci95_lower"]) > 0 else "fail",
        },
        {
            "gate": "entity_BAcc_lower_95CI_gt_0_5",
            "value": frozen_uncertainty["balanced_accuracy"]["entity_ci95_lower"],
            "rule": ">0.5",
            "status": "pass" if float(frozen_uncertainty["balanced_accuracy"]["entity_ci95_lower"]) > 0.5 else "fail",
        },
        {
            "gate": "cluster_MCC_lower_95CI_gt_0",
            "value": frozen_uncertainty["mcc"]["cluster_ci95_lower"],
            "rule": ">0",
            "status": "pass" if float(frozen_uncertainty["mcc"]["cluster_ci95_lower"]) > 0 else "fail",
        },
        {
            "gate": "cluster_BAcc_lower_95CI_gt_0_5",
            "value": frozen_uncertainty["balanced_accuracy"]["cluster_ci95_lower"],
            "rule": ">0.5",
            "status": "pass" if float(frozen_uncertainty["balanced_accuracy"]["cluster_ci95_lower"]) > 0.5 else "fail",
        },
        {
            "gate": "entity_ROC_AUC_lower_95CI_gt_0_5",
            "value": frozen_uncertainty["roc_auc"]["entity_ci95_lower"],
            "rule": ">0.5",
            "status": "pass" if float(frozen_uncertainty["roc_auc"]["entity_ci95_lower"]) > 0.5 else "fail",
        },
        {
            "gate": "leave_one_cluster_out_min_BAcc_gt_0_5",
            "value": finite(min(leave_out_bacc)),
            "rule": ">0.5",
            "status": "pass" if min(leave_out_bacc) > 0.5 else "fail",
        },
        {
            "gate": "nonabstain_coverage_at_least_0_5",
            "value": finite(nonabstain_coverage),
            "rule": ">=0.5",
            "status": "pass" if nonabstain_coverage >= 0.5 else "fail",
        },
        {
            "gate": "sensitivity_and_specificity_at_least_0_5",
            "value": f"sensitivity={frozen_point['sensitivity']:.8f};specificity={frozen_point['specificity']:.8f}",
            "rule": "both>=0.5",
            "status": "pass" if frozen_point["sensitivity"] >= 0.5 and frozen_point["specificity"] >= 0.5 else "fail",
        },
    ]
    validated = all(row["status"] == "pass" for row in gate_rows)

    write_csv(PREDICTIONS_PATH, predictions)
    write_csv(METRICS_PATH, metric_rows)
    write_csv(UNCERTAINTY_PATH, uncertainty_rows)
    write_csv(PAIRED_PATH, paired_rows)
    write_csv(CLUSTER_PATH, cluster_rows)
    write_csv(ERRORS_PATH, errors, fieldnames=list(predictions[0]))

    summary = {
        "run_date": "2026-08-03",
        "protocol_lock_path": str(LOCK_PATH),
        "protocol_lock_sha256": sha256_file(LOCK_PATH),
        "one_time_run_status": "completed_without_refit_recalibration_or_threshold_change",
        "frozen_model_sha256": sha256_file(MODEL_PATH),
        "strict_pool_v4_sha256": sha256_file(POOL_PATH),
        "n": len(predictions),
        "class_counts": dict(Counter(row["directional_label"] for row in predictions)),
        "upstream_experiment_clusters": len(unique_clusters),
        "multi_entity_clusters": {
            cluster: int(count)
            for cluster, count in Counter(clusters.tolist()).items()
            if count > 1
        },
        "applicability_counts": dict(Counter(row["applicability_band"] for row in predictions)),
        "nonabstain_n": nonabstain_n,
        "nonabstain_coverage": finite(nonabstain_coverage),
        "frozen_primary_metrics": {
            key: finite(value) if isinstance(value, float) else value
            for key, value in frozen_point.items()
        },
        "frozen_primary_uncertainty": {
            metric: {
                "value": row["value"],
                "entity_ci95": [row["entity_ci95_lower"], row["entity_ci95_upper"]],
                "cluster_ci95": [row["cluster_ci95_lower"], row["cluster_ci95_upper"]],
            }
            for metric, row in frozen_uncertainty.items()
        },
        "error_counts": dict(Counter(row["error_type"] for row in predictions)),
        "claim_gates": gate_rows,
        "validated_binary_predictor_claim_status": "pass" if validated else "fail",
        "action": (
            "retain frozen model and seek untouched ADC release-species validation"
            if validated
            else "do_not_tune_on_v4; retain bounded screener and use frozen errors to design assay-aware v0.2"
        ),
        "claim_boundary": "ordinary-drug direct-assay blind external evaluation; not intact-ADC or ADC release-species validation",
        "output_hashes": {},
    }
    for path in [PREDICTIONS_PATH, METRICS_PATH, UNCERTAINTY_PATH, PAIRED_PATH, CLUSTER_PATH, ERRORS_PATH]:
        summary["output_hashes"][str(path.relative_to(PROJECT))] = sha256_file(path)
    SUMMARY_PATH.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
