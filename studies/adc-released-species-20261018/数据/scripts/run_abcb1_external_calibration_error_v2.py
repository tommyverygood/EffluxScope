#!/usr/bin/env python3
"""Deep external calibration and error diagnostics for frozen ABCB1 scores."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, roc_auc_score
from sklearn.model_selection import GroupKFold

from run_abcb1_external_risk_experiments_v1 import (
    CLUSTER_COL,
    HIGH,
    LOW,
    N_BOOT,
    P_COL,
    PROJECT,
    PROCESSED,
    SEED,
    Y_COL,
    clean_frame,
    metrics_for_frame,
    round_frame,
    wilson_interval,
)


INPUT = PROCESSED / "abcb1_quantity_biogen_analysis_rows_v1.csv"
MODEL = PROJECT / "模型" / "ABCB1_selective_screener_v0_4_dev" / "ABCB1_selective_uncertainty_v0_4_dev.joblib"
CHANNEL_PROTOCOL = PROCESSED / "abcb1_external_evidence_channel_ablation_protocol_v1.json"
PROTOCOL = PROCESSED / "abcb1_external_calibration_error_protocol_v2.json"

PREDICTIONS = PROCESSED / "abcb1_external_crossfitted_calibration_predictions_v2.csv"
CALIBRATION_METRICS = PROCESSED / "abcb1_external_calibration_metrics_v2.csv"
CALIBRATION_BOOTSTRAP = PROCESSED / "abcb1_external_calibration_bootstrap_v2.csv"
CALIBRATION_PAIRWISE = PROCESSED / "abcb1_external_calibration_pairwise_bootstrap_v2.csv"
CALIBRATION_BINS = PROCESSED / "abcb1_external_calibration_bins_v2.csv"
BRIER_DECOMPOSITION = PROCESSED / "abcb1_external_brier_decomposition_v2.csv"
ERROR_COUNTS = PROCESSED / "abcb1_external_error_counts_v2.csv"
ERROR_FEATURES = PROCESSED / "abcb1_external_error_feature_comparisons_v2.csv"
ERROR_STRATA = PROCESSED / "abcb1_external_error_strata_v2.csv"
SCAFFOLD_HOTSPOTS = PROCESSED / "abcb1_external_error_scaffold_hotspots_v2.csv"
SUMMARY = PROCESSED / "abcb1_external_calibration_error_summary_v2.json"

METHOD_COLUMNS = {
    "raw_frozen_score": "raw_probability",
    "cross_fitted_platt": "platt_probability",
    "cross_fitted_isotonic": "isotonic_probability",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fit_calibration_line(y: np.ndarray, p: np.ndarray) -> tuple[float, float]:
    p = np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)
    y = np.asarray(y, dtype=float)
    x = np.column_stack([np.ones(len(p)), np.log(p / (1.0 - p))])
    beta = np.zeros(2, dtype=float)
    for _ in range(50):
        linear = np.clip(x @ beta, -30, 30)
        mu = 1.0 / (1.0 + np.exp(-linear))
        weights = np.clip(mu * (1.0 - mu), 1e-8, None)
        gradient = x.T @ (y - mu)
        hessian = x.T @ (weights[:, None] * x)
        try:
            step = np.linalg.solve(hessian, gradient)
        except np.linalg.LinAlgError:
            return np.nan, np.nan
        beta += step
        if np.max(np.abs(step)) < 1e-9:
            break
    return float(beta[0]), float(beta[1])


def bin_assignments(p: np.ndarray, strategy: str, n_bins: int = 10) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    if strategy == "equal_width":
        edges = np.linspace(0, 1, n_bins + 1)
    elif strategy == "adaptive":
        edges = np.unique(np.quantile(p, np.linspace(0, 1, n_bins + 1)))
        if len(edges) < 2:
            return np.zeros(len(p), dtype=int)
        edges[0] = min(edges[0], 0.0)
        edges[-1] = max(edges[-1], 1.0)
    else:
        raise ValueError(strategy)
    return np.clip(np.digitize(p, edges[1:-1], right=True), 0, len(edges) - 2)


def expected_calibration_error(y: np.ndarray, p: np.ndarray, strategy: str) -> float:
    bins = bin_assignments(p, strategy)
    total = len(y)
    value = 0.0
    for bin_id in np.unique(bins):
        mask = bins == bin_id
        value += mask.sum() / total * abs(float(np.mean(p[mask])) - float(np.mean(y[mask])))
    return float(value)


def calibration_metric_values(y: np.ndarray, p: np.ndarray) -> dict[str, float]:
    p = np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)
    y = np.asarray(y, dtype=int)
    intercept, slope = fit_calibration_line(y, p)
    return {
        "brier_score": float(brier_score_loss(y, p)),
        "log_loss": float(log_loss(y, p, labels=[0, 1])),
        "ece_equal_width_10": expected_calibration_error(y, p, "equal_width"),
        "ece_adaptive_10": expected_calibration_error(y, p, "adaptive"),
        "calibration_intercept": intercept,
        "calibration_slope": slope,
        "roc_auc": float(roc_auc_score(y, p)),
        "average_precision": float(average_precision_score(y, p)),
    }


def cross_fitted_predictions(work: pd.DataFrame, population: str, folds: int) -> pd.DataFrame:
    frame = work.reset_index(drop=True).copy()
    y = frame[Y_COL].to_numpy(dtype=int)
    raw_p = np.clip(frame[P_COL].to_numpy(dtype=float), 1e-6, 1 - 1e-6)
    groups = frame[CLUSTER_COL].astype(str).to_numpy()
    platt = np.full(len(frame), np.nan)
    isotonic = np.full(len(frame), np.nan)
    fold_id = np.full(len(frame), -1, dtype=int)
    splitter = GroupKFold(n_splits=folds, shuffle=True, random_state=SEED)
    for fold, (train, test) in enumerate(splitter.split(raw_p, y, groups=groups), start=1):
        train_logit = np.log(raw_p[train] / (1.0 - raw_p[train])).reshape(-1, 1)
        test_logit = np.log(raw_p[test] / (1.0 - raw_p[test])).reshape(-1, 1)
        platt_model = LogisticRegression(C=np.inf, l1_ratio=0.0, solver="lbfgs", max_iter=5000)
        platt_model.fit(train_logit, y[train])
        platt[test] = platt_model.predict_proba(test_logit)[:, 1]
        iso_model = IsotonicRegression(out_of_bounds="clip", y_min=1e-6, y_max=1 - 1e-6)
        iso_model.fit(raw_p[train], y[train])
        isotonic[test] = iso_model.predict(raw_p[test])
        fold_id[test] = fold
    assert np.isfinite(platt).all() and np.isfinite(isotonic).all() and (fold_id > 0).all()
    output = pd.DataFrame({
        "population": population,
        "external_row_id": frame["external_row_id"].astype(str),
        "standardized_inchikey": frame["standardized_inchikey"].astype(str),
        "bootstrap_cluster": groups,
        "fold_id": fold_id,
        "frozen_er2_binary_label": y,
        "raw_probability": raw_p,
        "platt_probability": np.clip(platt, 1e-6, 1 - 1e-6),
        "isotonic_probability": np.clip(isotonic, 1e-6, 1 - 1e-6),
    })
    return output


def calibration_bins(predictions: pd.DataFrame, population: str, method: str, p_col: str, strategy: str) -> list[dict[str, object]]:
    y = predictions[Y_COL].to_numpy(dtype=int)
    p = predictions[p_col].to_numpy(dtype=float)
    assignments = bin_assignments(p, strategy)
    rows = []
    for bin_id in sorted(np.unique(assignments)):
        mask = assignments == bin_id
        n = int(mask.sum())
        successes = int(y[mask].sum())
        low, high = wilson_interval(successes, n)
        rows.append({
            "population": population,
            "method": method,
            "strategy": strategy,
            "bin_id": int(bin_id),
            "n": n,
            "probability_min": float(np.min(p[mask])),
            "probability_max": float(np.max(p[mask])),
            "mean_predicted_probability": float(np.mean(p[mask])),
            "observed_positive_fraction": float(np.mean(y[mask])),
            "observed_positive_wilson95_low": low,
            "observed_positive_wilson95_high": high,
            "absolute_calibration_gap": abs(float(np.mean(p[mask])) - float(np.mean(y[mask]))),
        })
    return rows


def brier_decomposition(y: np.ndarray, p: np.ndarray) -> dict[str, float]:
    assignments = bin_assignments(p, "adaptive")
    prevalence = float(np.mean(y))
    reliability = 0.0
    resolution = 0.0
    for bin_id in np.unique(assignments):
        mask = assignments == bin_id
        weight = mask.mean()
        mean_p = float(np.mean(p[mask]))
        mean_y = float(np.mean(y[mask]))
        reliability += weight * (mean_p - mean_y) ** 2
        resolution += weight * (mean_y - prevalence) ** 2
    uncertainty = prevalence * (1.0 - prevalence)
    brier = float(np.mean((p - y) ** 2))
    residual = brier - (reliability - resolution + uncertainty)
    return {
        "brier_score": brier,
        "reliability": reliability,
        "resolution": resolution,
        "uncertainty": uncertainty,
        "binning_residual": residual,
    }


def calibration_analysis(predictions: pd.DataFrame, population: str) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    metric_rows = []
    bin_rows = []
    decomposition_rows = []
    y = predictions[Y_COL].to_numpy(dtype=int)
    for method, p_col in METHOD_COLUMNS.items():
        p = predictions[p_col].to_numpy(dtype=float)
        metric_rows.append({
            "population": population,
            "method": method,
            "n": int(len(y)),
            "n_positive": int(y.sum()),
            "n_negative": int(len(y) - y.sum()),
            "observed_prevalence": float(y.mean()),
            "mean_predicted_probability": float(p.mean()),
            **calibration_metric_values(y, p),
        })
        for strategy in ["equal_width", "adaptive"]:
            bin_rows.extend(calibration_bins(predictions, population, method, p_col, strategy))
        decomposition_rows.append({"population": population, "method": method, **brier_decomposition(y, p)})
    return metric_rows, bin_rows, decomposition_rows


def calibration_bootstrap(predictions: pd.DataFrame, population: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    frame = predictions.reset_index(drop=True)
    groups = [np.asarray(list(indices), dtype=int) for indices in frame.groupby(CLUSTER_COL, sort=False).groups.values()]
    rng = np.random.default_rng(SEED + 20000 + 100 * (population == "scaffold_clean"))
    values = {
        (method, metric): []
        for method in METHOD_COLUMNS
        for metric in ["brier_score", "log_loss", "ece_equal_width_10", "ece_adaptive_10", "calibration_intercept", "calibration_slope", "roc_auc", "average_precision"]
    }
    pair_values = {
        (method, metric): []
        for method in ["cross_fitted_platt", "cross_fitted_isotonic"]
        for metric in ["brier_score", "log_loss", "ece_equal_width_10", "ece_adaptive_10", "calibration_intercept", "calibration_slope", "roc_auc", "average_precision"]
    }
    for _ in range(N_BOOT):
        selected_groups = rng.integers(0, len(groups), size=len(groups))
        sample_indices = np.concatenate([groups[index] for index in selected_groups])
        sample = frame.iloc[sample_indices]
        y = sample[Y_COL].to_numpy(dtype=int)
        if len(np.unique(y)) < 2:
            continue
        replicate_metrics: dict[str, dict[str, float]] = {}
        for method, p_col in METHOD_COLUMNS.items():
            metrics = calibration_metric_values(y, sample[p_col].to_numpy(dtype=float))
            replicate_metrics[method] = metrics
            for metric, value in metrics.items():
                if np.isfinite(value):
                    values[(method, metric)].append(value)
        raw_metrics = replicate_metrics["raw_frozen_score"]
        for method in ["cross_fitted_platt", "cross_fitted_isotonic"]:
            calibrated_metrics = replicate_metrics[method]
            for metric in ["brier_score", "log_loss", "ece_equal_width_10", "ece_adaptive_10"]:
                favorable = raw_metrics[metric] - calibrated_metrics[metric]
                if np.isfinite(favorable):
                    pair_values[(method, metric)].append(favorable)
            intercept_gain = abs(raw_metrics["calibration_intercept"]) - abs(calibrated_metrics["calibration_intercept"])
            slope_gain = abs(raw_metrics["calibration_slope"] - 1.0) - abs(calibrated_metrics["calibration_slope"] - 1.0)
            auc_gain = calibrated_metrics["roc_auc"] - raw_metrics["roc_auc"]
            ap_gain = calibrated_metrics["average_precision"] - raw_metrics["average_precision"]
            for metric, favorable in [
                ("calibration_intercept", intercept_gain),
                ("calibration_slope", slope_gain),
                ("roc_auc", auc_gain),
                ("average_precision", ap_gain),
            ]:
                if np.isfinite(favorable):
                    pair_values[(method, metric)].append(favorable)
    point_lookup = {
        method: calibration_metric_values(frame[Y_COL].to_numpy(dtype=int), frame[p_col].to_numpy(dtype=float))
        for method, p_col in METHOD_COLUMNS.items()
    }
    rows = []
    for (method, metric), metric_values in values.items():
        array = np.asarray(metric_values, dtype=float)
        rows.append({
            "population": population,
            "method": method,
            "metric": metric,
            "point_estimate": point_lookup[method][metric],
            "ci95_low": float(np.quantile(array, 0.025)),
            "ci95_high": float(np.quantile(array, 0.975)),
            "valid_bootstrap_replicates": int(len(array)),
            "n_bootstrap_clusters": int(len(groups)),
        })
    pair_rows = []
    for (method, metric), metric_values in pair_values.items():
        array = np.asarray(metric_values, dtype=float)
        raw_point = point_lookup["raw_frozen_score"][metric]
        calibrated_point = point_lookup[method][metric]
        if metric in ["brier_score", "log_loss", "ece_equal_width_10", "ece_adaptive_10"]:
            point = raw_point - calibrated_point
            contrast = "raw_minus_calibrated_lower_is_better"
        elif metric == "calibration_intercept":
            point = abs(raw_point) - abs(calibrated_point)
            contrast = "raw_minus_calibrated_absolute_intercept_error"
        elif metric == "calibration_slope":
            point = abs(raw_point - 1.0) - abs(calibrated_point - 1.0)
            contrast = "raw_minus_calibrated_absolute_slope_error"
        else:
            point = calibrated_point - raw_point
            contrast = "calibrated_minus_raw_higher_is_better"
        pair_rows.append({
            "population": population,
            "calibrated_method": method,
            "metric": metric,
            "contrast": contrast,
            "raw_value": raw_point,
            "calibrated_value": calibrated_point,
            "favorable_improvement": point,
            "ci95_low": float(np.quantile(array, 0.025)),
            "ci95_high": float(np.quantile(array, 0.975)),
            "bootstrap_probability_favorable_gt_0": float(np.mean(array > 0)),
            "valid_bootstrap_replicates": int(len(array)),
        })
    return rows, pair_rows


def prepare_error_frame(work: pd.DataFrame, protocol: dict[str, object]) -> pd.DataFrame:
    frame = work.copy().reset_index(drop=True)
    frame["raw_probability"] = frame[P_COL]
    frame["forced_prediction"] = (frame[P_COL] >= 0.5).astype(int)
    frame["is_error"] = (frame["forced_prediction"] != frame[Y_COL]).astype(int)
    frame["error_outcome"] = np.select(
        [
            (frame[Y_COL] == 1) & (frame["forced_prediction"] == 1),
            (frame[Y_COL] == 0) & (frame["forced_prediction"] == 0),
            (frame[Y_COL] == 0) & (frame["forced_prediction"] == 1),
            (frame[Y_COL] == 1) & (frame["forced_prediction"] == 0),
        ],
        ["true_positive", "true_negative", "false_positive", "false_negative"],
        default="unresolved",
    )
    frame["selection_layer"] = np.where(frame["risk_call"].isin([HIGH, LOW]), "v0_4_accepted", "v0_4_rejected")
    thresholds = protocol["error_taxonomy"]["support_thresholds"]
    primary = pd.to_numeric(frame["primary_max_tanimoto"], errors="coerce").ge(float(thresholds["primary_max_tanimoto"]))
    direct = pd.to_numeric(frame["direct_max_tanimoto"], errors="coerce").ge(float(thresholds["direct_max_tanimoto"]))
    frame["support_combination"] = np.select(
        [primary & direct, primary & ~direct, ~primary & direct, ~primary & ~direct],
        ["both_supported", "primary_only", "direct_only", "neither_supported"],
        default="unresolved",
    )
    edges = protocol["error_taxonomy"]["probability_band_edges"]
    labels = [f"[{edges[index]:.2f},{edges[index + 1]:.2f}{']' if index == len(edges) - 2 else ')'}" for index in range(len(edges) - 1)]
    frame["probability_band"] = pd.cut(frame[P_COL], bins=edges, labels=labels, include_lowest=True, right=False)
    frame.loc[frame[P_COL].eq(1.0), "probability_band"] = labels[-1]

    descriptor_names = ["molecular_weight", "xlogp", "tpsa", "hbond_donors", "hbond_acceptors", "rotatable_bonds", "formal_charge", "heavy_atoms"]
    descriptor_rows = []
    for value in frame["physicochemical_descriptors_json"].fillna("{}"):
        try:
            parsed = json.loads(value)
        except (json.JSONDecodeError, TypeError):
            parsed = {}
        descriptor_rows.append({name: parsed.get(name, np.nan) for name in descriptor_names})
    descriptors = pd.DataFrame(descriptor_rows, index=frame.index)
    for name in descriptor_names:
        frame[name] = pd.to_numeric(descriptors[name], errors="coerce")
    return frame


def error_count_rows(frame: pd.DataFrame, population: str) -> list[dict[str, object]]:
    rows = []
    for outcome, group in frame.groupby("error_outcome", sort=False):
        rows.append({
            "population": population,
            "error_outcome": outcome,
            "n": int(len(group)),
            "fraction_of_population": float(len(group) / len(frame)),
            "median_raw_probability": float(group[P_COL].median()),
            "median_log10_efflux_ratio": float(group["log10_mdr1_mdck_er"].median()),
            "v0_4_accepted_n": int(group["risk_call"].isin([HIGH, LOW]).sum()),
            "v0_4_accepted_fraction": float(group["risk_call"].isin([HIGH, LOW]).mean()),
        })
    return rows


def error_strata_rows(frame: pd.DataFrame, population: str) -> list[dict[str, object]]:
    rows = []
    columns = {
        "selection_layer": "selection_layer",
        "risk_call": "risk_call",
        "applicability_status": "applicability_status",
        "uncertainty_status": "uncertainty_status",
        "support_combination": "support_combination",
        "probability_band": "probability_band",
    }
    for stratum_type, column in columns.items():
        for level, group in frame.groupby(column, dropna=False, observed=True):
            metrics = metrics_for_frame(group)
            errors = int(group["is_error"].sum())
            low, high = wilson_interval(errors, len(group))
            rows.append({
                "population": population,
                "stratum_type": stratum_type,
                "stratum_level": str(level),
                "n": int(len(group)),
                "fraction_of_population": float(len(group) / len(frame)),
                "error_n": errors,
                "error_rate": float(errors / len(group)),
                "error_rate_wilson95_low": low,
                "error_rate_wilson95_high": high,
                "false_positive_n": int((group["error_outcome"] == "false_positive").sum()),
                "false_negative_n": int((group["error_outcome"] == "false_negative").sum()),
                "balanced_accuracy": metrics["balanced_accuracy"],
                "mcc": metrics["mcc"],
                "sensitivity": metrics["sensitivity"],
                "specificity": metrics["specificity"],
            })
    return rows


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    result = mannwhitneyu(a, b, alternative="two-sided")
    delta = 2.0 * float(result.statistic) / (len(a) * len(b)) - 1.0
    return delta, float(result.pvalue)


def standardized_mean_difference(a: np.ndarray, b: np.ndarray) -> float:
    var_a = np.var(a, ddof=1) if len(a) > 1 else 0.0
    var_b = np.var(b, ddof=1) if len(b) > 1 else 0.0
    pooled = math.sqrt(((len(a) - 1) * var_a + (len(b) - 1) * var_b) / max(len(a) + len(b) - 2, 1))
    return float((np.mean(a) - np.mean(b)) / pooled) if pooled > 0 else np.nan


def bh_adjust(pvalues: pd.Series) -> pd.Series:
    values = pvalues.to_numpy(dtype=float)
    finite = np.isfinite(values)
    output = np.full(len(values), np.nan, dtype=float)
    finite_values = values[finite]
    if len(finite_values) == 0:
        return pd.Series(output, index=pvalues.index)
    order = np.argsort(finite_values)
    ranked = finite_values[order]
    adjusted = np.minimum.accumulate((ranked * len(finite_values) / np.arange(1, len(finite_values) + 1))[::-1])[::-1]
    finite_output = np.empty(len(finite_values), dtype=float)
    finite_output[order] = np.clip(adjusted, 0, 1)
    output[np.where(finite)[0]] = finite_output
    return pd.Series(output, index=pvalues.index)


def error_feature_rows(frame: pd.DataFrame, population: str) -> list[dict[str, object]]:
    feature_map = {
        "raw_probability": P_COL,
        "log10_efflux_ratio": "log10_mdr1_mdck_er",
        "bootstrap_p10": "bootstrap_probability_p10",
        "bootstrap_p90": "bootstrap_probability_p90",
        "bootstrap_interval_width": "bootstrap_interval_width",
        "direct_probability": "direct_panel_probability",
        "primary_max_tanimoto": "primary_max_tanimoto",
        "direct_max_tanimoto": "direct_max_tanimoto",
        "primary_local_positive_fraction": "primary_local_positive_fraction",
        "direct_local_positive_fraction": "direct_local_positive_fraction",
        "primary_local_conflict_index": "primary_local_conflict_index",
        "direct_local_conflict_index": "direct_local_conflict_index",
        "physchem_max_abs_robust_z": "physchem_max_abs_robust_z",
        "molecular_weight": "molecular_weight",
        "xlogp": "xlogp",
        "tpsa": "tpsa",
        "hbond_donors": "hbond_donors",
        "hbond_acceptors": "hbond_acceptors",
        "rotatable_bonds": "rotatable_bonds",
        "formal_charge": "formal_charge",
        "heavy_atoms": "heavy_atoms",
    }
    comparisons = {
        "all_error_vs_correct": (frame["is_error"].eq(1), frame["is_error"].eq(0), "error", "correct"),
        "false_positive_vs_true_positive": (frame["error_outcome"].eq("false_positive"), frame["error_outcome"].eq("true_positive"), "false_positive", "true_positive"),
        "false_negative_vs_true_negative": (frame["error_outcome"].eq("false_negative"), frame["error_outcome"].eq("true_negative"), "false_negative", "true_negative"),
    }
    rows = []
    for comparison, (mask_a, mask_b, label_a, label_b) in comparisons.items():
        for feature, column in feature_map.items():
            a = pd.to_numeric(frame.loc[mask_a, column], errors="coerce").dropna().to_numpy(dtype=float)
            b = pd.to_numeric(frame.loc[mask_b, column], errors="coerce").dropna().to_numpy(dtype=float)
            if len(a) < 2 or len(b) < 2:
                continue
            delta, pvalue = cliffs_delta(a, b)
            rows.append({
                "population": population,
                "comparison": comparison,
                "group_a": label_a,
                "group_b": label_b,
                "feature": feature,
                "group_a_n": int(len(a)),
                "group_b_n": int(len(b)),
                "group_a_mean": float(np.mean(a)),
                "group_b_mean": float(np.mean(b)),
                "group_a_median": float(np.median(a)),
                "group_b_median": float(np.median(b)),
                "median_difference_a_minus_b": float(np.median(a) - np.median(b)),
                "standardized_mean_difference": standardized_mean_difference(a, b),
                "cliffs_delta": delta,
                "mannwhitney_p": pvalue,
            })
    output = pd.DataFrame(rows)
    if not output.empty:
        output["fdr_bh_within_comparison"] = output.groupby(["population", "comparison"])["mannwhitney_p"].transform(bh_adjust)
    return output.to_dict(orient="records")


def scaffold_hotspot_rows(frame: pd.DataFrame, population: str, minimum_n: int) -> list[dict[str, object]]:
    overall_error = float(frame["is_error"].mean())
    prior_strength = 10.0
    alpha = overall_error * prior_strength
    beta = (1.0 - overall_error) * prior_strength
    rows = []
    for cluster, group in frame.groupby(CLUSTER_COL, sort=False):
        if len(group) < minimum_n:
            continue
        errors = int(group["is_error"].sum())
        low, high = wilson_interval(errors, len(group))
        rows.append({
            "population": population,
            "bootstrap_cluster": cluster,
            "n": int(len(group)),
            "error_n": errors,
            "raw_error_rate": float(errors / len(group)),
            "error_rate_wilson95_low": low,
            "error_rate_wilson95_high": high,
            "shrunken_error_rate": float((errors + alpha) / (len(group) + alpha + beta)),
            "false_positive_n": int((group["error_outcome"] == "false_positive").sum()),
            "false_negative_n": int((group["error_outcome"] == "false_negative").sum()),
            "true_positive_n": int((group["error_outcome"] == "true_positive").sum()),
            "true_negative_n": int((group["error_outcome"] == "true_negative").sum()),
            "median_raw_probability": float(group[P_COL].median()),
            "median_log10_efflux_ratio": float(group["log10_mdr1_mdck_er"].median()),
            "dominant_support_combination": str(group["support_combination"].mode().iloc[0]),
            "v0_4_accepted_n": int(group["risk_call"].isin([HIGH, LOW]).sum()),
        })
    return rows


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text())
    assert protocol["status"].startswith("locked_posthoc_deep_diagnostic")
    for entry in protocol["frozen_inputs"].values():
        path = PROJECT / entry["path"]
        assert path.exists(), path
        assert sha256(path) == entry["sha256"], f"Hash mismatch: {path}"
    assert protocol["calibration_diagnostics"]["cluster_bootstrap"]["replicates"] == N_BOOT
    assert protocol["calibration_diagnostics"]["cluster_bootstrap"]["seed"] == SEED

    raw = pd.read_csv(INPUT, low_memory=False)
    prediction_frames = []
    metric_rows = []
    bootstrap_rows = []
    calibration_pairwise_rows = []
    bin_rows = []
    decomposition_rows = []
    error_count_output = []
    error_strata_output = []
    error_feature_output = []
    hotspot_output = []
    error_frames: dict[str, pd.DataFrame] = {}
    for population in protocol["analysis_populations"]:
        work = clean_frame(raw, population)
        predictions = cross_fitted_predictions(work, population, protocol["calibration_diagnostics"]["cross_fitting"]["folds"])
        prediction_frames.append(predictions)
        metrics, bins, decompositions = calibration_analysis(predictions, population)
        metric_rows.extend(metrics)
        bin_rows.extend(bins)
        decomposition_rows.extend(decompositions)
        population_bootstrap, population_pairwise = calibration_bootstrap(predictions, population)
        bootstrap_rows.extend(population_bootstrap)
        calibration_pairwise_rows.extend(population_pairwise)

        error_frame = prepare_error_frame(work, protocol)
        error_frames[population] = error_frame
        error_count_output.extend(error_count_rows(error_frame, population))
        error_strata_output.extend(error_strata_rows(error_frame, population))
        error_feature_output.extend(error_feature_rows(error_frame, population))
        hotspot_output.extend(scaffold_hotspot_rows(
            error_frame,
            population,
            int(protocol["error_taxonomy"]["scaffold_hotspots"]["minimum_cluster_n"]),
        ))

    outputs = {
        PREDICTIONS: round_frame(pd.concat(prediction_frames, ignore_index=True)),
        CALIBRATION_METRICS: round_frame(pd.DataFrame(metric_rows)),
        CALIBRATION_BOOTSTRAP: round_frame(pd.DataFrame(bootstrap_rows)),
        CALIBRATION_PAIRWISE: round_frame(pd.DataFrame(calibration_pairwise_rows)),
        CALIBRATION_BINS: round_frame(pd.DataFrame(bin_rows)),
        BRIER_DECOMPOSITION: round_frame(pd.DataFrame(decomposition_rows)),
        ERROR_COUNTS: round_frame(pd.DataFrame(error_count_output)),
        ERROR_FEATURES: round_frame(pd.DataFrame(error_feature_output)),
        ERROR_STRATA: round_frame(pd.DataFrame(error_strata_output)),
        SCAFFOLD_HOTSPOTS: round_frame(pd.DataFrame(hotspot_output).sort_values(["population", "shrunken_error_rate"], ascending=[True, False])),
    }
    for path, frame in outputs.items():
        frame.to_csv(path, index=False)

    metrics_table = outputs[CALIBRATION_METRICS]
    identity_metrics = metrics_table[metrics_table["population"].eq("identity_clean")].set_index("method")
    raw_metrics = identity_metrics.loc["raw_frozen_score"]
    platt_metrics = identity_metrics.loc["cross_fitted_platt"]
    isotonic_metrics = identity_metrics.loc["cross_fitted_isotonic"]
    counts = outputs[ERROR_COUNTS]
    identity_counts = counts[counts["population"].eq("identity_clean")].set_index("error_outcome")
    strata = outputs[ERROR_STRATA]
    identity_strata = strata[(strata["population"] == "identity_clean") & (strata["n"] >= 50)].sort_values("error_rate", ascending=False)
    features = outputs[ERROR_FEATURES]
    identity_features = features[(features["population"] == "identity_clean") & (features["feature"] != "log10_efflux_ratio") & (features["fdr_bh_within_comparison"] <= 0.05)].copy()
    identity_features["absolute_cliffs_delta"] = identity_features["cliffs_delta"].abs()
    top_features = identity_features.sort_values("absolute_cliffs_delta", ascending=False).head(12)
    hotspots = outputs[SCAFFOLD_HOTSPOTS]
    identity_hotspots = hotspots[hotspots["population"].eq("identity_clean")].head(10)

    bootstrap_table = outputs[CALIBRATION_BOOTSTRAP]
    raw_slope_boot = bootstrap_table[(bootstrap_table["population"] == "identity_clean") & (bootstrap_table["method"] == "raw_frozen_score") & (bootstrap_table["metric"] == "calibration_slope")].iloc[0]
    raw_intercept_boot = bootstrap_table[(bootstrap_table["population"] == "identity_clean") & (bootstrap_table["method"] == "raw_frozen_score") & (bootstrap_table["metric"] == "calibration_intercept")].iloc[0]
    raw_ece_boot = bootstrap_table[(bootstrap_table["population"] == "identity_clean") & (bootstrap_table["method"] == "raw_frozen_score") & (bootstrap_table["metric"] == "ece_equal_width_10")].iloc[0]
    calibration_pairwise_table = outputs[CALIBRATION_PAIRWISE]
    platt_pairwise = calibration_pairwise_table[(calibration_pairwise_table["population"] == "identity_clean") & (calibration_pairwise_table["calibrated_method"] == "cross_fitted_platt")].set_index("metric")

    summary = {
        "analysis_version": "ABCB1_EXTERNAL_CALIBRATION_ERROR_V2",
        "analysis_date": "2026-08-05",
        "status": "completed_frozen_model_posthoc_deep_calibration_and_error_diagnostic",
        "protocol_path": str(PROTOCOL.relative_to(PROJECT)),
        "protocol_sha256": sha256(PROTOCOL),
        "identity_clean_raw_calibration": {
            "n": int(raw_metrics["n"]),
            "observed_prevalence": raw_metrics["observed_prevalence"],
            "mean_predicted_probability": raw_metrics["mean_predicted_probability"],
            "brier_score": raw_metrics["brier_score"],
            "log_loss": raw_metrics["log_loss"],
            "ece_equal_width_10": raw_metrics["ece_equal_width_10"],
            "ece_ci95": [raw_ece_boot["ci95_low"], raw_ece_boot["ci95_high"]],
            "calibration_intercept": raw_metrics["calibration_intercept"],
            "calibration_intercept_ci95": [raw_intercept_boot["ci95_low"], raw_intercept_boot["ci95_high"]],
            "calibration_slope": raw_metrics["calibration_slope"],
            "calibration_slope_ci95": [raw_slope_boot["ci95_low"], raw_slope_boot["ci95_high"]],
            "roc_auc": raw_metrics["roc_auc"],
            "average_precision": raw_metrics["average_precision"],
        },
        "within_source_crossfitted_recalibration": {
            "platt": {
                "brier_score": platt_metrics["brier_score"],
                "log_loss": platt_metrics["log_loss"],
                "ece_equal_width_10": platt_metrics["ece_equal_width_10"],
                "calibration_intercept": platt_metrics["calibration_intercept"],
                "calibration_slope": platt_metrics["calibration_slope"],
                "brier_favorable_improvement_ci95": [platt_pairwise.loc["brier_score", "ci95_low"], platt_pairwise.loc["brier_score", "ci95_high"]],
                "ece_favorable_improvement_ci95": [platt_pairwise.loc["ece_equal_width_10", "ci95_low"], platt_pairwise.loc["ece_equal_width_10", "ci95_high"]],
            },
            "isotonic": {
                "brier_score": isotonic_metrics["brier_score"],
                "log_loss": isotonic_metrics["log_loss"],
                "ece_equal_width_10": isotonic_metrics["ece_equal_width_10"],
                "calibration_intercept": isotonic_metrics["calibration_intercept"],
                "calibration_slope": isotonic_metrics["calibration_slope"],
            },
            "warning": "Same-source group-cross-fitted diagnostic only; not an independent calibrated model.",
        },
        "identity_clean_error_counts": {
            outcome: int(identity_counts.loc[outcome, "n"])
            for outcome in ["true_positive", "true_negative", "false_positive", "false_negative"]
        },
        "highest_error_strata_n_ge_50": identity_strata.head(10)[["stratum_type", "stratum_level", "n", "error_rate", "false_positive_n", "false_negative_n"]].to_dict(orient="records"),
        "largest_fdr_significant_feature_differences": top_features[["comparison", "feature", "group_a_median", "group_b_median", "cliffs_delta", "fdr_bh_within_comparison"]].to_dict(orient="records"),
        "top_shrunken_error_scaffold_hotspots": identity_hotspots[["bootstrap_cluster", "n", "raw_error_rate", "shrunken_error_rate", "false_positive_n", "false_negative_n"]].to_dict(orient="records"),
        "interpretation": protocol["interpretation_rules"],
        "output_hashes": {str(path.relative_to(PROJECT)): sha256(path) for path in outputs},
    }
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=lambda value: value.item() if isinstance(value, np.generic) else value) + "\n")
    print(json.dumps({
        "identity_clean_raw_calibration": summary["identity_clean_raw_calibration"],
        "within_source_crossfitted_recalibration": summary["within_source_crossfitted_recalibration"],
        "identity_clean_error_counts": summary["identity_clean_error_counts"],
        "highest_error_strata_n_ge_50": summary["highest_error_strata_n_ge_50"][:6],
        "outputs": [str(path.relative_to(PROJECT)) for path in [*outputs, SUMMARY]],
    }, ensure_ascii=False, indent=2, default=lambda value: value.item() if isinstance(value, np.generic) else value))


if __name__ == "__main__":
    main()
