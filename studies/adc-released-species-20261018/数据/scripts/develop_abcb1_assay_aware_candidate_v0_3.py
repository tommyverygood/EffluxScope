#!/usr/bin/env python3
"""Internally develop an assay-aware ABCB1 screening candidate.

All evaluation in this script is single-source development. The existing
external benchmark is intentionally not read, scored for selection, or used to
change thresholds.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import (
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    mean_absolute_error,
    mean_squared_error,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from abcb1_screener_core_v0_1 import (  # noqa: E402
    FP_GENERATOR,
    descriptor_array,
    load_package,
    morgan_array,
    predict_smiles,
)


PROTOCOL_PATH = PROJECT / "数据/processed/abcb1_assay_aware_development_protocol_lock_v1.json"
PANEL_PATH = PROJECT / "数据/processed/soskuti_2024_mdr1_species_panel_v1.csv"
FROZEN_MODEL_PATH = PROJECT / "模型/ABCB1_screener_v0_1/ABCB1_screener_v0_1.joblib"
OUTPUT_DIR = PROJECT / "数据/processed"
MODEL_DIR = PROJECT / "模型/ABCB1_assay_aware_candidate_v0_3_dev"
MODEL_PATH = MODEL_DIR / "ABCB1_assay_aware_candidate_v0_3_dev.joblib"

CLASSIFIERS = [
    "dummy_prior",
    "one_nearest_neighbor_morgan_tanimoto",
    "physchem_logistic",
    "morgan_logistic",
    "morgan_plus_physchem_logistic",
]
LEARNED_CLASSIFIERS = CLASSIFIERS[2:]
REGRESSORS = ["physchem_ridge", "morgan_ridge", "morgan_plus_physchem_ridge"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def safe_auc(y: np.ndarray, p: np.ndarray) -> float:
    return float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else float("nan")


def binary_metrics(y: np.ndarray, p: np.ndarray, threshold: float = 0.5) -> dict:
    pred = (p >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    return {
        "MCC": float(matthews_corrcoef(y, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "ROC_AUC": safe_auc(y, p),
        "Brier": float(brier_score_loss(y, p)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else float("nan"),
        "specificity": float(tn / (tn + fp)) if tn + fp else float("nan"),
        "F1": float(f1_score(y, pred, zero_division=0)),
        "TN": int(tn),
        "FP": int(fp),
        "FN": int(fn),
        "TP": int(tp),
    }


def continuous_metrics(y: np.ndarray, pred: np.ndarray) -> dict:
    rho = float(spearmanr(y, pred).statistic) if len(y) >= 3 else float("nan")
    return {
        "Spearman_rho": rho,
        "MAE_log2_net_ER": float(mean_absolute_error(y, pred)),
        "RMSE_log2_net_ER": float(math.sqrt(mean_squared_error(y, pred))),
    }


def scaffold_for_smiles(smiles: str, connectivity_block: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"cannot parse standardized structure: {smiles}")
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    return scaffold or f"ACYCLIC::{connectivity_block}"


def make_dataset(rows: list[dict], label_field: str) -> dict:
    selected = [row for row in rows if row.get(label_field, "") != ""]
    ids, names, smiles, inchikeys, groups = [], [], [], [], []
    labels, continuous, fps, desc, mols = [], [], [], [], []
    for row in selected:
        mol = Chem.MolFromSmiles(row["standardized_isomeric_smiles"])
        if mol is None:
            raise ValueError(row["panel_entity_id"])
        block = row["model_standardized_connectivity_block"]
        ids.append(row["panel_entity_id"])
        names.append(row["entity_name_source"])
        smiles.append(row["standardized_isomeric_smiles"])
        inchikeys.append(row["standardized_inchikey"])
        groups.append(scaffold_for_smiles(row["standardized_isomeric_smiles"], block))
        labels.append(int(row[label_field]))
        continuous.append(math.log2(float(row["human_mdr1_net_efflux_ratio_vs_mock"])))
        fps.append(morgan_array(mol))
        desc.append(descriptor_array(mol))
        mols.append(mol)
    return {
        "ids": np.asarray(ids, dtype=object),
        "names": np.asarray(names, dtype=object),
        "smiles": np.asarray(smiles, dtype=object),
        "inchikeys": np.asarray(inchikeys, dtype=object),
        "groups": np.asarray(groups, dtype=object),
        "y": np.asarray(labels, dtype=int),
        "continuous_y": np.asarray(continuous, dtype=float),
        "fps": np.asarray(fps, dtype=np.float32),
        "desc": np.asarray(desc, dtype=np.float64),
        "mols": mols,
    }


def make_features(kind: str, fps: np.ndarray, desc: np.ndarray, train_idx: np.ndarray, test_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray, StandardScaler | None]:
    scaler: StandardScaler | None = None
    if kind.startswith("physchem"):
        scaler = StandardScaler().fit(desc[train_idx])
        return scaler.transform(desc[train_idx]), scaler.transform(desc[test_idx]), scaler
    if kind.startswith("morgan_plus_physchem"):
        scaler = StandardScaler().fit(desc[train_idx])
        train = np.hstack([fps[train_idx], scaler.transform(desc[train_idx])])
        test = np.hstack([fps[test_idx], scaler.transform(desc[test_idx])])
        return train, test, scaler
    if kind.startswith("morgan"):
        return fps[train_idx], fps[test_idx], None
    raise ValueError(kind)


def inner_splits(y: np.ndarray, groups: np.ndarray, seed: int) -> list[tuple[np.ndarray, np.ndarray]]:
    for n_splits in (4, 3, 2):
        try:
            splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
            splits = list(splitter.split(np.zeros(len(y)), y, groups))
            if all(len(np.unique(y[train])) == 2 and len(np.unique(y[test])) == 2 for train, test in splits):
                return splits
        except ValueError:
            continue
    return []


def select_logistic_c(kind: str, dataset: dict, outer_train: np.ndarray, seed: int, c_grid: list[float]) -> float:
    y = dataset["y"][outer_train]
    groups = dataset["groups"][outer_train]
    splits = inner_splits(y, groups, seed)
    if not splits:
        return 0.1
    candidates = []
    for c in c_grid:
        mccs, baccs = [], []
        for inner_train_rel, inner_test_rel in splits:
            inner_train = outer_train[inner_train_rel]
            inner_test = outer_train[inner_test_rel]
            x_train, x_test, _ = make_features(kind, dataset["fps"], dataset["desc"], inner_train, inner_test)
            model = LogisticRegression(C=c, class_weight="balanced", solver="liblinear", max_iter=5000, random_state=seed)
            model.fit(x_train, dataset["y"][inner_train])
            pred = model.predict(x_test)
            mccs.append(matthews_corrcoef(dataset["y"][inner_test], pred))
            baccs.append(balanced_accuracy_score(dataset["y"][inner_test], pred))
        candidates.append((float(np.mean(mccs)), float(np.mean(baccs)), -c, c))
    return float(max(candidates)[-1])


def select_ridge_alpha(kind: str, dataset: dict, outer_train: np.ndarray, seed: int, alpha_grid: list[float]) -> float:
    y_class = dataset["y"][outer_train]
    groups = dataset["groups"][outer_train]
    splits = inner_splits(y_class, groups, seed)
    if not splits:
        return 10.0
    candidates = []
    for alpha in alpha_grid:
        maes, rhos = [], []
        for inner_train_rel, inner_test_rel in splits:
            inner_train = outer_train[inner_train_rel]
            inner_test = outer_train[inner_test_rel]
            x_train, x_test, _ = make_features(kind, dataset["fps"], dataset["desc"], inner_train, inner_test)
            model = Ridge(alpha=alpha)
            model.fit(x_train, dataset["continuous_y"][inner_train])
            pred = model.predict(x_test)
            maes.append(mean_absolute_error(dataset["continuous_y"][inner_test], pred))
            rho = spearmanr(dataset["continuous_y"][inner_test], pred).statistic
            rhos.append(float(rho) if np.isfinite(rho) else -1.0)
        # lower MAE, higher rho, and stronger regularization win ties
        candidates.append((-float(np.mean(maes)), float(np.mean(rhos)), alpha, alpha))
    return float(max(candidates)[-1])


def nearest_neighbor_probabilities(dataset: dict, train_idx: np.ndarray, test_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    train_fps = [FP_GENERATOR.GetFingerprint(dataset["mols"][i]) for i in train_idx]
    probs, sims = [], []
    for test_i in test_idx:
        fp = FP_GENERATOR.GetFingerprint(dataset["mols"][test_i])
        similarities = np.asarray(DataStructs.BulkTanimotoSimilarity(fp, train_fps), dtype=float)
        maximum = float(np.max(similarities))
        tied = np.isclose(similarities, maximum, rtol=0, atol=1e-12)
        probs.append(float(np.mean(dataset["y"][train_idx][tied])))
        sims.append(maximum)
    return np.asarray(probs), np.asarray(sims)


def repeated_classification(dataset: dict, protocol: dict, analysis_set: str) -> tuple[list[dict], list[dict], dict[str, np.ndarray]]:
    y = dataset["y"]
    prediction_lists: dict[str, list[list[float]]] = {model: [[] for _ in range(len(y))] for model in CLASSIFIERS}
    fold_rows: list[dict] = []
    detail_rows: list[dict] = []
    seeds = protocol["outer_evaluation"]["random_seeds"]
    c_grid = [float(value) for value in protocol["inner_selection"]["logistic_C_grid"]]
    for repeat_index, seed in enumerate(seeds, start=1):
        splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
        for fold_index, (train_idx, test_idx) in enumerate(splitter.split(np.zeros(len(y)), y, dataset["groups"]), start=1):
            fold_predictions: dict[str, np.ndarray] = {}
            selected_parameters: dict[str, str] = {}
            dummy = DummyClassifier(strategy="prior").fit(np.zeros((len(train_idx), 1)), y[train_idx])
            fold_predictions["dummy_prior"] = dummy.predict_proba(np.zeros((len(test_idx), 1)))[:, list(dummy.classes_).index(1)]
            nn_prob, nn_sim = nearest_neighbor_probabilities(dataset, train_idx, test_idx)
            fold_predictions["one_nearest_neighbor_morgan_tanimoto"] = nn_prob
            selected_parameters["one_nearest_neighbor_morgan_tanimoto"] = "tied_max_similarity_labels_averaged"
            for kind in LEARNED_CLASSIFIERS:
                c = select_logistic_c(kind, dataset, train_idx, seed * 100 + fold_index, c_grid)
                x_train, x_test, _ = make_features(kind, dataset["fps"], dataset["desc"], train_idx, test_idx)
                model = LogisticRegression(C=c, class_weight="balanced", solver="liblinear", max_iter=5000, random_state=seed)
                model.fit(x_train, y[train_idx])
                fold_predictions[kind] = model.predict_proba(x_test)[:, 1]
                selected_parameters[kind] = f"C={c:g}"
            for model_name, probabilities in fold_predictions.items():
                metrics = binary_metrics(y[test_idx], probabilities)
                fold_rows.append(
                    {
                        "analysis_set": analysis_set,
                        "repeat": repeat_index,
                        "seed": seed,
                        "fold": fold_index,
                        "model": model_name,
                        "n_train": len(train_idx),
                        "n_test": len(test_idx),
                        "n_test_positive": int(np.sum(y[test_idx] == 1)),
                        "n_test_negative": int(np.sum(y[test_idx] == 0)),
                        "selected_parameter": selected_parameters.get(model_name, "fold_prior"),
                        **metrics,
                    }
                )
                for local_position, entity_idx in enumerate(test_idx):
                    probability = float(probabilities[local_position])
                    prediction_lists[model_name][entity_idx].append(probability)
                    detail_rows.append(
                        {
                            "analysis_set": analysis_set,
                            "panel_entity_id": dataset["ids"][entity_idx],
                            "entity_name": dataset["names"][entity_idx],
                            "true_label": int(y[entity_idx]),
                            "scaffold_group": dataset["groups"][entity_idx],
                            "repeat": repeat_index,
                            "fold": fold_index,
                            "model": model_name,
                            "probability": round(probability, 8),
                            "binary_prediction": int(probability >= 0.5),
                            "nearest_training_tanimoto": round(float(nn_sim[local_position]), 6) if model_name == "one_nearest_neighbor_morgan_tanimoto" else "",
                        }
                    )
    aggregated = {model: np.asarray([np.mean(values) for values in per_entity], dtype=float) for model, per_entity in prediction_lists.items()}
    return fold_rows, detail_rows, aggregated


def repeated_regression(dataset: dict, protocol: dict) -> tuple[list[dict], list[dict], dict[str, np.ndarray]]:
    y_class = dataset["y"]
    y = dataset["continuous_y"]
    prediction_lists: dict[str, list[list[float]]] = {model: [[] for _ in range(len(y))] for model in REGRESSORS}
    fold_rows, detail_rows = [], []
    alpha_grid = [float(value) for value in protocol["inner_selection"]["ridge_alpha_grid"]]
    for repeat_index, seed in enumerate(protocol["outer_evaluation"]["random_seeds"], start=1):
        splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
        for fold_index, (train_idx, test_idx) in enumerate(splitter.split(np.zeros(len(y)), y_class, dataset["groups"]), start=1):
            for kind in REGRESSORS:
                alpha = select_ridge_alpha(kind, dataset, train_idx, seed * 100 + fold_index, alpha_grid)
                x_train, x_test, _ = make_features(kind, dataset["fps"], dataset["desc"], train_idx, test_idx)
                model = Ridge(alpha=alpha).fit(x_train, y[train_idx])
                pred = model.predict(x_test)
                fold_rows.append(
                    {
                        "analysis_set": "point_qc_eligible_continuous",
                        "repeat": repeat_index,
                        "seed": seed,
                        "fold": fold_index,
                        "model": kind,
                        "n_train": len(train_idx),
                        "n_test": len(test_idx),
                        "selected_parameter": f"alpha={alpha:g}",
                        **continuous_metrics(y[test_idx], pred),
                    }
                )
                for local_position, entity_idx in enumerate(test_idx):
                    value = float(pred[local_position])
                    prediction_lists[kind][entity_idx].append(value)
                    detail_rows.append(
                        {
                            "panel_entity_id": dataset["ids"][entity_idx],
                            "entity_name": dataset["names"][entity_idx],
                            "true_log2_net_ER": round(float(y[entity_idx]), 8),
                            "true_net_ER": round(float(2 ** y[entity_idx]), 8),
                            "scaffold_group": dataset["groups"][entity_idx],
                            "repeat": repeat_index,
                            "fold": fold_index,
                            "model": kind,
                            "predicted_log2_net_ER": round(value, 8),
                            "predicted_net_ER": round(float(2**value), 8),
                        }
                    )
    aggregated = {model: np.asarray([np.mean(values) for values in per_entity], dtype=float) for model, per_entity in prediction_lists.items()}
    return fold_rows, detail_rows, aggregated


def cluster_bootstrap_binary(y: np.ndarray, predictions: dict[str, np.ndarray], groups: np.ndarray, n_boot: int = 2000) -> tuple[list[dict], list[dict]]:
    unique_groups = np.unique(groups)
    group_indices = {group: np.flatnonzero(groups == group) for group in unique_groups}
    rng = np.random.default_rng(20260803)
    distributions: dict[tuple[str, str], list[float]] = defaultdict(list)
    pair_distributions: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    metrics_to_bootstrap = ["MCC", "balanced_accuracy", "ROC_AUC", "Brier"]
    successful = 0
    for _ in range(n_boot * 3):
        sampled_groups = rng.choice(unique_groups, size=len(unique_groups), replace=True)
        indices = np.concatenate([group_indices[group] for group in sampled_groups])
        if len(np.unique(y[indices])) < 2:
            continue
        boot_metrics = {model: binary_metrics(y[indices], pred[indices]) for model, pred in predictions.items()}
        for model, metric_values in boot_metrics.items():
            for metric in metrics_to_bootstrap:
                distributions[(model, metric)].append(metric_values[metric])
        for learned in LEARNED_CLASSIFIERS:
            for comparator in ("dummy_prior", "one_nearest_neighbor_morgan_tanimoto"):
                for metric in metrics_to_bootstrap:
                    pair_distributions[(learned, comparator, metric)].append(boot_metrics[learned][metric] - boot_metrics[comparator][metric])
        successful += 1
        if successful == n_boot:
            break
    summary_rows = []
    for model, pred in predictions.items():
        point = binary_metrics(y, pred)
        for metric, value in point.items():
            if metric in metrics_to_bootstrap:
                values = np.asarray(distributions[(model, metric)], dtype=float)
                summary_rows.append(
                    {
                        "model": model,
                        "metric": metric,
                        "point_estimate": value,
                        "bootstrap_median": float(np.nanmedian(values)),
                        "ci_low": float(np.nanpercentile(values, 2.5)),
                        "ci_high": float(np.nanpercentile(values, 97.5)),
                        "bootstrap_unit": "Murcko_scaffold_group",
                        "bootstrap_replicates": len(values),
                    }
                )
            else:
                summary_rows.append(
                    {
                        "model": model,
                        "metric": metric,
                        "point_estimate": value,
                        "bootstrap_median": "",
                        "ci_low": "",
                        "ci_high": "",
                        "bootstrap_unit": "not_bootstrapped_secondary_or_confusion_count",
                        "bootstrap_replicates": 0,
                    }
                )
    pair_rows = []
    for (learned, comparator, metric), values_list in pair_distributions.items():
        values = np.asarray(values_list, dtype=float)
        point = binary_metrics(y, predictions[learned])[metric] - binary_metrics(y, predictions[comparator])[metric]
        pair_rows.append(
            {
                "candidate_model": learned,
                "comparator_model": comparator,
                "metric": metric,
                "candidate_minus_comparator_point": point,
                "ci_low": float(np.nanpercentile(values, 2.5)),
                "ci_high": float(np.nanpercentile(values, 97.5)),
                "bootstrap_probability_difference_gt_0": float(np.mean(values > 0)),
                "interpretation_note": "For Brier, a negative difference favors the candidate; for other metrics, a positive difference favors the candidate.",
            }
        )
    return summary_rows, pair_rows


def select_final_classifier(metric_rows: list[dict]) -> str:
    lookup = {(row["model"], row["metric"]): row for row in metric_rows}
    ranked = []
    for model in LEARNED_CLASSIFIERS:
        ranked.append(
            (
                float(lookup[(model, "MCC")]["bootstrap_median"]),
                float(lookup[(model, "balanced_accuracy")]["bootstrap_median"]),
                -float(lookup[(model, "Brier")]["bootstrap_median"]),
                -LEARNED_CLASSIFIERS.index(model),
                model,
            )
        )
    return max(ranked)[-1]


def fit_full_classifier(kind: str, dataset: dict, protocol: dict) -> dict:
    all_idx = np.arange(len(dataset["y"]))
    c = select_logistic_c(kind, dataset, all_idx, 20260803, [float(v) for v in protocol["inner_selection"]["logistic_C_grid"]])
    x_train, _, scaler = make_features(kind, dataset["fps"], dataset["desc"], all_idx, all_idx)
    estimator = LogisticRegression(C=c, class_weight="balanced", solver="liblinear", max_iter=5000, random_state=20260803)
    estimator.fit(x_train, dataset["y"])
    training_fps = [FP_GENERATOR.GetFingerprint(mol) for mol in dataset["mols"]]
    self_nearest = []
    for i, fp in enumerate(training_fps):
        other = training_fps[:i] + training_fps[i + 1 :]
        self_nearest.append(max(DataStructs.BulkTanimotoSimilarity(fp, other)) if other else 0.0)
    return {
        "model_version": "ABCB1_assay_aware_candidate_v0.3_dev",
        "status": "single_source_internal_development_candidate_not_externally_validated",
        "feature_type": kind.removesuffix("_logistic"),
        "estimator": estimator,
        "descriptor_scaler": scaler,
        "selected_C": c,
        "binary_threshold": 0.5,
        "training_fingerprints": training_fps,
        "training_rows": [
            {
                "panel_entity_id": dataset["ids"][i],
                "entity_name": dataset["names"][i],
                "standardized_inchikey": dataset["inchikeys"][i],
                "standardized_isomeric_smiles": dataset["smiles"][i],
                "target_label": int(dataset["y"][i]),
                "log2_net_efflux_ratio": float(dataset["continuous_y"][i]),
                "scaffold_group": dataset["groups"][i],
            }
            for i in range(len(dataset["y"]))
        ],
        "applicability_thresholds": {
            "self_nearest_q10": float(np.quantile(self_nearest, 0.10)),
            "self_nearest_q25": float(np.quantile(self_nearest, 0.25)),
            "derivation": "direct-panel leave-one-out nearest-neighbor Tanimoto distribution",
        },
        "protocol_sha256": sha256(PROTOCOL_PATH),
        "development_dataset_sha256": sha256(PANEL_PATH),
        "claim_boundary": "internal single-source matched-parental screening candidate; no externally validated universal ABCB1 substrate or ADC outcome claim",
    }


def predict_full_package(package: dict, smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"candidate_probability": "", "candidate_binary_prediction": "", "candidate_max_training_tanimoto": ""}
    fps = morgan_array(mol).reshape(1, -1)
    desc = descriptor_array(mol).reshape(1, -1)
    feature_type = package["feature_type"]
    if feature_type == "physchem":
        x = package["descriptor_scaler"].transform(desc)
    elif feature_type == "morgan":
        x = fps
    else:
        x = np.hstack([fps, package["descriptor_scaler"].transform(desc)])
    probability = float(package["estimator"].predict_proba(x)[0, 1])
    fp = FP_GENERATOR.GetFingerprint(mol)
    similarity = float(max(DataStructs.BulkTanimotoSimilarity(fp, package["training_fingerprints"])))
    return {
        "candidate_probability": round(probability, 8),
        "candidate_binary_prediction": int(probability >= 0.5),
        "candidate_max_training_tanimoto": round(similarity, 6),
    }


def main() -> None:
    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    allowed_statuses = {
        "locked_before_model_result_generation",
        "locked_modeling_rules_with_posthoc_audit_only_identity_correction_2026-08-03",
    }
    if protocol["status"] not in allowed_statuses:
        raise ValueError("development protocol is not locked or transparently audit-corrected")
    if sha256(PANEL_PATH) != protocol["development_dataset"]["sha256"]:
        raise ValueError("development panel hash does not match locked protocol")
    panel_rows = read_csv(PANEL_PATH)
    primary = make_dataset(panel_rows, "conservative_human_mdr1_screen_binary")
    sensitivity = make_dataset(panel_rows, "derived_human_mdr1_screen_binary")

    primary_folds, primary_details, primary_agg = repeated_classification(primary, protocol, "primary_conservative_interval")
    sensitivity_folds, sensitivity_details, sensitivity_agg = repeated_classification(sensitivity, protocol, "sensitivity_point_threshold")
    regression_folds, regression_details, regression_agg = repeated_regression(sensitivity, protocol)

    primary_metrics, pairwise_rows = cluster_bootstrap_binary(primary["y"], primary_agg, primary["groups"])
    sensitivity_metrics, _ = cluster_bootstrap_binary(sensitivity["y"], sensitivity_agg, sensitivity["groups"])
    for row in primary_metrics:
        row["analysis_set"] = "primary_conservative_interval"
        row["n_entities"] = len(primary["y"])
        row["n_positive"] = int(np.sum(primary["y"] == 1))
        row["n_negative"] = int(np.sum(primary["y"] == 0))
    for row in sensitivity_metrics:
        row["analysis_set"] = "sensitivity_point_threshold"
        row["n_entities"] = len(sensitivity["y"])
        row["n_positive"] = int(np.sum(sensitivity["y"] == 1))
        row["n_negative"] = int(np.sum(sensitivity["y"] == 0))

    regression_metric_rows = []
    for model, pred in regression_agg.items():
        for metric, value in continuous_metrics(sensitivity["continuous_y"], pred).items():
            regression_metric_rows.append(
                {
                    "analysis_set": "point_qc_eligible_continuous",
                    "model": model,
                    "metric": metric,
                    "point_estimate": value,
                    "n_entities": len(sensitivity["y"]),
                    "claim_boundary": "repeated scaffold-grouped internal development; single source",
                }
            )

    aggregate_rows = []
    for dataset_name, dataset, aggregates in (
        ("primary_conservative_interval", primary, primary_agg),
        ("sensitivity_point_threshold", sensitivity, sensitivity_agg),
    ):
        for i in range(len(dataset["y"])):
            for model, values in aggregates.items():
                aggregate_rows.append(
                    {
                        "analysis_set": dataset_name,
                        "panel_entity_id": dataset["ids"][i],
                        "entity_name": dataset["names"][i],
                        "standardized_inchikey": dataset["inchikeys"][i],
                        "scaffold_group": dataset["groups"][i],
                        "true_label": int(dataset["y"][i]),
                        "model": model,
                        "mean_outer_probability": round(float(values[i]), 8),
                        "outer_binary_prediction": int(values[i] >= 0.5),
                        "outer_repeat_count": 10,
                    }
                )

    regression_aggregate_rows = []
    for i in range(len(sensitivity["y"])):
        for model, values in regression_agg.items():
            regression_aggregate_rows.append(
                {
                    "panel_entity_id": sensitivity["ids"][i],
                    "entity_name": sensitivity["names"][i],
                    "standardized_inchikey": sensitivity["inchikeys"][i],
                    "scaffold_group": sensitivity["groups"][i],
                    "true_log2_net_ER": round(float(sensitivity["continuous_y"][i]), 8),
                    "true_net_ER": round(float(2 ** sensitivity["continuous_y"][i]), 8),
                    "model": model,
                    "mean_outer_predicted_log2_net_ER": round(float(values[i]), 8),
                    "mean_outer_predicted_net_ER": round(float(2 ** values[i]), 8),
                    "outer_repeat_count": 10,
                }
            )

    selected_model = select_final_classifier(primary_metrics)
    package = fit_full_classifier(selected_model, primary, protocol)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    joblib.dump(package, MODEL_PATH)
    frozen_package = load_package(FROZEN_MODEL_PATH)
    all_prediction_rows = []
    primary_id_to_oof = {model: {primary["ids"][i]: float(values[i]) for i in range(len(primary["y"]))} for model, values in primary_agg.items()}
    for row in panel_rows:
        candidate = predict_full_package(package, row["standardized_isomeric_smiles"])
        frozen = predict_smiles(row["standardized_isomeric_smiles"], frozen_package)
        entity_id = row["panel_entity_id"]
        all_prediction_rows.append(
            {
                "panel_entity_id": entity_id,
                "entity_name": row["entity_name_source"],
                "low_recovery_flag": row["low_recovery_flag"],
                "permeability_control_flag": row["permeability_control_flag"],
                "direct_point_label": row["derived_human_mdr1_screen_binary"],
                "direct_conservative_label": row["conservative_human_mdr1_screen_binary"],
                "net_efflux_ratio": row["human_mdr1_net_efflux_ratio_vs_mock"],
                "threshold_certainty_status": row["threshold_certainty_status"],
                "frozen_training_overlap": row["overlap_with_frozen_training"],
                "weak_vs_direct_screen_label_status": row["weak_vs_direct_screen_label_status"],
                "selected_candidate_model": selected_model,
                **candidate,
                "selected_candidate_primary_mean_outer_probability": round(primary_id_to_oof[selected_model][entity_id], 8) if entity_id in primary_id_to_oof[selected_model] else "",
                "frozen_v0_1_probability": frozen.get("calibrated_probability", ""),
                "frozen_v0_1_binary_prediction": frozen.get("binary_prediction", ""),
                "frozen_v0_1_max_training_tanimoto": frozen.get("max_training_morgan_tanimoto", ""),
                "standardized_inchikey": row["standardized_inchikey"],
                "standardized_isomeric_smiles": row["standardized_isomeric_smiles"],
                "claim_boundary": "candidate full-fit score is development use only; only repeated outer prediction is used for internal performance",
            }
        )

    paths = {
        "classification_metric_summary": OUTPUT_DIR / "abcb1_assay_aware_classification_metrics_v1.csv",
        "classification_fold_metrics": OUTPUT_DIR / "abcb1_assay_aware_classification_fold_metrics_v1.csv",
        "classification_oof_details": OUTPUT_DIR / "abcb1_assay_aware_classification_oof_details_v1.csv",
        "classification_oof_aggregate": OUTPUT_DIR / "abcb1_assay_aware_classification_oof_aggregate_v1.csv",
        "pairwise_comparison": OUTPUT_DIR / "abcb1_assay_aware_pairwise_comparison_v1.csv",
        "continuous_metric_summary": OUTPUT_DIR / "abcb1_assay_aware_continuous_metrics_v1.csv",
        "continuous_fold_metrics": OUTPUT_DIR / "abcb1_assay_aware_continuous_fold_metrics_v1.csv",
        "continuous_oof_details": OUTPUT_DIR / "abcb1_assay_aware_continuous_oof_details_v1.csv",
        "continuous_oof_aggregate": OUTPUT_DIR / "abcb1_assay_aware_continuous_oof_aggregate_v1.csv",
        "candidate_predictions": OUTPUT_DIR / "abcb1_assay_aware_candidate_predictions_v1.csv",
    }
    write_csv(paths["classification_metric_summary"], primary_metrics + sensitivity_metrics)
    write_csv(paths["classification_fold_metrics"], primary_folds + sensitivity_folds)
    write_csv(paths["classification_oof_details"], primary_details + sensitivity_details)
    write_csv(paths["classification_oof_aggregate"], aggregate_rows)
    write_csv(paths["pairwise_comparison"], pairwise_rows)
    write_csv(paths["continuous_metric_summary"], regression_metric_rows)
    write_csv(paths["continuous_fold_metrics"], regression_folds)
    write_csv(paths["continuous_oof_details"], regression_details)
    write_csv(paths["continuous_oof_aggregate"], regression_aggregate_rows)
    write_csv(paths["candidate_predictions"], all_prediction_rows)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model_version": package["model_version"],
        "status": package["status"],
        "selected_internal_candidate": selected_model,
        "selected_C": package["selected_C"],
        "model_file": MODEL_PATH.name,
        "model_sha256": sha256(MODEL_PATH),
        "protocol_path": str(PROTOCOL_PATH.relative_to(PROJECT)),
        "protocol_sha256": sha256(PROTOCOL_PATH),
        "development_panel_path": str(PANEL_PATH.relative_to(PROJECT)),
        "development_panel_sha256": sha256(PANEL_PATH),
        "primary_n": len(primary["y"]),
        "primary_positive": int(np.sum(primary["y"] == 1)),
        "primary_negative": int(np.sum(primary["y"] == 0)),
        "sensitivity_n": len(sensitivity["y"]),
        "claim_boundary": package["claim_boundary"],
        "outputs": {name: str(path.relative_to(PROJECT)) for name, path in paths.items()},
    }
    (MODEL_DIR / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = {
        **manifest,
        "primary_metric_rows": [row for row in primary_metrics if row["model"] == selected_model],
        "pairwise_rows_selected_candidate": [row for row in pairwise_rows if row["candidate_model"] == selected_model],
        "continuous_metric_rows": regression_metric_rows,
        "validation_status": "not_externally_validated",
        "next_required_action": "collect and freeze a new untouched human ABCB1 matched-parental source or prospective panel before external scoring",
    }
    summary_path = OUTPUT_DIR / "abcb1_assay_aware_development_summary_v1.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
