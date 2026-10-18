#!/usr/bin/env python3
"""Train and audit the frozen 443-row ABCB1 research screener v0.1.

The script refuses to run if the frozen input hash changes or if a direct
structure-label conflict is detected. Model selection uses only repeated
scaffold-grouped internal cross-validation. External ADC and Roth records are
pressure tests and never enter fitting, calibration, or threshold selection.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import joblib
import numpy as np
import rdkit
import scipy
import sklearn
from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
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
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler

from abcb1_screener_core_v0_1 import (
    FP_GENERATOR,
    MORGAN_BITS,
    descriptor_array,
    morgan_array,
    predict_smiles,
    standardize_smiles,
)


DATA = Path(__file__).resolve().parents[1]
PROJECT = DATA.parent
PROCESSED = DATA / "processed"
INPUT = PROCESSED / "binary_supplement_training_interface_v1.csv"
EXPECTED_INPUT_SHA256 = "a1787c8ca0d02a9f142f3ca49a07c3434a2f07cd1977b8f754a55170c2da3bc5"
MODEL_DIR = PROJECT / "模型" / "ABCB1_screener_v0_1"
OUTPUT_DIR = PROJECT / "outputs" / "20260803_abcb1_screener_v0_1"
MODEL_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RUN_DATE = "2026-08-03"
SEEDS = [20260803, 20260804, 20260805, 20260806, 20260807]
N_SPLITS = 5
FIXED_BINARY_THRESHOLD = 0.50
GATE_MCC = 0.50
GATE_BALANCED_ACCURACY = 0.75


@dataclass(frozen=True)
class ModelSpec:
    name: str
    feature_type: str
    family: str
    c: float | None = None


MODEL_SPECS = [
    ModelSpec("dummy_prior", "morgan", "dummy"),
    ModelSpec("physchem_logistic_C1_balanced", "physchem", "logistic", 1.0),
    ModelSpec("morgan_logistic_C0.03_balanced", "morgan", "logistic", 0.03),
    ModelSpec("morgan_logistic_C0.1_balanced", "morgan", "logistic", 0.1),
    ModelSpec("morgan_logistic_C0.3_balanced", "morgan", "logistic", 0.3),
    ModelSpec("morgan_logistic_C1_balanced", "morgan", "logistic", 1.0),
    ModelSpec("morgan_physchem_logistic_C0.1_balanced", "morgan_physchem", "logistic", 0.1),
    ModelSpec("extra_trees_morgan_balanced", "morgan", "extra_trees"),
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
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def json_dump(path: Path, payload: dict | list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")


def finite_or_blank(value: float) -> float | str:
    return round(float(value), 8) if math.isfinite(float(value)) else ""


def ece_score(y: np.ndarray, probability: np.ndarray, n_bins: int = 10) -> float:
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    for lower, upper in zip(edges[:-1], edges[1:]):
        if upper == 1.0:
            mask = (probability >= lower) & (probability <= upper)
        else:
            mask = (probability >= lower) & (probability < upper)
        if not np.any(mask):
            continue
        ece += float(mask.mean()) * abs(float(probability[mask].mean()) - float(y[mask].mean()))
    return ece


def binary_metrics(y: np.ndarray, probability: np.ndarray, threshold: float = 0.5) -> dict:
    prediction = (probability >= threshold).astype(int)
    both = len(set(y.tolist())) == 2
    tn = int(((y == 0) & (prediction == 0)).sum())
    fp = int(((y == 0) & (prediction == 1)).sum())
    fn = int(((y == 1) & (prediction == 0)).sum())
    tp = int(((y == 1) & (prediction == 1)).sum())
    specificity = tn / (tn + fp) if tn + fp else math.nan
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
        "f1": float(f1_score(y, prediction, zero_division=0)),
        "roc_auc": float(roc_auc_score(y, probability)) if both else math.nan,
        "average_precision": float(average_precision_score(y, probability)) if both else math.nan,
        "brier": float(brier_score_loss(y, probability)),
        "log_loss": float(log_loss(y, np.clip(probability, 1e-8, 1 - 1e-8), labels=[0, 1])),
        "ece_10bin": ece_score(y, probability),
    }


def model_features(spec: ModelSpec, x_fp: np.ndarray, x_desc: np.ndarray, scaler: StandardScaler | None) -> np.ndarray:
    if spec.feature_type == "morgan":
        return x_fp
    if spec.feature_type == "physchem":
        if scaler is None:
            raise ValueError("descriptor scaler required")
        return scaler.transform(x_desc)
    if spec.feature_type == "morgan_physchem":
        if scaler is None:
            raise ValueError("descriptor scaler required")
        return np.hstack([x_fp, scaler.transform(x_desc)])
    raise ValueError(spec.feature_type)


def fit_spec(
    spec: ModelSpec,
    x_fp: np.ndarray,
    x_desc: np.ndarray,
    y: np.ndarray,
    train_indices: np.ndarray,
    seed: int,
) -> tuple[object, StandardScaler | None]:
    scaler = None
    if spec.feature_type in {"physchem", "morgan_physchem"}:
        scaler = StandardScaler().fit(x_desc[train_indices])
    x_train = model_features(spec, x_fp[train_indices], x_desc[train_indices], scaler)
    if spec.family == "dummy":
        estimator = DummyClassifier(strategy="prior")
    elif spec.family == "logistic":
        estimator = LogisticRegression(
            C=float(spec.c),
            class_weight="balanced",
            max_iter=5000,
            solver="liblinear",
            random_state=seed,
        )
    elif spec.family == "extra_trees":
        estimator = ExtraTreesClassifier(
            n_estimators=300,
            min_samples_leaf=2,
            max_features="sqrt",
            class_weight="balanced",
            n_jobs=-1,
            random_state=seed,
        )
    else:
        raise ValueError(spec.family)
    estimator.fit(x_train, y[train_indices])
    return estimator, scaler


def predict_spec(
    spec: ModelSpec,
    estimator: object,
    scaler: StandardScaler | None,
    x_fp: np.ndarray,
    x_desc: np.ndarray,
    indices: np.ndarray,
) -> np.ndarray:
    x = model_features(spec, x_fp[indices], x_desc[indices], scaler)
    return estimator.predict_proba(x)[:, 1]


def derived_generic_scaffold(mol: Chem.Mol, structure_key: str) -> str:
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumHeavyAtoms() == 0:
        return f"ACYCLIC::{structure_key}"
    generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    return Chem.MolToSmiles(generic, canonical=True, isomericSmiles=False)


def conflict_rows(records: list[dict]) -> list[dict]:
    checks = [
        ("standardized_inchikey", "hard_exact_structure_conflict"),
        ("derived_isomeric_smiles", "isomeric_structure_conflict"),
        ("derived_connectivity_smiles", "connectivity_conflict"),
        ("preferred_name_normalized", "preferred_name_conflict"),
    ]
    out: list[dict] = []
    for field, issue in checks:
        grouped: dict[str, list[dict]] = defaultdict(list)
        for row in records:
            key = row.get(field, "")
            if key:
                grouped[str(key)].append(row)
        for key, group in grouped.items():
            labels = sorted({str(item["target_label"]) for item in group})
            if len(labels) <= 1:
                continue
            out.append(
                {
                    "issue_type": issue,
                    "matching_field": field,
                    "matching_key": key,
                    "labels": " | ".join(labels),
                    "record_count": len(group),
                    "interface_ids": " | ".join(sorted(item["interface_id"] for item in group)),
                    "model_action": "exclude_and_adjudicate",
                }
            )
    return out


def self_nearest_similarities(fps: list) -> np.ndarray:
    maxima = []
    for index, fp in enumerate(fps):
        others = fps[:index] + fps[index + 1 :]
        maxima.append(max(DataStructs.BulkTanimotoSimilarity(fp, others)) if others else 0.0)
    return np.asarray(maxima, dtype=float)


def max_similarity_to_train(fp, train_fps: list) -> float:
    return float(max(DataStructs.BulkTanimotoSimilarity(fp, train_fps))) if train_fps else 0.0


def run_grouped_cv(
    records: list[dict],
    x_fp: np.ndarray,
    x_desc: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    split_name: str,
    specs: list[ModelSpec],
) -> tuple[list[dict], dict[str, np.ndarray], dict[str, np.ndarray], list[dict]]:
    fold_rows: list[dict] = []
    probability_sum = {spec.name: np.zeros(len(y), dtype=float) for spec in specs}
    probability_count = {spec.name: np.zeros(len(y), dtype=int) for spec in specs}
    support_sum = np.zeros(len(y), dtype=float)
    support_count = np.zeros(len(y), dtype=int)
    fps = [FP_GENERATOR.GetFingerprint(Chem.MolFromSmiles(row["derived_isomeric_smiles"])) for row in records]

    for repeat, seed in enumerate(SEEDS, start=1):
        splitter = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        for fold, (train_idx, test_idx) in enumerate(splitter.split(x_fp, y, groups), start=1):
            train_fps = [fps[index] for index in train_idx]
            for index in test_idx:
                support_sum[index] += max_similarity_to_train(fps[index], train_fps)
                support_count[index] += 1
            for spec in specs:
                estimator, scaler = fit_spec(spec, x_fp, x_desc, y, train_idx, seed)
                probability = predict_spec(spec, estimator, scaler, x_fp, x_desc, test_idx)
                probability_sum[spec.name][test_idx] += probability
                probability_count[spec.name][test_idx] += 1
                metric = binary_metrics(y[test_idx], probability, FIXED_BINARY_THRESHOLD)
                fold_rows.append(
                    {
                        "split_name": split_name,
                        "repeat": repeat,
                        "seed": seed,
                        "fold": fold,
                        "model": spec.name,
                        "feature_type": spec.feature_type,
                        "n_train": len(train_idx),
                        "n_test": len(test_idx),
                        "train_group_count": len(set(groups[train_idx].tolist())),
                        "test_group_count": len(set(groups[test_idx].tolist())),
                        **{key: finite_or_blank(value) if isinstance(value, float) else value for key, value in metric.items()},
                    }
                )

    oof = {name: probability_sum[name] / probability_count[name] for name in probability_sum}
    support = support_sum / support_count
    summary: list[dict] = []
    for spec in specs:
        relevant = [row for row in fold_rows if row["model"] == spec.name]
        aggregated = binary_metrics(y, oof[spec.name], FIXED_BINARY_THRESHOLD)
        summary.append(
            {
                "split_name": split_name,
                "model": spec.name,
                "feature_type": spec.feature_type,
                "valid_folds": len(relevant),
                "mean_fold_mcc": round(float(np.mean([float(row["mcc"]) for row in relevant])), 8),
                "sd_fold_mcc": round(float(np.std([float(row["mcc"]) for row in relevant])), 8),
                "mean_fold_balanced_accuracy": round(float(np.mean([float(row["balanced_accuracy"]) for row in relevant])), 8),
                "sd_fold_balanced_accuracy": round(float(np.std([float(row["balanced_accuracy"]) for row in relevant])), 8),
                "mean_fold_roc_auc": round(float(np.mean([float(row["roc_auc"]) for row in relevant])), 8),
                "sd_fold_roc_auc": round(float(np.std([float(row["roc_auc"]) for row in relevant])), 8),
                **{f"repeated_oof_{key}": finite_or_blank(value) if isinstance(value, float) else value for key, value in aggregated.items()},
                "gate_rule": f"mean fold MCC >= {GATE_MCC:.2f} and mean fold balanced accuracy >= {GATE_BALANCED_ACCURACY:.2f}",
                "gate_status": "pass" if np.mean([float(row["mcc"]) for row in relevant]) >= GATE_MCC and np.mean([float(row["balanced_accuracy"]) for row in relevant]) >= GATE_BALANCED_ACCURACY else "fail",
            }
        )
    support_rows = [
        {
            "interface_id": record["interface_id"],
            "mean_oof_max_training_tanimoto": round(float(support[index]), 6),
        }
        for index, record in enumerate(records)
    ]
    return fold_rows, oof, {"support": support}, summary


def run_random_reference_cv(
    x_fp: np.ndarray,
    x_desc: np.ndarray,
    y: np.ndarray,
    selected_spec: ModelSpec,
) -> dict:
    probabilities = np.zeros(len(y), dtype=float)
    counts = np.zeros(len(y), dtype=int)
    fold_metrics = []
    for seed in SEEDS:
        splitter = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        for train_idx, test_idx in splitter.split(x_fp, y):
            estimator, scaler = fit_spec(selected_spec, x_fp, x_desc, y, train_idx, seed)
            probability = predict_spec(selected_spec, estimator, scaler, x_fp, x_desc, test_idx)
            probabilities[test_idx] += probability
            counts[test_idx] += 1
            fold_metrics.append(binary_metrics(y[test_idx], probability, FIXED_BINARY_THRESHOLD))
    oof = probabilities / counts
    metric = binary_metrics(y, oof, FIXED_BINARY_THRESHOLD)
    return {
        "split_name": "optimistic_random_stratified_reference",
        "model": selected_spec.name,
        "valid_folds": len(fold_metrics),
        "mean_fold_mcc": round(float(np.mean([row["mcc"] for row in fold_metrics])), 8),
        "mean_fold_balanced_accuracy": round(float(np.mean([row["balanced_accuracy"] for row in fold_metrics])), 8),
        "mean_fold_roc_auc": round(float(np.mean([row["roc_auc"] for row in fold_metrics])), 8),
        **{f"repeated_oof_{key}": finite_or_blank(value) if isinstance(value, float) else value for key, value in metric.items()},
        "claim_boundary": "optimistic comparator only; not used for model selection or claim",
    }


def fit_calibrator(y: np.ndarray, probability: np.ndarray) -> LogisticRegression:
    clipped = np.clip(probability, 1e-6, 1 - 1e-6)
    logits = np.log(clipped / (1.0 - clipped)).reshape(-1, 1)
    calibrator = LogisticRegression(C=1e6, solver="lbfgs", max_iter=5000)
    calibrator.fit(logits, y)
    return calibrator


def apply_calibrator(calibrator: LogisticRegression, probability: np.ndarray) -> np.ndarray:
    clipped = np.clip(probability, 1e-6, 1 - 1e-6)
    logits = np.log(clipped / (1.0 - clipped)).reshape(-1, 1)
    return calibrator.predict_proba(logits)[:, 1]


def threshold_sweep(y: np.ndarray, probability: np.ndarray) -> list[dict]:
    rows = []
    for threshold in np.round(np.arange(0.20, 0.801, 0.01), 2):
        metric = binary_metrics(y, probability, float(threshold))
        rows.append({key: finite_or_blank(value) if isinstance(value, float) else value for key, value in metric.items()})
    return rows


def predict_external(package: dict) -> tuple[list[dict], dict]:
    output: list[dict] = []
    anchor_path = PROCESSED / "adc_direct_numeric_binary_anchors_v1.csv"
    if anchor_path.exists():
        for row in read_csv(anchor_path):
            if row.get("transporter_gene") != "ABCB1":
                continue
            prediction = predict_smiles(row.get("standardized_isomeric_smiles", ""), package)
            output.append(
                {
                    "evaluation_set": "adc_direct_numeric_anchor",
                    "record_id": row.get("anchor_id", ""),
                    "entity_name": row.get("entity_name", ""),
                    "entity_type": row.get("entity_type", ""),
                    "source_endpoint": "bidirectional_transport_numeric",
                    "source_value": row.get("no_inhibitor_efflux_ratio", ""),
                    "source_binary_label": row.get("derived_binary_label", ""),
                    "truth_set_eligible": row.get("truth_set_eligible", ""),
                    "source_url": row.get("source_urls", ""),
                    **{key: prediction.get(key, "") for key in ["standardized_inchikey", "calibrated_probability", "binary_prediction", "max_training_morgan_tanimoto", "applicability_band", "screening_call"]},
                    "claim_boundary": "external ADC pressure test only; never used in fitting or threshold selection",
                }
            )

    roth_path = PROCESSED / "external_truth_set_v2_pgp_bcrp_continuous_identity_audited.csv"
    registry_path = PROCESSED / "compound_registry_standardized_v1.csv"
    roth_probabilities = []
    roth_ratios = []
    if roth_path.exists() and registry_path.exists():
        registry = {row.get("reference_id", ""): row for row in read_csv(registry_path)}
        for row in read_csv(roth_path):
            if row.get("transporter_gene") != "ABCB1" or row.get("structure_based_external_metric_use") != "include":
                continue
            structure = registry.get(row.get("compound_id", ""), {})
            prediction = predict_smiles(structure.get("standardized_isomeric_smiles", ""), package)
            rr = float(row["resistance_ratio"]) if row.get("resistance_ratio", "") else math.nan
            if prediction.get("calibrated_probability", "") != "" and math.isfinite(rr):
                roth_probabilities.append(float(prediction["calibrated_probability"]))
                roth_ratios.append(rr)
            output.append(
                {
                    "evaluation_set": "roth_continuous_rr_identity_safe",
                    "record_id": row.get("assay_record_id", ""),
                    "entity_name": row.get("paper_compound_name", ""),
                    "entity_type": "free_payload",
                    "source_endpoint": "transfected_cell_cytotoxicity_resistance_ratio",
                    "source_value": row.get("resistance_ratio", ""),
                    "source_binary_label": "",
                    "truth_set_eligible": "continuous_directional_stress_test_only",
                    "source_url": row.get("source_url", ""),
                    **{key: prediction.get(key, "") for key in ["standardized_inchikey", "calibrated_probability", "binary_prediction", "max_training_morgan_tanimoto", "applicability_band", "screening_call"]},
                    "claim_boundary": "continuous functional RR is not a binary direct-transport truth label",
                }
            )
    if len(roth_probabilities) >= 3:
        rho, pvalue = spearmanr(roth_probabilities, roth_ratios)
        roth_summary = {
            "n": len(roth_probabilities),
            "spearman_probability_vs_resistance_ratio": finite_or_blank(float(rho)),
            "spearman_p_value": finite_or_blank(float(pvalue)),
            "expected_direction": "positive",
            "direction_gate": "pass" if rho > 0 else "fail",
            "claim_boundary": "directional pressure test only; assay endpoint differs from direct transport",
        }
    else:
        roth_summary = {"n": len(roth_probabilities), "status": "insufficient_records"}
    return output, roth_summary


def main() -> None:
    input_sha = sha256_file(INPUT)
    if input_sha != EXPECTED_INPUT_SHA256:
        raise RuntimeError(f"frozen input hash mismatch: expected {EXPECTED_INPUT_SHA256}, observed {input_sha}")

    source_rows = read_csv(INPUT)
    source_rows = [
        row
        for row in source_rows
        if row.get("transporter_gene") == "ABCB1"
        and row.get("target_role") == "candidate_supervised"
        and row.get("target_label") in {"0", "1"}
    ]
    if len(source_rows) != 443:
        raise RuntimeError(f"expected 443 frozen ABCB1 rows, observed {len(source_rows)}")

    records: list[dict] = []
    molecules: list[Chem.Mol] = []
    for row in source_rows:
        frozen_smiles = row.get("standardized_isomeric_smiles", "")
        mol = Chem.MolFromSmiles(frozen_smiles)
        if mol is None:
            raise RuntimeError(f"frozen standardized structure parse failure: {row['interface_id']}")
        reproduced = standardize_smiles(frozen_smiles)
        frozen_inchikey = row.get("standardized_inchikey", "") or Chem.MolToInchiKey(mol)
        frozen_isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        frozen_connectivity = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        structure_key = frozen_inchikey or frozen_connectivity
        record = {
            **row,
            "preferred_name_normalized": " ".join(row.get("preferred_name", "").upper().split()),
            "input_standardized_inchikey": row.get("standardized_inchikey", ""),
            "derived_standardized_inchikey": frozen_inchikey,
            "derived_isomeric_smiles": frozen_isomeric,
            "derived_connectivity_smiles": frozen_connectivity,
            "standardizer_reproduced_inchikey": reproduced.get("standardized_inchikey", ""),
            "inchikey_reproduction_status": "match" if reproduced.get("standardized_inchikey", "") == frozen_inchikey else "non_idempotent_tautomer_review_frozen_structure_retained",
            "exact_scaffold_group": row.get("murcko_scaffold_smiles", "") or f"ACYCLIC::{structure_key}",
            "generic_scaffold_group": derived_generic_scaffold(mol, structure_key),
        }
        records.append(record)
        molecules.append(mol)

    conflicts = conflict_rows(records)
    write_csv(PROCESSED / "abcb1_screener_v0_1_conflict_audit.csv", conflicts, ["issue_type", "matching_field", "matching_key", "labels", "record_count", "interface_ids", "model_action"])
    if conflicts:
        raise RuntimeError(f"direct structure-label conflicts detected: {len(conflicts)}")

    y = np.asarray([int(row["target_label"]) for row in records], dtype=int)
    x_fp = np.vstack([morgan_array(mol) for mol in molecules])
    x_desc = np.vstack([descriptor_array(mol) for mol in molecules])
    fps = [FP_GENERATOR.GetFingerprint(mol) for mol in molecules]
    exact_groups = np.asarray([row["exact_scaffold_group"] for row in records], dtype=object)
    generic_groups = np.asarray([row["generic_scaffold_group"] for row in records], dtype=object)

    fold_rows, exact_oof, exact_support_map, exact_summary = run_grouped_cv(
        records, x_fp, x_desc, y, exact_groups, "primary_exact_murcko_grouped", MODEL_SPECS
    )
    non_dummy = [row for row in exact_summary if row["model"] != "dummy_prior"]
    selected_summary = sorted(
        non_dummy,
        key=lambda row: (
            float(row["mean_fold_mcc"]),
            float(row["mean_fold_balanced_accuracy"]),
            float(row["mean_fold_roc_auc"]),
            -float(row["repeated_oof_brier"]),
        ),
        reverse=True,
    )[0]
    selected_spec = next(spec for spec in MODEL_SPECS if spec.name == selected_summary["model"])

    generic_fold_rows, generic_oof, _, generic_summary = run_grouped_cv(
        records, x_fp, x_desc, y, generic_groups, "secondary_generic_murcko_grouped", [selected_spec]
    )
    fold_rows.extend(generic_fold_rows)
    all_summary = exact_summary + generic_summary
    random_reference = run_random_reference_cv(x_fp, x_desc, y, selected_spec)

    selected_raw_oof = exact_oof[selected_spec.name]
    calibrator = fit_calibrator(y, selected_raw_oof)
    calibrated_oof = apply_calibrator(calibrator, selected_raw_oof)
    calibrated_metrics = binary_metrics(y, calibrated_oof, FIXED_BINARY_THRESHOLD)
    sweep = threshold_sweep(y, calibrated_oof)
    best_internal = sorted(
        sweep,
        key=lambda row: (float(row["mcc"]), float(row["balanced_accuracy"]), -abs(float(row["threshold"]) - 0.5)),
        reverse=True,
    )[0]
    for row in sweep:
        row["selected_for_deployment"] = "yes_fixed_preregistered" if float(row["threshold"]) == FIXED_BINARY_THRESHOLD else "no"
        row["best_internal_exploratory"] = "yes" if float(row["threshold"]) == float(best_internal["threshold"]) else "no"
        row["claim_boundary"] = "internal OOF diagnostic; the exploratory optimum is not deployed"

    self_nearest = self_nearest_similarities(fps)
    applicability = {
        "similarity_metric": "Morgan radius 2, 2048 bits, no chirality, Tanimoto",
        "self_nearest_q10": float(np.quantile(self_nearest, 0.10)),
        "self_nearest_q25": float(np.quantile(self_nearest, 0.25)),
        "self_nearest_median": float(np.quantile(self_nearest, 0.50)),
        "outside_rule": "query maximum training similarity < self-nearest q10 -> abstain",
        "caution_rule": "q10 <= similarity < q25 -> caution",
        "within_rule": "similarity >= q25 -> within empirical support; not proof of endpoint transferability",
    }

    all_indices = np.arange(len(y))
    final_estimator, final_scaler = fit_spec(selected_spec, x_fp, x_desc, y, all_indices, SEEDS[0])
    package = {
        "package_version": "ABCB1_screener_v0.1",
        "run_date": RUN_DATE,
        "transporter_gene": "ABCB1",
        "prediction_object": "standardized released or free small-molecule structure",
        "selected_model": selected_spec.name,
        "feature_type": selected_spec.feature_type,
        "estimator": final_estimator,
        "descriptor_scaler": final_scaler,
        "calibrator": calibrator,
        "binary_threshold": FIXED_BINARY_THRESHOLD,
        "binary_threshold_status": "fixed_not_external_calibrated",
        "applicability_thresholds": applicability,
        "training_fingerprints": fps,
        "training_structure_keys": [row["derived_standardized_inchikey"] for row in records],
        "training_structure_map": {row["derived_standardized_inchikey"]: row["derived_isomeric_smiles"] for row in records},
        "training_input_sha256": input_sha,
        "training_n": len(records),
        "training_n_positive": int(y.sum()),
        "training_n_negative": int((1 - y).sum()),
        "source_scope": sorted({row.get("source_database", "") for row in records}),
        "claim_status": "development_research_screener_gate_passed" if selected_summary["gate_status"] == "pass" else "development_research_screener_gate_failed",
        "claim_boundary": "research screening only; not an intact-ADC, exposure, resistance, efficacy, or universal assay-context predictor",
        "runtime": {
            "python": platform.python_version(),
            "rdkit": rdkit.__version__,
            "scikit_learn": sklearn.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "joblib": joblib.__version__,
        },
    }
    model_path = MODEL_DIR / "ABCB1_screener_v0_1.joblib"
    joblib.dump(package, model_path, compress=3)

    support = exact_support_map["support"]
    q10 = applicability["self_nearest_q10"]
    q25 = applicability["self_nearest_q25"]
    oof_rows: list[dict] = []
    for index, row in enumerate(records):
        calibrated = float(calibrated_oof[index])
        prediction = int(calibrated >= FIXED_BINARY_THRESHOLD)
        similarity = float(support[index])
        if similarity < q10:
            support_band = "outside_empirical_support_in_grouped_cv"
        elif similarity < q25:
            support_band = "low_empirical_support_in_grouped_cv"
        else:
            support_band = "within_empirical_support_in_grouped_cv"
        oof_rows.append(
            {
                "interface_id": row["interface_id"],
                "preferred_name": row["preferred_name"],
                "target_label": int(row["target_label"]),
                "source_database": row["source_database"],
                "derived_standardized_inchikey": row["derived_standardized_inchikey"],
                "exact_scaffold_group": row["exact_scaffold_group"],
                "generic_scaffold_group": row["generic_scaffold_group"],
                "selected_model": selected_spec.name,
                "repeated_exact_scaffold_oof_raw_probability": round(float(selected_raw_oof[index]), 8),
                "repeated_exact_scaffold_oof_calibrated_probability": round(calibrated, 8),
                "binary_threshold": FIXED_BINARY_THRESHOLD,
                "oof_binary_prediction": prediction,
                "oof_correct": "yes" if prediction == int(row["target_label"]) else "no",
                "mean_oof_max_training_tanimoto": round(similarity, 6),
                "oof_support_band": support_band,
                "claim_boundary": "internal repeated scaffold-grouped OOF only",
            }
        )

    external_rows, roth_summary = predict_external(package)

    model_summary_rows = []
    for row in all_summary:
        row = dict(row)
        row["selected_model"] = "yes" if row["split_name"] == "primary_exact_murcko_grouped" and row["model"] == selected_spec.name else "no"
        row["selection_rule"] = "highest primary mean fold MCC; ties by balanced accuracy, ROC-AUC, then lower Brier"
        row["claim_boundary"] = "internal cross-validation only"
        model_summary_rows.append(row)
    model_summary_rows.append(random_reference)

    conflict_summary = {
        "hard_exact_structure_conflict_groups": sum(row["issue_type"] == "hard_exact_structure_conflict" for row in conflicts),
        "isomeric_structure_conflict_groups": sum(row["issue_type"] == "isomeric_structure_conflict" for row in conflicts),
        "connectivity_conflict_groups": sum(row["issue_type"] == "connectivity_conflict" for row in conflicts),
        "preferred_name_conflict_groups": sum(row["issue_type"] == "preferred_name_conflict" for row in conflicts),
    }
    selected_gate_pass = selected_summary["gate_status"] == "pass"
    final_status = "GO_bounded_research_screener" if selected_gate_pass else "NO_GO_for_validated_binary_claim_GO_for_bounded_screening"
    summary = {
        "run_date": RUN_DATE,
        "status": final_status,
        "training_input": str(INPUT),
        "training_input_sha256": input_sha,
        "training_rows": len(records),
        "label_counts": {"0": int((y == 0).sum()), "1": int((y == 1).sum())},
        "unique_input_inchikeys": len({row["input_standardized_inchikey"] for row in records}),
        "unique_derived_inchikeys": len({row["derived_standardized_inchikey"] for row in records}),
        "exact_scaffold_groups": len(set(exact_groups.tolist())),
        "generic_scaffold_groups": len(set(generic_groups.tolist())),
        "single_source_scope": sorted({row["source_database"] for row in records}),
        "source_holdout_status": "not_identifiable_all_443_rows_from_one_source",
        "conflict_summary": conflict_summary,
        "selected_model": selected_spec.name,
        "selected_primary_cv": selected_summary,
        "selected_calibrated_repeated_oof_metrics_fixed_0_5": {key: finite_or_blank(value) if isinstance(value, float) else value for key, value in calibrated_metrics.items()},
        "exploratory_best_internal_threshold": best_internal,
        "deployed_threshold": FIXED_BINARY_THRESHOLD,
        "deployed_threshold_reason": "preserve fixed binary interface; exploratory OOF optimum not deployed",
        "applicability": applicability,
        "roth_external_continuous_pressure_test": roth_summary,
        "claim_boundary": package["claim_boundary"],
        "model_path": str(model_path),
    }

    write_csv(PROCESSED / "abcb1_screener_v0_1_training_audit.csv", records)
    write_csv(PROCESSED / "abcb1_screener_v0_1_cv_fold_metrics.csv", fold_rows)
    write_csv(PROCESSED / "abcb1_screener_v0_1_model_summary.csv", model_summary_rows)
    write_csv(PROCESSED / "abcb1_screener_v0_1_oof_predictions.csv", oof_rows)
    write_csv(PROCESSED / "abcb1_screener_v0_1_threshold_sweep.csv", sweep)
    write_csv(PROCESSED / "abcb1_screener_v0_1_external_pressure_predictions.csv", external_rows)
    json_dump(PROCESSED / "abcb1_screener_v0_1_summary.json", summary)

    artifact_paths = [
        model_path,
        PROCESSED / "abcb1_screener_v0_1_conflict_audit.csv",
        PROCESSED / "abcb1_screener_v0_1_training_audit.csv",
        PROCESSED / "abcb1_screener_v0_1_cv_fold_metrics.csv",
        PROCESSED / "abcb1_screener_v0_1_model_summary.csv",
        PROCESSED / "abcb1_screener_v0_1_oof_predictions.csv",
        PROCESSED / "abcb1_screener_v0_1_threshold_sweep.csv",
        PROCESSED / "abcb1_screener_v0_1_external_pressure_predictions.csv",
        PROCESSED / "abcb1_screener_v0_1_summary.json",
    ]
    manifest = {
        "package_version": "ABCB1_screener_v0.1",
        "generated_on": RUN_DATE,
        "artifacts": [
            {"path": str(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in artifact_paths
        ],
    }
    json_dump(MODEL_DIR / "manifest.json", manifest)

    model_card = f"""# ABCB1 substrate research screener v0.1

日期：{RUN_DATE}

## 结论

本版本是可运行的研究筛选器，不是已验证的 ADC 专用二分类器。训练对象为冻结接口中 443 个 ABCB1 候选监督结构（291 阳性、152 阴性），输入 SHA256 为 `{input_sha}`。训练集内部直接结构标签冲突为 0。

内部选择模型：`{selected_spec.name}`。主评估为 5 次重复、每次 5 折的 exact-Murcko 分组交叉验证。主门槛是 MCC ≥ {GATE_MCC:.2f} 且 balanced accuracy ≥ {GATE_BALANCED_ACCURACY:.2f}；本轮状态为 **{selected_summary['gate_status']}**。固定 0.5 阈值的校准 OOF 指标：MCC {calibrated_metrics['mcc']:.3f}、BAcc {calibrated_metrics['balanced_accuracy']:.3f}、ROC-AUC {calibrated_metrics['roc_auc']:.3f}、Brier {calibrated_metrics['brier']:.3f}。

## 模型判断层

输入是游离载荷或实际释放小分子的 SMILES。输出是 ABCB1 结构性底物倾向概率、固定 0.5 二分类、最大训练相似性和适用域状态。它不判断完整 ADC、释放效率、胞内游离暴露、耐药倍数、旁观者效应或临床疗效。

## 适用域

- 最大训练 Tanimoto < {applicability['self_nearest_q10']:.3f}：拒判；
- {applicability['self_nearest_q10']:.3f}–{applicability['self_nearest_q25']:.3f}：谨慎；
- ≥ {applicability['self_nearest_q25']:.3f}：仅表示处于经验结构支持域，不表示实验端点可迁移。

## 仍未解决的限制

443 条记录全部来自 ChEMBL/Daood 同一整理来源，无法执行真正的来源留出验证；标签仍是作者整理候选而非逐行 Tier-A 双向转运真值。Roth 连续 RR 只作方向压力测试，结果见 summary，不能转成二元真值。任何 ADC 对象即使给出高概率，只要处于适用域外，都必须按 `abstain_outside_domain` 处理。

## 复现

```bash
./.venv/bin/python 数据/scripts/train_abcb1_screener_v0_1.py
./.venv/bin/python 模型/ABCB1_screener_v0_1/predict_abcb1.py --smiles 'CCO'
```
"""
    (MODEL_DIR / "MODEL_CARD.md").write_text(model_card, encoding="utf-8")

    reproducibility_paths = artifact_paths + [
        MODEL_DIR / "MODEL_CARD.md",
        MODEL_DIR / "predict_abcb1.py",
        Path(__file__).resolve(),
        Path(__file__).resolve().parent / "abcb1_screener_core_v0_1.py",
    ]
    reproducibility_manifest = {
        "package_version": "ABCB1_screener_v0.1",
        "generated_on": RUN_DATE,
        "artifacts": [
            {"path": str(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in reproducibility_paths
        ],
    }
    json_dump(MODEL_DIR / "manifest.json", reproducibility_manifest)

    print(json.dumps(summary, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
