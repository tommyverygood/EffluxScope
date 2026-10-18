#!/usr/bin/env python3
"""Run the pre-locked source-aware external benchmark for ABCB1 v0.1.

This script never refits the frozen model and never tunes the 0.5 threshold.
Direct-source calls are used only as directional external diagnostics; they
are not written back into the frozen training interface as universal labels.
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
from scipy.stats import binomtest, rankdata, spearmanr
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)


DATA = Path(__file__).resolve().parents[1]
PROJECT = DATA.parent
PROCESSED = DATA / "processed"
MODEL_DIR = PROJECT / "模型" / "ABCB1_screener_v0_1"
MODEL_PATH = MODEL_DIR / "ABCB1_screener_v0_1.joblib"
LOCK_PATH = PROCESSED / "abcb1_q1_external_benchmark_protocol_lock_v1.json"
sys.path.insert(0, str(Path(__file__).resolve().parent))

from abcb1_screener_core_v0_1 import FP_GENERATOR, predict_smiles  # noqa: E402


RUN_DATE = "2026-08-03"
BOOTSTRAP_ITERATIONS = 2000
SEED = 20260803
EXPECTED_MODEL_SHA256 = "fbc92dfa99583386ba124141fed5e77c6ecf8a59850301a6533db3883b298ac3"

SOURCE_FILES = [
    "langthaler_2024_human_abcb1_reference_panel_observations_v1.csv",
    "pmda_regulatory_human_transporter_contrast_panel_observations_v1.csv",
    "independent_human_transporter_validation_panel_observations_v1.csv",
    "wegler_2021_human_abcb1_crispr_panel_observations_v1.csv",
    "independent_human_transporter_negative_contrast_panel_observations_v1.csv",
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


def safe_float(value) -> float | None:
    try:
        if value is None or str(value).strip() == "":
            return None
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def finite_or_blank(value: float) -> float | str:
    return round(float(value), 8) if math.isfinite(float(value)) else ""


def normalize_source_call(call: str) -> tuple[str, str]:
    text = (call or "").strip()
    if "human_ABCB1_positive_reference" in text or "human_ABCB1_positive_low_affinity_context" in text or text == "likely":
        return "positive_like", "1"
    if "human_ABCB1_negative_boundary" in text or "human_ABCB1_weak_or_not_supported" in text or text == "unlikely":
        return "negative_like", "0"
    return "context_only", ""


def source_url_registry() -> dict[str, str]:
    registry: dict[str, str] = {}
    for filename in SOURCE_FILES:
        for row in read_csv(PROCESSED / filename):
            source_id = row.get("source_id", "")
            url = row.get("source_url", "")
            if source_id and url:
                registry[source_id] = url
    return registry


def observation_index() -> tuple[dict[str, str], dict[tuple[str, str], list[tuple[float, str, str]]]]:
    smiles_by_inchikey: dict[str, str] = {}
    numeric_by_entity_source: dict[tuple[str, str], list[tuple[float, str, str]]] = defaultdict(list)

    for filename in SOURCE_FILES:
        rows = read_csv(PROCESSED / filename)
        for row in rows:
            if row.get("transporter_gene") != "ABCB1":
                continue
            inchikey = row.get("standardized_inchikey", "")
            smiles = row.get("standardized_isomeric_smiles", "")
            source_id = row.get("source_id", "")
            if inchikey and smiles:
                smiles_by_inchikey.setdefault(inchikey, smiles)
            if not inchikey or not source_id:
                continue

            inhibitor = (row.get("inhibitor", "") or "").strip().lower()
            is_transporter_specific_validation = (
                filename == "independent_human_transporter_validation_panel_observations_v1.csv"
                and safe_float(row.get("transporter_specific_efflux_ratio")) is not None
            )
            is_baseline = inhibitor in {"", "none", "none_required_matched_canine_mdr1_ko_control"}
            # The Zhang validation paper reports a transporter-specific ratio derived
            # from the matched tariquidar condition.  This is the locked priority
            # endpoint, not an inhibited baseline measurement, so retain it.
            is_baseline = is_baseline or is_transporter_specific_validation
            if not is_baseline:
                continue

            candidates: list[tuple[str, float | None]] = []
            if filename == "independent_human_transporter_validation_panel_observations_v1.csv":
                candidates = [
                    ("transporter_specific_efflux_ratio", safe_float(row.get("transporter_specific_efflux_ratio"))),
                    ("efflux_ratio_mean", safe_float(row.get("efflux_ratio_mean"))),
                ]
            elif filename == "pmda_regulatory_human_transporter_contrast_panel_observations_v1.csv":
                candidates = [
                    ("corrected_efflux_ratio_mean", safe_float(row.get("corrected_efflux_ratio_mean"))),
                    ("efflux_ratio_mean", safe_float(row.get("efflux_ratio_mean"))),
                ]
            elif filename == "independent_human_transporter_negative_contrast_panel_observations_v1.csv":
                candidates = [
                    ("net_flux_ratio", safe_float(row.get("net_flux_ratio"))),
                    ("efflux_ratio", safe_float(row.get("efflux_ratio"))),
                ]
            else:
                candidates = [("efflux_ratio_mean", safe_float(row.get("efflux_ratio_mean")))]

            selected = next(((name, value) for name, value in candidates if value is not None), None)
            if selected is not None:
                numeric_by_entity_source[(inchikey, source_id)].append((selected[1], selected[0], row.get("observation_id", "")))

    return smiles_by_inchikey, numeric_by_entity_source


def training_baseline(package: dict) -> tuple[list[int], list]:
    audit = read_csv(PROCESSED / "abcb1_screener_v0_1_training_audit.csv")
    if len(audit) != package["training_n"]:
        raise RuntimeError("training audit length does not match frozen package")
    keys = [row.get("derived_standardized_inchikey", "") for row in audit]
    if keys != package["training_structure_keys"]:
        raise RuntimeError("training audit order/key mismatch")
    labels = [int(row["target_label"]) for row in audit]
    return labels, package["training_fingerprints"]


def similarity_baselines(mol: Chem.Mol, training_labels: list[int], training_fps: list) -> dict:
    fp = FP_GENERATOR.GetFingerprint(mol)
    similarities = np.asarray(DataStructs.BulkTanimotoSimilarity(fp, training_fps), dtype=float)
    order = np.argsort(-similarities)
    nearest = int(order[0])
    top = order[:5]
    weights = np.maximum(similarities[top], 1e-8)
    return {
        "nearest_training_index": nearest,
        "nearest_training_tanimoto": round(float(similarities[nearest]), 6),
        "one_nn_probability": float(training_labels[nearest]),
        "knn5_weighted_probability": float(np.average(np.asarray(training_labels)[top], weights=weights)),
    }


def metric_dict(y: np.ndarray, probability: np.ndarray, threshold: float = 0.5) -> dict:
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
        "roc_auc": float(roc_auc_score(y, probability)) if both else math.nan,
        "average_precision": float(average_precision_score(y, probability)) if both else math.nan,
        "brier": float(brier_score_loss(y, probability)),
    }


def bootstrap_ci(y: np.ndarray, probability: np.ndarray, rng: np.random.Generator) -> dict[str, tuple[float, float]]:
    negative = np.flatnonzero(y == 0)
    positive = np.flatnonzero(y == 1)
    storage: dict[str, list[float]] = defaultdict(list)
    for _ in range(BOOTSTRAP_ITERATIONS):
        sampled = np.concatenate(
            [
                rng.choice(negative, size=len(negative), replace=True),
                rng.choice(positive, size=len(positive), replace=True),
            ]
        )
        metric = metric_dict(y[sampled], probability[sampled])
        for name in ["mcc", "balanced_accuracy", "sensitivity", "specificity", "roc_auc", "average_precision", "brier"]:
            value = metric[name]
            if math.isfinite(value):
                storage[name].append(value)
    return {
        name: (float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)))
        for name, values in storage.items()
        if values
    }


def paired_model_comparison(rows: list[dict]) -> list[dict]:
    """Paired entity bootstrap for frozen v0.1 versus the transparent 1-NN rule."""
    eligible = [row for row in rows if row.get("directional_label") in {"0", "1"}]
    y = np.asarray([int(row["directional_label"]) for row in eligible], dtype=int)
    frozen = np.asarray([float(row["frozen_probability"]) for row in eligible], dtype=float)
    one_nn = np.asarray([float(row["one_nn_probability"]) for row in eligible], dtype=float)
    if len(y) == 0 or len(set(y.tolist())) < 2:
        return []

    point_frozen = metric_dict(y, frozen)
    point_one_nn = metric_dict(y, one_nn)
    metric_names = ["mcc", "balanced_accuracy", "roc_auc", "brier"]
    storage: dict[str, list[float]] = defaultdict(list)
    negative = np.flatnonzero(y == 0)
    positive = np.flatnonzero(y == 1)
    rng = np.random.default_rng(SEED + 17)
    for _ in range(BOOTSTRAP_ITERATIONS):
        sampled = np.concatenate(
            [
                rng.choice(negative, size=len(negative), replace=True),
                rng.choice(positive, size=len(positive), replace=True),
            ]
        )
        frozen_metric = metric_dict(y[sampled], frozen[sampled])
        one_nn_metric = metric_dict(y[sampled], one_nn[sampled])
        for name in metric_names:
            difference = frozen_metric[name] - one_nn_metric[name]
            if math.isfinite(difference):
                storage[name].append(float(difference))

    output: list[dict] = []
    for name in metric_names:
        values = storage[name]
        lower, upper = (float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)))
        output.append(
            {
                "comparison": "frozen_ABCB1_v0_1_minus_one_nearest_neighbor",
                "metric": name,
                "frozen_value": finite_or_blank(point_frozen[name]),
                "one_nn_value": finite_or_blank(point_one_nn[name]),
                "paired_difference": finite_or_blank(point_frozen[name] - point_one_nn[name]),
                "ci95_lower": finite_or_blank(lower),
                "ci95_upper": finite_or_blank(upper),
                "n": len(y),
                "bootstrap_iterations": BOOTSTRAP_ITERATIONS,
                "direction_note": "positive favors frozen" if name != "brier" else "negative favors frozen because lower Brier is better",
                "claim_boundary": "paired historical external diagnostic; no model selection or tuning on these results",
            }
        )

    frozen_correct = (frozen >= 0.5).astype(int) == y
    one_nn_correct = (one_nn >= 0.5).astype(int) == y
    frozen_only_correct = int((frozen_correct & ~one_nn_correct).sum())
    one_nn_only_correct = int((~frozen_correct & one_nn_correct).sum())
    discordant = frozen_only_correct + one_nn_only_correct
    pvalue = float(binomtest(frozen_only_correct, n=discordant, p=0.5, alternative="two-sided").pvalue) if discordant else 1.0
    output.append(
        {
            "comparison": "frozen_ABCB1_v0_1_minus_one_nearest_neighbor",
            "metric": "mcnemar_exact_discordant_accuracy",
            "frozen_value": frozen_only_correct,
            "one_nn_value": one_nn_only_correct,
            "paired_difference": frozen_only_correct - one_nn_only_correct,
            "ci95_lower": "",
            "ci95_upper": "",
            "n": len(y),
            "bootstrap_iterations": 0,
            "direction_note": f"two-sided_exact_p={pvalue:.8f}; counts are frozen-only-correct versus 1NN-only-correct",
            "claim_boundary": "paired historical external diagnostic; no model selection or tuning on these results",
        }
    )
    return output


def metric_rows(
    rows: list[dict],
    benchmark_name: str,
    subset_name: str,
    probability_fields: dict[str, str],
    rng: np.random.Generator,
    add_bootstrap: bool = True,
) -> list[dict]:
    eligible = [row for row in rows if row.get("directional_label") in {"0", "1"}]
    y = np.asarray([int(row["directional_label"]) for row in eligible], dtype=int)
    output: list[dict] = []
    if not len(y):
        return output
    for model_name, field in probability_fields.items():
        probability = np.asarray([float(row[field]) for row in eligible], dtype=float)
        point = metric_dict(y, probability)
        intervals = bootstrap_ci(y, probability, rng) if add_bootstrap and len(set(y.tolist())) == 2 else {}
        for metric_name in ["mcc", "balanced_accuracy", "sensitivity", "specificity", "roc_auc", "average_precision", "brier"]:
            lower, upper = intervals.get(metric_name, (math.nan, math.nan))
            output.append(
                {
                    "benchmark": benchmark_name,
                    "subset": subset_name,
                    "model": model_name,
                    "probability_field": field,
                    "metric": metric_name,
                    "value": finite_or_blank(point[metric_name]),
                    "ci95_lower": finite_or_blank(lower),
                    "ci95_upper": finite_or_blank(upper),
                    "n": point["n"],
                    "n_negative": point["n_negative"],
                    "n_positive": point["n_positive"],
                    "tn": point["tn"],
                    "fp": point["fp"],
                    "fn": point["fn"],
                    "tp": point["tp"],
                    "threshold": 0.5,
                    "bootstrap_iterations": BOOTSTRAP_ITERATIONS if intervals else 0,
                    "claim_boundary": "source-call directional diagnostic; not universal binary substrate truth",
                }
            )
    return output


def build_direct_benchmark(package: dict, training_labels: list[int], training_fps: list) -> tuple[list[dict], list[dict], list[dict], list[dict], list[dict]]:
    overlay = read_csv(PROCESSED / "human_abcb1_external_entity_overlay_v1.csv")
    smiles_by_inchikey, numeric_by_entity_source = observation_index()
    source_urls = source_url_registry()
    training_keys = set(package["training_structure_keys"])
    training_connectivity = {key.split("-")[0] for key in training_keys if key}
    prior = package["training_n_positive"] / package["training_n"]
    rows: list[dict] = []

    for row in overlay:
        inchikey = row.get("standardized_inchikey", "")
        connectivity = row.get("inchikey_connectivity_block", "") or (inchikey.split("-")[0] if inchikey else "")
        smiles = smiles_by_inchikey.get(inchikey, "")
        source_ids = [value.strip() for value in row.get("source_ids", "").split("|") if value.strip()]
        source_id = source_ids[0] if len(source_ids) == 1 else " | ".join(source_ids)
        direct = int(row.get("human_direct_rows", "0") or 0) > 0
        if inchikey in training_keys:
            source_overlap = "exact_training_inchikey"
        elif connectivity in training_connectivity:
            source_overlap = "training_connectivity_block"
        else:
            source_overlap = "no_exact_or_connectivity_training_overlap"
        source_call_class, directional_label = normalize_source_call(row.get("source_assigned_calls", ""))
        base = {
            "benchmark_entity_id": f"ABCB1_EXT_{len(rows)+1:03d}",
            "entity_name": row.get("entity_names", ""),
            "standardized_inchikey": inchikey,
            "inchikey_connectivity_block": connectivity,
            "standardized_isomeric_smiles": smiles,
            "source_id": source_id,
            "source_url": " | ".join(source_urls.get(value, "") for value in source_ids if source_urls.get(value, "")),
            "source_assigned_call": row.get("source_assigned_calls", ""),
            "source_call_class": source_call_class,
            "directional_label": directional_label,
            "human_direct_rows": row.get("human_direct_rows", ""),
            "evidence_class": row.get("human_external_comparability_class", ""),
            "source_training_overlap_status": source_overlap,
            "training_overlap_status": source_overlap,
            "model_evaluation_eligible": "pending_model_standardization" if direct and smiles else "no",
            "exclusion_reason": "" if direct and smiles else ("not_human_direct" if not direct else "structure_missing"),
        }
        if smiles:
            prediction = predict_smiles(smiles, package)
            mol = Chem.MolFromSmiles(prediction.get("standardized_isomeric_smiles", ""))
            if mol is None:
                raise RuntimeError(f"prediction standardized molecule unresolved for {inchikey}")
            baselines = similarity_baselines(mol, training_labels, training_fps)
            model_inchikey = prediction.get("standardized_inchikey", "")
            model_connectivity = model_inchikey.split("-")[0] if model_inchikey else ""
            if model_inchikey in training_keys:
                model_overlap = "exact_training_inchikey_after_model_standardization"
            elif model_connectivity in training_connectivity:
                model_overlap = "training_connectivity_block_after_model_standardization"
            else:
                model_overlap = "no_exact_or_connectivity_training_overlap"
            model_eligible = direct and model_overlap == "no_exact_or_connectivity_training_overlap"
            base.update(
                {
                    "model_standardized_inchikey": model_inchikey,
                    "model_standardized_connectivity_block": model_connectivity,
                    "training_overlap_status": model_overlap,
                    "model_evaluation_eligible": "yes" if model_eligible else "no",
                    "exclusion_reason": "" if model_eligible else ("not_human_direct" if not direct else model_overlap),
                    "frozen_probability": prediction["calibrated_probability"],
                    "frozen_binary_prediction": prediction["binary_prediction"],
                    "max_training_tanimoto": prediction["max_training_morgan_tanimoto"],
                    "applicability_band": prediction["applicability_band"],
                    "screening_call": prediction["screening_call"],
                    "prior_probability": prior,
                    **baselines,
                }
            )
        else:
            base.update(
                {
                    "model_standardized_inchikey": "",
                    "model_standardized_connectivity_block": "",
                    "frozen_probability": "",
                    "frozen_binary_prediction": "",
                    "max_training_tanimoto": "",
                    "applicability_band": "unresolved_structure",
                    "screening_call": "unresolved_structure",
                    "prior_probability": prior,
                    "nearest_training_index": "",
                    "nearest_training_tanimoto": "",
                    "one_nn_probability": "",
                    "knn5_weighted_probability": "",
                }
            )

        numeric_values: list[float] = []
        numeric_fields: list[str] = []
        observation_ids: list[str] = []
        for sid in source_ids:
            for value, field, observation_id in numeric_by_entity_source.get((inchikey, sid), []):
                numeric_values.append(value)
                numeric_fields.append(field)
                observation_ids.append(observation_id)
        base.update(
            {
                "continuous_endpoint_median": round(float(np.median(numeric_values)), 8) if numeric_values else "",
                "continuous_endpoint_min": round(float(min(numeric_values)), 8) if numeric_values else "",
                "continuous_endpoint_max": round(float(max(numeric_values)), 8) if numeric_values else "",
                "continuous_endpoint_value_count": len(numeric_values),
                "continuous_endpoint_fields": " | ".join(sorted(set(numeric_fields))),
                "continuous_observation_ids": " | ".join(sorted(set(observation_ids))),
                "continuous_cross_source_pooling_allowed": "no_within_source_only",
                "claim_boundary": "direct-source directional/continuous external diagnostic only; no frozen label transfer",
            }
        )
        rows.append(base)

    eligible = [row for row in rows if row["model_evaluation_eligible"] == "yes"]
    call_rows = [row for row in eligible if row["directional_label"] in {"0", "1"}]
    rng = np.random.default_rng(SEED)
    probability_fields = {
        "frozen_ABCB1_v0_1": "frozen_probability",
        "training_class_prior": "prior_probability",
        "one_nearest_neighbor": "one_nn_probability",
        "five_nearest_neighbor_weighted": "knn5_weighted_probability",
    }
    metrics = metric_rows(call_rows, "direct_human_ABCB1_source_calls", "all_nonoverlap_explicit_calls", probability_fields, rng)
    for band in ["within_empirical_support", "low_empirical_support", "outside_empirical_support"]:
        subset = [row for row in call_rows if row["applicability_band"] == band]
        if subset:
            metrics.extend(metric_rows(subset, "direct_human_ABCB1_source_calls", band, probability_fields, rng, add_bootstrap=False))
    nonabstain = [row for row in call_rows if row["applicability_band"] != "outside_empirical_support"]
    if nonabstain:
        metrics.extend(metric_rows(nonabstain, "direct_human_ABCB1_source_calls", "within_or_caution_nonabstain", probability_fields, rng, add_bootstrap=False))
    for cutoff in [0.8, 0.5]:
        subset = [row for row in call_rows if float(row["max_training_tanimoto"]) < cutoff]
        if subset:
            metrics.extend(
                metric_rows(
                    subset,
                    "direct_human_ABCB1_source_calls",
                    f"max_training_tanimoto_below_{str(cutoff).replace('.', '_')}",
                    probability_fields,
                    rng,
                    add_bootstrap=(cutoff == 0.8),
                )
            )
    near_duplicate = [row for row in call_rows if float(row["max_training_tanimoto"]) >= 0.8]
    if near_duplicate:
        metrics.extend(
            metric_rows(
                near_duplicate,
                "direct_human_ABCB1_source_calls",
                "max_training_tanimoto_at_least_0_8",
                probability_fields,
                rng,
                add_bootstrap=False,
            )
        )

    source_sensitivity: list[dict] = []
    sources = sorted({row["source_id"] for row in call_rows})
    for source in sources:
        source_only = [row for row in call_rows if row["source_id"] == source]
        leave_out = [row for row in call_rows if row["source_id"] != source]
        for mode, subset in [("source_only", source_only), ("leave_one_source_out", leave_out)]:
            if not subset:
                continue
            y = np.asarray([int(row["directional_label"]) for row in subset], dtype=int)
            probability = np.asarray([float(row["frozen_probability"]) for row in subset], dtype=float)
            metric = metric_dict(y, probability)
            source_sensitivity.append(
                {
                    "mode": mode,
                    "source_id": source,
                    **{key: finite_or_blank(value) if isinstance(value, float) else value for key, value in metric.items()},
                    "claim_boundary": "source sensitivity; single-class source metrics are descriptive only",
                }
            )

    continuous_rows: list[dict] = []
    pooled_rank_records = []
    for source in sorted({row["source_id"] for row in eligible}):
        subset = [row for row in eligible if row["source_id"] == source and row["continuous_endpoint_median"] != ""]
        if len(subset) < 3 or len({float(row["continuous_endpoint_median"]) for row in subset}) < 3:
            continuous_rows.append(
                {
                    "source_id": source,
                    "n_entities": len(subset),
                    "endpoint_fields": " | ".join(sorted({row["continuous_endpoint_fields"] for row in subset})),
                    "model": "frozen_ABCB1_v0_1",
                    "spearman_rho": "",
                    "p_value": "",
                    "metric_status": "not_estimable_fewer_than_3_entities_or_values",
                    "claim_boundary": "within-source continuous endpoint only",
                }
            )
            continue
        endpoint = np.asarray([float(row["continuous_endpoint_median"]) for row in subset], dtype=float)
        endpoint_rank = rankdata(endpoint, method="average") / len(endpoint)
        for row, rank in zip(subset, endpoint_rank):
            pooled_rank_records.append((row, float(rank)))
        for model_name, field in probability_fields.items():
            probability = np.asarray([float(row[field]) for row in subset], dtype=float)
            rho, pvalue = spearmanr(probability, endpoint)
            continuous_rows.append(
                {
                    "source_id": source,
                    "n_entities": len(subset),
                    "endpoint_fields": " | ".join(sorted({row["continuous_endpoint_fields"] for row in subset})),
                    "model": model_name,
                    "spearman_rho": finite_or_blank(float(rho)),
                    "p_value": finite_or_blank(float(pvalue)),
                    "metric_status": "within_source_spearman",
                    "claim_boundary": "within-source continuous endpoint only; raw values are not pooled across sources",
                }
            )
    if pooled_rank_records:
        for model_name, field in probability_fields.items():
            probability = np.asarray([float(row[field]) for row, _ in pooled_rank_records], dtype=float)
            ranks = np.asarray([rank for _, rank in pooled_rank_records], dtype=float)
            rho, pvalue = spearmanr(probability, ranks)
            continuous_rows.append(
                {
                    "source_id": "POOLED_WITHIN_SOURCE_RANKS",
                    "n_entities": len(pooled_rank_records),
                    "endpoint_fields": "source-specific ranks only",
                    "model": model_name,
                    "spearman_rho": finite_or_blank(float(rho)),
                    "p_value": finite_or_blank(float(pvalue)),
                    "metric_status": "pooled_within_source_rank_sensitivity",
                    "claim_boundary": "source-normalized rank sensitivity; not pooled raw endpoint magnitude",
                }
            )

    model_comparison = paired_model_comparison(call_rows)
    return rows, metrics, source_sensitivity, continuous_rows, model_comparison


def build_legacy_benchmark(package: dict, training_labels: list[int], training_fps: list) -> tuple[list[dict], list[dict]]:
    rows = read_csv(PROCESSED / "public_independent_strict_binary_candidate_view_v1.csv")
    prior = package["training_n_positive"] / package["training_n"]
    output: list[dict] = []
    seen: set[str] = set()
    for row in rows:
        label = safe_float(row.get("binary_label"))
        if row.get("transporter_gene") != "ABCB1" or label not in {0.0, 1.0}:
            continue
        if row.get("overlap_type") != "no_daood_overlap" or row.get("source_structure_status") == "failed":
            continue
        key = row.get("standardized_inchikey", "") or row.get("standardized_smiles", "")
        if not key or key in seen:
            continue
        seen.add(key)
        smiles = row.get("standardized_isomeric_smiles", "") or row.get("standardized_smiles", "") or row.get("source_smiles", "")
        prediction = predict_smiles(smiles, package)
        mol = Chem.MolFromSmiles(prediction.get("standardized_isomeric_smiles", ""))
        if mol is None:
            continue
        baselines = similarity_baselines(mol, training_labels, training_fps)
        output.append(
            {
                "independent_record_id": row.get("independent_record_id", ""),
                "compound_name": row.get("compound_name", ""),
                "directional_label": str(int(label)),
                "source_dataset_split": row.get("source_dataset_split", ""),
                "source_doi": row.get("source_doi", ""),
                "source_url": row.get("source_url", ""),
                "evidence_grade": row.get("evidence_grade", ""),
                "standardized_inchikey": row.get("standardized_inchikey", ""),
                "frozen_probability": prediction["calibrated_probability"],
                "frozen_binary_prediction": prediction["binary_prediction"],
                "max_training_tanimoto": prediction["max_training_morgan_tanimoto"],
                "applicability_band": prediction["applicability_band"],
                "screening_call": prediction["screening_call"],
                "prior_probability": prior,
                **baselines,
                "claim_boundary": "single-source curated literature benchmark; not primary direct-assay validation",
            }
        )

    probability_fields = {
        "frozen_ABCB1_v0_1": "frozen_probability",
        "training_class_prior": "prior_probability",
        "one_nearest_neighbor": "one_nn_probability",
        "five_nearest_neighbor_weighted": "knn5_weighted_probability",
    }
    rng = np.random.default_rng(SEED + 101)
    metrics = metric_rows(output, "legacy_Bikadi_2011_single_source", "all_splits", probability_fields, rng)
    for split in ["Training", "Validation", "Test"]:
        subset = [row for row in output if row["source_dataset_split"] == split]
        metrics.extend(metric_rows(subset, "legacy_Bikadi_2011_single_source", split, probability_fields, rng, add_bootstrap=False))
    for band in ["within_empirical_support", "low_empirical_support", "outside_empirical_support"]:
        subset = [row for row in output if row["applicability_band"] == band]
        if subset:
            metrics.extend(metric_rows(subset, "legacy_Bikadi_2011_single_source", band, probability_fields, rng, add_bootstrap=False))
    return output, metrics


def metric_lookup(metrics: list[dict], benchmark: str, subset: str, model: str, metric: str) -> dict:
    return next(
        row
        for row in metrics
        if row["benchmark"] == benchmark and row["subset"] == subset and row["model"] == model and row["metric"] == metric
    )


def main() -> None:
    protocol = json.loads(LOCK_PATH.read_text(encoding="utf-8"))
    if sha256_file(MODEL_PATH) != EXPECTED_MODEL_SHA256:
        raise RuntimeError("frozen model hash mismatch")
    for entry in protocol["direct_human_abcb1_sources"]:
        path = PROJECT / entry["path"]
        if sha256_file(path) != entry["sha256"]:
            raise RuntimeError(f"source hash mismatch: {path}")
    legacy = protocol["legacy_binary_source"]
    if sha256_file(PROJECT / legacy["path"]) != legacy["sha256"]:
        raise RuntimeError("legacy source hash mismatch")

    package = joblib.load(MODEL_PATH)
    training_labels, training_fps = training_baseline(package)
    direct_rows, direct_metrics, source_sensitivity, continuous_metrics, model_comparison = build_direct_benchmark(package, training_labels, training_fps)
    legacy_rows, legacy_metrics = build_legacy_benchmark(package, training_labels, training_fps)
    all_metrics = direct_metrics + legacy_metrics

    eligible_direct = [row for row in direct_rows if row["model_evaluation_eligible"] == "yes"]
    call_direct = [row for row in eligible_direct if row["directional_label"] in {"0", "1"}]
    direct_mcc = metric_lookup(all_metrics, "direct_human_ABCB1_source_calls", "all_nonoverlap_explicit_calls", "frozen_ABCB1_v0_1", "mcc")
    direct_bacc = metric_lookup(all_metrics, "direct_human_ABCB1_source_calls", "all_nonoverlap_explicit_calls", "frozen_ABCB1_v0_1", "balanced_accuracy")
    direct_auc = metric_lookup(all_metrics, "direct_human_ABCB1_source_calls", "all_nonoverlap_explicit_calls", "frozen_ABCB1_v0_1", "roc_auc")
    estimable_sources = [
        row
        for row in continuous_metrics
        if row["model"] == "frozen_ABCB1_v0_1" and row["metric_status"] == "within_source_spearman" and row["source_id"] != "POOLED_WITHIN_SOURCE_RANKS"
    ]
    positive_continuous_sources = sum(float(row["spearman_rho"]) > 0 for row in estimable_sources if row["spearman_rho"] != "")
    loso = [row for row in source_sensitivity if row["mode"] == "leave_one_source_out" and row["balanced_accuracy"] != ""]
    min_loso_bacc = min(float(row["balanced_accuracy"]) for row in loso) if loso else math.nan
    within_n = sum(row["applicability_band"] == "within_empirical_support" for row in call_direct)
    nonabstain_n = sum(row["applicability_band"] != "outside_empirical_support" for row in call_direct)

    gates = [
        {
            "gate": "direct_MCC_lower_95CI_gt_0",
            "value": direct_mcc["ci95_lower"],
            "rule": "> 0",
            "status": "pass" if direct_mcc["ci95_lower"] != "" and float(direct_mcc["ci95_lower"]) > 0 else "fail",
            "interpretation": "directional source-call discrimination beyond chance correlation",
        },
        {
            "gate": "direct_BAcc_lower_95CI_gt_0_5",
            "value": direct_bacc["ci95_lower"],
            "rule": "> 0.5",
            "status": "pass" if direct_bacc["ci95_lower"] != "" and float(direct_bacc["ci95_lower"]) > 0.5 else "fail",
            "interpretation": "balanced accuracy lower bound exceeds chance",
        },
        {
            "gate": "positive_continuous_direction_in_two_sources",
            "value": positive_continuous_sources,
            "rule": ">= 2 independent sources with n>=3",
            "status": "pass" if positive_continuous_sources >= 2 else "fail",
            "interpretation": "continuous endpoint direction is not carried by one study",
        },
        {
            "gate": "leave_one_source_out_min_BAcc_gt_0_5",
            "value": finite_or_blank(min_loso_bacc),
            "rule": "> 0.5",
            "status": "pass" if math.isfinite(min_loso_bacc) and min_loso_bacc > 0.5 else "fail",
            "interpretation": "no single source determines the result",
        },
        {
            "gate": "nonabstain_coverage_at_least_50pct",
            "value": round(nonabstain_n / len(call_direct), 8) if call_direct else "",
            "rule": ">= 0.5",
            "status": "pass" if call_direct and nonabstain_n / len(call_direct) >= 0.5 else "fail",
            "interpretation": "applicability rejection does not hide most explicit calls",
        },
    ]

    summary = {
        "run_date": RUN_DATE,
        "protocol_lock_path": str(LOCK_PATH),
        "protocol_lock_sha256": sha256_file(LOCK_PATH),
        "benchmark_status": protocol["benchmark_status"],
        "frozen_model_sha256": sha256_file(MODEL_PATH),
        "direct_overlay_entities": len(direct_rows),
        "direct_evaluation_eligible_nonoverlap": len(eligible_direct),
        "direct_explicit_call_entities": len(call_direct),
        "direct_explicit_call_counts": dict(Counter(row["source_call_class"] for row in call_direct)),
        "direct_source_count": len({row["source_id"] for row in eligible_direct}),
        "direct_applicability_counts": dict(Counter(row["applicability_band"] for row in call_direct)),
        "direct_within_support_n": within_n,
        "direct_nonabstain_n": nonabstain_n,
        "direct_v0_1_MCC": direct_mcc,
        "direct_v0_1_balanced_accuracy": direct_bacc,
        "direct_v0_1_ROC_AUC": direct_auc,
        "continuous_estimable_sources": len(estimable_sources),
        "continuous_positive_direction_sources": positive_continuous_sources,
        "leave_one_source_out_min_balanced_accuracy": finite_or_blank(min_loso_bacc),
        "paired_frozen_vs_one_nn": model_comparison,
        "legacy_bikadi_rows": len(legacy_rows),
        "legacy_bikadi_label_counts": dict(Counter(row["directional_label"] for row in legacy_rows)),
        "q1_predictor_claim_gates": gates,
        "q1_predictor_claim_status": "pass" if all(row["status"] == "pass" for row in gates) else "fail",
        "action": "do_not_tune_on_external_results; retain bounded screener and reserve later untouched temporal/prospective set",
        "claim_boundary": "historical development external benchmark; not pristine prospective validation and no source-call label transfer to training",
    }

    write_csv(PROCESSED / "abcb1_q1_direct_external_entity_predictions_v1.csv", direct_rows)
    write_csv(PROCESSED / "abcb1_q1_external_metric_summary_v1.csv", all_metrics)
    write_csv(PROCESSED / "abcb1_q1_external_source_sensitivity_v1.csv", source_sensitivity)
    write_csv(PROCESSED / "abcb1_q1_external_continuous_metrics_v1.csv", continuous_metrics)
    write_csv(PROCESSED / "abcb1_q1_model_pairwise_comparison_v1.csv", model_comparison)
    write_csv(PROCESSED / "abcb1_q1_legacy_bikadi_predictions_v1.csv", legacy_rows)
    write_csv(PROCESSED / "abcb1_q1_predictor_claim_gates_v1.csv", gates)
    (PROCESSED / "abcb1_q1_external_benchmark_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
