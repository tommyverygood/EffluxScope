#!/usr/bin/env python3
"""Project frozen public-model outputs onto frozen project-reference OOD strata.

This is a retrospective, no-refit, no-threshold-change analysis.  The structural
strata describe distance from this project's ABCB1 reference chemistry only; they
must not be described as OOD relative to a public model's proprietary/development
training space.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import balanced_accuracy_score, matthews_corrcoef, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "数据" / "processed"

PROTOCOL = PROCESSED / "abcb1_public_predictor_ood_comparison_protocol_v1.json"
OOD_ASSIGNMENTS = PROCESSED / "abcb1_ood_structural_assignments_v1.csv"
OOD_LOCK = PROCESSED / "abcb1_ood_structure_assignment_lock_v1.json"
PUBLIC_LONG = PROCESSED / "abcb1_public_predictor_predictions_long_v1.csv"
PUBLIC_PROTOCOL = PROCESSED / "abcb1_public_predictor_benchmark_protocol_v1.json"
PUBLIC_SUMMARY = PROCESSED / "abcb1_public_predictor_benchmark_summary_v1.json"

ANALYSIS_ROWS = PROCESSED / "abcb1_public_predictor_ood_analysis_rows_v1.csv"
MODEL_STATUS = PROCESSED / "abcb1_public_predictor_ood_model_status_v1.csv"
METRICS = PROCESSED / "abcb1_public_predictor_ood_stratum_metrics_v1.csv"
CONTRASTS = PROCESSED / "abcb1_public_predictor_ood_contrasts_v1.csv"
SUMMARY = PROCESSED / "abcb1_public_predictor_ood_comparison_summary_v1.json"

EXPECTED_N = 2624
NATIVE_COVERAGE = 453 / 2624

MODEL_SPECS = {
    "astrazeneca_gnn_mtl_mdck_er_2025": {
        "display_name": "AZ GNN-MTL MDCK ER",
        "endpoint": "MDR1-MDCK LogER",
        "threshold": math.log10(2),
        "score_direction": "higher_is_more_efflux",
        "substrate_ranking_eligible": True,
        "model_output_status": "available_complete",
        "training_overlap_status": "not_auditable_proprietary_AZ_training_structures",
        "endpoint_limit": "same assay family; AZ training identities unavailable",
    },
    "astrazeneca_gnn_mtl_nih_mdck_er_2025": {
        "display_name": "AZ GNN-MTL NIH-MDCK ER",
        "endpoint": "NIH MDR1-MDCK LogER",
        "threshold": math.log10(2),
        "score_direction": "higher_is_more_efflux",
        "substrate_ranking_eligible": True,
        "model_output_status": "available_complete",
        "training_overlap_status": "not_auditable_proprietary_AZ_training_structures",
        "endpoint_limit": "transporter-aligned continuous endpoint; assay differs from Biogen panel",
    },
    "deeppk_pgp_substrate_2026_08_05": {
        "display_name": "Deep-PK P-gp substrate",
        "endpoint": "P-gp substrate probability",
        "threshold": 0.5,
        "score_direction": "higher_is_more_substrate_like",
        "substrate_ranking_eligible": True,
        "model_output_status": "available_complete",
        "training_overlap_status": "exact_development_identity_overlap_audited_previously",
        "endpoint_limit": "heterogeneous substrate label, not continuous MDR1-MDCK ER",
    },
    "admet_ai_v2_pgp_broccatelli_inhibition": {
        "display_name": "ADMET-AI P-gp inhibition",
        "endpoint": "P-gp inhibition probability",
        "threshold": 0.5,
        "score_direction": "higher_is_more_inhibitory",
        "substrate_ranking_eligible": False,
        "model_output_status": "available_endpoint_mismatch_negative_control",
        "training_overlap_status": "not_ranked_due_to_endpoint_mismatch",
        "endpoint_limit": "inhibition is not substrate transport; excluded from substrate ranking",
    },
}

INPUT_PATHS = [OOD_ASSIGNMENTS, OOD_LOCK, PUBLIC_LONG, PUBLIC_PROTOCOL, PUBLIC_SUMMARY]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        result = float(value)
        return result if math.isfinite(result) else None
    if pd.isna(value):
        return None
    return value


def verify_protocol_hashes(protocol: dict) -> None:
    expected = protocol["inputs_sha256"]
    for path in INPUT_PATHS:
        relative = str(path.relative_to(ROOT))
        observed = sha256(path)
        if expected.get(relative) != observed:
            raise ValueError(
                f"Frozen input hash mismatch for {relative}: "
                f"expected {expected.get(relative)}, observed {observed}"
            )


def model_status_frame() -> pd.DataFrame:
    rows = []
    for comparator_id, spec in MODEL_SPECS.items():
        rows.append({"comparator_id": comparator_id, "included_in_metrics": True, **spec})
    rows.append(
        {
            "comparator_id": "admetlab_3_pgp_substrate",
            "included_in_metrics": False,
            "display_name": "ADMETlab 3.0 P-gp substrate",
            "endpoint": "P-gp substrate probability",
            "threshold": math.nan,
            "score_direction": "higher_expected_to_be_more_substrate_like",
            "substrate_ranking_eligible": True,
            "model_output_status": "unavailable_no_predictions",
            "training_overlap_status": "not_evaluable",
            "endpoint_limit": (
                "official server-side failure recorded on 2026-08-05; no output was imputed "
                "and no alternate ADMETlab version was substituted"
            ),
        }
    )
    return pd.DataFrame(rows)


def load_and_merge() -> pd.DataFrame:
    assignments = pd.read_csv(OOD_ASSIGNMENTS)
    public = pd.read_csv(PUBLIC_LONG)
    public = public.loc[public["comparator_id"].isin(MODEL_SPECS)].copy()

    if len(assignments) != EXPECTED_N or not assignments["external_row_id"].is_unique:
        raise ValueError("Frozen OOD assignment table must contain 2,624 unique external_row_id rows")
    for comparator_id in MODEL_SPECS:
        group = public.loc[public["comparator_id"].eq(comparator_id)]
        if len(group) != EXPECTED_N or group["external_row_id"].nunique() != EXPECTED_N:
            raise ValueError(f"{comparator_id} does not have exactly 2,624 unique predictions")
        expected_threshold = MODEL_SPECS[comparator_id]["threshold"]
        observed_thresholds = group["binary_threshold"].dropna().unique()
        if len(observed_thresholds) != 1 or not np.isclose(
            observed_thresholds[0], expected_threshold, rtol=0.0, atol=1e-12
        ):
            raise ValueError(f"Frozen threshold mismatch for {comparator_id}")

    structural_columns = [
        "external_row_id",
        "standardized_inchikey",
        "connectivity_block",
        "joint_model_d5",
        "external_d5_quintile",
        "dual_ood",
        "kmeans_ood20",
        "strict_lohi_ood",
        "butina_cluster_id",
    ]
    merged = public.merge(
        assignments[structural_columns],
        on="external_row_id",
        how="left",
        validate="many_to_one",
        suffixes=("_public", "_ood"),
    )
    if merged["joint_model_d5"].isna().any():
        raise ValueError("At least one public-model row did not merge to a frozen OOD assignment")
    if not merged["standardized_inchikey_public"].eq(
        merged["standardized_inchikey_ood"]
    ).all():
        raise ValueError("Structure identity mismatch after external_row_id merge")

    merged["model_threshold_verified"] = True
    merged["substrate_ranking_eligible"] = merged["comparator_id"].map(
        lambda value: MODEL_SPECS[value]["substrate_ranking_eligible"]
    )
    merged["ood_reference_scope"] = "this_project_ABCB1_joint_model_reference"
    merged["public_model_training_space_ood_status"] = "not_established_by_this_analysis"
    merged = merged.sort_values(["comparator_id", "external_row_id"]).reset_index(drop=True)
    return merged


def safe_div(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator else math.nan


def calculate_metrics(frame: pd.DataFrame) -> dict:
    y = frame["observed_er2_label"].astype(int).to_numpy()
    score = frame["score"].astype(float).to_numpy()
    pred = frame["predicted_label"].astype(int).to_numpy()
    n = len(frame)
    n_positive = int((y == 1).sum())
    n_negative = int((y == 0).sum())
    tp = int(((y == 1) & (pred == 1)).sum())
    tn = int(((y == 0) & (pred == 0)).sum())
    fp = int(((y == 0) & (pred == 1)).sum())
    fn = int(((y == 1) & (pred == 0)).sum())
    both_classes = n_positive > 0 and n_negative > 0

    if n >= 2 and np.unique(score).size > 1:
        rho = spearmanr(score, frame["observed_log10_mdr1_mdck_er"].astype(float)).statistic
    else:
        rho = math.nan
    adequate = n >= 50 and n_positive >= 20 and n_negative >= 20
    return {
        "n": int(n),
        "n_positive": n_positive,
        "n_negative": n_negative,
        "positive_prevalence": safe_div(n_positive, n),
        "predicted_positive_n": int((pred == 1).sum()),
        "balanced_accuracy": (
            float(balanced_accuracy_score(y, pred)) if both_classes else math.nan
        ),
        "mcc": float(matthews_corrcoef(y, pred)) if both_classes else math.nan,
        "roc_auc": float(roc_auc_score(y, score)) if both_classes else math.nan,
        "spearman_rho": float(rho),
        "sensitivity": safe_div(tp, tp + fn),
        "specificity": safe_div(tn, tn + fp),
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "support_status": "adequate_descriptive_support" if adequate else "sparse_descriptive_only",
    }


def stratum_masks(group: pd.DataFrame) -> list[tuple[str, str, pd.Series]]:
    masks: list[tuple[str, str, pd.Series]] = [
        ("overall", "all", pd.Series(True, index=group.index))
    ]
    quintiles = ["Q1_nearest", "Q2", "Q3", "Q4", "Q5_farthest"]
    masks.extend(
        ("external_d5_quintile", value, group["external_d5_quintile"].eq(value))
        for value in quintiles
    )
    for variable in ["dual_ood", "kmeans_ood20", "strict_lohi_ood"]:
        masks.extend((variable, value, group[variable].eq(value)) for value in ["no", "yes"])
    return masks


def top_native_margin(frame: pd.DataFrame) -> pd.DataFrame:
    n_keep = max(1, int(math.floor(len(frame) * NATIVE_COVERAGE + 0.5)))
    return frame.sort_values(
        ["native_margin", "external_row_id"], ascending=[False, True]
    ).head(n_keep)


def analyze_strata(merged: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for comparator_id, group in merged.groupby("comparator_id", sort=False):
        spec = MODEL_SPECS[comparator_id]
        for stratum_type, stratum_value, mask in stratum_masks(group):
            full = group.loc[mask].copy()
            if full.empty:
                raise ValueError(f"Empty frozen stratum: {stratum_type}={stratum_value}")
            native = top_native_margin(full)
            for layer, selected in [
                ("full_coverage", full),
                ("model_native_margin_top17_3pct_within_stratum", native),
            ]:
                rows.append(
                    {
                        "comparator_id": comparator_id,
                        "display_name": spec["display_name"],
                        "endpoint": spec["endpoint"],
                        "substrate_ranking_eligible": spec["substrate_ranking_eligible"],
                        "training_overlap_status": spec["training_overlap_status"],
                        "stratum_type": stratum_type,
                        "stratum_value": stratum_value,
                        "evaluation_layer": layer,
                        "stratum_n": int(len(full)),
                        "selected_n": int(len(selected)),
                        "coverage_within_stratum": len(selected) / len(full),
                        "binary_threshold": spec["threshold"],
                        "selection_uses_outcomes": False,
                        **calculate_metrics(selected),
                    }
                )
    return pd.DataFrame(rows)


def contrast_table(metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    contrast_pairs = [
        ("external_d5_quintile", "Q5_farthest", "Q1_nearest", "Q5_minus_Q1"),
        ("dual_ood", "yes", "no", "yes_minus_no"),
        ("kmeans_ood20", "yes", "no", "yes_minus_no"),
        ("strict_lohi_ood", "yes", "no", "yes_minus_no"),
    ]
    measures = ["balanced_accuracy", "mcc", "roc_auc", "spearman_rho"]
    for comparator_id, model_frame in metrics.groupby("comparator_id", sort=False):
        for layer, layer_frame in model_frame.groupby("evaluation_layer", sort=False):
            for stratum_type, high, low, contrast_name in contrast_pairs:
                high_row = layer_frame.loc[
                    layer_frame["stratum_type"].eq(stratum_type)
                    & layer_frame["stratum_value"].eq(high)
                ]
                low_row = layer_frame.loc[
                    layer_frame["stratum_type"].eq(stratum_type)
                    & layer_frame["stratum_value"].eq(low)
                ]
                if len(high_row) != 1 or len(low_row) != 1:
                    raise ValueError(f"Contrast lookup failed for {comparator_id}/{stratum_type}")
                row = {
                    "comparator_id": comparator_id,
                    "display_name": high_row.iloc[0]["display_name"],
                    "substrate_ranking_eligible": bool(
                        high_row.iloc[0]["substrate_ranking_eligible"]
                    ),
                    "evaluation_layer": layer,
                    "stratum_type": stratum_type,
                    "contrast": contrast_name,
                    "higher_distance_or_pressure_value": high,
                    "lower_distance_or_reference_value": low,
                    "higher_n": int(high_row.iloc[0]["n"]),
                    "lower_n": int(low_row.iloc[0]["n"]),
                }
                for measure in measures:
                    row[f"{measure}_difference"] = (
                        float(high_row.iloc[0][measure]) - float(low_row.iloc[0][measure])
                    )
                rows.append(row)
    return pd.DataFrame(rows)


def write_summary(
    protocol: dict,
    merged: pd.DataFrame,
    status: pd.DataFrame,
    metrics: pd.DataFrame,
    contrasts: pd.DataFrame,
) -> dict:
    primary = metrics.loc[
        metrics["substrate_ranking_eligible"]
        & metrics["evaluation_layer"].eq("full_coverage")
    ].copy()
    native = metrics.loc[
        metrics["substrate_ranking_eligible"]
        & metrics["evaluation_layer"].eq(
            "model_native_margin_top17_3pct_within_stratum"
        )
    ].copy()
    summary = {
        "summary_id": "ABCB1_PUBLIC_PREDICTOR_OOD_COMPARISON_SUMMARY_V1",
        "status": "completed_available_models_admetlab_missing",
        "protocol": str(PROTOCOL.relative_to(ROOT)),
        "protocol_sha256": sha256(PROTOCOL),
        "n_external_rows": EXPECTED_N,
        "n_models_with_outputs": len(MODEL_SPECS),
        "n_substrate_ranking_eligible_models": int(
            status["substrate_ranking_eligible"].fillna(False).astype(bool).mul(
                status["included_in_metrics"]
            ).sum()
        ),
        "merge_validation": {
            "external_row_id_many_to_one_merge": "passed",
            "standardized_inchikey_identity_match": "passed",
            "frozen_model_thresholds": "passed",
            "outcome_based_selection": False,
        },
        "frozen_stratum_counts": {
            "external_d5_quintile": merged.drop_duplicates("external_row_id")[
                "external_d5_quintile"
            ].value_counts().sort_index().to_dict(),
            "dual_ood": merged.drop_duplicates("external_row_id")["dual_ood"].value_counts().to_dict(),
            "kmeans_ood20": merged.drop_duplicates("external_row_id")[
                "kmeans_ood20"
            ].value_counts().to_dict(),
            "strict_lohi_ood": merged.drop_duplicates("external_row_id")[
                "strict_lohi_ood"
            ].value_counts().to_dict(),
        },
        "model_availability": status.to_dict(orient="records"),
        "primary_full_coverage_metrics": primary[
            [
                "comparator_id",
                "stratum_type",
                "stratum_value",
                "n",
                "balanced_accuracy",
                "mcc",
                "roc_auc",
                "spearman_rho",
                "support_status",
            ]
        ].to_dict(orient="records"),
        "primary_native_margin_17_3pct_metrics": native[
            [
                "comparator_id",
                "stratum_type",
                "stratum_value",
                "n",
                "coverage_within_stratum",
                "balanced_accuracy",
                "mcc",
                "roc_auc",
                "spearman_rho",
                "support_status",
            ]
        ].to_dict(orient="records"),
        "claim_boundaries": protocol["prohibitions"],
        "interpretation": (
            "All OOD labels here are distances or pressure sets relative to this project's frozen "
            "ABCB1 reference chemistry. They do not establish OOD relative to AstraZeneca, Deep-PK, "
            "ADMET-AI, or ADMETlab training chemistry. ADMET-AI is endpoint-mismatched and ADMETlab "
            "has no available predictions."
        ),
        "outputs": {},
    }
    for path in [ANALYSIS_ROWS, MODEL_STATUS, METRICS, CONTRASTS]:
        summary["outputs"][str(path.relative_to(ROOT))] = sha256(path)
    SUMMARY.write_text(
        json.dumps(json_safe(summary), ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    )
    return summary


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text())
    verify_protocol_hashes(protocol)

    status = model_status_frame()
    merged = load_and_merge()
    metrics = analyze_strata(merged)
    contrasts = contrast_table(metrics)

    status.to_csv(MODEL_STATUS, index=False)
    merged.to_csv(ANALYSIS_ROWS, index=False)
    metrics.to_csv(METRICS, index=False)
    contrasts.to_csv(CONTRASTS, index=False)
    summary = write_summary(protocol, merged, status, metrics, contrasts)

    print(
        json.dumps(
            {
                "status": summary["status"],
                "n_external_rows": summary["n_external_rows"],
                "n_models_with_outputs": summary["n_models_with_outputs"],
                "outputs": summary["outputs"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
