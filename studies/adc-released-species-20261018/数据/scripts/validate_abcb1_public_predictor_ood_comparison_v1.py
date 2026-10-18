#!/usr/bin/env python3
"""Validate the frozen public-predictor/OOD projection outputs."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "数据" / "processed"
PROTOCOL = PROCESSED / "abcb1_public_predictor_ood_comparison_protocol_v1.json"
ANALYSIS_ROWS = PROCESSED / "abcb1_public_predictor_ood_analysis_rows_v1.csv"
MODEL_STATUS = PROCESSED / "abcb1_public_predictor_ood_model_status_v1.csv"
METRICS = PROCESSED / "abcb1_public_predictor_ood_stratum_metrics_v1.csv"
CONTRASTS = PROCESSED / "abcb1_public_predictor_ood_contrasts_v1.csv"
SUMMARY = PROCESSED / "abcb1_public_predictor_ood_comparison_summary_v1.json"
PRIOR_METRICS = PROCESSED / "abcb1_public_predictor_metrics_v1.csv"
VALIDATION = PROCESSED / "abcb1_public_predictor_ood_comparison_validation_v1.json"

EXPECTED_MODELS = {
    "astrazeneca_gnn_mtl_mdck_er_2025": math.log10(2),
    "astrazeneca_gnn_mtl_nih_mdck_er_2025": math.log10(2),
    "deeppk_pgp_substrate_2026_08_05": 0.5,
    "admet_ai_v2_pgp_broccatelli_inhibition": 0.5,
}
EXPECTED_N = 2624
NATIVE_RATE = 453 / 2624


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text())
    summary = json.loads(SUMMARY.read_text())
    analysis = pd.read_csv(ANALYSIS_ROWS)
    status = pd.read_csv(MODEL_STATUS)
    metrics = pd.read_csv(METRICS)
    contrasts = pd.read_csv(CONTRASTS)
    prior = pd.read_csv(PRIOR_METRICS)
    checks: list[dict] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"check": name, "passed": bool(condition), "detail": detail})
        if not condition:
            raise AssertionError(f"{name}: {detail}")

    for relative, expected in protocol["inputs_sha256"].items():
        observed = sha256(ROOT / relative)
        check(
            f"input_hash::{relative}",
            observed == expected,
            f"expected={expected}; observed={observed}",
        )

    check("analysis_row_count", len(analysis) == EXPECTED_N * 4, f"n={len(analysis)}")
    check(
        "analysis_key_unique",
        not analysis.duplicated(["comparator_id", "external_row_id"]).any(),
        "comparator_id/external_row_id is unique",
    )
    check(
        "model_set",
        set(analysis["comparator_id"]) == set(EXPECTED_MODELS),
        str(sorted(analysis["comparator_id"].unique())),
    )
    check(
        "structure_identity_complete",
        analysis["standardized_inchikey_public"].eq(
            analysis["standardized_inchikey_ood"]
        ).all(),
        "all InChIKey identities agree after external_row_id merge",
    )
    check(
        "reference_scope_boundary",
        analysis["public_model_training_space_ood_status"].eq(
            "not_established_by_this_analysis"
        ).all(),
        "public-model training-space OOD is never asserted",
    )
    for comparator_id, threshold in EXPECTED_MODELS.items():
        values = analysis.loc[
            analysis["comparator_id"].eq(comparator_id), "binary_threshold"
        ].unique()
        check(
            f"threshold::{comparator_id}",
            len(values) == 1 and np.isclose(values[0], threshold, rtol=0.0, atol=1e-12),
            f"observed={values.tolist()}; expected={threshold}",
        )

    check("model_status_rows", len(status) == 5, f"n={len(status)}")
    admetlab = status.loc[status["comparator_id"].eq("admetlab_3_pgp_substrate")]
    check(
        "admetlab_missing_not_imputed",
        len(admetlab) == 1
        and not bool(admetlab.iloc[0]["included_in_metrics"])
        and admetlab.iloc[0]["model_output_status"] == "unavailable_no_predictions",
        admetlab.to_dict(orient="records").__str__(),
    )
    admet_ai = status.loc[
        status["comparator_id"].eq("admet_ai_v2_pgp_broccatelli_inhibition")
    ]
    check(
        "admet_ai_endpoint_mismatch",
        len(admet_ai) == 1 and not bool(admet_ai.iloc[0]["substrate_ranking_eligible"]),
        admet_ai.to_dict(orient="records").__str__(),
    )

    check("metric_row_count", len(metrics) == 96, f"n={len(metrics)}")
    per_model = metrics.groupby("comparator_id").size().to_dict()
    check(
        "metric_rows_per_model",
        all(per_model.get(model) == 24 for model in EXPECTED_MODELS),
        str(per_model),
    )
    check(
        "no_outcome_based_selection",
        not metrics["selection_uses_outcomes"].astype(bool).any(),
        "all selection_uses_outcomes values are false",
    )

    full = metrics.loc[metrics["evaluation_layer"].eq("full_coverage")]
    check(
        "full_coverage_counts",
        full["selected_n"].eq(full["stratum_n"]).all()
        and np.allclose(full["coverage_within_stratum"], 1.0),
        "full rows retain each complete stratum",
    )
    native = metrics.loc[
        metrics["evaluation_layer"].eq("model_native_margin_top17_3pct_within_stratum")
    ].copy()
    expected_native_n = native["stratum_n"].map(
        lambda value: max(1, int(math.floor(value * NATIVE_RATE + 0.5)))
    )
    check(
        "native_coverage_counts",
        native["selected_n"].eq(expected_native_n).all() and native["n"].eq(expected_native_n).all(),
        "native selection uses the frozen 453/2624 rate with deterministic rounding",
    )

    unique_external = analysis.drop_duplicates("external_row_id")
    expected_counts = {
        ("external_d5_quintile", "Q1_nearest"): 525,
        ("external_d5_quintile", "Q2"): 525,
        ("external_d5_quintile", "Q3"): 524,
        ("external_d5_quintile", "Q4"): 525,
        ("external_d5_quintile", "Q5_farthest"): 525,
        ("dual_ood", "yes"): 609,
        ("dual_ood", "no"): 2015,
        ("kmeans_ood20", "yes"): 544,
        ("kmeans_ood20", "no"): 2080,
        ("strict_lohi_ood", "yes"): 2581,
        ("strict_lohi_ood", "no"): 43,
    }
    for (column, value), expected in expected_counts.items():
        observed = int(unique_external[column].eq(value).sum())
        check(f"stratum_count::{column}::{value}", observed == expected, f"n={observed}")

    prior_full = prior.loc[prior["population"].eq("all_identity_clean")].set_index(
        "comparator_id"
    )
    current_full = metrics.loc[
        metrics["evaluation_layer"].eq("full_coverage")
        & metrics["stratum_type"].eq("overall")
    ].set_index("comparator_id")
    for comparator_id in EXPECTED_MODELS:
        for metric in ["balanced_accuracy", "mcc", "roc_auc", "spearman_rho"]:
            observed = float(current_full.loc[comparator_id, metric])
            expected = float(prior_full.loc[comparator_id, metric])
            check(
                f"prior_anchor::{comparator_id}::{metric}",
                np.isclose(observed, expected, rtol=0.0, atol=1e-12),
                f"current={observed}; prior={expected}",
            )

    check("contrast_row_count", len(contrasts) == 32, f"n={len(contrasts)}")
    for relative, expected in summary["outputs"].items():
        observed = sha256(ROOT / relative)
        check(
            f"output_hash::{relative}",
            observed == expected,
            f"expected={expected}; observed={observed}",
        )

    validation = {
        "validation_id": "ABCB1_PUBLIC_PREDICTOR_OOD_COMPARISON_VALIDATION_V1",
        "status": "passed",
        "n_checks": len(checks),
        "all_passed": all(item["passed"] for item in checks),
        "checks": checks,
        "validated_outputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [ANALYSIS_ROWS, MODEL_STATUS, METRICS, CONTRASTS, SUMMARY]
        },
        "boundary": (
            "Validation confirms computational consistency only. It does not convert the "
            "retrospective project-reference distance audit into prospective or public-model-training-space OOD validation."
        ),
    }
    VALIDATION.write_text(json.dumps(validation, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "passed", "n_checks": len(checks)}, indent=2))


if __name__ == "__main__":
    main()
