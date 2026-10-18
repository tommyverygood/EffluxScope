#!/usr/bin/env python3
"""Validate the public-predictor rejection benchmark artifacts and anchors."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "数据" / "processed"
PROTOCOL = PROCESSED / "abcb1_public_predictor_benchmark_protocol_v1.json"
INPUT = PROCESSED / "abcb1_public_predictor_benchmark_input_v1.csv"
SUMMARY = PROCESSED / "abcb1_public_predictor_benchmark_summary_v1.json"
PREDICTIONS = PROCESSED / "abcb1_public_predictor_predictions_long_v1.csv"
METRICS = PROCESSED / "abcb1_public_predictor_metrics_v1.csv"
COVERAGE = PROCESSED / "abcb1_public_predictor_coverage_curve_v1.csv"
BOOTSTRAP = PROCESSED / "abcb1_public_predictor_cluster_bootstrap_v1.csv"
CONTRASTS = PROCESSED / "abcb1_public_predictor_rejection_contrasts_v1.csv"
OVERLAP = PROCESSED / "abcb1_public_predictor_deeppk_overlap_audit_v1.csv"
FIGURE_MANIFEST = PROCESSED / "abcb1_public_predictor_figure_manifest_v1.csv"
REPORT = ROOT / "ABCB1公开预测器拒判比较报告_v1_20260805.md"
VALIDATION = PROCESSED / "abcb1_public_predictor_benchmark_validation_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def close(value: float, expected: float, tolerance: float = 5e-4) -> bool:
    return math.isclose(float(value), expected, rel_tol=0, abs_tol=tolerance)


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text())
    summary = json.loads(SUMMARY.read_text())
    input_data = pd.read_csv(INPUT)
    predictions = pd.read_csv(PREDICTIONS)
    metrics = pd.read_csv(METRICS)
    coverage = pd.read_csv(COVERAGE)
    bootstrap = pd.read_csv(BOOTSTRAP)
    contrasts = pd.read_csv(CONTRASTS)
    overlap = pd.read_csv(OVERLAP)
    figure_manifest = pd.read_csv(FIGURE_MANIFEST)
    checks: list[dict] = []

    def check(name: str, condition: bool, evidence) -> None:
        checks.append({"check": name, "pass": bool(condition), "evidence": evidence})

    check("protocol_truthful_posthoc", "preliminary_astrazeneca" in protocol["timing"], protocol["timing"])
    check("protocol_input_hash", protocol["benchmark"]["input_sha256"] == sha256(INPUT), sha256(INPUT))
    check("input_n", len(input_data) == 2624, len(input_data))
    check("input_unique_ids", input_data["external_row_id"].nunique() == 2624, input_data["external_row_id"].nunique())
    check("input_accepted_n", int(input_data["v0_4_accepted"].sum()) == 453, int(input_data["v0_4_accepted"].sum()))

    check("summary_protocol_hash", summary["protocol_sha256"] == sha256(PROTOCOL), sha256(PROTOCOL))
    check("summary_status", summary["status"].startswith("completed_available"), summary["status"])
    check("admetlab_pending", summary["admetlab_status"]["status"] == "pending_official_server_failure", summary["admetlab_status"])
    check("admetlab_error_preserved", "BSEP" in summary["admetlab_status"]["observed_error"], summary["admetlab_status"]["observed_error"])

    expected_models = {
        "internal_v0_4_primary_score",
        "astrazeneca_gnn_mtl_mdck_er_2025",
        "astrazeneca_gnn_mtl_nih_mdck_er_2025",
        "deeppk_pgp_substrate_2026_08_05",
        "admet_ai_v2_pgp_broccatelli_inhibition",
    }
    check("prediction_model_set", set(predictions["comparator_id"]) == expected_models, sorted(predictions["comparator_id"].unique()))
    check("prediction_row_count", len(predictions) == 2624 * 5, len(predictions))
    check("prediction_each_model_n", predictions.groupby("comparator_id").size().eq(2624).all(), predictions.groupby("comparator_id").size().to_dict())
    check("prediction_scores_complete", predictions["score"].notna().all(), int(predictions["score"].isna().sum()))
    check("prediction_labels_binary", set(predictions["predicted_label"]) <= {0, 1}, sorted(predictions["predicted_label"].unique()))
    check("admet_ai_endpoint_mismatch", not predictions.loc[predictions["comparator_id"].str.startswith("admet_ai"), "eligible_primary_comparison"].any(), "ineligible")

    check("deeppk_overlap_rows", len(overlap) == 2624, len(overlap))
    check("deeppk_train_overlap_3", int(overlap["deeppk_exact_train_overlap"].sum()) == 3, int(overlap["deeppk_exact_train_overlap"].sum()))
    check("deeppk_validation_overlap_0", int(overlap["deeppk_exact_validation_overlap"].sum()) == 0, int(overlap["deeppk_exact_validation_overlap"].sum()))
    check("deeppk_test_overlap_0", int(overlap["deeppk_exact_test_overlap"].sum()) == 0, int(overlap["deeppk_exact_test_overlap"].sum()))
    check("deeppk_development_clean_n", int(overlap["deeppk_development_identity_clean"].sum()) == 2621, int(overlap["deeppk_development_identity_clean"].sum()))

    def metric(model: str, population: str) -> pd.Series:
        rows = metrics.loc[metrics["comparator_id"].eq(model) & metrics["population"].eq(population)]
        if len(rows) != 1:
            raise ValueError(f"Expected one metric row for {model}/{population}")
        return rows.iloc[0]

    anchors = {
        ("internal_v0_4_primary_score", "all_identity_clean"): (0.58169541, 0.24718176),
        ("internal_v0_4_primary_score", "v0_4_accepted_mask"): (0.73510345, 0.4792),
        ("astrazeneca_gnn_mtl_mdck_er_2025", "all_identity_clean"): (0.74325, 0.63007406),
        ("astrazeneca_gnn_mtl_mdck_er_2025", "v0_4_accepted_mask"): (0.7936, 0.6796),
        ("astrazeneca_gnn_mtl_mdck_er_2025", "native_margin_matched_453"): (0.9017, 0.6841),
        ("astrazeneca_gnn_mtl_nih_mdck_er_2025", "all_identity_clean"): (0.6986, 0.63546600),
        ("deeppk_pgp_substrate_2026_08_05", "all_identity_clean"): (0.5776, 0.4253),
        ("deeppk_pgp_substrate_2026_08_05", "v0_4_accepted_mask"): (0.6275, 0.5438),
        ("deeppk_pgp_substrate_2026_08_05", "native_margin_matched_453"): (0.5291, 0.1098),
        ("admet_ai_v2_pgp_broccatelli_inhibition", "all_identity_clean"): (0.5440, 0.1064),
    }
    for (model, population), (expected_bacc, expected_rho) in anchors.items():
        row = metric(model, population)
        check(f"anchor_bacc::{model}::{population}", close(row["balanced_accuracy"], expected_bacc), row["balanced_accuracy"])
        check(f"anchor_rho::{model}::{population}", close(row["spearman_rho"], expected_rho, 8e-4), row["spearman_rho"])

    deep_all = metric("deeppk_pgp_substrate_2026_08_05", "all_identity_clean")
    check("deeppk_threshold_asymmetry", deep_all["sensitivity"] < 0.25 and deep_all["specificity"] > 0.90, {"sensitivity": deep_all["sensitivity"], "specificity": deep_all["specificity"]})
    deep_clean = metric("deeppk_pgp_substrate_2026_08_05", "deeppk_development_identity_clean")
    check("deeppk_clean_result_stable", abs(deep_clean["balanced_accuracy"] - deep_all["balanced_accuracy"]) < 0.002, [deep_clean["balanced_accuracy"], deep_all["balanced_accuracy"]])

    check("coverage_rows", len(coverage) == 5 * 7, len(coverage))
    check("coverage_model_set", set(coverage["comparator_id"]) == expected_models, sorted(coverage["comparator_id"].unique()))
    check("bootstrap_rows", len(bootstrap) == 5 * 2000 * 4, len(bootstrap))
    check("bootstrap_replicates", bootstrap.groupby(["comparator_id", "population"])["replicate"].nunique().eq(2000).all(), bootstrap.groupby(["comparator_id", "population"])["replicate"].nunique().min())
    check("bootstrap_model_set", set(bootstrap["comparator_id"]) == expected_models, sorted(bootstrap["comparator_id"].unique()))

    def contrast(model: str, name: str) -> pd.Series:
        rows = contrasts.loc[contrasts["comparator_id"].eq(model) & contrasts["contrast"].eq(name)]
        if len(rows) != 1:
            raise ValueError(f"Expected one contrast row for {model}/{name}")
        return rows.iloc[0]

    deep_gain = contrast("deeppk_pgp_substrate_2026_08_05", "v0_4_mask_minus_all")
    check("deeppk_mask_gain_positive", deep_gain["ci95_low"] > 0, deep_gain.to_dict())
    deep_vs_margin = contrast("deeppk_pgp_substrate_2026_08_05", "v0_4_mask_minus_native_margin")
    check("deeppk_mask_beats_native_margin", deep_vs_margin["ci95_low"] > 0, deep_vs_margin.to_dict())
    az_gain = contrast("astrazeneca_gnn_mtl_mdck_er_2025", "v0_4_mask_minus_all")
    check("az_mdck_mask_gain_positive", az_gain["ci95_low"] > 0, az_gain.to_dict())
    az_vs_margin = contrast("astrazeneca_gnn_mtl_mdck_er_2025", "v0_4_mask_minus_native_margin")
    check("az_native_margin_beats_mask", az_vs_margin["ci95_high"] < 0, az_vs_margin.to_dict())
    check("contrast_rows", len(contrasts) == 10, len(contrasts))

    check("figure_manifest_rows", len(figure_manifest) == 3, len(figure_manifest))
    for _, row in figure_manifest.iterrows():
        path = ROOT / row["file"]
        check(f"figure_exists::{path.name}", path.exists() and path.stat().st_size > 1000, path.stat().st_size if path.exists() else 0)
        check(f"figure_hash::{path.name}", path.exists() and sha256(path) == row["sha256"], row["sha256"])

    report_text = REPORT.read_text()
    for token in [
        "AZ MDCK ER",
        "0.902",
        "Deep-PK",
        "0.628",
        "3/2,624",
        "P-gp inhibition",
        "KeyError: ['BSEP']",
        "事后比较",
    ]:
        check(f"report_token::{token}", token in report_text, token)
    check("report_substantial", REPORT.stat().st_size > 8000, REPORT.stat().st_size)

    failed = [item for item in checks if not item["pass"]]
    result = {
        "validation_version": "ABCB1_PUBLIC_PREDICTOR_BENCHMARK_VALIDATION_V1",
        "validation_date": "2026-08-05",
        "overall_status": "pass" if not failed else "fail",
        "checks_passed": len(checks) - len(failed),
        "checks_total": len(checks),
        "checks": checks,
    }
    VALIDATION.write_text(json.dumps(result, ensure_ascii=False, indent=2, default=str) + "\n")
    print(json.dumps({k: result[k] for k in ["overall_status", "checks_passed", "checks_total"]}, ensure_ascii=False, indent=2))
    if failed:
        print(json.dumps(failed, ensure_ascii=False, indent=2, default=str))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
