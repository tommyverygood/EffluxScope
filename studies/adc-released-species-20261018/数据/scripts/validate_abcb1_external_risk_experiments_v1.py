#!/usr/bin/env python3
"""Validate frozen-input integrity and outputs of ABCB1 external risk experiments."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"
PROTOCOL = PROCESSED / "abcb1_external_risk_experiment_protocol_v1.json"
SUMMARY = PROCESSED / "abcb1_external_risk_experiment_summary_v1.json"
VALIDATION = PROCESSED / "abcb1_external_risk_experiment_validation_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text())
    summary = json.loads(SUMMARY.read_text())
    curve = pd.read_csv(PROCESSED / "abcb1_external_risk_coverage_curve_v1.csv")
    operating = pd.read_csv(PROCESSED / "abcb1_external_risk_operating_points_v1.csv")
    cumulative = pd.read_csv(PROCESSED / "abcb1_external_rejector_cumulative_ablation_v1.csv")
    loo = pd.read_csv(PROCESSED / "abcb1_external_rejector_leave_one_out_v1.csv")
    calibration_bins = pd.read_csv(PROCESSED / "abcb1_external_calibration_bins_v1.csv")
    calibration_summary = pd.read_csv(PROCESSED / "abcb1_external_calibration_summary_v1.csv")
    error_strata = pd.read_csv(PROCESSED / "abcb1_external_error_strata_v1.csv")
    enrichment = pd.read_csv(PROCESSED / "abcb1_external_selective_enrichment_bootstrap_v1.csv")
    manifest = pd.read_csv(PROCESSED / "abcb1_external_risk_figure_manifest_v1.csv")

    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, evidence: object) -> None:
        checks.append({"check": name, "pass": bool(condition), "evidence": evidence})

    check("protocol_status_truthful_posthoc", protocol["status"].startswith("locked_posthoc_secondary_analysis"), protocol["status"])
    check("summary_completed", summary["status"] == "completed_frozen_model_posthoc_secondary_external_risk_experiments", summary["status"])
    check("protocol_hash_recorded", summary["protocol_sha256"] == sha256(PROTOCOL), summary["protocol_sha256"])
    for label, entry in protocol["frozen_inputs"].items():
        path = PROJECT / entry["path"]
        check(f"frozen_input_exists_{label}", path.exists(), entry["path"])
        check(f"frozen_input_hash_{label}", path.exists() and sha256(path) == entry["sha256"], entry["sha256"])

    for relative, expected_hash in summary["output_hashes"].items():
        path = PROJECT / relative
        check(f"output_exists_{path.name}", path.exists(), relative)
        check(f"output_hash_{path.name}", path.exists() and sha256(path) == expected_hash, expected_hash)

    check("identity_clean_n", summary["analysis_populations"]["identity_clean_n"] == 2624, summary["analysis_populations"])
    check("identity_clean_class_counts", summary["analysis_populations"]["identity_clean_positive_n"] == 1145 and summary["analysis_populations"]["identity_clean_negative_n"] == 1479, summary["analysis_populations"])
    check("scaffold_clean_n", summary["analysis_populations"]["scaffold_clean_n"] == 2562, summary["analysis_populations"])
    reproduction = summary["frozen_rule_reconstruction"]
    check("full_rule_reconstruction_exact", reproduction["exactly_reproduced"] and reproduction["symmetric_difference_n"] == 0, reproduction)
    check("full_rule_accepted_n", reproduction["frozen_accepted_n"] == 453 and reproduction["reconstructed_accepted_n"] == 453, reproduction)

    identity_points = operating[operating["population"].eq("identity_clean")].set_index("label")
    required_points = {"forced_binary_all", "frozen_probability_margin_only", "frozen_full_v0_4", "forced_binary_on_v0_4_rejected_layer"}
    check("all_operating_points_present", required_points.issubset(identity_points.index), sorted(identity_points.index))
    forced = identity_points.loc["forced_binary_all"]
    selective = identity_points.loc["frozen_full_v0_4"]
    rejected = identity_points.loc["forced_binary_on_v0_4_rejected_layer"]
    check("forced_binary_metrics_reproduced", np.isclose(forced["balanced_accuracy"], 0.58169541) and np.isclose(forced["mcc"], 0.16768795), forced[["balanced_accuracy", "mcc"]].to_dict())
    check("selective_metrics_reproduced", int(selective["accepted_n"]) == 453 and np.isclose(selective["balanced_accuracy"], 0.73510345) and np.isclose(selective["mcc"], 0.49857011), selective[["accepted_n", "balanced_accuracy", "mcc"]].to_dict())
    check("rejected_layer_metrics_reproduced", int(rejected["accepted_n"]) == 2171 and np.isclose(rejected["balanced_accuracy"], 0.54746152), rejected[["accepted_n", "balanced_accuracy", "mcc"]].to_dict())
    check("selective_coverage_reproduced", np.isclose(selective["coverage"], 0.17263720), float(selective["coverage"]))
    check("bootstrap_replicates_operating_points", bool((operating["bootstrap_valid_replicates"] == 2000).all()), sorted(operating["bootstrap_valid_replicates"].unique().tolist()))

    expected_curve_rows = 2 * len(protocol["experiments"]["probability_margin_coverage_reliability"]["target_coverages"])
    check("coverage_curve_row_count", len(curve) == expected_curve_rows, len(curve))
    curve_ok = True
    for _, group in curve.groupby("population"):
        ordered = group.sort_values("target_coverage", ascending=False)
        curve_ok &= bool(np.all(np.diff(ordered["coverage"].to_numpy()) <= 1e-12))
        curve_ok &= bool((np.abs(ordered["coverage"] - ordered["target_coverage"]) <= 1.0 / ordered["denominator_n"] + 1e-12).all())
    check("coverage_curve_targets_respected", curve_ok, "coverage follows locked target grid within one row")
    check("coverage_curve_metric_bounds", bool(curve["balanced_accuracy"].between(0, 1).all() and curve["mcc"].between(-1, 1).all()), "BAcc in [0,1], MCC in [-1,1]")

    identity_cumulative = cumulative[cumulative["population"].eq("identity_clean")].sort_values("stage_index")
    check("cumulative_stage_count", len(identity_cumulative) == 7, identity_cumulative[["stage_index", "rule_added"]].to_dict(orient="records"))
    check("cumulative_acceptance_nonincreasing", bool(np.all(np.diff(identity_cumulative["accepted_n"].to_numpy()) <= 0)), identity_cumulative["accepted_n"].tolist())
    check("cumulative_final_matches_frozen", int(identity_cumulative.iloc[-1]["accepted_n"]) == 453 and np.isclose(identity_cumulative.iloc[-1]["balanced_accuracy"], selective["balanced_accuracy"]), identity_cumulative.iloc[-1][["accepted_n", "balanced_accuracy"]].to_dict())

    identity_loo = loo[loo["population"].eq("identity_clean")].set_index("label")
    check("leave_one_out_complete", len(identity_loo) == 7, sorted(identity_loo.index))
    check("local_rule_external_failure_detected", identity_loo.loc["without_local_neighborhood_agreement", "balanced_accuracy"] > identity_loo.loc["full_rule_set", "balanced_accuracy"], identity_loo.loc[["full_rule_set", "without_local_neighborhood_agreement"], ["coverage", "balanced_accuracy", "mcc"]].to_dict(orient="index"))
    check("model_agreement_rule_external_value_detected", identity_loo.loc["without_supported_model_direction_agreement", "balanced_accuracy"] < identity_loo.loc["full_rule_set", "balanced_accuracy"], identity_loo.loc[["full_rule_set", "without_supported_model_direction_agreement"], ["coverage", "balanced_accuracy", "mcc"]].to_dict(orient="index"))

    identity_enrichment = enrichment[enrichment["population"].eq("identity_clean")]
    check("enrichment_metrics_complete", set(identity_enrichment["metric"]) == {"accuracy", "balanced_accuracy", "mcc", "error_rate"}, sorted(identity_enrichment["metric"].tolist()))
    check("enrichment_ci_excludes_zero", bool((identity_enrichment["ci95_low"] > 0).all()), identity_enrichment[["metric", "favorable_difference", "ci95_low", "ci95_high"]].to_dict(orient="records"))
    check("enrichment_bootstrap_replicates", bool((identity_enrichment["valid_bootstrap_replicates"] == 2000).all()), sorted(identity_enrichment["valid_bootstrap_replicates"].unique().tolist()))

    calibration_all = calibration_summary[calibration_summary["population"].eq("identity_clean_all")].iloc[0]
    check("calibration_population_complete", int(calibration_all["n"]) == 2624, calibration_all[["n", "n_positive", "n_negative"]].to_dict())
    calibration_values = pd.to_numeric(
        calibration_all[["brier_score", "log_loss", "ece_equal_width_10", "calibration_intercept", "calibration_slope", "roc_auc", "average_precision"]],
        errors="coerce",
    ).to_numpy(dtype=float)
    check("calibration_metrics_finite", bool(np.isfinite(calibration_values).all()), calibration_all.to_dict())
    check("calibration_probability_bounds", bool(calibration_bins["mean_predicted_probability"].between(0, 1).all() and calibration_bins["observed_positive_fraction"].between(0, 1).all()), "all calibration points in [0,1]")
    check("calibration_bin_strategies", set(calibration_bins["strategy"]) == {"equal_width", "adaptive"}, sorted(calibration_bins["strategy"].unique().tolist()))

    required_strata = {"selection_layer", "risk_call", "applicability_status", "uncertainty_status", "abstention_reason"}
    check("error_strata_complete", required_strata.issubset(set(error_strata["stratum_type"])), sorted(error_strata["stratum_type"].unique().tolist()))
    check("error_metrics_bounded", bool(error_strata["accuracy"].between(0, 1).all() and error_strata["balanced_accuracy"].between(0, 1).all()), "accuracy and BAcc bounded")

    check("figure_manifest_count", len(manifest) == 9, len(manifest))
    check("figure_formats_complete", set(manifest["format"]) == {"png", "svg", "pdf"} and manifest["figure"].nunique() == 3, manifest.groupby("figure")["format"].apply(list).to_dict())
    for _, row in manifest.iterrows():
        path = PROJECT / row["path"]
        check(f"figure_exists_{path.name}", path.exists(), row["path"])
        check(f"figure_hash_{path.name}", path.exists() and sha256(path) == row["sha256"], row["sha256"])
        check(f"figure_nontrivial_{path.name}", path.exists() and path.stat().st_size >= 20_000, int(path.stat().st_size) if path.exists() else 0)

    result = {
        "validation_version": "ABCB1_EXTERNAL_RISK_EXPERIMENT_VALIDATION_V1",
        "validation_date": "2026-08-05",
        "overall_status": "pass" if all(item["pass"] for item in checks) else "fail",
        "checks_passed": int(sum(item["pass"] for item in checks)),
        "checks_total": int(len(checks)),
        "checks": checks,
    }
    VALIDATION.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ["overall_status", "checks_passed", "checks_total"]}, ensure_ascii=False, indent=2))
    if result["overall_status"] != "pass":
        failed = [item for item in checks if not item["pass"]]
        raise SystemExit(json.dumps(failed, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
