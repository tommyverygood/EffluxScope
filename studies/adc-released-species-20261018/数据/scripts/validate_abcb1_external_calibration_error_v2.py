#!/usr/bin/env python3
"""Validate deep external calibration and error diagnostics for ABCB1."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"
PROTOCOL = PROCESSED / "abcb1_external_calibration_error_protocol_v2.json"
SUMMARY = PROCESSED / "abcb1_external_calibration_error_summary_v2.json"
VALIDATION = PROCESSED / "abcb1_external_calibration_error_validation_v2.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text())
    summary = json.loads(SUMMARY.read_text())
    predictions = pd.read_csv(PROCESSED / "abcb1_external_crossfitted_calibration_predictions_v2.csv")
    metrics = pd.read_csv(PROCESSED / "abcb1_external_calibration_metrics_v2.csv")
    bootstrap = pd.read_csv(PROCESSED / "abcb1_external_calibration_bootstrap_v2.csv")
    pairwise = pd.read_csv(PROCESSED / "abcb1_external_calibration_pairwise_bootstrap_v2.csv")
    bins = pd.read_csv(PROCESSED / "abcb1_external_calibration_bins_v2.csv")
    brier = pd.read_csv(PROCESSED / "abcb1_external_brier_decomposition_v2.csv")
    counts = pd.read_csv(PROCESSED / "abcb1_external_error_counts_v2.csv")
    features = pd.read_csv(PROCESSED / "abcb1_external_error_feature_comparisons_v2.csv")
    strata = pd.read_csv(PROCESSED / "abcb1_external_error_strata_v2.csv")
    hotspots = pd.read_csv(PROCESSED / "abcb1_external_error_scaffold_hotspots_v2.csv")
    figures = pd.read_csv(PROCESSED / "abcb1_external_calibration_error_figure_manifest_v2.csv")

    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, evidence: object) -> None:
        checks.append({"check": name, "pass": bool(condition), "evidence": evidence})

    check("protocol_truthful_posthoc", protocol["status"].startswith("locked_posthoc_deep_diagnostic"), protocol["status"])
    check("summary_completed", summary["status"] == "completed_frozen_model_posthoc_deep_calibration_and_error_diagnostic", summary["status"])
    check("protocol_hash_recorded", summary["protocol_sha256"] == sha256(PROTOCOL), summary["protocol_sha256"])
    for label, entry in protocol["frozen_inputs"].items():
        path = PROJECT / entry["path"]
        check(f"frozen_input_exists_{label}", path.exists(), entry["path"])
        check(f"frozen_input_hash_{label}", path.exists() and sha256(path) == entry["sha256"], entry["sha256"])
    for relative, expected_hash in summary["output_hashes"].items():
        path = PROJECT / relative
        check(f"output_exists_{path.name}", path.exists(), relative)
        check(f"output_hash_{path.name}", path.exists() and sha256(path) == expected_hash, expected_hash)

    check("prediction_row_count", len(predictions) == 2624 + 2562, len(predictions))
    check("prediction_population_counts", predictions["population"].value_counts().to_dict() == {"identity_clean": 2624, "scaffold_clean": 2562}, predictions["population"].value_counts().to_dict())
    check("five_folds_each_population", all(set(frame["fold_id"]) == {1, 2, 3, 4, 5} for _, frame in predictions.groupby("population")), predictions.groupby("population")["fold_id"].unique().apply(list).to_dict())
    cluster_fold_max = predictions.groupby(["population", "bootstrap_cluster"])["fold_id"].nunique().max()
    check("group_crossfit_no_cluster_split", int(cluster_fold_max) == 1, int(cluster_fold_max))
    probability_columns = ["raw_probability", "platt_probability", "isotonic_probability"]
    check("calibration_probabilities_finite_bounded", bool(np.isfinite(predictions[probability_columns].to_numpy(dtype=float)).all() and predictions[probability_columns].apply(lambda column: column.between(0, 1).all()).all()), predictions[probability_columns].describe().to_dict())

    check("calibration_metric_rows", len(metrics) == 6 and set(metrics["method"]) == {"raw_frozen_score", "cross_fitted_platt", "cross_fitted_isotonic"}, metrics[["population", "method"]].to_dict(orient="records"))
    identity = metrics[metrics["population"].eq("identity_clean")].set_index("method")
    raw = identity.loc["raw_frozen_score"]
    platt = identity.loc["cross_fitted_platt"]
    isotonic = identity.loc["cross_fitted_isotonic"]
    check("raw_calibration_anchors", np.isclose(raw["brier_score"], 0.25296644) and np.isclose(raw["ece_equal_width_10"], 0.12074011) and np.isclose(raw["calibration_slope"], 0.58381812), raw.to_dict())
    check("platt_calibration_anchors", np.isclose(platt["brier_score"], 0.23478098) and np.isclose(platt["ece_equal_width_10"], 0.01180836) and abs(platt["calibration_intercept"]) < 0.01 and abs(platt["calibration_slope"] - 1) < 0.03, platt.to_dict())
    check("isotonic_not_preferred_over_platt", isotonic["log_loss"] > platt["log_loss"] and isotonic["roc_auc"] < platt["roc_auc"], identity[["brier_score", "log_loss", "ece_equal_width_10", "roc_auc"]].to_dict(orient="index"))

    check("bootstrap_row_count", len(bootstrap) == 2 * 3 * 8, len(bootstrap))
    check("bootstrap_complete", bool((bootstrap["valid_bootstrap_replicates"] == 2000).all()), sorted(bootstrap["valid_bootstrap_replicates"].unique().tolist()))
    raw_slope = bootstrap[(bootstrap["population"] == "identity_clean") & (bootstrap["method"] == "raw_frozen_score") & (bootstrap["metric"] == "calibration_slope")].iloc[0]
    raw_intercept = bootstrap[(bootstrap["population"] == "identity_clean") & (bootstrap["method"] == "raw_frozen_score") & (bootstrap["metric"] == "calibration_intercept")].iloc[0]
    raw_ece = bootstrap[(bootstrap["population"] == "identity_clean") & (bootstrap["method"] == "raw_frozen_score") & (bootstrap["metric"] == "ece_equal_width_10")].iloc[0]
    check("raw_slope_ci_excludes_one", raw_slope["ci95_high"] < 1, raw_slope.to_dict())
    check("raw_intercept_ci_excludes_zero", raw_intercept["ci95_high"] < 0, raw_intercept.to_dict())
    check("raw_ece_ci_positive", raw_ece["ci95_low"] > 0.10, raw_ece.to_dict())

    check("pairwise_row_count", len(pairwise) == 2 * 2 * 8, len(pairwise))
    check("pairwise_bootstrap_complete", bool((pairwise["valid_bootstrap_replicates"] == 2000).all()), sorted(pairwise["valid_bootstrap_replicates"].unique().tolist()))
    platt_pair = pairwise[(pairwise["population"] == "identity_clean") & (pairwise["calibrated_method"] == "cross_fitted_platt")].set_index("metric")
    check("platt_brier_improvement_ci", platt_pair.loc["brier_score", "ci95_low"] > 0, platt_pair.loc["brier_score"].to_dict())
    check("platt_ece_improvement_ci", platt_pair.loc["ece_equal_width_10", "ci95_low"] > 0, platt_pair.loc["ece_equal_width_10"].to_dict())
    check("platt_does_not_improve_discrimination", platt_pair.loc["roc_auc", "ci95_high"] < 0 and platt_pair.loc["average_precision", "ci95_high"] < 0, platt_pair.loc[["roc_auc", "average_precision"]].to_dict(orient="index"))

    check("calibration_bins_complete", set(bins["strategy"]) == {"equal_width", "adaptive"} and set(bins["method"]) == {"raw_frozen_score", "cross_fitted_platt", "cross_fitted_isotonic"}, {"rows": len(bins)})
    check("calibration_bin_bounds", bool(bins["mean_predicted_probability"].between(0, 1).all() and bins["observed_positive_fraction"].between(0, 1).all()), "all bin values bounded")
    check("brier_decomposition_rows", len(brier) == 6 and np.isfinite(brier[["reliability", "resolution", "uncertainty", "binning_residual"]].to_numpy(dtype=float)).all(), brier.to_dict(orient="records"))

    identity_counts = counts[counts["population"].eq("identity_clean")].set_index("error_outcome")
    expected_counts = {"true_positive": 825, "true_negative": 655, "false_positive": 824, "false_negative": 320}
    check("error_taxonomy_counts", {key: int(identity_counts.loc[key, "n"]) for key in expected_counts} == expected_counts, identity_counts["n"].to_dict())
    check("error_counts_sum", int(identity_counts["n"].sum()) == 2624, int(identity_counts["n"].sum()))
    check("false_positives_dominate_errors", int(identity_counts.loc["false_positive", "n"]) > 2 * int(identity_counts.loc["false_negative", "n"]), expected_counts)

    identity_strata = strata[strata["population"].eq("identity_clean")]
    accepted = identity_strata[(identity_strata["stratum_type"] == "selection_layer") & (identity_strata["stratum_level"] == "v0_4_accepted")].iloc[0]
    rejected = identity_strata[(identity_strata["stratum_type"] == "selection_layer") & (identity_strata["stratum_level"] == "v0_4_rejected")].iloc[0]
    band = identity_strata[(identity_strata["stratum_type"] == "probability_band") & (identity_strata["stratum_level"] == "[0.50,0.70)")].iloc[0]
    neither = identity_strata[(identity_strata["stratum_type"] == "support_combination") & (identity_strata["stratum_level"] == "neither_supported")].iloc[0]
    check("selection_layer_error_gap", accepted["error_rate"] < rejected["error_rate"] and np.isclose(accepted["error_rate"], 0.24944812) and np.isclose(rejected["error_rate"], 0.47489636), {"accepted": accepted["error_rate"], "rejected": rejected["error_rate"]})
    check("gray_high_band_error_peak", np.isclose(band["error_rate"], 0.53619048) and int(band["false_positive_n"]) == 563, band.to_dict())
    check("outside_support_error_high", np.isclose(neither["error_rate"], 0.51815182), neither.to_dict())

    check("feature_comparisons_present", set(features["comparison"]) == {"all_error_vs_correct", "false_positive_vs_true_positive", "false_negative_vs_true_negative"}, features["comparison"].value_counts().to_dict())
    check("feature_fdr_computed", int(features["fdr_bh_within_comparison"].notna().sum()) >= len(features) - 4, {"finite": int(features["fdr_bh_within_comparison"].notna().sum()), "total": len(features)})
    fp_direct = features[(features["population"] == "identity_clean") & (features["comparison"] == "false_positive_vs_true_positive") & (features["feature"] == "direct_probability")].iloc[0]
    fn_heavy = features[(features["population"] == "identity_clean") & (features["comparison"] == "false_negative_vs_true_negative") & (features["feature"] == "heavy_atoms")].iloc[0]
    check("false_positive_direct_score_effect", fp_direct["cliffs_delta"] < -0.45 and fp_direct["fdr_bh_within_comparison"] < 0.05, fp_direct.to_dict())
    check("false_negative_size_effect", fn_heavy["cliffs_delta"] > 0.45 and fn_heavy["fdr_bh_within_comparison"] < 0.05, fn_heavy.to_dict())

    check("scaffold_hotspots_minimum_size", bool((hotspots["n"] >= 4).all()), hotspots["n"].describe().to_dict())
    check("scaffold_hotspot_rates_bounded", bool(hotspots["raw_error_rate"].between(0, 1).all() and hotspots["shrunken_error_rate"].between(0, 1).all()), "raw and shrunken rates bounded")
    check("scaffold_shrinkage_operates", bool((np.abs(hotspots["raw_error_rate"] - hotspots["shrunken_error_rate"]) > 0).all()), "all retained hotspots shrunk toward population mean")

    check("figure_manifest_complete", len(figures) == 6 and set(figures["figure"]) == {"Figure S5", "Figure S6"} and set(figures["format"]) == {"png", "svg", "pdf"}, figures[["figure", "format", "path"]].to_dict(orient="records"))
    for _, row in figures.iterrows():
        path = PROJECT / row["path"]
        check(f"figure_exists_{path.name}", path.exists(), row["path"])
        check(f"figure_hash_{path.name}", path.exists() and sha256(path) == row["sha256"], row["sha256"])
        check(f"figure_nontrivial_{path.name}", path.exists() and path.stat().st_size >= 20_000, int(path.stat().st_size) if path.exists() else 0)

    result = {
        "validation_version": "ABCB1_EXTERNAL_CALIBRATION_ERROR_VALIDATION_V2",
        "validation_date": "2026-08-05",
        "overall_status": "pass" if all(item["pass"] for item in checks) else "fail",
        "checks_passed": int(sum(item["pass"] for item in checks)),
        "checks_total": int(len(checks)),
        "checks": checks,
    }
    serializer = lambda value: value.item() if isinstance(value, np.generic) else value
    VALIDATION.write_text(json.dumps(result, ensure_ascii=False, indent=2, default=serializer) + "\n")
    print(json.dumps({key: result[key] for key in ["overall_status", "checks_passed", "checks_total"]}, ensure_ascii=False, indent=2))
    if result["overall_status"] != "pass":
        raise SystemExit(json.dumps([item for item in checks if not item["pass"]], ensure_ascii=False, indent=2, default=serializer))


if __name__ == "__main__":
    main()
