#!/usr/bin/env python3
"""Validate the source-traceable ABCB1 selective-screener manuscript package."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"
PLAN = PROJECT / "论文计划" / "论文一_ABCB1不确定性感知选择性筛查_主稿实施方案_20260805.md"
OUTPUT = PROCESSED / "abcb1_manuscript_package_validation_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def close(value: float, expected: float, tolerance: float = 1e-7) -> bool:
    return abs(float(value) - expected) <= tolerance


def main() -> None:
    results_path = PROCESSED / "abcb1_manuscript_primary_results_v1.csv"
    claims_path = PROCESSED / "abcb1_manuscript_claim_ledger_v1.csv"
    figures_path = PROCESSED / "abcb1_manuscript_figure_manifest_v1.csv"
    tables_path = PROCESSED / "abcb1_manuscript_table_manifest_v1.csv"
    readiness_path = PROCESSED / "abcb1_manuscript_readiness_summary_v1.json"
    results = pd.read_csv(results_path)
    claims = pd.read_csv(claims_path)
    figures = pd.read_csv(figures_path)
    tables = pd.read_csv(tables_path)
    readiness = json.loads(readiness_path.read_text())
    plan = PLAN.read_text()

    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: object) -> None:
        checks.append({"check": name, "status": "PASS" if condition else "FAIL", "detail": detail})

    check("primary_results_row_count", len(results) == 18, len(results))
    check("primary_result_ids_unique", results["result_id"].nunique() == len(results), results["result_id"].tolist())
    binary = results[results["result_id"].eq("R1_BINARY_BASELINE_IDENTITY_CLEAN")].iloc[0]
    selective = results[results["result_id"].eq("R2_SELECTIVE_PRIMARY_IDENTITY_CLEAN")].iloc[0]
    forced = results[results["result_id"].eq("R3_FORCED_REJECTED_LAYER_IDENTITY_CLEAN")].iloc[0]
    continuous = results[results["result_id"].eq("R6_CONTINUOUS_DIRECTION_IDENTITY_CLEAN")].iloc[0]
    daood = results[results["result_id"].eq("R8_DAOOD_LIMITATION")].iloc[0]
    enrichment = results[results["result_id"].eq("R12_SELECTIVE_ENRICHMENT_DIFFERENCE_POSTHOC")].iloc[0]
    calibration = results[results["result_id"].eq("R13_EXTERNAL_CALIBRATION_POSTHOC")].iloc[0]
    local_ablation = results[results["result_id"].eq("R14_LOCAL_RULE_ABLATION_HYPOTHESIS")].iloc[0]
    direct_channel = results[results["result_id"].eq("R15_DIRECT_EVIDENCE_CHANNEL_CONTRIBUTION")].iloc[0]
    raw_calibration = results[results["result_id"].eq("R16_RAW_EXTERNAL_CALIBRATION_LIMITATION")].iloc[0]
    platt_diagnostic = results[results["result_id"].eq("R17_WITHIN_SOURCE_PLATT_DIAGNOSTIC")].iloc[0]
    error_concentration = results[results["result_id"].eq("R18_FALSE_POSITIVE_ERROR_CONCENTRATION")].iloc[0]
    check("binary_baseline_anchor", int(binary["n_evaluated"]) == 2624 and close(binary["primary_value"], 0.58169541) and close(binary["mcc"], 0.16768795), binary.to_dict())
    check("selective_anchor", int(selective["n_evaluated"]) == 453 and close(selective["coverage"], 0.1726372) and close(selective["primary_value"], 0.73510345) and close(selective["mcc"], 0.49857011), selective.to_dict())
    check("forced_layer_anchor", int(forced["n_evaluated"]) == 2171 and close(forced["primary_value"], 0.54746152) and close(forced["mcc"], 0.09621794), forced.to_dict())
    check("continuous_anchor", int(continuous["n_evaluated"]) == 2624 and close(continuous["primary_value"], 0.24718176) and continuous["ci95_low"] > 0, continuous.to_dict())
    check("selective_outperforms_binary", selective["primary_value"] > binary["primary_value"] and selective["mcc"] > binary["mcc"] and selective["coverage"] < binary["coverage"], [binary["primary_value"], selective["primary_value"], binary["coverage"], selective["coverage"]])
    check("forced_layer_is_weaker", forced["primary_value"] < binary["primary_value"] and forced["mcc"] < binary["mcc"], [forced["primary_value"], forced["mcc"]])
    check("daood_limitation_preserved", int(daood["n_total"]) == 156 and "only 4 positive" in daood["analysis_subset"], daood[["n_total", "analysis_subset"]].to_dict())
    check("selective_enrichment_anchor", close(enrichment["primary_value"], 0.18764193) and enrichment["ci95_low"] > 0 and close(enrichment["mcc"], 0.40235217), enrichment.to_dict())
    check("external_calibration_anchor", close(calibration["primary_value"], 0.58367794) and "not an absolute probability" in calibration["allowed_interpretation"], calibration.to_dict())
    check("local_ablation_hypothesis_preserved", int(local_ablation["n_evaluated"]) == 551 and close(local_ablation["primary_value"], 0.76329591) and "requires new frozen validation" in local_ablation["evidence_timing"], local_ablation.to_dict())
    check("direct_channel_contribution_anchor", close(direct_channel["primary_value"], 0.05702906) and direct_channel["ci95_low"] > 0 and close(direct_channel["mcc"], 0.13951992), direct_channel.to_dict())
    check("raw_calibration_anchor_v2", close(raw_calibration["primary_value"], 0.58381812) and raw_calibration["ci95_high"] < 1 and "do not interpret" in raw_calibration["allowed_interpretation"], raw_calibration.to_dict())
    check("platt_diagnostic_anchor", close(platt_diagnostic["primary_value"], 0.01818546) and platt_diagnostic["ci95_low"] > 0 and "not independent validation" in platt_diagnostic["allowed_interpretation"], platt_diagnostic.to_dict())
    check("error_concentration_anchor", int(error_concentration["n_evaluated"]) == 1050 and close(error_concentration["primary_value"], 0.53619048) and "824 false positives" in error_concentration["allowed_interpretation"], error_concentration.to_dict())

    check("claim_ledger_has_allowed_and_forbidden", {"allowed", "allowed_with_boundary", "forbidden"}.issubset(set(claims["status"])), claims["status"].value_counts().to_dict())
    check("abstention_not_biological_class", claims["manuscript_wording"].str.contains("third biological class", case=False).any(), claims[claims["manuscript_wording"].str.contains("third biological class", case=False)]["claim_id"].tolist())
    check("adc_accuracy_claim_forbidden", claims[(claims["claim_id"] == "C9") & (claims["status"] == "forbidden")].shape[0] == 1, claims[claims["claim_id"] == "C9"].to_dict("records"))
    check("figure_manifest_complete", len(figures) == 11 and figures["figure_id"].str.startswith("Figure").all(), figures[["figure_id", "status"]].to_dict("records"))
    check("main_figure_data_ready", figures[figures["role"].isin(["primary", "required"])]["status"].isin(["data_ready", "ready_for_drawing", "figure_ready"]).all(), figures[figures["role"].isin(["primary", "required"])][["figure_id", "status"]].to_dict("records"))
    check("external_risk_figures_ready", set(figures.loc[(figures["status"].eq("figure_ready")) & (figures["role"].eq("primary")), "figure_id"]) == {"Figure 2", "Figure 3", "Figure 4"}, figures[figures["status"].eq("figure_ready")][["figure_id", "source_artifacts"]].to_dict("records"))
    check("evidence_channel_figure_ready", figures[(figures["figure_id"] == "Figure S4") & (figures["status"] == "figure_ready")].shape[0] == 1, figures[figures["figure_id"] == "Figure S4"].to_dict("records"))
    check("deep_calibration_figures_ready", set(figures.loc[(figures["status"] == "figure_ready") & (figures["figure_id"].isin(["Figure S5", "Figure S6"])), "figure_id"]) == {"Figure S5", "Figure S6"}, figures[figures["figure_id"].isin(["Figure S5", "Figure S6"])].to_dict("records"))
    check("table_manifest_complete", len(tables) == 12, len(tables))

    expected_gates = {
        "large_external_continuous_direction": "pass",
        "identity_clean_analysis": "pass",
        "scaffold_clean_sensitivity": "pass",
        "selective_reliability_improvement": "pass",
        "posthoc_external_risk_experiments": "pass",
        "posthoc_evidence_channel_ablation": "pass",
        "posthoc_deep_calibration_error_analysis": "pass",
        "binary_full_coverage_high_discrimination": "fail",
        "future_strict_matched_parental_or_inhibitor_validation": "pending",
        "adc_specific_external_validation": "fail",
    }
    check("readiness_gates_exact", readiness["gates"] == expected_gates, readiness["gates"])
    check("bounded_draft_go", readiness["manuscript_readiness"] == "go_for_draft_under_bounded_methods_claim", readiness["manuscript_readiness"])
    check("primary_method_is_selective", "v0.4 selective" in readiness["main_method"], readiness["main_method"])
    check("binary_is_baseline", "baseline" not in readiness["main_method"].lower() and "weak-label" in readiness["baseline"], readiness["baseline"])

    for relative, expected_hash in readiness["outputs"].items():
        path = PROJECT / relative
        check(f"readiness_output_hash::{relative}", path.exists() and sha256(path) == expected_hash, expected_hash)

    required_plan_tokens = [
        "An uncertainty-aware selective screening framework for ABCB1 efflux risk",
        "2,624",
        "0.582",
        "0.735",
        "17.3%",
        "2,171",
        "0.547",
        "拒判是决策状态",
        "不是第三种生物学类别",
        "ADC 释放物仅作风险分诊案例",
    ]
    for token in required_plan_tokens:
        check(f"plan_token::{token}", token in plan, token)
    check("plan_exists_and_substantial", PLAN.exists() and PLAN.stat().st_size > 10000, {"path": str(PLAN), "bytes": PLAN.stat().st_size})

    status = "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL"
    result = {
        "validation_version": "ABCB1_SELECTIVE_MANUSCRIPT_PACKAGE_VALIDATION_V1",
        "status": status,
        "checks_passed": sum(item["status"] == "PASS" for item in checks),
        "checks_total": len(checks),
        "manuscript_readiness": readiness["manuscript_readiness"],
        "plan_sha256": sha256(PLAN),
        "checks": checks,
    }
    OUTPUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
