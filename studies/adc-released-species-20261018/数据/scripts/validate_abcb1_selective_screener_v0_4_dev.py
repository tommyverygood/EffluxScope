#!/usr/bin/env python3
"""Validate v0.4-dev artifacts, boundaries, and deterministic behavior."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path


DATA = Path(__file__).resolve().parents[1]
PROJECT = DATA.parent
PROCESSED = DATA / "processed"
MODEL_DIR = PROJECT / "模型" / "ABCB1_selective_screener_v0_4_dev"
sys.path.insert(0, str(MODEL_DIR))

from predict_abcb1_selective import load_packages, predict_one  # noqa: E402


SELECTIVE = MODEL_DIR / "ABCB1_selective_uncertainty_v0_4_dev.joblib"
PRIMARY = PROJECT / "模型" / "ABCB1_screener_v0_1" / "ABCB1_screener_v0_1.joblib"
DIRECT = PROJECT / "模型" / "ABCB1_assay_aware_candidate_v0_3_dev" / "ABCB1_assay_aware_candidate_v0_3_dev.joblib"
PROTOCOL = PROCESSED / "abcb1_selective_screener_v0_4_protocol_lock_v1.json"
BUILD_SUMMARY = PROCESSED / "abcb1_selective_screener_v0_4_build_summary_v1.json"
DIAGNOSTIC = PROCESSED / "abcb1_selective_v0_4_posthoc_diagnostic_summary_v1.json"
ADC_SUMMARY = PROCESSED / "abcb1_selective_v0_4_adc_release_registry_summary_v1.json"
ADC_PREDICTIONS = PROCESSED / "abcb1_selective_v0_4_adc_release_registry_predictions_v1.csv"
VALIDATION = PROCESSED / "abcb1_selective_screener_v0_4_validation_v1.json"
MANIFEST = MODEL_DIR / "manifest.json"

EXPECTED_PRIMARY_HASH = "fbc92dfa99583386ba124141fed5e77c6ecf8a59850301a6533db3883b298ac3"
EXPECTED_DIRECT_HASH = "ebe6d1a8d5f06afc5c2d0b306dfd5266af56fd28fefd130147d20c5e015b907a"
EXPECTED_PROTOCOL_HASH = "716d5285bed29d7caaca192b2fd76ea4dc645fe973ee3ec64d889187ac51a641"
EXPECTED_PACKAGE_HASH = "c56fb6b8f6f8bf13472a1a2d433b05d4c6a35eb2490ba84ce9371395b8bbff8b"

HIGH_EXAMPLE = "O=C(NCc1c(F)cc(F)cc1F)c1cn2c(c(O)c1=O)C(=O)N1[C@H]3CC[C@H](C3)O[C@@H]1C2"
LOW_EXAMPLE = "CC(=O)Nc1ccc(O)cc1"
ABSTAIN_EXAMPLE = "CCO"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def check(name: str, passed: bool, detail: str) -> dict:
    if not passed:
        raise RuntimeError(f"validation failed: {name}: {detail}")
    return {"check": name, "status": "pass", "detail": detail}


def close(a: float, b: float, tolerance: float = 1e-8) -> bool:
    return math.isclose(float(a), float(b), rel_tol=tolerance, abs_tol=tolerance)


def main() -> None:
    checks: list[dict] = []
    checks.append(check("primary_model_hash", sha256_file(PRIMARY) == EXPECTED_PRIMARY_HASH, sha256_file(PRIMARY)))
    checks.append(check("direct_model_hash", sha256_file(DIRECT) == EXPECTED_DIRECT_HASH, sha256_file(DIRECT)))
    checks.append(check("protocol_hash", sha256_file(PROTOCOL) == EXPECTED_PROTOCOL_HASH, sha256_file(PROTOCOL)))
    checks.append(check("selective_package_hash", sha256_file(SELECTIVE) == EXPECTED_PACKAGE_HASH, sha256_file(SELECTIVE)))

    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    build = json.loads(BUILD_SUMMARY.read_text(encoding="utf-8"))
    diagnostic = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    adc = json.loads(ADC_SUMMARY.read_text(encoding="utf-8"))
    selective, primary, direct = load_packages(SELECTIVE, PRIMARY, DIRECT)

    checks.append(check("model_version", selective["model_version"] == "ABCB1_selective_screener_v0.4_dev", selective["model_version"]))
    checks.append(check("bootstrap_estimators", len(selective["bootstrap_models"]) == 64, str(len(selective["bootstrap_models"]))))
    checks.append(check("thresholds_frozen", selective["primary_probability_low_threshold"] == 0.35 and selective["primary_probability_high_threshold"] == 0.70, "low=0.35 high=0.70"))
    checks.append(check("protocol_reject_not_biological_class", protocol["three_output_semantics"]["abstain"].startswith("insufficient"), protocol["three_output_semantics"]["abstain"]))
    checks.append(check("no_historical_external_rule_tuning", build["historical_external_outcomes_used_for_rule_selection"] is False and diagnostic["rule_tuned_on_historical_external_outcomes"] is False, diagnostic["diagnostic_status"]))
    checks.append(check("future_confirmation_required", build["future_source_confirmation_required"] is True, build["status"]))

    internal = build["probability_margin_internal_metrics"]
    checks.append(check("internal_n", internal["n"] == 443, str(internal["n"])))
    checks.append(check("internal_probability_margin_coverage", close(internal["coverage"], 0.6117381489841986), str(internal["coverage"])))
    checks.append(check("internal_probability_margin_accuracy", close(internal["accepted_accuracy"], 0.8523985239852399), str(internal["accepted_accuracy"])))

    high = predict_one(HIGH_EXAMPLE, selective, primary, direct)
    low = predict_one(LOW_EXAMPLE, selective, primary, direct)
    abstain = predict_one(ABSTAIN_EXAMPLE, selective, primary, direct)
    invalid = predict_one("not_a_smiles", selective, primary, direct)
    intact_adc = predict_one(HIGH_EXAMPLE, selective, primary, direct, chemical_object_type="intact_ADC")
    clinical_ddi = predict_one(HIGH_EXAMPLE, selective, primary, direct, requested_context="clinical_DDI")
    checks.append(check("high_output_reachable", high["risk_call"] == "high_ABCB1_efflux_risk", high["risk_call"]))
    checks.append(check("low_output_reachable", low["risk_call"] == "low_ABCB1_efflux_risk", low["risk_call"]))
    checks.append(check("abstain_output_reachable", abstain["risk_call"] == "abstain" and bool(abstain["abstention_reasons"]), json.dumps(abstain["abstention_reasons"])))
    checks.append(check("invalid_structure_abstains", invalid["abstention_reasons"] == ["unresolved_structure"], json.dumps(invalid["abstention_reasons"])))
    checks.append(check("intact_adc_context_abstains", intact_adc["risk_call"] == "abstain" and "object_scope_incompatible" in intact_adc["abstention_reasons"], json.dumps(intact_adc["abstention_reasons"])))
    checks.append(check("clinical_ddi_context_abstains", clinical_ddi["risk_call"] == "abstain" and "requested_context_incompatible" in clinical_ddi["abstention_reasons"], json.dumps(clinical_ddi["abstention_reasons"])))
    checks.append(check("deterministic_single_prediction", json.dumps(high, sort_keys=True) == json.dumps(predict_one(HIGH_EXAMPLE, selective, primary, direct), sort_keys=True), "identical repeated JSON"))

    reg = diagnostic["regulatory_v4"]
    asap = diagnostic["asap_2025_weak_binary"]
    conrad = diagnostic["conrad_figure5"]
    checks.append(check("regulatory_posthoc_counts", reg["n"] == 33 and reg["accepted_n"] == 13 and reg["abstain_n"] == 20, json.dumps(reg)))
    checks.append(check("asap_posthoc_counts", asap["n"] == 62 and asap["accepted_n"] == 5 and asap["low_n"] == 0, json.dumps(asap)))
    checks.append(check("conrad_negative_result_preserved", conrad["identity_clean_pair_n"] == 11 and conrad["pair_direction_correct_n"] == 6, json.dumps(conrad)))

    adc_rows = read_csv(ADC_PREDICTIONS)
    allowed_calls = {"high_ABCB1_efflux_risk", "low_ABCB1_efflux_risk", "abstain"}
    checks.append(check("adc_registry_rows", len(adc_rows) == 10 and adc["explicit_structure_rows"] == 9, f"rows={len(adc_rows)} explicit={adc['explicit_structure_rows']}"))
    checks.append(check("adc_calls_allowed", {row["risk_call"] for row in adc_rows}.issubset(allowed_calls), str(sorted({row["risk_call"] for row in adc_rows}))))
    checks.append(check("accepted_rows_have_no_reasons", all(json.loads(row["abstention_reasons_json"]) == [] for row in adc_rows if row["risk_call"] != "abstain"), "all accepted rows clean"))
    checks.append(check("abstain_rows_have_reasons", all(bool(json.loads(row["abstention_reasons_json"])) for row in adc_rows if row["risk_call"] == "abstain"), "all abstain rows explained"))
    checks.append(check("adc_call_counts", adc["risk_call_counts"] == {"abstain": 7, "high_ABCB1_efflux_risk": 3}, json.dumps(adc["risk_call_counts"])))

    required_artifacts = [
        SELECTIVE,
        MODEL_DIR / "predict_abcb1_selective.py",
        MODEL_DIR / "README.md",
        MODEL_DIR / "MODEL_CARD.md",
        MODEL_DIR / "runtime-lock.txt",
        MODEL_DIR / "test_smoke.py",
        MODEL_DIR / "example_queries.csv",
        MODEL_DIR / "example_predictions.csv",
        PROJECT / "ABCB1三结果风险筛选器_v0.4开发报告_20260803.md",
        DATA / "scripts" / "build_abcb1_selective_screener_v0_4_dev.py",
        DATA / "scripts" / "diagnose_abcb1_selective_screener_v0_4_dev.py",
        DATA / "scripts" / "score_abcb1_selective_adc_release_registry_v0_4_dev.py",
        DATA / "scripts" / "validate_abcb1_selective_screener_v0_4_dev.py",
        PROTOCOL,
        BUILD_SUMMARY,
        DIAGNOSTIC,
        ADC_PREDICTIONS,
        ADC_SUMMARY,
        PROCESSED / "abcb1_selective_internal_metrics_v1.csv",
        PROCESSED / "abcb1_selective_probability_threshold_sweep_v1.csv",
        PROCESSED / "abcb1_selective_v0_4_posthoc_diagnostic_metrics_v1.csv",
        PROCESSED / "abcb1_selective_v0_4_regulatory_v4_posthoc_predictions_v1.csv",
        PROCESSED / "abcb1_selective_v0_4_asap2025_posthoc_predictions_v1.csv",
        PROCESSED / "abcb1_selective_v0_4_conrad_posthoc_predictions_v1.csv",
        PROCESSED / "abcb1_selective_v0_4_conrad_pairwise_posthoc_v1.csv",
    ]
    checks.append(check("required_artifacts_present", all(path.exists() and path.stat().st_size > 0 for path in required_artifacts), f"n={len(required_artifacts)}"))

    validation = {
        "validation_status": "pass",
        "checks_passed": len(checks),
        "checks_total": len(checks),
        "checks": checks,
        "model_version": selective["model_version"],
        "model_package_sha256": sha256_file(SELECTIVE),
        "operational_status": "three_output_development_screener_operational",
        "external_validation_status": "not_future_source_validated",
        "claim_boundary": selective["claim_boundary"],
    }
    VALIDATION.write_text(json.dumps(validation, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")

    artifacts = {
        str(path.relative_to(PROJECT)): sha256_file(path)
        for path in [*required_artifacts, VALIDATION]
    }
    manifest = {
        "model_version": selective["model_version"],
        "status": "development_selective_risk_screener_not_future_source_validated",
        "validated_by": str(VALIDATION.relative_to(PROJECT)),
        "validation_sha256": sha256_file(VALIDATION),
        "artifacts": artifacts,
        "entrypoint": "模型/ABCB1_selective_screener_v0_4_dev/predict_abcb1_selective.py",
        "claim_boundary": selective["claim_boundary"],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
