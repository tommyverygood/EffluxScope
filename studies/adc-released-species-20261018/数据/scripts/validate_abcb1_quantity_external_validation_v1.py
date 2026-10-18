#!/usr/bin/env python3
"""Validate the ABCB1 quantity-oriented external validation deliverables."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"
OUTPUT = PROJECT / "outputs" / "20260805_abcb1_quantity_external_validation_v1"
VALIDATION_OUT = PROCESSED / "abcb1_quantity_external_validation_validation_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def close(value: float, expected: float, tolerance: float = 1e-7) -> bool:
    return abs(float(value) - expected) <= tolerance


def main() -> None:
    summary_path = PROCESSED / "abcb1_quantity_external_validation_summary_v1.json"
    protocol_path = PROCESSED / "abcb1_quantity_external_validation_protocol_lock_v1.json"
    summary = json.loads(summary_path.read_text())
    protocol = json.loads(protocol_path.read_text())
    biogen = pd.read_csv(PROCESSED / "abcb1_quantity_biogen_analysis_rows_v1.csv", low_memory=False)
    daood = pd.read_csv(PROCESSED / "abcb1_quantity_daood_analysis_rows_v1.csv", low_memory=False)
    metrics = pd.read_csv(PROCESSED / "abcb1_quantity_external_metrics_v1.csv")
    workbook_validation = json.loads((OUTPUT / "workbook_validation_v1.json").read_text())

    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: object) -> None:
        checks.append({"check": name, "status": "PASS" if condition else "FAIL", "detail": detail})

    check("protocol_precedes_scoring_status", protocol["status"] == "locked_before_new_external_outcome_scoring", protocol["status"])
    check("protocol_hash_matches_summary", sha256(protocol_path) == summary["protocol_sha256"], summary["protocol_sha256"])
    check("model_hash_frozen", protocol["frozen_model"]["sha256"] == summary["model_sha256"], summary["model_sha256"])
    check("raw_scored_rows", summary["new_data_quantity"]["raw_scored_rows"] == 3249, summary["new_data_quantity"]["raw_scored_rows"])
    check("collapsed_unique_rows", len(biogen) + len(daood) == 3248, {"biogen": len(biogen), "daood": len(daood)})
    check("biogen_identity_clean_n", int(biogen["identity_clean"].eq("yes").sum()) == 2624, int(biogen["identity_clean"].eq("yes").sum()))
    check("biogen_scaffold_clean_n", int(biogen["scaffold_clean"].eq("yes").sum()) == 2562, int(biogen["scaffold_clean"].eq("yes").sum()))
    check("daood_identity_clean_n", int(daood["identity_clean"].eq("yes").sum()) == 156, int(daood["identity_clean"].eq("yes").sum()))
    check("daood_identity_clean_positive_n", int((daood.loc[daood["identity_clean"].eq("yes"), "author_pgp_substrate_label"] == 1).sum()) == 4, int((daood.loc[daood["identity_clean"].eq("yes"), "author_pgp_substrate_label"] == 1).sum()))

    for name, frame in [("biogen", biogen), ("daood", daood)]:
        clean = frame[frame["identity_clean"].eq("yes")]
        overlap = (
            clean["exact_connectivity_overlap_primary_443"].eq("yes")
            | clean["exact_connectivity_overlap_direct_84"].eq("yes")
        )
        check(f"{name}_identity_clean_has_no_exact_reference_overlap", not overlap.any(), int(overlap.sum()))
        scaffold_clean = frame[frame["scaffold_clean"].eq("yes")]
        scaffold_overlap = (
            scaffold_clean["murcko_overlap_primary_443"].eq("yes")
            | scaffold_clean["murcko_overlap_direct_84"].eq("yes")
        )
        check(f"{name}_scaffold_clean_has_no_reference_scaffold_overlap", not scaffold_overlap.any(), int(scaffold_overlap.sum()))

    primary = metrics[
        (metrics["source"] == "Biogen_Fang_2023")
        & (metrics["subset"] == "identity_clean")
        & (metrics["stratum"] == "all")
    ].iloc[0]
    scaffold = metrics[
        (metrics["source"] == "Biogen_Fang_2023")
        & (metrics["subset"] == "scaffold_clean")
        & (metrics["stratum"] == "all")
    ].iloc[0]
    check("primary_spearman", close(primary["spearman_probability_vs_endpoint"], 0.24718176), primary["spearman_probability_vs_endpoint"])
    check("primary_spearman_ci_excludes_zero", primary["spearman_ci95_low"] > 0 and primary["spearman_ci95_high"] > primary["spearman_ci95_low"], [primary["spearman_ci95_low"], primary["spearman_ci95_high"]])
    check("primary_coverage", close(primary["coverage"], 0.1726372), primary["coverage"])
    check("primary_accepted_n", int(primary["accepted_n"]) == 453, int(primary["accepted_n"]))
    check("primary_balanced_accuracy", close(primary["balanced_accuracy"], 0.73510345), primary["balanced_accuracy"])
    check("primary_mcc", close(primary["mcc"], 0.49857011), primary["mcc"])
    check("primary_high_endpoint_exceeds_low", primary["accepted_high_endpoint_median"] > primary["accepted_low_endpoint_median"], [primary["accepted_high_endpoint_median"], primary["accepted_low_endpoint_median"]])
    check("scaffold_clean_direction_robust", scaffold["spearman_ci95_low"] > 0 and scaffold["mcc"] > 0, [scaffold["spearman_ci95_low"], scaffold["mcc"]])
    check("predeclared_quantity_gate_pass", summary["predeclared_quantity_gate"]["quantity_support_status"] == "pass", summary["predeclared_quantity_gate"])
    check("strict_upgrade_remains_closed", summary["strict_predictor_upgrade"]["status"] == "not_permitted", summary["strict_predictor_upgrade"])

    for relative, expected_hash in summary["output_hashes"].items():
        path = PROJECT / relative
        check(f"output_hash::{relative}", path.exists() and sha256(path) == expected_hash, expected_hash)

    report = PROJECT / "ABCB1数量型外部风险验证报告_v1_20260805.md"
    workbook = OUTPUT / "ABCB1数量型外部风险验证_20260805.xlsx"
    check("report_exists", report.exists() and report.stat().st_size > 1000, str(report))
    check("workbook_exists", workbook.exists() and workbook.stat().st_size > 10000, str(workbook))
    check("workbook_post_export_validation", workbook_validation["status"] == "PASS", workbook_validation)

    status = "PASS" if all(row["status"] == "PASS" for row in checks) else "FAIL"
    result = {
        "validation_version": "ABCB1_QUANTITY_EXTERNAL_VALIDATION_OUTPUT_CHECK_V1",
        "status": status,
        "checks_passed": sum(row["status"] == "PASS" for row in checks),
        "checks_total": len(checks),
        "checks": checks,
    }
    VALIDATION_OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
