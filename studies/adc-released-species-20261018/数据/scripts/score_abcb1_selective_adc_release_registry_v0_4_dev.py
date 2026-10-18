#!/usr/bin/env python3
"""Apply v0.4-dev to the project's explicit ADC release-species registry."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path


DATA = Path(__file__).resolve().parents[1]
PROJECT = DATA.parent
PROCESSED = DATA / "processed"
MODEL_DIR = PROJECT / "模型" / "ABCB1_selective_screener_v0_4_dev"
sys.path.insert(0, str(MODEL_DIR))

from predict_abcb1_selective import flatten_batch_result, load_packages, predict_one  # noqa: E402


REGISTRY = PROCESSED / "adc_release_species_registry_standardized_v1.csv"
SELECTIVE = MODEL_DIR / "ABCB1_selective_uncertainty_v0_4_dev.joblib"
PRIMARY = PROJECT / "模型" / "ABCB1_screener_v0_1" / "ABCB1_screener_v0_1.joblib"
DIRECT = PROJECT / "模型" / "ABCB1_assay_aware_candidate_v0_3_dev" / "ABCB1_assay_aware_candidate_v0_3_dev.joblib"
OUTPUT = PROCESSED / "abcb1_selective_v0_4_adc_release_registry_predictions_v1.csv"
SUMMARY = PROCESSED / "abcb1_selective_v0_4_adc_release_registry_summary_v1.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    selective, primary, direct = load_packages(SELECTIVE, PRIMARY, DIRECT)
    rows = read_csv(REGISTRY)
    predictions: list[dict] = []
    for row in rows:
        smiles = row.get("standardized_isomeric_smiles", "")
        result = flatten_batch_result(
            predict_one(
                smiles,
                selective,
                primary,
                direct,
                chemical_object_type=row.get("reference_object_type", "release_species"),
                requested_context="structure_only_screening",
            )
        )
        predictions.append(
            {
                **row,
                **result,
                "application_role": "ADC_release_species_research_prioritization_not_validation",
            }
        )
    write_csv(OUTPUT, predictions)

    calls = Counter(row["risk_call"] for row in predictions)
    reasons = Counter(
        reason
        for row in predictions
        for reason in json.loads(row["abstention_reasons_json"])
    )
    summary = {
        "model_version": selective["model_version"],
        "application_status": "scored_explicit_ADC_release_species_registry_not_validation",
        "registry_path": str(REGISTRY.relative_to(PROJECT)),
        "registry_sha256": sha256_file(REGISTRY),
        "model_package_sha256": sha256_file(SELECTIVE),
        "registry_rows": len(rows),
        "explicit_structure_rows": sum(bool(row.get("standardized_isomeric_smiles", "")) for row in rows),
        "risk_call_counts": dict(calls),
        "abstention_reason_counts": dict(reasons),
        "output": str(OUTPUT.relative_to(PROJECT)),
        "interpretation": (
            "These calls demonstrate the intended ADC release-species interface and prioritize review or experiments. "
            "They are not accuracy estimates because the registry is not an ABCB1 truth set."
        ),
        "claim_boundary": selective["claim_boundary"],
    }
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
