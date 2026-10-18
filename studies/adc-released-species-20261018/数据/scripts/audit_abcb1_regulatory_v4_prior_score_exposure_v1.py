#!/usr/bin/env python3
"""Audit whether strict-pool v4 structures were scored before its blind run.

The audit searches every prior prediction CSV for exact InChIKey or connectivity
matches.  Exact normalized-name matches are retained as weaker alerts when the
historical file does not expose a structure.  It does not run the model.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT = SCRIPT_DIR.parents[1]
PROCESSED = PROJECT / "数据/processed"
sys.path.insert(0, str(SCRIPT_DIR))

from abcb1_screener_core_v0_1 import standardize_smiles  # noqa: E402


POOL = PROCESSED / "abcb1_regulatory_strict_external_pool_v4.csv"
OUT_CSV = PROCESSED / "abcb1_regulatory_v4_prior_score_exposure_audit_v1.csv"
OUT_JSON = PROCESSED / "abcb1_regulatory_v4_prior_score_exposure_summary_v1.json"

KEY_FIELDS = [
    "model_standardized_inchikey",
    "standardized_inchikey",
    "derived_standardized_inchikey",
    "inchikey",
]
SMILES_FIELDS = [
    "model_standardized_isomeric_smiles",
    "standardized_isomeric_smiles",
    "standardized_smiles",
    "input_smiles",
    "source_smiles",
    "smiles",
]
NAME_FIELDS = ["entity_name", "compound_name", "preferred_name", "name"]
SCORE_FIELDS = [
    "frozen_probability",
    "calibrated_probability",
    "predicted_probability",
    "raw_probability",
    "probability",
    "binary_prediction",
    "predicted_label",
]


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0]) if rows else [
        "strict_candidate_id",
        "entity_name",
        "historical_file",
        "historical_row_index",
        "match_type",
        "historical_name",
        "historical_inchikey",
        "historical_connectivity_block",
        "score_fields_exposed",
        "score_values_exposed",
        "adjudication",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def normalize_name(value: str) -> str:
    text = re.sub(r"[^a-z0-9]+", " ", (value or "").lower()).strip()
    text = re.sub(r"\b(db75|cpd0801|ro0*681133|ro0*731519|m1|m3)\b", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def historical_identity(row: dict) -> tuple[str, str, str]:
    key = next((row.get(field, "").strip() for field in KEY_FIELDS if row.get(field, "").strip()), "")
    if key:
        return key, key.split("-")[0], "reported_inchikey"
    for field in SMILES_FIELDS:
        smiles = row.get(field, "").strip()
        if not smiles:
            continue
        standardized = standardize_smiles(smiles)
        key = standardized.get("standardized_inchikey", "")
        if key:
            return key, key.split("-")[0], f"derived_from_{field}"
    return "", "", "unresolved"


def main() -> None:
    pool = read_csv(POOL)
    by_exact = {row["standardized_inchikey"]: row for row in pool if row.get("standardized_inchikey")}
    by_connectivity = {row["standardized_inchikey"].split("-")[0]: row for row in pool if row.get("standardized_inchikey")}
    by_name = {normalize_name(row["entity_name"]): row for row in pool if normalize_name(row["entity_name"])}

    matches = []
    files = sorted({*PROCESSED.glob("*prediction*.csv"), *PROCESSED.glob("*predictions*.csv")})
    for path in files:
        for row_index, row in enumerate(read_csv(path), start=2):
            transporter = (row.get("transporter_gene", "") or "").strip()
            if transporter and transporter != "ABCB1":
                continue
            if not transporter and "abcb1" not in path.name.lower() and "pgp" not in path.name.lower():
                continue
            key, block, identity_route = historical_identity(row)
            candidate = by_exact.get(key) or by_connectivity.get(block)
            match_type = ""
            if candidate is not None:
                match_type = "exact_inchikey" if key == candidate["standardized_inchikey"] else "connectivity_block"
            else:
                historical_name = next((row.get(field, "") for field in NAME_FIELDS if row.get(field, "")), "")
                candidate = by_name.get(normalize_name(historical_name))
                if candidate is not None:
                    match_type = "normalized_name_only"
            if candidate is None:
                continue
            exposed = [(field, row.get(field, "")) for field in SCORE_FIELDS if row.get(field, "") not in {"", None}]
            historical_name = next((row.get(field, "") for field in NAME_FIELDS if row.get(field, "")), "")
            matches.append(
                {
                    "strict_candidate_id": candidate["strict_candidate_id"],
                    "entity_name": candidate["entity_name"],
                    "historical_file": str(path.relative_to(PROJECT)),
                    "historical_row_index": row_index,
                    "match_type": match_type,
                    "historical_name": historical_name,
                    "historical_inchikey": key,
                    "historical_connectivity_block": block,
                    "identity_route": identity_route,
                    "score_fields_exposed": " | ".join(field for field, _ in exposed),
                    "score_values_exposed": " | ".join(str(value) for _, value in exposed),
                    "adjudication": (
                        "prior_score_exposure_confirmed"
                        if match_type in {"exact_inchikey", "connectivity_block"} and exposed
                        else "name_only_alert_manual_review"
                    ),
                }
            )

    write_csv(OUT_CSV, matches)
    confirmed_ids = {
        row["strict_candidate_id"]
        for row in matches
        if row["adjudication"] == "prior_score_exposure_confirmed"
    }
    name_only_ids = {
        row["strict_candidate_id"]
        for row in matches
        if row["adjudication"] == "name_only_alert_manual_review"
    }
    summary = {
        "pool_path": str(POOL.relative_to(PROJECT)),
        "historical_prediction_files_scanned": len(files),
        "strict_pool_entities": len(pool),
        "match_rows": len(matches),
        "confirmed_prior_score_exposure_entities": len(confirmed_ids),
        "confirmed_prior_score_exposure_ids": sorted(confirmed_ids),
        "name_only_alert_entities": len(name_only_ids),
        "name_only_alert_ids": sorted(name_only_ids),
        "primary_blind_eligible_entities": len(pool) - len(confirmed_ids),
        "adjudication_counts": dict(Counter(row["adjudication"] for row in matches)),
        "model_scoring_performed_by_this_audit": False,
    }
    OUT_JSON.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
