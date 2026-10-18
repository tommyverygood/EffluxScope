#!/usr/bin/env python3
"""Prepare frozen-source inputs for ABCB1 quantity-oriented external validation."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "数据" / "raw" / "abcb1_quantity_external_validation_20260805"
PROCESSED = ROOT / "数据" / "processed"

BIOGEN_SOURCE = RAW_DIR / "ADME_public_set_3521.csv"
DAOOD_SOURCE = PROCESSED / "daood_2025_pgp_substrate_external_candidates_v1.csv"
BIOGEN_OUTPUT = PROCESSED / "abcb1_quantity_biogen_model_input_v1.csv"
DAOOD_OUTPUT = PROCESSED / "abcb1_quantity_daood_model_input_v1.csv"
SOURCE_REGISTRY_OUTPUT = PROCESSED / "abcb1_quantity_external_source_registry_v1.csv"
SUMMARY_OUTPUT = PROCESSED / "abcb1_quantity_external_preparation_summary_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def prepare_biogen() -> tuple[list[dict[str, object]], dict[str, object]]:
    rows = read_csv(BIOGEN_SOURCE)
    endpoint = "LOG MDR1-MDCK ER (B-A/A-B)"
    output: list[dict[str, object]] = []
    for row_index, row in enumerate(rows, start=1):
        raw_value = str(row.get(endpoint, "")).strip()
        smiles = str(row.get("SMILES", "")).strip()
        if not raw_value or not smiles:
            continue
        try:
            log10_er = float(raw_value)
        except ValueError:
            continue
        if not math.isfinite(log10_er):
            continue
        output.append(
            {
                "external_row_id": f"BIOGEN_MDR1_{row_index:04d}",
                "source_id": "biogen_fang_2023_mdr1_mdck",
                "source_record_id": row.get("Internal ID", ""),
                "vendor_id": row.get("Vendor ID", ""),
                "collection_name": row.get("CollectionName", ""),
                "smiles": smiles,
                "log10_mdr1_mdck_er": f"{log10_er:.12g}",
                "mdr1_mdck_er": f"{10 ** log10_er:.12g}",
                "frozen_er2_binary_label": int(log10_er >= math.log10(2.0)),
                "chemical_object_type": "released_small_molecule",
                "requested_context": "structure_only_screening",
                "source_doi": "10.1021/acs.jcim.3c00160",
                "source_url": "https://github.com/molecularinformatics/Computational-ADME/blob/main/ADME_public_set_3521.csv",
                "endpoint_claim_boundary": "continuous MDR1-MDCK efflux-risk benchmark; no matched-parental or inhibitor-rescue truth",
            }
        )
    summary = {
        "source_rows": len(rows),
        "eligible_endpoint_rows": len(output),
        "frozen_er2_label_counts": dict(Counter(int(row["frozen_er2_binary_label"]) for row in output)),
    }
    return output, summary


def prepare_daood() -> tuple[list[dict[str, object]], dict[str, object]]:
    rows = read_csv(DAOOD_SOURCE)
    output: list[dict[str, object]] = []
    for row_index, row in enumerate(rows, start=1):
        label_text = str(row.get("author_pgp_substrate_label", "")).strip()
        smiles = str(row.get("smiles", "")).strip()
        if label_text not in {"0", "1"} or not smiles:
            continue
        output.append(
            {
                "external_row_id": f"DAOOD_PGP_{row_index:04d}",
                "source_id": "daood_2025_curated_pgp_substrate",
                "source_record_id": row.get("source_record_id", ""),
                "source_database": row.get("source_database", ""),
                "preferred_name": row.get("preferred_name", ""),
                "smiles": smiles,
                "author_pgp_substrate_label": int(label_text),
                "chemical_object_type": "released_small_molecule",
                "requested_context": "structure_only_screening",
                "source_doi": row.get("source_doi", "10.1021/acs.molpharmaceut.5c01065"),
                "source_url": row.get("source_url", "https://pmc.ncbi.nlm.nih.gov/articles/PMC12587445/"),
                "endpoint_claim_boundary": "heterogeneous author-curated P-gp substrate class; retrospective corroboration only",
            }
        )
    summary = {
        "source_rows": len(rows),
        "eligible_binary_rows": len(output),
        "author_label_counts": dict(Counter(int(row["author_pgp_substrate_label"]) for row in output)),
        "source_database_counts": dict(Counter(str(row["source_database"]) for row in output)),
    }
    return output, summary


def main() -> None:
    biogen_rows, biogen_summary = prepare_biogen()
    daood_rows, daood_summary = prepare_daood()
    write_csv(BIOGEN_OUTPUT, biogen_rows)
    write_csv(DAOOD_OUTPUT, daood_rows)

    source_registry = [
        {
            "source_id": "biogen_fang_2023_mdr1_mdck",
            "source_file": str(BIOGEN_SOURCE.relative_to(ROOT)),
            "sha256": sha256(BIOGEN_SOURCE),
            "eligible_rows": len(biogen_rows),
            "endpoint": "log10 MDR1-MDCK efflux ratio",
            "evidence_role": "primary_large_continuous_external_risk_benchmark",
            "source_doi": "10.1021/acs.jcim.3c00160",
            "source_url": "https://github.com/molecularinformatics/Computational-ADME/blob/main/ADME_public_set_3521.csv",
        },
        {
            "source_id": "daood_2025_curated_pgp_substrate",
            "source_file": str(DAOOD_SOURCE.relative_to(ROOT)),
            "sha256": sha256(DAOOD_SOURCE),
            "eligible_rows": len(daood_rows),
            "endpoint": "author-curated P-gp substrate binary class",
            "evidence_role": "secondary_large_heterogeneous_binary_corroboration",
            "source_doi": "10.1021/acs.molpharmaceut.5c01065",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12587445/",
        },
    ]
    write_csv(SOURCE_REGISTRY_OUTPUT, source_registry)

    summary = {
        "status": "prepared_after_protocol_lock_before_model_scoring",
        "protocol_path": "数据/processed/abcb1_quantity_external_validation_protocol_lock_v1.json",
        "biogen": biogen_summary,
        "daood": daood_summary,
        "outputs": {
            "biogen_model_input": str(BIOGEN_OUTPUT.relative_to(ROOT)),
            "daood_model_input": str(DAOOD_OUTPUT.relative_to(ROOT)),
            "source_registry": str(SOURCE_REGISTRY_OUTPUT.relative_to(ROOT)),
        },
    }
    SUMMARY_OUTPUT.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
