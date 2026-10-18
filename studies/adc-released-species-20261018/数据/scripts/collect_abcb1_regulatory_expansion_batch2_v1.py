#!/usr/bin/env python3
"""Collect a second, evidence-tiered ABCB1 regulatory expansion batch.

This batch prioritizes non-substrate evidence and keeps strict matched-parental
rows separate from quantitative transfectant-only, qualitative regulatory, and
Caco-2 evidence. It standardizes structures and audits identity overlap, but it
does not score or refit any model.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from abcb1_screener_core_v0_1 import standardize_smiles  # noqa: E402
from collect_abcb1_regulatory_reuse_panel_v1 import (  # noqa: E402
    fetch_pubchem,
    identity_index,
    read_csv,
    sha256,
    write_csv,
)


PROCESSED = PROJECT / "数据/processed"
RAW_DIR = PROJECT / "数据/raw/abcb1_regulatory_expansion_batch2"

PRIOR_DATASETS = {
    "frozen_443_training": (PROCESSED / "abcb1_screener_v0_1_training_audit.csv", ["standardized_inchikey"]),
    "soskuti2024_development": (PROCESSED / "soskuti_2024_mdr1_species_panel_v1.csv", ["standardized_inchikey"]),
    "q1_external_development": (
        PROCESSED / "abcb1_q1_direct_external_entity_predictions_v1.csv",
        ["model_standardized_inchikey", "standardized_inchikey"],
    ),
    "kelly2022_postfreeze_stress": (
        PROCESSED / "abcb1_kelly2022_table7_external_stress_v1.csv",
        ["standardized_inchikey"],
    ),
    "broccatelli2012_candidate_pool": (
        PROCESSED / "broccatelli_2012_abcb1_entities_v1.csv",
        ["standardized_inchikey"],
    ),
    "regulatory_reuse_v1": (
        PROCESSED / "abcb1_regulatory_reuse_entities_v1.csv",
        ["standardized_inchikey"],
    ),
}


SOURCES = {
    "PMDA_LORLATINIB": {
        "url": "https://www.pmda.go.jp/files/000228981.pdf",
        "local_pdf": "lorlatinib_pmda_000228981.pdf",
        "local_text": "lorlatinib_pmda_000228981.txt",
        "location": "Section 4.5.3 Transporters; PDF page 24 / printed review page 21; text lines 1228-1237",
    },
    "PMDA_GLYCOPYRRONIUM": {
        "url": "https://www.pmda.go.jp/files/000240532.pdf",
        "local_pdf": "glycopyrronium_pmda_000240532.pdf",
        "local_text": "glycopyrronium_pmda_000240532.txt",
        "location": "Non-clinical pharmacokinetics, transporter paragraph; PDF page 21 / printed review page 17; text lines 983-984",
    },
    "PMDA_DELAMANID": {
        "url": "https://www.pmda.go.jp/files/000209024.pdf",
        "local_pdf": "delamanid_pmda_000209024.pdf",
        "local_text": "delamanid_pmda_000209024.txt",
        "location": "Section 3(ii).A.(5).3; PDF pages 26-27 / printed pages 26-27; text lines 1582-1593",
    },
    "PMDA_SOFOSBUVIR": {
        "url": "https://www.pmda.go.jp/files/000210800.pdf",
        "local_pdf": "sofosbuvir_pmda_000210800.pdf",
        "local_text": "sofosbuvir_pmda_000210800.txt",
        "location": "Section 3(ii).A.(5).2; PDF page 24 / printed review page 24; text lines 1239-1247",
    },
    "PMDA_TENAPANOR": {
        "url": "https://www.pmda.go.jp/files/000271700.pdf",
        "local_pdf": "tenapanor_pmda_000271700.pdf",
        "local_text": "tenapanor_pmda_000271700.txt",
        "location": "Section 4.5.3; PDF pages 21-22 / printed review pages 18-19; text lines 1054-1085",
    },
    "PMDA_ACENEURAMIC_ACID": {
        "url": "https://www.pmda.go.jp/files/000273006.pdf",
        "local_pdf": "aceneuramic_acid_pmda_000273006.pdf",
        "local_text": "aceneuramic_acid_pmda_000273006.txt",
        "location": "Section (d) Transport by drug transporters; PDF page 22 / printed review page 19; text lines 1137-1144",
    },
    "PMDA_FORODESINE": {
        "url": "https://www.pmda.go.jp/files/000226817.pdf",
        "local_pdf": "forodesine_pmda_000226817.pdf",
        "local_text": "forodesine_pmda_000226817.txt",
        "location": "Section 4.5.3 Transporters; PDF page 17 / printed review page 14; text lines 813-822",
    },
    "FDA_4_AMINOPYRIDINE": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2010/022250Orig1s000ClinPharmR.pdf",
        "local_pdf": "fampridine_fda_022250_clinpharm.pdf",
        "local_text": "fampridine_fda_022250_clinpharm.txt",
        "location": "Study 8ACORP1; PDF pages 143-144 / NDA review pages 149-150; text lines 4851-4892",
    },
}


ENTITIES = [
    {
        "entity_name": "Lorlatinib",
        "pubchem_query_names": ["Lorlatinib"],
        "source_key": "PMDA_LORLATINIB",
        "label": 0,
        "assay_system": "human P-gp-expressing MDCKII; matched parent values not exposed",
        "concentration_uM": 2,
        "er_transporter": 1.8,
        "er_parental": "",
        "net_er": "",
        "inhibitor": "not reported for substrate assessment",
        "inhibitor_er": "",
        "positive_control": "not exposed",
        "source_call": "not a P-gp substrate",
        "raw_evidence_tier": "A2_quantitative_human_transfectant_no_parent",
        "strict_matched_parent_candidate": "no",
        "qc_note": "numeric ER below 2 and explicit regulatory call; no parent or substrate-assay inhibitor row",
    },
    {
        "entity_name": "Glycopyrronium",
        "pubchem_query_names": ["Glycopyrronium", "Glycopyrrolate"],
        "source_key": "PMDA_GLYCOPYRRONIUM",
        "label": 0,
        "assay_system": "human P-gp-expressing MDCK; parent values not exposed",
        "concentration_uM": "0.1-100",
        "er_transporter": "",
        "er_parental": "",
        "net_er": "",
        "inhibitor": "not exposed",
        "inhibitor_er": "",
        "positive_control": "not exposed",
        "source_call": "not a P-gp substrate",
        "raw_evidence_tier": "B1_official_directional_call_numeric_range_only",
        "strict_matched_parent_candidate": "no",
        "qc_note": "wide concentration range but no ER, matched parent, or inhibitor values",
    },
    {
        "entity_name": "Delamanid",
        "pubchem_query_names": ["Delamanid"],
        "source_key": "PMDA_DELAMANID",
        "label": 0,
        "assay_system": "human MDR1-expressing LLC-PK1; parent values not exposed",
        "concentration_uM": 5,
        "er_transporter": "",
        "er_parental": "",
        "net_er": "",
        "inhibitor": "not exposed",
        "inhibitor_er": "",
        "positive_control": "not exposed",
        "source_call": "not a P-gp substrate",
        "raw_evidence_tier": "B1_official_directional_call_concentration_only",
        "strict_matched_parent_candidate": "no",
        "qc_note": "human MDR1 system and explicit call, but no quantitative endpoint is disclosed",
    },
    {
        "entity_name": "Sofosbuvir",
        "pubchem_query_names": ["Sofosbuvir"],
        "source_key": "PMDA_SOFOSBUVIR",
        "label": 1,
        "assay_system": "Caco-2 endogenous P-gp context",
        "concentration_uM": "not exposed in substrate paragraph",
        "er_transporter": 12.4,
        "er_parental": "",
        "net_er": "",
        "inhibitor": "verapamil or cyclosporine A",
        "inhibitor_er": "3.66 or 4.66",
        "positive_control": "not exposed",
        "source_call": "P-gp substrate",
        "raw_evidence_tier": "B2_quantitative_Caco2_inhibitor_rescue",
        "strict_matched_parent_candidate": "no",
        "qc_note": "strong inhibitor reduction but Caco-2 is a multi-transporter system without matched parent",
    },
    {
        "entity_name": "GS-331007",
        "pubchem_query_names": ["GS-331007", "2'-Deoxy-2'-fluoro-2'-C-methyluridine"],
        "source_key": "PMDA_SOFOSBUVIR",
        "label": 0,
        "assay_system": "Caco-2 endogenous P-gp context",
        "concentration_uM": "not exposed in substrate paragraph",
        "er_transporter": "not exposed",
        "er_parental": "",
        "net_er": "",
        "inhibitor": "human P-gp inhibitor; identity not specified in sentence",
        "inhibitor_er": "not altered",
        "positive_control": "paired with inhibitor-responsive sofosbuvir",
        "source_call": "not a P-gp substrate",
        "raw_evidence_tier": "B2_paired_Caco2_inhibitor_nonresponse",
        "strict_matched_parent_candidate": "no",
        "qc_note": "useful parent-metabolite contrast; exact baseline ER is not disclosed",
    },
    {
        "entity_name": "Tenapanor",
        "pubchem_query_names": ["Tenapanor"],
        "source_key": "PMDA_TENAPANOR",
        "label": 0,
        "assay_system": "human MDR1-MDCKII versus untransfected MDCKII",
        "concentration_uM": 10,
        "er_transporter": ">10.7",
        "er_parental": 14.5,
        "net_er": "background-corrected regulatory negative; exact upper ER not exposed",
        "inhibitor": "zosuquidar 10 uM",
        "inhibitor_er": "1.91 in MDR1 / 2.39 in parent",
        "positive_control": "parallel P-gp probe-substrate system documented in same section",
        "source_call": "not a P-gp substrate",
        "raw_evidence_tier": "A1_matched_parent_and_inhibitor_background_adjudication",
        "strict_matched_parent_candidate": "yes",
        "qc_note": "high raw ER is a cell-background artifact; matched parent and inhibitor changes are similar",
    },
    {
        "entity_name": "Aceneuramic acid",
        "pubchem_query_names": ["Aceneuramic acid", "N-Acetylneuraminic acid"],
        "source_key": "PMDA_ACENEURAMIC_ACID",
        "label": 0,
        "assay_system": "human P-gp-expressing MDCKII; parent values not exposed",
        "concentration_uM": 100,
        "er_transporter": 0.914,
        "er_parental": "",
        "net_er": "",
        "inhibitor": "valspodar",
        "inhibitor_er": 0.983,
        "positive_control": "not exposed",
        "source_call": "not a P-gp substrate",
        "raw_evidence_tier": "A2_quantitative_human_transfectant_inhibitor_nonresponse",
        "strict_matched_parent_candidate": "no",
        "qc_note": "numeric transfectant ER and inhibitor condition; no matched parent values",
    },
    {
        "entity_name": "Forodesine",
        "pubchem_query_names": ["Forodesine"],
        "source_key": "PMDA_FORODESINE",
        "label": 0,
        "assay_system": "Caco-2 endogenous P-gp context",
        "concentration_uM": 5,
        "er_transporter": 1.41,
        "er_parental": "",
        "net_er": "",
        "inhibitor": "verapamil 100 uM",
        "inhibitor_er": 1.71,
        "positive_control": "not exposed",
        "source_call": "not a P-gp substrate",
        "raw_evidence_tier": "B2_quantitative_Caco2_inhibitor_nonresponse",
        "strict_matched_parent_candidate": "no",
        "qc_note": "raw ER below 2 and no inhibitor reduction; multi-transporter Caco-2 context",
    },
    {
        "entity_name": "4-Aminopyridine",
        "pubchem_query_names": ["4-Aminopyridine", "Fampridine"],
        "source_key": "FDA_4_AMINOPYRIDINE",
        "label": 0,
        "assay_system": "human MDR1-MDCK versus MDCK parent",
        "concentration_uM": "0.7,7,70",
        "er_transporter": "1.0,0.8,1.0",
        "er_parental": "measured in parallel; exact table not repeated in summary",
        "net_er": "0.8-1.0",
        "inhibitor": "not required for substrate classification",
        "inhibitor_er": "",
        "positive_control": "digoxin ER 123 in same MDR1-MDCK cell batch",
        "source_call": "not a P-gp substrate",
        "raw_evidence_tier": "A1_matched_parent_multiconcentration_positive_control",
        "strict_matched_parent_candidate": "yes",
        "qc_note": "three concentrations, triplicates, recovery, matched parent, and strong positive control",
    },
]


def fetch_pubchem_with_fallback(names: list[str]) -> tuple[dict, str, str, str]:
    errors: list[str] = []
    for name in names:
        record, status, error = fetch_pubchem(name)
        if record:
            return record, status, error, name
        errors.append(f"{name}: {status}: {error}")
    return {}, "pubchem_all_queries_failed", " | ".join(errors), ""


def common_strict_row(row: dict, origin: str) -> dict:
    return {
        "strict_candidate_id": row.get("regulatory_entity_id") or row.get("batch2_entity_id"),
        "entity_name": row["entity_name"],
        "directional_label": row["label"],
        "collection_origin": origin,
        "evidence_tier": row.get("evidence_tier", ""),
        "assay_system": row.get("assay_system", ""),
        "concentration_uM": row.get("concentration_uM", ""),
        "er_transporter": row.get("er_transporter", ""),
        "er_parental": row.get("er_parental", ""),
        "net_er": row.get("net_er", ""),
        "inhibitor": row.get("inhibitor", ""),
        "inhibitor_er": row.get("inhibitor_er", ""),
        "positive_control": row.get("positive_control", ""),
        "standardized_inchikey": row.get("standardized_inchikey", ""),
        "standardized_isomeric_smiles": row.get("standardized_isomeric_smiles", ""),
        "source_id": row.get("source_id", ""),
        "source_url": row.get("source_url", ""),
        "source_location": row.get("source_location", ""),
        "identity_status": row.get("identity_status", ""),
        "claim_boundary": row.get("claim_boundary", ""),
    }


def main() -> None:
    prior_indexes = {
        name: identity_index(path, columns)
        for name, (path, columns) in PRIOR_DATASETS.items()
    }

    entity_rows: list[dict] = []
    for index, manual in enumerate(ENTITIES, start=1):
        source = SOURCES[manual["source_key"]]
        pubchem, pubchem_status, pubchem_error, pubchem_query_used = fetch_pubchem_with_fallback(
            manual["pubchem_query_names"]
        )
        standard = standardize_smiles(pubchem.get("pubchem_smiles", ""))
        key = standard.get("standardized_inchikey", "")
        block = key.split("-")[0] if key else ""
        exact_hits = [name for name, (exact, _) in prior_indexes.items() if key and key in exact]
        connectivity_hits = [name for name, (_, conn) in prior_indexes.items() if block and block in conn]
        if standard.get("standardization_status") != "ok":
            identity_status = "exclude_unresolved_structure"
        elif exact_hits:
            identity_status = "prior_exact_identity_overlap"
        elif connectivity_hits:
            identity_status = "prior_connectivity_identity_overlap"
        else:
            identity_status = "no_exact_or_connectivity_prior_overlap"

        strict = manual["strict_matched_parent_candidate"] == "yes"
        candidate_eligible = strict and identity_status == "no_exact_or_connectivity_prior_overlap"
        local_pdf = RAW_DIR / source["local_pdf"]
        local_text = RAW_DIR / source["local_text"]
        manual_without_queries = {k: v for k, v in manual.items() if k != "pubchem_query_names"}
        entity_rows.append(
            {
                "batch2_entity_id": f"ABCB1_REG_BATCH2_{index:03d}",
                **manual_without_queries,
                **pubchem,
                "pubchem_query_used": pubchem_query_used,
                "pubchem_resolution_status": pubchem_status,
                "pubchem_resolution_error": pubchem_error,
                **standard,
                "inchikey_connectivity_block": block,
                "overlap_exact_prior_sets": " | ".join(exact_hits),
                "overlap_connectivity_prior_sets": " | ".join(connectivity_hits),
                "identity_status": identity_status,
                "strict_external_candidate_eligible": "yes" if candidate_eligible else "no",
                "source_id": manual["source_key"],
                "source_url": source["url"],
                "source_location": source["location"],
                "local_pdf": str(local_pdf.relative_to(PROJECT)),
                "local_pdf_sha256": sha256(local_pdf) if local_pdf.exists() else "missing",
                "local_text": str(local_text.relative_to(PROJECT)),
                "evidence_tier": manual["raw_evidence_tier"],
                "claim_boundary": "external regulatory evidence candidate; no model scoring or formal validation claim before frozen source-lineage audit",
            }
        )

    write_csv(PROCESSED / "abcb1_regulatory_expansion_batch2_entities_v1.csv", entity_rows)

    source_rows: list[dict] = []
    for source_id, source in SOURCES.items():
        local_pdf = RAW_DIR / source["local_pdf"]
        source_entities = [row for row in entity_rows if row["source_id"] == source_id]
        source_rows.append(
            {
                "source_id": source_id,
                "source_url": source["url"],
                "source_location": source["location"],
                "local_pdf": str(local_pdf.relative_to(PROJECT)),
                "local_pdf_sha256": sha256(local_pdf) if local_pdf.exists() else "missing",
                "entity_count": len(source_entities),
                "a1_raw_count": sum(row["strict_matched_parent_candidate"] == "yes" for row in source_entities),
                "strict_identity_clean_count": sum(row["strict_external_candidate_eligible"] == "yes" for row in source_entities),
                "allowed_use": "ABCB1 evidence-tiered external candidate collection; scoring only after separate frozen protocol",
            }
        )
    write_csv(PROCESSED / "abcb1_regulatory_expansion_batch2_source_registry_v1.csv", source_rows)

    v1_rows = read_csv(PROCESSED / "abcb1_regulatory_reuse_entities_v1.csv")
    tiered_rows: list[dict] = []
    for row in v1_rows:
        tiered_rows.append({"entity_record_id": row["regulatory_entity_id"], "collection_batch": "regulatory_reuse_v1", **row})
    for row in entity_rows:
        tiered_rows.append({"entity_record_id": row["batch2_entity_id"], "collection_batch": "regulatory_expansion_batch2", **row})
    write_csv(PROCESSED / "abcb1_regulatory_evidence_tiered_pool_v2.csv", tiered_rows)

    v1_strict = [row for row in v1_rows if row.get("strict_external_candidate_eligible") == "yes"]
    batch2_strict = [row for row in entity_rows if row.get("strict_external_candidate_eligible") == "yes"]
    strict_pool = [common_strict_row(row, "regulatory_reuse_v1") for row in v1_strict]
    strict_pool += [common_strict_row(row, "regulatory_expansion_batch2") for row in batch2_strict]
    write_csv(PROCESSED / "abcb1_regulatory_strict_external_pool_v2.csv", strict_pool)

    candidate_ledger_v1 = read_csv(PROCESSED / "abcb1_external_collection_candidate_ledger_v1.csv")
    candidate_ledger_v2 = list(candidate_ledger_v1)
    for row in entity_rows:
        is_strict = row["strict_external_candidate_eligible"] == "yes"
        inhibitor_context = row.get("inhibitor", "")
        if row.get("inhibitor_er"):
            inhibitor_context = f"{inhibitor_context}; inhibitor ER {row['inhibitor_er']}"
        candidate_ledger_v2.append(
            {
                "collection_candidate_id": f"EXTCOLL::{row['batch2_entity_id']}",
                "entity_name": row["entity_name"],
                "directional_label": row["label"],
                "evidence_route": (
                    "regulatory_direct_human_transfectant_matched_control"
                    if is_strict
                    else "regulatory_evidence_tiered_secondary"
                ),
                "formal_role_now": (
                    "pre_scoring_strict_external_candidate"
                    if is_strict
                    else "regulatory_secondary_evidence_only"
                ),
                "quantitative_endpoint": row.get("net_er") or row.get("er_transporter") or row.get("source_call", ""),
                "endpoint_definition": (
                    f"{row.get('assay_system', '')}; source call: {row.get('source_call', '')}"
                ),
                "matched_parental_control": "yes" if is_strict else "no_or_not_numeric",
                "selective_inhibitor_rescue": inhibitor_context,
                "source_id": row["source_id"],
                "source_cluster_count": 1,
                "source_clusters": row["source_id"],
                "standardized_inchikey": row.get("standardized_inchikey", ""),
                "inchikey_connectivity_block": row.get("inchikey_connectivity_block", ""),
                "standardized_isomeric_smiles": row.get("standardized_isomeric_smiles", ""),
                "source_url": row.get("source_url", ""),
                "source_location": row.get("source_location", ""),
                "evidence_tier": row.get("evidence_tier", ""),
                "remaining_gate": (
                    "freeze entity list and score once; retain source-clustered uncertainty"
                    if is_strict
                    else "retain outside formal external metric; use only for tier-stratified sensitivity analysis"
                ),
                "claim_boundary": row.get("claim_boundary", ""),
            }
        )
    write_csv(PROCESSED / "abcb1_external_collection_candidate_ledger_v2.csv", candidate_ledger_v2)

    queue_v1 = read_csv(PROCESSED / "abcb1_regulatory_expansion_search_queue_v1.csv")
    extracted_names = {
        "Lorlatinib",
        "Glycopyrronium (GP)",
        "Delamanid",
        "Sofosbuvir; GS-331007",
    }
    queue_v2 = []
    for row in queue_v1:
        new_row = dict(row)
        if row["candidate_names"] in extracted_names:
            new_row["status"] = "extracted_and_tiered_batch2"
        queue_v2.append(new_row)
    queue_v2.extend(
        [
            {
                "priority": 1,
                "candidate_names": "Tenapanor",
                "source_reported_class": "negative",
                "source_summary": "matched MDR1-MDCKII and untransfected MDCKII plus inhibitor conditions; background-driven high raw ER",
                "source_url": SOURCES["PMDA_TENAPANOR"]["url"],
                "route": "PMDA reverse search",
                "planned_role": "A1 strict negative candidate",
                "status": "extracted_and_tiered_batch2",
            },
            {
                "priority": 2,
                "candidate_names": "Aceneuramic acid",
                "source_reported_class": "negative",
                "source_summary": "human P-gp MDCKII ER 0.914 and valspodar ER 0.983",
                "source_url": SOURCES["PMDA_ACENEURAMIC_ACID"]["url"],
                "route": "PMDA reverse search",
                "planned_role": "A2 quantitative negative",
                "status": "extracted_and_tiered_batch2",
            },
            {
                "priority": 2,
                "candidate_names": "Forodesine",
                "source_reported_class": "negative",
                "source_summary": "Caco-2 ER 1.41 without and 1.71 with verapamil",
                "source_url": SOURCES["PMDA_FORODESINE"]["url"],
                "route": "PMDA reverse search",
                "planned_role": "B2 supportive negative",
                "status": "extracted_and_tiered_batch2",
            },
            {
                "priority": 1,
                "candidate_names": "4-Aminopyridine",
                "source_reported_class": "negative",
                "source_summary": "MDR1-MDCK versus MDCK net ER 0.8-1.0 at three concentrations; digoxin ER 123 positive control",
                "source_url": SOURCES["FDA_4_AMINOPYRIDINE"]["url"],
                "route": "FDA reverse search",
                "planned_role": "A1 strict negative candidate",
                "status": "extracted_and_tiered_batch2",
            },
        ]
    )
    write_csv(PROCESSED / "abcb1_regulatory_expansion_search_queue_v2.csv", queue_v2)

    class_counts = Counter(str(row["directional_label"]) for row in strict_pool)
    tier_counts = Counter(row["evidence_tier"] for row in entity_rows)
    summary = {
        "batch2_source_documents": len(source_rows),
        "batch2_entity_rows": len(entity_rows),
        "batch2_standardization_ok": sum(row["standardization_status"] == "ok" for row in entity_rows),
        "batch2_identity_clean_a1_candidates": len(batch2_strict),
        "batch2_identity_clean_a1_names": [row["entity_name"] for row in batch2_strict],
        "batch2_tier_counts": dict(tier_counts),
        "regulatory_strict_pool_v2_rows": len(strict_pool),
        "regulatory_strict_pool_v2_class_counts": dict(class_counts),
        "external_candidate_ledger_v2_rows": len(candidate_ledger_v2),
        "external_candidate_ledger_v2_role_counts": dict(
            Counter(row["formal_role_now"] for row in candidate_ledger_v2)
        ),
        "formal_external_shortfall_if_minimum_15_per_class": {
            "positive": max(0, 15 - class_counts.get("1", 0)),
            "negative": max(0, 15 - class_counts.get("0", 0)),
        },
        "model_scoring_performed": False,
        "formal_status": "evidence_tiered_candidate_collection_complete_pending_source_lineage_lock_and_scoring_protocol",
    }
    (PROCESSED / "abcb1_regulatory_expansion_batch2_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
