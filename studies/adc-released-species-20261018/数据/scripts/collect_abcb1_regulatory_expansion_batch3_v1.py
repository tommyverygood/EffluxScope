#!/usr/bin/env python3
"""Collect the predeclared third ABCB1 regulatory expansion batch.

The batch completes official PDF extraction, structure standardization, identity
overlap auditing, evidence-tier adjudication, and source-lineage registration.
It intentionally performs no model scoring or model/threshold modification.
"""

from __future__ import annotations

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
RAW_DIR = PROJECT / "数据/raw/abcb1_regulatory_expansion_batch3"

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
    "regulatory_evidence_pool_v2": (
        PROCESSED / "abcb1_regulatory_evidence_tiered_pool_v2.csv",
        ["standardized_inchikey"],
    ),
    "external_candidate_ledger_v2": (
        PROCESSED / "abcb1_external_collection_candidate_ledger_v2.csv",
        ["standardized_inchikey"],
    ),
}


SOURCES = {
    "FDA_ATRASENTAN": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2025/219208Orig1s000IntegratedR.pdf",
        "local_pdf": "atrasentan_fda_219208_integrated.pdf",
        "local_text": "atrasentan_fda_219208_integrated.txt",
        "location": "Section Atrasentan as Substrate of Gastrointestinal Efflux Transporters; PDF page 138 / printed page 125; text lines 6207-6220",
        "study_ids": "NC-24-0001 | Drug Metabolism Memo No. 32",
        "lineage_unit": "FDA_NDA219208_NC240001_ABCB1",
    },
    "PMDA_RIOCIGUAT": {
        "url": "https://www.pmda.go.jp/files/000207309.pdf",
        "local_pdf": "riociguat_pmda_000207309.pdf",
        "local_text": "riociguat_pmda_000207309.txt",
        "location": "Section 4.(ii).A.(1).4(b); PDF and printed page 38; text lines 2067-2073",
        "study_ids": "Attached document 4.2.2.7.2",
        "lineage_unit": "PMDA_RIOCIGUAT_CTD42272_ABCB1",
    },
    "PMDA_FUTIBATINIB": {
        "url": "https://www.pmda.go.jp/files/000267944.pdf",
        "local_pdf": "futibatinib_pmda_000267944.pdf",
        "local_text": "futibatinib_pmda_000267944.txt",
        "location": "Section 4.5.3 Transporters; PDF page 17 / printed page 14; text lines 809-817",
        "study_ids": "CTD transporter study summarized in Section 4.5.3; identifier not exposed",
        "lineage_unit": "PMDA_FUTIBATINIB_CTD453_ABCB1",
    },
    "PMDA_PIRTOBRUTINIB": {
        "url": "https://www.pmda.go.jp/files/000275223.pdf",
        "local_pdf": "pirtobrutinib_pmda_000275223.pdf",
        "local_text": "pirtobrutinib_pmda_000275223.txt",
        "location": "Section 4.5.3 Transporters; PDF page 21 / printed page 18; text lines 1035-1044",
        "study_ids": "CTD transporter study summarized in Section 4.5.3; identifier not exposed",
        "lineage_unit": "PMDA_PIRTOBRUTINIB_CTD453_ABCB1",
    },
    "FDA_GEPOTIDACIN": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2025/218230Orig1s000IntegratedR.pdf",
        "local_pdf": "gepotidacin_fda_218230_integrated.pdf",
        "local_text": "gepotidacin_fda_218230_integrated.txt",
        "location": "Section 14.1.4 and Table 110; PDF page 152 / printed page 150; text lines 7417-7483",
        "study_ids": "2015N235182 | 2020N436056 | 2015N235181 | 2013N159502 | 2013N158519 | 2011N124110",
        "lineage_unit": "FDA_NDA218230_GEPOTIDACIN_TRANSPORTER_REPORT_CLUSTER",
    },
    "FDA_BRENSOCATIB": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2025/217673Orig1s000IntegratedR.pdf",
        "local_pdf": "brensocatib_fda_217673_integrated.pdf",
        "local_text": "brensocatib_fda_217673_integrated.txt",
        "location": "In vitro DDI transporter summary table; PDF page 144 / printed page 126; text lines 6763-6783",
        "study_ids": "Insmed-01 | BJAA-0002-DV-PB",
        "lineage_unit": "FDA_NDA217673_INSMED01_BJAA0002_ABCB1",
    },
}


ENTITIES = [
    {
        "entity_name": "Atrasentan",
        "source_key": "FDA_ATRASENTAN",
        "label": 1,
        "assay_system": "human MDR1-transfected MDCKII versus wild-type; confirmatory Abcb1KO-MDCKII-MDR1 versus Mock-LV",
        "concentration_uM": "1",
        "er_transporter": "not separately exposed",
        "er_parental": "measured in net-ER denominator",
        "net_er": "60.7",
        "inhibitor": "LY335979 5 uM",
        "inhibitor_er": "1.41",
        "positive_control": "independent confirmatory ABCB1 study in same regulatory review",
        "source_call": "ABCB1/P-gp substrate",
        "evidence_tier": "A1_matched_background_inhibitor_confirmed",
        "strict_matched_parent_candidate": "yes",
        "qc_note": "knockout background and inhibitor collapse; BCRP evidence is separately adjudicated and not used for ABCB1 label",
    },
    {
        "entity_name": "Riociguat",
        "source_key": "PMDA_RIOCIGUAT",
        "label": 1,
        "assay_system": "human P-gp-expressing LLC-PK1 versus non-expressing LLC-PK1",
        "concentration_uM": "1",
        "er_transporter": "4.3-fold relative to parent",
        "er_parental": "matched non-expressing LLC-PK1",
        "net_er": "4.3",
        "inhibitor": "ivermectin 5 uM",
        "inhibitor_er": "decreased to parent-cell level",
        "positive_control": "dipyridamole net ratio 9.5",
        "source_call": "ABCB1/P-gp substrate",
        "evidence_tier": "A1_matched_parent_inhibitor_positive_control",
        "strict_matched_parent_candidate": "yes",
        "qc_note": "matched parent, inhibitor reversal, and same-page positive control",
    },
    {
        "entity_name": "Futibatinib",
        "source_key": "PMDA_FUTIBATINIB",
        "label": 1,
        "assay_system": "human P-gp-expressing LLC-PK1 versus non-expressing LLC-PK1",
        "concentration_uM": "0.3-30",
        "er_transporter": "relative efflux-ratio range 2.73-6.85",
        "er_parental": "matched non-expressing LLC-PK1 denominator",
        "net_er": "2.73-6.85",
        "inhibitor": "not exposed for substrate assessment",
        "inhibitor_er": "",
        "positive_control": "not exposed",
        "source_call": "ABCB1/P-gp substrate",
        "evidence_tier": "A1_matched_parent_multiconcentration_net_er",
        "strict_matched_parent_candidate": "yes",
        "qc_note": "all disclosed parent-corrected ratios exceed 2; inhibitor not required by locked rule",
    },
    {
        "entity_name": "Pirtobrutinib",
        "source_key": "PMDA_PIRTOBRUTINIB",
        "label": 1,
        "assay_system": "human P-gp-expressing MDCK versus non-expressing MDCK",
        "concentration_uM": "1-10",
        "er_transporter": "relative efflux-ratio range 17.6-22.3",
        "er_parental": "matched non-expressing MDCK denominator",
        "net_er": "17.6-22.3",
        "inhibitor": "elacridar 3 uM",
        "inhibitor_er": "1.17",
        "positive_control": "not exposed",
        "source_call": "ABCB1/P-gp substrate",
        "evidence_tier": "A1_matched_parent_inhibitor_collapse",
        "strict_matched_parent_candidate": "yes",
        "qc_note": "large matched-parent ratio and inhibitor collapse below 2",
    },
    {
        "entity_name": "Gepotidacin",
        "source_key": "FDA_GEPOTIDACIN",
        "label": 1,
        "assay_system": "human P-gp transfected bidirectional assay summarized across MDCK2 transporter study reports",
        "concentration_uM": "not separately exposed in Table 110",
        "er_transporter": "max flux/rate ratio 15.6",
        "er_parental": "net efflux denominator used; exact raw parent value not exposed",
        "net_er": "15.6",
        "inhibitor": "P-gp probe inhibitor; identity not exposed in summary",
        "inhibitor_er": ">50% reduction reported",
        "positive_control": "reviewer analysis against 2024 M12 criteria",
        "source_call": "ABCB1/P-gp substrate",
        "evidence_tier": "A1_net_er_inhibitor_reduction_reviewer_adjudicated",
        "strict_matched_parent_candidate": "yes",
        "qc_note": "review report explicitly states net efflux ratio >=2 and >50% probe-inhibitor reduction; study mapping retained as a cluster",
    },
    {
        "entity_name": "Brensocatib",
        "source_key": "FDA_BRENSOCATIB",
        "label": 1,
        "assay_system": "Abcb1KO-MDCKII-MDR1 qualitative result plus Caco-2 inhibitor-rescue experiment",
        "concentration_uM": "10 in Caco-2; transfectant concentration not exposed",
        "er_transporter": "3.86 in Caco-2; P-gp transfectant numeric endpoint not exposed",
        "er_parental": "Abcb1KO-MDCKII background used but numeric value not exposed",
        "net_er": "not exposed for P-gp transfectant",
        "inhibitor": "verapamil 25 uM | itraconazole 10 uM",
        "inhibitor_er": "54% | 38% reduction in Caco-2",
        "positive_control": "not exposed",
        "source_call": "likely MDR1 substrate; Caco-2 result called P-gp substrate",
        "evidence_tier": "A2_transfectant_directional_plus_B2_Caco2_rescue",
        "strict_matched_parent_candidate": "no",
        "qc_note": "direction is credible but strict P-gp transfectant numeric endpoint is absent; retain outside formal metric",
    },
]


def fetch_pubchem_with_fallback(name: str) -> tuple[dict, str, str]:
    record, status, error = fetch_pubchem(name)
    return record, status, error


def strict_row(row: dict) -> dict:
    return {
        "strict_candidate_id": row["batch3_entity_id"],
        "entity_name": row["entity_name"],
        "directional_label": row["label"],
        "collection_origin": "regulatory_expansion_batch3",
        "evidence_tier": row["evidence_tier"],
        "assay_system": row["assay_system"],
        "concentration_uM": row["concentration_uM"],
        "er_transporter": row["er_transporter"],
        "er_parental": row["er_parental"],
        "net_er": row["net_er"],
        "inhibitor": row["inhibitor"],
        "inhibitor_er": row["inhibitor_er"],
        "positive_control": row["positive_control"],
        "standardized_inchikey": row["standardized_inchikey"],
        "standardized_isomeric_smiles": row["standardized_isomeric_smiles"],
        "source_id": row["source_id"],
        "source_url": row["source_url"],
        "source_location": row["source_location"],
        "identity_status": row["identity_status"],
        "source_study_ids": row["source_study_ids"],
        "upstream_experiment_cluster": row["upstream_experiment_cluster"],
        "claim_boundary": row["claim_boundary"],
    }


def ledger_row(row: dict) -> dict:
    is_strict = row["strict_external_candidate_eligible"] == "yes"
    inhibitor_context = row["inhibitor"]
    if row["inhibitor_er"]:
        inhibitor_context += f"; {row['inhibitor_er']}"
    return {
        "collection_candidate_id": f"EXTCOLL::{row['batch3_entity_id']}",
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
        "quantitative_endpoint": row["net_er"] or row["er_transporter"],
        "endpoint_definition": f"{row['assay_system']}; source call: {row['source_call']}",
        "matched_parental_control": "yes" if is_strict else "numeric_value_not_exposed",
        "selective_inhibitor_rescue": inhibitor_context,
        "source_id": row["source_id"],
        "source_cluster_count": 1,
        "source_clusters": row["upstream_experiment_cluster"],
        "standardized_inchikey": row["standardized_inchikey"],
        "inchikey_connectivity_block": row["inchikey_connectivity_block"],
        "standardized_isomeric_smiles": row["standardized_isomeric_smiles"],
        "source_url": row["source_url"],
        "source_location": row["source_location"],
        "evidence_tier": row["evidence_tier"],
        "remaining_gate": (
            "complete source-lineage/vendor audit and score once under frozen protocol"
            if is_strict
            else "retain outside formal external metric; use only in tier-stratified sensitivity analysis"
        ),
        "claim_boundary": row["claim_boundary"],
    }


def main() -> None:
    prior_indexes = {name: identity_index(path, columns) for name, (path, columns) in PRIOR_DATASETS.items()}
    rows: list[dict] = []
    for index, manual in enumerate(ENTITIES, start=1):
        source = SOURCES[manual["source_key"]]
        pubchem, pubchem_status, pubchem_error = fetch_pubchem_with_fallback(manual["entity_name"])
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
        eligible = strict and identity_status == "no_exact_or_connectivity_prior_overlap"
        local_pdf = RAW_DIR / source["local_pdf"]
        local_text = RAW_DIR / source["local_text"]
        rows.append(
            {
                "batch3_entity_id": f"ABCB1_REG_BATCH3_{index:03d}",
                **manual,
                **pubchem,
                "pubchem_resolution_status": pubchem_status,
                "pubchem_resolution_error": pubchem_error,
                **standard,
                "inchikey_connectivity_block": block,
                "overlap_exact_prior_sets": " | ".join(exact_hits),
                "overlap_connectivity_prior_sets": " | ".join(connectivity_hits),
                "identity_status": identity_status,
                "strict_external_candidate_eligible": "yes" if eligible else "no",
                "source_id": manual["source_key"],
                "source_url": source["url"],
                "source_location": source["location"],
                "source_study_ids": source["study_ids"],
                "upstream_experiment_cluster": source["lineage_unit"],
                "source_lineage_status": "regulatory_report_and_sponsor_study_ids_locked_vendor_independence_pending",
                "local_pdf": str(local_pdf.relative_to(PROJECT)),
                "local_pdf_sha256": sha256(local_pdf) if local_pdf.exists() else "missing",
                "local_text": str(local_text.relative_to(PROJECT)),
                "claim_boundary": "pre-scoring regulatory candidate; not a validation result and not evidence for the ADC release-species domain",
            }
        )

    write_csv(PROCESSED / "abcb1_regulatory_expansion_batch3_entities_v1.csv", rows)

    source_rows = []
    lineage_rows = []
    for source_id, source in SOURCES.items():
        local_pdf = RAW_DIR / source["local_pdf"]
        source_entities = [row for row in rows if row["source_id"] == source_id]
        source_rows.append(
            {
                "source_id": source_id,
                "source_url": source["url"],
                "source_location": source["location"],
                "local_pdf": str(local_pdf.relative_to(PROJECT)),
                "local_pdf_sha256": sha256(local_pdf) if local_pdf.exists() else "missing",
                "entity_count": len(source_entities),
                "strict_identity_clean_count": sum(row["strict_external_candidate_eligible"] == "yes" for row in source_entities),
                "allowed_use": "pre-scoring evidence-tiered ABCB1 external candidate collection",
            }
        )
        lineage_rows.append(
            {
                "source_id": source_id,
                "entity_names": " | ".join(row["entity_name"] for row in source_entities),
                "regulatory_report_url": source["url"],
                "source_study_ids": source["study_ids"],
                "upstream_experiment_cluster": source["lineage_unit"],
                "within_batch_cluster_size": len(source_entities),
                "cross_report_shared_study_id": "none_observed_in_batch3",
                "assay_vendor": "not_exposed_in_regulatory_summary",
                "independence_status": "distinct_sponsor_experiment_cluster_pending_vendor_and_publication_backtrace",
                "analysis_rule": "cluster bootstrap/resampling by upstream_experiment_cluster; do not treat report count as proof of laboratory independence",
            }
        )
    write_csv(PROCESSED / "abcb1_regulatory_expansion_batch3_source_registry_v1.csv", source_rows)
    write_csv(PROCESSED / "abcb1_regulatory_expansion_batch3_source_lineage_v1.csv", lineage_rows)

    tiered_v2 = read_csv(PROCESSED / "abcb1_regulatory_evidence_tiered_pool_v2.csv")
    tiered_v3 = list(tiered_v2)
    for row in rows:
        tiered_v3.append({"entity_record_id": row["batch3_entity_id"], "collection_batch": "regulatory_expansion_batch3", **row})
    write_csv(PROCESSED / "abcb1_regulatory_evidence_tiered_pool_v3.csv", tiered_v3)

    strict_v2 = read_csv(PROCESSED / "abcb1_regulatory_strict_external_pool_v2.csv")
    batch3_strict = [row for row in rows if row["strict_external_candidate_eligible"] == "yes"]
    strict_v3 = list(strict_v2) + [strict_row(row) for row in batch3_strict]
    write_csv(PROCESSED / "abcb1_regulatory_strict_external_pool_v3.csv", strict_v3)

    ledger_v2 = read_csv(PROCESSED / "abcb1_external_collection_candidate_ledger_v2.csv")
    identity_clean_rows = [row for row in rows if row["identity_status"] == "no_exact_or_connectivity_prior_overlap"]
    ledger_v3 = list(ledger_v2) + [ledger_row(row) for row in identity_clean_rows]
    write_csv(PROCESSED / "abcb1_external_collection_candidate_ledger_v3.csv", ledger_v3)

    queue_v2 = read_csv(PROCESSED / "abcb1_regulatory_expansion_search_queue_v2.csv")
    extracted_names = {row["entity_name"] for row in rows}
    queue_v3 = []
    for row in queue_v2:
        new_row = dict(row)
        if row["candidate_names"] in extracted_names:
            new_row["status"] = "extracted_and_tiered_batch3"
        queue_v3.append(new_row)
    write_csv(PROCESSED / "abcb1_regulatory_expansion_search_queue_v3.csv", queue_v3)

    strict_counts = Counter(str(row["directional_label"]) for row in strict_v3)
    role_counts = Counter(row["formal_role_now"] for row in ledger_v3)
    summary = {
        "batch3_source_documents": len(source_rows),
        "batch3_entity_rows": len(rows),
        "batch3_standardization_ok": sum(row["standardization_status"] == "ok" for row in rows),
        "batch3_identity_clean_rows": len(identity_clean_rows),
        "batch3_strict_identity_clean_candidates": len(batch3_strict),
        "batch3_strict_names": [row["entity_name"] for row in batch3_strict],
        "batch3_secondary_names": [row["entity_name"] for row in rows if row not in batch3_strict],
        "batch3_identity_status_counts": dict(Counter(row["identity_status"] for row in rows)),
        "regulatory_strict_pool_v3_rows": len(strict_v3),
        "regulatory_strict_pool_v3_class_counts": dict(strict_counts),
        "external_candidate_ledger_v3_rows": len(ledger_v3),
        "external_candidate_ledger_v3_role_counts": dict(role_counts),
        "formal_external_shortfall_if_minimum_15_per_class": {
            "positive": max(0, 15 - strict_counts.get("1", 0)),
            "negative": max(0, 15 - strict_counts.get("0", 0)),
        },
        "source_lineage_cluster_count_batch3": len({row["upstream_experiment_cluster"] for row in rows}),
        "model_scoring_performed": False,
        "formal_status": "positive_side_candidate_count_gate_reached_if_identity_clean; negative_side_and_vendor_lineage_audit_pending",
    }
    (PROCESSED / "abcb1_regulatory_expansion_batch3_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
