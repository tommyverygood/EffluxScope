#!/usr/bin/env python3
"""Build the eighth gap-targeted ADC payload evidence batch.

This batch adds SG3199 as an exact released tesirine warhead, preserves the
FDA-versus-Roth ABCG2 conflict without adjudicating a project label, and adds a
separate exploratory PAMPA membrane-permeability layer from Khera et al.
Figure 4B.  It deliberately does not add binary model labels or promote PAMPA
values to transporter truth.
"""

from __future__ import annotations

import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


DATA = Path(__file__).resolve().parents[1]
PROCESSED = DATA / "processed"
RAW = DATA / "raw" / "eighth_batch_gap_targeted"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    if not rows and not fieldnames:
        raise ValueError(f"Cannot infer columns for empty table: {path}")
    columns = fieldnames or list(rows[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_unique(rows: list[dict[str, str]], key: str) -> None:
    values = [row[key] for row in rows]
    duplicates = [value for value, count in Counter(values).items() if count > 1]
    if duplicates:
        raise ValueError(f"Duplicate {key}: {duplicates}")


def build_sources() -> list[dict[str, object]]:
    return [
        {
            "source_id": "SRC_B8_001",
            "title": "Cellular-Resolution Imaging of Bystander Payload Tissue Penetration from Antibody-Drug Conjugates",
            "year": 2022,
            "source_type": "primary_article_figure_digitization",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8832507/",
            "use_in_batch": "PAMPA passive permeability means for nine ADC payloads/release species",
            "source_location": "Methods: PAMPA Permeability Assay; Figure 4B; article page 315",
            "primary_or_secondary": "primary_article",
            "notes": "Figure-only means were digitized at 220 dpi; error bars were not digitized. Values are exploratory and truth-set ineligible.",
        },
        {
            "source_id": "SRC_B8_002",
            "title": "Zynlonta BLA 761196 Multidiscipline Review",
            "year": 2021,
            "source_type": "primary_regulatory_summary",
            "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/761196Orig1s000MultidisciplineR.pdf",
            "use_in_batch": "source-assigned SG3199 P-gp substrate and BCRP non-substrate calls",
            "source_location": "PDF physical page 46; printed page 45",
            "primary_or_secondary": "primary_regulatory_summary",
            "notes": "The review does not reproduce assay system, concentration, Papp or ER. Preserve calls without creating a model label.",
        },
        {
            "source_id": "SRC_B8_003",
            "title": "Pre-clinical pharmacology and mechanism of action of SG3199, the PBD dimer warhead component of ADC payload tesirine",
            "year": 2018,
            "source_type": "primary_article_release_identity",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6041317/",
            "use_in_batch": "establish SG3199 as the released warhead component of tesirine",
            "source_location": "Abstract; Introduction; Figure 1",
            "primary_or_secondary": "primary_article",
            "notes": "DOI 10.1038/s41598-018-28533-4. The source explicitly uses the phrase released warhead component.",
        },
        {
            "source_id": "SRC_B8_004",
            "title": "PubChem compound record for SG3199",
            "year": 2026,
            "source_type": "primary_database_registry",
            "url": "https://pubchem.ncbi.nlm.nih.gov/compound/90132565",
            "use_in_batch": "machine-readable SG3199 structure crosswalk",
            "source_location": "CID 90132565 property record retrieved 2026-07-28",
            "primary_or_secondary": "primary_database",
            "notes": "CID, SMILES, InChIKey, formula and descriptors agree with the existing payload registry row ADC_PAYLOAD_032.",
        },
        {
            "source_id": "SRC_B8_005",
            "title": "Elahere NDA 761310 Multidiscipline Review",
            "year": 2022,
            "source_type": "primary_regulatory_recheck",
            "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2022/761310Orig1s000MultidisciplineR.pdf",
            "use_in_batch": "targeted recheck for missing numeric S-methyl-DM4 transport data",
            "source_location": "transporter summary",
            "primary_or_secondary": "primary_regulatory_summary",
            "notes": "Only qualitative DM4/S-methyl-DM4 P-gp calls were found; no new numeric Papp/ER was extracted.",
        },
    ]


def build_raw_log() -> list[dict[str, object]]:
    specs = [
        ("Khera_2022_ADC_payload_bystander_PAMPA_article.pdf", "SRC_B8_001", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8832507/", "downloaded_primary_article_pdf_page6_rendered_and_visually_checked", "pdf"),
        ("PMC8832507_BioC_fulltext.json", "SRC_B8_001", "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC8832507/unicode", "downloaded_primary_fulltext_bioc_json", "json"),
        ("FDA_Zynlonta_multidiscipline_review_2021.pdf", "SRC_B8_002", "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/761196Orig1s000MultidisciplineR.pdf", "downloaded_official_pdf_pages45_46_rendered_and_visually_checked", "pdf"),
        ("PMC6041317_SG3199_pharmacology_fulltext.xml", "SRC_B8_003", "https://pmc.ncbi.nlm.nih.gov/articles/PMC6041317/", "downloaded_primary_fulltext_bioc_xml", "xml"),
        ("PubChem_SG3199_properties_20260728.json", "SRC_B8_004", "https://pubchem.ncbi.nlm.nih.gov/compound/90132565", "downloaded_pubchem_property_json", "json"),
        ("FDA_Elahere_multidiscipline_review_2022.pdf", "SRC_B8_005", "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2022/761310Orig1s000MultidisciplineR.pdf", "downloaded_official_pdf_targeted_recheck_no_new_numeric", "pdf"),
        ("EXCLUDED_PMDA_GS5734_remdesivir_false_search_hit.pdf", "EXCLUDED_B8_001", "https://www.pmda.go.jp/drugs/2020/P20200518001/340018000_30200AMX00434000_H100_1.pdf", "excluded_false_search_hit_content_is_GS5734_remdesivir_not_Enhertu", "pdf"),
    ]
    rows = []
    for filename, source_id, url, status, kind in specs:
        path = RAW / filename
        if not path.exists():
            raise FileNotFoundError(path)
        signature = "no"
        if kind == "pdf":
            signature = "yes" if path.read_bytes()[:5] == b"%PDF-" else "no"
        elif kind == "json":
            json.loads(path.read_text(encoding="utf-8"))
            signature = "yes"
        elif kind == "xml":
            ET.parse(path)
            signature = "yes"
        rows.append(
            {
                "filename": f"raw/eighth_batch_gap_targeted/{filename}",
                "source_id": source_id,
                "source_url": url,
                "download_status": status,
                "signature_valid": signature,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
                "acquired_or_verified_at": "2026-07-28",
            }
        )
    return rows


def build_release_registry() -> list[dict[str, str]]:
    rows = read_csv(PROCESSED / "adc_release_species_registry_v3.csv")
    header = list(rows[0].keys())
    if not any(row["release_species_id"] == "ADC_REL_013" for row in rows):
        compound = next(row for row in read_csv(PROCESSED / "compound_registry_v1.csv") if row["compound_id"] == "ADC_PAYLOAD_032")
        new = {column: "" for column in header}
        new.update(
            {
                "release_species_id": "ADC_REL_013",
                "release_species_name": "SG3199",
                "aliases": "SG3199-PBD | released tesirine warhead",
                "payload_core": "SG3199",
                "residual_group": "none",
                "species_class": "free_payload",
                "representative_adc_or_format": "loncastuximab tesirine | camidanlumab tesirine | rovalpituzumab tesirine",
                "linker_class": "enzyme-cleavable tesirine payload",
                "release_identity_status": "confirmed_released_warhead",
                "release_evidence_method": "primary article explicit released-warhead assignment",
                "source_ids": "SRC_B8_002 | SRC_B8_003",
                "structure_query_names": "SG3199 | SG3199-PBD",
                "existing_payload_compound_id": "ADC_PAYLOAD_032",
                "modeling_status": "already_in_free_payload_registry",
                "pubchem_query_name": "SG3199",
                "pubchem_cid": compound["pubchem_cid"],
                "pubchem_title": compound["pubchem_title"],
                "isomeric_smiles": compound["isomeric_smiles"],
                "connectivity_smiles": compound["connectivity_smiles"],
                "inchikey": compound["inchikey"],
                "molecular_formula": compound["molecular_formula"],
                "molecular_weight": compound["molecular_weight"],
                "xlogp": compound["xlogp"],
                "tpsa": compound["tpsa"],
                "formal_charge": compound["formal_charge"],
                "structure_status": "resolved",
                "structure_note": "Exact SG3199 structure cross-linked to existing payload registry ADC_PAYLOAD_032 and PubChem CID 90132565.",
                "structure_source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/90132565",
            }
        )
        rows.append(new)
    assert_unique(rows, "release_species_id")
    write_csv(PROCESSED / "adc_release_species_registry_v4.csv", rows, header)
    return rows


def build_release_transport_evidence() -> list[dict[str, str]]:
    rows = read_csv(PROCESSED / "adc_release_species_transporter_evidence_long_v5.csv")
    header = list(rows[0].keys())

    def blank() -> dict[str, str]:
        return {column: "" for column in header}

    additions: list[dict[str, str]] = []
    for record_id, gene, protein, source_call, source_label, note in [
        ("ADC_TRANS_B8_FDA_SG3199_ABCB1", "ABCB1", "P-gp", "substrate", "1", "FDA source call is directionally concordant with Roth RR 120, but the FDA review reproduces no assay details or numeric result."),
        ("ADC_TRANS_B8_FDA_SG3199_ABCG2", "ABCG2", "BCRP", "not_substrate", "0", "FDA source call conflicts with Roth HEK293-R5 cytotoxicity RR 20; preserve both contexts and leave the project label blank."),
    ]:
        row = blank()
        row.update(
            {
                "transport_record_id": record_id,
                "release_species_id": "ADC_REL_013",
                "transporter_gene": gene,
                "transporter_protein": protein,
                "species": "Homo sapiens",
                "assay_system": "human in vitro transporter studies summarized by FDA",
                "test_concentrations": "not reproduced in review",
                "source_assigned_call": source_call,
                "source_binary_label": source_label,
                "model_binary_label": "",
                "model_label_status": "not_assigned_regulatory_call_preserved_conflict_audited",
                "evidence_tier": "B_regulatory_call_no_numeric",
                "truth_set_candidate": "call_only_cross_endpoint_conflict_audit",
                "source_id": "SRC_B8_002",
                "source_location": "PDF physical page 46; printed page 45",
                "curation_note": note,
                "release_species_name": "SG3199",
                "source_url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/761196Orig1s000MultidisciplineR.pdf",
                "assay_endpoint": "regulatory substrate call",
            }
        )
        additions.append(row)

    for record_id, gene, protein, cell, control, control_sd, treated, treated_sd, ratio in [
        ("ADC_TRANS_B8_ROTH_SG3199_ABCB1", "ABCB1", "P-gp", "HEK293-MDR19", "0.045", "0.011", "5.6", "1.6", "120"),
        ("ADC_TRANS_B8_ROTH_SG3199_ABCG2", "ABCG2", "BCRP", "HEK293-R5", "0.045", "0.011", "0.90", "0.31", "20"),
        ("ADC_TRANS_B8_ROTH_SG3199_ABCC1", "ABCC1", "MRP1", "HEK293-MRP1", "0.045", "0.011", "0.21", "0.026", "4.6"),
    ]:
        row = blank()
        row.update(
            {
                "transport_record_id": record_id,
                "release_species_id": "ADC_REL_013",
                "transporter_gene": gene,
                "transporter_protein": protein,
                "species": "Homo sapiens",
                "assay_system": "HEK293 stable transfectant cytotoxicity",
                "test_concentrations": "concentration-response; IC50 in nM",
                "source_assigned_call": "",
                "source_binary_label": "",
                "model_binary_label": "",
                "model_label_status": "not_assigned_continuous_rr_preserved",
                "evidence_tier": "B_functional_transfectant_resistance_ratio",
                "truth_set_candidate": "yes_functional_only_no_binary_cutoff",
                "source_id": "SRC_B6_002",
                "source_location": "Table 3",
                "curation_note": f"Exact SG3199 is now linked to an actual released tesirine warhead. IC50 SD values are {control_sd} nM control and {treated_sd} nM transporter condition; functional RR is not an intrinsic transport rate.",
                "release_species_name": "SG3199",
                "source_doi": "10.20517/cdr.2025.151",
                "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12883344/",
                "assay_endpoint": "IC50 cytotoxicity",
                "control_condition": "HEK293-pcDNA",
                "control_value": f"{control} nM",
                "transporter_condition": cell,
                "transporter_value": f"{treated} nM",
                "effect_ratio": ratio,
                "effect_ratio_definition": f"IC50({cell}) / IC50(HEK293-pcDNA)",
            }
        )
        additions.append(row)

    existing_ids = {row["transport_record_id"] for row in rows}
    rows.extend(row for row in additions if row["transport_record_id"] not in existing_ids)
    assert_unique(rows, "transport_record_id")
    if any(row["model_binary_label"].strip() for row in rows):
        raise ValueError("model_binary_label must remain blank")
    write_csv(PROCESSED / "adc_release_species_transporter_evidence_long_v6.csv", rows, header)
    return rows


def build_conflict_table() -> list[dict[str, object]]:
    rows = [
        {
            "conflict_record_id": "SG3199_CTX_001",
            "release_species_id": "ADC_REL_013",
            "entity_name": "SG3199",
            "transporter_gene": "ABCB1",
            "regulatory_source_call": "substrate",
            "regulatory_numeric_result": "not reproduced",
            "functional_source_endpoint": "HEK293-MDR19 / HEK293-pcDNA cytotoxicity resistance ratio",
            "functional_source_value": 120,
            "directional_relationship": "concordant_positive",
            "adjudication": "preserve_both_no_project_binary_label",
            "source_ids": "SRC_B8_002 | SRC_B6_002",
            "interpretation": "Same exact released warhead, but endpoint and assay details are not aligned; this is corroboration, not M0/M1 direct transport replication.",
        },
        {
            "conflict_record_id": "SG3199_CTX_002",
            "release_species_id": "ADC_REL_013",
            "entity_name": "SG3199",
            "transporter_gene": "ABCG2",
            "regulatory_source_call": "not_substrate",
            "regulatory_numeric_result": "not reproduced",
            "functional_source_endpoint": "HEK293-R5 / HEK293-pcDNA cytotoxicity resistance ratio",
            "functional_source_value": 20,
            "directional_relationship": "discordant_context_dependent",
            "adjudication": "explicit_conflict_preserve_both_no_project_binary_label",
            "source_ids": "SRC_B8_002 | SRC_B6_002",
            "interpretation": "The discordance may reflect assay sensitivity, concentration, endpoint or unreported regulatory method; it cannot be resolved from the available documents.",
        },
    ]
    write_csv(PROCESSED / "sg3199_cross_endpoint_conflict_v1.csv", rows)
    return rows


def build_matched_comparisons() -> list[dict[str, str]]:
    rows = read_csv(PROCESSED / "matched_context_comparisons_v3.csv")
    header = list(rows[0].keys())

    def make(comparison_id: str, gene: str, fda_call: str, rr: str, direction: str, limitation: str) -> dict[str, str]:
        row = {column: "" for column in header}
        row.update(
            {
                "comparison_id": comparison_id,
                "payload_core": "SG3199",
                "transporter_gene": gene,
                "free_entity_id": "ADC_REL_013",
                "free_entity_name": "SG3199 in FDA Zynlonta review",
                "released_entity_id": "ADC_REL_013",
                "released_entity_name": "SG3199 in Roth 2026",
                "free_evidence": f"FDA regulatory source call: {fda_call}; assay details and numeric result not reproduced",
                "released_evidence": f"Roth HEK293 stable-transfectant cytotoxicity resistance ratio {rr}",
                "matched_grade": "M4",
                "endpoint_alignment": "regulatory_qualitative_call_vs_functional_cytotoxicity_RR",
                "cell_background_alignment": "FDA_background_not_reproduced_vs_HEK293_transfectant",
                "concentration_alignment": "FDA_concentration_not_reproduced_vs_concentration_response",
                "identity_alignment": "same_exact_SG3199_released_warhead",
                "directional_result": direction,
                "cross_assay_fold_computed": "no",
                "training_use": "cross_endpoint_context_audit_only",
                "source_ids": "SRC_B8_002 | SRC_B6_002",
                "limitation": limitation,
            }
        )
        return row

    additions = [
        make("MATCH_B8_001", "ABCB1", "P-gp substrate", "120", "concordant positive direction", "FDA gives no assay system, concentration, Papp or ER; do not upgrade to direct independent transport replication."),
        make("MATCH_B8_002", "ABCG2", "not a BCRP substrate", "20", "discordant source contexts", "The available documents cannot adjudicate whether the conflict reflects assay endpoint, sensitivity, concentration or method; no project binary label is assigned."),
    ]
    existing_ids = {row["comparison_id"] for row in rows}
    rows.extend(row for row in additions if row["comparison_id"] not in existing_ids)
    assert_unique(rows, "comparison_id")
    write_csv(PROCESSED / "matched_context_comparisons_v4.csv", rows, header)
    return rows


def build_passive_permeability_view() -> list[dict[str, str]]:
    digitized = read_csv(PROCESSED / "khera_figure4b_pampa_digitization_v1.csv")
    rows = []
    for row in digitized:
        if not row["release_species_id"]:
            continue
        rows.append(
            {
                "release_species_id": row["release_species_id"],
                "release_species_name": row["normalized_entity_name"],
                "assay": row["assay"],
                "working_concentration_value": row["working_concentration_value"],
                "working_concentration_unit": row["working_concentration_unit"],
                "pampa_mean_digitized_1e6_cm_s": row["pampa_consensus_1e6_cm_s"],
                "one_pixel_uncertainty_1e6_cm_s": row["one_pixel_value_uncertainty_1e6_cm_s"],
                "error_bar_numeric_value": "",
                "numeric_origin": row["numeric_origin"],
                "truth_set_eligible": "no",
                "model_binary_label": "",
                "source_id": "SRC_B8_001",
                "source_location": row["source_location"],
                "interpretation_layer": "passive_membrane_permeability_not_transporter_substrate_truth",
            }
        )
    assert_unique(rows, "release_species_id")
    write_csv(PROCESSED / "adc_release_species_passive_permeability_v1.csv", rows)
    return rows


def update_issue_log() -> list[dict[str, str]]:
    path = DATA / "exclusion_and_conflict_log.csv"
    rows = read_csv(path)
    header = list(rows[0].keys())
    additions = [
        {
            "issue_id": "ISSUE-027",
            "entity_or_source": "SG3199 ABCG2",
            "issue_type": "cross_endpoint_source_call_conflict",
            "current_action": "preserve_FDA_non_substrate_call_and_Roth_RR20_separately_no_project_label",
            "reason": "FDA states SG3199 is not a BCRP substrate without reproducing numeric methods, while Roth reports HEK293-R5 cytotoxicity RR 20 for exact SG3199.",
            "next_evidence_needed": "Independent single-transporter bidirectional Papp or accumulation assay with concentration, recovery, controls and inhibitor rescue",
            "status": "open",
        },
        {
            "issue_id": "ISSUE-028",
            "entity_or_source": "Khera 2022 Figure 4B PAMPA",
            "issue_type": "single_figure_digitization_no_independent_curator",
            "current_action": "retain_nine_mean_bars_as_exploratory_membrane_layer_exclude_from_truth_set",
            "reason": "The article shows PAMPA means only as bars; error bars and analyte-specific 2.5 h versus 5 h assignments are not tabulated in the accessible main text.",
            "next_evidence_needed": "Author numerical table or independent second-curator digitization plus supplementary-method recovery",
            "status": "partially_resolved",
        },
        {
            "issue_id": "ISSUE-029",
            "entity_or_source": "PMDA search hit labeled as Enhertu CTD",
            "issue_type": "search_index_false_positive_wrong_compound",
            "current_action": "exclude_and_rename_raw_file_as_GS5734_remdesivir",
            "reason": "Visual/text verification showed the downloaded official PDF concerns GS-5734/remdesivir, not DXd or trastuzumab deruxtecan.",
            "next_evidence_needed": "Correct product-specific CTD source with explicit DXd single-transporter numeric tables",
            "status": "resolved_exclusion",
        },
    ]
    existing = {row["issue_id"] for row in rows}
    rows.extend(row for row in additions if row["issue_id"] not in existing)
    assert_unique(rows, "issue_id")
    write_csv(path, rows, header)
    return rows


def main() -> None:
    # Digitization must be run first so the passive-permeability view can be built.
    pampa_path = PROCESSED / "khera_figure4b_pampa_digitization_v1.csv"
    if not pampa_path.exists():
        raise FileNotFoundError(f"Run digitize_khera_figure4b.py first: {pampa_path}")

    sources = build_sources()
    raw_log = build_raw_log()
    registry = build_release_registry()
    transport = build_release_transport_evidence()
    conflicts = build_conflict_table()
    matched = build_matched_comparisons()
    passive = build_passive_permeability_view()
    issues = update_issue_log()

    write_csv(PROCESSED / "eighth_batch_source_registry.csv", sources)
    write_csv(PROCESSED / "eighth_batch_raw_download_log.csv", raw_log)

    pampa = read_csv(pampa_path)
    direct_numeric = read_csv(PROCESSED / "direct_bidirectional_transport_numeric_v3.csv")
    grade_counts = dict(sorted(Counter(row["matched_grade"] for row in matched).items()))
    open_or_partial = sum(row["status"] not in {"resolved", "resolved_exclusion"} for row in issues)

    validation = {
        "generated_at": "2026-07-28",
        "checks": {
            "pampa_rows_equal_9": len(pampa) == 9,
            "pampa_release_species_rows_equal_6": len(passive) == 6,
            "pampa_all_truth_set_ineligible": all(row["truth_set_eligible"] == "no" for row in pampa),
            "release_registry_entities_equal_13": len(registry) == 13,
            "release_transport_records_equal_39": len(transport) == 39,
            "release_transport_model_labels_all_blank": all(not row["model_binary_label"].strip() for row in transport),
            "sg3199_conflict_rows_equal_2": len(conflicts) == 2,
            "sg3199_abcg2_conflict_preserved": any(row["transporter_gene"] == "ABCG2" and row["directional_relationship"] == "discordant_context_dependent" for row in conflicts),
            "matched_m0_still_zero": grade_counts.get("M0", 0) == 0,
            "raw_signatures_all_valid": all(row["signature_valid"] == "yes" for row in raw_log),
            "direct_bidirectional_numeric_rows_unchanged_61": len(direct_numeric) == 61,
        },
    }
    validation["all_checks_passed"] = all(validation["checks"].values())
    if not validation["all_checks_passed"]:
        raise ValueError(f"Validation failed: {validation}")
    (PROCESSED / "eighth_batch_validation.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    summary = {
        "generated_at": "2026-07-28",
        "eighth_batch_relevant_sources": len(sources),
        "excluded_false_search_hits": 1,
        "raw_source_files": len(raw_log),
        "new_pampa_digitized_rows": len(pampa),
        "registered_release_species_with_pampa_rows": len(passive),
        "release_species_registry_v4_entities": len(registry),
        "release_transporter_records_v6": len(transport),
        "new_sg3199_release_transport_records": 5,
        "direct_bidirectional_transport_numeric_rows_v3": len(direct_numeric),
        "released_species_with_author_tabulated_bidirectional_papp_count": 5,
        "matched_comparisons_v4": len(matched),
        "matched_grade_counts": grade_counts,
        "m0_exact_same_assay_pairs": grade_counts.get("M0", 0),
        "m1_same_entity_independent_direct_pairs": grade_counts.get("M1", 0),
        "explicit_cross_endpoint_conflict_rows": len(conflicts),
        "model_binary_labels_created": 0,
        "open_or_partial_issue_count": open_or_partial,
        "key_results": [
            "SG3199 is added as ADC_REL_013 using a primary explicit released-warhead assignment and the existing exact PubChem-resolved structure.",
            "Five SG3199 release-evidence rows link two FDA source calls and three Roth continuous functional resistance ratios.",
            "FDA calls SG3199 a P-gp substrate but not a BCRP substrate; Roth reports RR 120 for ABCB1 and RR 20 for ABCG2. The ABCG2 discordance is preserved as a context conflict.",
            "Khera Figure 4B adds nine exploratory PAMPA mean bars; six map to registered actual release species. These values are not transporter labels.",
            "No new M0 matched free-payload/released-species pair and no new gene-isolated bidirectional Papp table were found.",
            "A PMDA search-index hit was verified as remdesivir/GS-5734 rather than Enhertu and excluded.",
        ],
    }
    (PROCESSED / "eighth_batch_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
