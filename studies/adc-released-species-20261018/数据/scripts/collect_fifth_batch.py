#!/usr/bin/env python3
"""Build fifth-batch matched transport evidence tables for the ADC predictor.

The batch prioritizes direct numeric transporter experiments and explicitly
grades free-payload/released-species comparability.  No cross-assay fold change
or binary model label is created.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "processed"
RAW5 = ROOT / "raw" / "fifth_batch_matched_pairs"
PROCESSED.mkdir(parents=True, exist_ok=True)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    fields = fieldnames or (list(rows[0].keys()) if rows else [])
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


SOURCES = [
    {
        "source_id": "SRC_B5_001",
        "title": "Kadcyla CTD 2.6.4: pharmacokinetics written summary",
        "year": "2013",
        "source_type": "regulatory_primary_submission_document",
        "url": "https://www.pmda.go.jp/drugs/2013/P201300137/450045000_22500AMX01816_I100_2.pdf",
        "use_in_batch": "direct numeric bidirectional transport of free DM1 in MDCKII and MDCKII-MDR1 with PSC833/verapamil rescue",
        "source_location": "PDF pages 33-34; section 2.6.4.7.2",
        "primary_or_secondary": "primary_regulatory",
        "notes": "DM1 0.5 uM: parental ER 0.89, MDR1 ER 6.65, background-corrected ER 7.47; separate inhibitor run ER 5.90, 1.10 with PSC833, and 1.76 with verapamil.",
    },
    {
        "source_id": "SRC_B5_002",
        "title": "Target-Site Drug Metabolism and Transport",
        "year": "2015",
        "source_type": "primary_article_figure",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11024933/",
        "use_in_batch": "Lys-MCC-DM1 direct bidirectional permeability Figure 8 and source conclusion",
        "source_location": "Figure 8",
        "primary_or_secondary": "primary_article",
        "notes": "DOI 10.1124/dmd.115.064576. Figure means were algorithmically digitized, but no independent second curator is available; numeric rows remain exploratory.",
    },
    {
        "source_id": "SRC_B5_003",
        "title": "Antibody-maytansinoid conjugates designed to bypass multidrug resistance",
        "year": "2010",
        "source_type": "primary_article_reused",
        "url": "https://pubmed.ncbi.nlm.nih.gov/20197459/",
        "use_in_batch": "same-system SMCC-DM1 versus PEG4Mal-DM1 functional and intracellular-metabolite contrasts",
        "source_location": "Table 1; Figures 3-6; Results and Discussion",
        "primary_or_secondary": "primary_article",
        "notes": "DOI 10.1158/0008-5472.CAN-09-3546. These are mechanistic same-system contrasts, not isolated bidirectional transport of both catabolites.",
    },
    {
        "source_id": "SRC_B5_004",
        "title": "FDA Enhertu multidisciplinary review",
        "year": "2019",
        "source_type": "regulatory_review",
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2019/761139Orig1s000MultidisciplineR.pdf",
        "use_in_batch": "regulatory substrate calls for released DXd/MAAA-1181a",
        "source_location": "PDF page 85 (printed page 84); clinical pharmacology sections",
        "primary_or_secondary": "primary_regulatory_summary",
        "notes": "MAAA-1181a is identified as a P-gp and BCRP substrate, among other transporters; no numeric P-gp/BCRP transport output is reproduced.",
    },
    {
        "source_id": "SRC_B5_005",
        "title": "Enhertu CTD 2.6.4: pharmacokinetics written summary",
        "year": "2020",
        "source_type": "regulatory_primary_submission_document",
        "url": "https://www.pmda.go.jp/drugs/2020/P20200420002/430574000_30200AMX00425_I100_1.pdf",
        "use_in_batch": "confirm MAAA-1181a identity as released DXd and primary excreted component",
        "source_location": "2.6.4 summary and tabulated sections",
        "primary_or_secondary": "primary_regulatory",
        "notes": "Supports released-species identity; this CTD file does not reproduce a numeric P-gp/BCRP substrate table.",
    },
    {
        "source_id": "SRC_B5_006",
        "title": "US7750116B1: Methods for treating drug-resistant cancer",
        "year": "2010",
        "source_type": "patent_primary_experiment",
        "url": "https://patents.google.com/patent/US7750116B1/en",
        "use_in_batch": "MMAF and vcMMAF-ADC cytotoxicity with verapamil in NCI/ADR-RES cells",
        "source_location": "PDF page 37; printed pages 61-62; Table 7",
        "primary_or_secondary": "primary_patent",
        "notes": "Table 7 reports +10 uM verapamil, while the immediately preceding methods text says 40 uM. This internal concentration conflict is retained.",
    },
    {
        "source_id": "SRC_B5_007",
        "title": "FDA Kadcyla pharmacology review",
        "year": "2013",
        "source_type": "regulatory_review",
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2013/125427Orig1s000PharmR.pdf",
        "use_in_batch": "qualitative regulatory corroboration of DM1 P-gp study",
        "source_location": "PDF page 83 (printed page 73)",
        "primary_or_secondary": "primary_regulatory_summary",
        "notes": "Review conclusion: DM1 is probably a P-gp substrate and is not an inhibitor; numeric values are taken from the PMDA CTD, not this review.",
    },
    {
        "source_id": "SRC_B5_008",
        "title": "FDA Elahere multidisciplinary review",
        "year": "2022",
        "source_type": "regulatory_review_reused",
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2022/761310Orig1s000MultidisciplineR.pdf",
        "use_in_batch": "same-source qualitative calls for DM4 and S-methyl-DM4",
        "source_location": "PDF page 77 (printed page 76)",
        "primary_or_secondary": "primary_regulatory_summary",
        "notes": "Both DM4 and S-methyl-DM4 are stated to be MDR1/P-gp substrates; no numeric transport output is reproduced.",
    },
    {
        "source_id": "SRC_B5_009",
        "title": "PMDA belantamab mafodotin report",
        "year": "2025",
        "source_type": "regulatory_review_reused",
        "url": "https://www.pmda.go.jp/files/000280862.pdf",
        "use_in_batch": "direct numeric cys-mcMMAF P-gp and BCRP transport",
        "source_location": "Section 4.5.3, page 32",
        "primary_or_secondary": "primary_regulatory_summary",
        "notes": "cys-mcMMAF 100 nM: P-gp ER 34 and 1.1 with verapamil 60 uM; BCRP ER <2 at 0.72, 7.2 and 14.4 nM.",
    },
]


def direct_numeric_rows() -> list[dict]:
    fields = [
        "transport_record_id", "entity_type", "entity_id", "entity_name", "payload_core",
        "transporter_gene", "cell_background", "transporter_system", "assay_endpoint",
        "concentration_value", "concentration_unit", "inhibitor", "inhibitor_concentration_value",
        "inhibitor_concentration_unit", "papp_ab_mean", "papp_ab_sd", "papp_ba_mean", "papp_ba_sd",
        "efflux_ratio_mean", "efflux_ratio_sd", "background_corrected_efflux_ratio",
        "source_assigned_call", "evidence_tier", "numeric_origin", "digitization_status",
        "truth_set_eligible", "comparability_group_id", "source_id", "source_url", "source_location", "notes",
    ]

    def row(**kwargs):
        return {field: kwargs.get(field, "") for field in fields}

    rows = [
        row(transport_record_id="B5_DM1_01", entity_type="free_payload", entity_id="ADC_PAYLOAD_008", entity_name="DM1", payload_core="DM1", transporter_gene="ABCB1", cell_background="MDCKII parental", transporter_system="parental control", assay_endpoint="bidirectional apparent permeability; ER=Papp B-A/Papp A-B", concentration_value="0.5", concentration_unit="uM", papp_ab_mean="3.02", papp_ab_sd="0.15", papp_ba_mean="2.69", papp_ba_sd="0.14", efflux_ratio_mean="0.89", efflux_ratio_sd="0.07", source_assigned_call="control_background", evidence_tier="A_direct_numeric", numeric_origin="source_table", truth_set_eligible="yes_continuous_no_binary_cutoff", comparability_group_id="DM1_MDCKII_0.5uM", source_id="SRC_B5_001", source_url=SOURCES[0]["url"], source_location="section 2.6.4.7.2", notes="Mean +/- SD, n=3."),
        row(transport_record_id="B5_DM1_02", entity_type="free_payload", entity_id="ADC_PAYLOAD_008", entity_name="DM1", payload_core="DM1", transporter_gene="ABCB1", cell_background="MDCKII-MDR1", transporter_system="human MDR1 transfectant", assay_endpoint="bidirectional apparent permeability; ER=Papp B-A/Papp A-B", concentration_value="0.5", concentration_unit="uM", papp_ab_mean="1.68", papp_ab_sd="0.38", papp_ba_mean="11.17", papp_ba_sd="0.50", efflux_ratio_mean="6.65", efflux_ratio_sd="0.23", background_corrected_efflux_ratio="7.47", source_assigned_call="substrate", evidence_tier="A_direct_numeric", numeric_origin="source_table", truth_set_eligible="yes_continuous_no_binary_cutoff", comparability_group_id="DM1_MDCKII_0.5uM", source_id="SRC_B5_001", source_url=SOURCES[0]["url"], source_location="section 2.6.4.7.2", notes="Mean +/- SD, n=3; net/background-corrected ER 7.47."),
        row(transport_record_id="B5_DM1_03", entity_type="free_payload", entity_id="ADC_PAYLOAD_008", entity_name="DM1", payload_core="DM1", transporter_gene="ABCB1", cell_background="MDCKII-MDR1", transporter_system="human MDR1 transfectant; inhibitor experiment", assay_endpoint="bidirectional apparent permeability; ER=Papp B-A/Papp A-B", concentration_value="0.5", concentration_unit="uM", papp_ab_mean="1.27", papp_ab_sd="0.18", papp_ba_mean="7.39", papp_ba_sd="0.15", efflux_ratio_mean="5.90", efflux_ratio_sd="0.77", source_assigned_call="substrate_no_inhibitor_run", evidence_tier="A_direct_numeric_inhibitor_rescue", numeric_origin="source_table", truth_set_eligible="yes_continuous_no_binary_cutoff", comparability_group_id="DM1_MDCKII_inhibitor_run", source_id="SRC_B5_001", source_url=SOURCES[0]["url"], source_location="section 2.6.4.7.2", notes="Separate inhibitor experiment; do not pool with ER 6.65 run."),
        row(transport_record_id="B5_DM1_04", entity_type="free_payload", entity_id="ADC_PAYLOAD_008", entity_name="DM1", payload_core="DM1", transporter_gene="ABCB1", cell_background="MDCKII-MDR1", transporter_system="human MDR1 transfectant; inhibitor experiment", assay_endpoint="bidirectional apparent permeability; ER=Papp B-A/Papp A-B", concentration_value="0.5", concentration_unit="uM", inhibitor="PSC833", inhibitor_concentration_value="10", inhibitor_concentration_unit="uM", papp_ab_mean="1.46", papp_ab_sd="0.24", papp_ba_mean="1.57", papp_ba_sd="0.10", efflux_ratio_mean="1.10", efflux_ratio_sd="0.20", source_assigned_call="inhibitor_rescue", evidence_tier="A_direct_numeric_inhibitor_rescue", numeric_origin="source_table", truth_set_eligible="yes_continuous_no_binary_cutoff", comparability_group_id="DM1_MDCKII_inhibitor_run", source_id="SRC_B5_001", source_url=SOURCES[0]["url"], source_location="section 2.6.4.7.2", notes="Mean +/- SD, n=3."),
        row(transport_record_id="B5_DM1_05", entity_type="free_payload", entity_id="ADC_PAYLOAD_008", entity_name="DM1", payload_core="DM1", transporter_gene="ABCB1", cell_background="MDCKII-MDR1", transporter_system="human MDR1 transfectant; inhibitor experiment", assay_endpoint="bidirectional apparent permeability; ER=Papp B-A/Papp A-B", concentration_value="0.5", concentration_unit="uM", inhibitor="verapamil", inhibitor_concentration_value="60", inhibitor_concentration_unit="uM", papp_ab_mean="1.90", papp_ab_sd="0.32", papp_ba_mean="3.19", papp_ba_sd="0.83", efflux_ratio_mean="1.76", efflux_ratio_sd="0.76", source_assigned_call="partial_inhibitor_rescue", evidence_tier="A_direct_numeric_inhibitor_rescue", numeric_origin="source_table", truth_set_eligible="yes_continuous_no_binary_cutoff", comparability_group_id="DM1_MDCKII_inhibitor_run", source_id="SRC_B5_001", source_url=SOURCES[0]["url"], source_location="section 2.6.4.7.2", notes="Mean +/- SD, n=3."),
        row(transport_record_id="B5_CYSMC_01", entity_type="released_species", entity_id="ADC_REL_006", entity_name="cys-mcMMAF", payload_core="MMAF", transporter_gene="ABCB1", cell_background="MDCKII human P-gp", transporter_system="human P-gp transfectant", assay_endpoint="bidirectional efflux ratio", concentration_value="100", concentration_unit="nM", efflux_ratio_mean="34", source_assigned_call="substrate", evidence_tier="A_direct_numeric_inhibitor_rescue", numeric_origin="regulatory_summary_numeric", truth_set_eligible="yes_continuous_no_binary_cutoff", comparability_group_id="CYSMCMMAF_MDCKII_100nM", source_id="SRC_B5_009", source_url=SOURCES[8]["url"], source_location="Section 4.5.3, page 32", notes="Papp values not reproduced in the review summary."),
        row(transport_record_id="B5_CYSMC_02", entity_type="released_species", entity_id="ADC_REL_006", entity_name="cys-mcMMAF", payload_core="MMAF", transporter_gene="ABCB1", cell_background="MDCKII human P-gp", transporter_system="human P-gp transfectant", assay_endpoint="bidirectional efflux ratio", concentration_value="100", concentration_unit="nM", inhibitor="verapamil", inhibitor_concentration_value="60", inhibitor_concentration_unit="uM", efflux_ratio_mean="1.1", source_assigned_call="inhibitor_rescue", evidence_tier="A_direct_numeric_inhibitor_rescue", numeric_origin="regulatory_summary_numeric", truth_set_eligible="yes_continuous_no_binary_cutoff", comparability_group_id="CYSMCMMAF_MDCKII_100nM", source_id="SRC_B5_009", source_url=SOURCES[8]["url"], source_location="Section 4.5.3, page 32", notes="Papp values not reproduced in the review summary."),
    ]

    foti_rows = read_csv(PROCESSED / "foti_figure8_digitization_v1.csv")
    for item in foti_rows:
        is_lmd = item["analyte"] == "Lys-MCC-DM1"
        rows.append(
            row(
                transport_record_id=f"B5_{item['digitization_record_id']}",
                entity_type="released_species" if is_lmd else "assay_control",
                entity_id="ADC_REL_001" if is_lmd else "",
                entity_name=item["analyte"],
                payload_core="DM1" if is_lmd else "",
                transporter_gene="ABCB1",
                cell_background=item["cell_system"],
                transporter_system="human P-gp transfectant" if "human P-gp" in item["cell_system"] else "parental control",
                assay_endpoint="efflux ratio (Papp B-A/Papp A-B)",
                concentration_value=item["concentration_value"],
                concentration_unit=item["concentration_unit"],
                inhibitor="elacridar" if item["inhibitor"] else "",
                efflux_ratio_mean=item["efflux_ratio_consensus"],
                source_assigned_call="no_evidence_as_pgp_substrate_low_passive_permeability" if is_lmd else "assay_control",
                evidence_tier="A_direct_numeric_digitized_exploratory" if is_lmd else "control_digitized_exploratory",
                numeric_origin=item["numeric_origin"],
                digitization_status=item["digitization_status"],
                truth_set_eligible="no_pending_independent_digitization",
                comparability_group_id=f"FOTI_F8_{item['figure_group']}",
                source_id="SRC_B5_002",
                source_url=SOURCES[1]["url"],
                source_location="Figure 8",
                notes=f"Bar mean only; one-pixel value uncertainty {item['one_pixel_value_uncertainty']}; error bar not digitized.",
            )
        )
    return rows


GRADE_DEFINITIONS = [
    {"matched_grade": "M0", "definition": "Same exact chemical entities compared in the same transporter cell background, assay endpoint, concentration framework and perturbation design.", "numeric_fold_allowed": "yes_within_prespecified_endpoint", "current_record_count": "0"},
    {"matched_grade": "M1", "definition": "Same exact chemical entity and transporter with independent direct assays, but not the same experimental run or background.", "numeric_fold_allowed": "no_cross_study_fold", "current_record_count": "0"},
    {"matched_grade": "M1H", "definition": "Potential same-entity independent support, held because exact material identity or reproduced direct numeric output is incomplete.", "numeric_fold_allowed": "no", "current_record_count": "2"},
    {"matched_grade": "M2", "definition": "Free payload and released species use the same endpoint family, but cell background, concentration or laboratory differs.", "numeric_fold_allowed": "no", "current_record_count": "1"},
    {"matched_grade": "M3", "definition": "Same cellular/ADC system provides functional or intracellular-retention contrast, but not isolated direct transport for both chemicals.", "numeric_fold_allowed": "no", "current_record_count": "1"},
    {"matched_grade": "M4", "definition": "Cross-assay directional or same-source qualitative comparison only.", "numeric_fold_allowed": "no", "current_record_count": "2"},
]


def matched_comparisons() -> list[dict]:
    return [
        {
            "comparison_id": "MATCH_B5_001", "payload_core": "DM1", "transporter_gene": "ABCB1",
            "free_entity_id": "ADC_PAYLOAD_008", "free_entity_name": "DM1", "released_entity_id": "ADC_REL_001", "released_entity_name": "Lys-MCC-DM1",
            "free_evidence": "PMDA direct MDCKII-MDR1 ER 6.65 at 0.5 uM; ER 5.90 falls to 1.10 with PSC833 and 1.76 with verapamil",
            "released_evidence": "Foti Figure 8 LLC-PK1 human P-gp ER approximately 1.27, 1.20, 0.75 and 0.59 at 0.1, 1, 5 and 10 uM; source concludes no evidence as P-gp substrate",
            "matched_grade": "M2", "endpoint_alignment": "same_bidirectional_efflux_ratio_family", "cell_background_alignment": "different_MDCKII_vs_LLC-PK1", "concentration_alignment": "different_0.5_vs_0.1_to_10_uM",
            "identity_alignment": "distinct_free_payload_vs_lysine_linker_catabolite", "directional_result": "residual lysine-MCC group attenuates P-gp signal in independent direct assays",
            "cross_assay_fold_computed": "no", "training_use": "external_mechanistic_validation_only", "source_ids": "SRC_B5_001 | SRC_B5_002",
            "limitation": "Strongest near-match found, but not M0: different cell background/lab/concentration, and Figure 8 values lack an independent second curator.",
        },
        {
            "comparison_id": "MATCH_B5_002", "payload_core": "MMAF", "transporter_gene": "ABCB1",
            "free_entity_id": "ADC_PAYLOAD_018", "free_entity_name": "MMAF", "released_entity_id": "ADC_REL_006", "released_entity_name": "cys-mcMMAF",
            "free_evidence": "NCI/ADR-RES patent experiment: MMAF IC50 0.45 to 0.07 with verapamil, fold 6.4; Roth HEK293 RR 2.0",
            "released_evidence": "PMDA direct MDCKII P-gp ER 34 at 100 nM to 1.1 with verapamil 60 uM",
            "matched_grade": "M4", "endpoint_alignment": "functional_cytotoxicity_vs_direct_transport", "cell_background_alignment": "different_NCI-ADR-RES_vs_MDCKII", "concentration_alignment": "different_and_patent_verapamil_conflict",
            "identity_alignment": "distinct_free_payload_vs_cysteine_linker_catabolite", "directional_result": "both show P-gp susceptibility, with stronger direct signal for cys-mcMMAF",
            "cross_assay_fold_computed": "no", "training_use": "hypothesis_support_only", "source_ids": "SRC_B5_006 | SRC_B5_009 | Roth 2026 Table 3",
            "limitation": "Patent methods text says 40 uM verapamil whereas Table 7 says 10 uM; endpoints and systems are incompatible.",
        },
        {
            "comparison_id": "MATCH_B5_003", "payload_core": "DM1", "transporter_gene": "ABCB1",
            "free_entity_id": "", "free_entity_name": "Lys-SMCC-DM1 from ADC catabolism", "released_entity_id": "ADC_REL_002", "released_entity_name": "Lys-PEG4Mal-DM1 from ADC catabolism",
            "free_evidence": "Kovtun same-cell experiments show lower intracellular metabolite and greater cyclosporin sensitivity for SMCC-DM1 ADC",
            "released_evidence": "PEG4Mal-DM1 ADC is 7-fold more potent in UO-31 and 5-fold more potent in COLO205MDR; metabolite retention is higher",
            "matched_grade": "M3", "endpoint_alignment": "same_cell_functional_and_intracellular_metabolite", "cell_background_alignment": "same_within_each_comparison", "concentration_alignment": "same_within_reported_comparison",
            "identity_alignment": "different_catabolites_generated_by_matched_ADC_formats", "directional_result": "PEG4Mal residual group reduces MDR1-associated loss relative to SMCC catabolite",
            "cross_assay_fold_computed": "no", "training_use": "mechanistic_context_only", "source_ids": "SRC_B5_003",
            "limitation": "The study compares processed ADC formats; it does not report isolated bidirectional transport for both catabolites.",
        },
        {
            "comparison_id": "MATCH_B5_004", "payload_core": "DM4", "transporter_gene": "ABCB1",
            "free_entity_id": "ADC_PAYLOAD_009", "free_entity_name": "DM4", "released_entity_id": "ADC_REL_004", "released_entity_name": "S-methyl-DM4",
            "free_evidence": "FDA review calls DM4 an MDR1/P-gp substrate", "released_evidence": "same FDA review calls S-methyl-DM4 an MDR1/P-gp substrate",
            "matched_grade": "M4", "endpoint_alignment": "same_source_qualitative_call", "cell_background_alignment": "not_reported", "concentration_alignment": "not_reported",
            "identity_alignment": "distinct_free_payload_vs_S-methyl_catabolite", "directional_result": "same positive direction", "cross_assay_fold_computed": "no",
            "training_use": "call_only_external_support", "source_ids": "SRC_B5_008", "limitation": "No numeric transport data are reproduced.",
        },
        {
            "comparison_id": "MATCH_B5_005", "payload_core": "DXd", "transporter_gene": "ABCB1",
            "free_entity_id": "ADC_PAYLOAD_029", "free_entity_name": "Deruxtecan (Roth material name)", "released_entity_id": "ADC_REL_010", "released_entity_name": "DXd / MAAA-1181a",
            "free_evidence": "Roth HEK293 MDR19 cytotoxicity RR 20", "released_evidence": "FDA Enhertu review identifies MAAA-1181a as a P-gp substrate",
            "matched_grade": "M1H", "endpoint_alignment": "functional_RR_vs_regulatory_substrate_call", "cell_background_alignment": "not_same_or_not_reproduced", "concentration_alignment": "not_same_or_not_reproduced",
            "identity_alignment": "likely_same_DXd_core_but_Roth_exact_test_material_unresolved", "directional_result": "same positive direction", "cross_assay_fold_computed": "no",
            "training_use": "independent_directional_validation_with_identity_hold", "source_ids": "SRC_B5_004 | Roth 2026 Table 3", "limitation": "Roth uses the name Deruxtecan; exact tested structure/material must be confirmed before treating this as an exact-entity replicate.",
        },
        {
            "comparison_id": "MATCH_B5_006", "payload_core": "DXd", "transporter_gene": "ABCG2",
            "free_entity_id": "ADC_PAYLOAD_029", "free_entity_name": "Deruxtecan (Roth material name)", "released_entity_id": "ADC_REL_010", "released_entity_name": "DXd / MAAA-1181a",
            "free_evidence": "Roth HEK293 R5 cytotoxicity RR 240", "released_evidence": "FDA Enhertu review identifies MAAA-1181a as a BCRP substrate",
            "matched_grade": "M1H", "endpoint_alignment": "functional_RR_vs_regulatory_substrate_call", "cell_background_alignment": "not_same_or_not_reproduced", "concentration_alignment": "not_same_or_not_reproduced",
            "identity_alignment": "likely_same_DXd_core_but_Roth_exact_test_material_unresolved", "directional_result": "same strong positive direction", "cross_assay_fold_computed": "no",
            "training_use": "independent_directional_validation_with_identity_hold", "source_ids": "SRC_B5_004 | Roth 2026 Table 3", "limitation": "No numeric BCRP transport value is reproduced in the FDA review; Roth material identity remains unresolved.",
        },
    ]


def mechanistic_rows() -> list[dict]:
    base = {"comparison_family": "Kovtun_2010_SMCC_vs_PEG4Mal_DM1", "transporter_gene": "ABCB1", "source_id": "SRC_B5_003", "source_url": SOURCES[2]["url"], "evidence_class": "M3_same_system_mechanistic", "training_use": "no_mechanistic_context_only"}
    entries = [
        ("MECH_B5_001", "UO-31", "cytotoxicity potency", "anti-EGFR-SMCC-DM1", "reference", "anti-EGFR-PEG4Mal-DM1", "7-fold more potent", "Results/Figure 4", "Cyclosporin restored SMCC-DM1 activity to the PEG4Mal-DM1 level and had little effect on PEG4Mal-DM1."),
        ("MECH_B5_002", "COLO205MDR", "cytotoxicity potency", "anti-CanAg-SMCC-DM1", "reference", "anti-CanAg-PEG4Mal-DM1", "5-fold more potent", "Results/Figure 5", "Matched ADC-format functional comparison in MDR1-expressing cells."),
        ("MECH_B5_003", "COLO205MDR", "intracellular catabolite amount", "Lys-SMCC-DM1", "approximately 50% lower", "Lys-PEG4Mal-DM1", "reference", "Results/Figure 6", "Parental COLO205 cells showed similar intracellular metabolite amounts."),
        ("MECH_B5_004", "COLO205MDR", "G2/M fraction", "anti-EpCAM-SMCC-DM1", "37", "anti-EpCAM-PEG4Mal-DM1", "67", "Results", "Percent of cells in G2/M; untreated 15%."),
        ("MECH_B5_005", "COLO205 parental", "G2/M fraction", "anti-EpCAM-SMCC-DM1", "85", "anti-EpCAM-PEG4Mal-DM1", "88", "Results", "Percent of cells in G2/M; untreated 27%."),
        ("MECH_B5_006", "COLO205MDR/HCT-15/UO-31", "MDR1 inhibitor response", "SMCC-DM1 ADC", "significant cyclosporin effect", "PEG4Mal-DM1 ADC", "marginal/no material effect", "Figures 3-5", "Supports differential MDR1 susceptibility of the generated catabolites; no isolated transport rate."),
    ]
    rows = []
    for rid, cell, endpoint, a, av, b, bv, loc, note in entries:
        rows.append({"mechanistic_record_id": rid, **base, "cell_system": cell, "endpoint": endpoint, "comparator_a": a, "comparator_a_value": av, "comparator_b": b, "comparator_b_value": bv, "source_location": loc, "note": note})
    return rows


def patent_functional_rows() -> list[dict]:
    values = [
        ("PAT_B5_001", "free_drug", "Doxorubicin", "10", "0.56", "17.9", "mM"),
        ("PAT_B5_002", "free_drug", "Cisplatin", "2.9", "3.3", "0.88", "mM"),
        ("PAT_B5_003", "free_payload", "MMAF", "0.45", "0.07", "6.4", "mM"),
        ("PAT_B5_004", "adc", "cBR96-vcMMAF", "1.9", "0.93", "2", "mg/mL"),
        ("PAT_B5_005", "adc", "mOKT9-vcMMAF", "0.27", "0.12", "2.2", "mg/mL"),
    ]
    return [
        {
            "functional_record_id": rid, "entity_type": entity_type, "entity_name": name,
            "cell_system": "NCI/ADR-RES", "assay_endpoint": "96 h resazurin IC50",
            "ic50_no_verapamil": no_ver, "ic50_with_verapamil": plus_ver, "ic50_unit": unit,
            "fold_difference": fold, "fold_definition": "IC50(no verapamil) / IC50(+verapamil)",
            "verapamil_concentration_table_header": "10 uM", "verapamil_concentration_methods_text": "40 uM",
            "concentration_conflict_status": "internal_source_conflict_unresolved", "evidence_tier": "P_patent_functional",
            "training_use": "no_functional_context_only", "source_id": "SRC_B5_006", "source_location": "Table 7; PDF page 37 (printed pages 61-62)",
            "source_url": SOURCES[5]["url"],
        }
        for rid, entity_type, name, no_ver, plus_ver, fold, unit in values
    ]


def release_transport_v3() -> list[dict]:
    rows = read_csv(PROCESSED / "adc_release_species_transporter_evidence_long_v2.csv")
    fields = list(rows[0].keys())

    def blank(**kwargs):
        return {field: kwargs.get(field, "") for field in fields}

    rows.extend([
        blank(transport_record_id="ADC_TRANS_B5_001", release_species_id="ADC_REL_010", transporter_gene="ABCB1", transporter_protein="P-gp", species="Homo sapiens", assay_system="human biomaterial transporter studies summarized by FDA", test_concentrations="not reproduced", source_assigned_call="substrate", source_binary_label="1", model_label_status="not_assigned_regulatory_call_only", evidence_tier="B_regulatory_call_no_numeric", truth_set_candidate="yes_call_only", source_id="SRC_B5_004", source_location="PDF page 85 (printed page 84)", curation_note="No numeric P-gp transport result is reproduced; keep separate from Roth functional RR.", release_species_name="DXd", source_url=SOURCES[3]["url"], assay_endpoint="regulatory substrate call"),
        blank(transport_record_id="ADC_TRANS_B5_002", release_species_id="ADC_REL_010", transporter_gene="ABCG2", transporter_protein="BCRP", species="Homo sapiens", assay_system="human biomaterial transporter studies summarized by FDA", test_concentrations="not reproduced", source_assigned_call="substrate", source_binary_label="1", model_label_status="not_assigned_regulatory_call_only", evidence_tier="B_regulatory_call_no_numeric", truth_set_candidate="yes_call_only", source_id="SRC_B5_004", source_location="PDF page 85 (printed page 84)", curation_note="No numeric BCRP transport result is reproduced; keep separate from Roth functional RR.", release_species_name="DXd", source_url=SOURCES[3]["url"], assay_endpoint="regulatory substrate call"),
    ])
    for item in read_csv(PROCESSED / "foti_figure8_digitization_v1.csv"):
        if item["analyte"] != "Lys-MCC-DM1" or "human P-gp" not in item["cell_system"]:
            continue
        rows.append(blank(
            transport_record_id=f"ADC_TRANS_B5_{int(item['digitization_record_id'][-2:]) + 2:03d}", release_species_id="ADC_REL_001", transporter_gene="ABCB1", transporter_protein="P-gp", species="Homo sapiens",
            assay_system="human MDR1-transfected LLC-PK1 bidirectional permeability", test_concentrations=f"{item['concentration_value']} {item['concentration_unit']}",
            efflux_ratio_no_inhibitor=item["efflux_ratio_consensus"] if not item["inhibitor"] else "",
            efflux_ratio_with_inhibitor=item["efflux_ratio_consensus"] if item["inhibitor"] else "", inhibitor=item["inhibitor"], source_assigned_call="no_evidence_as_pgp_substrate_low_potential", source_binary_label="0",
            model_label_status="hold_until_independent_second_digitization", evidence_tier="A_direct_numeric_digitized_exploratory", truth_set_candidate="pending_independent_curator", source_id="SRC_B5_002", source_location="Figure 8",
            curation_note="Algorithmic mean-bar digitization only; error bars not digitized; no independent second curator.", release_species_name="Lys-MCC-DM1", source_doi="10.1124/dmd.115.064576", source_url=SOURCES[1]["url"], assay_endpoint="efflux ratio (Papp B-A/Papp A-B)",
        ))
    return rows


def pair_evidence_v2() -> list[dict]:
    rows = read_csv(PROCESSED / "free_payload_release_pair_evidence_v1.csv")
    for row in rows:
        if row["pair_id"] == "ADC_PAIR_002":
            row["free_endpoint"] = "MDCKII-MDR1 bidirectional efflux ratio with inhibitor rescue"
            row["free_effect_value"] = "6.65 at 0.5 uM; 5.90 to 1.10 PSC833 or 1.76 verapamil in separate run"
            row["free_effect_definition"] = "Papp B-to-A / Papp A-to-B"
            row["release_endpoint"] = "LLC-PK1 human P-gp bidirectional efflux ratio; algorithmically digitized"
            row["release_effect_value"] = "1.27 | 1.20 | 0.75 | 0.59 at 0.1 | 1 | 5 | 10 uM; 0.46 at 5 uM + elacridar"
            row["release_effect_definition"] = "Papp B-to-A / Papp A-to-B; exploratory bar means"
            row["directional_comparison"] = "discordant_direct_assay_direction_free_DM1_positive_LysMCCDM1_low"
            row["endpoint_compatibility"] = "same_endpoint_family_different_cell_background_lab_and_concentration"
            row["free_source_id"] = "SRC_B5_001"
            row["release_source_ids"] = "SRC_B5_002"
            row["comparison_note"] = "Upgraded from cytotoxicity-vs-qualitative to direct bidirectional assays on both entities, but remains M2 rather than M0. No cross-study fold change; Foti values await independent digitization."
        if row["pair_id"] == "ADC_PAIR_001":
            row["free_endpoint"] = "NCI/ADR-RES inhibitor-rescue cytotoxicity plus Roth HEK293 transporter/control RR"
            row["free_effect_value"] = "patent fold 6.4; Roth RR 2.0"
            row["free_effect_definition"] = "IC50 no inhibitor / +verapamil; IC50(MDR19) / IC50(pcDNA)"
            row["free_source_id"] = "SRC_B5_006 | Roth 2026 Table 3"
            row["comparison_note"] = "Free MMAF now has a second positive functional-rescue source, but patent Table 7/methods disagree on verapamil concentration and the released-species evidence is still a different direct assay."

    fields = list(rows[0].keys())
    additions = [
        {
            "pair_id": "ADC_PAIR_008", "payload_core": "DXd", "free_payload_compound_id": "ADC_PAYLOAD_029", "free_payload_name": "Deruxtecan (Roth material name)", "release_species_id": "ADC_REL_010", "release_species_name": "DXd / MAAA-1181a", "transporter_gene": "ABCB1",
            "free_endpoint": "HEK293 transporter/control cytotoxicity RR", "free_effect_value": "20", "free_effect_definition": "IC50(MDR19) / IC50(pcDNA)", "release_endpoint": "FDA regulatory substrate call", "release_effect_value": "P-gp substrate; numeric transport value not reproduced", "release_effect_definition": "source call",
            "directional_comparison": "same_positive_direction_identity_pending", "endpoint_compatibility": "not_comparable_cross_assay", "causal_interpretation": "independent_directional_validation_only", "free_source_id": "Roth 2026 Table 3", "release_source_ids": "SRC_B5_004", "comparison_note": "Likely same released DXd core, but exact Roth test material remains unresolved; do not merge the C52 Deruxtecan registry structure with C26 DXd without source confirmation.",
        },
        {
            "pair_id": "ADC_PAIR_009", "payload_core": "DXd", "free_payload_compound_id": "ADC_PAYLOAD_029", "free_payload_name": "Deruxtecan (Roth material name)", "release_species_id": "ADC_REL_010", "release_species_name": "DXd / MAAA-1181a", "transporter_gene": "ABCG2",
            "free_endpoint": "HEK293 transporter/control cytotoxicity RR", "free_effect_value": "240", "free_effect_definition": "IC50(R5) / IC50(pcDNA)", "release_endpoint": "FDA regulatory substrate call", "release_effect_value": "BCRP substrate; numeric transport value not reproduced", "release_effect_definition": "source call",
            "directional_comparison": "same_strong_positive_direction_identity_pending", "endpoint_compatibility": "not_comparable_cross_assay", "causal_interpretation": "independent_directional_validation_only", "free_source_id": "Roth 2026 Table 3", "release_source_ids": "SRC_B5_004", "comparison_note": "Strong independent directional agreement, but no reproduced FDA numeric value and exact Roth material identity remains unresolved.",
        },
    ]
    for addition in additions:
        rows.append({field: addition.get(field, "") for field in fields})
    return rows


def raw_log() -> list[dict]:
    specs = [
        (RAW5 / "PMDA_Kadcyla_2013_Nonclinical_PK.pdf", SOURCES[0]["url"], "SRC_B5_001"),
        (ROOT / "raw" / "second_batch_release_species" / "Foti_2015_Figure8_LysMCCDM1_Pgp.jpg", SOURCES[1]["url"], "SRC_B5_002"),
        (ROOT / "raw" / "second_batch_release_species" / "Kovtun_2010_article_copy.pdf", SOURCES[2]["url"], "SRC_B5_003"),
        (RAW5 / "FDA_Enhertu_2019_MultidisciplineReview.pdf", SOURCES[3]["url"], "SRC_B5_004"),
        (RAW5 / "PMDA_Enhertu_2020_Pharmacokinetics.pdf", SOURCES[4]["url"], "SRC_B5_005"),
        (RAW5 / "US7750116B1_Auristatin_Patent.pdf", "https://patentimages.storage.googleapis.com/e2/0c/11/3d27c2b2f98574/US7750116.pdf", "SRC_B5_006"),
        (RAW5 / "FDA_Kadcyla_2013_Pharmacology_Review.pdf", SOURCES[6]["url"], "SRC_B5_007"),
        (ROOT / "raw" / "third_batch_pair_resolution" / "FDA_Elahere_2022_MultidisciplineReview.pdf", SOURCES[7]["url"], "SRC_B5_008"),
        (ROOT / "raw" / "second_batch_release_species" / "PMDA_belantamab_mafodotin_English_report.pdf", SOURCES[8]["url"], "SRC_B5_009"),
    ]
    rows = []
    for path, url, source_id in specs:
        data = path.read_bytes()
        valid = data.startswith(b"%PDF") if path.suffix.lower() == ".pdf" else data.startswith(b"\xff\xd8")
        rows.append({
            "filename": str(path.relative_to(ROOT)), "source_id": source_id, "source_url": url,
            "download_status": "reused_existing_validated_signature" if "fifth_batch" not in str(path) else "downloaded_validated_signature",
            "signature_valid": "yes" if valid else "no", "size_bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(), "acquired_or_verified_at": "2026-07-28",
        })
    failed = ROOT / "tmp" / "failed_downloads" / "fifth_batch" / "Foti_2015_PMC11024933_BioC_error.html"
    if failed.exists():
        data = failed.read_bytes()
        rows.append({"filename": str(failed.relative_to(ROOT)), "source_id": "SRC_B5_002", "source_url": "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC11024933/unicode", "download_status": "failed_endpoint_html_error_excluded_from_raw_sources", "signature_valid": "no", "size_bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(), "acquired_or_verified_at": "2026-07-28"})
    return rows


def main() -> None:
    direct = direct_numeric_rows()
    matched = matched_comparisons()
    mechanistic = mechanistic_rows()
    patent_rows = patent_functional_rows()
    release_v3 = release_transport_v3()
    pairs_v2 = pair_evidence_v2()
    downloads = raw_log()

    write_csv(PROCESSED / "direct_bidirectional_transport_numeric_v1.csv", direct)
    write_csv(PROCESSED / "matched_context_comparisons_v1.csv", matched)
    write_csv(PROCESSED / "matched_evidence_grade_definitions_v1.csv", GRADE_DEFINITIONS)
    write_csv(PROCESSED / "paired_release_species_mechanistic_evidence_v1.csv", mechanistic)
    write_csv(PROCESSED / "patent_mmaf_verapamil_functional_v1.csv", patent_rows)
    write_csv(PROCESSED / "adc_release_species_transporter_evidence_long_v3.csv", release_v3)
    write_csv(PROCESSED / "free_payload_release_pair_evidence_v2.csv", pairs_v2)
    write_csv(PROCESSED / "fifth_batch_source_registry.csv", SOURCES)
    write_csv(PROCESSED / "fifth_batch_raw_download_log.csv", downloads)

    grade_counts = Counter(row["matched_grade"] for row in matched)
    summary = {
        "generated_at": "2026-07-28",
        "fifth_batch_sources": len(SOURCES),
        "raw_source_files_or_attempts": len(downloads),
        "direct_numeric_rows": len(direct),
        "direct_numeric_target_rows": sum(row["entity_type"] != "assay_control" for row in direct),
        "foti_digitized_bars": 16,
        "foti_digitization_truth_set_eligible": 0,
        "matched_comparisons": len(matched),
        "matched_grade_counts": dict(sorted(grade_counts.items())),
        "m0_exact_same_assay_pairs": grade_counts.get("M0", 0),
        "same_system_mechanistic_rows": len(mechanistic),
        "patent_functional_rows": len(patent_rows),
        "release_transporter_records_v3": len(release_v3),
        "free_release_pairs_v2": len(pairs_v2),
        "model_binary_labels_created": 0,
        "key_results": [
            "Free DM1 now has direct numeric MDCKII-MDR1 transport with inhibitor rescue: ER 6.65; separate run 5.90 to 1.10 with PSC833 and 1.76 with verapamil.",
            "Lys-MCC-DM1 Figure 8 bar means were reproducibly digitized, supporting low P-gp signal, but remain outside the truth set pending an independent curator.",
            "DXd/MAAA-1181a gained FDA P-gp and BCRP substrate calls, directionally agreeing with Roth RR 20 and 240, subject to exact Roth material identity resolution.",
            "No M0 free-payload/released-species pair in the same cell background, endpoint and concentration framework was found.",
        ],
        "feasibility_conclusion": "Evidence quality improved and supports an ADC-specific evidence-layer predictor, but the matched validation set is still too small and heterogeneous for a credible universal substrate classifier. Continue with endpoint-stratified prediction plus abstention, not binary label transfer.",
    }
    (PROCESSED / "fifth_batch_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
