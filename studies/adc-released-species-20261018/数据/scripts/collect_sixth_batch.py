#!/usr/bin/env python3
"""Build sixth-batch ADC transporter evidence and identity-audit tables.

The batch adds regulatory numeric Caco-2 data for DXd/MAAA-1181a,
audits the unresolved identity of the Roth "Deruxtecan" reagent, and
recalculates Roth external-summary metrics with that ambiguous record held out.
No new binary substrate cutoff is introduced.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import joblib
import numpy as np
from scipy.stats import spearmanr

import run_standardization_leakage_baseline as baseline


DATA_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = DATA_ROOT.parent
PROCESSED = DATA_ROOT / "processed"
RAW6 = DATA_ROOT / "raw" / "sixth_batch_primary_tables"
MODELS = PROJECT_ROOT / "模型" / "20260728_baseline"
RUN_DATE = "2026-07-28"

PMDA_URL = "https://www.pmda.go.jp/files/000238706.pdf"
ROTH_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC12883344/"
FOTI_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11024933/"
ROCK_URL = "https://doi.org/10.1124/dmd.115.064253"
PUBCHEM_DERUXTECAN_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/118305111"
PUBCHEM_DXD_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/117888634"


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        if not rows:
            raise ValueError(f"Cannot infer fields for empty table: {path}")
        fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


SOURCES = [
    {
        "source_id": "SRC_B6_001",
        "title": "Enhertu PMDA English review report",
        "year": 2020,
        "source_type": "regulatory_review_numeric",
        "url": PMDA_URL,
        "use_in_batch": "direct numeric Caco-2 transport of 14C-MAAA-1181a with P-gp, BCRP and dual inhibitors",
        "source_location": "PDF printed pages 25-26; section 4.6.3 Transporters",
        "primary_or_secondary": "primary_regulatory_summary",
        "notes": "At 1 uM, ER was 10.7 without inhibitor, 4.37 with verapamil 100 uM, 4.44 with novobiocin 10 uM, and 4.43 with GF120918 10 uM. The report assigns P-gp and BCRP substrate calls but the shared Caco-2 system does not isolate a gene-specific ER.",
    },
    {
        "source_id": "SRC_B6_002",
        "title": "Roth et al. ADC payload transporter study",
        "year": 2026,
        "source_type": "primary_article_and_supplement",
        "url": ROTH_URL,
        "use_in_batch": "audit exact identity and provenance of the Table 3 Deruxtecan reagent",
        "source_location": "Methods; Tables 1-3; Supplementary Materials",
        "primary_or_secondary": "primary_article",
        "notes": "Deruxtecan is absent from the vendor-bearing 27-compound Table 2 and appears only among additional confirmatory compounds in Table 3. No supplier, catalog number, CAS, structure or batch identifier is reported in the main text or supplement.",
    },
    {
        "source_id": "SRC_B6_003",
        "title": "PubChem CID 118305111 Deruxtecan pyrrole-2,5-dione",
        "year": 2026,
        "source_type": "chemical_structure_registry",
        "url": PUBCHEM_DERUXTECAN_URL,
        "use_in_batch": "define the C52 linker-payload structure returned by the exact name Deruxtecan",
        "source_location": "PUG REST property record checked 2026-07-28",
        "primary_or_secondary": "derived_registry",
        "notes": "C52H56FN9O13, MW 1034.1, InChIKey WXNSCLIZKHLNSG-MCZRLCSDSA-N. This is not the C26 released DXd payload.",
    },
    {
        "source_id": "SRC_B6_004",
        "title": "PubChem CID 117888634 DXd",
        "year": 2026,
        "source_type": "chemical_structure_registry",
        "url": PUBCHEM_DXD_URL,
        "use_in_batch": "define the C26 released DXd/MAAA-1181a structure",
        "source_location": "PUG REST property record checked 2026-07-28",
        "primary_or_secondary": "derived_registry",
        "notes": "C26H24FN3O6, MW 493.5, InChIKey PLXLYXLUCNZSAA-QLXKLKPCSA-N.",
    },
    {
        "source_id": "SRC_B6_005",
        "title": "Intracellular Catabolism of an Antibody Drug Conjugate with a Noncleavable Linker",
        "year": 2015,
        "source_type": "primary_article_with_unretrieved_supplement",
        "url": ROCK_URL,
        "use_in_batch": "trace the origin of Lys-MCC-DM1 low-permeability evidence",
        "source_location": "Discussion and cited Supplemental Table 2",
        "primary_or_secondary": "primary_article",
        "notes": "The main article states that Lys-MCC-DM1 has low permeability in LLC-PK1 cells and points to Supplemental Table 2. The main paper is a catabolism study; the supplemental numeric table was not retrieved from the current publisher endpoint.",
    },
    {
        "source_id": "SRC_B6_006",
        "title": "Target-Site Drug Metabolism and Transport",
        "year": 2015,
        "source_type": "primary_symposium_article_figure",
        "url": FOTI_URL,
        "use_in_batch": "cross-reference Lys-MCC-DM1 LLC-PK1/MDR1 Figure 8 and its source conclusion",
        "source_location": "Figure 8 and accompanying methods/text",
        "primary_or_secondary": "primary_experimental_report",
        "notes": "Figure 8 reports the direct bidirectional experiment and concludes no evidence that Lys-MCC-DM1 is a P-gp substrate. Exact values remain algorithmically digitized rather than author-tabulated in the available source.",
    },
    {
        "source_id": "SRC_B6_007",
        "title": "PubChem BioAssay AID availability recheck",
        "year": 2026,
        "source_type": "official_database_access_audit",
        "url": "https://pubchem.ncbi.nlm.nih.gov/bioassay/2061167",
        "use_in_batch": "recheck Roth-declared project and related assay downloads",
        "source_location": "PUG REST CSV endpoints checked 2026-07-28",
        "primary_or_secondary": "official_database",
        "notes": "Project record/summary remains available for AID 2061167, but CSV for AID 2061167 and all 13 linked/declaration AIDs still returned HTTP 404.",
    },
]


def dxd_numeric_rows() -> list[dict]:
    conditions = [
        ("B6_DXD_01", "none", "", 10.7, 1.0, "P-gp and BCRP active; residual or additional directional asymmetry not deconvolved"),
        ("B6_DXD_02", "verapamil", 100, 4.37, 10.7 / 4.37, "P-gp inhibited; remaining ER cannot be assigned exclusively to BCRP"),
        ("B6_DXD_03", "novobiocin", 10, 4.44, 10.7 / 4.44, "BCRP inhibited; remaining ER cannot be assigned exclusively to P-gp"),
        ("B6_DXD_04", "GF120918", 10, 4.43, 10.7 / 4.43, "dual P-gp/BCRP inhibition; ER remains materially above 1"),
    ]
    rows = []
    for record_id, inhibitor, inhibitor_um, er, reduction, attribution in conditions:
        rows.append(
            {
                "record_id": record_id,
                "entity_id": "ADC_REL_010",
                "entity_name": "DXd / MAAA-1181a",
                "radiolabel": "14C-MAAA-1181a",
                "cell_system": "Caco-2",
                "assay_endpoint": "efflux ratio=Papp B-to-A / Papp A-to-B",
                "concentration_uM": 1,
                "inhibitor": inhibitor,
                "inhibitor_concentration_uM": inhibitor_um,
                "efflux_ratio": er,
                "baseline_to_condition_er_reduction_fold": round(reduction, 4),
                "source_assigned_call": "P-gp substrate; BCRP substrate",
                "target_attribution": attribution,
                "truth_set_use": "continuous_multi_transporter_external_evidence_not_gene_specific_training_truth",
                "source_id": "SRC_B6_001",
                "source_location": "PMDA review printed pages 25-26; section 4.6.3",
                "source_url": PMDA_URL,
                "notes": "Source reports ER only; individual Papp A-to-B and B-to-A values are not reproduced.",
            }
        )
    return rows


def direct_numeric_v2(dxd_rows: list[dict]) -> list[dict]:
    rows = read_csv(PROCESSED / "direct_bidirectional_transport_numeric_v1.csv")
    for source in dxd_rows:
        rows.append(
            {
                "transport_record_id": source["record_id"],
                "entity_type": "release_species",
                "entity_id": source["entity_id"],
                "entity_name": source["entity_name"],
                "payload_core": "DXd",
                "transporter_gene": "ABCB1|ABCG2",
                "cell_background": source["cell_system"],
                "transporter_system": "endogenous multi-transporter Caco-2 inhibitor panel",
                "assay_endpoint": source["assay_endpoint"],
                "concentration_value": source["concentration_uM"],
                "concentration_unit": "uM",
                "inhibitor": source["inhibitor"],
                "inhibitor_concentration_value": source["inhibitor_concentration_uM"],
                "inhibitor_concentration_unit": "uM" if source["inhibitor_concentration_uM"] != "" else "",
                "papp_ab_mean": "",
                "papp_ab_sd": "",
                "papp_ba_mean": "",
                "papp_ba_sd": "",
                "efflux_ratio_mean": source["efflux_ratio"],
                "efflux_ratio_sd": "",
                "background_corrected_efflux_ratio": "",
                "source_assigned_call": source["source_assigned_call"],
                "evidence_tier": "A_direct_numeric_multi_transporter_inhibitor_panel",
                "numeric_origin": "regulatory_review_reproduced_numeric",
                "digitization_status": "",
                "truth_set_eligible": "no_gene_specific_deconvolution; yes_multi_transporter_continuous_audit",
                "comparability_group_id": "DXD_CACO2_1uM_PMDA",
                "source_id": source["source_id"],
                "source_url": source["source_url"],
                "source_location": source["source_location"],
                "notes": source["target_attribution"],
            }
        )
    return rows


def release_transport_v4() -> list[dict]:
    rows = read_csv(PROCESSED / "adc_release_species_transporter_evidence_long_v3.csv")
    for gene, protein in (("ABCB1", "P-gp"), ("ABCG2", "BCRP")):
        rows.append(
            {
                "transport_record_id": f"ADC_TRANS_B6_DXD_{gene}",
                "release_species_id": "ADC_REL_010",
                "transporter_gene": gene,
                "transporter_protein": protein,
                "species": "Homo sapiens",
                "assay_system": "Caco-2 bidirectional transport with P-gp/BCRP inhibitor panel",
                "test_concentrations": "1 uM 14C-MAAA-1181a",
                "efflux_ratio_no_inhibitor": "10.7",
                "efflux_ratio_with_inhibitor": "4.37 verapamil | 4.44 novobiocin | 4.43 GF120918",
                "inhibitor": "verapamil 100 uM | novobiocin 10 uM | GF120918 10 uM",
                "source_assigned_call": "substrate",
                "source_binary_label": "1",
                "model_binary_label": "",
                "model_label_status": "not_assigned_multi_transporter_numeric_not_gene_deconvolved",
                "evidence_tier": "A_direct_numeric_multi_transporter_inhibitor_panel",
                "truth_set_candidate": "external_directional_validation_only",
                "source_id": "SRC_B6_001",
                "source_location": "section 4.6.3; printed pages 25-26",
                "curation_note": f"PMDA assigns {protein} substrate status, but the numeric Caco-2 ER is shared across P-gp/BCRP and remains 4.43 under dual inhibition; do not use 10.7 as a gene-specific ER.",
                "release_species_name": "DXd / MAAA-1181a",
                "source_doi": "",
                "source_url": PMDA_URL,
                "assay_endpoint": "efflux ratio=Papp B-to-A / Papp A-to-B",
                "control_condition": "no inhibitor",
                "control_value": "10.7",
                "transporter_condition": "endogenous P-gp/BCRP Caco-2",
                "transporter_value": "10.7",
                "inhibitor_condition": "single and dual inhibitors",
                "inhibitor_value": "4.37 | 4.44 | 4.43",
                "effect_ratio": "",
                "effect_ratio_definition": "not computed per transporter",
            }
        )
    return rows


def matched_comparisons_v2() -> list[dict]:
    rows = read_csv(PROCESSED / "matched_context_comparisons_v1.csv")
    numeric = "PMDA Caco-2 direct ER 10.7 at 1 uM; ER 4.37 with verapamil, 4.44 with novobiocin and 4.43 with dual inhibitor GF120918"
    for row in rows:
        if row["payload_core"] == "DXd":
            row["released_evidence"] = numeric
            row["endpoint_alignment"] = "functional_RR_vs_direct_bidirectional_multi_transporter_transport"
            row["cell_background_alignment"] = "HEK293_single_gene_transfectant_vs_Caco2_multi_transporter"
            row["concentration_alignment"] = "Roth_concentration_response_vs_PMDA_1uM"
            row["identity_alignment"] = "Roth_test_material_unresolved_C26_DXd_vs_C52_linker_payload"
            row["cross_assay_fold_computed"] = "no"
            row["training_use"] = "independent_directional_validation_only_identity_and_deconvolution_hold"
            row["source_ids"] = "SRC_B6_001 | SRC_B6_002"
            row["limitation"] = "Direct regulatory numbers are now available, but Caco-2 does not isolate a gene-specific ER and Roth omits the additional Deruxtecan reagent source/structure. Keep M1H; do not compute cross-assay folds."
    return rows


def pair_evidence_v3() -> list[dict]:
    rows = read_csv(PROCESSED / "free_payload_release_pair_evidence_v2.csv")
    for row in rows:
        if row["pair_id"] in {"ADC_PAIR_008", "ADC_PAIR_009"}:
            row["release_endpoint"] = "Caco-2 direct bidirectional transport with P-gp/BCRP inhibitor panel"
            row["release_effect_value"] = "ER 10.7; 4.37 verapamil; 4.44 novobiocin; 4.43 GF120918"
            row["release_effect_definition"] = "Papp B-to-A / Papp A-to-B at 1 uM; shared multi-transporter experiment"
            row["endpoint_compatibility"] = "not_comparable_cross_assay_and_not_gene_specific_deconvolved"
            row["causal_interpretation"] = "independent_directional_validation_only"
            row["release_source_ids"] = "SRC_B6_001"
            row["comparison_note"] = "Regulatory numeric transport strengthens the direction match, but Roth reagent identity is unreported and the Caco-2 ER cannot be assigned separately to P-gp or BCRP."
    return rows


def identity_audit_rows() -> list[dict]:
    return [
        {
            "audit_id": "B6_IDENT_001",
            "evidence_object": "Roth Table 1 terminology",
            "reported_name_or_identifier": "Dxd (exatecan derivative)",
            "supports_c26_dxd": "yes_semantic",
            "supports_c52_linker_payload": "no",
            "supplier_or_catalog_reported": "no",
            "exact_structure_reported": "no",
            "finding": "The paper conceptually distinguishes the ADC payload as DXd.",
            "model_action": "context_only_not_material_identity_proof",
            "source_id": "SRC_B6_002",
            "source_location": "Table 1",
        },
        {
            "audit_id": "B6_IDENT_002",
            "evidence_object": "Roth vendor-bearing screen panel",
            "reported_name_or_identifier": "27 compounds in Table 2",
            "supports_c26_dxd": "no",
            "supports_c52_linker_payload": "no",
            "supplier_or_catalog_reported": "for_table2_compounds_only",
            "exact_structure_reported": "no",
            "finding": "Deruxtecan and Exatecan are absent from Table 2; Table 2 cannot establish the source of the Table 3 Deruxtecan reagent.",
            "model_action": "do_not_impute_vendor_from_table2",
            "source_id": "SRC_B6_002",
            "source_location": "Methods and Table 2",
        },
        {
            "audit_id": "B6_IDENT_003",
            "evidence_object": "Roth confirmatory assay",
            "reported_name_or_identifier": "Deruxtecan; Table 3 RR 20 ABCB1 and 240 ABCG2",
            "supports_c26_dxd": "ambiguous",
            "supports_c52_linker_payload": "ambiguous",
            "supplier_or_catalog_reported": "no",
            "exact_structure_reported": "no",
            "finding": "The compound is an additional confirmatory payload, but its vial identity is not documented.",
            "model_action": "retain_continuous_rr_exclude_structure_based_metric_until_resolved",
            "source_id": "SRC_B6_002",
            "source_location": "Figure 2 and Table 3",
        },
        {
            "audit_id": "B6_IDENT_004",
            "evidence_object": "Roth supplementary material",
            "reported_name_or_identifier": "Supplementary Table 1 and figures",
            "supports_c26_dxd": "no",
            "supports_c52_linker_payload": "no",
            "supplier_or_catalog_reported": "no",
            "exact_structure_reported": "no",
            "finding": "The supplement documents cell lines and response figures but does not add compound procurement metadata.",
            "model_action": "identity_remains_open",
            "source_id": "SRC_B6_002",
            "source_location": "Supplementary Materials PDF",
        },
        {
            "audit_id": "B6_IDENT_005",
            "evidence_object": "PubChem exact-name resolution",
            "reported_name_or_identifier": "Deruxtecan -> CID 118305111",
            "supports_c26_dxd": "no",
            "supports_c52_linker_payload": "yes_registry",
            "supplier_or_catalog_reported": "not_applicable",
            "exact_structure_reported": "yes_registry",
            "finding": "Exact-name query returns a C52 maleimide-GGFG-linker-DXd species, MW 1034.1.",
            "model_action": "one_identity_hypothesis_only_not_safe_default_for_roth",
            "source_id": "SRC_B6_003",
            "source_location": "PUG REST property record",
        },
        {
            "audit_id": "B6_IDENT_006",
            "evidence_object": "PubChem payload synonym resolution",
            "reported_name_or_identifier": "DXd -> CID 117888634",
            "supports_c26_dxd": "yes_registry",
            "supports_c52_linker_payload": "no",
            "supplier_or_catalog_reported": "not_applicable",
            "exact_structure_reported": "yes_registry",
            "finding": "DXd query returns the C26 released payload, MW 493.5.",
            "model_action": "use_for_ADC_REL_010_release_species_only",
            "source_id": "SRC_B6_004",
            "source_location": "PUG REST property record",
        },
        {
            "audit_id": "B6_IDENT_007",
            "evidence_object": "Final identity judgment",
            "reported_name_or_identifier": "Roth Deruxtecan",
            "supports_c26_dxd": "plausible",
            "supports_c52_linker_payload": "plausible",
            "supplier_or_catalog_reported": "no",
            "exact_structure_reported": "no",
            "finding": "The paper's language favors a DXd interpretation, while the exact chemical name favors a marketed C52 linker-payload interpretation; neither resolves the vial.",
            "model_action": "mandatory_identity_abstention_and_dual_structure_sensitivity_only",
            "source_id": "SRC_B6_002 | SRC_B6_003 | SRC_B6_004",
            "source_location": "cross-source audit",
        },
    ]


def lys_mcc_primary_table_audit_rows() -> list[dict]:
    return [
        {
            "audit_id": "B6_LMD_001",
            "source_id": "SRC_B6_005",
            "document": "Rock et al. 2015 dmd.115.064253",
            "evidence_location": "main text Discussion -> Supplemental Table 2",
            "what_is_directly_supported": "Lys-MCC-DM1 was obtained from ImmunoGen; low permeability was measured in LLC-PK1 cells.",
            "numeric_table_retrieved": "no",
            "relationship_to_foti_figure8": "same analyte and LLC-PK1 family; exact dataset identity not proven",
            "current_action": "retain_as_primary_low_permeability_source_but_do_not_invent_papp_values",
            "access_or_quality_note": "Current ASPET supplement endpoint returned HTTP 403; the accessible four-page main article does not include the supplement.",
        },
        {
            "audit_id": "B6_LMD_002",
            "source_id": "SRC_B6_006",
            "document": "Foti et al. 2015 dmd.115.064576",
            "evidence_location": "Figure 8 and accompanying methods",
            "what_is_directly_supported": "14C-Lys-MCC-DM1 bidirectional transport in LLC-PK1 and human MDR1-transfected cells at 0.1, 1, 5 and 10 uM; source concludes low passive permeability and no evidence of P-gp substrate behavior.",
            "numeric_table_retrieved": "no_author_table; figure_only",
            "relationship_to_foti_figure8": "direct source",
            "current_action": "retain_algorithmic_bar_digitization_as_exploratory_only",
            "access_or_quality_note": "No machine-readable values or error-bar numbers are available in the accessible article representation.",
        },
        {
            "audit_id": "B6_LMD_003",
            "source_id": "SRC_B6_005 | SRC_B6_006",
            "document": "cross-source provenance judgment",
            "evidence_location": "Rock Supplemental Table 2 citation versus Foti Figure 8",
            "what_is_directly_supported": "Both sources support low permeability; only Foti explicitly reports the MDR1/P-gp inhibitor experiment in the accessible text.",
            "numeric_table_retrieved": "no",
            "relationship_to_foti_figure8": "possibly related experiment; not safe to equate without the supplement",
            "current_action": "do_not_replace_figure_digitization_with_main_article_claim",
            "access_or_quality_note": "A second independent digitization or the original supplement/author data remains necessary before truth-set promotion.",
        },
    ]


def pubchem_recheck_rows() -> list[dict]:
    roles = {
        2061167: "project_summary",
        2202241: "paper_declared_second_aid",
    }
    aids = [2061167, 2061304, 2061305, 2061306, 2061317, 2061318, 2061319, 2061320, 2061321, 2061322, 2061323, 2061324, 2202239, 2202241]
    rows = []
    for aid in aids:
        rows.append(
            {
                "aid": aid,
                "role": roles.get(aid, "related_assay"),
                "pug_rest_csv_http_status": 404,
                "checked_at": RUN_DATE,
                "availability_change_since_2026_07_26": "none",
                "model_or_data_action": "retain_access_log_no_raw_assay_rows",
                "url": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/CSV",
            }
        )
    return rows


def external_truth_v2() -> list[dict]:
    rows = read_csv(PROCESSED / "external_truth_set_v1_pgp_bcrp_continuous_26.csv")
    output = []
    for row in rows:
        ambiguous = row["compound_id"] == "ADC_PAYLOAD_029"
        output.append(
            {
                **row,
                "structure_identity_status": "unresolved_C26_DXd_vs_C52_linker_payload" if ambiguous else "registered_structure_retained",
                "structure_based_external_metric_use": "exclude_until_reagent_identity_resolved" if ambiguous else "include",
                "identity_audit_source": "SRC_B6_002 | SRC_B6_003 | SRC_B6_004" if ambiguous else "",
            }
        )
    return output


def audited_roth_predictions_v2() -> list[dict]:
    rows = read_csv(PROCESSED / "baseline_roth_external_predictions_v1.csv")
    output = []
    for row in rows:
        ambiguous = row["compound_id"] == "ADC_PAYLOAD_029"
        output.append(
            {
                **row,
                "structure_identity_status": "unresolved_C26_DXd_vs_C52_linker_payload" if ambiguous else "registered_structure_retained",
                "prediction_validity": "invalid_for_material_specific_interpretation" if ambiguous else "unchanged",
                "identity_safe_external_metric_include": "no" if ambiguous else "yes",
                "validity_reason": "Roth does not report supplier/catalog/structure for the additional Deruxtecan reagent; the v1 probability used PubChem C52 by exact-name resolution." if ambiguous else "",
            }
        )
    return output


def identity_safe_metrics(audited_predictions: list[dict]) -> list[dict]:
    metrics = []
    for gene in ("ABCB1", "ABCG2"):
        for model_name in ("dummy_prior", "physchem_logistic", "morgan_logistic"):
            rows = [
                row for row in audited_predictions
                if row["transporter_gene"] == gene
                and row["model_name"] == model_name
                and row["identity_safe_external_metric_include"] == "yes"
                and row["predicted_probability"] != ""
            ]
            ratios = np.array([float(row["resistance_ratio"]) for row in rows])
            probabilities = np.array([float(row["predicted_probability"]) for row in rows])
            rho_result = spearmanr(ratios, probabilities)
            rho = float(rho_result.statistic)
            p_value = float(rho_result.pvalue)
            metrics.append(
                {
                    "transporter_gene": gene,
                    "model_name": model_name,
                    "n_compounds": len(rows),
                    "excluded_identity_ambiguous_compounds": "Deruxtecan (Roth Table 3)",
                    "spearman_probability_vs_continuous_rr": round(rho, 6) if math.isfinite(rho) else "",
                    "spearman_p_value": round(p_value, 6) if math.isfinite(p_value) else "",
                    "outside_similarity_support_count": sum(row["abstain_recommendation"] == "yes" for row in rows),
                    "low_similarity_support_count": sum(row["abstain_recommendation"] == "yes_caution" for row in rows),
                    "binary_threshold_used": "no",
                    "interpretation": "identity_safe_directional_association_only_small_external_panel",
                }
            )
    return metrics


def gate_summary_v2(metrics: list[dict]) -> list[dict]:
    rows = read_csv(PROCESSED / "baseline_model_gate_summary_v1.csv")
    lookup = {(row["transporter_gene"], row["model_name"]): row for row in metrics}
    for row in rows:
        metric = lookup[(row["transporter_gene"], row["model_name"])]
        row["roth_external_n"] = metric["n_compounds"]
        row["roth_spearman_probability_vs_rr"] = metric["spearman_probability_vs_continuous_rr"]
        rho = metric["spearman_probability_vs_continuous_rr"]
        row["roth_directional_signal"] = "positive" if rho != "" and float(rho) > 0 else "negative" if rho != "" and float(rho) < 0 else "none_or_undefined"
        row["final_predictor_claim_status"] = "not_eligible_identity_safe_external_panel_small_and_no_independent_matched_release_validation"
    return rows


def deruxtecan_structure_sensitivity() -> tuple[list[dict], list[dict]]:
    free_rows = read_csv(PROCESSED / "compound_registry_standardized_v1.csv")
    release_rows = read_csv(PROCESSED / "adc_release_species_registry_standardized_v1.csv")
    decisions = read_csv(PROCESSED / "daood_training_leakage_decisions_v1.csv")
    c52 = next(row for row in free_rows if row["reference_id"] == "ADC_PAYLOAD_029")
    c26 = next(row for row in release_rows if row["reference_id"] == "ADC_REL_010")
    hypotheses = [
        ("C52_PubChem_Deruxtecan_linker_payload", c52, "SRC_B6_003"),
        ("C26_PubChem_DXd_released_payload", c26, "SRC_B6_004"),
    ]
    rr = {"ABCB1": 20.0, "ABCG2": 240.0}
    rows = []
    for gene in ("ABCB1", "ABCG2"):
        training_rows = [
            row for row in decisions
            if row["transporter_gene"] == gene and row["leakage_model_action"] == "include_leakage_clean_training"
        ]
        for model_name in ("dummy_prior", "physchem_logistic", "morgan_logistic"):
            bundle = joblib.load(MODELS / f"{gene}_{model_name}_baseline.joblib")
            model = bundle["model"]
            thresholds = bundle["similarity_support_thresholds"]
            for hypothesis, ref, source_id in hypotheses:
                probability = baseline.predict_probability(model, model_name, ref)
                maximum, band, abstain = baseline.support_for(ref, training_rows, thresholds)
                rows.append(
                    {
                        "transporter_gene": gene,
                        "model_name": model_name,
                        "roth_resistance_ratio": rr[gene],
                        "structure_hypothesis": hypothesis,
                        "structure_reference_id": ref["reference_id"],
                        "pubchem_cid": ref["source_pubchem_cid"],
                        "molecular_formula": "C52H56FN9O13" if hypothesis.startswith("C52") else "C26H24FN3O6",
                        "predicted_probability": round(float(probability), 8),
                        "max_clean_training_morgan_tanimoto": round(float(maximum), 6),
                        "similarity_support_band": band,
                        "abstain_recommendation": "yes_identity_unresolved" if abstain == "no_similarity_abstention" else f"{abstain}_plus_identity_unresolved",
                        "interpretation": "dual_structure_sensitivity_only_not_material_specific_prediction",
                        "source_id": source_id,
                    }
                )
    lookup = {(row["transporter_gene"], row["model_name"], row["structure_hypothesis"]): row for row in rows}
    deltas = []
    for gene in ("ABCB1", "ABCG2"):
        for model_name in ("dummy_prior", "physchem_logistic", "morgan_logistic"):
            c52_row = lookup[(gene, model_name, "C52_PubChem_Deruxtecan_linker_payload")]
            c26_row = lookup[(gene, model_name, "C26_PubChem_DXd_released_payload")]
            deltas.append(
                {
                    "transporter_gene": gene,
                    "model_name": model_name,
                    "c52_probability": c52_row["predicted_probability"],
                    "c26_probability": c26_row["predicted_probability"],
                    "c26_minus_c52_probability": round(float(c26_row["predicted_probability"]) - float(c52_row["predicted_probability"]), 8),
                    "c52_max_similarity": c52_row["max_clean_training_morgan_tanimoto"],
                    "c26_max_similarity": c26_row["max_clean_training_morgan_tanimoto"],
                    "c26_minus_c52_similarity": round(float(c26_row["max_clean_training_morgan_tanimoto"]) - float(c52_row["max_clean_training_morgan_tanimoto"]), 6),
                    "material_identity_status": "unresolved",
                    "allowed_use": "sensitivity_analysis_only",
                }
            )
    return rows, deltas


def update_issue_log() -> list[dict]:
    path = DATA_ROOT / "exclusion_and_conflict_log.csv"
    rows = read_csv(path)
    by_id = {row["issue_id"]: row for row in rows}
    by_id["ISSUE-003"].update(
        {
            "issue_type": "test_material_source_and_structure_omitted",
            "current_action": "exclude_Roth_Deruxtecan_from_structure_based_external_metrics; retain_RR_and_dual_structure_sensitivity",
            "reason": "Roth Table 3 adds Deruxtecan but the vendor-bearing Table 2 and supplement omit it. PubChem exact-name Deruxtecan is C52 whereas DXd/MAAA-1181a is C26; the tested vial cannot be identified.",
            "next_evidence_needed": "Author confirmation, purchase record, catalog/CAS, exact mass or structure for the Roth Table 3 reagent",
            "status": "partially_resolved",
        }
    )
    by_id["ISSUE-006"].update(
        {
            "reason": "AID 2061167 project record remains readable, but all project/linked/declaration assay CSV endpoints still returned HTTP 404 on 2026-07-28.",
            "next_evidence_needed": "Future PubChem release or author-provided raw files",
            "status": "open",
        }
    )
    by_id["ISSUE-014"].update(
        {
            "reason": "Rock 2015 confirms an LLC-PK1 permeability Supplemental Table 2 exists, but the current supplement endpoint was not retrievable; the accessible Foti source still exposes only Figure 8 for exact ER values.",
            "next_evidence_needed": "Retrieve Rock Supplemental Table 2 and verify whether it is identical to Foti Figure 8, or obtain author data/independent digitization",
            "status": "open",
        }
    )
    by_id["ISSUE-016"].update(
        {
            "next_evidence_needed": "Independent second curator or author-provided numerical table; Rock Supplemental Table 2 is a targeted retrieval route but identity with Figure 8 is not yet proven",
            "status": "partially_resolved",
        }
    )
    additions = [
        {
            "issue_id": "ISSUE-019",
            "entity_or_source": "DXd/MAAA-1181a PMDA Caco-2 transport",
            "issue_type": "multi_transporter_numeric_not_gene_deconvolved",
            "current_action": "retain_four_numeric_conditions_for_external_audit; exclude_10.7_from_gene_specific_training_truth",
            "reason": "P-gp and BCRP are tested in one Caco-2 inhibitor panel and ER remains 4.43 under dual inhibition, so the baseline ER cannot be uniquely assigned to either transporter.",
            "next_evidence_needed": "Matched single-transporter MDCK/LLC-PK1 assays or knockout/transfectant experiments reporting individual Papp values",
            "status": "open",
        },
        {
            "issue_id": "ISSUE-020",
            "entity_or_source": "Baseline Roth Deruxtecan predictions v1",
            "issue_type": "structure_identity_invalidates_material_specific_interpretation",
            "current_action": "supersede_with_identity_audited_v2_and_remove_from_external_metric",
            "reason": "The v1 baseline used the C52 PubChem exact-name structure for ADC_PAYLOAD_029, but Roth did not establish that this was the tested material.",
            "next_evidence_needed": "Resolve Roth material identity, then rerun with the confirmed structure",
            "status": "partially_resolved",
        },
    ]
    existing_ids = set(by_id)
    for addition in additions:
        if addition["issue_id"] not in existing_ids:
            rows.append(addition)
    write_csv(path, rows)
    return rows


def raw_download_log() -> list[dict]:
    specs = [
        (RAW6 / "PMDA_Enhertu_2020_English_review.pdf", "SRC_B6_001", PMDA_URL, "downloaded_validated_signature"),
        (DATA_ROOT / "raw" / "roth_2026" / "PMC12883344_fulltext.xml", "SRC_B6_002", ROTH_URL, "reused_existing_xml"),
        (DATA_ROOT / "raw" / "roth_2026" / "cdr-9-2-SupplementaryMaterials.pdf", "SRC_B6_002", ROTH_URL, "reused_existing_validated_signature"),
        (RAW6 / "PubChem_Deruxtecan_properties_20260728.json", "SRC_B6_003", PUBCHEM_DERUXTECAN_URL, "downloaded_official_json"),
        (RAW6 / "PubChem_DXd_properties_20260728.json", "SRC_B6_004", PUBCHEM_DXD_URL, "downloaded_official_json"),
        (DATA_ROOT / "raw" / "second_batch_release_species" / "Foti_2015_Figure8_LysMCCDM1_Pgp.jpg", "SRC_B6_006", FOTI_URL, "reused_existing_validated_signature"),
    ]
    rows = []
    for path, source_id, url, status in specs:
        if not path.exists():
            continue
        rows.append(
            {
                "filename": str(path.relative_to(DATA_ROOT)),
                "source_id": source_id,
                "source_url": url,
                "download_status": status,
                "signature_valid": "yes",
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "acquired_or_verified_at": RUN_DATE,
            }
        )
    rows.append(
        {
            "filename": "",
            "source_id": "SRC_B6_005",
            "source_url": "https://dmd.aspetjournals.org/content/suppl/2015/06/22/dmd.115.064253.DC1",
            "download_status": "publisher_supplement_endpoint_http_403_no_file_saved",
            "signature_valid": "not_applicable",
            "size_bytes": "",
            "sha256": "",
            "acquired_or_verified_at": RUN_DATE,
        }
    )
    return rows


def main() -> None:
    dxd = dxd_numeric_rows()
    direct = direct_numeric_v2(dxd)
    release = release_transport_v4()
    matched = matched_comparisons_v2()
    pairs = pair_evidence_v3()
    identity = identity_audit_rows()
    lys_audit = lys_mcc_primary_table_audit_rows()
    pubchem = pubchem_recheck_rows()
    external_truth = external_truth_v2()
    predictions = audited_roth_predictions_v2()
    metrics = identity_safe_metrics(predictions)
    gate = gate_summary_v2(metrics)
    sensitivity, sensitivity_delta = deruxtecan_structure_sensitivity()
    issue_rows = update_issue_log()
    downloads = raw_download_log()

    write_csv(PROCESSED / "dxd_maaa1181a_caco2_transport_numeric_v1.csv", dxd)
    write_csv(PROCESSED / "direct_bidirectional_transport_numeric_v2.csv", direct)
    write_csv(PROCESSED / "adc_release_species_transporter_evidence_long_v4.csv", release)
    write_csv(PROCESSED / "matched_context_comparisons_v2.csv", matched)
    write_csv(PROCESSED / "free_payload_release_pair_evidence_v3.csv", pairs)
    write_csv(PROCESSED / "roth_deruxtecan_identity_audit_v1.csv", identity)
    write_csv(PROCESSED / "lys_mcc_dm1_primary_table_audit_v1.csv", lys_audit)
    write_csv(PROCESSED / "pubchem_assay_availability_recheck_20260728.csv", pubchem)
    write_csv(PROCESSED / "external_truth_set_v2_pgp_bcrp_continuous_identity_audited.csv", external_truth)
    write_csv(PROCESSED / "baseline_roth_external_predictions_v2_identity_audited.csv", predictions)
    write_csv(PROCESSED / "baseline_roth_external_metrics_v2_identity_safe.csv", metrics)
    write_csv(PROCESSED / "baseline_model_gate_summary_v2_identity_safe.csv", gate)
    write_csv(PROCESSED / "roth_deruxtecan_structure_sensitivity_v1.csv", sensitivity)
    write_csv(PROCESSED / "roth_deruxtecan_structure_sensitivity_delta_v1.csv", sensitivity_delta)
    write_csv(PROCESSED / "sixth_batch_source_registry.csv", SOURCES)
    write_csv(PROCESSED / "sixth_batch_raw_download_log.csv", downloads)

    grade_counts = Counter(row["matched_grade"] for row in matched)
    metric_lookup = {(row["transporter_gene"], row["model_name"]): row for row in metrics}
    delta_lookup = {(row["transporter_gene"], row["model_name"]): row for row in sensitivity_delta}
    summary = {
        "generated_at": RUN_DATE,
        "sixth_batch_sources": len(SOURCES),
        "raw_source_files_or_attempts": len(downloads),
        "new_dxd_numeric_conditions": len(dxd),
        "direct_numeric_rows_v2": len(direct),
        "release_transporter_records_v4": len(release),
        "matched_comparisons_v2": len(matched),
        "matched_grade_counts": dict(sorted(grade_counts.items())),
        "m0_exact_same_assay_pairs": grade_counts.get("M0", 0),
        "roth_identity_audit_rows": len(identity),
        "pubchem_aids_rechecked": len(pubchem),
        "pubchem_csv_rows_available": 0,
        "roth_identity_safe_external_n_per_gene": 12,
        "identity_safe_spearman": {
            gene: {
                model: metric_lookup[(gene, model)]["spearman_probability_vs_continuous_rr"]
                for model in ("dummy_prior", "physchem_logistic", "morgan_logistic")
            }
            for gene in ("ABCB1", "ABCG2")
        },
        "deruxtecan_structure_probability_deltas_c26_minus_c52": {
            gene: {
                model: delta_lookup[(gene, model)]["c26_minus_c52_probability"]
                for model in ("dummy_prior", "physchem_logistic", "morgan_logistic")
            }
            for gene in ("ABCB1", "ABCG2")
        },
        "open_or_partial_issue_count": sum(row["status"] in {"open", "partially_resolved"} for row in issue_rows),
        "model_binary_labels_created": 0,
        "key_results": [
            "PMDA supplies numeric DXd/MAAA-1181a Caco-2 ER values: 10.7 without inhibitor, 4.37 with verapamil, 4.44 with novobiocin, and 4.43 with GF120918 at 1 uM.",
            "The Caco-2 numeric panel supports P-gp/BCRP substrate direction but cannot be converted into separate ABCB1 and ABCG2 ER truth values; residual dual-inhibitor ER is 4.43.",
            "Roth Deruxtecan reagent identity remains unresolved because the additional Table 3 compound is absent from the vendor-bearing Table 2 and no catalog, CAS or structure appears in the supplement.",
            "The v1 structure-based Roth predictions for Deruxtecan are invalid for material-specific interpretation and are excluded from identity-safe external metrics.",
            "Rock 2015 points to an inaccessible Supplemental Table 2 for LLC-PK1 permeability, so Foti Figure 8 digitization remains exploratory rather than upgraded truth.",
            "PubChem project and related assay CSV endpoints remained unavailable on 2026-07-28.",
        ],
        "feasibility_conclusion": "The ADC-specific evidence layer is stronger, particularly for DXd, but the new data reinforce the need for target-specific deconvolution and material-identity abstention. No M0 pair was found and the universal predictor claim gate remains closed.",
    }
    (PROCESSED / "sixth_batch_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
