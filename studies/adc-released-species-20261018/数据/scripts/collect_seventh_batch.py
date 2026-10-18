#!/usr/bin/env python3
"""Build seventh-batch direct transport and same-system ADC evidence tables.

This batch prioritizes author-tabulated bidirectional Papp/efflux data for
actual ADC payloads or named released metabolites. It preserves Caco-2 target
attribution limits, keeps complete-ADC resistance separate from payload
transport truth, and does not create new model binary labels.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


DATA_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = DATA_ROOT.parent
PROCESSED = DATA_ROOT / "processed"
RAW7 = DATA_ROOT / "raw" / "seventh_batch_direct_transport"
RUN_DATE = "2026-07-28"

OBI_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11791482/"
ZD_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11063767/"
PMDA_PADCEV_URL = "https://www.pmda.go.jp/drugs/2021/P20211011002/800126000_30300AMX00454000_I100_1.pdf"
PMDA_ADCETRIS_URL = "https://www.pmda.go.jp/files/000229787.pdf"
MMAE_MDCK_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6650582/"
BAY716_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC7709048/"
KSP_LINKER_URL = "https://doi.org/10.1021/acs.bioconjchem.0c00357"
SME_DM1_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC4324968/"


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
        "source_id": "SRC_B7_001",
        "title": "OBI-992, a Novel TROP2-Targeted Antibody-Drug Conjugate, Demonstrates Antitumor Activity in Multiple Cancer Models",
        "year": 2025,
        "source_type": "primary_article_and_supplement_numeric",
        "url": OBI_URL,
        "use_in_batch": "Caco-2 Papp/ER for SN-38, DXd and exatecan; matched HCC827 payload/ADC resistance and inhibitor rescue",
        "source_location": "Methods; Supplementary Tables S5-S7",
        "primary_or_secondary": "primary_article",
        "notes": "All Caco-2 compounds were tested at 10 uM for 2 h. The paper maps SN-38, MMAE, DXd and exatecan to Acros Organics, BroadPharm, MedChemExpress and Seedchem, respectively; catalog and lot are not reported. Endogenous Caco-2 transport is inhibitor-attributed, not gene-isolated.",
    },
    {
        "source_id": "SRC_B7_002",
        "title": "Design and Evaluation of ZD06519, a Novel Camptothecin Payload for Antibody Drug Conjugates",
        "year": 2024,
        "source_type": "primary_article_supplement_qualitative",
        "url": ZD_URL,
        "use_in_batch": "source-assigned Caco-2 MDR1 and BCRP calls for DXd, exatecan, SN-38, FD1 and FD5",
        "source_location": "Supplementary Table S3",
        "primary_or_secondary": "primary_article",
        "notes": "The table supplies check/X/marginal calls but no Papp or ER values; retain as qualitative corroboration only.",
    },
    {
        "source_id": "SRC_B7_003",
        "title": "Padcev PMDA CTD pharmacokinetic study tables",
        "year": 2021,
        "source_type": "primary_regulatory_study_tables_numeric",
        "url": PMDA_PADCEV_URL,
        "use_in_batch": "author-tabulated Caco-2 Papp and ER for MMAE in Study 8004 and RPT-01709",
        "source_location": "CTD sections 2.6.5.15.1-2; PDF pages 76-77; printed pages 39-40",
        "primary_or_secondary": "primary_regulatory_table",
        "notes": "Study 8004 reports mean(SD); RPT-01709 reports mean +/- SEM, n=3. The latter uses a source-defined ER>2 substrate criterion.",
    },
    {
        "source_id": "SRC_B7_004",
        "title": "Lysosomal P-gp-MDR1 Confers Drug Resistance of Brentuximab Vedotin and Its Cytotoxic Payload Monomethyl Auristatin E in Tumor Cells",
        "year": 2019,
        "source_type": "primary_article_direct_transport",
        "url": MMAE_MDCK_URL,
        "use_in_batch": "independent MDCK-WT/MDCK-MDR1 direct MMAE efflux ratios and inhibitor rescue",
        "source_location": "Methods; Results P-gp-MDR1 Mediated MMAE Efflux; Figure 2A",
        "primary_or_secondary": "primary_article",
        "notes": "At 1 uM MMAE, exact text ERs are 13.6 in MDCK-WT and 44.5 in MDCK-MDR1; quinidine or ketoconazole 10 uM reduces ER below 1.5. Exact Papp values are figure-only and are not digitized here.",
    },
    {
        "source_id": "SRC_B7_005",
        "title": "IL3RA-Targeting Antibody-Drug Conjugate BAY-943 with a Kinesin Spindle Protein Inhibitor Payload",
        "year": 2020,
        "source_type": "primary_article_and_supplement_direct_transport",
        "url": BAY716_URL,
        "use_in_batch": "Caco-2 Papp and ER for the author-named active payload metabolite BAY-716",
        "source_location": "main Results and Supplementary Methods Caco-2 permeability assay",
        "primary_or_secondary": "primary_article",
        "notes": "BAY-716 was tested at 2 uM for 2 h in Caco-2 ACC 169 monolayers; exact Papp A-B 1.8 nm/s, Papp B-A 2.7 nm/s and ER 1.6. Replicate count is not reported.",
    },
    {
        "source_id": "SRC_B7_006",
        "title": "Tailored Linker Chemistries for the Efficient and Selective Activation of ADCs with KSPi Payloads",
        "year": 2020,
        "source_type": "primary_chemistry_article",
        "url": KSP_LINKER_URL,
        "use_in_batch": "chemical-structure context for Bayer legumain-cleavable KSPi active metabolite 8",
        "source_location": "Figure 1; Scheme 1",
        "primary_or_secondary": "primary_article",
        "notes": "The active metabolite structure is shown, but the paper never labels it BAY-716. Do not assert an exact crosswalk without explicit source confirmation.",
    },
    {
        "source_id": "SRC_B7_007",
        "title": "Adcetris PMDA English review report",
        "year": 2019,
        "source_type": "regulatory_review_summary",
        "url": PMDA_ADCETRIS_URL,
        "use_in_batch": "independent regulatory summary supporting MMAE as a P-gp substrate",
        "source_location": "nonclinical pharmacokinetics transporter section",
        "primary_or_secondary": "primary_regulatory_summary",
        "notes": "Use as direction-level regulatory corroboration; do not duplicate a study as an independent experiment when study identifiers overlap.",
    },
    {
        "source_id": "SRC_B7_008",
        "title": "High-Affinity Accumulation of a Maytansinoid in Cells via Weak Tubulin Interaction",
        "year": 2015,
        "source_type": "primary_article_cellular_efflux_mechanism",
        "url": SME_DM1_URL,
        "use_in_batch": "quantitative S-methyl-DM1 cellular retention and rebinding mechanism",
        "source_location": "Methods; Results; Figure 1E",
        "primary_or_secondary": "primary_article",
        "notes": "After loading MCF7 cells with 3.3 nM radiolabeled S-methyl-DM1, 65 +/- 3% remained after 2 h; excess unlabeled compound increased efflux three-fold. This is not a transporter-specific assay.",
    },
]


def obi_caco2_rows() -> list[dict]:
    values = [
        ("estrone-3-sulfate", "assay_control", "", "BCRP", "", 0.19, 0.01, 4.68, 0.89, 24.6),
        ("estrone-3-sulfate", "assay_control", "", "BCRP", "Ko143", 1.54, 0.69, 2.59, 0.14, 1.7),
        ("digoxin", "assay_control", "", "P-gp", "", 0.08, 0.06, 0.42, 0.32, 5.6),
        ("digoxin", "assay_control", "", "P-gp", "verapamil", 0.31, 0.15, 0.87, 0.25, 2.8),
        ("SN-38", "released_payload", "ADC_REL_009", "BCRP", "", 0.26, 0.07, 2.56, 0.68, 9.8),
        ("SN-38", "released_payload", "ADC_REL_009", "BCRP", "Ko143", 2.32, 0.28, 4.04, 0.92, 1.7),
        ("SN-38", "released_payload", "ADC_REL_009", "P-gp", "", 0.32, 0.04, 6.45, 0.36, 20.0),
        ("SN-38", "released_payload", "ADC_REL_009", "P-gp", "verapamil", 0.57, 0.15, 5.90, 0.28, 10.5),
        ("DXd", "released_payload", "ADC_REL_010", "BCRP", "", 0.25, 0.08, 4.42, 0.87, 17.7),
        ("DXd", "released_payload", "ADC_REL_010", "BCRP", "Ko143", 0.34, 0.18, 3.11, 0.40, 9.1),
        ("DXd", "released_payload", "ADC_REL_010", "P-gp", "", 0.12, 0.10, 2.61, 0.16, 22.1),
        ("DXd", "released_payload", "ADC_REL_010", "P-gp", "verapamil", 0.22, 0.07, 2.96, 0.52, 13.8),
        ("Exatecan", "released_payload", "ADC_REL_011", "BCRP", "", 2.60, 0.13, 6.30, 0.48, 2.4),
        ("Exatecan", "released_payload", "ADC_REL_011", "BCRP", "Ko143", 2.74, 1.42, 3.28, 0.71, 1.2),
        ("Exatecan", "released_payload", "ADC_REL_011", "P-gp", "", 1.32, 0.35, 3.59, 0.75, 2.7),
        ("Exatecan", "released_payload", "ADC_REL_011", "P-gp", "verapamil", 2.43, 0.79, 4.22, 0.13, 1.7),
    ]
    baseline = {(name, arm): er for name, _, _, arm, inhibitor, *_, er in values if inhibitor == ""}
    rows = []
    for idx, (name, entity_type, entity_id, arm, inhibitor, ab, ab_sd, ba, ba_sd, er) in enumerate(values, 1):
        base = baseline[(name, arm)]
        gene = "ABCG2" if arm == "BCRP" else "ABCB1"
        if entity_type == "assay_control":
            call = "reference_substrate_control" if not inhibitor else "reference_inhibitor_response"
        elif name in {"SN-38", "DXd"}:
            call = f"{arm}_involved_in_efflux"
        else:
            call = f"lower_{arm}_efflux_than_SN-38_and_DXd"
        rows.append(
            {
                "record_id": f"B7_OBI_S5_{idx:02d}",
                "compound": name,
                "entity_type": entity_type,
                "entity_id": entity_id,
                "transporter_gene": gene,
                "transporter_arm": arm,
                "cell_system": "Caco-2 C2BBe1 (BCRC 60182)",
                "test_concentration_uM": 10,
                "incubation_hours": 2,
                "inhibitor": inhibitor or "none",
                "inhibitor_concentration_uM": "",
                "papp_ab_mean_1e6_cm_s": ab,
                "papp_ab_sd": ab_sd,
                "papp_ba_mean_1e6_cm_s": ba,
                "papp_ba_sd": ba_sd,
                "efflux_ratio": er,
                "baseline_to_condition_er_reduction_fold": round(base / er, 4),
                "uncertainty_type": "SD",
                "source_assigned_call": call,
                "truth_set_use": "continuous_target_inhibitor_attribution_audit_not_gene_isolated_training_truth",
                "source_id": "SRC_B7_001",
                "source_location": "Supplementary Table S5",
                "notes": "Inhibitor concentration is not reported in Table S5 or the accessible permeability methods.",
            }
        )
    return rows


def pmda_mmae_rows() -> list[dict]:
    rows = []
    study_8004 = [
        ("MMAE", 1, "none", "", 0.5, 0.0, 8.1, 1.1, 15.0, 0.8, "ABCB1", "test_article"),
        ("MMAE", 10, "none", "", 0.2, 0.0, 6.8, 1.6, 34.0, 1.1, "ABCB1", "test_article"),
        ("MMAE", 100, "none", "", 0.1, 0.0, 6.7, 0.4, 49.0, 0.3, "ABCB1", "test_article"),
        ("MMAE", 10, "none", "", 0.3, 0.0, 8.4, 2.1, 28.0, 0.3, "ABCB1", "test_article_separate_inhibitor_run"),
        ("MMAE", 10, "PSC833", 10, 1.7, 0.2, 2.3, 0.6, 1.3, 0.3, "ABCB1", "test_article_inhibitor"),
        ("MMAE", 10, "verapamil", 60, 1.9, 0.4, 5.5, 0.3, 2.9, 0.2, "ABCB1", "test_article_inhibitor"),
        ("digoxin", 10, "none", "", 1.3, 0.1, 19.0, 1.4, 14.0, 0.1, "ABCB1", "assay_control"),
        ("digoxin", 10, "verapamil", 60, 4.2, 1.2, 6.3, 0.4, 1.5, 0.3, "ABCB1", "assay_control"),
        ("digoxin", 10, "MMAE", 50, 1.5, 0.1, 17.0, 1.4, 11.0, 0.1, "ABCB1", "inhibition_control"),
    ]
    for idx, values in enumerate(study_8004, 1):
        compound, conc, inhibitor, inh_conc, ab, ab_u, ba, ba_u, er, er_u, gene, role = values
        rows.append(
            {
                "record_id": f"B7_PMDA_8004_{idx:02d}",
                "study_id": "8004",
                "compound": compound,
                "entity_type": "released_payload" if compound == "MMAE" else "assay_control",
                "entity_id": "ADC_REL_007" if compound == "MMAE" else "",
                "transporter_gene": gene,
                "cell_system": "Caco-2 monolayer",
                "test_concentration_uM": conc,
                "inhibitor": inhibitor,
                "inhibitor_concentration_uM": inh_conc,
                "papp_ab_mean_1e6_cm_s": ab,
                "papp_ab_uncertainty": ab_u,
                "papp_ba_mean_1e6_cm_s": ba,
                "papp_ba_uncertainty": ba_u,
                "efflux_ratio": er,
                "efflux_ratio_uncertainty": er_u,
                "uncertainty_type": "SD",
                "replicates": "",
                "condition_role": role,
                "source_assigned_call": "P-gp_substrate" if compound == "MMAE" else "control",
                "truth_set_use": "continuous_direct_ABCB1_evidence_no_new_model_binary_label" if compound == "MMAE" else "assay_control",
                "source_id": "SRC_B7_003",
                "source_location": "CTD 2.6.5.15.1; Study 8004; PDF page 76",
                "notes": "The 10 uM no-inhibitor row belongs to a separate inhibitor experiment and is not pooled with the concentration-series 10 uM row.",
            }
        )
    rpt = [
        ("MMAE", 10, "none", "", 0.10, 0.02, 0.91, 0.04, 8.8, "test_article"),
        ("MMAE", 10, "LY335979", 10, 0.13, 0.03, 0.20, 0.01, 1.6, "P-gp_inhibitor"),
        ("MMAE", 10, "Ko143", 5, 0.15, 0.05, 0.93, 0.06, 6.0, "BCRP_inhibitor"),
        ("MMAE", 10, "indomethacin", 100, 0.25, 0.10, 1.10, 0.02, 4.3, "MRP2_inhibitor"),
        ("propranolol", 10, "none", "", 17.60, 2.90, 27.02, 0.27, 1.5, "transcellular_control"),
        ("nadolol", 10, "none", "", 0.37, 0.02, 0.90, 0.07, 2.4, "paracellular_control"),
    ]
    for idx, values in enumerate(rpt, 1):
        compound, conc, inhibitor, inh_conc, ab, ab_u, ba, ba_u, er, role = values
        rows.append(
            {
                "record_id": f"B7_PMDA_RPT01709_{idx:02d}",
                "study_id": "RPT-01709",
                "compound": compound,
                "entity_type": "released_payload" if compound == "MMAE" else "assay_control",
                "entity_id": "ADC_REL_007" if compound == "MMAE" else "",
                "transporter_gene": "ABCB1|ABCG2|ABCC2" if compound == "MMAE" else "",
                "cell_system": "Caco-2 monolayer",
                "test_concentration_uM": conc,
                "inhibitor": inhibitor,
                "inhibitor_concentration_uM": inh_conc,
                "papp_ab_mean_1e6_cm_s": ab,
                "papp_ab_uncertainty": ab_u,
                "papp_ba_mean_1e6_cm_s": ba,
                "papp_ba_uncertainty": ba_u,
                "efflux_ratio": er,
                "efflux_ratio_uncertainty": "",
                "uncertainty_type": "SEM",
                "replicates": 3,
                "condition_role": role,
                "source_assigned_call": "source_rule_ER_gt_2_indicates_efflux_pump_substrate" if compound == "MMAE" else "control",
                "truth_set_use": "continuous_multi_transporter_inhibitor_panel_not_gene_specific_training_truth" if compound == "MMAE" else "assay_control",
                "source_id": "SRC_B7_003",
                "source_location": "CTD 2.6.5.15.2; Study RPT-01709; PDF page 77",
                "notes": "The source states an ER>2 criterion, but the project does not convert it into a new model binary label; inhibitor selectivity is interpreted by arm only.",
            }
        )
    return rows


def mmae_mdck_rows() -> list[dict]:
    return [
        {
            "record_id": "B7_MMAE_MDCK_01",
            "compound": "MMAE",
            "entity_id": "ADC_REL_007",
            "cell_system": "MDCK-WT",
            "transporter_gene": "ABCB1",
            "test_concentration_uM": 1,
            "efflux_ratio": 13.6,
            "background_corrected_ratio_of_ratios": "",
            "inhibitor_condition": "quinidine 10 uM or ketoconazole 10 uM reduces ER below 1.5",
            "source_assigned_call": "parental_background_directional_efflux",
            "numeric_origin": "exact_author_text_ER; Papp bars not digitized",
            "truth_set_use": "independent_direct_continuous_evidence",
            "source_id": "SRC_B7_004",
            "source_location": "Results P-gp-MDR1 Mediated MMAE Efflux; Figure 2A",
        },
        {
            "record_id": "B7_MMAE_MDCK_02",
            "compound": "MMAE",
            "entity_id": "ADC_REL_007",
            "cell_system": "MDCK-MDR1",
            "transporter_gene": "ABCB1",
            "test_concentration_uM": 1,
            "efflux_ratio": 44.5,
            "background_corrected_ratio_of_ratios": 3.3,
            "inhibitor_condition": "quinidine 10 uM or ketoconazole 10 uM reduces ER below 1.5",
            "source_assigned_call": "P-gp-MDR1_substrate",
            "numeric_origin": "exact_author_text_ER; Papp bars not digitized",
            "truth_set_use": "independent_direct_continuous_ABCB1_evidence",
            "source_id": "SRC_B7_004",
            "source_location": "Results P-gp-MDR1 Mediated MMAE Efflux; Figure 2A",
        },
    ]


def bay716_rows() -> list[dict]:
    return [
        {
            "record_id": "B7_BAY716_01",
            "compound": "BAY-716",
            "entity_id": "ADC_REL_012",
            "entity_type": "author_named_active_payload_metabolite",
            "cell_system": "Caco-2 ACC 169",
            "transporter_gene": "ABCB1",
            "test_concentration_uM": 2,
            "incubation_hours": 2,
            "papp_ab_mean_nm_s": 1.8,
            "papp_ba_mean_nm_s": 2.7,
            "papp_ab_mean_1e6_cm_s": 0.18,
            "papp_ba_mean_1e6_cm_s": 0.27,
            "efflux_ratio": 1.6,
            "replicates": "not_reported",
            "source_assigned_call": "poor_P-gp_substrate_no_active_efflux_in_Caco-2",
            "truth_set_use": "continuous_source_interpretation_not_gene_isolated_training_truth",
            "source_id": "SRC_B7_005",
            "source_location": "main Results and Supplementary Methods Caco-2 permeability assay",
            "notes": "Converted 1 nm/s to 0.1 x10^-6 cm/s for harmonized columns. Exact chemical structure/crosswalk remains unresolved.",
        }
    ]


def zd_calls() -> list[dict]:
    calls = {
        "FD1": {"ABCB1": "substrate", "ABCG2": "marginal"},
        "FD5": {"ABCB1": "not_substrate", "ABCG2": "not_substrate"},
        "DXd": {"ABCB1": "substrate", "ABCG2": "marginal"},
        "Exatecan": {"ABCB1": "not_substrate", "ABCG2": "not_substrate"},
        "SN-38": {"ABCB1": "substrate", "ABCG2": "substrate"},
    }
    ids = {"DXd": "ADC_REL_010", "Exatecan": "ADC_REL_011", "SN-38": "ADC_REL_009"}
    rows = []
    for compound, per_gene in calls.items():
        for gene, call in per_gene.items():
            rows.append(
                {
                    "record_id": f"B7_ZD_{compound.upper().replace('-', '')}_{gene}",
                    "compound": compound,
                    "release_species_id": ids.get(compound, ""),
                    "transporter_gene": gene,
                    "assay_system": "Caco-2 source-assigned MDR1/BCRP call",
                    "reported_symbol": "marginal" if call == "marginal" else "check" if call == "substrate" else "X",
                    "normalized_source_call": call,
                    "numeric_papp_or_er_available": "no",
                    "model_binary_label": "",
                    "truth_set_use": "qualitative_corroboration_only",
                    "source_id": "SRC_B7_002",
                    "source_location": "Supplementary Table S3",
                    "notes": "Do not infer a project threshold or numerical ER from the check/X/marginal symbol.",
                }
            )
    return rows


def obi_hcc827_rr_rows() -> list[dict]:
    values = [
        ("MMAE", "payload_control", "ADC_REL_007", 0.34, 596.1, 0.13, 1753.2, 0.4),
        ("SN-38", "payload", "ADC_REL_009", 4.75, 5.1, 48.6, 1.1, 10.2),
        ("DXd", "payload", "ADC_REL_010", 4.27, 14.8, 39.6, 3.5, 9.3),
        ("Exatecan", "payload", "ADC_REL_011", 0.49, 0.26, 2.1, 0.53, 4.3),
        ("Sacituzumab govitecan", "complete_ADC", "", 0.58, 1.05, 19.5, 1.8, 33.6),
        ("Datopotamab deruxtecan", "complete_ADC", "", 12.3, 370.3, 582.7, 30.1, 47.4),
        ("OBI-992", "complete_ADC", "", 77.1, 122.2, 238.6, 1.6, 3.1),
    ]
    rows = []
    for idx, values_row in enumerate(values, 1):
        name, entity_type, release_id, parent, abcb1, abcg2, rr1, rr2 = values_row
        rows.append(
            {
                "record_id": f"B7_OBI_S6_{idx:02d}",
                "entity_name": name,
                "entity_type": entity_type,
                "release_species_id": release_id,
                "cell_background": "HCC827 parental and stable overexpression derivatives",
                "parental_ic50_nM": parent,
                "abcb1_ic50_nM": abcb1,
                "abcg2_ic50_nM": abcg2,
                "abcb1_resistance_ratio": rr1,
                "abcg2_resistance_ratio": rr2,
                "endpoint": "6-day cell viability IC50",
                "training_use": "payload_functional_context_only" if "payload" in entity_type else "complete_ADC_context_not_payload_transport_truth",
                "source_id": "SRC_B7_001",
                "source_location": "Supplementary Table S6",
                "notes": "RR is IC50(overexpression)/IC50(parental); complete ADC rows reflect the whole delivery system.",
            }
        )
    return rows


def obi_inhibitor_rescue_rows() -> list[dict]:
    values = [
        ("Sacituzumab govitecan", "ABCB1", 0.58, 1.05, 1.14, 0.92, "verapamil", 40),
        ("Sacituzumab govitecan", "ABCG2", 0.58, 19.47, 0.44, 44.25, "Ko143", 10),
        ("Datopotamab deruxtecan", "ABCB1", 12.3, 370.3, 45.39, 8.16, "verapamil", 40),
        ("Datopotamab deruxtecan", "ABCG2", 12.3, 582.7, 0.50, 1165.4, "Ko143", 10),
        ("OBI-992", "ABCB1", 77.1, 122.2, 110.7, 1.10, "verapamil", 40),
        ("OBI-992", "ABCG2", 77.1, 238.6, 84.55, 2.82, "Ko143", 10),
    ]
    rows = []
    for idx, values_row in enumerate(values, 1):
        adc, gene, parent, no_inh, with_inh, fold, inhibitor, inhibitor_conc = values_row
        rows.append(
            {
                "record_id": f"B7_OBI_S7_{idx:02d}",
                "adc_name": adc,
                "transporter_gene": gene,
                "parental_ic50_nM": parent,
                "overexpression_ic50_no_inhibitor_nM": no_inh,
                "overexpression_ic50_with_inhibitor_nM": with_inh,
                "restored_sensitivity_fold": fold,
                "inhibitor": inhibitor,
                "inhibitor_concentration_uM": inhibitor_conc,
                "endpoint": "6-day cell viability IC50 inhibitor rescue",
                "training_use": "complete_ADC_functional_context_not_payload_transport_truth",
                "source_id": "SRC_B7_001",
                "source_location": "Supplementary Table S7",
                "notes": "Fold is IC50(no inhibitor)/IC50(with inhibitor) as defined by the source.",
            }
        )
    return rows


def same_system_payload_adc_rows(rr_rows: list[dict]) -> list[dict]:
    by_name = {row["entity_name"]: row for row in rr_rows}
    mappings = [
        ("SN-38", "Sacituzumab govitecan", "SG"),
        ("DXd", "Datopotamab deruxtecan", "Dato-DXd"),
        ("Exatecan", "OBI-992", "OBI-992"),
    ]
    rows = []
    for payload, adc, platform in mappings:
        for gene, field in (("ABCB1", "abcb1_resistance_ratio"), ("ABCG2", "abcg2_resistance_ratio")):
            payload_rr = float(by_name[payload][field])
            adc_rr = float(by_name[adc][field])
            rows.append(
                {
                    "comparison_id": f"B7_CTX_{platform.upper().replace('-', '')}_{gene}",
                    "platform": platform,
                    "transporter_gene": gene,
                    "payload_name": payload,
                    "payload_release_species_id": by_name[payload]["release_species_id"],
                    "payload_resistance_ratio": payload_rr,
                    "complete_adc_name": adc,
                    "complete_adc_resistance_ratio": adc_rr,
                    "same_cell_background": "yes_HCC827_parental_and_same_overexpression_derivatives",
                    "same_endpoint": "yes_6-day_IC50_RR",
                    "adc_to_payload_rr_ratio": round(adc_rr / payload_rr, 4) if payload_rr else "",
                    "allowed_interpretation": "same_system_context_link_only",
                    "training_use": "no_complete_ADC_is_not_payload_transport_truth",
                    "source_id": "SRC_B7_001",
                    "source_location": "Supplementary Table S6",
                    "limitation": "ADC uptake, trafficking, release, antigen expression and payload efflux all contribute; do not interpret the ratio as transporter kinetics.",
                }
            )
    return rows


def membrane_retention_rows() -> list[dict]:
    return [
        {
            "record_id": "B7_MEM_001",
            "entity_id": "ADC_REL_005",
            "entity_name": "S-methyl-DM1",
            "cell_system": "MCF7 cells",
            "loading_condition": "3.3 nM radiolabeled S-methyl-DM1 for 2 h at 37 C",
            "washout_condition": "fresh medium for 2 h at 37 C",
            "retained_fraction_mean_percent": 65,
            "retained_fraction_uncertainty_percent": 3,
            "uncertainty_type": "source reports +/-; type not explicitly stated in result sentence",
            "replicates": 2,
            "perturbation": "0.33 uM unlabeled S-methyl-DM1 during washout",
            "perturbation_effect": "efflux_rate_increased_3_fold",
            "mechanistic_interpretation": "intracellular tubulin rebinding contributes to retention after passive entry",
            "transporter_specific": "no",
            "training_use": "membrane_distribution_and_retention_mechanistic_layer_only",
            "source_id": "SRC_B7_008",
            "source_location": "Methods; Results; Figure 1E",
        }
    ]


def release_registry_v3() -> list[dict]:
    rows = read_csv(PROCESSED / "adc_release_species_registry_v2.csv")
    compounds = {row["compound_id"]: row for row in read_csv(PROCESSED / "compound_registry_v1.csv")}
    exatecan = compounds["ADC_PAYLOAD_028"]
    rows.extend(
        [
            {
                "release_species_id": "ADC_REL_011",
                "release_species_name": "Exatecan",
                "aliases": "free exatecan payload",
                "payload_core": "Exatecan",
                "residual_group": "none",
                "species_class": "free_payload",
                "representative_adc_or_format": "OBI-992",
                "linker_class": "maleimide-PEG24-VA-PAB cleavable linker",
                "release_identity_status": "author_assigned_payload_not_direct_catabolite_MS",
                "release_evidence_method": "ADC design and matched payload activity",
                "source_ids": "SRC_B7_001",
                "structure_query_names": "Exatecan",
                "existing_payload_compound_id": "ADC_PAYLOAD_028",
                "modeling_status": "already_in_free_payload_registry",
                "pubchem_query_name": exatecan["pubchem_query_name"],
                "pubchem_cid": exatecan["pubchem_cid"],
                "pubchem_title": exatecan["pubchem_title"],
                "isomeric_smiles": exatecan["isomeric_smiles"],
                "connectivity_smiles": exatecan["connectivity_smiles"],
                "inchikey": exatecan["inchikey"],
                "molecular_formula": exatecan["molecular_formula"],
                "molecular_weight": exatecan["molecular_weight"],
                "xlogp": exatecan["xlogp"],
                "tpsa": exatecan["tpsa"],
                "formal_charge": exatecan["formal_charge"],
                "structure_status": "resolved",
                "structure_note": "Exact free exatecan structure is resolved; intracellular catabolite identity for OBI-992 is author-assigned rather than directly measured in the cited paper.",
                "structure_source_url": exatecan["structure_source_url"],
            },
            {
                "release_species_id": "ADC_REL_012",
                "release_species_name": "BAY-716",
                "aliases": "active toxophore metabolite of BAY-943",
                "payload_core": "proprietary pyrrole KSP inhibitor",
                "residual_group": "physicochemical modifier retained after legumain cleavage",
                "species_class": "author_named_active_payload_metabolite",
                "representative_adc_or_format": "BAY-943 IL3RA-ADC",
                "linker_class": "lysine-linked legumain-cleavable peptide",
                "release_identity_status": "author_named_active_metabolite_structure_crosswalk_unresolved",
                "release_evidence_method": "ADC design, biological activity and author naming",
                "source_ids": "SRC_B7_005 | SRC_B7_006",
                "structure_query_names": "BAY-716",
                "existing_payload_compound_id": "",
                "modeling_status": "exclude_until_exact_structure_crosswalk_confirmed",
                "pubchem_query_name": "",
                "pubchem_cid": "",
                "pubchem_title": "",
                "isomeric_smiles": "",
                "connectivity_smiles": "",
                "inchikey": "",
                "molecular_formula": "",
                "molecular_weight": "",
                "xlogp": "",
                "tpsa": "",
                "formal_charge": "",
                "structure_status": "primary_structure_context_but_exact_name_crosswalk_unresolved",
                "structure_note": "Lerchen 2020 shows active KSPi metabolite 8, but neither source explicitly states that compound 8 is BAY-716; do not model by inferred structure.",
                "structure_source_url": f"{BAY716_URL} | {KSP_LINKER_URL}",
            },
        ]
    )
    return rows


def direct_numeric_v3(obi: list[dict], pmda: list[dict], mdck: list[dict], bay: list[dict]) -> list[dict]:
    rows = read_csv(PROCESSED / "direct_bidirectional_transport_numeric_v2.csv")
    for source in obi:
        rows.append(
            {
                "transport_record_id": source["record_id"],
                "entity_type": "released_species" if source["entity_id"] else "assay_control",
                "entity_id": source["entity_id"],
                "entity_name": source["compound"],
                "payload_core": source["compound"] if source["entity_id"] else "",
                "transporter_gene": source["transporter_gene"],
                "cell_background": source["cell_system"],
                "transporter_system": "endogenous multi-transporter Caco-2 with target-selective inhibitor arm",
                "assay_endpoint": "bidirectional apparent permeability; ER=Papp B-A/Papp A-B",
                "concentration_value": source["test_concentration_uM"],
                "concentration_unit": "uM",
                "inhibitor": source["inhibitor"],
                "inhibitor_concentration_value": source["inhibitor_concentration_uM"],
                "inhibitor_concentration_unit": "",
                "papp_ab_mean": source["papp_ab_mean_1e6_cm_s"],
                "papp_ab_sd": source["papp_ab_sd"],
                "papp_ba_mean": source["papp_ba_mean_1e6_cm_s"],
                "papp_ba_sd": source["papp_ba_sd"],
                "efflux_ratio_mean": source["efflux_ratio"],
                "efflux_ratio_sd": "",
                "background_corrected_efflux_ratio": "",
                "source_assigned_call": source["source_assigned_call"],
                "evidence_tier": "A_direct_numeric_target_inhibitor_attributed_Caco2",
                "numeric_origin": "source_supplement_table",
                "digitization_status": "",
                "truth_set_eligible": "no_gene_isolated_training; yes_continuous_external_audit",
                "comparability_group_id": f"OBI_CACO2_{source['compound'].upper().replace('-', '')}_{source['transporter_gene']}_10uM",
                "source_id": source["source_id"],
                "source_url": OBI_URL,
                "source_location": source["source_location"],
                "notes": source["notes"],
            }
        )
    for source in pmda:
        rows.append(
            {
                "transport_record_id": source["record_id"],
                "entity_type": "released_species" if source["entity_id"] else "assay_control",
                "entity_id": source["entity_id"],
                "entity_name": source["compound"],
                "payload_core": "MMAE" if source["entity_id"] else "",
                "transporter_gene": source["transporter_gene"],
                "cell_background": source["cell_system"],
                "transporter_system": "endogenous Caco-2 inhibitor panel",
                "assay_endpoint": "bidirectional apparent permeability; ER=Papp B-A/Papp A-B",
                "concentration_value": source["test_concentration_uM"],
                "concentration_unit": "uM",
                "inhibitor": source["inhibitor"],
                "inhibitor_concentration_value": source["inhibitor_concentration_uM"],
                "inhibitor_concentration_unit": "uM" if source["inhibitor_concentration_uM"] != "" else "",
                "papp_ab_mean": source["papp_ab_mean_1e6_cm_s"],
                "papp_ab_sd": source["papp_ab_uncertainty"] if source["uncertainty_type"] == "SD" else "",
                "papp_ba_mean": source["papp_ba_mean_1e6_cm_s"],
                "papp_ba_sd": source["papp_ba_uncertainty"] if source["uncertainty_type"] == "SD" else "",
                "efflux_ratio_mean": source["efflux_ratio"],
                "efflux_ratio_sd": source["efflux_ratio_uncertainty"] if source["uncertainty_type"] == "SD" else "",
                "background_corrected_efflux_ratio": "",
                "source_assigned_call": source["source_assigned_call"],
                "evidence_tier": "A_direct_numeric_regulatory_table",
                "numeric_origin": f"source_table_{source['uncertainty_type']}",
                "digitization_status": "",
                "truth_set_eligible": source["truth_set_use"],
                "comparability_group_id": f"MMAE_PMDA_{source['study_id'].replace('-', '')}",
                "source_id": source["source_id"],
                "source_url": PMDA_PADCEV_URL,
                "source_location": source["source_location"],
                "notes": source["notes"] + (" SEM is retained in pmda_mmae_caco2_transport_numeric_v1.csv and intentionally not written into inherited *_sd columns." if source["uncertainty_type"] == "SEM" else ""),
            }
        )
    for source in mdck:
        rows.append(
            {
                "transport_record_id": source["record_id"],
                "entity_type": "released_species",
                "entity_id": source["entity_id"],
                "entity_name": source["compound"],
                "payload_core": "MMAE",
                "transporter_gene": source["transporter_gene"],
                "cell_background": source["cell_system"],
                "transporter_system": "parental or human MDR1-transfected MDCKII",
                "assay_endpoint": "bidirectional apparent permeability; ER=Papp B-A/Papp A-B",
                "concentration_value": source["test_concentration_uM"],
                "concentration_unit": "uM",
                "inhibitor": "none",
                "inhibitor_concentration_value": "",
                "inhibitor_concentration_unit": "",
                "papp_ab_mean": "",
                "papp_ab_sd": "",
                "papp_ba_mean": "",
                "papp_ba_sd": "",
                "efflux_ratio_mean": source["efflux_ratio"],
                "efflux_ratio_sd": "",
                "background_corrected_efflux_ratio": source["background_corrected_ratio_of_ratios"],
                "source_assigned_call": source["source_assigned_call"],
                "evidence_tier": "A_direct_numeric_single_transporter_text_ER",
                "numeric_origin": source["numeric_origin"],
                "digitization_status": "Papp_not_digitized",
                "truth_set_eligible": source["truth_set_use"],
                "comparability_group_id": "MMAE_MDCK_MDR1_1uM_LIU2019",
                "source_id": source["source_id"],
                "source_url": MMAE_MDCK_URL,
                "source_location": source["source_location"],
                "notes": source["inhibitor_condition"],
            }
        )
    for source in bay:
        rows.append(
            {
                "transport_record_id": source["record_id"],
                "entity_type": "released_species",
                "entity_id": source["entity_id"],
                "entity_name": source["compound"],
                "payload_core": "proprietary pyrrole KSP inhibitor",
                "transporter_gene": "ABCB1",
                "cell_background": source["cell_system"],
                "transporter_system": "endogenous Caco-2 source interpretation",
                "assay_endpoint": "bidirectional apparent permeability; ER=Papp B-A/Papp A-B",
                "concentration_value": source["test_concentration_uM"],
                "concentration_unit": "uM",
                "inhibitor": "none",
                "inhibitor_concentration_value": "",
                "inhibitor_concentration_unit": "",
                "papp_ab_mean": source["papp_ab_mean_1e6_cm_s"],
                "papp_ab_sd": "",
                "papp_ba_mean": source["papp_ba_mean_1e6_cm_s"],
                "papp_ba_sd": "",
                "efflux_ratio_mean": source["efflux_ratio"],
                "efflux_ratio_sd": "",
                "background_corrected_efflux_ratio": "",
                "source_assigned_call": source["source_assigned_call"],
                "evidence_tier": "A_direct_numeric_Caco2_structure_unresolved",
                "numeric_origin": "source_main_text_plus_supplement_methods",
                "digitization_status": "",
                "truth_set_eligible": source["truth_set_use"],
                "comparability_group_id": "BAY716_CACO2_2uM",
                "source_id": source["source_id"],
                "source_url": BAY716_URL,
                "source_location": source["source_location"],
                "notes": source["notes"],
            }
        )
    return rows


def release_transport_v5(obi: list[dict], pmda: list[dict], mdck: list[dict], bay: list[dict], zd: list[dict]) -> list[dict]:
    rows = read_csv(PROCESSED / "adc_release_species_transporter_evidence_long_v4.csv")
    names = {"ADC_REL_007": "MMAE", "ADC_REL_009": "SN-38", "ADC_REL_010": "DXd", "ADC_REL_011": "Exatecan", "ADC_REL_012": "BAY-716"}
    protein = {"ABCB1": "P-gp", "ABCG2": "BCRP"}
    for entity_id in ("ADC_REL_009", "ADC_REL_010", "ADC_REL_011"):
        for gene in ("ABCB1", "ABCG2"):
            selected = [row for row in obi if row["entity_id"] == entity_id and row["transporter_gene"] == gene]
            baseline = next(row for row in selected if row["inhibitor"] == "none")
            inhibited = next(row for row in selected if row["inhibitor"] != "none")
            rows.append(
                {
                    "transport_record_id": f"ADC_TRANS_B7_OBI_{entity_id}_{gene}",
                    "release_species_id": entity_id,
                    "transporter_gene": gene,
                    "transporter_protein": protein[gene],
                    "species": "Homo sapiens",
                    "assay_system": "Caco-2 C2BBe1 bidirectional transport with target-selective inhibitor arm",
                    "test_concentrations": "10 uM",
                    "efflux_ratio_no_inhibitor": baseline["efflux_ratio"],
                    "efflux_ratio_with_inhibitor": inhibited["efflux_ratio"],
                    "inhibitor": inhibited["inhibitor"],
                    "source_assigned_call": baseline["source_assigned_call"],
                    "source_binary_label": "",
                    "model_binary_label": "",
                    "model_label_status": "not_assigned_Caco2_inhibitor_attribution_not_gene_isolation",
                    "evidence_tier": "A_direct_numeric_target_inhibitor_attributed_Caco2",
                    "truth_set_candidate": "external_continuous_target_attribution_only",
                    "source_id": "SRC_B7_001",
                    "source_location": "Supplementary Table S5",
                    "curation_note": "Complete Papp and SD are available, but endogenous Caco-2 is not a single-gene transporter system and inhibitor concentrations are not reported.",
                    "release_species_name": names[entity_id],
                    "source_doi": "10.1158/1535-7163.MCT-24-0588",
                    "source_url": OBI_URL,
                    "assay_endpoint": "Papp A-B, Papp B-A and efflux ratio",
                    "control_condition": "no inhibitor",
                    "control_value": baseline["efflux_ratio"],
                    "transporter_condition": "endogenous Caco-2 target arm",
                    "transporter_value": baseline["efflux_ratio"],
                    "inhibitor_condition": inhibited["inhibitor"],
                    "inhibitor_value": inhibited["efflux_ratio"],
                    "effect_ratio": round(float(baseline["efflux_ratio"]) / float(inhibited["efflux_ratio"]), 4),
                    "effect_ratio_definition": "ER(no inhibitor)/ER(target inhibitor)",
                }
            )
    rows.append(
        {
            "transport_record_id": "ADC_TRANS_B7_PMDA_MMAE_ABCB1",
            "release_species_id": "ADC_REL_007",
            "transporter_gene": "ABCB1",
            "transporter_protein": "P-gp",
            "species": "Homo sapiens",
            "assay_system": "Caco-2 Study 8004 direct bidirectional transport",
            "test_concentrations": "1 | 10 | 100 uM; inhibitor run at 10 uM",
            "efflux_ratio_no_inhibitor": "15 | 34 | 49; separate inhibitor-run baseline 28",
            "efflux_ratio_with_inhibitor": "1.3 PSC833 | 2.9 verapamil",
            "inhibitor": "PSC833 10 uM | verapamil 60 uM",
            "source_assigned_call": "P-gp_substrate",
            "source_binary_label": "1",
            "model_binary_label": "",
            "model_label_status": "not_assigned_project_preserves_continuous_source_data",
            "evidence_tier": "A_direct_numeric_regulatory_table",
            "truth_set_candidate": "yes_continuous_independent_direct_ABCB1_evidence",
            "source_id": "SRC_B7_003",
            "source_location": "CTD 2.6.5.15.1; Study 8004",
            "curation_note": "The two 10 uM no-inhibitor rows are distinct source runs and are not pooled.",
            "release_species_name": "MMAE",
            "source_doi": "",
            "source_url": PMDA_PADCEV_URL,
            "assay_endpoint": "Papp A-B, Papp B-A and efflux ratio",
            "control_condition": "10 uM MMAE no inhibitor in inhibitor run",
            "control_value": 28,
            "transporter_condition": "endogenous P-gp Caco-2",
            "transporter_value": 28,
            "inhibitor_condition": "PSC833 | verapamil",
            "inhibitor_value": "1.3 | 2.9",
            "effect_ratio": "21.54 | 9.66",
            "effect_ratio_definition": "ER(no inhibitor)/ER(inhibitor)",
        }
    )
    rows.append(
        {
            "transport_record_id": "ADC_TRANS_B7_PMDA_MMAE_MULTI",
            "release_species_id": "ADC_REL_007",
            "transporter_gene": "ABCB1|ABCG2|ABCC2",
            "transporter_protein": "P-gp|BCRP|MRP2",
            "species": "Homo sapiens",
            "assay_system": "Caco-2 RPT-01709 multi-transporter inhibitor panel",
            "test_concentrations": "10 uM",
            "efflux_ratio_no_inhibitor": 8.8,
            "efflux_ratio_with_inhibitor": "1.6 LY335979 | 6.0 Ko143 | 4.3 indomethacin",
            "inhibitor": "LY335979 10 uM | Ko143 5 uM | indomethacin 100 uM",
            "source_assigned_call": "source_rule_ER_gt_2_indicates_efflux_pump_substrate",
            "source_binary_label": "",
            "model_binary_label": "",
            "model_label_status": "not_assigned_multi_transporter_numeric_not_gene_deconvolved",
            "evidence_tier": "A_direct_numeric_multi_transporter_inhibitor_panel",
            "truth_set_candidate": "external_continuous_target_attribution_only",
            "source_id": "SRC_B7_003",
            "source_location": "CTD 2.6.5.15.2; Study RPT-01709",
            "curation_note": "Selective P-gp inhibition has the dominant effect; Ko143 and indomethacin effects remain partial. Do not create three gene-specific baseline ER labels.",
            "release_species_name": "MMAE",
            "source_doi": "",
            "source_url": PMDA_PADCEV_URL,
            "assay_endpoint": "Papp A-B, Papp B-A and efflux ratio",
            "control_condition": "no inhibitor",
            "control_value": 8.8,
            "transporter_condition": "endogenous multi-transporter Caco-2",
            "transporter_value": 8.8,
            "inhibitor_condition": "three inhibitor arms",
            "inhibitor_value": "1.6 | 6.0 | 4.3",
            "effect_ratio": "5.5 | 1.47 | 2.05",
            "effect_ratio_definition": "ER(no inhibitor)/ER(inhibitor)",
        }
    )
    rows.append(
        {
            "transport_record_id": "ADC_TRANS_B7_MMAE_MDCK_ABCB1",
            "release_species_id": "ADC_REL_007",
            "transporter_gene": "ABCB1",
            "transporter_protein": "P-gp-MDR1",
            "species": "Homo sapiens",
            "assay_system": "MDCK-WT and human MDR1-transfected MDCKII",
            "test_concentrations": "1 uM",
            "efflux_ratio_no_inhibitor": "13.6 WT | 44.5 MDR1",
            "efflux_ratio_with_inhibitor": "<1.5",
            "inhibitor": "quinidine 10 uM | ketoconazole 10 uM",
            "source_assigned_call": "substrate",
            "source_binary_label": "1",
            "model_binary_label": "",
            "model_label_status": "not_assigned_project_preserves_continuous_source_data",
            "evidence_tier": "A_direct_numeric_single_transporter_inhibitor_rescue",
            "truth_set_candidate": "yes_continuous_independent_direct_ABCB1_evidence",
            "source_id": "SRC_B7_004",
            "source_location": "Results; Figure 2A",
            "curation_note": "Exact text ERs and ratio-of-ratios 3.3 are retained; figure-only Papp bars are not digitized.",
            "release_species_name": "MMAE",
            "source_doi": "10.3389/fphar.2019.00749",
            "source_url": MMAE_MDCK_URL,
            "assay_endpoint": "efflux ratio and MDR1/WT ratio-of-ratios",
            "control_condition": "MDCK-WT",
            "control_value": 13.6,
            "transporter_condition": "MDCK-MDR1",
            "transporter_value": 44.5,
            "inhibitor_condition": "quinidine or ketoconazole",
            "inhibitor_value": "<1.5",
            "effect_ratio": 3.3,
            "effect_ratio_definition": "ER(MDCK-MDR1)/ER(MDCK-WT)",
        }
    )
    source = bay[0]
    rows.append(
        {
            "transport_record_id": "ADC_TRANS_B7_BAY716_ABCB1",
            "release_species_id": "ADC_REL_012",
            "transporter_gene": "ABCB1",
            "transporter_protein": "P-gp",
            "species": "Homo sapiens",
            "assay_system": "Caco-2 ACC 169 bidirectional transport",
            "test_concentrations": "2 uM",
            "efflux_ratio_no_inhibitor": source["efflux_ratio"],
            "efflux_ratio_with_inhibitor": "",
            "inhibitor": "none",
            "source_assigned_call": "poor_P-gp_substrate_no_active_efflux",
            "source_binary_label": "0_source_interpretation",
            "model_binary_label": "",
            "model_label_status": "not_assigned_structure_unresolved_and_Caco2_not_gene_isolated",
            "evidence_tier": "A_direct_numeric_Caco2_structure_unresolved",
            "truth_set_candidate": "no_until_structure_and_target_attribution_resolved",
            "source_id": "SRC_B7_005",
            "source_location": source["source_location"],
            "curation_note": source["notes"],
            "release_species_name": "BAY-716",
            "source_doi": "10.3390/cancers12113464",
            "source_url": BAY716_URL,
            "assay_endpoint": "Papp A-B, Papp B-A and efflux ratio",
            "control_condition": "no inhibitor",
            "control_value": 1.6,
            "transporter_condition": "endogenous Caco-2",
            "transporter_value": 1.6,
            "inhibitor_condition": "not tested",
            "inhibitor_value": "",
            "effect_ratio": "",
            "effect_ratio_definition": "",
        }
    )
    for call in zd:
        if call["release_species_id"] == "":
            continue
        rows.append(
            {
                "transport_record_id": call["record_id"],
                "release_species_id": call["release_species_id"],
                "transporter_gene": call["transporter_gene"],
                "transporter_protein": protein[call["transporter_gene"]],
                "species": "Homo sapiens",
                "assay_system": call["assay_system"],
                "test_concentrations": "not_reported",
                "efflux_ratio_no_inhibitor": "",
                "efflux_ratio_with_inhibitor": "",
                "inhibitor": "",
                "source_assigned_call": call["normalized_source_call"],
                "source_binary_label": "",
                "model_binary_label": "",
                "model_label_status": "not_assigned_source_symbol_only",
                "evidence_tier": "B_primary_source_qualitative_no_numeric",
                "truth_set_candidate": "qualitative_corroboration_only",
                "source_id": "SRC_B7_002",
                "source_location": "Supplementary Table S3",
                "curation_note": call["notes"],
                "release_species_name": names[call["release_species_id"]],
                "source_doi": "10.1158/1535-7163.MCT-23-0822",
                "source_url": ZD_URL,
                "assay_endpoint": "source-assigned Caco-2 transporter call",
                "control_condition": "",
                "control_value": "",
                "transporter_condition": call["transporter_gene"],
                "transporter_value": call["normalized_source_call"],
                "inhibitor_condition": "",
                "inhibitor_value": "",
                "effect_ratio": "",
                "effect_ratio_definition": "",
            }
        )
    return rows


def matched_comparisons_v3() -> list[dict]:
    rows = read_csv(PROCESSED / "matched_context_comparisons_v2.csv")
    rows.extend(
        [
            {
                "comparison_id": "MATCH_B7_001",
                "payload_core": "MMAE",
                "transporter_gene": "ABCB1",
                "free_entity_id": "ADC_REL_007",
                "free_entity_name": "MMAE in PMDA Padcev Caco-2 studies",
                "released_entity_id": "ADC_REL_007",
                "released_entity_name": "MMAE in independent MDCK-MDR1 study",
                "free_evidence": "PMDA Study 8004 Caco-2 ER 15, 34 and 49 at 1, 10 and 100 uM; separate 10 uM run ER 28 to 1.3 PSC833 and 2.9 verapamil",
                "released_evidence": "Liu-Kreyche 2019 at 1 uM: MDCK-WT ER 13.6, MDCK-MDR1 ER 44.5, ratio-of-ratios 3.3; quinidine or ketoconazole lowers ER below 1.5",
                "matched_grade": "M1",
                "endpoint_alignment": "same_bidirectional_efflux_ratio_family_with_inhibitor_rescue",
                "cell_background_alignment": "different_Caco2_vs_MDCK_WT_MDR1",
                "concentration_alignment": "overlap_at_1uM_but_different_designs",
                "identity_alignment": "same_exact_MMAE_entity_and_actual_ADC_release_payload",
                "directional_result": "independent_direct_assays_concordantly_support_ABCB1_efflux",
                "cross_assay_fold_computed": "no",
                "training_use": "independent_direct_continuous_validation_not_M0",
                "source_ids": "SRC_B7_003 | SRC_B7_004",
                "limitation": "Independent studies are not the same run/background; PMDA Caco-2 contains endogenous transporters, whereas the MDCK study provides the single-transporter contrast.",
            },
            {
                "comparison_id": "MATCH_B7_002",
                "payload_core": "DXd",
                "transporter_gene": "ABCB1",
                "free_entity_id": "ADC_REL_010",
                "free_entity_name": "DXd / MAAA-1181a PMDA material",
                "released_entity_id": "ADC_REL_010",
                "released_entity_name": "DXd purchased from MedChemExpress in OBI-992 study",
                "free_evidence": "PMDA 1 uM Caco-2 ER 10.7; 4.37 verapamil; 4.44 novobiocin; 4.43 dual inhibitor",
                "released_evidence": "OBI-992 study 10 uM P-gp arm ER 22.1 without inhibitor and 13.8 with verapamil, with complete Papp/SD",
                "matched_grade": "M1H",
                "endpoint_alignment": "same_bidirectional_Caco2_efflux_family",
                "cell_background_alignment": "Caco2_in_both_independent_labs",
                "concentration_alignment": "different_1_vs_10uM",
                "identity_alignment": "same_DXd_name_but_OBI_catalog_and_lot_omitted",
                "directional_result": "same_positive_direction",
                "cross_assay_fold_computed": "no",
                "training_use": "independent_directional_validation_multi_transporter_hold",
                "source_ids": "SRC_B6_001 | SRC_B7_001",
                "limitation": "Both are endogenous Caco-2 panels and neither yields a clean single-gene baseline ER; OBI reports vendor but not catalog/lot.",
            },
            {
                "comparison_id": "MATCH_B7_003",
                "payload_core": "DXd",
                "transporter_gene": "ABCG2",
                "free_entity_id": "ADC_REL_010",
                "free_entity_name": "DXd / MAAA-1181a PMDA material",
                "released_entity_id": "ADC_REL_010",
                "released_entity_name": "DXd purchased from MedChemExpress in OBI-992 study",
                "free_evidence": "PMDA 1 uM Caco-2 ER 10.7; 4.37 verapamil; 4.44 novobiocin; 4.43 dual inhibitor",
                "released_evidence": "OBI-992 study 10 uM BCRP arm ER 17.7 without inhibitor and 9.1 with Ko143, with complete Papp/SD",
                "matched_grade": "M1H",
                "endpoint_alignment": "same_bidirectional_Caco2_efflux_family",
                "cell_background_alignment": "Caco2_in_both_independent_labs",
                "concentration_alignment": "different_1_vs_10uM",
                "identity_alignment": "same_DXd_name_but_OBI_catalog_and_lot_omitted",
                "directional_result": "same_positive_direction",
                "cross_assay_fold_computed": "no",
                "training_use": "independent_directional_validation_multi_transporter_hold",
                "source_ids": "SRC_B6_001 | SRC_B7_001",
                "limitation": "Both are endogenous Caco-2 panels and neither yields a clean single-gene baseline ER; OBI reports vendor but not catalog/lot.",
            },
        ]
    )
    return rows


def update_issue_log() -> list[dict]:
    path = DATA_ROOT / "exclusion_and_conflict_log.csv"
    rows = read_csv(path)
    additions = [
        {
            "issue_id": "ISSUE-021",
            "entity_or_source": "OBI-992 Supplementary Table S5",
            "issue_type": "Caco2_inhibitor_attribution_not_gene_isolation",
            "current_action": "retain_complete_Papp_SD_and_ER_as_continuous_target_attribution; exclude_from_gene_isolated_training_truth",
            "reason": "P-gp and BCRP arms use the same endogenous Caco-2 background; inhibitor concentrations are not reported and residual ER remains substantial for several compounds.",
            "next_evidence_needed": "Single-transporter transfectant or knockout assays with matched Papp and inhibitor conditions",
            "status": "open",
        },
        {
            "issue_id": "ISSUE-022",
            "entity_or_source": "OBI-992 Supplementary Tables S6-S7",
            "issue_type": "complete_ADC_context_not_payload_transport_truth",
            "current_action": "store_payload_and_ADC_RR_in_separate_same_system_context_table",
            "reason": "Whole-ADC IC50 reflects uptake, trafficking, linker cleavage, release and payload efflux, not just substrate transport.",
            "next_evidence_needed": "Direct intracellular released-species concentration or flux measured in the same HCC827 backgrounds",
            "status": "open",
        },
        {
            "issue_id": "ISSUE-023",
            "entity_or_source": "MMAE regulatory Caco-2 studies across vedotin products",
            "issue_type": "possible_study_reuse_or_duplicate_regulatory_summary",
            "current_action": "deduplicate_by_study_number_8004_and_RPT01709_not_by_product_name",
            "reason": "The same MMAE study may be summarized in more than one product review; counting product dossiers as independent experiments would inflate evidence.",
            "next_evidence_needed": "Original study reports or matching study-number audit across dossiers",
            "status": "partially_resolved",
        },
        {
            "issue_id": "ISSUE-024",
            "entity_or_source": "BAY-716",
            "issue_type": "active_metabolite_structure_crosswalk_unresolved",
            "current_action": "retain_numeric_transport_and_exclude_from_structure_model",
            "reason": "The BAY-943 paper names BAY-716 and reports Papp/ER, while a related Bayer paper shows active metabolite 8; neither explicitly equates the two identifiers.",
            "next_evidence_needed": "Author/patent confirmation mapping BAY-716 to an exact structure, formula, CAS or compound number",
            "status": "open",
        },
        {
            "issue_id": "ISSUE-025",
            "entity_or_source": "ZD06519 Supplementary Table S3",
            "issue_type": "source_symbol_call_without_numeric_transport",
            "current_action": "retain_qualitative_corroboration_only_no_model_label",
            "reason": "Check, X and marginal calls are reported without Papp, ER, concentration or inhibitor values.",
            "next_evidence_needed": "Underlying Caco-2 numeric table or author data",
            "status": "open",
        },
        {
            "issue_id": "ISSUE-026",
            "entity_or_source": "S-methyl-DM1 MCF7 efflux",
            "issue_type": "cellular_retention_not_transporter_specific",
            "current_action": "retain_in_membrane_retention_mechanistic_layer_only",
            "reason": "The experiment quantifies radiolabeled cellular washout and intracellular tubulin rebinding but does not perturb or measure P-gp/ABCG2.",
            "next_evidence_needed": "Matched ABCB1/ABCG2 transporter assay for S-methyl-DM1",
            "status": "open",
        },
    ]
    existing = {row["issue_id"] for row in rows}
    rows.extend(row for row in additions if row["issue_id"] not in existing)
    write_csv(path, rows)
    return rows


def raw_download_log() -> list[dict]:
    specs = [
        ("PMC11791482_fulltext.xml", "SRC_B7_001", OBI_URL, "downloaded_primary_fulltext_xml"),
        ("PMC11791482_supplementaryFiles.zip", "SRC_B7_001", OBI_URL, "downloaded_primary_supplement_zip"),
        ("mct-24-0588_supplementary_methods_and_tables_suppsmt.docx", "SRC_B7_001", OBI_URL, "extracted_primary_supplement_docx_rendered_and_visually_checked"),
        ("PMC11063767_fulltext.xml", "SRC_B7_002", ZD_URL, "downloaded_primary_fulltext_xml"),
        ("PMC11063767_supplementaryFiles.zip", "SRC_B7_002", ZD_URL, "downloaded_primary_supplement_zip"),
        ("mct-23-0822_supplementary_table_s3_suppst3.docx", "SRC_B7_002", ZD_URL, "extracted_primary_supplement_docx_rendered_and_visually_checked"),
        ("PMDA_Padcev_MMAE_transport_CTD_2021.pdf", "SRC_B7_003", PMDA_PADCEV_URL, "downloaded_official_pdf_rendered_pages75_77_visually_checked"),
        ("PMC6650582_fulltext.xml", "SRC_B7_004", MMAE_MDCK_URL, "downloaded_primary_fulltext_xml"),
        ("PMC6650582_supplementaryFiles.zip", "SRC_B7_004", MMAE_MDCK_URL, "downloaded_primary_article_assets_zip"),
        ("fphar-10-00749-g002.jpg", "SRC_B7_004", MMAE_MDCK_URL, "downloaded_primary_figure_visually_checked"),
        ("PMC7709048_fulltext.xml", "SRC_B7_005", BAY716_URL, "downloaded_primary_fulltext_xml"),
        ("PMC7709048_supplementaryFiles.zip", "SRC_B7_005", BAY716_URL, "downloaded_primary_supplement_zip"),
        ("cancers-12-03464-g001.jpg", "SRC_B7_005", BAY716_URL, "downloaded_primary_figure_visually_checked"),
        ("bay716_supp_nested/cancers-12-03464-s001.zip", "SRC_B7_005", BAY716_URL, "extracted_nested_supplement_zip"),
        ("bay716_supp_nested/Kirchhoff et al_Appendix A_Suppl methods.docx", "SRC_B7_005", BAY716_URL, "extracted_rendered_all_4_pages_visually_checked"),
        ("bay716_supp_nested/Kirchhoff et al_Appendix B_Suppl figs and tables_201113.docx", "SRC_B7_005", BAY716_URL, "extracted_rendered_all_6_pages_visually_checked"),
        ("Lerchen_2020_BioconjugateChemistry_KSPi_linkers_article.pdf", "SRC_B7_006", KSP_LINKER_URL, "downloaded_primary_article_rendered_all_6_pages_visually_checked"),
        ("PMDA_Adcetris_English_review_2019.pdf", "SRC_B7_007", PMDA_ADCETRIS_URL, "downloaded_official_pdf"),
        ("PMC4324968_fulltext.xml", "SRC_B7_008", SME_DM1_URL, "downloaded_primary_fulltext_xml"),
        ("PMC4324968_supplementaryFiles.zip", "SRC_B7_008", SME_DM1_URL, "downloaded_primary_article_assets_zip"),
        ("pone.0117523.g001.jpg", "SRC_B7_008", SME_DM1_URL, "downloaded_primary_figure_visually_checked"),
    ]
    rows = []
    for relative, source_id, url, status in specs:
        path = RAW7 / relative
        if not path.exists():
            raise FileNotFoundError(path)
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
    failed = RAW7 / "OBI992_MCT240588_supplementary_methods_tables_download_failed_recaptcha.html"
    rows.append(
        {
            "filename": str(failed.relative_to(DATA_ROOT)) if failed.exists() else "",
            "source_id": "SRC_B7_001",
            "source_url": "publisher direct supplementary endpoint",
            "download_status": "failed_recaptcha_html_not_used_as_source",
            "signature_valid": "no_html_challenge",
            "size_bytes": failed.stat().st_size if failed.exists() else "",
            "sha256": sha256_file(failed) if failed.exists() else "",
            "acquired_or_verified_at": RUN_DATE,
        }
    )
    return rows


def main() -> None:
    obi = obi_caco2_rows()
    pmda = pmda_mmae_rows()
    mdck = mmae_mdck_rows()
    bay = bay716_rows()
    zd = zd_calls()
    rr = obi_hcc827_rr_rows()
    rescue = obi_inhibitor_rescue_rows()
    payload_adc = same_system_payload_adc_rows(rr)
    membrane = membrane_retention_rows()
    registry = release_registry_v3()
    direct = direct_numeric_v3(obi, pmda, mdck, bay)
    release = release_transport_v5(obi, pmda, mdck, bay, zd)
    matched = matched_comparisons_v3()
    issues = update_issue_log()
    downloads = raw_download_log()

    write_csv(PROCESSED / "obi992_caco2_payload_transport_numeric_v1.csv", obi)
    write_csv(PROCESSED / "pmda_mmae_caco2_transport_numeric_v1.csv", pmda)
    write_csv(PROCESSED / "mmae_mdck_mdr1_text_numeric_v1.csv", mdck)
    write_csv(PROCESSED / "bay716_caco2_transport_numeric_v1.csv", bay)
    write_csv(PROCESSED / "zd06519_payload_transporter_calls_v1.csv", zd)
    write_csv(PROCESSED / "obi992_hcc827_payload_adc_resistance_v1.csv", rr)
    write_csv(PROCESSED / "obi992_adc_inhibitor_rescue_v1.csv", rescue)
    write_csv(PROCESSED / "same_system_payload_adc_context_v1.csv", payload_adc)
    write_csv(PROCESSED / "s_methyl_dm1_cellular_retention_v1.csv", membrane)
    write_csv(PROCESSED / "adc_release_species_registry_v3.csv", registry)
    write_csv(PROCESSED / "direct_bidirectional_transport_numeric_v3.csv", direct)
    write_csv(PROCESSED / "adc_release_species_transporter_evidence_long_v5.csv", release)
    write_csv(PROCESSED / "matched_context_comparisons_v3.csv", matched)
    write_csv(PROCESSED / "seventh_batch_source_registry.csv", SOURCES)
    write_csv(PROCESSED / "seventh_batch_raw_download_log.csv", downloads)

    id_checks = {
        "obi_record_ids_unique": len({row["record_id"] for row in obi}) == len(obi),
        "pmda_record_ids_unique": len({row["record_id"] for row in pmda}) == len(pmda),
        "direct_transport_ids_unique": len({row["transport_record_id"] for row in direct}) == len(direct),
        "release_transport_ids_unique": len({row["transport_record_id"] for row in release}) == len(release),
        "matched_comparison_ids_unique": len({row["comparison_id"] for row in matched}) == len(matched),
        "release_species_ids_unique": len({row["release_species_id"] for row in registry}) == len(registry),
    }
    model_label_checks = {
        "release_model_binary_labels_all_blank": all(row["model_binary_label"] == "" for row in release),
        "zd_model_binary_labels_all_blank": all(row["model_binary_label"] == "" for row in zd),
    }
    count_checks = {
        "direct_v3_equals_v2_plus_34": len(direct) == len(read_csv(PROCESSED / "direct_bidirectional_transport_numeric_v2.csv")) + 34,
        "obi_s5_has_16_rows": len(obi) == 16,
        "pmda_has_15_rows": len(pmda) == 15,
        "matched_has_one_confirmed_M1": sum(row["matched_grade"] == "M1" for row in matched) == 1,
        "matched_has_zero_M0": sum(row["matched_grade"] == "M0" for row in matched) == 0,
    }
    raw_checks = {
        "all_expected_raw_files_present": all(row["filename"] != "" for row in downloads[:-1]),
        "all_used_raw_files_have_sha256": all(row["sha256"] != "" for row in downloads[:-1]),
        "failed_recaptcha_not_used": downloads[-1]["download_status"] == "failed_recaptcha_html_not_used_as_source",
    }
    validation = {
        "generated_at": RUN_DATE,
        "id_checks": id_checks,
        "model_label_checks": model_label_checks,
        "count_checks": count_checks,
        "raw_checks": raw_checks,
    }
    validation["all_checks_pass"] = all(
        value
        for group in (id_checks, model_label_checks, count_checks, raw_checks)
        for value in group.values()
    )
    if not validation["all_checks_pass"]:
        raise AssertionError(validation)
    (PROCESSED / "seventh_batch_validation.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    released_author_papp_entities = sorted(
        {
            row["entity_id"]
            for row in direct
            if row["entity_type"] == "released_species"
            and row["entity_id"]
            and row["papp_ab_mean"] != ""
            and row["numeric_origin"] not in {"algorithmic_digitization_from_jpeg_bar_mean"}
        }
    )
    summary = {
        "generated_at": RUN_DATE,
        "seventh_batch_sources": len(SOURCES),
        "raw_source_files_or_attempts": len(downloads),
        "new_obi_caco2_numeric_rows": len(obi),
        "new_pmda_mmae_numeric_rows": len(pmda),
        "new_mmae_mdck_exact_text_er_rows": len(mdck),
        "new_bay716_numeric_rows": len(bay),
        "new_zd06519_qualitative_calls": len(zd),
        "direct_numeric_rows_v3": len(direct),
        "direct_numeric_rows_added_vs_v2": len(direct) - len(read_csv(PROCESSED / "direct_bidirectional_transport_numeric_v2.csv")),
        "release_transporter_records_v5": len(release),
        "release_species_registry_v3_entities": len(registry),
        "released_species_with_author_tabulated_papp_count": len(released_author_papp_entities),
        "released_species_with_author_tabulated_papp_ids": released_author_papp_entities,
        "same_system_payload_adc_context_rows": len(payload_adc),
        "matched_comparisons_v3": len(matched),
        "matched_grade_counts": dict(sorted(Counter(row["matched_grade"] for row in matched).items())),
        "m0_exact_same_assay_pairs": sum(row["matched_grade"] == "M0" for row in matched),
        "m1_same_entity_independent_direct_pairs": sum(row["matched_grade"] == "M1" for row in matched),
        "model_binary_labels_created": 0,
        "open_or_partial_issue_count": sum(row["status"] in {"open", "partially_resolved"} for row in issues),
        "key_results": [
            "OBI-992 Table S5 adds complete Papp means/SD and ER for SN-38, DXd and exatecan in P-gp and BCRP inhibitor arms at 10 uM.",
            "PMDA Padcev tables add 15 author-tabulated MMAE/control rows, including Study 8004 and the triplicate RPT-01709 inhibitor panel.",
            "MMAE now has the first confirmed M1 same-entity independent direct comparison: PMDA Caco-2 and Liu-Kreyche MDCK-WT/MDR1.",
            "BAY-716 adds a named active ADC metabolite with Papp A-B 1.8 nm/s, Papp B-A 2.7 nm/s and ER 1.6 at 2 uM, but exact structure remains unresolved.",
            "OBI-992 Tables S6-S7 provide same-system payload/complete-ADC functional context; complete ADC rows are not treated as payload transport truth.",
            "S-methyl-DM1 cellular retention supports a tubulin-rebinding layer for membrane/retention modeling but is not transporter-specific evidence.",
            "No strict M0 free-payload/released-species same-run pair was found; no new model binary labels were created.",
        ],
        "feasibility_conclusion": "The direct ADC-payload evidence layer improves materially and now contains a confirmed M1 independent direct replication for MMAE, but M0 remains absent and Caco-2 target attribution is not gene isolation. The universal ADC substrate-classifier claim gate remains closed; an evidence-layered continuous risk and abstention tool is increasingly feasible.",
    }
    (PROCESSED / "seventh_batch_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
