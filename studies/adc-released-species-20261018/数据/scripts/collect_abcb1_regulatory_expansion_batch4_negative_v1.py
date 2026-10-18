#!/usr/bin/env python3
"""Freeze the fourth, negative-focused ABCB1 external-candidate batch.

This collection uses two routes that were not the main route in earlier batches:
source-traceable multi-compound transporter studies and named metabolites in
regulatory reviews.  It performs structure standardization, connectivity
leakage auditing, evidence-tier adjudication, and source-lineage registration.
It intentionally performs no model scoring, refit, calibration, or threshold
modification.
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
RAW_DIR = PROJECT / "数据/raw/abcb1_regulatory_expansion_batch4_negative"

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
    "regulatory_evidence_pool_v3": (
        PROCESSED / "abcb1_regulatory_evidence_tiered_pool_v3.csv",
        ["standardized_inchikey"],
    ),
    "external_candidate_ledger_v3": (
        PROCESSED / "abcb1_external_collection_candidate_ledger_v3.csv",
        ["standardized_inchikey"],
    ),
}


SOURCES = {
    "FDA_SEPIAPTERIN": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2025/219666Orig1s000IntegratedR.pdf",
        "local_file": "sepiapterin_fda_219666_integrated.pdf",
        "location": "PDF page 188 / printed page 170; extracted text lines 8644-8650",
        "study_ids": "PTC923-2020-002",
        "cluster": "FDA_NDA219666_PTC9232020002_MDR1",
        "source_type": "FDA integrated review",
    },
    "YANG2014_DIAMIDINES": {
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068519/",
        "local_file": "diamidines_pmc4068519.html",
        "location": "Methods: BBB permeability and P-gp-mediated efflux assays; Table 1; Figure 1",
        "study_ids": "PMCID PMC4068519 | PMID 24798280 | DOI 10.1128/AAC.02605-14",
        "cluster": "YANG2014_WT_MDCK_MDR1_MDCK_TABLE1",
        "source_type": "primary multi-compound transporter study",
    },
    "SHEN2020_FREE_DRUG": {
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7138916/",
        "local_file": "free_drug_pmc7138916.html",
        "location": "Methods section 2.3 and Table 3 (n=4; 1 uM)",
        "study_ids": "PMCID PMC7138916 | DOI 10.1002/prp2.575",
        "cluster": "SHEN2020_PARENT_MDCK_MDR1_MDCK_TABLE3",
        "source_type": "primary multi-compound transporter study",
    },
    "FDA_NETUPITANT_METABOLITES": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2014/205718Orig1s000ClinPharmR.pdf",
        "local_file": "netupitant_fda_205718_clinpharm.pdf",
        "location": "PDF page 118 Table 32 / printed page 117; reviewer adjudication on PDF page 119",
        "study_ids": "NDA 205718 MDCKII-MDR1 substrate assessment; sponsor study identifier not exposed",
        "cluster": "FDA_NDA205718_MDCKII_MDR1_METABOLITES_M1_M3",
        "source_type": "FDA clinical pharmacology review",
    },
    "FDA_VANDETANIB": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2011/022405Orig1s000ClinPharmR.pdf",
        "local_file": "vandetanib_fda_022405_clinpharm.pdf",
        "location": "PDF pages 33-34 / printed pages 32-33; Tables 17-18",
        "study_ids": "KMN070",
        "cluster": "FDA_NDA022405_KMN070_MDR1",
        "source_type": "FDA clinical pharmacology review",
    },
    "FDA_KYGEVVI": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2025/219792Orig1s000IntegratedR.pdf",
        "local_file": "kygevvi_fda_219792_integrated.pdf",
        "location": "PDF pages 188-189; NCD3991 repeat transporter study",
        "study_ids": "NCD3991",
        "cluster": "FDA_NDA219792_NCD3991_MDR1",
        "source_type": "FDA integrated review",
    },
    "FDA_FOSDENOPTERIN": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/214018Orig1s000IntegratedR.pdf",
        "local_file": "fosdenopterin_fda_214018_integrated.pdf",
        "location": "PDF pages 91-92",
        "study_ids": "substrate study identifier not exposed",
        "cluster": "FDA_NDA214018_CPMP_MDR1",
        "source_type": "FDA integrated review",
    },
    "FDA_ZURANOLONE": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2023/217369Orig2s000IntegratedR.pdf",
        "local_file": "zuranolone_fda_217369_integrated.pdf",
        "location": "PDF page 167",
        "study_ids": "substrate study identifier not exposed",
        "cluster": "FDA_NDA217369_ZURANOLONE_MDR1",
        "source_type": "FDA integrated review",
    },
    "FDA_ALYFTREK": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2025/218730Orig1s000IntegratedR.pdf",
        "local_file": "alyftrek_fda_218730_integrated.pdf",
        "location": "PDF page 137 / Table 63",
        "study_ids": "substrate study identifier not exposed",
        "cluster": "FDA_NDA218730_VANZACAFTOR_MDR1",
        "source_type": "FDA integrated review",
    },
    "FDA_CEFTOLOZANE": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2014/206829Orig1s000ClinPharmR.pdf",
        "local_file": "ceftolozane_fda_206829_clinpharm.pdf",
        "location": "Caco-2 substrate assessment",
        "study_ids": "substrate study identifier not exposed",
        "cluster": "FDA_NDA206829_CACO2_TRANSPORT",
        "source_type": "FDA clinical pharmacology review",
    },
    "FDA_AVANAFIL": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2012/202276Orig1s000ClinPharmR.pdf",
        "local_file": "avanafil_fda_202276_clinpharm.pdf",
        "location": "MDR1/wild-type substrate summary",
        "study_ids": "substrate study identifier not exposed",
        "cluster": "FDA_NDA202276_AVANAFIL_MDR1",
        "source_type": "FDA clinical pharmacology review",
    },
    "FDA_COPPER_HISTIDINATE": {
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2026/211241Orig1s000IntegratedR.pdf",
        "local_file": "copper_histidinate_fda_211241_integrated.pdf",
        "location": "paired Abcb1KO-MDCKII-MDR1/control substrate study",
        "study_ids": "substrate study identifier not exposed",
        "cluster": "FDA_NDA211241_COPPER_HISTIDINATE_MDR1",
        "source_type": "FDA integrated review",
    },
    "THC2021": {
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8626638/",
        "local_file": "thc_pmc8626638.html",
        "location": "Abcb1KO-MDCKII transporter assays",
        "study_ids": "PMCID PMC8626638",
        "cluster": "THC2021_ABCB1KO_MDCKII_MDR1",
        "source_type": "primary transporter study",
    },
}


DB820_SMILES = "C(N)(=N)C1=CC=C(C=C1)C1=CC=C(O1)C1=NC=C(C(=N)N)C=C1"
DB829_SMILES = "O1C(=CC=C1C1=NC=C(C(=N)N)C=C1)C1=NC=C(C(=N)N)C=C1"


STRICT_CANDIDATES = [
    {
        "entity_name": "Sepiapterin",
        "pubchem_query": "Sepiapterin",
        "source_id": "FDA_SEPIAPTERIN",
        "assay_system": "Abcb1KO-MDCKII-Mock-LV versus Abcb1KO-MDCKII-MDR1",
        "concentration_uM": "multiple; exact values not exposed in review summary",
        "er_transporter": "not separately exposed",
        "er_parental": "matched Mock-LV denominator",
        "net_er": "<2 at all concentrations",
        "positive_control": "digoxin net ER 35.4; valspodar reduced it to 0.98",
        "inhibitor": "valspodar competence condition on digoxin",
        "evidence_tier": "A1_knockout_background_mock_positive_control_inhibitor_confirmed",
        "qc_note": "review explicitly reports all sepiapterin net ER values below 2",
    },
    {
        "entity_name": "Pentamidine",
        "pubchem_query": "Pentamidine",
        "source_id": "YANG2014_DIAMIDINES",
        "assay_system": "WT MDCK versus human MDR1-MDCK; triplicate; TEER and mass-balance QC",
        "concentration_uM": "1 | 10",
        "er_transporter": "1.2 | 3.1",
        "er_parental": "2.2 | 4.9",
        "net_er": "0.55 | 0.63",
        "positive_control": "digoxin net flux ratio 3.7",
        "inhibitor": "not required; matched-parent same-batch positive control",
        "evidence_tier": "A1_matched_wildtype_multiconcentration_positive_control_qc",
        "qc_note": "mass balance 92-99%; TEER >130 ohm cm2",
    },
    {
        "entity_name": "Furamidine (DB75)",
        "pubchem_query": "Furamidine",
        "source_id": "YANG2014_DIAMIDINES",
        "assay_system": "WT MDCK versus human MDR1-MDCK; triplicate; TEER and mass-balance QC",
        "concentration_uM": "1 | 10",
        "er_transporter": "0.5 | 0.9",
        "er_parental": "0.6 | 2.4",
        "net_er": "0.83 | 0.38",
        "positive_control": "digoxin net flux ratio 3.7",
        "inhibitor": "not required; matched-parent same-batch positive control",
        "evidence_tier": "A1_matched_wildtype_multiconcentration_positive_control_qc",
        "qc_note": "mass balance 92-105%; TEER >130 ohm cm2",
    },
    {
        "entity_name": "DB820",
        "input_smiles": DB820_SMILES,
        "structure_source": "Yang 2014 Figure 1; systematic name in Ismail et al. 2003 DOI 10.1021/jm0302602",
        "source_id": "YANG2014_DIAMIDINES",
        "assay_system": "WT MDCK versus human MDR1-MDCK; triplicate; TEER and mass-balance QC",
        "concentration_uM": "1 | 10",
        "er_transporter": "0.8 | 0.7",
        "er_parental": "0.6 | 0.9",
        "net_er": "1.3 | 0.78",
        "positive_control": "digoxin net flux ratio 3.7",
        "inhibitor": "not required; matched-parent same-batch positive control",
        "evidence_tier": "A1_matched_wildtype_multiconcentration_positive_control_qc",
        "qc_note": "mass balance 96-102%; exact structure transcribed from primary figure and systematic name",
    },
    {
        "entity_name": "DB829 (CPD0801)",
        "input_smiles": DB829_SMILES,
        "structure_source": "Yang 2014 Figure 1; symmetric diaza-furamidine structure",
        "source_id": "YANG2014_DIAMIDINES",
        "assay_system": "WT MDCK versus human MDR1-MDCK; triplicate; TEER and mass-balance QC",
        "concentration_uM": "1 | 10",
        "er_transporter": "0.8 | 1.1",
        "er_parental": "0.95 | 1.4",
        "net_er": "0.84 | 0.62",
        "positive_control": "digoxin net flux ratio 3.7",
        "inhibitor": "not required; matched-parent same-batch positive control",
        "evidence_tier": "A1_matched_wildtype_multiconcentration_positive_control_qc",
        "qc_note": "mass balance 84-103%; exact structure transcribed from primary figure",
    },
    {
        "entity_name": "Acetaminophen",
        "pubchem_query": "Acetaminophen",
        "source_id": "SHEN2020_FREE_DRUG",
        "assay_system": "parental MDCK versus human MDR1-MDCK; n=4",
        "concentration_uM": "1",
        "er_transporter": "0.772",
        "er_parental": "0.850299 (14.2/16.7)",
        "net_er": "0.907887",
        "positive_control": "digoxin MDR1 ER 17.1; parent-corrected net ER 13.85",
        "inhibitor": "10 uM cyclosporin A used to suppress intrinsic canine P-gp in parental MDCK context",
        "evidence_tier": "A1_matched_parent_positive_control_numeric_papp",
        "qc_note": "net ER recomputed from published directional Papp values",
    },
    {
        "entity_name": "4-Aminoantipyrine",
        "pubchem_query": "4-Aminoantipyrine",
        "source_id": "SHEN2020_FREE_DRUG",
        "assay_system": "parental MDCK versus human MDR1-MDCK; n=4",
        "concentration_uM": "1",
        "er_transporter": "0.662",
        "er_parental": "0.645833 (27.9/43.2)",
        "net_er": "1.025032",
        "positive_control": "digoxin MDR1 ER 17.1; parent-corrected net ER 13.85",
        "inhibitor": "10 uM cyclosporin A used to suppress intrinsic canine P-gp in parental MDCK context",
        "evidence_tier": "A1_matched_parent_positive_control_numeric_papp",
        "qc_note": "net ER recomputed from published directional Papp values",
    },
    {
        "entity_name": "Theophylline",
        "pubchem_query": "Theophylline",
        "source_id": "SHEN2020_FREE_DRUG",
        "assay_system": "parental MDCK versus human MDR1-MDCK; n=4",
        "concentration_uM": "1",
        "er_transporter": "0.889",
        "er_parental": "0.733333 (16.5/22.5)",
        "net_er": "1.212273",
        "positive_control": "digoxin MDR1 ER 17.1; parent-corrected net ER 13.85",
        "inhibitor": "10 uM cyclosporin A used to suppress intrinsic canine P-gp in parental MDCK context",
        "evidence_tier": "A1_matched_parent_positive_control_numeric_papp",
        "qc_note": "net ER recomputed from published directional Papp values",
    },
    {
        "entity_name": "Ganciclovir",
        "pubchem_query": "Ganciclovir",
        "source_id": "SHEN2020_FREE_DRUG",
        "assay_system": "parental MDCK versus human MDR1-MDCK; n=4",
        "concentration_uM": "1",
        "er_transporter": "1.69",
        "er_parental": "10.627615 (2.54/0.239)",
        "net_er": "0.159021",
        "positive_control": "digoxin MDR1 ER 17.1; parent-corrected net ER 13.85",
        "inhibitor": "10 uM cyclosporin A used to suppress intrinsic canine P-gp in parental MDCK context",
        "evidence_tier": "A1_matched_parent_positive_control_numeric_papp_low_permeability",
        "qc_note": "very low A-to-B parent Papp and high parent ER; retain with low-permeability flag",
    },
    {
        "entity_name": "Pemetrexed",
        "pubchem_query": "Pemetrexed",
        "source_id": "SHEN2020_FREE_DRUG",
        "assay_system": "parental MDCK versus human MDR1-MDCK; n=4",
        "concentration_uM": "1",
        "er_transporter": "0.623",
        "er_parental": "1.754425 (0.793/0.452)",
        "net_er": "0.355101",
        "positive_control": "digoxin MDR1 ER 17.1; parent-corrected net ER 13.85",
        "inhibitor": "10 uM cyclosporin A used to suppress intrinsic canine P-gp in parental MDCK context",
        "evidence_tier": "A1_matched_parent_positive_control_numeric_papp_low_permeability",
        "qc_note": "low-permeability compound; net ER recomputed from published directional Papp values",
    },
    {
        "entity_name": "N-desmethyl netupitant (M1; RO0681133)",
        "pubchem_query": "N-desmethyl netupitant",
        "source_id": "FDA_NETUPITANT_METABOLITES",
        "assay_system": "parental MDCKII versus human MDCKII-MDR1",
        "concentration_uM": "3 | 10 | 30",
        "er_transporter": "0.13+/-0.59 | 0.25+/-0.28 | 2.37+/-0.31",
        "er_parental": "0.17+/-0.64 | 0.13+/-0.58 | 15.07+/-0.32",
        "net_er": "0.74+/-0.87 | 1.86+/-0.64 | 0.16+/-0.44",
        "positive_control": "digoxin net ER 4.00+/-1.03; PSC833 reduced it to 1.28+/-0.07",
        "inhibitor": "PSC833 competence condition on digoxin",
        "evidence_tier": "A1_matched_parent_multiconcentration_positive_control_inhibitor_confirmed",
        "qc_note": "FDA reviewer explicitly adjudicated M1 as non-substrate at all concentrations",
    },
    {
        "entity_name": "Monohydroxy netupitant (M3; RO0731519)",
        "pubchem_query": "Monohydroxy Netupitant",
        "source_id": "FDA_NETUPITANT_METABOLITES",
        "assay_system": "parental MDCKII versus human MDCKII-MDR1",
        "concentration_uM": "3 | 10 | 30",
        "er_transporter": "0.14+/-0.44 | 0.31+/-0.36 | 0.40+/-0.12",
        "er_parental": "5.92+/-0.73 | 0.25+/-0.27 | 0.42+/-0.12",
        "net_er": "0.02+/-0.85 | 1.27+/-0.45 | 0.96+/-0.17",
        "positive_control": "digoxin net ER 4.00+/-1.03; PSC833 reduced it to 1.28+/-0.07",
        "inhibitor": "PSC833 competence condition on digoxin",
        "evidence_tier": "A1_matched_parent_multiconcentration_positive_control_inhibitor_confirmed",
        "qc_note": "FDA reviewer explicitly adjudicated M3 as non-substrate at all concentrations",
    },
]


SECONDARY_CANDIDATES = [
    {
        "entity_name": "Vandetanib",
        "pubchem_query": "Vandetanib",
        "source_id": "FDA_VANDETANIB",
        "assay_system": "MDCKII-control versus human MDCKII-MDR1; study KMN070",
        "concentration_uM": "1 | 10 | 25 | 50",
        "er_transporter": "uncorrected ER 1.4-2.2",
        "er_parental": "matched MDCKII-control",
        "net_er": "1.13 | 1.06 | 0.99 | 0.75",
        "positive_control": "digoxin competence shown in separate study KMN096, not exposed for KMN070 batch",
        "inhibitor": "verapamil condition caused only minor change in vandetanib permeability",
        "evidence_tier": "A2_matched_parent_multiconcentration_missing_same_batch_competence",
        "qc_note": "kept secondary because the locked rule requires same-batch competence",
    }
]


FREE_DRUG_EXCLUSIONS = [
    ("Antipyrine", "0", "0.984727", "prior connectivity overlap: soskuti/q1/broccatelli"),
    ("Lamotrigine", "0", "1.085123", "prior connectivity overlap: soskuti/broccatelli"),
    ("Propranolol", "0", "1.454970", "prior connectivity overlap: soskuti/q1/broccatelli"),
    ("Carbamazepine", "0", "1.093041", "prior connectivity overlap: soskuti/broccatelli"),
    ("Ofloxacin", "0", "1.304310", "prior connectivity overlap: frozen443/soskuti/kelly"),
    ("Atenolol", "0", "0.184590", "prior connectivity overlap: frozen443/soskuti/q1/broccatelli"),
    ("Cimetidine", "0", "0.682278", "prior connectivity overlap: frozen443/soskuti/q1/broccatelli"),
    ("Citalopram", "ambiguous", "2.998165", "raw MDR1 ER <2 but parent-corrected net ER >=2; also prior overlap"),
    ("Diltiazem", "ambiguous", "2.541176", "raw MDR1 ER <2 but parent-corrected net ER >=2; also training overlap"),
    ("Quinidine", "1", "not needed", "same-study ABCB1 positive control; not a negative candidate"),
    ("Fexofenadine", "1", "not needed", "same-study ABCB1 substrate; not a negative candidate"),
    ("Digoxin", "1", "13.85", "same-study ABCB1 positive control; not a negative candidate"),
    ("Sumatriptan", "1", "not needed", "same-study ABCB1 substrate and prior external overlap"),
]


OTHER_EXCLUSIONS = [
    ("Doxecitine", "FDA_KYGEVVI", "secondary_transfectant_only", "repeat assay did not expose a matched parent/mock comparator"),
    ("Doxribtimine", "FDA_KYGEVVI", "secondary_transfectant_only", "repeat assay did not expose a matched parent/mock comparator"),
    ("cPMP (fosdenopterin active moiety)", "FDA_FOSDENOPTERIN", "secondary_missing_competence", "matched comparison exposed but same-batch positive-control/inhibitor competence not exposed"),
    ("Zuranolone", "FDA_ZURANOLONE", "secondary_qualitative_numeric_incomplete", "review narrative does not expose a complete interpretable matched-parent numeric endpoint"),
    ("Vanzacaftor", "FDA_ALYFTREK", "secondary_qualitative_only", "table reports P-gp substrate No without paired numeric assay"),
    ("Ceftolozane", "FDA_CEFTOLOZANE", "secondary_caco2_only", "locked rule excludes Caco-2-only substrate assessment"),
    ("Tazobactam", "FDA_CEFTOLOZANE", "secondary_caco2_only", "locked rule excludes Caco-2-only substrate assessment"),
    ("Avanafil", "FDA_AVANAFIL", "secondary_missing_competence", "single MDR1/wild-type net ER 1.8; same-batch competence not exposed"),
    ("Copper histidinate", "FDA_COPPER_HISTIDINATE", "secondary_structure_domain_unresolved", "paired assay is strong but the exact metal-complex structure is not representation-safe for the current model"),
    ("Secnidazole", "FDA_CEFTOLOZANE", "excluded_prior_external_overlap", "strict assay evidence exists but the structure was already used in Q1 external development"),
    ("Delta-9-THC", "THC2021", "secondary_transfectant_only", "Abcb1KO-MDCKII-MDR1 plus inhibitor but no matched mock/parent comparator"),
    ("11-OH-THC", "THC2021", "secondary_transfectant_only", "Abcb1KO-MDCKII-MDR1 plus inhibitor but no matched mock/parent comparator"),
    ("THC-COOH", "THC2021", "secondary_transfectant_only", "Abcb1KO-MDCKII-MDR1 plus inhibitor but no matched mock/parent comparator"),
]


def resolve_structure(candidate: dict) -> dict:
    if candidate.get("input_smiles"):
        return {
            "pubchem_cid": "",
            "pubchem_inchikey": "",
            "pubchem_smiles": "",
            "molecular_formula": "",
            "molecular_weight": "",
            "xlogp": "",
            "tpsa": "",
            "pubchem_cache_file": "",
            "pubchem_resolution_status": "manual_primary_figure_and_systematic_name",
            "pubchem_resolution_error": "",
            **standardize_smiles(candidate["input_smiles"]),
        }
    query = candidate["pubchem_query"]
    pubchem, status, error = fetch_pubchem(query)
    return {
        **pubchem,
        "pubchem_resolution_status": status,
        "pubchem_resolution_error": error,
        **standardize_smiles(pubchem.get("pubchem_smiles", "")),
    }


def strict_row(row: dict) -> dict:
    return {
        "strict_candidate_id": row["batch4_entity_id"],
        "entity_name": row["entity_name"],
        "directional_label": 0,
        "collection_origin": "regulatory_expansion_batch4_negative",
        "evidence_tier": row["evidence_tier"],
        "assay_system": row["assay_system"],
        "concentration_uM": row["concentration_uM"],
        "er_transporter": row["er_transporter"],
        "er_parental": row["er_parental"],
        "net_er": row["net_er"],
        "inhibitor": row["inhibitor"],
        "inhibitor_er": "",
        "positive_control": row["positive_control"],
        "standardized_inchikey": row["standardized_inchikey"],
        "standardized_isomeric_smiles": row["standardized_isomeric_smiles"],
        "source_id": row["source_id"],
        "source_url": row["source_url"],
        "source_location": row["source_location"],
        "identity_status": row["identity_status"],
        "claim_boundary": row["claim_boundary"],
        "source_study_ids": row["source_study_ids"],
        "upstream_experiment_cluster": row["upstream_experiment_cluster"],
    }


def ledger_row(row: dict) -> dict:
    strict = row["strict_external_candidate_eligible"] == "yes"
    return {
        "collection_candidate_id": f"EXTCOLL::{row['batch4_entity_id']}",
        "entity_name": row["entity_name"],
        "directional_label": 0,
        "evidence_route": "multi_compound_primary_or_regulatory_named_metabolite_matched_control",
        "formal_role_now": "pre_scoring_strict_external_candidate" if strict else "regulatory_secondary_evidence_only",
        "quantitative_endpoint": row["net_er"],
        "endpoint_definition": row["assay_system"],
        "matched_parental_control": "yes",
        "selective_inhibitor_rescue": row["inhibitor"],
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
            "score once under the separately frozen external-scoring protocol"
            if strict
            else "retain outside the formal metric"
        ),
        "claim_boundary": row["claim_boundary"],
    }


def main() -> None:
    prior_indexes = {name: identity_index(path, columns) for name, (path, columns) in PRIOR_DATASETS.items()}
    rows = []
    for index, manual in enumerate(STRICT_CANDIDATES + SECONDARY_CANDIDATES, start=1):
        source = SOURCES[manual["source_id"]]
        structure = resolve_structure(manual)
        key = structure.get("standardized_inchikey", "")
        block = key.split("-")[0] if key else ""
        exact_hits = [name for name, (exact, _) in prior_indexes.items() if key and key in exact]
        connectivity_hits = [name for name, (_, connectivity) in prior_indexes.items() if block and block in connectivity]
        if structure.get("standardization_status") != "ok":
            identity_status = "exclude_unresolved_structure"
        elif exact_hits:
            identity_status = "prior_exact_identity_overlap"
        elif connectivity_hits:
            identity_status = "prior_connectivity_identity_overlap"
        else:
            identity_status = "no_exact_or_connectivity_prior_overlap"
        predeclared_strict = manual in STRICT_CANDIDATES
        eligible = predeclared_strict and identity_status == "no_exact_or_connectivity_prior_overlap"
        local_file = RAW_DIR / source["local_file"]
        rows.append(
            {
                "batch4_entity_id": f"ABCB1_REG_BATCH4_{index:03d}",
                "collection_batch": "regulatory_expansion_batch4_negative",
                "directional_label": 0,
                **manual,
                **structure,
                "inchikey_connectivity_block": block,
                "overlap_exact_prior_sets": " | ".join(exact_hits),
                "overlap_connectivity_prior_sets": " | ".join(connectivity_hits),
                "identity_status": identity_status,
                "strict_external_candidate_eligible": "yes" if eligible else "no",
                "protocol_adjudication": "strict_A1" if eligible else "secondary_A2",
                "source_url": source["url"],
                "source_location": source["location"],
                "source_study_ids": source["study_ids"],
                "upstream_experiment_cluster": source["cluster"],
                "source_lineage_status": "single upstream experiment cluster; cluster-aware resampling required",
                "local_source_file": str(local_file.relative_to(PROJECT)),
                "local_source_sha256": sha256(local_file) if local_file.exists() else "missing",
                "claim_boundary": "pre-scoring ABCB1 non-substrate candidate; not a model-validation result and not evidence for ADC release species",
            }
        )

    write_csv(PROCESSED / "abcb1_regulatory_expansion_batch4_negative_entities_v1.csv", rows)

    source_rows = []
    for source_id, source in SOURCES.items():
        local_file = RAW_DIR / source["local_file"]
        source_entities = [row for row in rows if row["source_id"] == source_id]
        source_rows.append(
            {
                "source_id": source_id,
                "source_type": source["source_type"],
                "source_url": source["url"],
                "source_location": source["location"],
                "source_study_ids": source["study_ids"],
                "upstream_experiment_cluster": source["cluster"],
                "local_source_file": str(local_file.relative_to(PROJECT)),
                "local_source_sha256": sha256(local_file) if local_file.exists() else "missing",
                "entity_rows_in_structured_batch": len(source_entities),
                "strict_identity_clean_count": sum(row["strict_external_candidate_eligible"] == "yes" for row in source_entities),
                "allowed_use": "pre-scoring evidence-tiered ABCB1 external candidate collection",
            }
        )
    write_csv(PROCESSED / "abcb1_regulatory_expansion_batch4_negative_source_registry_v1.csv", source_rows)

    lineage_rows = []
    for source_id, source in SOURCES.items():
        source_entities = [row for row in rows if row["source_id"] == source_id]
        lineage_rows.append(
            {
                "source_id": source_id,
                "entity_names": " | ".join(row["entity_name"] for row in source_entities),
                "source_study_ids": source["study_ids"],
                "upstream_experiment_cluster": source["cluster"],
                "within_batch_cluster_size": len(source_entities),
                "assay_vendor": "not exposed unless stated in primary report",
                "independence_status": "one source cluster regardless of entity count",
                "analysis_rule": "cluster bootstrap/resampling by upstream_experiment_cluster; do not count entities from one table as independent laboratories",
            }
        )
    write_csv(PROCESSED / "abcb1_regulatory_expansion_batch4_negative_source_lineage_v1.csv", lineage_rows)

    exclusion_rows = []
    for name, source_label, net_er, reason in FREE_DRUG_EXCLUSIONS:
        exclusion_rows.append(
            {
                "entity_name": name,
                "source_id": "SHEN2020_FREE_DRUG",
                "reported_or_derived_label": source_label,
                "net_er_or_status": net_er,
                "adjudication": "excluded_from_batch4_strict_negative",
                "reason": reason,
            }
        )
    for name, source_id, adjudication, reason in OTHER_EXCLUSIONS:
        exclusion_rows.append(
            {
                "entity_name": name,
                "source_id": source_id,
                "reported_or_derived_label": "0_candidate",
                "net_er_or_status": "see primary source",
                "adjudication": adjudication,
                "reason": reason,
            }
        )
    exclusion_rows.append(
        {
            "entity_name": "Vandetanib",
            "source_id": "FDA_VANDETANIB",
            "reported_or_derived_label": "0_candidate",
            "net_er_or_status": "1.13 | 1.06 | 0.99 | 0.75",
            "adjudication": "secondary_missing_same_batch_competence",
            "reason": "KMN070 exposes matched-parent multiconcentration data, but the digoxin competence study is KMN096 rather than the same substrate batch",
        }
    )
    write_csv(PROCESSED / "abcb1_regulatory_expansion_batch4_negative_exclusion_log_v1.csv", exclusion_rows)

    tiered_v3 = read_csv(PROCESSED / "abcb1_regulatory_evidence_tiered_pool_v3.csv")
    tiered_v4 = list(tiered_v3) + [
        {"entity_record_id": row["batch4_entity_id"], **row}
        for row in rows
    ]
    write_csv(PROCESSED / "abcb1_regulatory_evidence_tiered_pool_v4.csv", tiered_v4)

    strict_v3 = read_csv(PROCESSED / "abcb1_regulatory_strict_external_pool_v3.csv")
    batch4_strict = [row for row in rows if row["strict_external_candidate_eligible"] == "yes"]
    strict_v4 = list(strict_v3) + [strict_row(row) for row in batch4_strict]
    write_csv(PROCESSED / "abcb1_regulatory_strict_external_pool_v4.csv", strict_v4)

    ledger_v3 = read_csv(PROCESSED / "abcb1_external_collection_candidate_ledger_v3.csv")
    identity_clean_rows = [row for row in rows if row["identity_status"] == "no_exact_or_connectivity_prior_overlap"]
    ledger_v4 = list(ledger_v3) + [ledger_row(row) for row in identity_clean_rows]
    write_csv(PROCESSED / "abcb1_external_collection_candidate_ledger_v4.csv", ledger_v4)

    strict_counts = Counter(str(row["directional_label"]) for row in strict_v4)
    summary = {
        "batch4_source_documents_or_html_snapshots": len(SOURCES),
        "batch4_structured_entity_rows": len(rows),
        "batch4_standardization_ok": sum(row["standardization_status"] == "ok" for row in rows),
        "batch4_identity_clean_rows": len(identity_clean_rows),
        "batch4_new_strict_negative_candidates": len(batch4_strict),
        "batch4_new_strict_negative_names": [row["entity_name"] for row in batch4_strict],
        "batch4_strict_upstream_cluster_count": len({row["upstream_experiment_cluster"] for row in batch4_strict}),
        "batch4_secondary_names": [row["entity_name"] for row in rows if row["strict_external_candidate_eligible"] != "yes"],
        "batch4_exclusion_log_rows": len(exclusion_rows),
        "regulatory_strict_pool_v4_rows": len(strict_v4),
        "regulatory_strict_pool_v4_class_counts": dict(strict_counts),
        "formal_external_shortfall_if_minimum_15_per_class": {
            "positive": max(0, 15 - strict_counts.get("1", 0)),
            "negative": max(0, 15 - strict_counts.get("0", 0)),
        },
        "stopping_rule_reached": strict_counts.get("0", 0) >= 15,
        "source_dependence_warning": "five Shen compounds share one experiment and four diamidines share one experiment; use source-cluster bootstrap and report entity-level plus cluster-level uncertainty",
        "model_scoring_performed": False,
        "model_refit_or_threshold_change_performed": False,
        "formal_status": "negative count gate reached; final source-lineage audit is complete for this batch; scoring remains locked",
    }
    (PROCESSED / "abcb1_regulatory_expansion_batch4_negative_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
