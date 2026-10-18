#!/usr/bin/env python3
"""Collect the second ADC efflux data batch: released species and direct transporter evidence."""

from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw" / "second_batch_release_species"
PROCESSED = ROOT / "processed"
RAW.mkdir(parents=True, exist_ok=True)
PROCESSED.mkdir(parents=True, exist_ok=True)

ACQUIRED_AT = "2026-07-28"
USER_AGENT = "ADC-release-species-curation/1.0"


SOURCES = [
    {
        "source_id": "SRC_B2_001",
        "pmid": "16618769",
        "doi": "10.1158/0008-5472.CAN-05-4489",
        "title": "Antibody-maytansinoid conjugates are activated in targeted cancer cells by lysosomal degradation and linker-dependent intracellular processing",
        "year": 2006,
        "source_type": "primary_article",
        "url": "https://pubmed.ncbi.nlm.nih.gov/16618769/",
        "use_in_batch": "maytansinoid catabolite identity",
    },
    {
        "source_id": "SRC_B2_002",
        "pmid": "20197459",
        "doi": "10.1158/0008-5472.CAN-09-3546",
        "title": "Antibody-maytansinoid conjugates designed to bypass multidrug resistance",
        "year": 2010,
        "source_type": "primary_article",
        "url": "https://pubmed.ncbi.nlm.nih.gov/20197459/",
        "use_in_batch": "Lys-PEG4Mal-DM1 identity and MDR1 comparative retention",
    },
    {
        "source_id": "SRC_B2_003",
        "pmid": "25986849",
        "pmcid": "PMC11024933",
        "doi": "10.1124/dmd.115.064576",
        "title": "Target-Site Drug Metabolism and Transport",
        "year": 2015,
        "source_type": "primary_methods_article",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11024933/",
        "use_in_batch": "Lys-MCC-DM1 lysosomal identity and direct P-gp assay",
    },
    {
        "source_id": "SRC_B2_004",
        "pmid": "26108745",
        "doi": "10.1124/dmd.115.064253",
        "title": "Intracellular Catabolism of an Antibody Drug Conjugate with a Noncleavable Linker",
        "year": 2015,
        "source_type": "primary_article",
        "url": "https://pubmed.ncbi.nlm.nih.gov/26108745/",
        "use_in_batch": "Lys-MCC-DM1 lysosomal LC-MS confirmation",
    },
    {
        "source_id": "SRC_B2_005",
        "pmid": "36741448",
        "pmcid": "PMC9896951",
        "doi": "10.1093/pnasnexus/pgac063",
        "title": "SLC46A3 is a lysosomal proton-coupled steroid conjugate and bile acid transporter involved in transport of active catabolites of T-DM1",
        "year": 2022,
        "source_type": "primary_article",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9896951/",
        "use_in_batch": "Lys-SMCC-DM1 active catabolite and lysosomal export context",
    },
    {
        "source_id": "SRC_B2_006",
        "pmid": "29592884",
        "doi": "10.1158/1535-7163.MCT-17-0672",
        "title": "Mechanistic Modeling of Antibody-Drug Conjugate Internalization at the Cellular Level Reveals Inefficient Processing Steps",
        "year": 2018,
        "source_type": "primary_article",
        "url": "https://pubmed.ncbi.nlm.nih.gov/29592884/",
        "use_in_batch": "quantitative MMAE and cys-mcMMAF intracellular release",
    },
    {
        "source_id": "SRC_B2_007",
        "pmid": "29046337",
        "doi": "10.1158/0008-5472.CAN-17-2391",
        "title": "Cathepsin B Is Dispensable for Cellular Processing of Cathepsin B-Cleavable Antibody-Drug Conjugates",
        "year": 2017,
        "source_type": "primary_article",
        "url": "https://pubmed.ncbi.nlm.nih.gov/29046337/",
        "use_in_batch": "cys-VC-MMAE processing intermediate",
    },
    {
        "source_id": "SRC_B2_008",
        "pmid": "36190986",
        "pmcid": "PMC9812893",
        "doi": "10.1158/1535-7163.MCT-22-0375",
        "title": "Antibody-Drug Conjugate Sacituzumab Govitecan Drives Efficient Tissue Penetration and Rapid Intracellular Drug Release",
        "year": 2023,
        "source_type": "primary_article",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9812893/",
        "use_in_batch": "SN-38 release identity and location",
    },
    {
        "source_id": "SRC_B2_009",
        "pmid": "38643548",
        "doi": "10.1016/j.dmpk.2024.101001",
        "title": "Unveiling the intra-tumor fate of trastuzumab deruxtecan in a xenograft model to support its mechanism of action",
        "year": 2024,
        "source_type": "primary_article",
        "url": "https://pubmed.ncbi.nlm.nih.gov/38643548/",
        "use_in_batch": "DXd release in lysosomal fractions",
    },
    {
        "source_id": "SRC_B2_010",
        "pmid": "41951857",
        "pmcid": "PMC13342044",
        "doi": "10.1007/s10565-026-10183-2",
        "title": "Mechanism of uptake and toxicity of a BCMA antibody drug conjugate with a MMAF payload by nonantigen expressing cells",
        "year": 2026,
        "source_type": "primary_article",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC13342044/",
        "use_in_batch": "quantitative intracellular cys-mcMMAF release",
    },
    {
        "source_id": "SRC_B2_011",
        "pmid": "",
        "doi": "",
        "title": "PMDA Report on the Deliberation Results: belantamab mafodotin",
        "year": 2025,
        "source_type": "regulatory_review_primary_submission_summary",
        "url": "https://www.pmda.go.jp/files/000280862.pdf",
        "use_in_batch": "direct cys-mcMMAF P-gp and BCRP bidirectional transport data",
    },
    {
        "source_id": "SRC_B2_012",
        "pmid": "22875610",
        "pmcid": "PMC3475867",
        "doi": "10.1208/s12248-012-9386-x",
        "title": "ADME of antibody-maytansinoid conjugates",
        "year": 2012,
        "source_type": "review_with_primary_data_context",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3475867/",
        "use_in_batch": "maytansinoid catabolite nomenclature and permeability context",
    },
    {
        "source_id": "SRC_B2_013",
        "pmid": "22829590",
        "doi": "10.1124/dmd.112.046169",
        "title": "In Vitro Characterization of the Drug-Drug Interaction Potential of Catabolites of Antibody-Maytansinoid Conjugates",
        "year": 2012,
        "source_type": "primary_article",
        "url": "https://pubmed.ncbi.nlm.nih.gov/22829590/",
        "use_in_batch": "S-methyl-DM1 and Lys-MCC-DM1 catabolite identity",
    },
]


RELEASE_SPECIES = [
    {
        "release_species_id": "ADC_REL_001",
        "release_species_name": "Lys-MCC-DM1",
        "aliases": "Lys-SMCC-DM1 | lysine-SMCC-DM1 | lysine-MCC-DM1",
        "payload_core": "DM1",
        "residual_group": "lysine-MCC/SMCC",
        "species_class": "amino_acid_linker_catabolite",
        "representative_adc_or_format": "trastuzumab emtansine / Ab-SMCC-DM1",
        "linker_class": "noncleavable thioether",
        "release_identity_status": "confirmed_major_catabolite",
        "release_evidence_method": "radiolabel, LC-MS, isolated lysosome",
        "source_ids": "SRC_B2_001 | SRC_B2_003 | SRC_B2_004 | SRC_B2_005",
        "structure_query_names": "Lysine-MCC-DM1 | Lys-SMCC-DM1",
        "structure_status_override": "manual_review",
        "structure_note_override": "PubChem resolves two records with the same formula but different stereochemical InChIKeys; preserve both aliases until exact experimental stereochemistry is reconciled.",
        "existing_payload_compound_id": "ADC_PAYLOAD_008",
        "modeling_status": "candidate_after_stereochemistry_review",
    },
    {
        "release_species_id": "ADC_REL_002",
        "release_species_name": "Lys-PEG4Mal-DM1",
        "aliases": "lysine-PEG4Mal-DM1",
        "payload_core": "DM1",
        "residual_group": "lysine-PEG4-maleimide",
        "species_class": "hydrophilic_amino_acid_linker_catabolite",
        "representative_adc_or_format": "Ab-PEG4Mal-DM1 experimental format",
        "linker_class": "noncleavable hydrophilic PEG4 thioether",
        "release_identity_status": "confirmed_catabolite",
        "release_evidence_method": "cellular catabolism and metabolite retention",
        "source_ids": "SRC_B2_002",
        "structure_query_names": "Lysine-PEG4Mal-DM1",
        "structure_status_override": "",
        "structure_note_override": "",
        "existing_payload_compound_id": "ADC_PAYLOAD_008",
        "modeling_status": "exclude_until_structure_resolved",
    },
    {
        "release_species_id": "ADC_REL_003",
        "release_species_name": "Lys-SPDB-DM4",
        "aliases": "lysine-SPDB-DM4",
        "payload_core": "DM4",
        "residual_group": "lysine-SPDB",
        "species_class": "reducible_linker_initial_catabolite",
        "representative_adc_or_format": "huC242-SPDB-DM4",
        "linker_class": "cleavable disulfide",
        "release_identity_status": "confirmed_intermediate",
        "release_evidence_method": "cellular catabolism",
        "source_ids": "SRC_B2_001",
        "structure_query_names": "Lys-SPDB-DM4",
        "structure_status_override": "",
        "structure_note_override": "",
        "existing_payload_compound_id": "ADC_PAYLOAD_009",
        "modeling_status": "exclude_until_structure_resolved",
    },
    {
        "release_species_id": "ADC_REL_004",
        "release_species_name": "S-methyl-DM4",
        "aliases": "DM4-SMe | S-methyl ravtansine",
        "payload_core": "DM4",
        "residual_group": "S-methyl",
        "species_class": "downstream_thiol_metabolite",
        "representative_adc_or_format": "huC242-SPDB-DM4",
        "linker_class": "cleavable disulfide",
        "release_identity_status": "confirmed_active_downstream_metabolite",
        "release_evidence_method": "cellular catabolism",
        "source_ids": "SRC_B2_001",
        "structure_query_names": "S-methyl-DM4 | DM4-SMe",
        "structure_status_override": "",
        "structure_note_override": "",
        "existing_payload_compound_id": "ADC_PAYLOAD_009",
        "modeling_status": "exclude_until_structure_resolved",
    },
    {
        "release_species_id": "ADC_REL_005",
        "release_species_name": "S-methyl-DM1",
        "aliases": "DM1-SMe",
        "payload_core": "DM1",
        "residual_group": "S-methyl",
        "species_class": "downstream_thiol_metabolite",
        "representative_adc_or_format": "cleavable antibody-DM1 formats",
        "linker_class": "cleavable disulfide",
        "release_identity_status": "reported_catabolite",
        "release_evidence_method": "in vivo/in vitro catabolite characterization",
        "source_ids": "SRC_B2_012 | SRC_B2_013",
        "structure_query_names": "S-methyl-DM1 | DM1-SMe",
        "structure_status_override": "",
        "structure_note_override": "",
        "existing_payload_compound_id": "ADC_PAYLOAD_008",
        "modeling_status": "candidate_after_source_specificity_review",
    },
    {
        "release_species_id": "ADC_REL_006",
        "release_species_name": "cys-mcMMAF",
        "aliases": "cysteine-maleimidocaproyl-MMAF",
        "payload_core": "MMAF",
        "residual_group": "cysteine-maleimidocaproyl",
        "species_class": "amino_acid_linker_catabolite",
        "representative_adc_or_format": "belantamab mafodotin",
        "linker_class": "noncleavable maleimidocaproyl",
        "release_identity_status": "confirmed_active_catabolite",
        "release_evidence_method": "LC-MS and regulatory bioanalysis",
        "source_ids": "SRC_B2_006 | SRC_B2_010 | SRC_B2_011",
        "structure_query_names": "cys-mcMMAF",
        "structure_status_override": "manual_review",
        "structure_note_override": "PubChem name query returns a small-molecule formula compatible with cys-mcMMAF but the title is Denintuzumab mafodotin; exact entity naming must be checked manually.",
        "existing_payload_compound_id": "ADC_PAYLOAD_018",
        "modeling_status": "candidate_after_identity_review",
    },
    {
        "release_species_id": "ADC_REL_007",
        "release_species_name": "MMAE",
        "aliases": "monomethyl auristatin E",
        "payload_core": "MMAE",
        "residual_group": "none",
        "species_class": "free_payload",
        "representative_adc_or_format": "vc-PABC-MMAE ADC",
        "linker_class": "protease-cleavable self-immolative",
        "release_identity_status": "confirmed_released_payload",
        "release_evidence_method": "quantitative intracellular mass spectrometry",
        "source_ids": "SRC_B2_006",
        "structure_query_names": "Monomethyl auristatin E | MMAE",
        "structure_status_override": "",
        "structure_note_override": "",
        "existing_payload_compound_id": "ADC_PAYLOAD_017",
        "modeling_status": "already_in_free_payload_registry",
    },
    {
        "release_species_id": "ADC_REL_008",
        "release_species_name": "cys-VC-MMAE",
        "aliases": "cysteine-valine-citrulline-MMAE",
        "payload_core": "MMAE",
        "residual_group": "cysteine-VC",
        "species_class": "processing_intermediate",
        "representative_adc_or_format": "anti-HER2-VC(R/S)-MMAE experimental ADC",
        "linker_class": "protease-cleavable peptide",
        "release_identity_status": "confirmed_processing_product",
        "release_evidence_method": "mass spectrometry",
        "source_ids": "SRC_B2_007",
        "structure_query_names": "cys-VC-MMAE",
        "structure_status_override": "",
        "structure_note_override": "",
        "existing_payload_compound_id": "ADC_PAYLOAD_017",
        "modeling_status": "exclude_until_structure_and_activity_role_resolved",
    },
    {
        "release_species_id": "ADC_REL_009",
        "release_species_name": "SN-38",
        "aliases": "7-ethyl-10-hydroxycamptothecin",
        "payload_core": "SN-38",
        "residual_group": "none",
        "species_class": "free_payload",
        "representative_adc_or_format": "sacituzumab govitecan",
        "linker_class": "hydrolysable CL2A",
        "release_identity_status": "confirmed_released_payload",
        "release_evidence_method": "dual-label imaging and payload pharmacodynamics",
        "source_ids": "SRC_B2_008",
        "structure_query_names": "SN-38",
        "structure_status_override": "",
        "structure_note_override": "",
        "existing_payload_compound_id": "ADC_PAYLOAD_024",
        "modeling_status": "already_in_free_payload_registry",
    },
    {
        "release_species_id": "ADC_REL_010",
        "release_species_name": "DXd",
        "aliases": "released deruxtecan payload",
        "payload_core": "DXd",
        "residual_group": "none",
        "species_class": "free_payload",
        "representative_adc_or_format": "trastuzumab deruxtecan",
        "linker_class": "GGFG peptide plus self-immolative spacer",
        "release_identity_status": "confirmed_released_payload",
        "release_evidence_method": "lysosomal fraction protease release",
        "source_ids": "SRC_B2_009",
        "structure_query_names": "DXd",
        "structure_status_override": "resolved_candidate_manual_crosscheck",
        "structure_note_override": "DXd name query resolves PubChem CID 117888634, unlike the ambiguous Deruxtecan name used in batch 1; retain as a derived correction candidate until exact Roth test material is confirmed.",
        "existing_payload_compound_id": "ADC_PAYLOAD_029",
        "modeling_status": "candidate_after_link_to_roth_material_review",
    },
]


RELEASE_OBSERVATIONS = [
    {
        "release_species_id": "ADC_REL_001",
        "adc_or_format": "huC242-SMCC-DM1",
        "cell_or_system": "COLO 205 target cells",
        "duration_h": "",
        "adc_concentration_value": "",
        "adc_concentration_unit": "",
        "measurement_operator": "",
        "measured_value": "",
        "measurement_unit": "",
        "method": "radiolabeled cellular catabolism",
        "observation": "Lys-SMCC-DM1 was the sole detected metabolite from the thioether-linked conjugate.",
        "source_id": "SRC_B2_001",
        "source_location": "abstract and catabolism results",
        "evidence_tier": "A_identity",
    },
    {
        "release_species_id": "ADC_REL_003",
        "adc_or_format": "huC242-SPDB-DM4",
        "cell_or_system": "COLO 205 target cells",
        "duration_h": "",
        "adc_concentration_value": "",
        "adc_concentration_unit": "",
        "measurement_operator": "",
        "measured_value": "",
        "measurement_unit": "",
        "method": "radiolabeled cellular catabolism",
        "observation": "Lys-SPDB-DM4 was generated first and then reduced and S-methylated.",
        "source_id": "SRC_B2_001",
        "source_location": "abstract and catabolism results",
        "evidence_tier": "A_identity",
    },
    {
        "release_species_id": "ADC_REL_004",
        "adc_or_format": "huC242-SPDB-DM4",
        "cell_or_system": "COLO 205 target cells",
        "duration_h": "",
        "adc_concentration_value": "",
        "adc_concentration_unit": "",
        "measurement_operator": "",
        "measured_value": "",
        "measurement_unit": "",
        "method": "radiolabeled cellular catabolism",
        "observation": "S-methyl-DM4 was identified as the lipophilic potent downstream metabolite.",
        "source_id": "SRC_B2_001",
        "source_location": "abstract and catabolism results",
        "evidence_tier": "A_identity",
    },
    {
        "release_species_id": "ADC_REL_002",
        "adc_or_format": "anti-EpCAM-PEG4Mal-DM1",
        "cell_or_system": "COLO 205 and COLO 205MDR",
        "duration_h": "21",
        "adc_concentration_value": "",
        "adc_concentration_unit": "",
        "measurement_operator": "",
        "measured_value": "",
        "measurement_unit": "",
        "method": "cellular processing and intracellular retention",
        "observation": "Lys-PEG4Mal-DM1 was the observed hydrophilic metabolite and was retained better in MDR1-expressing cells than Lys-SMCC-DM1.",
        "source_id": "SRC_B2_002",
        "source_location": "results and discussion; 3 h exposure plus 18 h chase",
        "evidence_tier": "B_identity_functional",
    },
    {
        "release_species_id": "ADC_REL_001",
        "adc_or_format": "example Ab-MCC-DM1",
        "cell_or_system": "HepG2 cells and isolated lysosomes",
        "duration_h": "24 | 48",
        "adc_concentration_value": "",
        "adc_concentration_unit": "",
        "measurement_operator": "",
        "measured_value": "",
        "measurement_unit": "",
        "method": "LC-MS/MS",
        "observation": "Lys-MCC-DM1 was the sole detectable catabolite; lysosomal localization was directly demonstrated.",
        "source_id": "SRC_B2_003",
        "source_location": "Figures 7-8 and accompanying text",
        "evidence_tier": "A_identity",
    },
    {
        "release_species_id": "ADC_REL_001",
        "adc_or_format": "anti-CD70 Ab-MCC-DM1",
        "cell_or_system": "target cells with lysosome enrichment",
        "duration_h": "",
        "adc_concentration_value": "",
        "adc_concentration_unit": "",
        "measurement_operator": "",
        "measured_value": "",
        "measurement_unit": "",
        "method": "lysosome enrichment plus mass spectrometry",
        "observation": "Catabolism to active Lys-MCC-DM1 in the lysosomal compartment was conclusively demonstrated.",
        "source_id": "SRC_B2_004",
        "source_location": "abstract and results",
        "evidence_tier": "A_identity",
    },
    {
        "release_species_id": "ADC_REL_007",
        "adc_or_format": "Ab033-vcMMAE",
        "cell_or_system": "A431 cells",
        "duration_h": "24",
        "adc_concentration_value": "100",
        "adc_concentration_unit": "nM",
        "measurement_operator": ">",
        "measured_value": "5000000",
        "measurement_unit": "drug molecules per cell",
        "method": "quantitative mass spectrometry",
        "observation": "Freed drug from the cleavable linker exceeded five million molecules per A431 cell after 24 h.",
        "source_id": "SRC_B2_006",
        "source_location": "abstract and Figure 4 context",
        "evidence_tier": "A_quantitative_release",
    },
    {
        "release_species_id": "ADC_REL_006",
        "adc_or_format": "Ab033 noncleavable MMAF ADC",
        "cell_or_system": "A431 and H441 cells",
        "duration_h": "0-24",
        "adc_concentration_value": "100",
        "adc_concentration_unit": "nM",
        "measurement_operator": "",
        "measured_value": "",
        "measurement_unit": "drug molecules per cell",
        "method": "quantitative mass spectrometry",
        "observation": "Released cys-mcMMAF was quantified and was lower than freed MMAE from the matched cleavable format.",
        "source_id": "SRC_B2_006",
        "source_location": "Figure 4 and abstract",
        "evidence_tier": "A_identity_quantitative_method",
    },
    {
        "release_species_id": "ADC_REL_008",
        "adc_or_format": "anti-HER2-VC(R/S)-MMAE",
        "cell_or_system": "parental and cathepsin-B-deleted cells",
        "duration_h": "0 | 4 | 24",
        "adc_concentration_value": "",
        "adc_concentration_unit": "",
        "measurement_operator": "",
        "measured_value": "",
        "measurement_unit": "",
        "method": "mass spectrometry in lysate and media",
        "observation": "cys-VC-MMAE processing product was measured after ADC treatment.",
        "source_id": "SRC_B2_007",
        "source_location": "catabolite processing results",
        "evidence_tier": "A_identity",
    },
    {
        "release_species_id": "ADC_REL_009",
        "adc_or_format": "sacituzumab govitecan",
        "cell_or_system": "NCI-N87 cells and xenografts",
        "duration_h": "",
        "adc_concentration_value": "",
        "adc_concentration_unit": "",
        "measurement_operator": "",
        "measured_value": "",
        "measurement_unit": "",
        "method": "dual-label fluorescence and gamma-H2AX pharmacodynamics",
        "observation": "Rapid intracellular release of SN-38 was supported; little to no extracellular tumor release was observed in the dual-label experiment.",
        "source_id": "SRC_B2_008",
        "source_location": "abstract and Figure 6 context",
        "evidence_tier": "B_release_location",
    },
    {
        "release_species_id": "ADC_REL_010",
        "adc_or_format": "trastuzumab deruxtecan",
        "cell_or_system": "tumor-cell lysosomal fractions",
        "duration_h": "",
        "adc_concentration_value": "",
        "adc_concentration_unit": "",
        "measurement_operator": "",
        "measured_value": "",
        "measurement_unit": "",
        "method": "lysosomal fraction cleavage assay",
        "observation": "DXd was released by cysteine proteases including cathepsins in a pH-dependent lysosomal-fraction assay.",
        "source_id": "SRC_B2_009",
        "source_location": "abstract",
        "evidence_tier": "A_identity",
    },
    {
        "release_species_id": "ADC_REL_006",
        "adc_or_format": "belantamab mafodotin",
        "cell_or_system": "live HCEC",
        "duration_h": "48",
        "adc_concentration_value": "10",
        "adc_concentration_unit": "ug/mL",
        "measurement_operator": "equal",
        "measured_value": "21.30",
        "measurement_unit": "fg/cell",
        "method": "LC-MS",
        "observation": "Mean intracellular cys-mcMMAF at the final time point.",
        "source_id": "SRC_B2_010",
        "source_location": "Results, LC-MS release section",
        "evidence_tier": "A_quantitative_release",
    },
    {
        "release_species_id": "ADC_REL_006",
        "adc_or_format": "belantamab mafodotin",
        "cell_or_system": "live HCEC",
        "duration_h": "48",
        "adc_concentration_value": "100",
        "adc_concentration_unit": "ug/mL",
        "measurement_operator": "equal",
        "measured_value": "102.72",
        "measurement_unit": "fg/cell",
        "method": "LC-MS",
        "observation": "Mean intracellular cys-mcMMAF at the final time point.",
        "source_id": "SRC_B2_010",
        "source_location": "Results, LC-MS release section",
        "evidence_tier": "A_quantitative_release",
    },
    {
        "release_species_id": "ADC_REL_006",
        "adc_or_format": "belantamab mafodotin",
        "cell_or_system": "live RPTEC",
        "duration_h": "48",
        "adc_concentration_value": "10",
        "adc_concentration_unit": "ug/mL",
        "measurement_operator": "equal",
        "measured_value": "10.99",
        "measurement_unit": "fg/cell",
        "method": "LC-MS",
        "observation": "Mean intracellular cys-mcMMAF at the final time point.",
        "source_id": "SRC_B2_010",
        "source_location": "Results, LC-MS release section",
        "evidence_tier": "A_quantitative_release",
    },
    {
        "release_species_id": "ADC_REL_006",
        "adc_or_format": "belantamab mafodotin",
        "cell_or_system": "live RPTEC",
        "duration_h": "48",
        "adc_concentration_value": "100",
        "adc_concentration_unit": "ug/mL",
        "measurement_operator": "equal",
        "measured_value": "36.91",
        "measurement_unit": "fg/cell",
        "method": "LC-MS",
        "observation": "Mean intracellular cys-mcMMAF at the final time point.",
        "source_id": "SRC_B2_010",
        "source_location": "Results, LC-MS release section",
        "evidence_tier": "A_quantitative_release",
    },
]


TRANSPORT_EVIDENCE = [
    {
        "transport_record_id": "ADC_TRANS_B2_001",
        "release_species_id": "ADC_REL_006",
        "transporter_gene": "ABCB1",
        "transporter_protein": "P-gp",
        "species": "Homo sapiens",
        "assay_system": "human P-gp-expressing MDCK II bidirectional transport",
        "test_concentrations": "100 nM",
        "efflux_ratio_no_inhibitor": "34",
        "efflux_ratio_with_inhibitor": "1.1",
        "inhibitor": "verapamil 60 uM",
        "source_assigned_call": "substrate",
        "source_binary_label": "1",
        "model_binary_label": "",
        "model_label_status": "not_assigned_source_call_preserved",
        "evidence_tier": "A_direct_numeric_inhibitor_rescue",
        "truth_set_candidate": "yes",
        "source_id": "SRC_B2_011",
        "source_location": "Section 4.5.3, page 32",
        "curation_note": "Regulatory summary of applicant-submitted direct transport experiment.",
    },
    {
        "transport_record_id": "ADC_TRANS_B2_002",
        "release_species_id": "ADC_REL_006",
        "transporter_gene": "ABCG2",
        "transporter_protein": "BCRP",
        "species": "Homo sapiens",
        "assay_system": "human BCRP-expressing MDCK II bidirectional transport",
        "test_concentrations": "0.72 | 7.2 | 14.4 nM",
        "efflux_ratio_no_inhibitor": "<2 at all concentrations",
        "efflux_ratio_with_inhibitor": "",
        "inhibitor": "",
        "source_assigned_call": "not_identified_as_substrate_under_tested_conditions",
        "source_binary_label": "0",
        "model_binary_label": "",
        "model_label_status": "not_assigned_source_call_preserved",
        "evidence_tier": "A_direct_numeric_no_rescue_reported",
        "truth_set_candidate": "yes",
        "source_id": "SRC_B2_011",
        "source_location": "Section 4.5.3, page 32",
        "curation_note": "Source lists cys-mcMMAF among P-gp/MRP/OATP substrates but not BCRP substrates; no BCRP inhibitor result was reported in this summary.",
    },
    {
        "transport_record_id": "ADC_TRANS_B2_003",
        "release_species_id": "ADC_REL_001",
        "transporter_gene": "ABCB1",
        "transporter_protein": "P-gp",
        "species": "Homo sapiens",
        "assay_system": "human MDR1-transfected LLC-PK1 bidirectional permeability",
        "test_concentrations": "0.1 | 1.0 | 5.0 | 10 uM",
        "efflux_ratio_no_inhibitor": "numeric values only in Figure 8",
        "efflux_ratio_with_inhibitor": "numeric values only in Figure 8",
        "inhibitor": "elacridar",
        "source_assigned_call": "no_evidence_as_pgp_substrate_low_potential",
        "source_binary_label": "0",
        "model_binary_label": "",
        "model_label_status": "hold_until_figure_values_dual_curated",
        "evidence_tier": "A_direct_qualitative_inhibitor_control",
        "truth_set_candidate": "pending_numeric_extraction",
        "source_id": "SRC_B2_003",
        "source_location": "Figure 8 and accompanying text",
        "curation_note": "Do not train as a binary record until Figure 8 numeric values are independently extracted and checked.",
    },
    {
        "transport_record_id": "ADC_TRANS_B2_004",
        "release_species_id": "ADC_REL_002",
        "transporter_gene": "ABCB1",
        "transporter_protein": "P-gp",
        "species": "Homo sapiens",
        "assay_system": "COLO 205MDR intracellular metabolite retention and ADC cytotoxicity",
        "test_concentrations": "ADC exposure sufficient for killing; exact value not transcribed",
        "efflux_ratio_no_inhibitor": "",
        "efflux_ratio_with_inhibitor": "",
        "inhibitor": "",
        "source_assigned_call": "better_retained_than_lys_smcc_dm1",
        "source_binary_label": "",
        "model_binary_label": "",
        "model_label_status": "comparative_only_not_transport_truth",
        "evidence_tier": "C_comparative_cellular_retention",
        "truth_set_candidate": "no",
        "source_id": "SRC_B2_002",
        "source_location": "results and discussion",
        "curation_note": "Direct incubation showed no MDR1 inhibition, but inhibition is not the same endpoint as substrate transport.",
    },
    {
        "transport_record_id": "ADC_TRANS_B2_005",
        "release_species_id": "ADC_REL_001",
        "transporter_gene": "ABCB1",
        "transporter_protein": "P-gp",
        "species": "Homo sapiens",
        "assay_system": "COLO 205MDR intracellular metabolite retention and ADC cytotoxicity",
        "test_concentrations": "ADC exposure sufficient for killing; exact value not transcribed",
        "efflux_ratio_no_inhibitor": "",
        "efflux_ratio_with_inhibitor": "",
        "inhibitor": "",
        "source_assigned_call": "less_retained_than_lys_peg4mal_dm1",
        "source_binary_label": "",
        "model_binary_label": "",
        "model_label_status": "comparative_only_not_transport_truth",
        "evidence_tier": "C_comparative_cellular_retention",
        "truth_set_candidate": "no",
        "source_id": "SRC_B2_002",
        "source_location": "results and discussion",
        "curation_note": "This comparison does not override the later direct LLC-PK1 assay; assay context and endpoints differ.",
    },
]


RAW_DOWNLOADS = [
    (
        "Foti_2015_PubMed_metadata.xml",
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=25986849&retmode=xml",
    ),
    (
        "Foti_2015_Figure8_LysMCCDM1_Pgp.jpg",
        "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/140c/11024933/5e41ac4d2471/dmd.115.064576f8.jpg",
    ),
    (
        "Tomabechi_2022_PMC9896951_fulltext.xml",
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9896951/fullTextXML",
    ),
    (
        "Kopp_2023_PubMed_metadata.xml",
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=36190986&retmode=xml",
    ),
    (
        "Newman_2026_PubMed_metadata.xml",
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=41951857&retmode=xml",
    ),
    (
        "PMDA_belantamab_mafodotin_English_report.pdf",
        "https://www.pmda.go.jp/files/000280862.pdf",
    ),
    (
        "Kovtun_2010_article_copy.pdf",
        "https://ptacts.uspto.gov/ptacts/public-informations/petitions/1550303/download-documents?artifactId=Hy520_sEocFVwBEmKcRVP3G3fkn4fsIMKz65ditngVzJTGwJneud85M",
    ),
]

OBSOLETE_GENERATED_FILES = [
    "Foti_2015_PMC11024933.html",
    "Tomabechi_2022_PMC9896951.html",
    "Kopp_2023_PMC9812893.html",
    "Newman_2026_PMC13342044.html",
]


def request_bytes(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fetch_source_metadata() -> list[dict]:
    rows = []
    for source in SOURCES:
        row = dict(source)
        row.setdefault("pmcid", "")
        query = source.get("pmid") or source.get("doi")
        if not query:
            row["metadata_status"] = "not_applicable_regulatory_source"
            rows.append(row)
            continue
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
            {"query": query, "format": "json"}
        )
        try:
            data = json.loads(request_bytes(url))
            results = data.get("resultList", {}).get("result", [])
            if results:
                result = results[0]
                row["metadata_status"] = "resolved"
                row["europe_pmc_id"] = result.get("id", "")
                row["author_string"] = result.get("authorString", "")
                row["journal_title"] = result.get("journalTitle", "")
                row["is_open_access"] = result.get("isOpenAccess", "")
                row["cited_by_count_at_collection"] = result.get("citedByCount", "")
            else:
                row["metadata_status"] = "not_found"
        except Exception as exc:  # noqa: BLE001
            row["metadata_status"] = f"error:{type(exc).__name__}:{exc}"
        rows.append(row)
        time.sleep(0.25)
    return rows


def fetch_pubchem_properties(query_name: str) -> tuple[dict | None, str, str]:
    properties = "Title,CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,MolecularWeight,XLogP,TPSA,Charge"
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        + urllib.parse.quote(query_name, safe="")
        + f"/property/{properties}/JSON"
    )
    try:
        data = json.loads(request_bytes(url))
        records = data.get("PropertyTable", {}).get("Properties", [])
        return (records[0] if records else None), "200", url
    except urllib.error.HTTPError as exc:
        return None, str(exc.code), url
    except Exception as exc:  # noqa: BLE001
        return None, f"error:{type(exc).__name__}:{exc}", url


def add_structures(species_rows: list[dict]) -> tuple[list[dict], list[dict]]:
    resolved_rows = []
    query_log = []
    for row in species_rows:
        result = None
        used_name = ""
        attempts = []
        final_status = ""
        final_url = ""
        for name in row["structure_query_names"].split(" | "):
            result, final_status, final_url = fetch_pubchem_properties(name)
            attempts.append(f"{name}:{final_status}")
            time.sleep(0.27)
            if result:
                used_name = name
                break

        output = dict(row)
        output["pubchem_query_name"] = used_name
        output["pubchem_cid"] = result.get("CID", "") if result else ""
        output["pubchem_title"] = result.get("Title", "") if result else ""
        output["isomeric_smiles"] = result.get("SMILES", result.get("IsomericSMILES", "")) if result else ""
        output["connectivity_smiles"] = (
            result.get("ConnectivitySMILES", result.get("CanonicalSMILES", "")) if result else ""
        )
        output["inchikey"] = result.get("InChIKey", "") if result else ""
        output["molecular_formula"] = result.get("MolecularFormula", "") if result else ""
        output["molecular_weight"] = result.get("MolecularWeight", "") if result else ""
        output["xlogp"] = result.get("XLogP", "") if result else ""
        output["tpsa"] = result.get("TPSA", "") if result else ""
        output["formal_charge"] = result.get("Charge", "") if result else ""
        output["structure_status"] = row["structure_status_override"] or ("resolved" if result else "unresolved")
        output["structure_note"] = row["structure_note_override"] or (
            "" if result else "No PubChem record resolved under the recorded aliases; manual structure curation is required."
        )
        output["structure_source_url"] = final_url
        output.pop("structure_status_override", None)
        output.pop("structure_note_override", None)
        resolved_rows.append(output)
        query_log.append(
            {
                "release_species_id": row["release_species_id"],
                "release_species_name": row["release_species_name"],
                "queries_attempted": " | ".join(attempts),
                "resolved_query_name": used_name,
                "pubchem_cid": output["pubchem_cid"],
                "final_http_status": final_status,
                "final_query_url": final_url,
                "queried_at": ACQUIRED_AT,
            }
        )
    return resolved_rows, query_log


def download_raw_sources() -> list[dict]:
    for filename in OBSOLETE_GENERATED_FILES:
        (RAW / filename).unlink(missing_ok=True)
    rows = []
    for filename, url in RAW_DOWNLOADS:
        path = RAW / filename
        try:
            payload = request_bytes(url)
            suffix = path.suffix.lower()
            valid = (
                (suffix == ".pdf" and payload.startswith(b"%PDF"))
                or (suffix in {".jpg", ".jpeg"} and payload.startswith(b"\xff\xd8"))
                or (suffix == ".xml" and payload.lstrip().startswith(b"<"))
            )
            if not valid:
                raise ValueError(f"unexpected file signature for {suffix}")
            path.write_bytes(payload)
            status = "downloaded_validated_signature"
            detail = ""
        except Exception as exc:  # noqa: BLE001
            payload = b""
            status = "failed"
            detail = f"{type(exc).__name__}:{exc}"
        rows.append(
            {
                "filename": filename,
                "source_url": url,
                "download_status": status,
                "size_bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest() if payload else "",
                "detail": detail,
                "acquired_at": ACQUIRED_AT,
            }
        )
        time.sleep(0.25)
    return rows


def join_reference_fields(
    rows: list[dict],
    species_lookup: dict[str, dict],
    source_lookup: dict[str, dict],
    generated_id_field: str | None = None,
) -> list[dict]:
    outputs = []
    for index, row in enumerate(rows, start=1):
        output = dict(row)
        if generated_id_field:
            output.setdefault(generated_id_field, f"ADC_RELEASE_OBS_{index:03d}")
        species = species_lookup[row["release_species_id"]]
        source = source_lookup[row["source_id"]]
        output["release_species_name"] = species["release_species_name"]
        output["source_doi"] = source.get("doi", "")
        output["source_url"] = source["url"]
        outputs.append(output)
    return outputs


def main() -> None:
    source_rows = fetch_source_metadata()
    source_lookup = {row["source_id"]: row for row in source_rows}
    species_rows, structure_log = add_structures(RELEASE_SPECIES)
    species_lookup = {row["release_species_id"]: row for row in species_rows}
    release_rows = join_reference_fields(
        RELEASE_OBSERVATIONS, species_lookup, source_lookup, generated_id_field="release_observation_id"
    )
    transport_rows = join_reference_fields(TRANSPORT_EVIDENCE, species_lookup, source_lookup)
    raw_log = download_raw_sources()

    write_csv(PROCESSED / "second_batch_source_registry.csv", source_rows)
    write_csv(PROCESSED / "adc_release_species_registry_v1.csv", species_rows)
    write_csv(PROCESSED / "adc_release_observations_long.csv", release_rows)
    write_csv(PROCESSED / "adc_release_species_transporter_evidence_long.csv", transport_rows)
    write_csv(PROCESSED / "second_batch_pubchem_structure_query_log.csv", structure_log)
    write_csv(PROCESSED / "second_batch_raw_download_log.csv", raw_log)

    summary = {
        "generated_at": ACQUIRED_AT,
        "second_batch_sources": len(source_rows),
        "release_species_entities": len(species_rows),
        "release_observation_records": len(release_rows),
        "direct_or_comparative_transporter_records": len(transport_rows),
        "truth_set_candidates": sum(row["truth_set_candidate"] in {"yes", "pending_numeric_extraction"} for row in transport_rows),
        "structures_resolved_or_reviewable": sum(bool(row["pubchem_cid"]) for row in species_rows),
        "structures_unresolved": sum(not row["pubchem_cid"] for row in species_rows),
        "raw_downloads_succeeded": sum(row["download_status"] == "downloaded_validated_signature" for row in raw_log),
        "key_result": "cys-mcMMAF P-gp efflux ratio 34 without inhibitor and 1.1 with verapamil; BCRP efflux ratio <2 at all three tested concentrations.",
        "warning": "Released species are distinct chemical entities from free payloads. Comparative retention and cytotoxicity records are not direct substrate truth.",
    }
    (PROCESSED / "second_batch_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
