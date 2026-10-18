#!/usr/bin/env python3
"""Collect the third ADC efflux batch: structure resolution and free-payload/released-species pairs.

The script deliberately preserves assay heterogeneity. It does not convert any
functional or regulatory source call into a model label, and it does not compute
numeric fold changes across incompatible assay systems.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import time
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "processed"
RAW = ROOT / "raw" / "third_batch_pair_resolution"
RAW.mkdir(parents=True, exist_ok=True)
PROCESSED.mkdir(parents=True, exist_ok=True)

ACQUIRED_AT = "2026-07-28"
USER_AGENT = "ADC-release-species-curation/1.1 (public-data research; contact via local project)"


SOURCES = [
    {
        "source_id": "SRC_B3_001",
        "title": "NCATS GSRS: DM4-SPDB-L-LYSINE",
        "year": "",
        "source_type": "regulatory_substance_registry",
        "url": "https://drugs.ncats.io/substance/6MDF4WB915",
        "use_in_batch": "Lys-SPDB-DM4 identity, formula, registry stereochemistry and identifiers",
        "source_location": "substance record 6MDF4WB915",
        "primary_or_secondary": "primary_registry",
        "notes": "Preferred name DM4-SPDB-L-LYSINE; UNII 6MDF4WB915; CAS 1280215-91-9; PubChem CID 131801096.",
    },
    {
        "source_id": "SRC_B3_002",
        "title": "PubChem CID 131801096: DM4-spdb-L-lysine",
        "year": "",
        "source_type": "chemical_structure_registry",
        "url": "https://pubchem.ncbi.nlm.nih.gov/compound/131801096",
        "use_in_batch": "machine-readable Lys-SPDB-DM4 structure and descriptors",
        "source_location": "PUG REST property record",
        "primary_or_secondary": "derived_registry",
        "notes": "Matched to the NCATS record by PubChem CID and formula.",
    },
    {
        "source_id": "SRC_B3_003",
        "title": "PMDA Polatuzumab Vedotin CTD: nonclinical pharmacokinetics",
        "year": "2021",
        "source_type": "regulatory_primary_submission_document",
        "url": "https://www.pmda.go.jp/drugs/2021/P20210322002/450045000_30300AMX00256_I100_1.pdf",
        "use_in_batch": "Cys-vc-MMAE chemical scheme, exact mass, retention time and excretion fraction",
        "source_location": "PDF page 30, Table 2.6.4.5.3.1-1",
        "primary_or_secondary": "primary_regulatory",
        "notes": "The table reports Cys-vc-MMAE mass 1436.79 Da and total excretion 2.60% dose; MMAE total 19.04% dose.",
    },
    {
        "source_id": "SRC_B3_004",
        "title": "PubChem CID 132471209: Cysteine-vc-mmae",
        "year": "",
        "source_type": "chemical_structure_registry",
        "url": "https://pubchem.ncbi.nlm.nih.gov/compound/132471209",
        "use_in_batch": "machine-readable cys-VC-MMAE structure and descriptors",
        "source_location": "PUG REST property record",
        "primary_or_secondary": "derived_registry",
        "notes": "Structure includes the maleimidocaproyl-valine-citrulline-PABC residual group; cross-checked against the PMDA chemical scheme.",
    },
    {
        "source_id": "SRC_B3_005",
        "title": "FDA Elahere multidisciplinary review",
        "year": "2022",
        "source_type": "regulatory_review",
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2022/761310Orig1s000MultidisciplineR.pdf",
        "use_in_batch": "regulatory source call that DM4 and S-methyl-DM4 are MDR1/P-gp substrates",
        "source_location": "PDF page 77 (printed page 76)",
        "primary_or_secondary": "primary_regulatory_summary",
        "notes": "No numeric transport result is reproduced in the review.",
    },
    {
        "source_id": "SRC_B3_006",
        "title": "NCATS GSRS: LYSINE-MCC-DM1",
        "year": "",
        "source_type": "regulatory_substance_registry",
        "url": "https://drugs.ncats.io/substance/6E9CC1P9UA",
        "use_in_batch": "resolve Lys-MCC-DM1 stereochemical registry conflict",
        "source_location": "substance record 6E9CC1P9UA",
        "primary_or_secondary": "primary_registry",
        "notes": "The material is explicitly EPIMERIC with 11 of 12 stereocenters defined; a single PubChem stereoisomer is not a literal representation of the registered material.",
    },
    {
        "source_id": "SRC_B3_007",
        "title": "PubChem CID 86278355 synonyms: Cys-McMMAF",
        "year": "",
        "source_type": "chemical_structure_registry",
        "url": "https://pubchem.ncbi.nlm.nih.gov/compound/86278355",
        "use_in_batch": "resolve cys-mcMMAF small-molecule identity despite whole-ADC title aliases",
        "source_location": "PUG REST synonym record",
        "primary_or_secondary": "derived_registry",
        "notes": "Synonyms include Cys-McMMAF, CYSTEINYL-MFFA and CAS 1160590-05-5; whole-ADC aliases remain a title-quality warning.",
    },
    {
        "source_id": "SRC_B3_008",
        "title": "Antibody-maytansinoid conjugates designed to bypass multidrug resistance",
        "year": "2010",
        "source_type": "primary_article_reused_from_batch2",
        "url": "https://pubmed.ncbi.nlm.nih.gov/20197459/",
        "use_in_batch": "Lys-PEG4Mal-DM1 primary chemical scheme and S-methyl-DM1 MDR1 functional data",
        "source_location": "Figure 2 and Table 1",
        "primary_or_secondary": "primary_article",
        "notes": "Figure 2 gives a chemical scheme but no machine identifier for Lys-PEG4Mal-DM1; Table 1 reports S-methyl-DM1 IC50 ratios of 8, 7 and 6.",
    },
]


RAW_DOWNLOADS = [
    (
        "FDA_Elahere_2022_MultidisciplineReview.pdf",
        "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2022/761310Orig1s000MultidisciplineR.pdf",
        "pdf",
    ),
    (
        "PMDA_Polatuzumab_2021_CTD.pdf",
        "https://www.pmda.go.jp/drugs/2021/P20210322002/450045000_30300AMX00256_I100_1.pdf",
        "pdf",
    ),
    (
        "PubChem_CID131801096_properties.json",
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/131801096/property/Title,CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,MolecularWeight,XLogP,TPSA,Charge/JSON",
        "json",
    ),
    (
        "PubChem_CID132471209_properties.json",
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/132471209/property/Title,CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,MolecularWeight,XLogP,TPSA,Charge/JSON",
        "json",
    ),
    (
        "PubChem_CID86278355_synonyms.json",
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/86278355/synonyms/JSON",
        "json",
    ),
    (
        "NCATS_6MDF4WB915_LysSPDBDM4.html",
        "https://drugs.ncats.io/substance/6MDF4WB915",
        "html",
    ),
    (
        "NCATS_6E9CC1P9UA_LysMCCDM1.html",
        "https://drugs.ncats.io/substance/6E9CC1P9UA",
        "html",
    ),
]


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def request_bytes(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=90) as response:
        return response.read()


def validate_content(data: bytes, expected: str) -> tuple[bool, str]:
    if expected == "pdf":
        return data.startswith(b"%PDF"), "pdf_signature"
    if expected == "json":
        try:
            json.loads(data.decode("utf-8"))
            return True, "json_parse"
        except (UnicodeDecodeError, json.JSONDecodeError):
            return False, "json_parse_failed"
    if expected == "html":
        lower = data[:4096].lower()
        return b"<html" in lower or b"<!doctype html" in lower, "html_signature"
    return bool(data), "nonempty"


def download_raw() -> list[dict]:
    rows = []
    for filename, url, expected in RAW_DOWNLOADS:
        target = RAW / filename
        status = "failed"
        detail = ""
        data = b""
        for attempt in range(3):
            try:
                data = request_bytes(url)
                valid, method = validate_content(data, expected)
                if not valid:
                    raise ValueError(f"invalid {expected}: {method}")
                target.write_bytes(data)
                status = "downloaded_validated_signature"
                detail = method
                break
            except (urllib.error.URLError, TimeoutError, ValueError) as exc:
                detail = f"attempt {attempt + 1}: {type(exc).__name__}: {exc}"
                time.sleep(1.5 * (attempt + 1))
        rows.append(
            {
                "filename": filename,
                "source_url": url,
                "download_status": status,
                "size_bytes": len(data) if status.startswith("downloaded") else "",
                "sha256": hashlib.sha256(data).hexdigest() if status.startswith("downloaded") else "",
                "detail": detail,
                "acquired_at": ACQUIRED_AT,
            }
        )

    reused = ROOT / "raw" / "second_batch_release_species" / "Kovtun_2010_article_copy.pdf"
    if reused.exists():
        data = reused.read_bytes()
        rows.append(
            {
                "filename": "../second_batch_release_species/Kovtun_2010_article_copy.pdf",
                "source_url": "https://pubmed.ncbi.nlm.nih.gov/20197459/",
                "download_status": "reused_existing_validated_signature",
                "size_bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
                "detail": "PDF signature validated in batch 2; reused without modification",
                "acquired_at": ACQUIRED_AT,
            }
        )
    else:
        rows.append(
            {
                "filename": "../second_batch_release_species/Kovtun_2010_article_copy.pdf",
                "source_url": "https://pubmed.ncbi.nlm.nih.gov/20197459/",
                "download_status": "missing_expected_batch2_file",
                "size_bytes": "",
                "sha256": "",
                "detail": "Run collect_second_batch.py first",
                "acquired_at": ACQUIRED_AT,
            }
        )
    return rows


def pubchem_property(filename: str) -> dict:
    payload = json.loads((RAW / filename).read_text(encoding="utf-8"))
    return payload["PropertyTable"]["Properties"][0]


def build_registry_v2() -> tuple[list[dict], list[dict]]:
    rows = read_csv(PROCESSED / "adc_release_species_registry_v1.csv")
    by_id = {row["release_species_id"]: row for row in rows}

    lys_spdb = pubchem_property("PubChem_CID131801096_properties.json")
    cys_vc = pubchem_property("PubChem_CID132471209_properties.json")

    lys_mcc = by_id["ADC_REL_001"]
    lys_mcc["modeling_status"] = "candidate_after_epimer_handling_rule"
    lys_mcc["structure_status"] = "regulatory_epimeric_mixture"
    lys_mcc["structure_note"] = (
        "NCATS GSRS identifies LYSINE-MCC-DM1 (UNII 6E9CC1P9UA) as EPIMERIC with 11/12 "
        "stereocenters defined. The retained PubChem SMILES/InChIKey is a representative single "
        "epimer, not a literal representation of the registered mixture; use connectivity only "
        "after a preregistered epimer rule."
    )
    lys_mcc["structure_source_url"] = "https://drugs.ncats.io/substance/6E9CC1P9UA | https://pubchem.ncbi.nlm.nih.gov/compound/71587969"

    lys_peg = by_id["ADC_REL_002"]
    lys_peg["structure_status"] = "primary_scheme_found_machine_structure_unresolved"
    lys_peg["structure_note"] = (
        "Kovtun 2010 Figure 2A provides a primary chemical scheme and MS-confirmed synthetic metabolite, "
        "but no verified identifier, formula or machine-readable structure was found. Do not manually "
        "transcribe the drawing into a modeling structure without independent chemical review."
    )
    lys_peg["structure_source_url"] = "https://pubmed.ncbi.nlm.nih.gov/20197459/"
    lys_peg["modeling_status"] = "exclude_until_machine_structure_independently_curated"

    lys_spdb_row = by_id["ADC_REL_003"]
    lys_spdb_row.update(
        {
            "pubchem_query_name": "CID 131801096",
            "pubchem_cid": str(lys_spdb.get("CID", 131801096)),
            "pubchem_title": lys_spdb.get("Title", "DM4-spdb-L-lysine"),
            "isomeric_smiles": lys_spdb.get("SMILES", lys_spdb.get("IsomericSMILES", "")),
            "connectivity_smiles": lys_spdb.get("ConnectivitySMILES", lys_spdb.get("CanonicalSMILES", "")),
            "inchikey": lys_spdb.get("InChIKey", ""),
            "molecular_formula": lys_spdb.get("MolecularFormula", ""),
            "molecular_weight": lys_spdb.get("MolecularWeight", ""),
            "xlogp": lys_spdb.get("XLogP", ""),
            "tpsa": lys_spdb.get("TPSA", ""),
            "formal_charge": lys_spdb.get("Charge", ""),
            "structure_status": "resolved_regulatory_registry",
            "structure_note": "NCATS GSRS UNII 6MDF4WB915 maps DM4-SPDB-L-LYSINE to PubChem CID 131801096; structure and formula agree across registries.",
            "structure_source_url": "https://drugs.ncats.io/substance/6MDF4WB915 | https://pubchem.ncbi.nlm.nih.gov/compound/131801096",
            "modeling_status": "candidate_after_standardization_qc",
        }
    )

    cys_mcmmaf = by_id["ADC_REL_006"]
    cys_mcmmaf["structure_status"] = "resolved_with_synonym_title_conflict"
    cys_mcmmaf["structure_note"] = (
        "PubChem CID 86278355 synonyms include Cys-McMMAF, CYSTEINYL-MFFA and CAS 1160590-05-5, "
        "supporting the small-molecule identity. The preferred PubChem title still conflates the "
        "record with whole-ADC aliases, so the title conflict is retained as a provenance flag."
    )
    cys_mcmmaf["structure_source_url"] = "https://pubchem.ncbi.nlm.nih.gov/compound/86278355"
    cys_mcmmaf["modeling_status"] = "candidate_after_synonym_identity_qc"

    cys_vc_row = by_id["ADC_REL_008"]
    cys_vc_row["aliases"] = "cysteine-valine-citrulline-MMAE | Cys-MC-VC-PAB-MMAE | Cysteine-vc-mmae"
    cys_vc_row["residual_group"] = "cysteine-maleimidocaproyl-valine-citrulline-PABC"
    cys_vc_row.update(
        {
            "pubchem_query_name": "CID 132471209",
            "pubchem_cid": str(cys_vc.get("CID", 132471209)),
            "pubchem_title": cys_vc.get("Title", "Cysteine-vc-mmae"),
            "isomeric_smiles": cys_vc.get("SMILES", cys_vc.get("IsomericSMILES", "")),
            "connectivity_smiles": cys_vc.get("ConnectivitySMILES", cys_vc.get("CanonicalSMILES", "")),
            "inchikey": cys_vc.get("InChIKey", ""),
            "molecular_formula": cys_vc.get("MolecularFormula", ""),
            "molecular_weight": cys_vc.get("MolecularWeight", ""),
            "xlogp": cys_vc.get("XLogP", ""),
            "tpsa": cys_vc.get("TPSA", ""),
            "formal_charge": cys_vc.get("Charge", ""),
            "structure_status": "resolved_candidate_primary_scheme_plus_pubchem",
            "structure_note": "PubChem CID 132471209 (Cysteine-vc-mmae) matches the PMDA-reported mass/formula and chemical scheme; keep the cys-VC shorthand mapping explicit because the structure retains mc-VC-PABC.",
            "structure_source_url": "https://www.pmda.go.jp/drugs/2021/P20210322002/450045000_30300AMX00256_I100_1.pdf | https://pubchem.ncbi.nlm.nih.gov/compound/132471209",
            "modeling_status": "candidate_after_source_structure_crosscheck",
        }
    )

    updates = [
        {
            "release_species_id": "ADC_REL_001",
            "release_species_name": "Lys-MCC-DM1",
            "prior_structure_status": "manual_review",
            "new_structure_status": "regulatory_epimeric_mixture",
            "resolution_outcome": "stereochemical_conflict_explained_not_collapsed",
            "registry_or_pubchem_id": "UNII 6E9CC1P9UA | PubChem CID 71587969",
            "molecular_formula": "C53H75ClN6O15S",
            "molecular_weight": "1103.712",
            "xlogp": "-0.2",
            "tpsa": "312",
            "machine_structure_available": "representative_single_epimer_only",
            "modeling_action": "use only after epimer/connectivity rule",
            "primary_support": "NCATS GSRS explicitly marks EPIMERIC; 11/12 stereocenters defined",
            "source_ids": "SRC_B3_006",
        },
        {
            "release_species_id": "ADC_REL_002",
            "release_species_name": "Lys-PEG4Mal-DM1",
            "prior_structure_status": "unresolved",
            "new_structure_status": "primary_scheme_found_machine_structure_unresolved",
            "resolution_outcome": "scheme_found_but_not_machine_resolved",
            "registry_or_pubchem_id": "",
            "molecular_formula": "",
            "molecular_weight": "",
            "xlogp": "",
            "tpsa": "",
            "machine_structure_available": "no",
            "modeling_action": "exclude",
            "primary_support": "Kovtun 2010 Figure 2A; identity confirmed by electrospray MS",
            "source_ids": "SRC_B3_008",
        },
        {
            "release_species_id": "ADC_REL_003",
            "release_species_name": "Lys-SPDB-DM4",
            "prior_structure_status": "unresolved",
            "new_structure_status": "resolved_regulatory_registry",
            "resolution_outcome": "resolved",
            "registry_or_pubchem_id": "UNII 6MDF4WB915 | PubChem CID 131801096",
            "molecular_formula": lys_spdb.get("MolecularFormula", ""),
            "molecular_weight": lys_spdb.get("MolecularWeight", ""),
            "xlogp": lys_spdb.get("XLogP", ""),
            "tpsa": lys_spdb.get("TPSA", ""),
            "machine_structure_available": "yes",
            "modeling_action": "candidate after standardization QC",
            "primary_support": "NCATS GSRS identity plus PubChem structure",
            "source_ids": "SRC_B3_001 | SRC_B3_002",
        },
        {
            "release_species_id": "ADC_REL_006",
            "release_species_name": "cys-mcMMAF",
            "prior_structure_status": "manual_review",
            "new_structure_status": "resolved_with_synonym_title_conflict",
            "resolution_outcome": "small_molecule_identity_supported_title_flag_retained",
            "registry_or_pubchem_id": "PubChem CID 86278355 | CAS 1160590-05-5",
            "molecular_formula": cys_mcmmaf["molecular_formula"],
            "molecular_weight": cys_mcmmaf["molecular_weight"],
            "xlogp": cys_mcmmaf["xlogp"],
            "tpsa": cys_mcmmaf["tpsa"],
            "machine_structure_available": "yes",
            "modeling_action": "candidate after synonym identity QC",
            "primary_support": "PubChem synonym record contains Cys-McMMAF and CYSTEINYL-MFFA",
            "source_ids": "SRC_B3_007",
        },
        {
            "release_species_id": "ADC_REL_008",
            "release_species_name": "cys-VC-MMAE",
            "prior_structure_status": "unresolved",
            "new_structure_status": "resolved_candidate_primary_scheme_plus_pubchem",
            "resolution_outcome": "resolved_with_shorthand_mapping_flag",
            "registry_or_pubchem_id": "PubChem CID 132471209 | CAS 3051855-32-1",
            "molecular_formula": cys_vc.get("MolecularFormula", ""),
            "molecular_weight": cys_vc.get("MolecularWeight", ""),
            "xlogp": cys_vc.get("XLogP", ""),
            "tpsa": cys_vc.get("TPSA", ""),
            "machine_structure_available": "yes",
            "modeling_action": "candidate after source-structure crosscheck",
            "primary_support": "PMDA chemical scheme/mass plus PubChem structure",
            "source_ids": "SRC_B3_003 | SRC_B3_004",
        },
    ]
    return rows, updates


def build_observations_v2() -> list[dict]:
    rows = read_csv(PROCESSED / "adc_release_observations_long.csv")
    additions = [
        {
            "release_species_id": "ADC_REL_008",
            "adc_or_format": "radiolabeled polatuzumab vedotin",
            "cell_or_system": "rat urine and bile metabolite profile",
            "duration_h": "",
            "adc_concentration_value": "",
            "adc_concentration_unit": "",
            "measurement_operator": "equal",
            "measured_value": "2.60",
            "measurement_unit": "% administered radioactive dose",
            "method": "radiochromatographic metabolite profiling",
            "observation": "Cys-vc-MMAE: bile 2.60% dose, urine listed as NA in the table, total 2.60% dose; exact mass 1436.79 Da and retention time 34.3 min.",
            "source_id": "SRC_B3_003",
            "source_location": "PDF page 30, Table 2.6.4.5.3.1-1",
            "evidence_tier": "A_regulatory_quantitative_metabolite",
            "release_observation_id": "ADC_RELEASE_OBS_016",
            "release_species_name": "cys-VC-MMAE",
            "source_doi": "",
            "source_url": "https://www.pmda.go.jp/drugs/2021/P20210322002/450045000_30300AMX00256_I100_1.pdf",
        },
        {
            "release_species_id": "ADC_REL_007",
            "adc_or_format": "radiolabeled polatuzumab vedotin",
            "cell_or_system": "rat urine and bile metabolite profile",
            "duration_h": "",
            "adc_concentration_value": "",
            "adc_concentration_unit": "",
            "measurement_operator": "equal",
            "measured_value": "19.04",
            "measurement_unit": "% administered radioactive dose",
            "method": "radiochromatographic metabolite profiling",
            "observation": "MMAE: urine 2.56% dose and bile 16.48% dose, total 19.04% dose; exact mass 717.50 Da and retention time 21.6 min.",
            "source_id": "SRC_B3_003",
            "source_location": "PDF page 30, Table 2.6.4.5.3.1-1",
            "evidence_tier": "A_regulatory_quantitative_metabolite",
            "release_observation_id": "ADC_RELEASE_OBS_017",
            "release_species_name": "MMAE",
            "source_doi": "",
            "source_url": "https://www.pmda.go.jp/drugs/2021/P20210322002/450045000_30300AMX00256_I100_1.pdf",
        },
    ]
    return rows + additions


TRANSPORT_EXTRA_FIELDS = [
    "assay_endpoint",
    "control_condition",
    "control_value",
    "transporter_condition",
    "transporter_value",
    "inhibitor_condition",
    "inhibitor_value",
    "effect_ratio",
    "effect_ratio_definition",
]


def build_transport_v2() -> list[dict]:
    rows = read_csv(PROCESSED / "adc_release_species_transporter_evidence_long.csv")
    for row in rows:
        for field in TRANSPORT_EXTRA_FIELDS:
            row[field] = ""

    common = {
        "release_species_id": "ADC_REL_005",
        "transporter_gene": "ABCB1",
        "transporter_protein": "P-gp/MDR1",
        "species": "Homo sapiens",
        "test_concentrations": "cytotoxicity concentration-response",
        "efflux_ratio_no_inhibitor": "",
        "efflux_ratio_with_inhibitor": "",
        "inhibitor": "",
        "source_assigned_call": "MDR1_susceptibility_supported",
        "source_binary_label": "1",
        "model_binary_label": "",
        "model_label_status": "not_assigned_functional_endpoint_preserved",
        "evidence_tier": "B_functional_resistance_or_inhibitor_rescue",
        "truth_set_candidate": "yes_functional_only",
        "source_id": "SRC_B3_008",
        "source_location": "Table 1",
        "curation_note": "Functional cytotoxicity evidence; not a bidirectional permeability measurement and not numerically comparable with Roth Table 3.",
        "release_species_name": "S-methyl-DM1",
        "source_doi": "10.1158/0008-5472.CAN-09-3546",
        "source_url": "https://pubmed.ncbi.nlm.nih.gov/20197459/",
        "assay_endpoint": "IC50 cytotoxicity",
    }
    additions = []
    comparisons = [
        (
            "ADC_TRANS_B3_001",
            "COLO 205MDR versus parental COLO 205",
            "COLO 205 parental",
            "0.05 nM",
            "COLO 205MDR",
            "0.42 nM",
            "",
            "",
            "8",
            "IC50(COLO 205MDR) / IC50(COLO 205)",
        ),
        (
            "ADC_TRANS_B3_002",
            "HCT-15 with and without cyclosporin A",
            "HCT-15 without cyclosporin A",
            "0.14 nM",
            "",
            "",
            "HCT-15 plus cyclosporin A 1 uM",
            "0.02 nM",
            "7",
            "IC50(-cyclosporin A) / IC50(+cyclosporin A)",
        ),
        (
            "ADC_TRANS_B3_003",
            "UO-31 with and without cyclosporin A",
            "UO-31 without cyclosporin A",
            "0.40 nM",
            "",
            "",
            "UO-31 plus cyclosporin A 1 uM",
            "0.07 nM",
            "6",
            "IC50(-cyclosporin A) / IC50(+cyclosporin A)",
        ),
    ]
    for record_id, system, control, control_value, trans, trans_value, inhibitor, inhibitor_value, ratio, definition in comparisons:
        row = dict(common)
        row.update(
            {
                "transport_record_id": record_id,
                "assay_system": system,
                "control_condition": control,
                "control_value": control_value,
                "transporter_condition": trans,
                "transporter_value": trans_value,
                "inhibitor_condition": inhibitor,
                "inhibitor_value": inhibitor_value,
                "effect_ratio": ratio,
                "effect_ratio_definition": definition,
            }
        )
        additions.append(row)

    additions.append(
        {
            "transport_record_id": "ADC_TRANS_B3_004",
            "release_species_id": "ADC_REL_004",
            "transporter_gene": "ABCB1",
            "transporter_protein": "P-gp/MDR1",
            "species": "Homo sapiens",
            "assay_system": "in vitro transporter assessment summarized in FDA multidisciplinary review",
            "test_concentrations": "not reproduced in review",
            "efflux_ratio_no_inhibitor": "not reproduced in review",
            "efflux_ratio_with_inhibitor": "",
            "inhibitor": "",
            "source_assigned_call": "substrate",
            "source_binary_label": "1",
            "model_binary_label": "",
            "model_label_status": "not_assigned_regulatory_call_only",
            "evidence_tier": "B_regulatory_call_no_numeric",
            "truth_set_candidate": "yes_call_only",
            "source_id": "SRC_B3_005",
            "source_location": "PDF page 77 (printed page 76)",
            "curation_note": "FDA review states that DM4 and S-methyl-DM4 are MDR1/P-gp substrates and neither inhibits P-gp; numeric assay output is not shown.",
            "release_species_name": "S-methyl-DM4",
            "source_doi": "",
            "source_url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2022/761310Orig1s000MultidisciplineR.pdf",
            "assay_endpoint": "regulatory substrate call",
            "control_condition": "",
            "control_value": "",
            "transporter_condition": "",
            "transporter_value": "",
            "inhibitor_condition": "",
            "inhibitor_value": "",
            "effect_ratio": "",
            "effect_ratio_definition": "",
        }
    )
    return rows + additions


def numeric(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def delta(a: str, b: str) -> str:
    left, right = numeric(a), numeric(b)
    if left is None or right is None:
        return ""
    return f"{right - left:.2f}"


def build_pair_table(registry: list[dict]) -> list[dict]:
    payloads = {row["compound_id"]: row for row in read_csv(PROCESSED / "compound_registry_v1.csv")}
    released = {row["release_species_id"]: row for row in registry}
    specs = [
        {
            "pair_id": "ADC_PAIR_001",
            "payload_id": "ADC_PAYLOAD_018",
            "release_id": "ADC_REL_006",
            "transporter_gene": "ABCB1",
            "free_endpoint": "HEK293 transporter/control cytotoxicity RR",
            "free_effect_value": "2.0",
            "free_effect_definition": "IC50(MDR19) / IC50(pcDNA)",
            "release_endpoint": "MDCK-II bidirectional efflux ratio with inhibitor rescue",
            "release_effect_value": "34; 1.1 with verapamil",
            "release_effect_definition": "B-to-A / A-to-B; no inhibitor and +verapamil",
            "directional_comparison": "discordant_higher_efflux_signal_after_residual_group",
            "release_source_ids": "SRC_B2_011",
            "note": "Strongest current directionality example, but cross-assay and not a matched causal experiment.",
        },
        {
            "pair_id": "ADC_PAIR_002",
            "payload_id": "ADC_PAYLOAD_008",
            "release_id": "ADC_REL_001",
            "transporter_gene": "ABCB1",
            "free_endpoint": "HEK293 transporter/control cytotoxicity RR",
            "free_effect_value": "52",
            "free_effect_definition": "IC50(MDR19) / IC50(pcDNA)",
            "release_endpoint": "LLC-PK1 bidirectional permeability, qualitative conclusion",
            "release_effect_value": "no evidence as P-gp substrate; numeric Figure 8 pending",
            "release_effect_definition": "source conclusion with elacridar control",
            "directional_comparison": "discordant_candidate_lower_efflux_signal_after_residual_group",
            "release_source_ids": "SRC_B2_003",
            "note": "Do not convert the qualitative conclusion to a matched RR or digitize Figure 8 only once.",
        },
        {
            "pair_id": "ADC_PAIR_003",
            "payload_id": "ADC_PAYLOAD_008",
            "release_id": "ADC_REL_005",
            "transporter_gene": "ABCB1",
            "free_endpoint": "HEK293 transporter/control cytotoxicity RR",
            "free_effect_value": "52",
            "free_effect_definition": "IC50(MDR19) / IC50(pcDNA)",
            "release_endpoint": "resistant/parental or inhibitor-rescue cytotoxicity ratios",
            "release_effect_value": "8 | 7 | 6",
            "release_effect_definition": "COLO205MDR/parental; HCT-15 -/+cycA; UO-31 -/+cycA",
            "directional_comparison": "same_positive_direction_lower_functional_ratio_cross_assay",
            "release_source_ids": "SRC_B3_008",
            "note": "Supports retained MDR1 susceptibility; magnitudes are not comparable across cell lines and designs.",
        },
        {
            "pair_id": "ADC_PAIR_004",
            "payload_id": "ADC_PAYLOAD_008",
            "release_id": "ADC_REL_002",
            "transporter_gene": "ABCB1",
            "free_endpoint": "HEK293 transporter/control cytotoxicity RR",
            "free_effect_value": "52",
            "free_effect_definition": "IC50(MDR19) / IC50(pcDNA)",
            "release_endpoint": "comparative intracellular retention",
            "release_effect_value": "better retained than Lys-SMCC-DM1",
            "release_effect_definition": "relative cellular retention after ADC processing",
            "directional_comparison": "insufficient_for_substrate_direction",
            "release_source_ids": "SRC_B2_002 | SRC_B3_008",
            "note": "Comparative retention is mechanism support only; machine-readable structure remains unresolved.",
        },
        {
            "pair_id": "ADC_PAIR_005",
            "payload_id": "ADC_PAYLOAD_009",
            "release_id": "ADC_REL_004",
            "transporter_gene": "ABCB1",
            "free_endpoint": "HEK293 transporter/control cytotoxicity RR",
            "free_effect_value": "96",
            "free_effect_definition": "IC50(MDR19) / IC50(pcDNA)",
            "release_endpoint": "FDA regulatory substrate call",
            "release_effect_value": "P-gp substrate; no numeric value reproduced",
            "release_effect_definition": "source call",
            "directional_comparison": "same_positive_direction_no_released_numeric",
            "release_source_ids": "SRC_B3_005",
            "note": "Useful as call-only external evidence, not as a continuous endpoint.",
        },
        {
            "pair_id": "ADC_PAIR_006",
            "payload_id": "ADC_PAYLOAD_009",
            "release_id": "ADC_REL_003",
            "transporter_gene": "ABCB1",
            "free_endpoint": "HEK293 transporter/control cytotoxicity RR",
            "free_effect_value": "96",
            "free_effect_definition": "IC50(MDR19) / IC50(pcDNA)",
            "release_endpoint": "no released-species transporter evidence located",
            "release_effect_value": "",
            "release_effect_definition": "",
            "directional_comparison": "unknown",
            "release_source_ids": "",
            "note": "Structure is now resolved; transporter evidence remains a collection gap.",
        },
        {
            "pair_id": "ADC_PAIR_007",
            "payload_id": "ADC_PAYLOAD_017",
            "release_id": "ADC_REL_008",
            "transporter_gene": "ABCB1",
            "free_endpoint": "HEK293 transporter/control cytotoxicity RR",
            "free_effect_value": "400",
            "free_effect_definition": "IC50(MDR19) / IC50(pcDNA)",
            "release_endpoint": "no released-species transporter evidence located",
            "release_effect_value": "",
            "release_effect_definition": "",
            "directional_comparison": "unknown",
            "release_source_ids": "",
            "note": "Structure is now resolved as the larger Cys-MC-VC-PABC-MMAE species; do not copy the free MMAE label.",
        },
    ]

    rows = []
    for spec in specs:
        free = payloads[spec["payload_id"]]
        rel = released[spec["release_id"]]
        rows.append(
            {
                "pair_id": spec["pair_id"],
                "payload_core": free["paper_compound_name"],
                "free_payload_compound_id": spec["payload_id"],
                "free_payload_name": free["paper_compound_name"],
                "release_species_id": spec["release_id"],
                "release_species_name": rel["release_species_name"],
                "transporter_gene": spec["transporter_gene"],
                "free_endpoint": spec["free_endpoint"],
                "free_effect_value": spec["free_effect_value"],
                "free_effect_definition": spec["free_effect_definition"],
                "release_endpoint": spec["release_endpoint"],
                "release_effect_value": spec["release_effect_value"],
                "release_effect_definition": spec["release_effect_definition"],
                "free_mw": free["molecular_weight"],
                "release_mw": rel["molecular_weight"],
                "delta_mw_release_minus_free": delta(free["molecular_weight"], rel["molecular_weight"]),
                "free_xlogp": free["xlogp"],
                "release_xlogp": rel["xlogp"],
                "delta_xlogp_release_minus_free": delta(free["xlogp"], rel["xlogp"]),
                "free_tpsa": free["tpsa"],
                "release_tpsa": rel["tpsa"],
                "delta_tpsa_release_minus_free": delta(free["tpsa"], rel["tpsa"]),
                "directional_comparison": spec["directional_comparison"],
                "endpoint_compatibility": "not_comparable_cross_assay",
                "causal_interpretation": "hypothesis_only_not_causal",
                "free_source_id": "Roth 2026 Table 3",
                "release_source_ids": spec["release_source_ids"],
                "comparison_note": spec["note"] + " Descriptor changes alone do not define P-gp substrate behavior.",
            }
        )
    return rows


def update_issue_log() -> list[dict]:
    rows = read_csv(ROOT / "exclusion_and_conflict_log.csv")
    by_id = {row["issue_id"]: row for row in rows}

    by_id["ISSUE-009"].update(
        {
            "issue_type": "regulatory_epimeric_mixture",
            "current_action": "retain_pubchem_epimer_and_require_epimer_handling_rule",
            "reason": "NCATS GSRS explicitly identifies LYSINE-MCC-DM1 as EPIMERIC with 11 of 12 stereocenters defined; prior PubChem stereochemical divergence is therefore explainable rather than a simple registry error",
            "next_evidence_needed": "Preregister whether modeling uses connectivity-only, enumerated epimers or an explicit mixture representation",
            "status": "partially_resolved",
        }
    )
    by_id["ISSUE-010"].update(
        {
            "entity_or_source": "Lys-PEG4Mal-DM1 / Lys-SPDB-DM4 / cys-VC-MMAE",
            "issue_type": "release_species_structure_resolution_mixed",
            "current_action": "include_Lys-SPDB-DM4_and_cys-VC-MMAE_after_QC_exclude_Lys-PEG4Mal-DM1",
            "reason": "Lys-SPDB-DM4 and cys-VC-MMAE now have registry/PubChem structures cross-linked to primary evidence; Lys-PEG4Mal-DM1 has only a primary chemical scheme without a verified machine structure",
            "next_evidence_needed": "Independent cheminformatics curation of Lys-PEG4Mal-DM1 from the primary scheme",
            "status": "partially_resolved",
        }
    )
    by_id["ISSUE-011"].update(
        {
            "current_action": "retain_structure_and_flag_title_alias_conflict",
            "reason": "PubChem CID 86278355 synonyms explicitly include Cys-McMMAF and CYSTEINYL-MFFA, supporting the small-molecule structure, while whole-ADC aliases still contaminate the preferred title",
            "next_evidence_needed": "Independent registry or reference-standard stereochemistry crosscheck before final model freeze",
            "status": "partially_resolved",
        }
    )

    new_issue = {
        "issue_id": "ISSUE-015",
        "entity_or_source": "Free-payload versus released-species pair table",
        "issue_type": "cross_assay_endpoint_noncomparability",
        "current_action": "allow_directional_hypothesis_only_do_not_compute_cross_assay_fold_change",
        "reason": "Roth cytotoxicity RR, MDCK efflux ratio, inhibitor rescue, comparative retention and regulatory source calls are different endpoints and assay systems",
        "next_evidence_needed": "Matched free-payload/released-species bidirectional transport or accumulation assays in the same cell background",
        "status": "open",
    }
    if "ISSUE-015" not in by_id:
        rows.append(new_issue)
    else:
        by_id["ISSUE-015"].update(new_issue)
    return rows


def validate_dataset(
    raw_log: list[dict],
    registry: list[dict],
    observations: list[dict],
    transport: list[dict],
    pairs: list[dict],
    issues: list[dict],
) -> None:
    checks = [
        (len(registry) == 10, "expected 10 release-species entities"),
        (len(observations) == 17, "expected 17 release observations"),
        (len(transport) == 9, "expected 9 release-species transporter records"),
        (len(pairs) == 7, "expected 7 free-payload/released-species pairs"),
        (len({row["release_species_id"] for row in registry}) == len(registry), "duplicate release_species_id"),
        (len({row["release_observation_id"] for row in observations}) == len(observations), "duplicate release_observation_id"),
        (len({row["transport_record_id"] for row in transport}) == len(transport), "duplicate transport_record_id"),
        (len({row["pair_id"] for row in pairs}) == len(pairs), "duplicate pair_id"),
        (all(row["model_binary_label"] == "" for row in transport), "model_binary_label must remain blank"),
        (all(row["endpoint_compatibility"] == "not_comparable_cross_assay" for row in pairs), "pair endpoint boundary changed"),
        (all(row["causal_interpretation"] == "hypothesis_only_not_causal" for row in pairs), "pair causal boundary changed"),
        (
            sum("unresolved" in row["structure_status"] for row in registry) == 1,
            "expected exactly one machine-unresolved released-species structure",
        ),
        (len({row["issue_id"] for row in issues}) == len(issues), "duplicate issue_id"),
        (
            all(
                row["download_status"]
                in {"downloaded_validated_signature", "reused_existing_validated_signature"}
                for row in raw_log
            ),
            "one or more required raw files were not validated",
        ),
    ]
    failed = [message for passed, message in checks if not passed]
    if failed:
        raise RuntimeError("Third-batch validation failed: " + "; ".join(failed))


def main() -> None:
    raw_log = download_raw()
    failures = [row for row in raw_log if row["download_status"].startswith("failed")]
    required_json = {
        "PubChem_CID131801096_properties.json",
        "PubChem_CID132471209_properties.json",
        "PubChem_CID86278355_synonyms.json",
    }
    missing_required = [name for name in required_json if not (RAW / name).exists()]
    if missing_required:
        raise RuntimeError(f"Required raw JSON missing after download: {missing_required}")

    registry_v2, structure_updates = build_registry_v2()
    observations_v2 = build_observations_v2()
    transport_v2 = build_transport_v2()
    pairs = build_pair_table(registry_v2)
    issues = update_issue_log()
    validate_dataset(raw_log, registry_v2, observations_v2, transport_v2, pairs, issues)

    write_csv(PROCESSED / "third_batch_source_registry.csv", SOURCES)
    write_csv(PROCESSED / "third_batch_raw_download_log.csv", raw_log)
    write_csv(PROCESSED / "release_structure_resolution_updates_v1.csv", structure_updates)
    write_csv(PROCESSED / "adc_release_species_registry_v2.csv", registry_v2)
    write_csv(PROCESSED / "adc_release_observations_long_v2.csv", observations_v2)
    transport_fields = list(transport_v2[0].keys())
    write_csv(PROCESSED / "adc_release_species_transporter_evidence_long_v2.csv", transport_v2, transport_fields)
    write_csv(PROCESSED / "free_payload_release_pair_evidence_v1.csv", pairs)
    write_csv(ROOT / "exclusion_and_conflict_log.csv", issues)

    summary = {
        "generated_at": ACQUIRED_AT,
        "third_batch_sources": len(SOURCES),
        "third_batch_raw_records": len(raw_log),
        "third_batch_raw_failures": len(failures),
        "release_species_entities_v2": len(registry_v2),
        "release_observation_records_v2": len(observations_v2),
        "release_transporter_records_v2": len(transport_v2),
        "free_payload_release_pairs": len(pairs),
        "previously_unresolved_structures": 3,
        "newly_machine_resolved_structures": 2,
        "remaining_machine_unresolved_structures": 1,
        "new_functional_or_regulatory_transporter_records": 4,
        "model_binary_labels_assigned": 0,
        "key_results": [
            "Lys-SPDB-DM4 resolved through NCATS UNII 6MDF4WB915 and PubChem CID 131801096.",
            "cys-VC-MMAE resolved as PubChem CID 132471209 and cross-checked against the PMDA chemical scheme.",
            "Lys-PEG4Mal-DM1 remains excluded because only a primary scheme, not a verified machine structure, is available.",
            "S-methyl-DM1 has three MDR1 functional ratios (8, 7, 6); S-methyl-DM4 has an FDA P-gp substrate call without reproduced numeric output.",
        ],
        "warning": "Pair rows are cross-assay directional hypotheses. No cross-assay numeric fold change or model label was created.",
    }
    (PROCESSED / "third_batch_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
