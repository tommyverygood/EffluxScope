#!/usr/bin/env python3
"""Collect Langthaler et al. 2024 direct human ABCB1 reference panel.

The paper reports conventional 1 h bidirectional Papp/ER measurements in
MDCKII-MDR1 cells transfected with human MDR1.  It is kept as an independent
external source-heldout continuous layer.  The paper's substrate-class
column is retained as source provenance but is never copied into the frozen
binary target.
"""

from __future__ import annotations

import hashlib
import json
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw" / "external_abcb1"
PROCESSED = ROOT / "processed"
SOURCE_ID = "SRC_B23_LANGTHALER2024_HUMAN_ABCB1_REFERENCE_PANEL"
DOI = "10.1186/s12987-023-00495-4"
PMCID = "PMC10809502"
SOURCE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC10809502/"
BIOC_URL = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC10809502/unicode"

STRUCTURES = {
    "altanserin": ("3033677", "C1CN(CCC1C(=O)C2=CC=C(C=C2)F)CCN3C(=O)C4=CC=CC=C4NC3=S"),
    "antipyrine": ("2206", "CC1=CC(=O)N(N1C)C2=CC=CC=C2"),
    "atenolol": ("2249", "CC(C)NCC(COC1=CC=C(C=C1)CC(=O)N)O"),
    "buspirone": ("2477", "C1CCC2(C1)CC(=O)N(C(=O)C2)CCCCN3CCN(CC3)C4=NC=CC=N4"),
    "cimetidine": ("2756", "CC1=C(N=CN1)CSCCNC(=NC)NC#N"),
    "citalopram": ("2771", "CN(C)CCCC1(C2=C(CO1)C=C(C=C2)C#N)C3=CC=C(C=C3)F"),
    "diphenhydramine": ("3100", "CN(C)CCOC(C1=CC=CC=C1)C2=CC=CC=C2"),
    "doxepin": ("3158", "CN(C)CCC=C1C2=CC=CC=C2COC3=CC=CC=C31"),
    "fluoxetine": ("3386", "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F"),
    "gabapentin": ("3446", "C1CCC(CC1)(CC(=O)O)CN"),
    "indomethacin": ("3715", "CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=CC(=C2)OC)CC(=O)O"),
    "metoclopramide": ("4168", "CCN(CC)CCNC(=O)C1=CC(=C(C=C1OC)N)Cl"),
    "n-desmethylclozapine": ("135409468", "C1CN(CCN1)C2=NC3=C(C=CC(=C3)Cl)NC4=CC=CC=C42"),
    "propranolol": ("4946", "CC(C)NCC(COC1=CC=CC2=CC=CC=C21)O"),
    "risperidone": ("5073", "CC1=C(C(=O)N2CCCCC2=N1)CCN3CCC(CC3)C4=NOC5=C4C=CC(=C5)F"),
    "way-100635": ("5684", "COC1=CC=CC=C1N2CCN(CC2)CCN(C3=CC=CC=N3)C(=O)C4CCCCC4"),
}

# Table 1: mean +/- SD, Papp in 10^-6 cm/s; one triplicate test occasion.
VALUES = {
    "altanserin": (3.0, 0.6, 3.4, 0.3, 1.1, "unlikely"),
    "antipyrine": (39.0, 1.0, 48.0, 3.0, 1.2, "unlikely"),
    "atenolol": ("<LLOQ", "", "<LLOQ", "", "NA", "NA"),
    "buspirone": (21.0, 2.0, 20.0, 3.0, 0.9, "unlikely"),
    "cimetidine": ("<LLOQ", "", 2.5, 0.2, "NA", "possible"),
    "citalopram": (4.7, 0.4, 7.7, 0.3, 1.7, "possible"),
    "diphenhydramine": (9.8, 1.3, 16.0, 1.0, 1.6, "possible"),
    "doxepin": (2.0, 0.1, 2.2, 0.1, 1.1, "unlikely"),
    "fluoxetine": (0.4, 0.1, 0.4, 0.1, 1.0, "unlikely"),
    "gabapentin": (2.9, 0.7, 1.6, 0.2, 0.5, "unlikely"),
    "indomethacin": (3.7, 0.5, 2.7, 0.3, 0.7, "unlikely"),
    "metoclopramide": (24.0, 1.0, 49.0, 3.0, 2.0, "likely"),
    "n-desmethylclozapine": (2.1, 0.5, 4.4, 0.3, 2.1, "likely"),
    "propranolol": (4.5, 0.4, 3.7, 0.1, 0.8, "unlikely"),
    "risperidone": (8.5, 1.2, 37.0, 4.0, 4.4, "likely"),
    "way-100635": (13.0, 1.0, 17.0, 3.0, 1.3, "unlikely"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def identity(name: str) -> dict[str, str]:
    cid, smi = STRUCTURES[name]
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(name)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    ik = inchi.MolToInchiKey(mol)
    return {
        "pubchem_cid": cid,
        "standardized_isomeric_smiles": canonical,
        "standardized_inchikey": ik,
        "inchikey_connectivity_block": ik.split("-")[0],
        "source_structure_candidate_smiles": canonical,
        "structure_resolution_status": "resolved_pubchem_name_and_article_table_name",
    }


def observations() -> pd.DataFrame:
    rows = []
    for name, (ab, ab_sd, ba, ba_sd, er, call) in VALUES.items():
        ent = identity(name)
        rows.append({
            "observation_id": f"LANGTHALER2024_{name.upper().replace('-', '_')}_MDR1_0_5UM",
            "entity_name": name,
            "entity_id": f"PUBCHEM_CID_{ent['pubchem_cid']}",
            "transporter_gene": "ABCB1",
            "transporter_species": "human_ABCB1_transfectant",
            "transporter_system": "MDCKII-MDR1",
            "cell_background": "MDCKII",
            "assay_endpoint": "bidirectional_transwell_efflux_ratio",
            "substrate_concentration_value": 0.5,
            "substrate_concentration_unit": "uM",
            "duration_value": 60.0,
            "duration_unit": "min",
            "inhibitor": "none",
            "inhibitor_concentration_value": "",
            "inhibitor_concentration_unit": "",
            "papp_a_to_b_mean_x10^-6_cm_s": ab,
            "papp_a_to_b_sd": ab_sd,
            "papp_b_to_a_mean_x10^-6_cm_s": ba,
            "papp_b_to_a_sd": ba_sd,
            "efflux_ratio_mean": er,
            "efflux_ratio_error_type": "derived from triplicate means; source reports SD for Papp",
            "source_substrate_classification": call,
            "mass_recovery_boundary": ">=85% including filter/cell-associated amount; >=65% excluding it",
            "standardized_isomeric_smiles": ent["standardized_isomeric_smiles"],
            "standardized_inchikey": ent["standardized_inchikey"],
            "inchikey_connectivity_block": ent["inchikey_connectivity_block"],
            "pubchem_cid": ent["pubchem_cid"],
            "source_structure_candidate_smiles": ent["source_structure_candidate_smiles"],
            "structure_resolution_status": ent["structure_resolution_status"],
            "source_id": SOURCE_ID,
            "source_doi": DOI,
            "source_url": SOURCE_URL,
            "source_location": "Table 1; conventional bidirectional 1-h transport experiment",
            "evidence_tier": "A2_primary_direct_numeric_human_ABCB1",
            "source_assigned_call": call,
            "gene_specific_context": "human MDR1-transfected MDCKII; positive digoxin control; triplicate; recovery and monolayer controls",
            "truth_set_eligible": "no_binary_label_continuous_source_heldout_only",
            "model_binary_label": "",
            "claim_boundary": "direct human ABCB1 continuous external source; source call is not transferred into frozen binary task",
            "notes": "Final test concentration 0.5 uM; all compounds tested in triplicate on one occasion; ER >2 likely, >1.5 possible, <1.5 unlikely per source protocol.",
        })
    return pd.DataFrame(rows)


def overlap(out: pd.DataFrame) -> pd.DataFrame:
    daood = pd.read_csv(PROCESSED / "daood_training_records_standardized_v1.csv", dtype=str, keep_default_na=False)
    cand = pd.read_csv(PROCESSED / "transporter_truthset_AB_candidate_v1.csv", dtype=str, keep_default_na=False)
    rows = []
    for name, group in out.groupby("entity_name", sort=True):
        key = group.iloc[0]["standardized_inchikey"]
        dh = daood[(daood.transporter_gene == "ABCB1") & (daood.standardized_inchikey == key)]
        ch = cand[(cand.transporter_gene == "ABCB1") & (cand.standardized_inchikey == key)]
        rows.append({
            "entity_name": name,
            "standardized_inchikey": key,
            "daood_match_count": len(dh),
            "daood_record_ids": " | ".join(dh.training_record_id.tolist()),
            "daood_source_databases": " | ".join(dh.source_database.tolist()),
            "daood_published_labels": " | ".join(dh.published_binary_label.tolist()),
            "candidate_match_count": len(ch),
            "candidate_audit_ids": " | ".join(ch.audit_record_id.tolist()),
            "action": "retain_human_ABCB1_continuous_source_heldout_no_binary_transfer",
            "claim_boundary": "source-heldout direct ER/Papp evidence; existing labels remain frozen",
        })
    return pd.DataFrame(rows)


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    bioc_path = RAW / "langthaler_2024_human_abcb1_reference_panel_bioc.json"
    text_path = RAW / "langthaler_2024_human_abcb1_reference_panel_bioc.txt"
    if not bioc_path.exists():
        bioc_path.write_bytes(urllib.request.urlopen(BIOC_URL, timeout=60).read())
    data = json.loads(bioc_path.read_text())
    text_path.write_text("\n\n".join(f"[{p.get('infons', {}).get('type', '')}] {p.get('text', '')}" for p in data[0]["documents"][0]["passages"] if p.get("text")), encoding="utf-8")
    out = observations()
    out.to_csv(PROCESSED / "langthaler_2024_human_abcb1_reference_panel_observations_v1.csv", index=False, encoding="utf-8-sig")
    audit = overlap(out)
    audit.to_csv(PROCESSED / "langthaler_2024_human_abcb1_reference_panel_identity_overlap_audit_v1.csv", index=False, encoding="utf-8-sig")
    registry = pd.DataFrame([{
        "source_id": SOURCE_ID,
        "doi": DOI,
        "pmcid": PMCID,
        "source_url": SOURCE_URL,
        "bioc_url": BIOC_URL,
        "raw_relative_path": "raw/external_abcb1/langthaler_2024_human_abcb1_reference_panel_bioc.json",
        "raw_sha256": sha256(bioc_path),
        "parsed_observations": len(out),
        "unique_entities": out.entity_name.nunique(),
        "numeric_rows": int(out.efflux_ratio_mean.astype(str).str.match(r"^\d").sum()),
        "model_binary_label_nonblank": int(out.model_binary_label.astype(str).str.strip().ne("").sum()),
        "binary_label_transfer": "none",
        "allowed_use": "direct human ABCB1 source-heldout continuous benchmark; no binary label transfer",
    }])
    registry.to_csv(PROCESSED / "langthaler_2024_human_abcb1_reference_panel_source_registry_v1.csv", index=False, encoding="utf-8-sig")
    summary = {
        "batch": "langthaler_2024_human_abcb1_reference_panel_v1",
        "acquired_at": date.today().isoformat(),
        "raw_sha256": sha256(bioc_path),
        "observations": len(out),
        "entities": sorted(out.entity_name.unique().tolist()),
        "numeric_rows": int(out.efflux_ratio_mean.astype(str).str.match(r"^\d").sum()),
        "source_class_counts": out.source_substrate_classification.value_counts().to_dict(),
        "model_binary_label_nonblank": int(out.model_binary_label.astype(str).str.strip().ne("").sum()),
        "daood_exact_overlap_count": int(audit.daood_match_count.sum()),
        "candidate_exact_overlap_count": int(audit.candidate_match_count.sum()),
        "claim_boundary": "direct human ABCB1 ER/Papp source-heldout continuous evidence; no binary label transfer",
    }
    (PROCESSED / "langthaler_2024_human_abcb1_reference_panel_summary_v1.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
