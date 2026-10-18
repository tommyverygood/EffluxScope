#!/usr/bin/env python3
"""Validate independently retrieved S-Me-DM4 identity without standardization.

Accepted SMILES is copied from PubChem CID 89493732, never synthesized by
deleting a sulfur from the legacy structure. The graph edit below is only a
consistency diagnostic against the independently retrieved DM4 parent.
"""
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors

BASE = Path(__file__).resolve().parent
RAW = BASE / "raw"
ROOT = BASE.parents[2]

def props(name):
    return json.loads((RAW / name).read_text())["PropertyTable"]["Properties"]

def molecular_audit(row):
    mol = Chem.MolFromSmiles(row["SMILES"])
    if mol is None:
        raise ValueError("SMILES parse failed")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    sulfurs = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 16]
    return {
        "cid": row["CID"], "source_smiles": row["SMILES"],
        "rdkit_isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
        "rdkit_inchikey": Chem.MolToInchiKey(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "average_mw": Descriptors.MolWt(mol),
        "exact_mass": rdMolDescriptors.CalcExactMolWt(mol),
        "defined_atom_stereo": sum(1 for _, tag in Chem.FindMolChiralCenters(mol, includeUnassigned=True) if tag != "?"),
        "undefined_atom_stereo": sum(1 for _, tag in Chem.FindMolChiralCenters(mol, includeUnassigned=True) if tag == "?"),
        "double_bond_stereo": [str(b.GetStereo()) for b in mol.GetBonds() if b.GetStereo() != Chem.BondStereo.STEREONONE],
        "sulfur_count": len(sulfurs),
        "sulfur_sulfur_bonds": sum(1 for b in mol.GetBonds() if b.GetBeginAtom().GetAtomicNum() == 16 and b.GetEndAtom().GetAtomicNum() == 16),
        "formal_charge": Chem.GetFormalCharge(mol),
        "source_key_matches_rdkit": Chem.MolToInchiKey(mol) == row["InChIKey"],
        "source_formula_matches_rdkit": rdMolDescriptors.CalcMolFormula(mol) == row["MolecularFormula"],
    }, mol

def main():
    accepted = props("correct_cid_properties.json")[0]
    key_query = props("correct_key_properties.json")[0]
    candidates = props("s_me_dm4.json")
    old = props("old_cid_properties.json")[0]
    parent = props("dm4_parent_properties.json")[0]
    cas_candidates = props("cas_890654_03_2.json")
    audits = []
    accepted_audit, accepted_mol = molecular_audit(accepted)
    for row in candidates + [old, parent]:
        audit, _ = molecular_audit(row)
        audits.append(audit)
    sdf = Chem.SDMolSupplier(str(RAW / "correct_cid_sdf.sdf"), removeHs=True)[0]
    if sdf is None:
        raise ValueError("SDF parse failed")
    sdf_key = Chem.MolToInchiKey(sdf)

    # An independent identity was already retrieved above. De-methylation is a
    # downstream scaffold/stereo conservation check, never the identity source.
    terminal_methyl = []
    for atom in accepted_mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 3 and atom.GetNeighbors()[0].GetAtomicNum() == 16:
            terminal_methyl.append(atom.GetIdx())
    assert len(terminal_methyl) == 1
    demethylated = Chem.RWMol(accepted_mol)
    demethylated.RemoveAtom(terminal_methyl[0])
    dm4_consistency_mol = demethylated.GetMol()
    Chem.SanitizeMol(dm4_consistency_mol)
    dm4_consistency_key = Chem.MolToInchiKey(dm4_consistency_mol)

    checks = {
        "rdkit_version": rdBase.rdkitVersion,
        "validated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "raw_file_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(RAW.iterdir()) if p.is_file() and p.suffix not in [".txt"]},
        "accepted_cid": accepted["CID"],
        "cid_query_equals_key_query": accepted["InChIKey"] == key_query["InChIKey"] and accepted["SMILES"] == key_query["SMILES"],
        "api_smiles_equals_sdf_key": sdf_key == accepted["InChIKey"],
        "accepted_inchikey_matches": accepted_audit["source_key_matches_rdkit"],
        "accepted_formula_matches": accepted_audit["source_formula_matches_rdkit"],
        "single_sulfur_no_sulfur_sulfur_bond": accepted_audit["sulfur_count"] == 1 and accepted_audit["sulfur_sulfur_bonds"] == 0,
        "all_8_atom_centers_defined": accepted_audit["defined_atom_stereo"] == 8 and accepted_audit["undefined_atom_stereo"] == 0,
        "both_double_bonds_E": accepted_audit["double_bond_stereo"] == ["STEREOE", "STEREOE"],
        "parent_core_stereo_conserved": dm4_consistency_key == parent["InChIKey"],
        "demethylated_diagnostic_inchikey": dm4_consistency_key,
        "independent_parent_registry_inchikey": parent["InChIKey"],
        "cas_query_first_cid": cas_candidates[0]["CID"],
        "cas_query_first_inchikey": cas_candidates[0]["InChIKey"],
        "cas_query_first_is_accepted_identity": cas_candidates[0]["InChIKey"] == accepted["InChIKey"],
        "primary2010_figure1_manual_review": "Reviewed original PMC figure pixels: right panel has tertiary C-S-Me, no S-S bond; two drawn backbone double bonds are trans. Registry retains all natural-core stereocenters.",
        "primary2006_experimental_species_evidence": "Erickson 2006 DOI 10.1158/0008-5472.CAN-05-4489, Results/Fig. 3: terminal S-methyl-DM4 peak resistant to DTT/selenol; reported [M+Na]+816.4 and [M+K]+832.5. Original publisher indexed full-text/figure legends read; publisher direct full-text HTTP403, no PDF/SI saved.",
        "primary2022_mass_crosscheck": "Pouzin 2022 DOI 10.1007/s10928-021-09799-0, Table 2: DM4 molecular weight780 and MeDM4 molecular weight794.",
        "gsrs_exact_sme_results": json.loads((RAW / "gsrs_exact_sme_query.json").read_text()).get("total"),
        "gsrs_exact_inchikey_results": json.loads((RAW / "gsrs_sme_inchikey_query.json").read_text()).get("total"),
        "gsrs_independent_sme_identity_confirmed": False,
    }
    required = ["cid_query_equals_key_query", "api_smiles_equals_sdf_key", "accepted_inchikey_matches", "accepted_formula_matches", "single_sulfur_no_sulfur_sulfur_bond", "all_8_atom_centers_defined", "both_double_bonds_E", "parent_core_stereo_conserved"]
    assert all(checks[k] for k in required), checks
    obj = {
        "entity_id": "ADC_REL_004", "name": "S-methyl-DM4 (terminal methyl thioether)",
        "status": "accepted_independent_registry_identity_source_consistent_for_frozen_computational_scoring",
        "exact_isomeric_smiles": accepted["SMILES"],
        "exact_inchikey": accepted["InChIKey"], "inchi": accepted["InChI"],
        "formula": accepted["MolecularFormula"], "registry_average_mw": accepted["MolecularWeight"],
        "registry_exact_mass": accepted["ExactMass"], "pubchem_cid": accepted["CID"],
        "cas_crossreference": "890654-03-2", "cas_status": "ambiguous_registry_stereo_mapping_do_not_select_by_first_CAS_hit",
        "cas_not_present_in_accepted_cid_synonyms": "890654-03-2" not in json.dumps(json.loads((RAW / "correct_cid_synonyms.json").read_text())),
        "primary_structure_source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2954514/#F1",
        "primary_structure_locator": "Lopus et al 2010, Fig.1 right panel; Materials and Methods, Synthesis and Chemistry",
        "registry_structure_source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/89493732",
        "raw_machine_source": "raw/correct_cid_properties.json",
        "approved_for_frozen_scoring": "yes", "independently_reviewed": "yes",
        "review_scope": "Research-agent primary-figure/registry graph audit independent of pair-score builder; no new chemical sample or experimental identity validation.",
        "legacy_record_preserved": True,
    }
    (BASE / "accepted_structure.json").write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")
    (BASE / "identity_checks.json").write_text(json.dumps({"checks": checks, "audits": audits}, indent=2, ensure_ascii=False) + "\n")
    fields = ["entity_id", "exact_isomeric_smiles", "exact_inchikey", "primary_structure_source_url", "registry_structure_source_url", "approved_for_frozen_scoring", "independently_reviewed"]
    with (BASE / "approved_structures_v2.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader(); writer.writerow({f: obj[f] for f in fields})
    cas_wrong = cas_candidates[0]
    sensitivity = [
        {"sensitivity_id": "SMe_DM4_correct_EE_thioether", "role": "accepted_main_case", "structure": accepted},
        {"sensitivity_id": "SMe_DM4_legacy_EE_methyl_disulfide", "role": "wrong_identity_diagnostic_only", "structure": old},
        {"sensitivity_id": "SMe_DM4_CAS_first_ZZ_thioether", "role": "wrong_stereo_diagnostic_only", "structure": cas_wrong},
        {"sensitivity_id": "DM4_parent_EE_thiol", "role": "independent_parent_control", "structure": parent},
    ]
    fields = ["sensitivity_id", "role", "pubchem_cid", "exact_inchikey", "exact_isomeric_smiles", "formula", "molecular_weight"]
    with (BASE / "identity_sensitivity_inputs.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
        for item in sensitivity:
            r = item["structure"]
            w.writerow({"sensitivity_id": item["sensitivity_id"], "role": item["role"], "pubchem_cid": r["CID"], "exact_inchikey": r["InChIKey"], "exact_isomeric_smiles": r["SMILES"], "formula": r["MolecularFormula"], "molecular_weight": r["MolecularWeight"]})
    primary_text = (RAW / "primary2010_pmc.html").read_text()
    assert "Synthesis and Chemistry" in primary_text and "methyl iodide" in primary_text
    checks_path = BASE / "identity_checks.json"
    print(json.dumps({"accepted_cid": accepted["CID"], "accepted_inchikey": accepted["InChIKey"], "formula": accepted["MolecularFormula"], "checks_passed": required, "all_passed": True, "output": str(checks_path)}, indent=2))

if __name__ == "__main__":
    main()
