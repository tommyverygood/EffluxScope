#!/usr/bin/env python3
"""Standardize ADC structures, audit similarity leakage, and fit evidence-tiered baselines.

This pipeline never overwrites source SMILES and never assigns model labels to
Roth continuous resistance ratios or ADC released-species evidence. Models are
development baselines with scaffold-aware cross-validation and explicit
similarity-support flags, not final deployable predictors.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, pstdev

import joblib
import numpy as np
import scipy
import sklearn
from rdkit import Chem, DataStructs, RDLogger, rdBase
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


DATA_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = DATA_ROOT.parent
PROCESSED = DATA_ROOT / "processed"
OUTPUT_DIR = DATA_ROOT / "outputs" / "20260728_standardization_leakage_baseline"
MODELS_DIR = PROJECT_ROOT / "模型" / "20260728_baseline"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
MODELS_DIR.mkdir(parents=True, exist_ok=True)

RUN_DATE = "2026-07-28"
RANDOM_SEED = 20260728
MORGAN_RADIUS = 2
MORGAN_BITS = 2048
NEAR_DUPLICATE_THRESHOLD = 0.80
SCAFFOLD_FAMILY_THRESHOLD = 0.50
COMPLEX_SCAFFOLD_MIN_HEAVY_ATOMS = 10
CV_SPLITS = 5
INTERNAL_MCC_GATE = 0.50
INTERNAL_BALANCED_ACCURACY_GATE = 0.75

FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(
    radius=MORGAN_RADIUS,
    fpSize=MORGAN_BITS,
    includeChirality=False,
)
TAUTOMER_ENUMERATOR = rdMolStandardize.TautomerEnumerator()
UNCHARGER = rdMolStandardize.Uncharger(canonicalOrder=True)
# The enumerator status is captured explicitly below, so repeated limit warnings
# are suppressed without hiding parse failures or standardization exceptions.
RDLogger.DisableLog("rdApp.warning")

DESCRIPTOR_FIELDS = [
    "std_molecular_weight",
    "std_xlogp",
    "std_tpsa",
    "std_hbond_donor_count",
    "std_hbond_acceptor_count",
    "std_rotatable_bond_count",
    "std_formal_charge",
    "std_heavy_atom_count",
]

CONFLICT_FIELDS = [
    "training_record_id",
    "pid",
    "preferred_name",
    "transporter_gene",
    "issue_type",
    "detail",
    "model_action",
]


STANDARDIZATION_RULES = {
    "version": "v1",
    "run_date": RUN_DATE,
    "input_policy": "use source isomeric SMILES when present; otherwise source connectivity/published SMILES",
    "steps": [
        "RDKit MolFromSmiles with sanitization",
        "rdMolStandardize.Cleanup",
        "rdMolStandardize.FragmentParent with organic preference inherited from defaults",
        "rdMolStandardize.Uncharger canonical order; permanent charges may remain",
        "TautomerEnumerator canonical tautomer",
        "canonical isomeric and connectivity SMILES plus Standard InChIKey",
    ],
    "raw_structure_policy": "source SMILES are immutable and retained in source fields",
    "epimeric_policy": "registered epimeric mixtures use connectivity for audit/model input and retain an explicit representative-structure flag",
    "unresolved_policy": "blank or unparsable structures are excluded from structure models without imputation",
    "tautomer_qc_policy": (
        "enumerate then select the RDKit canonical tautomer; retain enumeration status/count and flag "
        "MaxTautomersReached, MaxTransformsReached, or Canceled as best-effort standardization QC"
    ),
    "tautomer_enumerator_limits": {
        "max_tautomers": TAUTOMER_ENUMERATOR.GetMaxTautomers(),
        "max_transforms": TAUTOMER_ENUMERATOR.GetMaxTransforms(),
    },
    "fingerprint": {
        "type": "Morgan bit fingerprint",
        "radius": MORGAN_RADIUS,
        "n_bits": MORGAN_BITS,
        "include_chirality": False,
        "similarity": "Tanimoto",
    },
    "leakage_rules": {
        "exact": "same standardized connectivity SMILES or same InChIKey connectivity block",
        "near_duplicate": f"Morgan Tanimoto >= {NEAR_DUPLICATE_THRESHOLD:.2f}",
        "complex_scaffold_family": (
            f"same exact Murcko scaffold with >= {COMPLEX_SCAFFOLD_MIN_HEAVY_ATOMS} heavy atoms "
            f"and Morgan Tanimoto >= {SCAFFOLD_FAMILY_THRESHOLD:.2f}"
        ),
        "model_action": "exclude exact, near-duplicate, and complex-scaffold-family training rows from external evaluation models",
        "common_scaffold_action": "flag only; do not automatically exclude simple/common scaffolds",
    },
    "model_rules": {
        "targets": "Daood published substrate labels only, ABCB1 and ABCG2 fitted separately",
        "duplicate_rule": "collapse identical standardized structure-transporter records when labels agree; exclude conflicts",
        "cross_validation": "5-fold StratifiedGroupKFold grouped by exact Murcko scaffold; acyclic structures grouped by structure ID",
        "baselines": ["class-prior dummy", "8-descriptor logistic regression", "Morgan logistic regression"],
        "Roth_role": "external continuous stress-test; evaluate Spearman association with resistance ratio without binarization",
        "released_species_role": "external evidence-tier audit using source labels only; no model label assignment",
        "claim_boundary": "no final predictor claim without independent matched released-species validation",
    },
}


EVIDENCE_TIER_MAP = {
    "A_direct_numeric_inhibitor_rescue": (
        "A1_direct_numeric_rescue",
        "bidirectional_transport",
        "external_binary_audit",
        "aggregate one entity-transporter record; retain numeric inhibitor rescue",
    ),
    "A_direct_numeric_no_rescue_reported": (
        "A2_direct_numeric_no_rescue",
        "bidirectional_transport",
        "external_binary_audit",
        "aggregate one entity-transporter record; retain absence of inhibitor control",
    ),
    "A_direct_qualitative_inhibitor_control": (
        "A3_direct_qualitative",
        "bidirectional_transport_qualitative",
        "external_binary_audit_with_numeric_hold",
        "aggregate one entity-transporter record; do not invent Figure values",
    ),
    "B_functional_resistance_or_inhibitor_rescue": (
        "B1_functional_rescue",
        "cytotoxicity_resistance_or_rescue",
        "external_functional_audit",
        "aggregate replicate/cell-system records at entity-transporter level",
    ),
    "B_regulatory_call_no_numeric": (
        "B2_regulatory_call",
        "regulatory_source_call",
        "external_call_only_audit",
        "one entity-transporter source call; numeric endpoint remains missing",
    ),
    "C_comparative_cellular_retention": (
        "C_comparative_only",
        "comparative_retention",
        "exclude_from_binary_evaluation",
        "retain as mechanism context only",
    ),
}


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def short_hash(text: str, length: int = 12) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:length]


def safe_float(value) -> float | None:
    try:
        if value is None or value == "":
            return None
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def connectivity_block(inchikey: str) -> str:
    return inchikey.split("-")[0] if inchikey else ""


def formal_charge(mol: Chem.Mol) -> int:
    return int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))


def standardize_smiles(source_smiles: str) -> dict:
    base = {
        "source_smiles": source_smiles or "",
        "standardization_status": "unresolved_blank_structure",
        "standardized_isomeric_smiles": "",
        "standardized_connectivity_smiles": "",
        "standardized_inchikey": "",
        "inchikey_connectivity_block": "",
        "murcko_scaffold_smiles": "",
        "generic_murcko_scaffold_smiles": "",
        "scaffold_heavy_atom_count": "",
        "fragment_count_source": "",
        "fragment_removed_flag": "",
        "source_formal_charge": "",
        "standardized_formal_charge": "",
        "tautomer_changed_flag": "",
        "tautomer_enumeration_status": "not_run",
        "tautomer_enumerated_count": "",
        "tautomer_limit_flag": "",
        "defined_stereocenter_count": "",
        "unassigned_stereocenter_count": "",
        "standardization_actions": "",
        **{field: "" for field in DESCRIPTOR_FIELDS},
    }
    if not source_smiles or not source_smiles.strip():
        return base

    try:
        source_mol = Chem.MolFromSmiles(source_smiles, sanitize=True)
        if source_mol is None:
            base["standardization_status"] = "parse_failed"
            return base
        source_fragments = len(Chem.GetMolFrags(source_mol))
        source_charge = formal_charge(source_mol)
        cleaned = rdMolStandardize.Cleanup(source_mol)
        parent = rdMolStandardize.FragmentParent(cleaned, skipStandardize=True)
        uncharged = UNCHARGER.uncharge(parent)
        pre_tautomer_connectivity = Chem.MolToSmiles(uncharged, canonical=True, isomericSmiles=False)
        tautomer_result = TAUTOMER_ENUMERATOR.Enumerate(uncharged)
        tautomer_status = getattr(tautomer_result.status, "name", str(tautomer_result.status))
        tautomer_count = len(tautomer_result)
        tautomer_limit_flag = "no" if tautomer_status == "Completed" else "yes"
        standardized = TAUTOMER_ENUMERATOR.PickCanonical(tautomer_result)
        Chem.AssignStereochemistry(standardized, cleanIt=True, force=True)

        iso = Chem.MolToSmiles(standardized, canonical=True, isomericSmiles=True)
        conn = Chem.MolToSmiles(standardized, canonical=True, isomericSmiles=False)
        try:
            inchikey = Chem.MolToInchiKey(standardized)
        except Exception:
            inchikey = ""

        scaffold_mol = MurckoScaffold.GetScaffoldForMol(standardized)
        scaffold = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False) if scaffold_mol.GetNumAtoms() else ""
        if scaffold_mol.GetNumAtoms():
            generic_scaffold = Chem.MolToSmiles(
                MurckoScaffold.MakeScaffoldGeneric(scaffold_mol),
                canonical=True,
                isomericSmiles=False,
            )
        else:
            generic_scaffold = ""

        chiral = Chem.FindMolChiralCenters(standardized, includeUnassigned=True, useLegacyImplementation=False)
        defined = sum(label != "?" for _, label in chiral)
        unassigned = sum(label == "?" for _, label in chiral)
        std_charge = formal_charge(standardized)
        actions = ["cleanup", "fragment_parent", "uncharge", "canonical_tautomer"]
        if source_fragments > 1:
            actions.append("fragment_removed")
        if source_charge != std_charge:
            actions.append("charge_changed")
        if pre_tautomer_connectivity != conn:
            actions.append("tautomer_changed")
        if tautomer_limit_flag == "yes":
            actions.append(f"tautomer_enumeration_{tautomer_status}")

        base.update(
            {
                "standardization_status": "standardized",
                "standardized_isomeric_smiles": iso,
                "standardized_connectivity_smiles": conn,
                "standardized_inchikey": inchikey,
                "inchikey_connectivity_block": connectivity_block(inchikey),
                "murcko_scaffold_smiles": scaffold,
                "generic_murcko_scaffold_smiles": generic_scaffold,
                "scaffold_heavy_atom_count": int(scaffold_mol.GetNumHeavyAtoms()),
                "fragment_count_source": source_fragments,
                "fragment_removed_flag": "yes" if source_fragments > 1 else "no",
                "source_formal_charge": source_charge,
                "standardized_formal_charge": std_charge,
                "tautomer_changed_flag": "yes" if pre_tautomer_connectivity != conn else "no",
                "tautomer_enumeration_status": tautomer_status,
                "tautomer_enumerated_count": tautomer_count,
                "tautomer_limit_flag": tautomer_limit_flag,
                "defined_stereocenter_count": defined,
                "unassigned_stereocenter_count": unassigned,
                "standardization_actions": " | ".join(actions),
                "std_molecular_weight": round(float(Descriptors.MolWt(standardized)), 4),
                "std_xlogp": round(float(Crippen.MolLogP(standardized)), 4),
                "std_tpsa": round(float(rdMolDescriptors.CalcTPSA(standardized)), 4),
                "std_hbond_donor_count": int(Lipinski.NumHDonors(standardized)),
                "std_hbond_acceptor_count": int(Lipinski.NumHAcceptors(standardized)),
                "std_rotatable_bond_count": int(Lipinski.NumRotatableBonds(standardized)),
                "std_formal_charge": std_charge,
                "std_heavy_atom_count": int(standardized.GetNumHeavyAtoms()),
            }
        )
        return base
    except Exception as exc:
        base["standardization_status"] = "standardization_failed"
        base["standardization_actions"] = f"{type(exc).__name__}: {exc}"
        return base


def mol_from_std(row: dict) -> Chem.Mol | None:
    smiles = row.get("standardized_isomeric_smiles", "")
    return Chem.MolFromSmiles(smiles) if smiles else None


def fingerprint(row: dict):
    mol = mol_from_std(row)
    return FP_GENERATOR.GetFingerprint(mol) if mol is not None else None


def fp_numpy(fp) -> np.ndarray:
    array = np.zeros((MORGAN_BITS,), dtype=np.float64)
    DataStructs.ConvertToNumpyArray(fp, array)
    return array


def descriptor_numpy(row: dict) -> np.ndarray:
    return np.array([float(row[field]) for field in DESCRIPTOR_FIELDS], dtype=np.float64)


def make_standardized_registries() -> tuple[list[dict], list[dict], list[dict]]:
    free_source = read_csv(PROCESSED / "compound_registry_v1.csv")
    release_source = read_csv(PROCESSED / "adc_release_species_registry_v2.csv")
    training_source = read_csv(PROCESSED / "daood_2025_pgp_bcrp_substrate_training_long.csv")

    free_rows = []
    for source in free_source:
        structure = source.get("isomeric_smiles") or source.get("connectivity_smiles") or ""
        std = standardize_smiles(structure)
        free_rows.append(
            {
                "reference_object_type": "free_payload",
                "reference_id": source["compound_id"],
                "reference_name": source["paper_compound_name"],
                "source_structure_status": source["structure_status"],
                "source_pubchem_cid": source["pubchem_cid"],
                "source_inchikey": source["inchikey"],
                "source_structure_url": source["structure_source_url"],
                "special_structure_flag": (
                    "salt_form_review" if source["structure_status"] == "salt_form_review" else ""
                ),
                **std,
            }
        )

    release_rows = []
    for source in release_source:
        structure = source.get("isomeric_smiles") or source.get("connectivity_smiles") or ""
        std = standardize_smiles(structure)
        special = ""
        if source["structure_status"] == "regulatory_epimeric_mixture":
            special = "epimeric_mixture_representative_connectivity_only"
        elif "unresolved" in source["structure_status"]:
            special = "machine_structure_unresolved"
        elif "candidate" in source["structure_status"] or "conflict" in source["structure_status"]:
            special = "source_crosscheck_flag_retained"
        release_rows.append(
            {
                "reference_object_type": "release_species",
                "reference_id": source["release_species_id"],
                "reference_name": source["release_species_name"],
                "payload_core": source["payload_core"],
                "residual_group": source["residual_group"],
                "source_structure_status": source["structure_status"],
                "source_pubchem_cid": source["pubchem_cid"],
                "source_inchikey": source["inchikey"],
                "source_structure_url": source["structure_source_url"],
                "special_structure_flag": special,
                **std,
            }
        )

    training_rows = []
    for source in training_source:
        std = standardize_smiles(source["smiles_as_published"])
        structure_id = (
            f"DAOOD_STD_{short_hash(std['standardized_connectivity_smiles'])}"
            if std["standardization_status"] == "standardized"
            else f"DAOOD_UNRESOLVED_{source['training_record_id']}"
        )
        training_rows.append(
            {
                "training_record_id": source["training_record_id"],
                "training_structure_id": structure_id,
                "pid": source["pid"],
                "source_database": source["source_database"],
                "preferred_name": source["preferred_name"],
                "transporter_gene": source["transporter_gene"],
                "transporter_protein": source["transporter_protein"],
                "published_binary_label": source["published_binary_label"],
                "source_processing_status": source["processing_status"],
                "source_doi": source["source_doi"],
                "source_url": source["source_url"],
                **std,
            }
        )
    return free_rows, release_rows, training_rows


def collapse_training(training_rows: list[dict]) -> tuple[list[dict], list[dict]]:
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    unresolved = []
    for row in training_rows:
        if row["standardization_status"] != "standardized":
            unresolved.append(
                {
                    "training_record_id": row["training_record_id"],
                    "pid": row["pid"],
                    "preferred_name": row["preferred_name"],
                    "transporter_gene": row["transporter_gene"],
                    "issue_type": "structure_not_standardized",
                    "detail": row["standardization_actions"],
                    "model_action": "exclude",
                }
            )
            continue
        grouped[(row["transporter_gene"], row["standardized_connectivity_smiles"])].append(row)

    consensus = []
    conflicts = list(unresolved)
    for (gene, smiles), members in sorted(grouped.items()):
        labels = sorted({member["published_binary_label"] for member in members})
        first = members[0]
        label_status = "consensus_agree" if len(labels) == 1 else "conflicting_labels"
        consensus_label = labels[0] if len(labels) == 1 else ""
        if len(labels) > 1:
            conflicts.append(
                {
                    "training_record_id": " | ".join(member["training_record_id"] for member in members),
                    "pid": " | ".join(sorted({member["pid"] for member in members})),
                    "preferred_name": " | ".join(sorted({member["preferred_name"] for member in members})),
                    "transporter_gene": gene,
                    "issue_type": "conflicting_labels_after_standardization",
                    "detail": f"published labels: {' | '.join(labels)}",
                    "model_action": "exclude",
                }
            )
        consensus.append(
            {
                "consensus_record_id": f"DAOOD_CONS_{gene}_{short_hash(smiles)}",
                "training_structure_id": first["training_structure_id"],
                "transporter_gene": gene,
                "transporter_protein": first["transporter_protein"],
                "consensus_binary_label": consensus_label,
                "label_status": label_status,
                "source_record_count": len(members),
                "source_training_record_ids": " | ".join(member["training_record_id"] for member in members),
                "source_pids": " | ".join(sorted({member["pid"] for member in members})),
                "source_databases": " | ".join(sorted({member["source_database"] for member in members})),
                "preferred_names": " | ".join(sorted({member["preferred_name"] for member in members})),
                **{
                    field: first[field]
                    for field in [
                        "standardization_status",
                        "standardized_isomeric_smiles",
                        "standardized_connectivity_smiles",
                        "standardized_inchikey",
                        "inchikey_connectivity_block",
                        "murcko_scaffold_smiles",
                        "generic_murcko_scaffold_smiles",
                        "scaffold_heavy_atom_count",
                        "tautomer_enumeration_status",
                        "tautomer_enumerated_count",
                        "tautomer_limit_flag",
                        *DESCRIPTOR_FIELDS,
                    ]
                },
            }
        )
    return consensus, conflicts


def audit_similarity(
    references: list[dict], consensus: list[dict]
) -> tuple[list[dict], list[dict], list[dict]]:
    valid_refs = [row for row in references if row["standardization_status"] == "standardized"]
    ref_fps = [fingerprint(row) for row in valid_refs]
    ref_by_smiles: dict[str, list[dict]] = defaultdict(list)
    ref_by_conn_block: dict[str, list[dict]] = defaultdict(list)
    for ref in valid_refs:
        ref_by_smiles[ref["standardized_connectivity_smiles"]].append(ref)
        if ref["inchikey_connectivity_block"]:
            ref_by_conn_block[ref["inchikey_connectivity_block"]].append(ref)

    unique_training: dict[str, dict] = {}
    genes_by_smiles: dict[str, set[str]] = defaultdict(set)
    names_by_smiles: dict[str, set[str]] = defaultdict(set)
    ids_by_smiles: dict[str, set[str]] = defaultdict(set)
    for row in consensus:
        smiles = row["standardized_connectivity_smiles"]
        unique_training.setdefault(smiles, row)
        genes_by_smiles[smiles].add(row["transporter_gene"])
        names_by_smiles[smiles].update(row["preferred_names"].split(" | "))
        ids_by_smiles[smiles].add(row["training_structure_id"])

    unique_rows = list(unique_training.values())
    unique_fps = [fingerprint(row) for row in unique_rows]

    flagged_pairs = []
    reference_summary = []
    for ref, ref_fp in zip(valid_refs, ref_fps):
        similarities = DataStructs.BulkTanimotoSimilarity(ref_fp, unique_fps)
        order = sorted(range(len(similarities)), key=lambda idx: (-similarities[idx], unique_rows[idx]["training_structure_id"]))
        exact_count = 0
        near_count = 0
        moderate_count = 0
        same_scaffold_count = 0
        complex_scaffold_count = 0
        for rank, idx in enumerate(order, start=1):
            train = unique_rows[idx]
            similarity = float(similarities[idx])
            exact_smiles = train["standardized_connectivity_smiles"] == ref["standardized_connectivity_smiles"]
            exact_block = bool(
                train["inchikey_connectivity_block"]
                and train["inchikey_connectivity_block"] == ref["inchikey_connectivity_block"]
            )
            exact = exact_smiles or exact_block
            same_scaffold = bool(
                ref["murcko_scaffold_smiles"]
                and train["murcko_scaffold_smiles"] == ref["murcko_scaffold_smiles"]
            )
            scaffold_heavy = int(ref["scaffold_heavy_atom_count"] or 0)
            complex_scaffold = (
                same_scaffold
                and scaffold_heavy >= COMPLEX_SCAFFOLD_MIN_HEAVY_ATOMS
                and similarity >= SCAFFOLD_FAMILY_THRESHOLD
            )
            exact_count += int(exact)
            near_count += int(similarity >= NEAR_DUPLICATE_THRESHOLD)
            moderate_count += int(similarity >= SCAFFOLD_FAMILY_THRESHOLD)
            same_scaffold_count += int(same_scaffold)
            complex_scaffold_count += int(complex_scaffold)
            if rank <= 5 or exact or similarity >= SCAFFOLD_FAMILY_THRESHOLD or complex_scaffold:
                flagged_pairs.append(
                    {
                        "reference_object_type": ref["reference_object_type"],
                        "reference_id": ref["reference_id"],
                        "reference_name": ref["reference_name"],
                        "training_structure_id": train["training_structure_id"],
                        "training_preferred_names": " | ".join(sorted(names_by_smiles[train["standardized_connectivity_smiles"]])),
                        "training_transporter_genes": " | ".join(sorted(genes_by_smiles[train["standardized_connectivity_smiles"]])),
                        "similarity_rank_for_reference": rank,
                        "morgan_tanimoto": round(similarity, 6),
                        "exact_standardized_connectivity": "yes" if exact else "no",
                        "same_murcko_scaffold": "yes" if same_scaffold else "no",
                        "complex_scaffold_family_flag": "yes" if complex_scaffold else "no",
                        "reference_murcko_scaffold": ref["murcko_scaffold_smiles"],
                        "training_murcko_scaffold": train["murcko_scaffold_smiles"],
                        "audit_action": (
                            "exclude_exact"
                            if exact
                            else "exclude_near_duplicate"
                            if similarity >= NEAR_DUPLICATE_THRESHOLD
                            else "exclude_complex_scaffold_family"
                            if complex_scaffold
                            else "flag_only"
                        ),
                    }
                )
        nearest = unique_rows[order[0]] if order else None
        reference_summary.append(
            {
                "reference_object_type": ref["reference_object_type"],
                "reference_id": ref["reference_id"],
                "reference_name": ref["reference_name"],
                "standardization_status": ref["standardization_status"],
                "special_structure_flag": ref.get("special_structure_flag", ""),
                "nearest_training_structure_id": nearest["training_structure_id"] if nearest else "",
                "nearest_training_names": (
                    " | ".join(sorted(names_by_smiles[nearest["standardized_connectivity_smiles"]])) if nearest else ""
                ),
                "max_morgan_tanimoto": round(float(similarities[order[0]]), 6) if order else "",
                "exact_match_training_structure_count": exact_count,
                "near_duplicate_training_structure_count_ge_080": near_count,
                "moderate_similarity_count_ge_050": moderate_count,
                "same_murcko_scaffold_count": same_scaffold_count,
                "complex_scaffold_family_count": complex_scaffold_count,
                "external_evaluation_leakage_status": (
                    "exact_overlap_present"
                    if exact_count
                    else "near_duplicate_present"
                    if near_count
                    else "complex_scaffold_family_present"
                    if complex_scaffold_count
                    else "no_exclusion_level_overlap"
                ),
            }
        )

    decisions = []
    for row in consensus:
        fp = fingerprint(row)
        similarities = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
        max_idx = int(np.argmax(similarities))
        max_ref = valid_refs[max_idx]
        max_similarity = float(similarities[max_idx])
        exact_refs = {
            ref["reference_id"]
            for ref in ref_by_smiles.get(row["standardized_connectivity_smiles"], [])
        }
        if row["inchikey_connectivity_block"]:
            exact_refs.update(
                ref["reference_id"]
                for ref in ref_by_conn_block.get(row["inchikey_connectivity_block"], [])
            )
        same_scaffold_refs = [
            ref
            for ref in valid_refs
            if row["murcko_scaffold_smiles"]
            and ref["murcko_scaffold_smiles"] == row["murcko_scaffold_smiles"]
        ]
        complex_scaffold_refs = [
            ref
            for ref in same_scaffold_refs
            if int(ref["scaffold_heavy_atom_count"] or 0) >= COMPLEX_SCAFFOLD_MIN_HEAVY_ATOMS
            and DataStructs.TanimotoSimilarity(fp, fingerprint(ref)) >= SCAFFOLD_FAMILY_THRESHOLD
        ]
        if row["label_status"] == "conflicting_labels":
            action = "exclude_label_conflict"
        elif exact_refs:
            action = "exclude_external_exact_match"
        elif max_similarity >= NEAR_DUPLICATE_THRESHOLD:
            action = "exclude_external_near_duplicate"
        elif complex_scaffold_refs:
            action = "exclude_external_complex_scaffold_family"
        else:
            action = "include_leakage_clean_training"
        decisions.append(
            {
                **row,
                "nearest_external_reference_type": max_ref["reference_object_type"],
                "nearest_external_reference_id": max_ref["reference_id"],
                "nearest_external_reference_name": max_ref["reference_name"],
                "max_external_morgan_tanimoto": round(max_similarity, 6),
                "exact_external_reference_ids": " | ".join(sorted(exact_refs)),
                "same_murcko_scaffold_reference_ids": " | ".join(sorted(ref["reference_id"] for ref in same_scaffold_refs)),
                "complex_scaffold_family_reference_ids": " | ".join(sorted(ref["reference_id"] for ref in complex_scaffold_refs)),
                "leakage_model_action": action,
            }
        )
    return reference_summary, flagged_pairs, decisions


def compute_metrics(y_true: np.ndarray, probabilities: np.ndarray) -> dict:
    predictions = (probabilities >= 0.5).astype(int)
    both_classes = len(set(y_true.tolist())) == 2
    return {
        "roc_auc": float(roc_auc_score(y_true, probabilities)) if both_classes else math.nan,
        "average_precision": float(average_precision_score(y_true, probabilities)) if both_classes else math.nan,
        "balanced_accuracy": float(balanced_accuracy_score(y_true, predictions)),
        "mcc": float(matthews_corrcoef(y_true, predictions)),
        "brier_score": float(brier_score_loss(y_true, probabilities)),
    }


def make_estimators() -> dict:
    return {
        "dummy_prior": DummyClassifier(strategy="prior"),
        "physchem_logistic": Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        C=1.0,
                        class_weight="balanced",
                        max_iter=3000,
                        solver="liblinear",
                        random_state=RANDOM_SEED,
                    ),
                ),
            ]
        ),
        "morgan_logistic": LogisticRegression(
            C=1.0,
            class_weight="balanced",
            max_iter=3000,
            solver="liblinear",
            random_state=RANDOM_SEED,
        ),
    }


def features_for(rows: list[dict], model_name: str) -> np.ndarray:
    if model_name == "physchem_logistic":
        return np.vstack([descriptor_numpy(row) for row in rows])
    if model_name == "morgan_logistic":
        return np.vstack([fp_numpy(fingerprint(row)) for row in rows])
    return np.zeros((len(rows), 1), dtype=np.float64)


def similarity_support_thresholds(rows: list[dict]) -> dict:
    fps = [fingerprint(row) for row in rows]
    maxima = []
    for idx, fp in enumerate(fps):
        others = fps[:idx] + fps[idx + 1 :]
        maxima.append(max(DataStructs.BulkTanimotoSimilarity(fp, others)) if others else 0.0)
    return {
        "q10": float(np.quantile(maxima, 0.10)),
        "q25": float(np.quantile(maxima, 0.25)),
        "median": float(np.quantile(maxima, 0.50)),
        "mean": float(np.mean(maxima)),
    }


def support_for(row: dict, training_rows: list[dict], thresholds: dict) -> tuple[float, str, str]:
    fp = fingerprint(row)
    train_fps = [fingerprint(item) for item in training_rows]
    maximum = max(DataStructs.BulkTanimotoSimilarity(fp, train_fps)) if train_fps else 0.0
    if maximum < thresholds["q10"]:
        band = "outside_training_similarity_support"
        abstain = "yes"
    elif maximum < thresholds["q25"]:
        band = "low_similarity_support"
        abstain = "yes_caution"
    else:
        band = "within_empirical_similarity_support"
        abstain = "no_similarity_abstention"
    return float(maximum), band, abstain


def fit_baselines(decisions: list[dict]):
    cv_metrics = []
    oof_predictions = []
    support_rows = []
    fitted_models = {}
    training_by_gene = {}
    thresholds_by_gene = {}

    for gene in ["ABCB1", "ABCG2"]:
        rows = [
            row
            for row in decisions
            if row["transporter_gene"] == gene
            and row["leakage_model_action"] == "include_leakage_clean_training"
            and row["consensus_binary_label"] in {"0", "1"}
        ]
        training_by_gene[gene] = rows
        thresholds = similarity_support_thresholds(rows)
        thresholds_by_gene[gene] = thresholds
        support_rows.append(
            {
                "transporter_gene": gene,
                "training_structure_count": len(rows),
                "self_nearest_similarity_q10": round(thresholds["q10"], 6),
                "self_nearest_similarity_q25": round(thresholds["q25"], 6),
                "self_nearest_similarity_median": round(thresholds["median"], 6),
                "self_nearest_similarity_mean": round(thresholds["mean"], 6),
                "abstain_rule": "max similarity below q10",
                "caution_rule": "q10 <= max similarity < q25",
            }
        )

        y = np.array([int(row["consensus_binary_label"]) for row in rows], dtype=int)
        groups = np.array(
            [
                row["murcko_scaffold_smiles"]
                if row["murcko_scaffold_smiles"]
                else f"ACYCLIC::{row['training_structure_id']}"
                for row in rows
            ]
        )
        splitter = StratifiedGroupKFold(n_splits=CV_SPLITS, shuffle=True, random_state=RANDOM_SEED)
        split_indices = list(splitter.split(np.zeros(len(rows)), y, groups))

        for model_name, estimator in make_estimators().items():
            fold_metric_rows = []
            for fold, (train_idx, test_idx) in enumerate(split_indices, start=1):
                train_rows = [rows[idx] for idx in train_idx]
                test_rows = [rows[idx] for idx in test_idx]
                x_train = features_for(train_rows, model_name)
                x_test = features_for(test_rows, model_name)
                y_train = y[train_idx]
                y_test = y[test_idx]
                model = sklearn.base.clone(estimator)
                model.fit(x_train, y_train)
                probabilities = model.predict_proba(x_test)[:, 1]
                metrics = compute_metrics(y_test, probabilities)
                metric_row = {
                    "transporter_gene": gene,
                    "model_name": model_name,
                    "fold": fold,
                    "n_train": len(train_idx),
                    "n_test": len(test_idx),
                    "test_positive_count": int(y_test.sum()),
                    "test_negative_count": int(len(y_test) - y_test.sum()),
                    **{key: round(value, 6) if math.isfinite(value) else "" for key, value in metrics.items()},
                    "split_strategy": "StratifiedGroupKFold_exact_Murcko_scaffold",
                }
                cv_metrics.append(metric_row)
                fold_metric_rows.append(metric_row)
                for local_idx, probability in zip(test_idx, probabilities):
                    row = rows[int(local_idx)]
                    oof_predictions.append(
                        {
                            "consensus_record_id": row["consensus_record_id"],
                            "training_structure_id": row["training_structure_id"],
                            "transporter_gene": gene,
                            "model_name": model_name,
                            "fold": fold,
                            "published_binary_label": row["consensus_binary_label"],
                            "predicted_probability": round(float(probability), 8),
                            "predicted_label_at_0_5": int(probability >= 0.5),
                            "murcko_scaffold_smiles": row["murcko_scaffold_smiles"],
                            "preferred_names": row["preferred_names"],
                        }
                    )

            for summary_name, function in [("mean", mean), ("sd", pstdev)]:
                summary = {
                    "transporter_gene": gene,
                    "model_name": model_name,
                    "fold": summary_name,
                    "n_train": "",
                    "n_test": sum(int(row["n_test"]) for row in fold_metric_rows),
                    "test_positive_count": sum(int(row["test_positive_count"]) for row in fold_metric_rows),
                    "test_negative_count": sum(int(row["test_negative_count"]) for row in fold_metric_rows),
                    "split_strategy": "StratifiedGroupKFold_exact_Murcko_scaffold",
                }
                for metric in ["roc_auc", "average_precision", "balanced_accuracy", "mcc", "brier_score"]:
                    values = [float(row[metric]) for row in fold_metric_rows if row[metric] != ""]
                    summary[metric] = round(function(values), 6) if values else ""
                cv_metrics.append(summary)

            final_model = sklearn.base.clone(estimator)
            final_model.fit(features_for(rows, model_name), y)
            fitted_models[(gene, model_name)] = final_model

    return (
        cv_metrics,
        oof_predictions,
        support_rows,
        fitted_models,
        training_by_gene,
        thresholds_by_gene,
    )


def predict_probability(model, model_name: str, row: dict) -> float:
    features = features_for([row], model_name)
    return float(model.predict_proba(features)[0, 1])


def external_roth_predictions(
    free_rows: list[dict],
    fitted_models: dict,
    training_by_gene: dict,
    thresholds_by_gene: dict,
) -> tuple[list[dict], list[dict]]:
    free_by_id = {row["reference_id"]: row for row in free_rows}
    source_rows = read_csv(PROCESSED / "external_truth_set_v1_pgp_bcrp_continuous_26.csv")
    predictions = []
    for source in source_rows:
        gene = source["transporter_gene"]
        ref = free_by_id[source["compound_id"]]
        for model_name in ["dummy_prior", "physchem_logistic", "morgan_logistic"]:
            model = fitted_models[(gene, model_name)]
            if ref["standardization_status"] != "standardized":
                probability = ""
                maximum = ""
                band = "structure_unresolved"
                abstain = "yes"
            else:
                probability = round(predict_probability(model, model_name, ref), 8)
                maximum_raw, band, abstain = support_for(
                    ref,
                    training_by_gene[gene],
                    thresholds_by_gene[gene],
                )
                maximum = round(maximum_raw, 6)
            predictions.append(
                {
                    "assay_record_id": source["assay_record_id"],
                    "compound_id": source["compound_id"],
                    "paper_compound_name": source["paper_compound_name"],
                    "transporter_gene": gene,
                    "model_name": model_name,
                    "resistance_ratio": source["resistance_ratio"],
                    "binary_substrate_label": "",
                    "predicted_probability": probability,
                    "max_clean_training_morgan_tanimoto": maximum,
                    "similarity_support_band": band,
                    "abstain_recommendation": abstain,
                    "evaluation_role": "continuous_external_stress_test_no_binarization",
                    "source_doi": source["source_doi"],
                    "source_url": source["source_url"],
                }
            )

    metrics = []
    for gene in ["ABCB1", "ABCG2"]:
        for model_name in ["dummy_prior", "physchem_logistic", "morgan_logistic"]:
            rows = [
                row
                for row in predictions
                if row["transporter_gene"] == gene
                and row["model_name"] == model_name
                and row["predicted_probability"] != ""
            ]
            ratios = np.array([float(row["resistance_ratio"]) for row in rows])
            probs = np.array([float(row["predicted_probability"]) for row in rows])
            rho_result = spearmanr(ratios, probs)
            outside = sum(row["abstain_recommendation"] == "yes" for row in rows)
            caution = sum(row["abstain_recommendation"] == "yes_caution" for row in rows)
            metrics.append(
                {
                    "transporter_gene": gene,
                    "model_name": model_name,
                    "n_compounds": len(rows),
                    "spearman_probability_vs_continuous_rr": (
                        round(float(rho_result.statistic), 6)
                        if math.isfinite(float(rho_result.statistic))
                        else ""
                    ),
                    "spearman_p_value": (
                        round(float(rho_result.pvalue), 6)
                        if math.isfinite(float(rho_result.pvalue))
                        else ""
                    ),
                    "outside_similarity_support_count": outside,
                    "low_similarity_support_count": caution,
                    "binary_threshold_used": "no",
                    "interpretation": "directional_association_only_small_external_panel",
                }
            )
    return predictions, metrics


def evidence_priority(tier: str) -> int:
    order = {
        "A_direct_numeric_inhibitor_rescue": 1,
        "A_direct_numeric_no_rescue_reported": 2,
        "A_direct_qualitative_inhibitor_control": 3,
        "B_functional_resistance_or_inhibitor_rescue": 4,
        "B_regulatory_call_no_numeric": 5,
        "C_comparative_cellular_retention": 6,
    }
    return order.get(tier, 99)


def released_evidence_predictions(
    release_rows: list[dict],
    fitted_models: dict,
    training_by_gene: dict,
    thresholds_by_gene: dict,
) -> tuple[list[dict], list[dict]]:
    release_by_id = {row["reference_id"]: row for row in release_rows}
    source_rows = read_csv(PROCESSED / "adc_release_species_transporter_evidence_long_v2.csv")
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for source in source_rows:
        if source["source_binary_label"] in {"0", "1"}:
            grouped[(source["release_species_id"], source["transporter_gene"])].append(source)

    predictions = []
    for (release_id, gene), members in sorted(grouped.items()):
        labels = sorted({row["source_binary_label"] for row in members})
        ref = release_by_id[release_id]
        chosen_tier = sorted((row["evidence_tier"] for row in members), key=evidence_priority)[0]
        stratum, endpoint_family, allowed_role, aggregation_rule = EVIDENCE_TIER_MAP[chosen_tier]
        for model_name in ["dummy_prior", "physchem_logistic", "morgan_logistic"]:
            if ref["standardization_status"] != "standardized" or len(labels) != 1:
                probability = ""
                maximum = ""
                band = "structure_or_label_conflict"
                abstain = "yes"
            else:
                model = fitted_models[(gene, model_name)]
                probability = round(predict_probability(model, model_name, ref), 8)
                maximum_raw, band, abstain = support_for(
                    ref,
                    training_by_gene[gene],
                    thresholds_by_gene[gene],
                )
                maximum = round(maximum_raw, 6)
            predictions.append(
                {
                    "release_species_id": release_id,
                    "release_species_name": ref["reference_name"],
                    "payload_core": ref.get("payload_core", ""),
                    "transporter_gene": gene,
                    "model_name": model_name,
                    "source_binary_label": labels[0] if len(labels) == 1 else "",
                    "model_binary_label": "",
                    "evidence_tier": chosen_tier,
                    "evidence_stratum": stratum,
                    "endpoint_family": endpoint_family,
                    "source_record_count_aggregated": len(members),
                    "source_record_ids": " | ".join(row["transport_record_id"] for row in members),
                    "source_effect_ratios": " | ".join(row["effect_ratio"] for row in members if row["effect_ratio"]),
                    "predicted_probability": probability,
                    "max_clean_training_morgan_tanimoto": maximum,
                    "similarity_support_band": band,
                    "abstain_recommendation": abstain,
                    "allowed_evaluation_role": allowed_role,
                    "aggregation_rule": aggregation_rule,
                    "interpretation": "descriptive_external_evidence_audit_not_model_label",
                }
            )

    metrics = []
    for gene in ["ABCB1", "ABCG2"]:
        for model_name in ["dummy_prior", "physchem_logistic", "morgan_logistic"]:
            rows = [
                row
                for row in predictions
                if row["transporter_gene"] == gene
                and row["model_name"] == model_name
                and row["predicted_probability"] != ""
            ]
            if not rows:
                continue
            y = np.array([int(row["source_binary_label"]) for row in rows])
            probs = np.array([float(row["predicted_probability"]) for row in rows])
            result = compute_metrics(y, probs)
            metrics.append(
                {
                    "transporter_gene": gene,
                    "model_name": model_name,
                    "evidence_stratum": "all_eligible_source_labels",
                    "n_entity_transporter_records": len(rows),
                    "positive_count": int(y.sum()),
                    "negative_count": int(len(y) - y.sum()),
                    **{key: round(value, 6) if math.isfinite(value) else "" for key, value in result.items()},
                    "metric_status": "descriptive_only_n_below_20",
                    "model_label_assignment": "none",
                }
            )
            for stratum in sorted({row["evidence_stratum"] for row in rows}):
                subset = [row for row in rows if row["evidence_stratum"] == stratum]
                metrics.append(
                    {
                        "transporter_gene": gene,
                        "model_name": model_name,
                        "evidence_stratum": stratum,
                        "n_entity_transporter_records": len(subset),
                        "positive_count": sum(int(row["source_binary_label"]) for row in subset),
                        "negative_count": len(subset) - sum(int(row["source_binary_label"]) for row in subset),
                        "roc_auc": "",
                        "average_precision": "",
                        "balanced_accuracy": "",
                        "mcc": "",
                        "brier_score": "",
                        "metric_status": "not_computed_within_tiny_evidence_stratum",
                        "model_label_assignment": "none",
                    }
                )
    return predictions, metrics


def pair_probability_shifts(
    free_rows: list[dict],
    release_rows: list[dict],
    fitted_models: dict,
    training_by_gene: dict,
    thresholds_by_gene: dict,
) -> list[dict]:
    free_by_id = {row["reference_id"]: row for row in free_rows}
    release_by_id = {row["reference_id"]: row for row in release_rows}
    source_pairs = read_csv(PROCESSED / "free_payload_release_pair_evidence_v1.csv")
    output = []
    gene = "ABCB1"
    for pair in source_pairs:
        free = free_by_id[pair["free_payload_compound_id"]]
        release = release_by_id[pair["release_species_id"]]
        for model_name in ["dummy_prior", "physchem_logistic", "morgan_logistic"]:
            model = fitted_models[(gene, model_name)]
            if free["standardization_status"] == "standardized":
                free_prob = predict_probability(model, model_name, free)
                free_support, free_band, free_abstain = support_for(
                    free, training_by_gene[gene], thresholds_by_gene[gene]
                )
            else:
                free_prob, free_support, free_band, free_abstain = None, None, "structure_unresolved", "yes"
            if release["standardization_status"] == "standardized":
                release_prob = predict_probability(model, model_name, release)
                release_support, release_band, release_abstain = support_for(
                    release, training_by_gene[gene], thresholds_by_gene[gene]
                )
            else:
                release_prob, release_support, release_band, release_abstain = None, None, "structure_unresolved", "yes"
            output.append(
                {
                    "pair_id": pair["pair_id"],
                    "payload_core": pair["payload_core"],
                    "free_payload_compound_id": pair["free_payload_compound_id"],
                    "free_payload_name": pair["free_payload_name"],
                    "release_species_id": pair["release_species_id"],
                    "release_species_name": pair["release_species_name"],
                    "transporter_gene": gene,
                    "model_name": model_name,
                    "free_predicted_probability": round(free_prob, 8) if free_prob is not None else "",
                    "release_predicted_probability": round(release_prob, 8) if release_prob is not None else "",
                    "probability_shift_release_minus_free": (
                        round(release_prob - free_prob, 8)
                        if release_prob is not None and free_prob is not None
                        else ""
                    ),
                    "free_max_training_similarity": round(free_support, 6) if free_support is not None else "",
                    "release_max_training_similarity": round(release_support, 6) if release_support is not None else "",
                    "free_similarity_support_band": free_band,
                    "release_similarity_support_band": release_band,
                    "free_abstain_recommendation": free_abstain,
                    "release_abstain_recommendation": release_abstain,
                    "observed_directional_comparison": pair["directional_comparison"],
                    "endpoint_compatibility": pair["endpoint_compatibility"],
                    "causal_interpretation": pair["causal_interpretation"],
                    "analysis_status": "model_structure_sensitivity_only_not_causal_validation",
                }
            )
    return output


def aggregate_cv_lookup(cv_metrics: list[dict]) -> dict[tuple[str, str], dict]:
    return {
        (row["transporter_gene"], row["model_name"]): row
        for row in cv_metrics
        if row["fold"] == "mean"
    }


def make_gate_summary(
    cv_metrics: list[dict],
    roth_metrics: list[dict],
    release_metrics: list[dict],
) -> list[dict]:
    cv_lookup = aggregate_cv_lookup(cv_metrics)
    roth_lookup = {(row["transporter_gene"], row["model_name"]): row for row in roth_metrics}
    release_lookup = {
        (row["transporter_gene"], row["model_name"]): row
        for row in release_metrics
        if row["evidence_stratum"] == "all_eligible_source_labels"
    }
    rows = []
    for gene in ["ABCB1", "ABCG2"]:
        for model_name in ["dummy_prior", "physchem_logistic", "morgan_logistic"]:
            cv = cv_lookup[(gene, model_name)]
            roth = roth_lookup[(gene, model_name)]
            release = release_lookup.get((gene, model_name), {})
            mcc = safe_float(cv["mcc"])
            bacc = safe_float(cv["balanced_accuracy"])
            internal_pass = bool(
                mcc is not None
                and bacc is not None
                and mcc >= INTERNAL_MCC_GATE
                and bacc >= INTERNAL_BALANCED_ACCURACY_GATE
            )
            rho = safe_float(roth["spearman_probability_vs_continuous_rr"])
            rows.append(
                {
                    "transporter_gene": gene,
                    "model_name": model_name,
                    "scaffold_cv_mean_mcc": cv["mcc"],
                    "scaffold_cv_mean_balanced_accuracy": cv["balanced_accuracy"],
                    "scaffold_cv_mean_roc_auc": cv["roc_auc"],
                    "internal_gate_rule": f"MCC >= {INTERNAL_MCC_GATE:.2f} AND balanced accuracy >= {INTERNAL_BALANCED_ACCURACY_GATE:.2f}",
                    "internal_gate_status": "pass_internal_only" if internal_pass else "fail_internal_gate",
                    "roth_external_n": roth["n_compounds"],
                    "roth_spearman_probability_vs_rr": roth["spearman_probability_vs_continuous_rr"],
                    "roth_directional_signal": (
                        "positive"
                        if rho is not None and rho > 0
                        else "negative"
                        if rho is not None and rho < 0
                        else "none_or_undefined"
                    ),
                    "release_evidence_n": release.get("n_entity_transporter_records", 0),
                    "release_metric_status": release.get("metric_status", "no_eligible_records"),
                    "recommended_use": (
                        "ranking_only_with_similarity_abstention"
                        if internal_pass
                        else "development_baseline_only"
                    ),
                    "final_predictor_claim_status": "not_eligible_no_independent_matched_release_validation",
                }
            )
    return rows


def save_models(
    fitted_models: dict,
    training_by_gene: dict,
    thresholds_by_gene: dict,
) -> list[dict]:
    manifest = []
    for (gene, model_name), model in sorted(fitted_models.items()):
        filename = f"{gene}_{model_name}_baseline.joblib"
        path = MODELS_DIR / filename
        bundle = {
            "transporter_gene": gene,
            "model_name": model_name,
            "model": model,
            "run_date": RUN_DATE,
            "random_seed": RANDOM_SEED,
            "feature_definition": (
                {"descriptor_fields": DESCRIPTOR_FIELDS}
                if model_name == "physchem_logistic"
                else {"morgan_radius": MORGAN_RADIUS, "morgan_bits": MORGAN_BITS, "include_chirality": False}
                if model_name == "morgan_logistic"
                else {"strategy": "class prior"}
            ),
            "standardization_rules": STANDARDIZATION_RULES,
            "similarity_support_thresholds": thresholds_by_gene[gene],
            "training_consensus_record_ids": [row["consensus_record_id"] for row in training_by_gene[gene]],
            "claim_boundary": "development baseline; not a deployable ADC released-species predictor",
        }
        joblib.dump(bundle, path, compress=3)
        manifest.append(
            {
                "filename": filename,
                "transporter_gene": gene,
                "model_name": model_name,
                "training_consensus_structure_count": len(training_by_gene[gene]),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "relative_path": str(path.relative_to(PROJECT_ROOT)),
                "claim_status": "development_baseline_only",
            }
        )
    return manifest


def evidence_mapping_rows() -> list[dict]:
    rows = []
    for source_tier, (stratum, endpoint, role, aggregation) in EVIDENCE_TIER_MAP.items():
        rows.append(
            {
                "source_evidence_tier": source_tier,
                "evidence_stratum": stratum,
                "endpoint_family": endpoint,
                "allowed_modeling_role": role,
                "source_binary_label_use": (
                    "external_audit_only" if "external" in role else "not_used_as_binary_target"
                ),
                "model_binary_label_assignment": "never_assigned_in_this_run",
                "aggregation_rule": aggregation,
                "rationale": "Preserve assay hierarchy and prevent direct transport, functional rescue, regulatory call and comparative retention from being treated as interchangeable truth.",
            }
        )
    return rows


def validate_outputs(
    free_rows: list[dict],
    release_rows: list[dict],
    training_rows: list[dict],
    decisions: list[dict],
    cv_metrics: list[dict],
    roth_predictions: list[dict],
    release_predictions: list[dict],
    gate_summary: list[dict],
) -> None:
    checks = [
        (len(free_rows) == 32, "expected 32 free-payload registry rows"),
        (len(release_rows) == 10, "expected 10 released-species registry rows"),
        (len(training_rows) == 1114, "expected 1114 Daood training records"),
        (
            sum(row["standardization_status"] == "standardized" for row in release_rows) == 9,
            "expected 9 standardized released-species structures",
        ),
        (
            all(row["model_binary_label"] == "" for row in release_predictions),
            "released-species model labels must remain blank",
        ),
        (
            all(
                not row["exact_external_reference_ids"]
                for row in decisions
                if row["leakage_model_action"] == "include_leakage_clean_training"
            ),
            "exact external overlap leaked into model training",
        ),
        (len([row for row in cv_metrics if row["fold"] == "mean"]) == 6, "missing CV aggregate rows"),
        (len(roth_predictions) == 26 * 3, "expected three model predictions for 26 Roth records"),
        (len(gate_summary) == 6, "expected six transporter-model gate rows"),
        (
            all(row["final_predictor_claim_status"].startswith("not_eligible") for row in gate_summary),
            "final predictor claim boundary changed",
        ),
    ]
    failed = [message for passed, message in checks if not passed]
    if failed:
        raise RuntimeError("Validation failed: " + "; ".join(failed))


def main() -> None:
    free_rows, release_rows, training_rows = make_standardized_registries()
    references = free_rows + release_rows
    consensus, conflicts = collapse_training(training_rows)
    reference_summary, flagged_pairs, decisions = audit_similarity(references, consensus)
    (
        cv_metrics,
        oof_predictions,
        support_rows,
        fitted_models,
        training_by_gene,
        thresholds_by_gene,
    ) = fit_baselines(decisions)
    roth_predictions, roth_metrics = external_roth_predictions(
        free_rows,
        fitted_models,
        training_by_gene,
        thresholds_by_gene,
    )
    release_predictions, release_metrics = released_evidence_predictions(
        release_rows,
        fitted_models,
        training_by_gene,
        thresholds_by_gene,
    )
    pair_shifts = pair_probability_shifts(
        free_rows,
        release_rows,
        fitted_models,
        training_by_gene,
        thresholds_by_gene,
    )
    gate_summary = make_gate_summary(cv_metrics, roth_metrics, release_metrics)
    manifest = save_models(fitted_models, training_by_gene, thresholds_by_gene)
    mapping = evidence_mapping_rows()

    validate_outputs(
        free_rows,
        release_rows,
        training_rows,
        decisions,
        cv_metrics,
        roth_predictions,
        release_predictions,
        gate_summary,
    )

    write_csv(PROCESSED / "compound_registry_standardized_v1.csv", free_rows)
    write_csv(PROCESSED / "adc_release_species_registry_standardized_v1.csv", release_rows)
    write_csv(PROCESSED / "daood_training_records_standardized_v1.csv", training_rows)
    write_csv(PROCESSED / "daood_training_consensus_v1.csv", consensus)
    write_csv(
        PROCESSED / "daood_structure_and_label_conflicts_v1.csv",
        conflicts,
        fieldnames=CONFLICT_FIELDS,
    )
    write_csv(PROCESSED / "similarity_leakage_reference_summary_v1.csv", reference_summary)
    write_csv(PROCESSED / "similarity_leakage_flagged_pairs_v1.csv", flagged_pairs)
    write_csv(PROCESSED / "daood_training_leakage_decisions_v1.csv", decisions)
    write_csv(PROCESSED / "evidence_tier_mapping_preregistered_v1.csv", mapping)
    write_csv(PROCESSED / "baseline_similarity_support_thresholds_v1.csv", support_rows)
    write_csv(PROCESSED / "baseline_scaffold_cv_metrics_v1.csv", cv_metrics)
    write_csv(PROCESSED / "baseline_scaffold_oof_predictions_v1.csv", oof_predictions)
    write_csv(PROCESSED / "baseline_roth_external_predictions_v1.csv", roth_predictions)
    write_csv(PROCESSED / "baseline_roth_external_metrics_v1.csv", roth_metrics)
    write_csv(PROCESSED / "baseline_release_evidence_predictions_v1.csv", release_predictions)
    write_csv(PROCESSED / "baseline_release_evidence_metrics_v1.csv", release_metrics)
    write_csv(PROCESSED / "baseline_pair_probability_shifts_v1.csv", pair_shifts)
    write_csv(PROCESSED / "baseline_model_gate_summary_v1.csv", gate_summary)
    write_csv(PROCESSED / "baseline_model_manifest_v1.csv", manifest)

    (PROCESSED / "structure_standardization_and_baseline_config_v1.json").write_text(
        json.dumps(STANDARDIZATION_RULES, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    counts_by_action = defaultdict(int)
    for row in decisions:
        counts_by_action[row["leakage_model_action"]] += 1
    counts_by_status = defaultdict(int)
    for row in references:
        counts_by_status[row["standardization_status"]] += 1
    reference_tautomer_counts = defaultdict(int)
    for row in references:
        reference_tautomer_counts[row["tautomer_enumeration_status"]] += 1
    training_tautomer_counts = defaultdict(int)
    for row in training_rows:
        training_tautomer_counts[row["tautomer_enumeration_status"]] += 1
    summary = {
        "run_date": RUN_DATE,
        "runtime": {
            "python": platform.python_version(),
            "rdkit": rdBase.rdkitVersion,
            "scikit_learn": sklearn.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "joblib": joblib.__version__,
        },
        "standardization": {
            "free_payload_rows": len(free_rows),
            "release_species_rows": len(release_rows),
            "reference_standardization_status_counts": dict(sorted(counts_by_status.items())),
            "reference_tautomer_enumeration_status_counts": dict(sorted(reference_tautomer_counts.items())),
            "daood_source_records": len(training_rows),
            "daood_tautomer_enumeration_status_counts": dict(sorted(training_tautomer_counts.items())),
            "daood_consensus_structure_transporter_records": len(consensus),
            "structure_or_label_conflict_rows": len(conflicts),
        },
        "leakage": {
            "reference_objects_audited": len(reference_summary),
            "flagged_similarity_pairs": len(flagged_pairs),
            "training_action_counts": dict(sorted(counts_by_action.items())),
            "thresholds": {
                "near_duplicate": NEAR_DUPLICATE_THRESHOLD,
                "complex_scaffold_similarity": SCAFFOLD_FAMILY_THRESHOLD,
                "complex_scaffold_min_heavy_atoms": COMPLEX_SCAFFOLD_MIN_HEAVY_ATOMS,
            },
        },
        "models": {
            "ABCB1_clean_training_records": len(training_by_gene["ABCB1"]),
            "ABCG2_clean_training_records": len(training_by_gene["ABCG2"]),
            "saved_model_artifacts": len(manifest),
            "Roth_external_prediction_rows": len(roth_predictions),
            "released_species_external_prediction_rows": len(release_predictions),
            "pair_probability_shift_rows": len(pair_shifts),
            "model_binary_labels_assigned": 0,
            "final_predictor_claims": 0,
        },
        "gate_summary": gate_summary,
        "warning": "Internal scaffold CV and small external audits are development evidence only. A credible ADC released-species predictor still requires independent matched transport validation.",
    }
    (PROCESSED / "standardization_leakage_baseline_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
