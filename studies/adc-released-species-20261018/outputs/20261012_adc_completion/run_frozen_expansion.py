#!/usr/bin/env python3
"""Reproduce the frozen panel and score newly source-resolved ADC species.

No training, truth imputation, threshold changes, or edits to older inputs.
Run with the recorded .venv-experiments interpreter. Writes into scoring/.
"""
from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "scoring"
MODEL = ROOT / "模型/ABCB1_selective_screener_v0_4_dev"
sys.path.insert(0, str(MODEL))
print("Importing frozen model modules", flush=True)
from predict_abcb1_selective import flatten_batch_result, load_packages, predict_one
from abcb1_screener_core_v0_1 import standardize_smiles
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import joblib
import numpy as np
import rdkit
import sklearn


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def main():
    OUT.mkdir(exist_ok=True)
    if (OUT / "completion.json").exists():
        raise RuntimeError("Completed results are immutable; use a separately versioned run.")
    prior_path = ROOT / "outputs/20261009_adc_pair_experiments/pairs/results_curated_v2/entities_scored_or_hold.csv"
    new_path = HERE / "family_expansion/new_machine_entities.csv"
    input_paths = [prior_path, new_path, HERE / "protocol.json", Path(__file__),
                   MODEL / "predict_abcb1_selective.py",
                   ROOT / "数据/scripts/abcb1_screener_core_v0_1.py",
                   MODEL / "ABCB1_selective_uncertainty_v0_4_dev.joblib",
                   ROOT / "模型/ABCB1_screener_v0_1/ABCB1_screener_v0_1.joblib",
                   ROOT / "模型/ABCB1_assay_aware_candidate_v0_3_dev/ABCB1_assay_aware_candidate_v0_3_dev.joblib"]
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in input_paths}
    dump(OUT / "prerun_inputs.json", {"created_utc": datetime.now(timezone.utc).isoformat(), "input_sha256": hashes})
    seen = set()
    entities = []
    for row in read(prior_path):
        if row["score_status"] != "scored_frozen_model" or row["unique_structure_id"] in seen:
            continue
        seen.add(row["unique_structure_id"])
        entities.append({"entity_id": row["entity_id"], "entity_name": row["entity_name"],
                         "source_smiles": row["source_input_smiles"], "source_inchikey": row["source_inchikey"],
                         "role": row["object_type"], "family": "previous_corrected_panel",
                         "input_source": row["primary_structure_source_url"], "prior": row,
                         "identity_limitation": row["identity_limitation"], "newly_resolved": False})
    for row in read(new_path):
        if row["identity_review_status"] != "resolved" or not row["primary_structure_source_url"]:
            raise RuntimeError("New structure requires explicit source-resolved identity.")
        entities.append({"entity_id": row["entity_id"], "entity_name": row["entity_name"],
                         "source_smiles": row["exact_isomeric_smiles"], "source_inchikey": row["exact_inchikey"],
                         "role": row["object_role"], "family": row["family"],
                         "input_source": row["primary_structure_source_url"], "prior": None,
                         "identity_limitation": row["identity_limitation"], "newly_resolved": True})
    print(f"Loading fixed packages for {len(entities)} unique source structures", flush=True)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        selective, primary, direct = load_packages(input_paths[-3], input_paths[-2], input_paths[-1])
        pk = set(primary["training_structure_keys"])
        dk = {r["standardized_inchikey"] for r in direct["training_rows"]}
        pconn, dconn = {k.split("-")[0] for k in pk}, {k.split("-")[0] for k in dk}
        results, reproduction = [], []
        for entity in entities:
            mol = Chem.MolFromSmiles(entity["source_smiles"])
            if mol is None or Chem.MolToInchiKey(mol) != entity["source_inchikey"]:
                raise RuntimeError(f"Source SMILES/key mismatch for {entity['entity_id']}")
            standard = standardize_smiles(entity["source_smiles"])
            result = predict_one(entity["source_smiles"], selective, primary, direct,
                                 chemical_object_type="released_or_free_small_molecule",
                                 requested_context="structure_only_screening")
            keys = {entity["source_inchikey"], result.get("standardized_inchikey", "")} - {""}
            conn = {k.split("-")[0] for k in keys}
            flat = flatten_batch_result(result)
            row = {k: v for k, v in entity.items() if k != "prior"}
            row.update(flat)
            row.update(source_formula=rdMolDescriptors.CalcMolFormula(mol),
                       source_key_confirmed=True, tautomer_status=standard["tautomer_status"],
                       enumeration_complete=standard["tautomer_status"] == "Completed",
                       diagnostic_enumeration_hold=standard["tautomer_status"] != "Completed",
                       primary_exact_identity=bool(keys & pk), direct_exact_identity=bool(keys & dk),
                       primary_connectivity_overlap=bool(conn & pconn), direct_connectivity_overlap=bool(conn & dconn),
                       identity_is_not_biological_validation=True)
            results.append(row)
            old = entity["prior"]
            if old is not None:
                for field in ["primary_calibrated_probability", "direct_panel_probability"]:
                    difference = abs(float(result[field]) - float(old[field]))
                    reproduction.append({"entity_id": entity["entity_id"], "field": field,
                                         "old": old[field], "current": result[field],
                                         "absolute_difference": difference, "passed": difference <= 1e-8})
                for field in ["risk_call", "standardized_inchikey"]:
                    reproduction.append({"entity_id": entity["entity_id"], "field": field,
                                         "old": old[field], "current": result[field],
                                         "absolute_difference": "", "passed": result[field] == old[field]})
            print(entity["entity_id"], result["primary_calibrated_probability"],
                  result["direct_panel_probability"], result["risk_call"], standard["tautomer_status"], flush=True)
    if not all(r["passed"] for r in reproduction):
        raise RuntimeError("Historical frozen-panel reproduction failed.")
    after = {str(p.relative_to(ROOT)): sha(p) for p in input_paths}
    if after != hashes:
        raise RuntimeError("Frozen input mutated during experiment.")
    write(OUT / "species_predictions.csv", results)
    write(OUT / "historical_reproduction.csv", reproduction)
    dump(OUT / "completion.json", {
        "status": "PASS_frozen_computation_only", "completed_utc": datetime.now(timezone.utc).isoformat(),
        "structures_scored": len(results), "prior_structures_reproduced": len(seen),
        "reproduction_checks": len(reproduction), "all_reproduction_checks_pass": True,
        "new_structures_scored": sum(r["newly_resolved"] for r in results),
        "new_identity_or_connectivity_development_overlap": [r["entity_id"] for r in results if r["newly_resolved"] and
            any(r[f] for f in ["primary_exact_identity", "direct_exact_identity", "primary_connectivity_overlap", "direct_connectivity_overlap"])],
        "new_structure_enumeration_holds": [r["entity_id"] for r in results if r["newly_resolved"] and r["diagnostic_enumeration_hold"]],
        "input_sha256": hashes, "frozen_inputs_unchanged": True, "no_refit_or_selector_change": True,
        "new_wet_lab_experiments": 0, "new_measured_parent_release_transport_pairs": 0,
        "runtime": {"executable": sys.executable, "python": platform.python_version(),
                    "numpy": np.__version__, "rdkit": rdkit.__version__, "sklearn": sklearn.__version__, "joblib": joblib.__version__},
        "warnings": sorted({f"{w.category.__name__}: {w.message}" for w in captured}),
        "claim_boundary": "Frozen component scores and reproducibility, not ADC accuracy, transport change, joint transport, synergy, or efficacy. New exact structure can still have limited assay-material mapping."})


if __name__ == "__main__":
    main()
