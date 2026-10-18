#!/usr/bin/env python3
"""Compare wrong legacy S2 representation with a separately approved S-methyl-DM4 structure.

The wrong structure is a diagnostic input, never a biological material or paired
transport observation. All results use the original frozen predictor packages.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
MODEL_DIR = PROJECT / "模型/ABCB1_selective_screener_v0_4_dev"
sys.path.insert(0, str(MODEL_DIR))
from predict_abcb1_selective import flatten_batch_result, load_packages, predict_one
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--approved-structures", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, default=HERE / "identity_sensitivity_v1")
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if not output_dir.is_relative_to(HERE):
        raise ValueError("Write only within the new pairs directory")
    output_dir.mkdir(parents=True, exist_ok=True)
    if (output_dir / "identity_sensitivity_summary.json").exists():
        raise ValueError("Existing result bundle is protected; choose a fresh --output-dir")
    corrected = next(r for r in read_csv(args.approved_structures) if r["entity_id"] == "ADC_REL_004")
    if corrected.get("approved_for_frozen_scoring") != "yes" or corrected.get("independently_reviewed") != "yes":
        raise ValueError("New representation lacks explicit independent review approval")
    legacy_registry = PROJECT / "outputs/20260810_adc_demo_reconstruction_v2/adc_application_demo_registry_reconstructed_v2.csv"
    legacy_predictions = PROJECT / "outputs/20260810_adc_demo_reconstruction_v2/abcb1_adc_application_demo_predictions_v2.csv"
    old = next(r for r in read_csv(legacy_registry) if r["entity_id"] == "ADC_REL_004")
    historical = next(r for r in read_csv(legacy_predictions) if r["entity_id"] == "ADC_REL_004")
    selective_path = MODEL_DIR / "ABCB1_selective_uncertainty_v0_4_dev.joblib"
    primary_path = PROJECT / "模型/ABCB1_screener_v0_1/ABCB1_screener_v0_1.joblib"
    direct_path = PROJECT / "模型/ABCB1_assay_aware_candidate_v0_3_dev/ABCB1_assay_aware_candidate_v0_3_dev.joblib"
    inputs = [args.approved_structures.resolve(), legacy_registry, legacy_predictions,
              selective_path, primary_path, direct_path, Path(__file__).resolve()]
    before = {str(p): digest(p) for p in inputs}
    selective, primary, direct = load_packages(selective_path, primary_path, direct_path)
    specs = [
        {"representation_id": "LEGACY_WRONG_S2", "material_status": "diagnostic_not_biological_material",
         "input_smiles": old["source_isomeric_smiles"], "source_url": old["source_structure_url"],
         "allowed_in_biological_pair_statistics": "no"},
        {"representation_id": "CURATED_S_METHYL_DM4", "material_status": "independently_reviewed_actual_metabolite_representation",
         "input_smiles": corrected["exact_isomeric_smiles"], "source_url": corrected["primary_structure_source_url"],
         "allowed_in_biological_pair_statistics": "structure_score_only_not_transport_truth"},
    ]
    predictions = []
    for spec in specs:
        mol = Chem.MolFromSmiles(spec["input_smiles"])
        if mol is None:
            raise ValueError("Approved or diagnostic input does not parse")
        result = predict_one(spec["input_smiles"], selective, primary, direct,
                             chemical_object_type="released_small_molecule",
                             requested_context="structure_only_screening")
        predictions.append({**spec, "source_inchikey": Chem.MolToInchiKey(mol),
                            "source_formula": rdMolDescriptors.CalcMolFormula(mol),
                            "source_molecular_weight": round(Descriptors.MolWt(mol), 6),
                            **flatten_batch_result(result), "biological_direction_truth": "not_available"})
    checks = [
        {"field": field, "fresh_wrong_structure_value": predictions[0][field],
         "legacy_value": historical[field], "pass": str(predictions[0][field]) == str(historical[field])}
        for field in ["primary_calibrated_probability", "direct_panel_probability", "risk_call", "abstention_reasons_json"]
    ]
    if not all(r["pass"] for r in checks):
        raise RuntimeError("The fresh diagnostic wrong representation does not reproduce legacy output")
    after = {str(p): digest(p) for p in inputs}
    if before != after:
        raise RuntimeError("Frozen input changed during sensitivity run")
    summary = {
        "status": "completed_frozen_representation_sensitivity_not_biological_experiment",
        "old_wrong_input_only_diagnostic_not_counted_as_pair": True,
        "legacy_wrong_structure_reproduction_checks": checks,
        "legacy_wrong_structure_reproduction_all_pass": all(r["pass"] for r in checks),
        "corrected_minus_wrong_v0_1_score": round(float(predictions[1]["primary_calibrated_probability"])
                                                - float(predictions[0]["primary_calibrated_probability"]), 8),
        "corrected_minus_wrong_v0_3_score": round(float(predictions[1]["direct_panel_probability"])
                                                - float(predictions[0]["direct_panel_probability"]), 8),
        "wrong_full_decision": predictions[0]["risk_call"], "corrected_full_decision": predictions[1]["risk_call"],
        "old_inputs_unchanged": before == after, "input_hashes": before,
        "claim_boundary": "The effect measured is a representation correction on frozen scores, not a chemistry intervention, transporter effect, or model accuracy estimate.",
    }
    write_csv(output_dir / "identity_sensitivity_scores.csv", predictions)
    (output_dir / "identity_sensitivity_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "input_hashes"}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
