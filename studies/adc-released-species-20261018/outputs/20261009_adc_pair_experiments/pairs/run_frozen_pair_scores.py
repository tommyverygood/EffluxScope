#!/usr/bin/env python3
"""Score ADC parent/release structures with the existing frozen v0.4 package.

Writes only alongside this script. No fitting, threshold changes, truth labels,
combined dual-payload score, or modification of processed data is performed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
MODEL_DIR = PROJECT / "模型/ABCB1_selective_screener_v0_4_dev"
sys.path.insert(0, str(MODEL_DIR))
from predict_abcb1_selective import flatten_batch_result, load_packages, predict_one
from abcb1_screener_core_v0_1 import standardize_smiles
from rdkit import Chem
import joblib
import numpy as np
import rdkit
import sklearn


AUDIT = PROJECT / "outputs/20261009_adc_central_pairing/local_pair_audit.csv"
FREE = PROJECT / "数据/processed/compound_registry_v1.csv"
RELEASE = PROJECT / "outputs/20260810_adc_demo_reconstruction_v2/adc_application_demo_registry_reconstructed_v2.csv"
SELECTIVE = MODEL_DIR / "ABCB1_selective_uncertainty_v0_4_dev.joblib"
PRIMARY = PROJECT / "模型/ABCB1_screener_v0_1/ABCB1_screener_v0_1.joblib"
DIRECT = PROJECT / "模型/ABCB1_assay_aware_candidate_v0_3_dev/ABCB1_assay_aware_candidate_v0_3_dev.joblib"


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h))


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def j(value) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def numeric(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--approved-structures", type=Path,
                        help="Optional independently reviewed correction CSV; strict admission fields required")
    parser.add_argument("--output-dir", type=Path, default=HERE / "results_v1",
                        help="Fresh result directory; existing result files are never overwritten")
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if not output_dir.is_relative_to(HERE):
        raise ValueError("Experiment result directory must be within the new pairs directory")
    output_dir.mkdir(parents=True, exist_ok=True)
    if (output_dir / "scoring_summary.json").exists():
        raise ValueError("Existing result bundle is protected; choose a fresh --output-dir")
    frozen_inputs = [AUDIT, FREE, RELEASE, SELECTIVE, PRIMARY, DIRECT,
                     MODEL_DIR / "predict_abcb1_selective.py", Path(__file__).resolve()]
    overrides = {}
    if args.approved_structures:
        frozen_inputs.append(args.approved_structures.resolve())
        for row in read_csv(args.approved_structures):
            if (row.get("approved_for_frozen_scoring") != "yes"
                    or row.get("independently_reviewed") != "yes"
                    or not row.get("exact_isomeric_smiles")
                    or not row.get("primary_structure_source_url")):
                raise ValueError("Correction CSV lacks explicit approval, independent review, exact structure, or primary source")
            overrides[row["entity_id"]] = row
    before = {str(p.relative_to(PROJECT)) if p.is_relative_to(PROJECT) else str(p): digest(p)
              for p in frozen_inputs}
    selective, primary, direct = load_packages(SELECTIVE, PRIMARY, DIRECT)
    audit = read_csv(AUDIT)
    free = {r["compound_id"]: r for r in read_csv(FREE)}
    release = {r["entity_id"]: r for r in read_csv(RELEASE)}
    entity_ids = list(dict.fromkeys(x for pair in audit for x in [pair["free_id"], pair["released_id"]]))
    entities = {}
    for entity_id in entity_ids:
        if entity_id in release:
            r = release[entity_id]
            entities[entity_id] = {
                "entity_id": entity_id, "entity_name": r["entity_name"],
                "source_input_smiles": r["source_isomeric_smiles"],
                "source_inchikey": r["source_inchikey"],
                "object_type": "release_species", "release_stage": r["release_stage"],
                "release_object_class": r["audited_object_class"],
                "identity_status": r["structure_identity_verdict"],
                "identity_limitation": r["prohibited_claim"],
                "primary_structure_source_url": r["source_structure_url"],
                "input_registry_path": str(RELEASE.relative_to(PROJECT)),
                "identity_hold": "",
            }
        else:
            r = free[entity_id]
            entities[entity_id] = {
                "entity_id": entity_id, "entity_name": r["paper_compound_name"],
                "source_input_smiles": r["isomeric_smiles"], "source_inchikey": r["inchikey"],
                "object_type": "free_payload", "release_stage": "free_parent_reference",
                "release_object_class": "free_payload", "identity_status": r["structure_status"],
                "identity_limitation": r["structure_note"],
                "primary_structure_source_url": r["structure_source_url"],
                "input_registry_path": str(FREE.relative_to(PROJECT)), "identity_hold": "",
            }
    # These are explicit chemistry holds, separate from a model abstention.
    holds = {
        "ADC_REL_002": "unresolved_Lys_PEG4Mal_DM1_machine_structure",
        "ADC_REL_004": "wrong_S_methyl_DM4_registry_maps_two_sulfur_methyl_disulfide",
        "ADC_REL_005": "S_methyl_DM1_stereochemistry_unspecified_not_exact_material",
        "ADC_PAYLOAD_029": "Roth_Deruxtecan_material_ambiguous_C26_DXd_vs_C52_linker_payload",
    }
    for entity_id, reason in holds.items():
        if entity_id in entities:
            entities[entity_id]["identity_hold"] = reason
    for entity_id, correction in overrides.items():
        if entity_id != "ADC_REL_004":
            raise ValueError("Only the explicitly curated S-methyl-DM4 correction is admitted by this experiment")
        entities[entity_id].update(
            source_input_smiles=correction["exact_isomeric_smiles"],
            source_inchikey=correction.get("exact_inchikey", ""),
            identity_status="independently_reviewed_corrected_exact_structure",
            identity_hold="", primary_structure_source_url=correction["primary_structure_source_url"],
            input_registry_path=str(args.approved_structures.resolve()),
            identity_limitation="Corrected representation is admitted only by explicit independent review; not a new experimental transporter label.",
        )
    # An identity-preserving control uses the actual released chemical identity.
    # It is NOT the ambiguous Roth Deruxtecan structure stored under PAYLOAD_029.
    control_parent_ids = {}
    for pair in audit:
        if pair["pair_type"] == "identity_preserving_control":
            parent_id = "CONTROL_PARENT::" + pair["released_id"]
            control_parent_ids[pair["audit_record_id"]] = parent_id
            entities[parent_id] = {
                **entities[pair["released_id"]], "entity_id": parent_id,
                "entity_name": pair["released_name"] + " free parent identity control",
                "object_type": "free_payload", "release_stage": "identity_preserving_control",
                "release_object_class": "free_payload",
            }
    entities["DUAL_MATERIAL::TRIPTOLIDE_UNVERIFIED"] = {
        "entity_id": "DUAL_MATERIAL::TRIPTOLIDE_UNVERIFIED", "entity_name": "Triptolide (dual-paper exact material unverified)",
        "source_input_smiles": "", "source_inchikey": "", "object_type": "free_payload",
        "release_stage": "dual_candidate_material_unverified", "release_object_class": "free_payload",
        "identity_status": "HOLD", "identity_limitation": "No independently verified source-specific structure admitted",
        "primary_structure_source_url": "", "input_registry_path": "",
        "identity_hold": "triptolide_dual_source_exact_material_not_independently_confirmed",
    }
    scores_by_structure = {}
    entity_rows = []
    unique_rows = []
    primary_keys = set(primary["training_structure_keys"])
    direct_keys = {r["standardized_inchikey"] for r in direct["training_rows"]}
    for entity in entities.values():
        base = {**entity, "score_status": "identity_HOLD_not_scored", "unique_structure_id": "",
                "effective_analysis_status": "identity_HOLD", "supports_independent_validation": "no"}
        if entity["identity_hold"] or not entity["source_input_smiles"]:
            entity_rows.append(base)
            continue
        normalized = standardize_smiles(entity["source_input_smiles"])
        key = normalized["standardized_inchikey"]
        if normalized["standardization_status"] != "ok" or not key:
            entity_rows.append({**base, "identity_hold": "standardization_failed", **normalized})
            continue
        if key not in scores_by_structure:
            result = predict_one(entity["source_input_smiles"], selective, primary, direct,
                                 chemical_object_type="released_or_free_small_molecule",
                                 requested_context="structure_only_screening")
            raw_mol = Chem.MolFromSmiles(entity["source_input_smiles"])
            raw_key = Chem.MolToInchiKey(raw_mol) if raw_mol is not None else ""
            query_key = result.get("standardized_inchikey", "")
            primary_exact = any(x in primary_keys for x in [raw_key, query_key, key] if x)
            direct_exact = any(x in direct_keys for x in [raw_key, query_key, key] if x)
            contamination = primary_exact or direct_exact or result["exact_reference_identity"] == "yes"
            ps = numeric(result.get("primary_max_tanimoto"))
            ds = numeric(result.get("direct_max_tanimoto"))
            primary_support = ps is not None and ps >= float(selective["primary_similarity_q10"])
            direct_support = ds is not None and ds >= float(selective["direct_similarity_q10"])
            mask = ("both" if primary_support and direct_support else "primary_only" if primary_support
                    else "direct_only" if direct_support else "neither")
            unique_id = "STRUCT_" + str(len(scores_by_structure) + 1).zfill(3)
            scores_by_structure[key] = {
                "unique_structure_id": unique_id, **flatten_batch_result(result),
                "v0_1_risk_score": result["primary_calibrated_probability"],
                "v0_3_risk_score": result["direct_panel_probability"],
                "primary_exact_development_identity": "yes" if primary_exact else "no",
                "direct_exact_development_identity": "yes" if direct_exact else "no",
                "any_exact_development_identity": "yes" if contamination else "no",
                "supported_mask": mask, "primary_supported": "yes" if primary_support else "no",
                "direct_supported": "yes" if direct_support else "no",
                "primary_similarity_support_threshold": selective["primary_similarity_q10"],
                "direct_similarity_support_threshold": selective["direct_similarity_q10"],
                "score_status": "scored_frozen_model",
                "effective_analysis_status": "development_identity_contaminated_not_validation" if contamination
                    else "model_abstain" if result["risk_call"] == "abstain" else "accepted_model_risk_screen",
                "supports_independent_validation": "no",
                "biological_direction_truth": "not_available_matched_M0_pair",
            }
            unique_rows.append({"first_entity_id": entity["entity_id"], "first_entity_name": entity["entity_name"],
                                "dedup_standardized_key": key, **scores_by_structure[key]})
        entity_rows.append({**base, **scores_by_structure[key]})
    by_id = {r["entity_id"]: r for r in entity_rows}
    contrasts = []
    identity_checks = []
    for pair in audit:
        parent_id = control_parent_ids.get(pair["audit_record_id"], pair["free_id"])
        a, b = by_id[parent_id], by_id[pair["released_id"]]
        row = {**pair, "scored_parent_entity_id": parent_id,
               "parent_score_status": a["score_status"], "release_score_status": b["score_status"],
               "parent_identity_hold": a["identity_hold"], "release_identity_hold": b["identity_hold"],
               "parent_unique_structure_id": a["unique_structure_id"], "release_unique_structure_id": b["unique_structure_id"],
               "parent_full_decision": a.get("risk_call", "identity_HOLD"),
               "release_full_decision": b.get("risk_call", "identity_HOLD"),
               "parent_effective_analysis_status": a["effective_analysis_status"],
               "release_effective_analysis_status": b["effective_analysis_status"],
               "parent_supported_mask": a.get("supported_mask", "not_scored"),
               "release_supported_mask": b.get("supported_mask", "not_scored"),
               "parent_reasons_json": a.get("abstention_reasons_json", j([a["identity_hold"]])),
               "release_reasons_json": b.get("abstention_reasons_json", j([b["identity_hold"]])),
               "paired_biological_delta_truth": "not_available_no_matched_M0_pair",
               "delta_interpretation": "score_change_only_not_delta_ER_or_validated_efflux_direction",
               "accuracy_or_causal_claim_allowed": "no"}
        for version in ["v0_1", "v0_3"]:
            x, y = numeric(a.get(version + "_risk_score")), numeric(b.get(version + "_risk_score"))
            row["parent_" + version + "_risk_score"] = "" if x is None else x
            row["release_" + version + "_risk_score"] = "" if y is None else y
            row[version + "_release_minus_parent_score"] = "" if x is None or y is None else round(y - x, 8)
        if pair["pair_type"] == "identity_preserving_control":
            passed = (a["unique_structure_id"] == b["unique_structure_id"]
                      and all(row[v + "_release_minus_parent_score"] == 0 for v in ["v0_1", "v0_3"])
                      and row["parent_full_decision"] == row["release_full_decision"])
            identity_checks.append({"audit_record_id": pair["audit_record_id"], "passed": passed})
            row["identity_control_exact_score_equality"] = "pass" if passed else "FAIL"
        contrasts.append(row)
    dual_specs = [
        ("DUAL_COMPONENT_001", "ADC_PAYLOAD_017", "ADC_PAYLOAD_018", "free_MMAE_plus_free_MMAF", ""),
        ("DUAL_COMPONENT_002", "ADC_PAYLOAD_017", "ADC_REL_006", "free_MMAE_plus_cys_mcMMAF", ""),
        ("DUAL_COMPONENT_003_HOLD", "ADC_REL_011", "DUAL_MATERIAL::TRIPTOLIDE_UNVERIFIED",
         "exatecan_plus_triptolide_source_material_hold", "source_specific_dual_material_structures_not_independently_confirmed"),
    ]
    duals = []
    for dual_id, a_id, b_id, design, hold in dual_specs:
        a, b = by_id[a_id], by_id[b_id]
        row = {"dual_component_comparison_id": dual_id, "design": design,
               "comparison_status": "identity_HOLD" if hold else "individual_frozen_scores_completed",
               "left_entity_id": a_id, "right_entity_id": b_id,
               "left_entity_name": a["entity_name"], "right_entity_name": b["entity_name"],
               "joint_or_synergy_model": "none", "joint_score": "", "synergy_truth": "not_available",
               "joint_ER_truth": "not_available", "dual_ADC_experiment_performed": "no",
               "interpretation": "separate_ABCB1_component_risk_screen_only_no_combined_prediction",
               "identity_hold": hold, "independent_validation_claim_allowed": "no"}
        for side, ent in [("left", a), ("right", b)]:
            for field in ["v0_1_risk_score", "v0_3_risk_score", "risk_call", "supported_mask",
                          "any_exact_development_identity", "abstention_reasons_json"]:
                row[side + "_" + field] = "" if hold else ent.get(field, "")
        duals.append(row)
    after = {str(p.relative_to(PROJECT)) if p.is_relative_to(PROJECT) else str(p): digest(p)
             for p in frozen_inputs}
    if before != after:
        raise RuntimeError("A frozen input changed during scoring; refuse to report a stable run")
    if not all(c["passed"] for c in identity_checks):
        raise RuntimeError("Identity-preserving control score equality failed")
    write_csv(output_dir / "entities_scored_or_hold.csv", entity_rows)
    write_csv(output_dir / "unique_structure_predictions.csv", unique_rows)
    write_csv(output_dir / "parent_release_score_contrasts.csv", contrasts)
    write_csv(output_dir / "dual_component_scores.csv", duals)
    summary = {
        "status": "completed_frozen_computational_experiment_not_wet_lab_validation",
        "model_version": selective["model_version"], "no_refit_or_threshold_tuning": True,
        "old_inputs_unchanged": before == after,
        "entity_rows_including_control_aliases": len(entity_rows),
        "unique_structures_scored": len(unique_rows),
        "identity_hold_entity_count": sum(r["score_status"] == "identity_HOLD_not_scored" for r in entity_rows),
        "unique_structure_full_decisions": dict(Counter(r["risk_call"] for r in unique_rows)),
        "unique_structure_exact_development_identities": sum(r["any_exact_development_identity"] == "yes" for r in unique_rows),
        "chemical_changed_pair_count": sum(r["pair_type"] == "chemical_changed_parent_to_catabolite" for r in contrasts),
        "chemical_changed_pairs_both_scored": sum(r["pair_type"] == "chemical_changed_parent_to_catabolite"
            and r["parent_score_status"] == r["release_score_status"] == "scored_frozen_model" for r in contrasts),
        "identity_preserving_control_checks": identity_checks,
        "identity_preserving_control_checks_all_pass": all(c["passed"] for c in identity_checks),
        "dual_individual_component_comparisons_completed": sum(r["comparison_status"] == "individual_frozen_scores_completed" for r in duals),
        "matched_M0_biological_transport_pairs": 0,
        "output_dir": str(output_dir),
        "claim_boundary": "Only frozen individual structure scores and abstention behavior were measured; paired score deltas are not transporter effects, ADC efficacy, or synergy.",
        "input_hashes": before,
        "runtime": {"executable": sys.executable, "python": platform.python_version(), "numpy": np.__version__,
                    "scikit_learn": sklearn.__version__, "joblib": joblib.__version__, "rdkit": rdkit.__version__},
    }
    (output_dir / "scoring_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ["input_hashes", "runtime"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
