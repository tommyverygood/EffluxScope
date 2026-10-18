#!/usr/bin/env python3
"""Frozen ADC structure/feature sensitivity audit; no fitting or entrypoint edits."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
import warnings
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "数据/scripts"))
import abcb1_screener_core_v0_1 as core
import joblib
import numpy as np
import rdkit
import sklearn
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

BASE = ROOT / "outputs/20261009_adc_pair_experiments/pairs/results_curated_v2"
PRIMARY = ROOT / "模型/ABCB1_screener_v0_1/ABCB1_screener_v0_1.joblib"
DIRECT = ROOT / "模型/ABCB1_assay_aware_candidate_v0_3_dev/ABCB1_assay_aware_candidate_v0_3_dev.joblib"
VARIANTS = ["source", "no_tautomer", "stereo_preserved_tautomer", "frozen_old", "source_stereo_stripped"]
DESC_NAMES = ["molecular_weight", "xlogp", "tpsa", "hbond_donors", "hbond_acceptors", "rotatable_bonds", "formal_charge", "heavy_atoms"]

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def js(x):
    return json.dumps(x, ensure_ascii=False, sort_keys=True, separators=(",", ":"))

def read_csv(p):
    with p.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def write_csv(name, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (HERE / name).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def mol_record(mol):
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return {"isomeric_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            "inchikey": Chem.MolToInchiKey(mol),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "formal_charge": sum(a.GetFormalCharge() for a in mol.GetAtoms()),
            "assigned_tetrahedral": sum(s != "?" for _, s in centers),
            "unassigned_tetrahedral": sum(s == "?" for _, s in centers),
            "assigned_double_bond_stereo": sum(b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS) for b in mol.GetBonds()),
            "atom_count": mol.GetNumAtoms(), "fragment_count": len(Chem.GetMolFrags(mol))}

def annotate_source(mol):
    out = Chem.Mol(mol)
    for a in out.GetAtoms():
        a.SetIntProp("_audit_source_atom", a.GetIdx())
    return out

def mapped_centers(mol):
    return {mol.GetAtomWithIdx(i).GetIntProp("_audit_source_atom"): s for i, s in Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            if mol.GetAtomWithIdx(i).HasProp("_audit_source_atom")}

def stages_and_variants(source):
    stages = [("source", Chem.Mol(source), "not_run")]
    clean = rdMolStandardize.Cleanup(source)
    stages.append(("cleanup", clean, "not_run"))
    parent = rdMolStandardize.FragmentParent(clean, skipStandardize=True)
    stages.append(("fragment_parent", parent, "not_run"))
    uncharged = core.UNCHARGER.uncharge(parent)
    stages.append(("uncharged", uncharged, "not_run"))
    old_enum = rdMolStandardize.TautomerEnumerator()
    old_t = old_enum.Enumerate(uncharged)
    old = old_enum.PickCanonical(old_t)
    Chem.AssignStereochemistry(old, cleanIt=True, force=True)
    stages.append(("frozen_old", old, str(old_t.status)))
    preserve_enum = rdMolStandardize.TautomerEnumerator()
    preserve_enum.SetRemoveSp3Stereo(False)
    preserve_enum.SetRemoveBondStereo(False)
    preserved_t = preserve_enum.Enumerate(uncharged)
    preserved = preserve_enum.PickCanonical(preserved_t)
    Chem.AssignStereochemistry(preserved, cleanIt=True, force=True)
    stages.append(("stereo_preserved_tautomer", preserved, str(preserved_t.status)))
    stripped = Chem.Mol(source)
    Chem.RemoveStereochemistry(stripped)
    variants = {"source": Chem.Mol(source), "no_tautomer": uncharged,
                "stereo_preserved_tautomer": preserved, "frozen_old": old,
                "source_stereo_stripped": stripped}
    status = {"frozen_old": str(old_t.status), "stereo_preserved_tautomer": str(preserved_t.status)}
    return stages, variants, status

def feature(package, mol):
    fp = core.morgan_array(mol).reshape(1, -1)
    raw_desc = core.descriptor_array(mol).reshape(1, -1)
    ft = package["feature_type"]
    if ft == "morgan":
        return fp
    if ft == "physchem":
        return package["descriptor_scaler"].transform(raw_desc)
    if ft in ("morgan_physchem", "morgan_plus_physchem"):
        return np.hstack([fp, package["descriptor_scaler"].transform(raw_desc)])
    raise ValueError("unsupported feature_type " + str(ft))

def score(package, mol, primary=False):
    x = feature(package, mol)
    if not np.isfinite(x).all():
        raise ValueError("non_finite_estimator_features")
    raw = float(package["estimator"].predict_proba(x)[0, 1])
    transformed = raw
    if primary and package.get("calibrator") is not None:
        clipped = float(np.clip(raw, 1e-6, 1-1e-6))
        logit = math.log(clipped / (1-clipped))
        transformed = float(package["calibrator"].predict_proba(np.asarray([[logit]]))[0, 1])
    return raw, transformed, x

def main():
    if (HERE / "run_manifest.json").exists():
        raise RuntimeError("Existing result bundle protected; create a fresh directory for another run")
    frozen = [BASE / x for x in ("entities_scored_or_hold.csv", "unique_structure_predictions.csv", "parent_release_score_contrasts.csv", "pair_release_stage_classification.csv")]
    frozen += [PRIMARY, DIRECT, Path(core.__file__), Path(__file__), HERE / "PROTOCOL.md"]
    before = {str(p.relative_to(ROOT)): digest(p) for p in frozen}
    (HERE / "prerun_freeze.json").write_text(json.dumps({"frozen_utc": datetime.now(timezone.utc).isoformat(), "input_sha256": before}, indent=2, ensure_ascii=False) + "\n")
    captured = []
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        primary, direct = joblib.load(PRIMARY), joblib.load(DIRECT)
        entities = read_csv(BASE / "entities_scored_or_hold.csv")
        old = {r["unique_structure_id"]: r for r in read_csv(BASE / "unique_structure_predictions.csv")}
        by_sid = {}
        for r in entities:
            if r.get("score_status") == "scored_frozen_model":
                by_sid.setdefault(r["unique_structure_id"], r)
        assert len(by_sid) == len(old) == 13
        stage_rows, atom_rows, score_rows, reproduction = [], [], [], []
        scores, mols_by_sid = {}, {}
        for sid, entity in by_sid.items():
            source = Chem.MolFromSmiles(entity["source_input_smiles"])
            assert source is not None
            source = annotate_source(source)
            source_rec = mol_record(source)
            stages, variants, enum_status = stages_and_variants(source)
            source_centers = mapped_centers(source)
            source_atoms = {a.GetIntProp("_audit_source_atom"): a for a in source.GetAtoms()}
            for stage, mol, enum in stages:
                rec = mol_record(mol)
                stage_rows.append({"unique_structure_id": sid, "entity_id": entity["entity_id"], "entity_name": entity["entity_name"],
                                   "release_stage": entity["release_stage"], "stage": stage, **rec,
                                   "full_key_changed_from_source": rec["inchikey"] != source_rec["inchikey"],
                                   "connectivity_block_changed_from_source": rec["inchikey"].split("-")[0] != source_rec["inchikey"].split("-")[0],
                                   "enumeration_status": enum})
                centers = mapped_centers(mol)
                mapped_atoms = {a.GetIntProp("_audit_source_atom"): a for a in mol.GetAtoms() if a.HasProp("_audit_source_atom")}
                assert len(mapped_atoms) == mol.GetNumAtoms()
                for src_idx, src_cip in source_centers.items():
                    a = mapped_atoms.get(src_idx)
                    stage_cip = centers.get(src_idx, "not_tetrahedral_center")
                    atom_rows.append({"unique_structure_id": sid, "entity_name": entity["entity_name"], "stage": stage,
                                      "source_atom_index_zero_based": src_idx, "source_assigned_cip": src_cip,
                                      "stage_cip": stage_cip,
                                      "source_atom_chiral_tag": str(source_atoms[src_idx].GetChiralTag()),
                                      "stage_atom_chiral_tag": str(a.GetChiralTag()) if a else "atom_removed",
                                      "stage_hybridization": str(a.GetHybridization()) if a else "atom_removed",
                                      "assigned_stereo_lost": src_cip != "?" and stage_cip in ("?", "not_tetrahedral_center"),
                                      "loss_type": "not_lost" if src_cip == "?" or stage_cip not in ("?", "not_tetrahedral_center") else "center_still_sp3_unassigned" if a and a.GetHybridization() == Chem.HybridizationType.SP3 else "chosen_tautomer_no_longer_sp3_center"})
            mols_by_sid[sid] = variants
            src_fp = core.morgan_array(source)
            src_desc = core.descriptor_array(source)
            old_fp = core.morgan_array(variants["frozen_old"])
            for variant, mol in variants.items():
                rec = mol_record(mol)
                hold = ""
                if variant in enum_status and enum_status[variant] != "Completed":
                    hold = "tautomer_enumeration_incomplete"
                if rec["formula"] != source_rec["formula"] or rec["formal_charge"] != source_rec["formal_charge"]:
                    hold = "formula_or_charge_changed_variant_not_like_for_like"
                desc = core.descriptor_array(mol)
                fp = core.morgan_array(mol)
                row = {"unique_structure_id": sid, "entity_id": entity["entity_id"], "entity_name": entity["entity_name"],
                       "release_stage": entity["release_stage"], "source_identity_limitation": entity["identity_limitation"],
                       "variant": variant, "diagnostic_status": "HOLD" if hold else "computed_representation_sensitivity_only",
                       "hold_reason": hold, **rec,
                       "morgan_bit_differences_vs_source": int(np.count_nonzero(fp != src_fp)),
                       "morgan_bit_differences_vs_frozen": int(np.count_nonzero(fp != old_fp)),
                       "descriptor_max_abs_difference_vs_source": float(np.abs(desc-src_desc).max()),
                       "descriptors_json": js(dict(zip(DESC_NAMES, map(float, desc)))),
                       "representation_contract": "training_frozen_preprocessing" if variant == "frozen_old" else "unvalidated_alternative_representation_diagnostic_only",
                       "is_independent_biological_validation": False,
                       "alternative_full_selector_decision": "not_computed"}
                for version, package, primary_flag in [("v0_1", primary, True), ("v0_3", direct, False)]:
                    row[version + "_raw_estimator_score"] = ""
                    row[version + "_risk_score"] = ""
                    if not hold:
                        raw_score, model_score, x = score(package, mol, primary_flag)
                        row[version + "_raw_estimator_score"] = raw_score
                        row[version + "_risk_score"] = model_score
                        row[version + "_binary_direction_at_0_5_diagnostic"] = int(model_score >= 0.5)
                        row[version + "_feature_array_sha256"] = hashlib.sha256(x.tobytes()).hexdigest()
                        row[version + "_max_training_tanimoto"] = core.max_training_similarity(mol, package["training_fingerprints"])
                        row[version + "_outside_frozen_similarity_q10"] = row[version + "_max_training_tanimoto"] < float(package["applicability_thresholds"]["self_nearest_q10"])
                score_rows.append(row)
                scores[(sid, variant)] = row
            # Check actual v0.1 frozen override, without invoking predict_smiles.
            frozen_primary_smiles = primary.get("training_structure_map", {}).get(source_rec["inchikey"], Chem.MolToSmiles(variants["frozen_old"], isomericSmiles=True))
            effective = Chem.MolFromSmiles(frozen_primary_smiles)
            for version, package, mol, prior_field in [("v0_1", primary, effective, "primary_calibrated_probability"),
                                                       ("v0_3", direct, variants["frozen_old"], "direct_panel_probability")]:
                _, val, _ = score(package, mol, version == "v0_1")
                previous = float(old[sid][prior_field])
                reproduction.append({"unique_structure_id": sid, "entity_name": entity["entity_name"], "version": version,
                                     "prior_score": previous, "recomputed_score": val,
                                     "absolute_difference": abs(val-previous), "passed_1e_8": abs(val-previous) <= 1e-8,
                                     "primary_training_structure_override": version == "v0_1" and source_rec["inchikey"] in primary.get("training_structure_map", {})})
        for row in score_rows:
            frozen_score = scores[(row["unique_structure_id"], "frozen_old")]
            source_score = scores[(row["unique_structure_id"], "source")]
            for v in ("v0_1", "v0_3"):
                a, b, c = [x[v+"_risk_score"] for x in (row, frozen_score, source_score)]
                row[v+"_score_minus_frozen"] = a-b if a != "" and b != "" else ""
                row[v+"_score_minus_source"] = a-c if a != "" and c != "" else ""
                old_field = "primary_calibrated_probability" if v == "v0_1" else "direct_panel_probability"
                historic = float(old[row["unique_structure_id"]][old_field])
                row[v+"_historical_frozen_score_not_newly_admitted"] = historic
                row[v+"_score_minus_historical_frozen_diagnostic"] = a-historic if a != "" else ""
        stages = {r["audit_record_id"]: r for r in read_csv(BASE / "pair_release_stage_classification.csv")}
        pair_rows, controls = [], []
        for pair in read_csv(BASE / "parent_release_score_contrasts.csv"):
            for variant in VARIANTS:
                psid, rsid = pair["parent_unique_structure_id"], pair["release_unique_structure_id"]
                a, b = scores.get((psid, variant)), scores.get((rsid, variant))
                stage = stages[pair["audit_record_id"]]
                row = {"audit_record_id": pair["audit_record_id"], "pair_type": pair["pair_type"], "parent_name": pair["free_name"],
                       "release_name": pair["released_name"], "release_stage_for_analysis": stage["release_stage_for_analysis"],
                       "may_count_as_terminal_chemical_change_pair": stage["may_count_as_terminal_chemical_change_pair"],
                       "variant": variant, "pair_status": "both_computed_diagnostic" if a and b and a["diagnostic_status"] != "HOLD" and b["diagnostic_status"] != "HOLD" else "representation_HOLD" if a and b else "identity_HOLD",
                       "parent_hold": pair["parent_identity_hold"], "release_hold": pair["release_identity_hold"],
                       "matched_M0_biological_truth": "not_available", "supports_accuracy_or_causal_claim": False}
                for v in ("v0_1", "v0_3"):
                    pa = a[v+"_risk_score"] if a else ""
                    rb = b[v+"_risk_score"] if b else ""
                    row["parent_"+v+"_risk_score"] = pa
                    row["release_"+v+"_risk_score"] = rb
                    row[v+"_release_minus_parent_score"] = rb-pa if rb != "" and pa != "" else ""
                    fpa, frb = scores.get((psid, "frozen_old")), scores.get((rsid, "frozen_old"))
                    row[v+"_pair_delta_change_vs_frozen"] = row[v+"_release_minus_parent_score"] - (frb[v+"_risk_score"]-fpa[v+"_risk_score"]) if row[v+"_release_minus_parent_score"] != "" and frb and fpa and frb[v+"_risk_score"] != "" and fpa[v+"_risk_score"] != "" else ""
                    historic = float(pair[v+"_release_minus_parent_score"]) if pair.get(v+"_release_minus_parent_score") else ""
                    row[v+"_historical_frozen_pair_delta_not_newly_admitted"] = historic
                    row[v+"_pair_delta_change_vs_historical_frozen_diagnostic"] = row[v+"_release_minus_parent_score"]-historic if row[v+"_release_minus_parent_score"] != "" and historic != "" else ""
                pair_rows.append(row)
                if pair["pair_type"] == "identity_preserving_control":
                    controls.append({"audit_record_id": pair["audit_record_id"], "variant": variant,
                                     "same_structure_id": psid == rsid,
                                     "same_scores": row["v0_1_release_minus_parent_score"] == 0 and row["v0_3_release_minus_parent_score"] == 0})
        stereo_checks = []
        for sid in by_sid:
            a, b = scores[(sid, "source")], scores[(sid, "source_stereo_stripped")]
            stereo_checks.append({"unique_structure_id": sid, "entity_name": a["entity_name"],
                                  "morgan_features_identical": a["morgan_bit_differences_vs_source"] == b["morgan_bit_differences_vs_source"] == 0,
                                  "descriptors_identical": b["descriptor_max_abs_difference_vs_source"] == 0,
                                  "v0_1_score_delta_strip_minus_source": b["v0_1_risk_score"]-a["v0_1_risk_score"],
                                  "v0_3_score_delta_strip_minus_source": b["v0_3_risk_score"]-a["v0_3_risk_score"]})
        captured = [{"category": w.category.__name__, "message": str(w.message), "filename": w.filename, "lineno": w.lineno} for w in ws]
    after = {str(p.relative_to(ROOT)): digest(p) for p in frozen}
    checks = {"frozen_inputs_unchanged": before == after, "prior_scores_reproduced": all(r["passed_1e_8"] for r in reproduction),
              "same_identity_controls_all_equal": all(r["same_structure_id"] and r["same_scores"] for r in controls),
              "stereo_only_features_scores_invariant": all(r["morgan_features_identical"] and r["descriptors_identical"] and r["v0_1_score_delta_strip_minus_source"] == r["v0_3_score_delta_strip_minus_source"] == 0 for r in stereo_checks)}
    write_csv("stage_structure_audit.csv", stage_rows)
    write_csv("source_atom_stereo_audit.csv", atom_rows)
    write_csv("variant_scores.csv", score_rows)
    write_csv("parent_release_variant_contrasts.csv", pair_rows)
    write_csv("same_identity_control_checks.csv", controls)
    write_csv("stereo_only_negative_control.csv", stereo_checks)
    write_csv("frozen_score_reproduction.csv", reproduction)
    summary = {"status": "completed" if all(checks.values()) else "verification_failed", "completed_utc": datetime.now(timezone.utc).isoformat(),
               "runtime": {"executable": sys.executable, "python": platform.python_version(), "rdkit": rdkit.__version__, "sklearn": sklearn.__version__, "numpy": np.__version__, "joblib": joblib.__version__},
               "models": {"v0_1_feature_type": primary["feature_type"], "v0_3_feature_type": direct["feature_type"], "morgan_include_chirality": False,
                          "default_remove_sp3_stereo": core.TAUTOMER_ENUMERATOR.GetRemoveSp3Stereo(), "default_remove_bond_stereo": core.TAUTOMER_ENUMERATOR.GetRemoveBondStereo()},
               "unique_source_structures": len(by_sid), "scored_variant_rows": len(score_rows), "identity_control_count": len(controls),
               "input_sha256": before, "checks": checks, "captured_python_warnings": captured,
               "rdkit_warning_policy": "core disables rdApp.warning; stderr retained separately; enumeration status recorded and incomplete variants HOLD",
               "no_new_truth_labels": True, "no_refit": True, "alternate_representation_validated": False,
               "source_vs_frozen_full_key_changed": sum(r["stage"] == "frozen_old" and r["full_key_changed_from_source"] for r in stage_rows),
               "source_vs_frozen_connectivity_block_changed": sum(r["stage"] == "frozen_old" and r["connectivity_block_changed_from_source"] for r in stage_rows),
               "stereo_loss_source_atom_count": sum(r["stage"] == "frozen_old" and r["assigned_stereo_lost"] for r in atom_rows)}
    for variant in VARIANTS:
        summary[variant+"_hold_count"] = sum(r["variant"] == variant and r["diagnostic_status"] == "HOLD" for r in score_rows)
    (HERE / "run_manifest.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": summary["status"], "checks": checks, "unique_structures": len(by_sid), "full_key_changes": summary["source_vs_frozen_full_key_changed"], "connectivity_changes": summary["source_vs_frozen_connectivity_block_changed"]}, ensure_ascii=False))
    if not all(checks.values()):
        raise RuntimeError("meaningful verification checks failed; inspect bundle")

if __name__ == "__main__":
    main()
