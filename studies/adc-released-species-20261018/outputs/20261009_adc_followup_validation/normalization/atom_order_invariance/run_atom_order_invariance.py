#!/usr/bin/env python3
from __future__ import annotations
import csv
import hashlib
import importlib.util
import json
import platform
import sys
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
NORMALIZATION = HERE.parent
ROOT = HERE.parents[3]
spec = importlib.util.spec_from_file_location("normalization_validation", NORMALIZATION / "run_normalization_validation.py")
nv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nv)
Chem, core, np = nv.Chem, nv.core, nv.np
SEED = 20261009
N_RANDOM = 16

def write_csv(name, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (HERE / name).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def main():
    if (HERE / "run_manifest.json").exists():
        raise RuntimeError("Existing bundle protected")
    inputs = [nv.BASE / name for name in ["entities_scored_or_hold.csv", "unique_structure_predictions.csv", "parent_release_score_contrasts.csv", "pair_release_stage_classification.csv"]]
    inputs += [nv.PRIMARY, nv.DIRECT, Path(core.__file__), Path(nv.__file__), Path(__file__), HERE / "PROTOCOL.md"]
    before = {str(p.relative_to(ROOT)): nv.digest(p) for p in inputs}
    (HERE / "prerun_freeze.json").write_text(json.dumps({"frozen_utc": datetime.now(timezone.utc).isoformat(), "seed": SEED, "input_sha256": before}, ensure_ascii=False, indent=2) + "\n")
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        primary, direct = nv.joblib.load(nv.PRIMARY), nv.joblib.load(nv.DIRECT)
        unique = {}
        for r in nv.read_csv(nv.BASE / "entities_scored_or_hold.csv"):
            if r["score_status"] == "scored_frozen_model":
                unique.setdefault(r["unique_structure_id"], r)
        prior = {r["unique_structure_id"]: r for r in nv.read_csv(nv.BASE / "unique_structure_predictions.csv")}
        rows, summaries, repro = [], [], []
        for entity_index, (sid, ent) in enumerate(unique.items()):
            source = Chem.MolFromSmiles(ent["source_input_smiles"])
            assert source is not None
            expected_key = Chem.MolToInchiKey(source)
            inputs_smiles = [ent["source_input_smiles"]] + list(Chem.MolToRandomSmilesVect(source, N_RANDOM, randomSeed=SEED+entity_index, isomericSmiles=True))
            entity_rows = []
            for rep_index, smiles in enumerate(inputs_smiles):
                raw = Chem.MolFromSmiles(smiles)
                raw_key = Chem.MolToInchiKey(raw) if raw is not None else ""
                row = {"unique_structure_id": sid, "entity_id": ent["entity_id"], "entity_name": ent["entity_name"],
                       "release_stage": ent["release_stage"], "representation_index": rep_index,
                       "input_kind": "original" if rep_index == 0 else "random_isomeric_smiles", "random_seed": SEED+entity_index,
                       "input_smiles": smiles, "source_raw_inchikey": expected_key, "input_raw_inchikey": raw_key,
                       "same_raw_identity": raw_key == expected_key, "status": "raw_identity_failure_excluded"}
                if raw_key != expected_key:
                    rows.append(row)
                    entity_rows.append(row)
                    continue
                normalized = core.standardize_smiles(smiles)
                row.update({k: normalized[k] for k in ["standardization_status", "standardized_isomeric_smiles", "standardized_inchikey", "tautomer_status", "fragment_removed_flag", "standardization_error"]})
                if normalized["standardization_status"] != "ok":
                    row["status"] = "standardization_failure_excluded"
                else:
                    row["status"] = "diagnostic_computed_cap_hit_not_newly_admitted" if normalized["tautomer_status"] != "Completed" else "diagnostic_computed"
                    mol = Chem.MolFromSmiles(normalized["standardized_isomeric_smiles"])
                    for version, package, primary_flag in [("v0_1", primary, True), ("v0_3", direct, False)]:
                        raw_score, score, x = nv.score(package, mol, primary_flag)
                        row[version+"_raw_estimator_score"] = raw_score
                        row[version+"_risk_score"] = score
                        row[version+"_feature_sha256"] = hashlib.sha256(x.tobytes()).hexdigest()
                        row[version+"_binary_direction_at_0_5_diagnostic"] = int(score >= 0.5)
                        if version == "v0_1":
                            row["primary_margin_category_diagnostic"] = "high" if score >= 0.70 else "low" if score <= 0.35 else "gray"
                        if rep_index == 0:
                            old_field = "primary_calibrated_probability" if version == "v0_1" else "direct_panel_probability"
                            old_score = float(prior[sid][old_field])
                            repro.append({"unique_structure_id": sid, "entity_name": ent["entity_name"], "version": version,
                                          "historical_score": old_score, "recomputed_original_score": score,
                                          "passed_1e_8": abs(score-old_score) <= 1e-8, "absolute_difference": abs(score-old_score)})
                rows.append(row)
                entity_rows.append(row)
            successful = [r for r in entity_rows if "v0_1_risk_score" in r]
            summary = {"unique_structure_id": sid, "entity_id": ent["entity_id"], "entity_name": ent["entity_name"], "release_stage": ent["release_stage"],
                       "source_raw_inchikey": expected_key, "represented_input_rows": len(entity_rows), "unique_input_smiles_count": len(set(inputs_smiles)),
                       "raw_identity_failures": sum(not r["same_raw_identity"] for r in entity_rows), "successfully_diagnosed_rows": len(successful),
                       "unique_standardized_key_count": len({r["standardized_inchikey"] for r in successful}),
                       "unique_standardized_keys_json": nv.js(sorted({r["standardized_inchikey"] for r in successful})),
                       "tautomer_status_counts_json": nv.js(dict(Counter(r["tautomer_status"] for r in successful))),
                       "any_enumeration_cap_hit": any(r["tautomer_status"] != "Completed" for r in successful)}
            for version in ["v0_1", "v0_3"]:
                vals = [r[version+"_risk_score"] for r in successful]
                summary[version+"_score_min"] = min(vals) if vals else ""
                summary[version+"_score_max"] = max(vals) if vals else ""
                summary[version+"_score_range"] = max(vals)-min(vals) if vals else ""
                summary[version+"_unique_rounded8_score_count"] = len({round(v, 8) for v in vals})
                summary[version+"_unique_feature_count"] = len({r[version+"_feature_sha256"] for r in successful})
                summary[version+"_unique_direction_count"] = len({r[version+"_binary_direction_at_0_5_diagnostic"] for r in successful})
            summary["primary_margin_unique_categories_json"] = nv.js(sorted({r["primary_margin_category_diagnostic"] for r in successful}))
            summary["primary_margin_category_count"] = len({r["primary_margin_category_diagnostic"] for r in successful})
            summary["key_varies_score_stays_same_rounded8"] = summary["unique_standardized_key_count"] > 1 and summary["v0_1_unique_rounded8_score_count"] == summary["v0_3_unique_rounded8_score_count"] == 1
            summaries.append(summary)
            print(nv.js({"entity": ent["entity_name"], "keys": summary["unique_standardized_key_count"], "v01_range": summary["v0_1_score_range"], "v03_range": summary["v0_3_score_range"], "cap_hit": summary["any_enumeration_cap_hit"]}), flush=True)
        warning_rows = [{"category": w.category.__name__, "message": str(w.message)} for w in ws]
    stage_map = {r["audit_record_id"]: r for r in nv.read_csv(nv.BASE / "pair_release_stage_classification.csv")}
    by_sid = {r["unique_structure_id"]: r for r in summaries}
    pair_ranges = []
    for p in nv.read_csv(nv.BASE / "parent_release_score_contrasts.csv"):
        a, b = by_sid.get(p["parent_unique_structure_id"]), by_sid.get(p["release_unique_structure_id"])
        row = {"audit_record_id": p["audit_record_id"], "pair_type": p["pair_type"], "parent_name": p["free_name"], "release_name": p["released_name"],
               "release_stage_for_analysis": stage_map[p["audit_record_id"]]["release_stage_for_analysis"],
               "may_count_as_terminal_chemical_change_pair": stage_map[p["audit_record_id"]]["may_count_as_terminal_chemical_change_pair"],
               "diagnostic_status": "both_computed_fixed_identity_variation" if a and b else "identity_HOLD",
               "range_is_biological_CI": False, "independent_biological_validation": False}
        for v in ["v0_1", "v0_3"]:
            if a and b:
                lo, hi = b[v+"_score_min"]-a[v+"_score_max"], b[v+"_score_max"]-a[v+"_score_min"]
                row[v+"_release_minus_parent_lower"] = lo
                row[v+"_release_minus_parent_upper"] = hi
                row[v+"_delta_observed_direction"] = "positive_all_observed" if lo > 1e-8 else "negative_all_observed" if hi < -1e-8 else "zero_or_spans_zero"
                old_delta = float(p[v+"_release_minus_parent_score"])
                row[v+"_historical_pair_delta"] = old_delta
                row[v+"_historical_delta_within_observed_range_tolerance_1e_8"] = lo-1e-8 <= old_delta <= hi+1e-8
            else:
                row[v+"_release_minus_parent_lower"] = row[v+"_release_minus_parent_upper"] = ""
        pair_ranges.append(row)
    after = {str(p.relative_to(ROOT)): nv.digest(p) for p in inputs}
    checks = {"input_hashes_unchanged": before == after, "all_raw_identities_unchanged": all(r["same_raw_identity"] for r in rows),
              "all_221_rows_retained": len(rows) == 13*17, "all_original_scores_reproduced_1e_8": all(r["passed_1e_8"] for r in repro),
              "same_identity_controls_range_contains_zero": all(r[v+"_release_minus_parent_lower"] <= 0 <= r[v+"_release_minus_parent_upper"] for r in pair_ranges if r["pair_type"] == "identity_preserving_control" for v in ["v0_1", "v0_3"])}
    write_csv("representation_rows.csv", rows)
    write_csv("structure_invariance_summary.csv", summaries)
    write_csv("pair_score_ranges.csv", pair_ranges)
    write_csv("historical_score_reproduction.csv", repro)
    manifest = {"status": "completed" if all(checks.values()) else "verification_failed", "completed_utc": datetime.now(timezone.utc).isoformat(),
                "runtime": {"executable": sys.executable, "python": platform.python_version(), "rdkit": nv.rdkit.__version__, "numpy": np.__version__, "sklearn": nv.sklearn.__version__, "joblib": nv.joblib.__version__},
                "seed": SEED, "random_representations_per_structure": N_RANDOM, "source_structure_count": len(summaries), "input_rows": len(rows),
                "default_max_tautomers": core.TAUTOMER_ENUMERATOR.GetMaxTautomers(), "default_max_transforms": core.TAUTOMER_ENUMERATOR.GetMaxTransforms(),
                "input_sha256": before, "checks": checks, "python_warnings": warning_rows,
                "cap_hit_structure_count": sum(r["any_enumeration_cap_hit"] for r in summaries),
                "key_variable_structure_count": sum(r["unique_standardized_key_count"] > 1 for r in summaries),
                "score_variable_structure_count_rounded8": sum(r["v0_1_unique_rounded8_score_count"] > 1 or r["v0_3_unique_rounded8_score_count"] > 1 for r in summaries),
                "threshold_direction_variable_structure_count": sum(r["v0_1_unique_direction_count"] > 1 or r["v0_3_unique_direction_count"] > 1 for r in summaries),
                "primary_margin_category_variable_structure_count": sum(r["primary_margin_category_count"] > 1 for r in summaries),
                "finite_test_proves_tautomer_completeness": False, "variant_hold_policy_unchanged": True, "no_refit": True}
    (HERE / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(nv.js({"status": manifest["status"], "checks": checks, "key_variable": manifest["key_variable_structure_count"], "score_variable": manifest["score_variable_structure_count_rounded8"]}), flush=True)
    if not all(checks.values()):
        raise RuntimeError("meaningful validation checks failed")

if __name__ == "__main__":
    main()
