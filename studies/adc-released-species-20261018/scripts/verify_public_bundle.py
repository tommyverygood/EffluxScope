#!/usr/bin/env python3
"""Offline verification after extracting the public release to a fresh directory.

Reads frozen scientific files only. Every generated file goes to --output-dir,
which must not already exist. Does not fit models, run bootstrap, or fetch data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
import traceback
from collections import Counter
from pathlib import Path


BIOGEN = "outputs/20260810_soskuti_identity_double_review_v1/reanalysis/external/analysis/abcb1_quantity_biogen_analysis_rows_identity_audited_v2.csv"
ADC = "outputs/20261012_adc_completion/scoring/species_predictions.csv"
SCOPE = "outputs/20261012_adc_scope_revision"
BENCHMARK = "outputs/20261009_adc_pair_experiments/benchmark"
REGIONS = ("MW_ge600", "TPSA_ge140", "amide_ge4_MW_ge500")


def read(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def truth(value):
    return str(value).lower() in {"yes", "true", "1"}


def close(observed, expected, tolerance=1e-10):
    assert math.isfinite(float(observed)) and math.isfinite(float(expected))
    assert math.isclose(float(observed), float(expected), rel_tol=0, abs_tol=tolerance), (observed, expected)


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields)
        writer.writeheader()
        writer.writerows(rows)


def verify_pairs_and_asap_aggregation(root, output):
    pair_dir = root / "outputs/20261009_adc_pair_experiments/pairs/results_curated_v2"
    pairs = read(pair_dir / "parent_release_score_contrasts.csv")
    stages = {row["audit_record_id"]: row for row in read(pair_dir / "pair_release_stage_classification.csv")}
    assert len(pairs) == len(stages) == 14
    assert {row["audit_record_id"] for row in pairs} == set(stages)
    changed = [row for row in pairs if row["pair_type"] == "chemical_changed_parent_to_catabolite"
               and row["parent_score_status"] == row["release_score_status"] == "scored_frozen_model"]
    controls = [row for row in pairs if row["pair_type"] == "identity_preserving_control"]
    assert len(changed) == 5 and {row["skeleton_family"] for row in changed} == {"auristatin", "maytansinoid"}
    assert len(controls) == 5
    results = []
    for row in pairs:
        stage = stages[row["audit_record_id"]]
        assert stage["matched_M0_biological_truth"] == "no"
        both_scored = row["parent_score_status"] == row["release_score_status"] == "scored_frozen_model"
        assert both_scored == truth(stage["both_sides_computationally_scored"])
        if row in changed:
            assert "abstain" in (row["parent_full_decision"], row["release_full_decision"])
        if row in changed or row in controls:
            for channel in ("v0_1", "v0_3"):
                parent = row[f"parent_{channel}_risk_score"]
                release = row[f"release_{channel}_risk_score"]
                delta = row[f"{channel}_release_minus_parent_score"]
                # Empty or nonnumeric values raise; never replace missing evidence by zero.
                close(float(release) - float(parent), delta, 2e-8)
                if row in controls:
                    close(parent, release)
                    close(delta, 0)
            if row in controls:
                assert row["identity_control_exact_score_equality"] == "pass"
        else:
            assert not both_scored
            assert row["v0_1_release_minus_parent_score"] == row["v0_3_release_minus_parent_score"] == ""
        results.append({"audit_record_id": row["audit_record_id"], "pair_type": row["pair_type"],
                        "family": row["skeleton_family"], "both_sides_scored": both_scored,
                        "changed_computational_pair": row in changed, "identity_control": row in controls,
                        "matched_transport_truth_available": False})
    write_csv(output / "pair_inventory_reproduction.csv", results)
    asap_dir = root / "outputs/20261009_adc_followup_validation/asap_secondary"
    observations = read(root / "数据/processed/abcb1_selective_v0_4_asap2025_posthoc_predictions_v1.csv")
    ledger = read(asap_dir / "compound_aggregation_ledger.csv")
    identities = read(asap_dir / "repetition_audit/unique_compound_identity_and_aggregation.csv")
    assert len(observations) == len({row["observation_id"] for row in observations}) == 62
    counts = Counter(row["standardized_inchikey"] for row in observations)
    assert len(counts) == len(ledger) == len(identities) == 56
    assert set(counts) == {row["id"] for row in ledger} == {row["standardized_inchikey"] for row in identities}
    assert sum(int(row["n_observation_rows"]) for row in ledger) == sum(int(row["observation_count"]) for row in identities) == 62
    for row in ledger:
        assert row["key"] == row["id"] and int(row["n_observation_rows"]) == counts[row["id"]]
    for row in identities:
        key = row["standardized_inchikey"]
        ids = json.loads(row["observation_ids"])
        expected_ids = {item["observation_id"] for item in observations if item["standardized_inchikey"] == key}
        assert int(row["observation_count"]) == len(ids) == counts[key]
        assert len(ids) == len(set(ids)) and set(ids) == expected_ids
    return {"status": "PASS", "candidate_inventory_rows": 14, "changed_computational_pairs": 5,
            "changed_pair_families": ["auristatin", "maytansinoid"], "identity_preserving_controls": 5,
            "changed_pairs_with_at_least_one_abstention": 5, "matched_M0_transport_truth_pairs": 0,
            "pair_delta_definition": "release_score_minus_parent_score; not delta_ER",
            "ASAP_observation_rows": 62, "ASAP_aggregated_unique_compounds": 56,
            "ASAP_aggregation_counts_and_observation_ids_match": True,
            "ASAP_replicate_independence_established": False}


def verify_matched_count(root, output):
    expected = read(root / BENCHMARK / "metrics.csv")
    sources = {
        "identity_audited_v2": root / BIOGEN,
        "historical_v1_sensitivity": root / "数据/processed/abcb1_quantity_biogen_analysis_rows_v1.csv",
    }
    results = []
    for version, path in sources.items():
        rows = sorted((row for row in read(path) if row["identity_clean"] == "yes" and row["model_usable"] == "yes"), key=lambda row: row["external_row_id"])
        assert len(rows) == len({row["connectivity_block"] for row in rows}) == 2624
        full = [i for i, row in enumerate(rows) if row["accepted_prediction"] == "yes"]
        no_local = [i for i, row in enumerate(rows) if not (set(json.loads(row["abstention_reasons_json"])) - {"primary_local_neighborhood_strongly_opposite", "direct_local_neighborhood_strongly_opposite"})]
        assert len(full) == 453 and len(no_local) == 551
        for policy, pool, score_field in [
            ("v01_margin", range(len(rows)), "primary_calibrated_probability"),
            ("v03_margin", range(len(rows)), "direct_panel_probability"),
            ("full_v04_pool_primary_margin", full, "primary_calibrated_probability"),
            ("no_local_pool_primary_margin", no_local, "primary_calibrated_probability"),
        ]:
            selected = sorted(pool, key=lambda i: (-abs(float(rows[i][score_field]) - .5), rows[i]["external_row_id"]))[:453]
            tp = tn = fp = fn = 0
            for index in selected:
                positive = int(rows[index]["frozen_er2_binary_label"]) == 1
                predicted = float(rows[index][score_field]) >= .5
                tp += positive and predicted
                tn += not positive and not predicted
                fp += not positive and predicted
                fn += positive and not predicted
            sensitivity, specificity = tp / (tp + fn), tn / (tn + fp)
            observed = {"selected_n": len(selected), "selected_positive_n": tp + fn, "selected_negative_n": tn + fp,
                        "predicted_high_n": tp + fp, "predicted_low_n": tn + fn, "tp": tp, "tn": tn, "fp": fp, "fn": fn,
                        "accuracy": (tp + tn) / len(selected), "balanced_accuracy": (sensitivity + specificity) / 2,
                        "sensitivity": sensitivity, "specificity": specificity, "ppv": tp / (tp + fp), "npv": tn / (tn + fn),
                        "mcc": (tp * tn - fp * fn) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)),
                        "coverage": len(selected) / len(rows),
                        "positive_class_coverage": (tp + fn) / sum(int(row["frozen_er2_binary_label"]) for row in rows),
                        "negative_class_coverage": (tn + fp) / sum(1 - int(row["frozen_er2_binary_label"]) for row in rows)}
            target = next(row for row in expected if row["version"] == version and row["policy"] == policy and row["target"] == "fixed_full_coverage")
            for key, value in observed.items():
                close(value, target[key])
            results.append({"version": version, "policy": policy, **observed, "point_estimates_match": True, "bootstrap_recomputed": False})
    write_csv(output / "matched_count_reproduction.csv", results)
    return {"status": "PASS", "comparison_rows": len(results), "selected_per_policy": 453,
            "point_estimates_checked": True, "bootstrap_and_CI_recomputed": False}


def verify_features_and_gates(root, output):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    features = read(root / SCOPE / "feasibility/features.csv")
    assert len(features) == 2694
    assert Counter(row["source"] for row in features) == {"Biogen2624": 2624, "ASAP56": 56, "ADC14": 14}
    assert len({(row["source"], row["uid"]) for row in features}) == 2694
    pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    for row in features:
        mol = Chem.MolFromSmiles(row["standardized_isomeric_smiles"])
        assert mol is not None and Chem.MolToInchiKey(mol) == row["standardized_inchikey"]
        mw, tpsa = Descriptors.MolWt(mol), rdMolDescriptors.CalcTPSA(mol)
        amides = len(mol.GetSubstructMatches(pattern, uniquify=True))
        close(mw, row["molecular_weight"], .002)
        close(tpsa, row["tpsa"], 1e-5)
        assert amides == int(row["amide_count"])
        flags = {"MW_ge600": mw >= 600, "TPSA_ge140": tpsa >= 140, "amide_ge4_MW_ge500": amides >= 4 and mw >= 500}
        assert all(flags[key] == truth(row[key]) for key in REGIONS), row["uid"]
    count_record = json.loads((root / SCOPE / "feasibility/count_feasibility.json").read_text())
    reproduced = []
    for expected in count_record["counts"]:
        source, region = expected["source"], expected["region"]
        rows = [row for row in features if row["source"] == source and (region == "all_identity_clean" or truth(row[region]))]
        accepted = [row for row in rows if truth(row["full_accepted"])]
        observed = {"unique_compound_n": len(rows), "positive_n": sum(row["er_ge2_weak_label"] == "1" for row in rows),
                    "negative_n": sum(row["er_ge2_weak_label"] == "0" for row in rows),
                    "scaffold_n": len({row["murcko_scaffold"] for row in rows}), "full_accepted_n": len(accepted),
                    "full_accepted_positive_n": sum(row["er_ge2_weak_label"] == "1" for row in accepted),
                    "full_accepted_negative_n": sum(row["er_ge2_weak_label"] == "0" for row in accepted),
                    "full_accepted_scaffold_n": len({row["murcko_scaffold"] for row in accepted}),
                    "adc14_region_n": sum(row["source"] == "ADC14" and (region == "all_identity_clean" or truth(row[region])) for row in features)}
        assert all(observed[key] == expected[key] for key in observed), (source, region)
        reproduced.append({"source": source, "region": region, **observed})
    protocol = json.loads((root / SCOPE / "protocol.json").read_text())
    status = json.loads((root / SCOPE / "regional_evaluation_status.json").read_text())
    assert len(status["source_region_rows"]) == 6 and status["qualifying_regions"] == 0
    null_checked = 0
    for row in status["source_region_rows"]:
        counts = next(item for item in reproduced if item["source"] == row["source"] and item["region"] == row["region"])
        population, selected = protocol["population_gate"], protocol["matched_count_gate"]
        tests = {"population_n": (counts["unique_compound_n"], population["minimum_n"]),
                 "population_positive_n": (counts["positive_n"], population["minimum_each_measured_class"]),
                 "population_negative_n": (counts["negative_n"], population["minimum_each_measured_class"]),
                 "population_scaffold_n": (counts["scaffold_n"], population["minimum_scaffolds"]),
                 "full_selected_n": (counts["full_accepted_n"], selected["minimum_k"]),
                 "full_selected_positive_n": (counts["full_accepted_positive_n"], selected["minimum_each_selected_measured_class"]),
                 "full_selected_negative_n": (counts["full_accepted_negative_n"], selected["minimum_each_selected_measured_class"]),
                 "full_selected_scaffold_n": (counts["full_accepted_scaffold_n"], selected["minimum_selected_scaffolds"])}
        failures = {key: {"observed": value, "required": minimum} for key, (value, minimum) in tests.items() if value < minimum}
        assert failures and failures == row["failed_report_gates"]
        assert row["regional_comparison_status"] == "not_evaluable" and row["bootstrap_executed"] is False
        for field in ("BAcc", "MCC", "paired_delta_BAcc", "paired_CI95"):
            assert row[field] is None
            null_checked += 1
    write_csv(output / "chemical_region_count_reproduction.csv", reproduced)
    return {"status": "PASS", "individual_feature_rows": len(features), "count_rows": len(reproduced),
            "regional_rows": 6, "explicit_null_values_checked": null_checked, "performance_or_bootstrap_computed": False}


def verify_adc(root, output):
    sys.path.insert(0, str(root / "模型/ABCB1_selective_screener_v0_4_dev"))
    from predict_abcb1_selective import load_packages, predict_one, flatten_batch_result
    model_paths = [root / "模型/ABCB1_selective_screener_v0_4_dev/ABCB1_selective_uncertainty_v0_4_dev.joblib",
                   root / "模型/ABCB1_screener_v0_1/ABCB1_screener_v0_1.joblib",
                   root / "模型/ABCB1_assay_aware_candidate_v0_3_dev/ABCB1_assay_aware_candidate_v0_3_dev.joblib"]
    packages = load_packages(*model_paths)
    expected_rows = read(root / ADC)
    assert len(expected_rows) == len({row["entity_id"] for row in expected_rows}) == 14
    numbers = ("primary_calibrated_probability", "direct_panel_probability", "bootstrap_probability_p10", "bootstrap_probability_median", "bootstrap_probability_p90", "physchem_max_abs_robust_z", "primary_max_tanimoto", "direct_max_tanimoto")
    texts = ("standardization_status", "standardized_inchikey", "risk_call", "accepted_prediction", "applicability_status", "uncertainty_status")
    results = []
    for index, expected in enumerate(expected_rows, 1):
        print(f"ADC {index}/14: {expected['entity_id']}", flush=True)
        result = flatten_batch_result(predict_one(expected["source_smiles"], *packages,
                     chemical_object_type=expected["chemical_object_type"], requested_context=expected["requested_context"]))
        for field in numbers:
            close(result[field], expected[field], 5e-7)
        for field in texts:
            assert result[field] == expected[field], (expected["entity_id"], field)
        assert json.loads(result["abstention_reasons_json"]) == json.loads(expected["abstention_reasons_json"])
        results.append({"entity_id": expected["entity_id"], "entity_name": expected["entity_name"], **result})
    assert Counter(row["risk_call"] for row in results) == {"high_ABCB1_efflux_risk": 5, "abstain": 9}
    write_csv(output / "adc14_reproduction.csv", results)
    return {"status": "PASS", "species_n": len(results), "numeric_fields_per_species": len(numbers),
            "categorical_fields_per_species": len(texts) + 1, "numeric_tolerance_absolute": 5e-7,
            "model_sha256": {str(path.relative_to(root)): sha(path) for path in model_paths}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="Extracted public bundle root; preserve its relative directory tree")
    parser.add_argument("--output-dir", required=True, type=Path, help="Fresh verification directory; may be outside the extracted tree")
    parser.add_argument("--skip-inference", action="store_true", help="Only verify published scores/counts/features, without loading frozen models")
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output_dir.resolve()
    if not root.is_dir():
        parser.error("--root must be an existing extracted bundle directory")
    if output.exists():
        parser.error("--output-dir must not already exist")
    output.mkdir(parents=True)
    protected = [path for path in root.rglob("*") if path.is_file() and not path.is_relative_to(output) and "__pycache__" not in path.parts and ".git" not in path.parts]
    before = {str(path.relative_to(root)): sha(path) for path in protected}
    receipt = {"status": "RUNNING", "python": platform.python_version(), "bundle_root": str(root), "output_dir": str(output),
               "no_training": True, "no_network": True, "no_bootstrap_recomputed": True}
    failed = False
    try:
        print("Checking five changed computational pairs and ASAP 62-to-56 aggregation", flush=True)
        receipt["pairs_and_asap_aggregation"] = verify_pairs_and_asap_aggregation(root, output)
        print("Checking matched-count point estimates", flush=True)
        receipt["matched_count"] = verify_matched_count(root, output)
        print("Checking 2694 molecular-feature rows and reporting gates", flush=True)
        receipt["features_and_gates"] = verify_features_and_gates(root, output)
        receipt["adc14"] = {"status": "SKIPPED_by_explicit_argument"} if args.skip_inference else verify_adc(root, output)
        changed = [name for name, digest in before.items() if not (root / name).is_file() or sha(root / name) != digest]
        assert not changed, changed
        receipt["original_bundle_files_unchanged"] = True
        receipt["protected_file_n"] = len(before)
        receipt["status"] = "PASS" if not args.skip_inference else "PASS_counts_and_features_only_inference_skipped"
    except Exception as exc:
        failed = True
        receipt["status"] = "FAIL"
        receipt["error"] = f"{type(exc).__name__}: {exc}"
        receipt["traceback"] = traceback.format_exc()
        print(receipt["traceback"], file=sys.stderr)
    (output / "verification.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": receipt["status"], "receipt": str(output / "verification.json")}, ensure_ascii=False), flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
