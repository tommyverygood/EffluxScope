#!/usr/bin/env python3
"""Independent graph-count and cross-agent count receipts, no performance metrics."""
import csv
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Fragments

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
FORMAL = HERE.parent / "protocol.json"
REVIEW = HERE.parent / "protocol_review/independent_region_feasibility.csv"

def read(path):
    with path.open(encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def literal_graph_amides(mol):
    # Literal equivalent of [CX3](=[OX1])[NX3], implemented without SMARTS.
    matches = set()
    for c in mol.GetAtoms():
        if c.GetAtomicNum() != 6 or c.GetIsAromatic() or c.GetTotalDegree() != 3:
            continue
        oxygens = [o for o in c.GetNeighbors() if o.GetAtomicNum() == 8 and not o.GetIsAromatic() and o.GetTotalDegree() == 1 and mol.GetBondBetweenAtoms(c.GetIdx(), o.GetIdx()).GetBondType() == Chem.BondType.DOUBLE]
        nitrogens = [n for n in c.GetNeighbors() if n.GetAtomicNum() == 7 and not n.GetIsAromatic() and n.GetTotalDegree() == 3 and mol.GetBondBetweenAtoms(c.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
        matches |= {(c.GetIdx(), o.GetIdx(), n.GetIdx()) for o in oxygens for n in nitrogens}
    return len(matches)

def dump(name, obj):
    (HERE / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")

rs = read(HERE / "features.csv")
manifest = json.loads((HERE / "completion.json").read_text())
feature_checks, fragment_definition_differences = [], []
for r in rs:
    mol = Chem.MolFromSmiles(r["standardized_isomeric_smiles"])
    mw = abs(Descriptors.MolWt(mol) - float(r["molecular_weight"])) < 1e-6
    tpsa = abs(rdMolDescriptors.CalcTPSA(mol) - float(r["tpsa"])) < 1e-6
    amide = literal_graph_amides(mol) == int(r["amide_count"])
    assert mw and tpsa and amide, r["uid"]
    feature_checks.append({"source": r["source"], "uid": r["uid"], "MW_match": mw, "TPSA_match": tpsa, "literal_graph_amide_match": amide})
    fragment = int(Fragments.fr_amide(mol))
    if fragment != int(r["amide_count"]):
        fragment_definition_differences.append({"source": r["source"], "uid": r["uid"], "formal_SMARTS_count": int(r["amide_count"]), "Fragments_fr_amide_count": fragment, "affects_amide4_MW500_membership": (fragment >= 4 and float(r["molecular_weight"]) >= 500) != (r["amide_ge4_MW_ge500"] == "True")})
assert not any(r["affects_amide4_MW500_membership"] for r in fragment_definition_differences)
assert all(sha(ROOT / p) == value for p, value in manifest["input_sha256"].items())

region_map = {"MW_ge_600": "MW_ge600", "TPSA_ge_140": "TPSA_ge140", "amide_ge_4_and_MW_ge_500": "amide_ge4_MW_ge500"}
comparison_checks = []
for independent in read(REVIEW):
    pop = [r for r in rs if r["source"] == independent["source"] and r[region_map[independent["region"]]] == "True"]
    k = sum(r["full_accepted"] == "True" for r in pop)
    if independent["policy"] == "population": selected = pop
    elif independent["policy"] == "full": selected = [r for r in pop if r["full_accepted"] == "True"]
    else:
        score_field = "primary_score" if independent["policy"] == "primary_margin" else "auxiliary_score"
        selected = sorted(pop, key=lambda r: (-abs(float(r[score_field]) - .5), r["uid"]))[:k]
    calculated = {"population_n": len(pop), "selected_n": len(selected),
                  "positive_n": sum(r["er_ge2_weak_label"] == "1" for r in selected),
                  "negative_n": sum(r["er_ge2_weak_label"] == "0" for r in selected),
                  "cluster_n": len({r["murcko_scaffold"] for r in selected}),
                  "coverage": len(selected) / len(pop) if pop else None}
    for field, value in calculated.items():
        expected = float(independent[field]) if independent[field] else None
        passed = expected is None and value is None or expected is not None and value is not None and abs(expected - value) < 1e-12
        comparison_checks.append({"source": independent["source"], "region": independent["region"], "policy": independent["policy"], "field": field, "this_run": value, "independent_review": expected, "passed": passed})
        assert passed, comparison_checks[-1]

protocol = json.loads(FORMAL.read_text())
summaries = json.loads((HERE / "count_feasibility.json").read_text())["counts"]
gate_receipt = []
for r in summaries:
    if r["region"] == "all_identity_clean": continue
    population_pass = r["unique_compound_n"] >= protocol["population_gate"]["minimum_n"] and min(r["positive_n"], r["negative_n"]) >= protocol["population_gate"]["minimum_each_measured_class"] and r["scaffold_n"] >= protocol["population_gate"]["minimum_scaffolds"]
    matched_pass = r["full_accepted_n"] >= protocol["matched_count_gate"]["minimum_k"] and min(r["full_accepted_positive_n"], r["full_accepted_negative_n"]) >= protocol["matched_count_gate"]["minimum_each_selected_measured_class"] and r["full_accepted_scaffold_n"] >= protocol["matched_count_gate"]["minimum_selected_scaffolds"]
    gate_receipt.append({"source": r["source"], "region": r["region"], "formal_population_gate_pass": population_pass, "formal_full_matched_gate_pass": matched_pass, "performance_comparison_allowed": population_pass and matched_pass})
assert not any(r["performance_comparison_allowed"] for r in gate_receipt)
dump("independent_feature_count_verification.json", {"status": "PASS", "feature_rows_checked": len(feature_checks), "MW_TPSA_recomputed_match_count": len(feature_checks), "literal_graph_formal_SMARTS_match_count": len(feature_checks), "Fragments_definition_differences": fragment_definition_differences, "no_difference_changes_region_membership": True, "all_frozen_source_hashes_reverified": True, "model_predictions_not_recomputed": True, "model_performance_not_computed": True})
dump("cross_agent_count_receipt.json", {"status": "PASS", "independent_review_sha256": sha(REVIEW), "cross_agent_rows_checked": len(read(REVIEW)), "numeric_field_checks": len(comparison_checks), "all_checks_pass": True, "checks": comparison_checks, "amide_definition_note": "Independent review uses Fragments.fr_amide, formal protocol uses [CX3](=[OX1])[NX3]; seven per-entity counts differ but no region membership differs"})
dump("formal_protocol_supersession.json", {"status": "RECORDED", "sole_performance_protocol": str(FORMAL.relative_to(ROOT)), "formal_protocol_sha256": sha(FORMAL), "superseded_count_only_practical_gate": "feasibility/protocol.json minimum_practical_gate_not_power_analysis (10 accepted/3 each class) is diagnostic only and has no authority for performance comparison", "formal_gate_receipt": gate_receipt, "three_regions_all_fail_formal_gate": True, "thresholds_unchanged": True, "no_performance_computed": True})
print(json.dumps({"features_verified": len(feature_checks), "cross_agent_rows": len(read(REVIEW)), "field_checks": len(comparison_checks), "formal_region_gate": gate_receipt, "motif_definition_differences": fragment_definition_differences}, ensure_ascii=False))
