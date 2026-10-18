#!/usr/bin/env python3
"""Count predefined chemical stress regions without computing model performance."""
from __future__ import annotations
import csv
import hashlib
import json
import platform
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import rdkit
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.warning")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INPUTS = {
    "biogen": ROOT / "outputs/20260810_soskuti_identity_double_review_v1/reanalysis/external/analysis/abcb1_quantity_biogen_analysis_rows_identity_audited_v2.csv",
    "asap56": ROOT / "outputs/20261009_adc_followup_validation/asap_secondary/compound_aggregation_ledger.csv",
    "asap_identity56": ROOT / "outputs/20261009_adc_followup_validation/asap_secondary/repetition_audit/unique_compound_identity_and_aggregation.csv",
    "asap62": ROOT / "数据/processed/abcb1_selective_v0_4_asap2025_posthoc_predictions_v1.csv",
    "adc14": ROOT / "outputs/20261012_adc_completion/scoring/species_predictions.csv",
    "primary443": ROOT / "数据/processed/abcb1_screener_v0_1_training_audit.csv",
    "direct84": ROOT / "数据/processed/soskuti_2024_mdr1_species_panel_identity_audited_v2.csv",
    "protocol": HERE / "protocol.json",
    "script": Path(__file__).resolve(),
}
REGIONS = ["all_identity_clean", "MW_ge600", "TPSA_ge140", "amide_ge4_MW_ge500"]
AMIDE = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")

def read(path):
    with path.open(encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def dump(name, obj):
    (HERE / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n")

def write(name, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (HERE / name).open("w", encoding="utf-8-sig", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields)
        w.writeheader(); w.writerows(rows)

def true(value):
    return str(value).lower() in {"true", "yes", "1"}

def make_feature(source, uid, r, source_smiles, label, er, scaffold, observations, source_key=""):
    descriptors = json.loads(r["physicochemical_descriptors_json"])
    mol = Chem.MolFromSmiles(r["standardized_isomeric_smiles"])
    if mol is None:
        raise ValueError(f"Cannot parse frozen model graph: {uid}")
    key = r["standardized_inchikey"]
    assert Chem.MolToInchiKey(mol) == key, (uid, key, Chem.MolToInchiKey(mol))
    source_mol = Chem.MolFromSmiles(source_smiles)
    assert source_mol is not None
    computed_source_key = Chem.MolToInchiKey(source_mol)
    if source_key:
        assert source_key == computed_source_key, (uid, source_key, computed_source_key)
    amides = len(mol.GetSubstructMatches(AMIDE, uniquify=True))
    mw, tpsa = float(descriptors["molecular_weight"]), float(descriptors["tpsa"])
    source_keys = {computed_source_key, key}
    conn = key.split("-")[0]
    feature = {
        "source": source, "uid": uid, "entity_name": r.get("entity_name", r.get("molecule_pref_name", "")),
        "standardized_inchikey": key, "connectivity_block": conn,
        "source_exact_inchikey_computed": computed_source_key,
        "standardized_isomeric_smiles": r["standardized_isomeric_smiles"],
        "source_smiles": source_smiles,
        "molecular_weight": mw, "tpsa": tpsa, "amide_count": amides,
        "murcko_scaffold": scaffold or "acyclic::" + conn,
        "observation_count": observations, "er_ge2_weak_label": label,
        "measured_er": er, "full_accepted": true(r["accepted_prediction"]),
        "primary_score": float(r["primary_calibrated_probability"]),
        "auxiliary_score": float(r["direct_panel_probability"]),
        "risk_call": r["risk_call"], "abstention_reasons_json": r["abstention_reasons_json"],
        "validation_role": "ADC_component_mapping_no_measured_truth" if source == "ADC14" else "retrospective_general_molecule_weak_ER_endpoint",
        "posthoc_aggregation": source == "ASAP56",
        "MW_ge600": mw >= 600, "TPSA_ge140": tpsa >= 140,
        "amide_ge4_MW_ge500": amides >= 4 and mw >= 500,
    }
    feature["development_standardized_exact_overlap"] = key in reference_keys
    feature["development_standardized_connectivity_overlap"] = conn in reference_conn
    feature["development_source_exact_overlap"] = bool(source_keys & reference_keys)
    feature["development_source_connectivity_overlap"] = bool({k.split("-")[0] for k in source_keys} & reference_conn)
    return feature

def counts(rows, source, region, adc):
    selected = [r for r in rows if r["source"] == source and (region == "all_identity_clean" or r[region])]
    accepted = [r for r in selected if r["full_accepted"]]
    out = {"source": source, "region": region, "unique_compound_n": len(selected),
           "positive_n": sum(r["er_ge2_weak_label"] == 1 for r in selected),
           "negative_n": sum(r["er_ge2_weak_label"] == 0 for r in selected),
           "scaffold_n": len({r["murcko_scaffold"] for r in selected}),
           "full_accepted_n": len(accepted),
           "full_accepted_positive_n": sum(r["er_ge2_weak_label"] == 1 for r in accepted),
           "full_accepted_negative_n": sum(r["er_ge2_weak_label"] == 0 for r in accepted),
           "full_accepted_scaffold_n": len({r["murcko_scaffold"] for r in accepted}),
           "adc14_region_n": 14 if region == "all_identity_clean" else sum(r[region] for r in adc),
           "adc14_region_fraction": 1.0 if region == "all_identity_clean" else sum(r[region] for r in adc) / 14,
           "missing_full_selected_class": not (any(r["er_ge2_weak_label"] == 1 for r in accepted) and any(r["er_ge2_weak_label"] == 0 for r in accepted)),
           "no_full_accepted": not accepted,
           "minimum_practical_count_gate_pass": False}
    gate = protocol["minimum_practical_gate_not_power_analysis"]
    mapped = {"population_n": "unique_compound_n", "positive_n": "positive_n", "negative_n": "negative_n",
              "scaffold_n": "scaffold_n", "full_accepted_n": "full_accepted_n", "full_accepted_positive_n": "full_accepted_positive_n",
              "full_accepted_negative_n": "full_accepted_negative_n"}
    out["minimum_practical_count_gate_pass"] = all(out[mapped[k]] >= v for k, v in gate.items())
    out["count_gate_shortfalls"] = {k: {"observed": out[mapped[k]], "required": v} for k, v in gate.items() if out[mapped[k]] < v}
    return out

if __name__ == "__main__":
    if (HERE / "completion.json").exists():
        raise RuntimeError("Do not overwrite completed feasibility run")
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in INPUTS.values()}
    protocol = json.loads(INPUTS["protocol"].read_text())
    dump("prerun_input_manifest.json", {"created_utc": datetime.now(timezone.utc).isoformat(), "sha256": hashes})
    primary = read(INPUTS["primary443"])
    direct = read(INPUTS["direct84"])
    assert len(primary) == 443 and len(direct) == 84
    reference_keys = {r["derived_standardized_inchikey"] for r in primary if r["derived_standardized_inchikey"]}
    reference_keys |= {r["standardized_inchikey"] for r in direct if r["standardized_inchikey"]}
    reference_conn = {k.split("-")[0] for k in reference_keys}
    features, excluded, audit = [], [], []
    biogen = read(INPUTS["biogen"])
    for r in biogen:
        if r["identity_clean"] != "yes" or r["model_usable"] != "yes":
            excluded.append({"source": "Biogen2624", "uid": r["external_row_id"], "reason": "historical_identity_clean_or_model_usable_exclusion", "standardized_inchikey": r["standardized_inchikey"]})
            continue
        feature = make_feature("Biogen2624", r["external_row_id"], r, r["input_smiles"], int(r["frozen_er2_binary_label"]), float(r["mdr1_mdck_er"]), r["bootstrap_cluster"], int(r["within_source_row_count"]))
        assert feature["er_ge2_weak_label"] == int(feature["measured_er"] >= 2)
        assert not feature["development_standardized_connectivity_overlap"]
        if feature["development_source_connectivity_overlap"]:
            excluded.append({"source": "Biogen2624", "uid": feature["uid"], "reason": "additional_source_or_standardized_reference_connectivity_overlap", "standardized_inchikey": feature["standardized_inchikey"]})
        else:
            features.append(feature)
    asap62 = read(INPUTS["asap62"])
    grouped = defaultdict(list)
    for r in asap62: grouped[r["standardized_inchikey"]].append(r)
    asap_identity = {r["standardized_inchikey"]: r for r in read(INPUTS["asap_identity56"])}
    asap56 = read(INPUTS["asap56"])
    assert len(asap62) == 62 and len(asap56) == len(grouped) == len(asap_identity) == 56
    for aggregate in asap56:
        key = aggregate["key"]
        rs, identity = grouped[key], asap_identity[key]
        r = rs[0]
        assert len(rs) == int(aggregate["n_observation_rows"])
        for field in ["standardized_isomeric_smiles", "physicochemical_descriptors_json", "primary_calibrated_probability", "direct_panel_probability", "accepted_prediction"]:
            assert len({x[field] for x in rs}) == 1, (key, field)
        assert true(aggregate["full"]) == true(r["accepted_prediction"])
        er = 2 ** float(aggregate["log2_er"])
        feature = make_feature("ASAP56", key, r, identity["source_canonical_smiles"], int(float(aggregate["log2_er"]) >= 1), er, aggregate["cluster"], len(rs), identity["source_exact_inchikey_computed"])
        assert feature["er_ge2_weak_label"] == int(identity["aggregated_weak_er_ge2_label"])
        if feature["development_source_connectivity_overlap"]:
            excluded.append({"source": "ASAP56", "uid": key, "reason": "source_or_standardized_reference_connectivity_overlap", "standardized_inchikey": key})
        else: features.append(feature)
    adc = []
    for r in read(INPUTS["adc14"]):
        feature = make_feature("ADC14", r["entity_id"], r, r["source_smiles"], None, None, "", 0, r["source_inchikey"])
        adc.append(feature)
    assert len(adc) == 14
    for source in ["Biogen2624", "ASAP56", "ADC14"]:
        rs = [r for r in features + adc if r["source"] == source]
        assert len(rs) == len({r["uid"] for r in rs})
        assert len(rs) == len({r["connectivity_block"] for r in rs}), source
    summaries = [counts(features, source, region, adc) for source in ["Biogen2624", "ASAP56"] for region in REGIONS]
    cross_source = sorted({r["connectivity_block"] for r in features if r["source"] == "Biogen2624"} & {r["connectivity_block"] for r in features if r["source"] == "ASAP56"})
    adc_general_overlap = {source: sorted({r["connectivity_block"] for r in features if r["source"] == source} & {r["connectivity_block"] for r in adc}) for source in ["Biogen2624", "ASAP56"]}
    write("features.csv", features + adc)
    write("exclusion_ledger.csv", excluded)
    write("count_feasibility.csv", [{k: v for k, v in r.items() if k != "count_gate_shortfalls"} for r in summaries])
    schema = {
        "feature_table": "features.csv", "primary_analysis_filter": "source in Biogen2624 or ASAP56; apply one prespecified region boolean; ADC14 has no measured labels",
        "unique_unit": "connectivity_block within source", "join_key": "source + uid", "model_graph": "standardized_isomeric_smiles",
        "chemical_stress_features": {"molecular_weight": "existing frozen descriptor JSON", "tpsa": "existing frozen descriptor JSON", "amide_count": protocol["amide_definition"]},
        "labels": {"er_ge2_weak_label": "0/1 for general transport sources, empty for ADC14", "measured_er": "Biogen original ER; ASAP existing median-log2ER derived ER; empty for ADC14"},
        "frozen_predictions": {"full_accepted": "bool existing full v0.4 mask", "primary_score": "fixed v0.1 channel", "auxiliary_score": "fixed direct v0.3 channel"},
        "scaffold_cluster": "murcko_scaffold; empty original scaffold falls back to entity connectivity",
        "status_boundaries": ["already observed labels; retrospective", "ASAP56 aggregation is known posthoc secondary sensitivity", "general molecule negative does not equal ADC negative", "stress regions are operational thresholds, not validated ADC applicability domains"],
        "source_field_paths": {"Biogen": str(INPUTS["biogen"].relative_to(ROOT)), "ASAP62": str(INPUTS["asap62"].relative_to(ROOT)), "ASAP56": str(INPUTS["asap56"].relative_to(ROOT)), "ADC14": str(INPUTS["adc14"].relative_to(ROOT))},
    }
    dump("input_schema.json", schema)
    report = {"status": "PASS_count_only_feasibility", "counts": summaries, "adc14": [{k: r[k] for k in ["uid", "entity_name", "molecular_weight", "tpsa", "amide_count", "MW_ge600", "TPSA_ge140", "amide_ge4_MW_ge500", "full_accepted"]} for r in adc], "excluded": excluded,
              "reference_unique_exact_key_n": len(reference_keys), "reference_unique_connectivity_n": len(reference_conn),
              "cross_source_connectivity_overlap": cross_source, "ADC_general_source_connectivity_overlap": adc_general_overlap,
              "all_source_graphs_parse_and_frozen_model_key_reproduces": True, "no_performance_computed": True,
              "all_inputs_unchanged": hashes == {str(p.relative_to(ROOT)): sha(p) for p in INPUTS.values()}}
    assert report["all_inputs_unchanged"]
    dump("count_feasibility.json", report)
    completion = {"status": "PASS", "completed_utc": datetime.now(timezone.utc).isoformat(), "input_sha256": hashes,
                  "runtime": {"executable": sys.executable, "python": platform.python_version(), "rdkit": rdkit.__version__},
                  "feature_rows": len(features + adc), "counts_by_source": dict(Counter(r["source"] for r in features + adc)),
                  "no_model_prediction_refit_or_performance_computation": True, "all_inputs_unchanged": True}
    dump("completion.json", completion)
    lines = ["# 预定义化学压力区：计数可行性", "", "本轮不新增检索、训练、预测或模型性能计算。三个区域按结构特征事先冻结；它们是操作性化学压力区，不能称为已验证的 ADC 领域。普通分子的 ER<2 记录不能称为 ADC 负例。", "",
             "Biogen 采用 identity_audited_v2 的 identity_clean=yes/model_usable=yes 既有唯一 connectivity 人群。ASAP 的 56 个唯一化合物来自已知 posthoc median-log2ER 聚合，属于既有结果后的敏感性复核。ER≥2 是来源特定弱端点，缺少共同亲本或抑制剂确认。", "",
             "|来源|结构区域|唯一分子|ER≥2|ER<2|scaffold|full接受|接受阳性|接受阴性|ADC14落入|计数门槛|", "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in summaries:
        lines.append(f"|{r['source']}|{r['region']}|{r['unique_compound_n']}|{r['positive_n']}|{r['negative_n']}|{r['scaffold_n']}|{r['full_accepted_n']}|{r['full_accepted_positive_n']}|{r['full_accepted_negative_n']}|{r['adc14_region_n']}/14|{'满足' if r['minimum_practical_count_gate_pass'] else '未满足'}|")
    lines += ["", "amide_count 使用已冻结标准化结构中的 [CX3](=[OX1])[NX3] 唯一匹配数量；这是酰胺基序计数，未验证为肽类分类器。MW/TPSA 均取既有模型 descriptor JSON，未重跑标准化。", "",
              "计数门槛仅用于有限回顾性同名额比较的实用性判断：总体 n≥30、正负各≥5、scaffold≥10，full接受≥10、接受正负各≥3。它不是统计功效分析，也不是性能或 ADC 泛化证明。所有三个区都应报告，不能依据性能选择某一区。", "",
              f"Biogen原表{len(biogen)}行，保留{sum(r['source']=='Biogen2624' for r in features)}个；ASAP56保留{sum(r['source']=='ASAP56' for r in features)}个。开发参考 exact/connectivity 重叠已剔除，细目见 exclusion_ledger.csv。Biogen/ASAP跨来源connectivity重叠{len(cross_source)}个，未合并为新增独立观测。", "",
              "ADC14仅用于压力区覆盖映射，没有运输标签，因此不能从其落区比例计算 ADC 准确率。features.csv 给出结构、化学特征、冻结分数、既有 full 接受名单和弱 ER 标签；input_schema.json 指定可复用字段和边界。", "", "原始输入 SHA256 和 Python/RDKit 环境见 completion.json。"]
    (HERE / "count_feasibility.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"counts": summaries, "source_n": completion["counts_by_source"], "additional_excluded": [r for r in excluded if r['reason'] != 'historical_identity_clean_or_model_usable_exclusion']}, ensure_ascii=False))
