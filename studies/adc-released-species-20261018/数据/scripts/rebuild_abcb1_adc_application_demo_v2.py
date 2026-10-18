#!/usr/bin/env python3
"""Rebuild the ADC-related application set and score it with frozen ABCB1 v0.4.

This workflow keeps three distinct states:
1) a release-related entity with a model-eligible structure;
2) a model abstention after valid structural scoring; and
3) an entity that is not scored because its chemical identity is unresolved.

The script never overwrites the historical 10-row registry or its predictions.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import sys
from collections import Counter
from datetime import date
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors


DATA = Path(__file__).resolve().parents[1]
PROJECT = DATA.parent
PROCESSED = DATA / "processed"
OUTPUT_DIR = PROJECT / "outputs" / "20260810_adc_demo_reconstruction_v2"
MODEL_DIR = PROJECT / "模型" / "ABCB1_selective_screener_v0_4_dev"
PRIMARY_MODEL = PROJECT / "模型" / "ABCB1_screener_v0_1" / "ABCB1_screener_v0_1.joblib"
DIRECT_MODEL = PROJECT / "模型" / "ABCB1_assay_aware_candidate_v0_3_dev" / "ABCB1_assay_aware_candidate_v0_3_dev.joblib"
SELECTIVE_MODEL = MODEL_DIR / "ABCB1_selective_uncertainty_v0_4_dev.joblib"

SOURCE_REGISTRY = PROCESSED / "adc_release_species_registry_v4.csv"
EVIDENCE_REVIEW = PROCESSED / "adc_demo_release_evidence_review_v2.csv"
TRANSPORT_CONTEXT = PROCESSED / "adc_release_species_transporter_evidence_long_v6.csv"

REGISTRY_OUT = PROCESSED / "adc_application_demo_registry_reconstructed_v2.csv"
PREDICTIONS_OUT = PROCESSED / "abcb1_adc_application_demo_predictions_v2.csv"
SUMMARY_OUT = PROCESSED / "abcb1_adc_application_demo_summary_v2.json"
SOURCES_OUT = PROCESSED / "adc_application_demo_sources_v2.csv"
CORRECTIONS_OUT = PROCESSED / "adc_application_demo_source_corrections_v2.csv"
TRANSPORT_OUT = PROCESSED / "adc_application_demo_abcb1_context_v2.csv"
VALIDATION_OUT = PROCESSED / "adc_application_demo_validation_v2.json"
MANIFEST_OUT = OUTPUT_DIR / "adc_application_demo_manifest_v2.json"
REPORT_OUT = OUTPUT_DIR / "ADC应用示范集重建与ABCB1重评分报告_20260810.md"
WORKBOOK_OUT = OUTPUT_DIR / "ADC应用示范对象分层与ABCB1重评分审计_20260810.xlsx"
MANUSCRIPT = PROJECT / "论文主稿" / "ABCB1化学空间迁移与选择性拒判_中文主稿_v0.1_20260809.md"

sys.path.insert(0, str(MODEL_DIR))
sys.path.insert(0, str(DATA / "scripts"))

from abcb1_screener_core_v0_1 import standardize_smiles  # noqa: E402
from predict_abcb1_selective import flatten_batch_result, load_packages, predict_one  # noqa: E402


RELEASE_GRADE_MAP = {
    "A1": "R3",
    "A2": "R2",
    "B1": "R2",
    "B2": "R2",
    "C": "R1",
}

IDENTITY_GRADES = {
    "ADC_REL_001": "ID2",
    "ADC_REL_002": "ID1",
    "ADC_REL_003": "ID2",
    "ADC_REL_004": "ID3",
    "ADC_REL_005": "ID2",
    "ADC_REL_006": "ID2",
    "ADC_REL_007": "ID3",
    "ADC_REL_008": "ID2",
    "ADC_REL_009": "ID3",
    "ADC_REL_010": "ID3",
    "ADC_REL_011": "ID3",
    "ADC_REL_012": "ID0",
    "ADC_REL_013": "ID3",
}

RELATIONSHIP_ROLE = {
    "ADC_REL_001": "analytically_confirmed_release",
    "ADC_REL_002": "analytically_confirmed_release",
    "ADC_REL_003": "analytically_confirmed_release",
    "ADC_REL_004": "analytically_confirmed_release",
    "ADC_REL_005": "source_assigned_release",
    "ADC_REL_006": "analytically_confirmed_release",
    "ADC_REL_007": "analytically_confirmed_release",
    "ADC_REL_008": "analytically_confirmed_release",
    "ADC_REL_009": "analytically_confirmed_release",
    "ADC_REL_010": "analytically_confirmed_release",
    "ADC_REL_011": "analytically_confirmed_release_systemic_unconjugated",
    "ADC_REL_012": "source_assigned_release",
    "ADC_REL_013": "source_assigned_release",
}

ABCB1_CONTEXT = {
    "ADC_REL_001": (
        "direct_qualitative_low_or_uncertain",
        "Human MDR1-LLC-PK1 bidirectional experiment was interpreted by the source as low P-gp potential; figure-level values remain a separate curation task.",
        "contextual_not_accuracy_grade",
    ),
    "ADC_REL_002": (
        "comparative_cellular_retention_only",
        "Relative retention in MDR1-expressing cells is not a direct transporter truth label.",
        "no",
    ),
    "ADC_REL_003": ("none_endpoint_aligned", "No endpoint-aligned exact-entity ABCB1 truth was used.", "no"),
    "ADC_REL_004": (
        "regulatory_substrate_call",
        "A regulatory review calls S-methyl-DM4 a P-gp substrate but does not reproduce numeric transport data.",
        "call_only",
    ),
    "ADC_REL_005": (
        "functional_resistance_positive",
        "MDR1-associated cytotoxicity shifts support susceptibility but are not bidirectional flux measurements.",
        "functional_only",
    ),
    "ADC_REL_006": (
        "direct_numeric_positive",
        "Human P-gp-MDCK ER was 34 and fell to 1.1 with verapamil in the regulatory summary.",
        "yes_direct_numeric",
    ),
    "ADC_REL_007": (
        "direct_numeric_positive",
        "Independent Caco-2 and MDR1-MDCK studies report high ER with inhibitor rescue.",
        "yes_direct_numeric",
    ),
    "ADC_REL_008": ("none_endpoint_aligned", "No exact-entity endpoint-aligned ABCB1 truth was used.", "no"),
    "ADC_REL_009": (
        "Caco2_inhibitor_attributed_positive",
        "Caco-2 experiments support P-gp involvement, but the system is not single-gene isolated.",
        "contextual_continuous",
    ),
    "ADC_REL_010": (
        "regulatory_plus_Caco2_positive",
        "Regulatory substrate calls and Caco-2 inhibitor-attributed ER are directionally positive but heterogeneous.",
        "contextual_continuous",
    ),
    "ADC_REL_011": (
        "context_dependent",
        "Caco-2 results are low or source-called non-substrate in one study; no clean universal binary truth is assigned.",
        "no_clean_binary",
    ),
    "ADC_REL_012": (
        "Caco2_low_structure_unresolved",
        "The source reports ER 1.6, but the exact structure crosswalk is unresolved and the system is not gene isolated.",
        "no",
    ),
    "ADC_REL_013": (
        "regulatory_plus_functional_positive",
        "FDA P-gp substrate call and MDR1-transfectant functional resistance are directionally concordant; neither is a reproduced direct-flux dataset here.",
        "supportive_not_accuracy_grade",
    ),
}

SOURCE_CORRECTIONS = [
    {
        "source_id": "SRC_B2_004",
        "field": "author_string_and_journal",
        "historical_value": "Wu W, Xu Z, Zhan C, Yin X, Yu S.; Environ Monit Assess",
        "corrected_value": "Rock BM, Tometsko ME, Patel SK, Hamblett KJ, Fanslow WC, Rock DA; Drug Metab Dispos",
        "unchanged_identifier": "PMID 26108745; DOI 10.1124/dmd.115.064253",
        "basis_url": "https://doi.org/10.1124/dmd.115.064253",
        "correction_effect": "metadata_only_release_evidence_unchanged",
    },
    {
        "source_id": "SRC_B2_013",
        "field": "pmid_author_string_and_journal",
        "historical_value": "PMID 22829590; Sakurai A, Matsuda M, Kiyokawa E.; J Biol Chem",
        "corrected_value": "PMID 22752008; Davis JA, Rock DA, Wienkers LC, Pearson JT; Drug Metab Dispos",
        "unchanged_identifier": "DOI 10.1124/dmd.112.046169",
        "basis_url": "https://pubmed.ncbi.nlm.nih.gov/22752008/",
        "correction_effect": "metadata_and_PubMed_crosswalk_repaired_release_grade_remains_R2",
    },
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, obj: dict | list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def chiral_counts(smiles: str) -> tuple[str, str]:
    if not smiles:
        return "", ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "", ""
    centers = Chem.FindMolChiralCenters(
        mol,
        includeUnassigned=True,
        includeCIP=True,
        useLegacyImplementation=False,
    )
    assigned = sum(label != "?" for _, label in centers)
    unassigned = sum(label == "?" for _, label in centers)
    return str(assigned), str(unassigned)


def descriptor_fields(standardized_smiles: str) -> dict[str, str]:
    mol = Chem.MolFromSmiles(standardized_smiles) if standardized_smiles else None
    if mol is None:
        return {
            "std_molecular_weight": "",
            "std_xlogp": "",
            "std_tpsa": "",
            "std_hbond_donor_count": "",
            "std_hbond_acceptor_count": "",
            "std_rotatable_bond_count": "",
            "std_formal_charge": "",
            "std_heavy_atom_count": "",
            "std_molecular_formula": "",
        }
    return {
        "std_molecular_weight": f"{Descriptors.MolWt(mol):.4f}",
        "std_xlogp": f"{Crippen.MolLogP(mol):.4f}",
        "std_tpsa": f"{rdMolDescriptors.CalcTPSA(mol):.4f}",
        "std_hbond_donor_count": str(Lipinski.NumHDonors(mol)),
        "std_hbond_acceptor_count": str(Lipinski.NumHAcceptors(mol)),
        "std_rotatable_bond_count": str(Lipinski.NumRotatableBonds(mol)),
        "std_formal_charge": str(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
        "std_heavy_atom_count": str(mol.GetNumHeavyAtoms()),
        "std_molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
    }


def model_input_type(audited_object_class: str) -> str:
    if audited_object_class == "bare_free_payload":
        return "free_payload"
    if "linker" in audited_object_class or "catabolite" in audited_object_class:
        return "released_linker_payload"
    return "released_small_molecule"


def summarize(label: str, rows: list[dict]) -> dict:
    evaluable = [row for row in rows if row["model_score_eligible"] == "yes"]
    accepted = [row for row in evaluable if row["accepted_prediction"] == "yes"]
    abstained = [row for row in evaluable if row["final_decision_state"] == "model_abstain"]
    not_scored = [row for row in rows if row["model_score_eligible"] != "yes"]
    coverage = len(accepted) / len(evaluable) if evaluable else None
    return {
        "stratum": label,
        "registry_entities": len(rows),
        "model_evaluable": len(evaluable),
        "accepted": len(accepted),
        "model_abstain": len(abstained),
        "not_scored_identity": len(not_scored),
        "high_risk": sum(row["risk_call"] == "high_ABCB1_efflux_risk" for row in accepted),
        "low_risk": sum(row["risk_call"] == "low_ABCB1_efflux_risk" for row in accepted),
        "coverage": round(coverage, 8) if coverage is not None else None,
        "accepted_entity_ids": [row["entity_id"] for row in accepted],
        "abstained_entity_ids": [row["entity_id"] for row in abstained],
        "not_scored_entity_ids": [row["entity_id"] for row in not_scored],
    }


def explode_sources(evidence_rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    seen: set[tuple[str, str, str]] = set()
    for row in evidence_rows:
        for role, field in (("primary_release_or_identity", "primary_source_url"), ("supporting", "supporting_source_url")):
            urls = [item.strip() for item in row.get(field, "").split("|") if item.strip()]
            for rank, url in enumerate(urls, start=1):
                key = (row["reference_id"], role, url)
                if key in seen:
                    continue
                seen.add(key)
                output.append(
                    {
                        "source_link_id": f"ADC_DEMO_SRC_{len(output) + 1:03d}",
                        "entity_id": row["reference_id"],
                        "entity_name": row["reference_name"],
                        "evidence_role": role,
                        "within_role_rank": rank,
                        "source_url": url,
                        "supports": "release_relationship_and_or_chemical_identity",
                        "evidence_semantics": row.get("evidence_semantics", ""),
                        "allowed_claim": row.get("allowed_claim", ""),
                        "prohibited_claim": row.get("prohibited_claim", ""),
                        "review_date": row.get("review_date", ""),
                    }
                )
    return output


def build_report(summary: dict, rows: list[dict]) -> str:
    main = summary["strata"]["R2_plus_R3_main"]
    strict = summary["strata"]["R3_strict"]
    extended = summary["strata"]["R2_extended"]
    table_lines = [
        "| 实体 | 释放证据 | 身份等级 | 分析对象 | 最终状态 | 主模型概率 | 主要限制 |",
        "|---|---|---|---|---|---:|---|",
    ]
    for row in rows:
        probability = row.get("primary_calibrated_probability", "")
        if isinstance(probability, float):
            probability = f"{probability:.3f}"
        elif probability not in ("", None):
            try:
                probability = f"{float(probability):.3f}"
            except (TypeError, ValueError):
                probability = str(probability)
        table_lines.append(
            "| {name} | {release} | {identity} | {obj} | {state} | {prob} | {limit} |".format(
                name=row["entity_name"],
                release=row["release_evidence_grade"],
                identity=row["identity_grade"],
                obj=row["audited_object_class"],
                state=row["final_decision_state"],
                prob=probability or "NA",
                limit=(row.get("ineligibility_reason") or row.get("abstention_reasons_json") or "无"),
            )
        )
    return f"""# ADC释放相关实体应用集重建与ABCB1 v0.4重评分报告

## 结论

数据锁定于{summary['data_lock_date']}。历史10条快照已重建为13个唯一ADC释放相关化学实体。13条均达到R2或R3释放关系证据，其中R3分析确认层{strict['registry_entities']}条、R2来源指认层{extended['registry_entities']}条。11条具有ID2或ID3机器结构并进入冻结v0.4，2条因身份或结构未闭合而记为`not_scored_identity`，不再并入模型拒判。

主应用集接受{main['accepted']}/{main['model_evaluable']}条，选择性覆盖率为{main['coverage'] * 100:.1f}%；接受结果全部为高ABCB1外排风险。另有{main['model_abstain']}条有效结构被模型拒判，主要原因仍是同时缺乏两套参考面板支持或超出理化包络。R3严格层接受{strict['accepted']}/{strict['model_evaluable']}条，覆盖率为{strict['coverage'] * 100:.1f}%。这些比例描述当前登记实体上的接口行为，不是ADC专用准确率，也不能估计ADC载荷总体风险流行率。

## 对历史结果的修正

旧表的“3条高风险、7条拒判”混合了两种状态。按相同口径重新表述应为：历史10个登记实体中9个可评分，3/9被接受，6/9为模型拒判，另1/10因结构未解析未评分。v2同时补入Exatecan、BAY-716和SG3199，并执行证据与身份门控。

## 分层定义

- R3：实体由LC-MS/MS、结构性MS、放射性代谢物分析或其他实体分辨方法直接确认；Exatecan在OBI-992给药后作为系统性未偶联载荷测得，但不据此声称其全部来自胞内裂解。
- R2：原始论文或监管文件明确指认释放物或活性分解物，但当前未获得命名ADC后的实体分辨释放测量。
- ID3：精确结构和身份交叉核验通过；ID2：连接结构明确但存在立体、混合物、互变异构或命名表示局限；ID1/ID0不输入模型。
- `model_abstain`只用于有效结构已进入模型但证据支持不足的对象；`not_scored_identity`用于结构或身份未闭合对象。

## 逐实体结果

{chr(10).join(table_lines)}

## 证据与模型边界

本应用集按ADC释放关系和化学身份构建，不读取ABCB1方向或模型输出决定纳入。已有ABCB1结果仅作为预测后案例注释，因其混合了直接双向转运、监管定性调用、Caco-2抑制剂归因和功能耐药端点，未合并计算准确率。模型对象为标准化的游离或释放后小分子，不是完整ADC，也不预测胞内游离暴露、疗效、耐药或临床DDI。

## 数据修复

- `SRC_B2_004`的标题、PMID和DOI正确，但历史作者与期刊元数据错配；已在v2修正表登记。
- `SRC_B2_013`的DOI正确，但历史PMID、作者与期刊错配；正确PMID为22752008。
- v4中Exatecan仍写为仅作者指认。2025和2026原始研究已用LC-MS/MS测得OBI-992后的游离/未偶联Exatecan，v2据此升级为R3，同时保留“非胞内来源特异”限制。
"""


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    source_rows = read_csv(SOURCE_REGISTRY)
    evidence_rows = read_csv(EVIDENCE_REVIEW)
    source_by_id = {row["release_species_id"]: row for row in source_rows}
    evidence_by_id = {row["reference_id"]: row for row in evidence_rows}

    expected_ids = [f"ADC_REL_{index:03d}" for index in range(1, 14)]
    if sorted(source_by_id) != expected_ids:
        raise AssertionError(f"Unexpected v4 registry IDs: {sorted(source_by_id)}")
    if sorted(evidence_by_id) != expected_ids:
        raise AssertionError(f"Evidence review must contain all 13 IDs: {sorted(evidence_by_id)}")

    selective, primary, direct = load_packages(SELECTIVE_MODEL, PRIMARY_MODEL, DIRECT_MODEL)
    model_hash_before = sha256_file(SELECTIVE_MODEL)
    legacy_hashes_before = {
        "legacy_registry_v1": sha256_file(PROCESSED / "adc_release_species_registry_standardized_v1.csv"),
        "legacy_predictions_v1": sha256_file(PROCESSED / "abcb1_selective_v0_4_adc_release_registry_predictions_v1.csv"),
        "legacy_summary_v1": sha256_file(PROCESSED / "abcb1_selective_v0_4_adc_release_registry_summary_v1.json"),
    }

    registry_rows: list[dict] = []
    prediction_rows: list[dict] = []
    for entity_id in expected_ids:
        source = source_by_id[entity_id]
        evidence = evidence_by_id[entity_id]
        source_smiles = source.get("isomeric_smiles", "").strip()
        standardized = standardize_smiles(source_smiles)
        identity_grade = IDENTITY_GRADES[entity_id]
        release_grade = RELEASE_GRADE_MAP[evidence["release_evidence_grade"]]
        model_eligible = identity_grade in {"ID2", "ID3"} and bool(source_smiles)
        source_assigned, source_unassigned = chiral_counts(source_smiles)
        std_assigned, std_unassigned = chiral_counts(standardized.get("standardized_isomeric_smiles", ""))
        source_block = source.get("inchikey", "").split("-")[0] if source.get("inchikey") else ""
        standardized_key = standardized.get("standardized_inchikey", "")
        standardized_block = standardized_key.split("-")[0] if standardized_key else ""
        abcb1_level, abcb1_summary, validation_eligible = ABCB1_CONTEXT[entity_id]
        ineligibility_reason = ""
        if entity_id == "ADC_REL_002":
            ineligibility_reason = "no_independently_verified_machine_readable_structure"
        elif entity_id == "ADC_REL_012":
            ineligibility_reason = "exact_name_to_structure_crosswalk_unresolved"

        registry = {
            "panel_record_id": f"ADC_APP_V2_{int(entity_id[-3:]):03d}",
            "entity_id": entity_id,
            "entity_name": source["release_species_name"],
            "entity_name_as_source": source["release_species_name"],
            "aliases": source.get("aliases", ""),
            "representative_adc_or_format": evidence.get("representative_adc_or_format", source.get("representative_adc_or_format", "")),
            "payload_core": source.get("payload_core", ""),
            "residual_group": source.get("residual_group", ""),
            "audited_object_class": evidence.get("audited_object_class", source.get("species_class", "")),
            "legacy_species_class": source.get("species_class", ""),
            "linker_class": source.get("linker_class", ""),
            "adc_relationship_role": RELATIONSHIP_ROLE[entity_id],
            "release_stage": evidence.get("terminality", ""),
            "release_evidence_grade": release_grade,
            "release_evidence_grade_source_scale": evidence.get("release_evidence_grade", ""),
            "release_identity_verdict": evidence.get("release_identity_verdict", ""),
            "release_evidence_method": source.get("release_evidence_method", ""),
            "evidence_semantics": evidence.get("evidence_semantics", ""),
            "primary_source_url": evidence.get("primary_source_url", ""),
            "supporting_source_url": evidence.get("supporting_source_url", ""),
            "identity_grade": identity_grade,
            "structure_evidence_grade_source_scale": evidence.get("structure_evidence_grade", ""),
            "structure_identity_verdict": evidence.get("structure_identity_verdict", ""),
            "structure_status": source.get("structure_status", ""),
            "identity_crosswalk_status": evidence.get("structure_identity_verdict", ""),
            "source_pubchem_cid": source.get("pubchem_cid", ""),
            "source_inchikey": source.get("inchikey", ""),
            "source_connectivity_block": source_block,
            "source_isomeric_smiles": source_smiles,
            "source_connectivity_smiles": source.get("connectivity_smiles", ""),
            "source_structure_url": source.get("structure_source_url", ""),
            "source_assigned_stereocenters": source_assigned,
            "source_unassigned_stereocenters": source_unassigned,
            "registry_standardization_status": standardized.get("standardization_status", ""),
            "standardized_isomeric_smiles": standardized.get("standardized_isomeric_smiles", ""),
            "standardized_connectivity_smiles": standardized.get("standardized_connectivity_smiles", ""),
            "standardized_inchikey": standardized_key,
            "standardized_connectivity_block": standardized_block,
            "standardized_assigned_stereocenters": std_assigned,
            "standardized_unassigned_stereocenters": std_unassigned,
            "source_to_standardized_connectivity_block_changed": yes_no(bool(source_block and standardized_block and source_block != standardized_block)),
            "tautomer_status": standardized.get("tautomer_status", ""),
            "fragment_removed_flag": standardized.get("fragment_removed_flag", ""),
            **descriptor_fields(standardized.get("standardized_isomeric_smiles", "")),
            "primary_demo_eligible": yes_no(release_grade in {"R2", "R3"}),
            "strict_r3_eligible": yes_no(release_grade == "R3"),
            "model_score_eligible": yes_no(model_eligible),
            "model_input_object_type": model_input_type(evidence.get("audited_object_class", "")),
            "ineligibility_reason": ineligibility_reason,
            "legacy_modeling_status": source.get("modeling_status", ""),
            "allowed_claim": evidence.get("allowed_claim", ""),
            "prohibited_claim": evidence.get("prohibited_claim", ""),
            "abcb1_evidence_level": abcb1_level,
            "abcb1_evidence_summary": abcb1_summary,
            "validation_eligible": validation_eligible,
            "data_lock_date": "2026-08-10",
            "registry_version": "ADC_application_demo_v2",
        }
        registry_rows.append(registry)

        if model_eligible:
            result = flatten_batch_result(
                predict_one(
                    source_smiles,
                    selective,
                    primary,
                    direct,
                    chemical_object_type=registry["model_input_object_type"],
                    requested_context="structure_only_screening",
                )
            )
            final_state = result["risk_call"] if result["accepted_prediction"] == "yes" else "model_abstain"
            prediction = {
                **registry,
                **result,
                "final_decision_state": final_state,
                "application_role": "ADC_release_related_entity_research_prioritization_not_accuracy_validation",
            }
        else:
            prediction = {
                **registry,
                "input_smiles": source_smiles,
                "chemical_object_type": registry["model_input_object_type"],
                "requested_context": "structure_only_screening",
                "context_compatibility": "not_evaluated_identity_gate",
                "risk_call": "not_scored",
                "decision_zh": "身份未解析，不予评分",
                "accepted_prediction": "no",
                "primary_calibrated_probability": "",
                "primary_low_threshold": "",
                "primary_high_threshold": "",
                "bootstrap_probability_p10": "",
                "bootstrap_probability_median": "",
                "bootstrap_probability_p90": "",
                "bootstrap_interval_width": "",
                "direct_panel_probability": "",
                "primary_max_tanimoto": "",
                "direct_max_tanimoto": "",
                "physchem_max_abs_robust_z": "",
                "physchem_threshold_q99": "",
                "applicability_status": "not_evaluated_identity_gate",
                "uncertainty_status": "not_evaluated_identity_gate",
                "exact_reference_identity": "not_evaluated",
                "stereo_blind_collision_warning": "not_evaluated",
                "model_version": selective["model_version"],
                "validation_status": selective["status"],
                "claim_boundary": selective["claim_boundary"],
                "abstention_reasons_json": json.dumps(["identity_or_structure_unresolved"], ensure_ascii=False),
                "abstention_reasons_zh_json": json.dumps(["化学身份或机器可读结构未闭合"], ensure_ascii=False),
                "nearest_primary_references_json": "[]",
                "nearest_direct_references_json": "[]",
                "physicochemical_descriptors_json": "{}",
                "final_decision_state": "not_scored_identity",
                "application_role": "ADC_release_related_entity_registry_context_only_not_model_scored",
            }
        prediction_rows.append(prediction)

    write_csv(REGISTRY_OUT, registry_rows)
    write_csv(PREDICTIONS_OUT, prediction_rows)
    source_links = explode_sources(evidence_rows)
    write_csv(SOURCES_OUT, source_links)
    write_csv(CORRECTIONS_OUT, SOURCE_CORRECTIONS)

    transport_rows = [
        row
        for row in read_csv(TRANSPORT_CONTEXT)
        if row.get("release_species_id", "") in set(expected_ids)
    ]
    write_csv(TRANSPORT_OUT, transport_rows)

    all_main = [row for row in prediction_rows if row["release_evidence_grade"] in {"R2", "R3"}]
    strict_r3 = [row for row in prediction_rows if row["release_evidence_grade"] == "R3"]
    extended_r2 = [row for row in prediction_rows if row["release_evidence_grade"] == "R2"]
    summary = {
        "dataset_version": "ADC_application_demo_v2",
        "data_lock_date": "2026-08-10",
        "model_version": selective["model_version"],
        "model_package_sha256": model_hash_before,
        "application_status": "stratified_ADC_release_related_entity_demonstration_not_accuracy_validation",
        "strata": {
            "R2_plus_R3_main": summarize("R2_plus_R3_main", all_main),
            "R3_strict": summarize("R3_strict", strict_r3),
            "R2_extended": summarize("R2_extended", extended_r2),
        },
        "release_grade_counts": dict(Counter(row["release_evidence_grade"] for row in prediction_rows)),
        "identity_grade_counts": dict(Counter(row["identity_grade"] for row in prediction_rows)),
        "decision_state_counts": dict(Counter(row["final_decision_state"] for row in prediction_rows)),
        "model_abstention_reason_counts": dict(
            Counter(
                reason
                for row in prediction_rows
                if row["final_decision_state"] == "model_abstain"
                for reason in json.loads(row["abstention_reasons_json"])
            )
        ),
        "historical_snapshot_reexpression": {
            "registry_entities": 10,
            "model_evaluable": 9,
            "accepted": 3,
            "model_abstain": 6,
            "not_scored_identity": 1,
            "coverage": round(3 / 9, 8),
            "note": "The historical phrase 3 high-risk and 7 abstentions mixed one unscored identity with six model abstentions.",
        },
        "claim_boundary": selective["claim_boundary"],
        "interpretation": (
            "This case series demonstrates identity gating, selective scoring and machine-readable abstention for ADC-related released entities. "
            "It is not an ADC-specific accuracy estimate, a prevalence sample, an intact-ADC transport model or a universal substrate truth set."
        ),
        "processed_outputs": {
            "registry": str(REGISTRY_OUT.relative_to(PROJECT)),
            "predictions": str(PREDICTIONS_OUT.relative_to(PROJECT)),
            "sources": str(SOURCES_OUT.relative_to(PROJECT)),
            "source_corrections": str(CORRECTIONS_OUT.relative_to(PROJECT)),
            "abcb1_context": str(TRANSPORT_OUT.relative_to(PROJECT)),
        },
    }
    write_json(SUMMARY_OUT, summary)

    expected_accepted = {"ADC_REL_007", "ADC_REL_009", "ADC_REL_010", "ADC_REL_011", "ADC_REL_013"}
    expected_abstained = {"ADC_REL_001", "ADC_REL_003", "ADC_REL_004", "ADC_REL_005", "ADC_REL_006", "ADC_REL_008"}
    expected_not_scored = {"ADC_REL_002", "ADC_REL_012"}
    actual_accepted = {row["entity_id"] for row in prediction_rows if row["accepted_prediction"] == "yes"}
    actual_abstained = {row["entity_id"] for row in prediction_rows if row["final_decision_state"] == "model_abstain"}
    actual_not_scored = {row["entity_id"] for row in prediction_rows if row["final_decision_state"] == "not_scored_identity"}

    checks: list[dict] = []

    def json_safe(value):
        if isinstance(value, set):
            return sorted(value)
        if isinstance(value, dict):
            return {key: json_safe(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [json_safe(item) for item in value]
        return value

    def check(check_id: str, observed, expected) -> None:
        checks.append(
            {
                "check_id": check_id,
                "status": "pass" if observed == expected else "fail",
                "observed": json_safe(observed),
                "expected": json_safe(expected),
            }
        )

    check("registry_13_unique", len(registry_rows), 13)
    check("registry_unique_ids", len({row["entity_id"] for row in registry_rows}), 13)
    check("release_grade_counts", summary["release_grade_counts"], {"R3": 10, "R2": 3})
    check("identity_grade_counts", summary["identity_grade_counts"], {"ID2": 5, "ID1": 1, "ID3": 6, "ID0": 1})
    check("model_evaluable_11", sum(row["model_score_eligible"] == "yes" for row in registry_rows), 11)
    check("accepted_entities", sorted(actual_accepted), sorted(expected_accepted))
    check("model_abstained_entities", sorted(actual_abstained), sorted(expected_abstained))
    check("not_scored_entities", sorted(actual_not_scored), sorted(expected_not_scored))
    check("all_accepted_are_high_risk", {row["risk_call"] for row in prediction_rows if row["accepted_prediction"] == "yes"}, {"high_ABCB1_efflux_risk"})
    check("strict_r3_registry", summary["strata"]["R3_strict"]["registry_entities"], 10)
    check("strict_r3_evaluable", summary["strata"]["R3_strict"]["model_evaluable"], 9)
    check("strict_r3_accepted", summary["strata"]["R3_strict"]["accepted"], 4)
    check("main_coverage", summary["strata"]["R2_plus_R3_main"]["coverage"], round(5 / 11, 8))
    check("strict_coverage", summary["strata"]["R3_strict"]["coverage"], round(4 / 9, 8))
    check("model_package_unchanged", sha256_file(SELECTIVE_MODEL), model_hash_before)
    check(
        "legacy_files_unchanged",
        {
            "legacy_registry_v1": sha256_file(PROCESSED / "adc_release_species_registry_standardized_v1.csv"),
            "legacy_predictions_v1": sha256_file(PROCESSED / "abcb1_selective_v0_4_adc_release_registry_predictions_v1.csv"),
            "legacy_summary_v1": sha256_file(PROCESSED / "abcb1_selective_v0_4_adc_release_registry_summary_v1.json"),
        },
        legacy_hashes_before,
    )
    check("source_corrections_count", len(SOURCE_CORRECTIONS), 2)
    check("source_links_nonempty", len(source_links) > 0, True)
    check("abcb1_context_nonempty", len(transport_rows) > 0, True)
    manuscript_text = MANUSCRIPT.read_text(encoding="utf-8")
    check(
        "manuscript_methods_synced",
        "### 2.10 ADC释放相关实体应用集的重建与分层评分" in manuscript_text,
        True,
    )
    check(
        "manuscript_results_synced",
        all(
            phrase in manuscript_text
            for phrase in [
                "注册表包含13个唯一ADC释放相关化学实体",
                "v0.4接受5/11个模型可评价实体",
                "选择性覆盖率为45.5%",
                "覆盖率为44.4%",
            ]
        ),
        True,
    )
    check(
        "historical_heading_removed",
        "### 3.10 ADC释放小分子的应用示范暴露了当前参考空间的不足" in manuscript_text,
        False,
    )
    check("workbook_export_present", WORKBOOK_OUT.exists() and WORKBOOK_OUT.stat().st_size > 10000, True)
    formula_scan = OUTPUT_DIR / "exported_workbook_formula_errors.ndjson"
    check(
        "workbook_formula_scan_clean",
        formula_scan.exists() and "matched 0 entries" in formula_scan.read_text(encoding="utf-8"),
        True,
    )

    failures = [row for row in checks if row["status"] != "pass"]
    validation = {
        "validation_version": "ADC_application_demo_validation_v2",
        "run_date": str(date.today()),
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failures),
        "checks_failed": len(failures),
        "status": "pass" if not failures else "fail",
        "checks": checks,
    }
    write_json(VALIDATION_OUT, validation)
    if failures:
        raise AssertionError(f"Validation failed: {failures}")

    report = build_report(summary, prediction_rows)
    REPORT_OUT.write_text(report, encoding="utf-8")

    packaged = [
        REGISTRY_OUT,
        PREDICTIONS_OUT,
        SUMMARY_OUT,
        SOURCES_OUT,
        CORRECTIONS_OUT,
        TRANSPORT_OUT,
        VALIDATION_OUT,
        EVIDENCE_REVIEW,
    ]
    for path in packaged:
        shutil.copy2(path, OUTPUT_DIR / path.name)

    manifest_files = [
        SOURCE_REGISTRY,
        EVIDENCE_REVIEW,
        SELECTIVE_MODEL,
        PRIMARY_MODEL,
        DIRECT_MODEL,
        *packaged,
        REPORT_OUT,
        WORKBOOK_OUT,
        MANUSCRIPT,
        OUTPUT_DIR / "workbook_summary_inspection.ndjson",
        OUTPUT_DIR / "exported_workbook_summary_inspection.ndjson",
        OUTPUT_DIR / "workbook_formula_errors.ndjson",
        OUTPUT_DIR / "exported_workbook_formula_errors.ndjson",
    ]
    manifest = {
        "manifest_version": "ADC_application_demo_manifest_v2",
        "data_lock_date": "2026-08-10",
        "model_runtime": {
            "python": "3.12.13",
            "numpy": "2.5.1",
            "scikit_learn": "1.9.0",
            "scipy": "1.18.0",
            "rdkit": "2026.03.4",
            "joblib": "1.5.3",
        },
        "files": [
            {
                "path": str(path.relative_to(PROJECT)) if path.is_relative_to(PROJECT) else str(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in manifest_files
        ],
    }
    write_json(MANIFEST_OUT, manifest)
    print(json.dumps({"summary": summary, "validation": validation}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
