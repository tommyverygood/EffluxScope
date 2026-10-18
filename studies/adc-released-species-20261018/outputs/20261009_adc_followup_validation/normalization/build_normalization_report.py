#!/usr/bin/env python3
"""Read-only result checks and report rendering, separate from frozen computation."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
from collections import Counter
from rdkit import Chem

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = ROOT / "outputs/20261009_adc_pair_experiments/pairs/results_curated_v2"

def read_csv(path):
    return list(csv.DictReader(path.open(encoding="utf-8-sig", newline="")))

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def num(x):
    return float(x)

def main():
    manifest = json.loads((HERE / "run_manifest.json").read_text())
    scores = read_csv(HERE / "variant_scores.csv")
    stage = read_csv(HERE / "stage_structure_audit.csv")
    atoms = read_csv(HERE / "source_atom_stereo_audit.csv")
    pairs = read_csv(HERE / "parent_release_variant_contrasts.csv")
    prior = {r["unique_structure_id"]: r for r in read_csv(BASE / "unique_structure_predictions.csv")}
    by = {(r["unique_structure_id"],r["variant"]):r for r in scores}
    src = [r for r in scores if r["variant"] == "source"]
    modern = []
    for row in stage:
        m = Chem.MolFromSmiles(row["isomeric_smiles"])
        centers = Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
        modern.append({"unique_structure_id": row["unique_structure_id"], "entity_name": row["entity_name"], "stage": row["stage"],
                       "assigned_tetrahedral_nonlegacy": sum(s != "?" for _,s in centers),
                       "unassigned_tetrahedral_nonlegacy": sum(s == "?" for _,s in centers),
                       "chiral_center_api": "FindMolChiralCenters(includeUnassigned=True,useLegacyImplementation=False)"})
    fields = list(modern[0])
    with (HERE / "nonlegacy_stereo_count_crosscheck.csv").open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(modern)
    extra = {"frozen_standardized_full_keys_match_saved_unique_predictions": all(r["inchikey"] == prior[r["unique_structure_id"]]["standardized_inchikey"] for r in stage if r["stage"] == "frozen_old"),
             "source_formula_and_charge_unchanged_after_no_tautomer_cleanup": all(by[(r["unique_structure_id"], "no_tautomer")]["formula"] == r["formula"] and by[(r["unique_structure_id"], "no_tautomer")]["formal_charge"] == r["formal_charge"] for r in src),
             "source_no_tautomer_full_keys_identical": all(by[(r["unique_structure_id"], "no_tautomer")]["inchikey"] == r["inchikey"] for r in src),
             "all_frozen_input_hashes_still_match_manifest": all(digest(ROOT / p) == h for p,h in manifest["input_sha256"].items())}
    (HERE / "additional_verification.json").write_text(json.dumps({"status":"pass" if all(extra.values()) else "FAIL","checks":extra},ensure_ascii=False,indent=2)+"\n")
    assert all(extra.values())
    lost = {(r["unique_structure_id"]):[a["source_atom_index_zero_based"] for a in atoms if a["unique_structure_id"] == r["unique_structure_id"] and a["stage"] == "frozen_old" and a["assigned_stereo_lost"] == "True"] for r in src}
    rows = ["# ADC 化学表示敏感性验证结果", "", "已完成 13 个独立 source 结构 × 5 种表示的冻结计算。模型、core 和输入 SHA256 前后不变；历史 26 个模型分数全部在 1e-8 内复现；5 个同身份控制 × 5 个 variant 的 25 项相等检查全部通过。原始化学身份与冻结模型输入表示继续分开保存。", "", "## 结果与机制", "", "10/13 个完整 InChIKey 改变，其中 2/13 connectivity block 改变。Cleanup、FragmentParent、Uncharger 在这批结构上不改变完整 key；变化来自互变异构步骤。默认枚举删除参与互变异构的立体标记，合计 27 个原 source 已定义四面体中心变为未定义/不再为该中心；保留立体设置恢复 source 的已定义中心数，但不会补充原本未定义的中心。", "", "两模型的 Morgan 指纹设置 includeChirality=False，8 个描述符也不编码立体信息。13/13 个仅去立体的阴性对照，feature 和两模型分数均精确相等。因此，保存正确 source 手性对化学身份、采购/分析对应很必要，但本模型不能据此区分立体异构体。不能把恢复 SMILES 的手性标记写成模型获得了立体化学识别能力。", "", "| 材料 | source 已定义中心 | 冻结后已定义中心 | source 零基原子编号的立体丢失 | 完整 key 改变 | connectivity block 改变 | 枚举状态 |", "|---|---:|---:|---|---|---|---|"]
    for r in src:
        f = by[(r["unique_structure_id"],"frozen_old")]
        t = next(t for t in stage if t["unique_structure_id"] == r["unique_structure_id"] and t["stage"] == "frozen_old")
        rows.append(f"| {r['entity_name']} | {r['assigned_tetrahedral']} | {f['assigned_tetrahedral']} | {', '.join(lost[r['unique_structure_id']]) or '无'} | {t['full_key_changed_from_source']} | {t['connectivity_block_changed_from_source']} | {t['enumeration_status']} |")
    rows += ["", "stage_structure_audit.csv 使用 RDKit 默认 legacy tetrahedral-center 检测。nonlegacy_stereo_count_crosscheck.csv 另以 useLegacyImplementation=False 对齐上一批审计；Lys-MCC-DM1 环上的潜在未定义中心计数会因该 API 选项不同而不同，已定义中心及本轮 27 项立体丢失结论不受影响。", "", "只有 cys-mcMMAF 和 Lys-MCC-DM1 改变键型/指纹/描述符与数值。两者旧标准化将 succinimide 表示变成 hydroxylated aromatic tautomer 的计算表示，分子式不变，但 HBD、TPSA、LogP 发生改变。其他 11 个结构 source 与冻结分数相等（历史 CSV 保留 8 位小数，差值小于 1e-8）。", "", "| 材料 | source v0.1 | 历史冻结 v0.1 | source−历史 v0.1 | source v0.3 | 历史冻结 v0.3 | source−历史 v0.3 | 指纹位差 |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in src:
        if int(r["morgan_bit_differences_vs_frozen"]) == 0:
            continue
        rows.append("| "+r["entity_name"]+" | "+" | ".join(f"{num(r[k]):.8f}" for k in ["v0_1_risk_score","v0_1_historical_frozen_score_not_newly_admitted","v0_1_score_minus_historical_frozen_diagnostic","v0_3_risk_score","v0_3_historical_frozen_score_not_newly_admitted","v0_3_score_minus_historical_frozen_diagnostic"])+" | "+r["morgan_bit_differences_vs_frozen"]+" |")
    rows += ["", "## 不完整枚举与 HOLD", "", "默认 MaxTautomers=1000、MaxTransforms=1000。cys-mcMMAF、Lys-MCC-DM1 均为 MaxTransformsReached，cys-VC-MMAE 为 MaxTautomersReached；默认与保留立体两种枚举都触发同样的 cap。本轮预冻结协议将这些 canonical variant 标 HOLD，空置新准入 score。另表 frozen_score_reproduction.csv 保留历史输出的数值复现，historical_frozen_* 列仅供诊断对照，不能绕过 HOLD。", "", "cap hit 本身不证明数值对输入顺序不稳定，也不证明输出结构错误；它表示枚举未完成，不能未经验证把该结果当作唯一充分 canonical 表示。原模型并未修改。source/no_tautomer 的 estimator 分数仅为替代表示敏感性诊断，未校准、未验证，不可替代原部署输入契约或完整 v0.4 判定。", "", "## 配对分数变化", "", "| 配对 | 阶段 | source Δv0.1 | 历史冻结 Δv0.1 | source Δv0.3 | 历史冻结 Δv0.3 |", "|---|---|---:|---:|---:|---:|"]
    for r in pairs:
        if r["variant"] != "source" or r["pair_type"] != "chemical_changed_parent_to_catabolite" or r["pair_status"] != "both_computed_diagnostic":
            continue
        label = "终末/下游" if r["may_count_as_terminal_chemical_change_pair"] == "True" else "处理中间体"
        rows.append("| "+r["parent_name"]+"→"+r["release_name"]+" | "+label+" | "+" | ".join(f"{num(r[k]):.8f}" for k in ["v0_1_release_minus_parent_score","v0_1_historical_frozen_pair_delta_not_newly_admitted","v0_3_release_minus_parent_score","v0_3_historical_frozen_pair_delta_not_newly_admitted"])+" |")
    rows += ["", "这 5 个计算配对的 source 诊断 delta 均为正，方向与历史结果相同；化学改变没有产生当前模型分数的高低风险翻转。MMAF→cys-mcMMAF 的 v0.1 delta 从 0.22285 变为 0.14840，提示报告定量差值时必须写明结构处理方式。分数差不是 ER 差、ABCB1因果效应或完整 ADC 效能差，仍无 matched M0 配对真值。S-methyl-DM1、Lys-PEG4Mal-DM1 与歧义 Roth 材料继续 HOLD。", "", "## 复跑与文件", "", "```bash", "./.venv-experiments/bin/python outputs/20261009_adc_followup_validation/normalization/run_normalization_validation.py", "```", "", "现有结果目录受保护，脚本不会覆盖 run_manifest.json；如需重新执行，应复制脚本及协议到新的实验目录并正确调整项目根路径。", "", "[冻结运行记录](run_manifest.json) · [variant scores](variant_scores.csv) · [配对 delta](parent_release_variant_contrasts.csv) · [原子立体审计](source_atom_stereo_audit.csv) · [补充验证](additional_verification.json)", "", "机制依据：[RDKit 官方 TautomerEnumerator 文档](https://www.rdkit.org/docs/cppapi/classRDKit_1_1MolStandardize_1_1TautomerEnumerator.html)。运行采用 RDKit 2026.03.4、scikit-learn 1.9.0、Python 3.12.14；记录了 joblib/NumPy 2.5 的数组 shape DeprecationWarning，无 sklearn 版本不一致警告。"]
    inv = HERE / "atom_order_invariance/run_manifest.json"
    if inv.exists():
        im = json.loads(inv.read_text())
        rows += ["", "## 固定身份的随机输入顺序追加验证", "", f"13 个 source 各保留原 SMILES 与 16 个固定 seed 随机 isomeric SMILES，共 {im['input_rows']} 行。raw InChIKey 全部与各 source 相同，历史 score 全部复现。有限随机化下，标准化 key 改变的结构数 {im['key_variable_structure_count']}，两模型 8位分数改变的结构数 {im['score_variable_structure_count_rounded8']}，0.5 方向改变数 {im['threshold_direction_variable_structure_count']}，主模型 0.35/0.70 margin 类别改变数 {im['primary_margin_category_variable_structure_count']}。", "", "该测试用于 cap-hit 稳定性证据；不能证明遍历所有原子顺序或所有 tautomers 完备。详细协议、逐行结构/分数、pair 保守范围保存在 [atom_order_invariance](atom_order_invariance/PROTOCOL.md) 与 [结构汇总](atom_order_invariance/structure_invariance_summary.csv)。"]
    (HERE / "化学表示敏感性验证结论.md").write_text("\n".join(rows)+"\n",encoding="utf-8")
    print(json.dumps({"extra_checks":extra,"report":str(HERE / "化学表示敏感性验证结论.md")},ensure_ascii=False))

if __name__ == "__main__":
    main()
