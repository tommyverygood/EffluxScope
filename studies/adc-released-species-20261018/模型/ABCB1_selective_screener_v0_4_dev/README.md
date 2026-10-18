# ABCB1 三结果风险筛选器 v0.4-dev

本工具面向明确的游离载荷或释放小分子结构，输出：

- `high_ABCB1_efflux_risk`：高 ABCB1 外排风险；
- `low_ABCB1_efflux_risk`：低 ABCB1 外排风险；
- `abstain`：拒判，并给出机器可读原因。

`abstain` 是决策策略，不是训练得到的第三种生物学类别。

## 单个结构

```bash
.venv/bin/python \
  模型/ABCB1_selective_screener_v0_4_dev/predict_abcb1_selective.py \
  --smiles 'CC(=O)Nc1ccc(O)cc1'
```

可选上下文参数：

```bash
  --chemical-object-type free_payload \
  --requested-context structure_only_screening
```

如果对象是 `intact_ADC`，或请求上下文是 `clinical_DDI`、`efficacy`、`resistance` 等超出范围的任务，工具会直接拒判。批量 CSV 可使用 `chemical_object_type` 和 `requested_context` 两列；缺省时按释放/游离小分子结构筛选处理。

## 批量结构

输入 CSV 至少包含 `smiles` 列：

```bash
.venv/bin/python \
  模型/ABCB1_selective_screener_v0_4_dev/predict_abcb1_selective.py \
  --input-csv query.csv \
  --smiles-column smiles \
  --output predictions.csv
```

## 主要字段

- `decision_zh`：高风险、低风险或拒判；
- `chemical_object_type`、`requested_context`、`context_compatibility`：对象与任务边界；
- `primary_calibrated_probability`：冻结 v0.1 主模型概率；
- `bootstrap_probability_p10/p90`：按 Murcko 骨架重采样的模型区间；
- `applicability_status`：两套结构参考面板和理化包络状态；
- `abstention_reasons_json`：批量结果中的拒判原因；
- `nearest_primary_references_json`、`nearest_direct_references_json`：可追溯近邻证据；
- `claim_boundary`：每行输出的声明边界。

## 当前边界

本版本是开发期风险分诊工具，不是成熟的 ADC 专用底物分类器。它不判断完整 ADC、释放效率、胞内游离浓度、耐药、疗效或临床 DDI。历史监管池和 ASAP 2025 已在规则锁定前被观察，因此只能作事后诊断；仍需要一个新的、结果未知的匹配亲本或抑制剂救援来源完成前瞻确认。
