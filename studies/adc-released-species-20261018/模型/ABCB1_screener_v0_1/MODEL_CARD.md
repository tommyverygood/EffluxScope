# ABCB1 substrate research screener v0.1

日期：2026-08-03

## 结论

本版本是可运行的研究筛选器，不是已验证的 ADC 专用二分类器。训练对象为冻结接口中 443 个 ABCB1 候选监督结构（291 阳性、152 阴性），输入 SHA256 为 `a1787c8ca0d02a9f142f3ca49a07c3434a2f07cd1977b8f754a55170c2da3bc5`。训练集内部直接结构标签冲突为 0。

内部选择模型：`morgan_physchem_logistic_C0.1_balanced`。主评估为 5 次重复、每次 5 折的 exact-Murcko 分组交叉验证。主门槛是 MCC ≥ 0.50 且 balanced accuracy ≥ 0.75；本轮状态为 **fail**。固定 0.5 阈值的校准 OOF 指标：MCC 0.495、BAcc 0.730、ROC-AUC 0.807、Brier 0.163。

## 模型判断层

输入是游离载荷或实际释放小分子的 SMILES。输出是 ABCB1 结构性底物倾向概率、固定 0.5 二分类、最大训练相似性和适用域状态。它不判断完整 ADC、释放效率、胞内游离暴露、耐药倍数、旁观者效应或临床疗效。

## 适用域

- 最大训练 Tanimoto < 0.230：拒判；
- 0.230–0.304：谨慎；
- ≥ 0.304：仅表示处于经验结构支持域，不表示实验端点可迁移。

## 仍未解决的限制

443 条记录全部来自 ChEMBL/Daood 同一整理来源，无法执行真正的来源留出验证；标签仍是作者整理候选而非逐行 Tier-A 双向转运真值。Roth 连续 RR 只作方向压力测试，结果见 summary，不能转成二元真值。任何 ADC 对象即使给出高概率，只要处于适用域外，都必须按 `abstain_outside_domain` 处理。

## 复现

```bash
./.venv/bin/python 数据/scripts/train_abcb1_screener_v0_1.py
./.venv/bin/python 模型/ABCB1_screener_v0_1/predict_abcb1.py --smiles 'CCO'
```
