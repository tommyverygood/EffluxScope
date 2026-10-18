# 冻结模型的化学表示敏感性验证

本协议先于本轮计算保存。问题是：上一批 13 个独立、已准入的 source 结构在冻结标准化后发生的身份表示变化，是否影响冻结 v0.1 / v0.3 estimator 的分数或 parent–release 分数差？本轮不使用新的生物学标签，不重新训练，不调整阈值，不改变生产入口或原始化学身份。

输入为 `20261009_adc_pair_experiments/pairs/results_curated_v2/` 中的 entities、unique predictions、pair contrasts 和 release-stage classification。已有 identity HOLD 仍保留为 HOLD，不补猜结构。Lys-MCC-DM1 只沿用已准入单一表示并保留其 epimeric-mixture 限制。

对每个独立 source 结构计算以下表示：

1. `source`：解析并保留 source 立体/键型，不调用 standardizer。
2. `no_tautomer`：同旧流程的 Cleanup → FragmentParent → Uncharger，但不做 tautomer enumeration。
3. `stereo_preserved_tautomer`：相同前处理；SetRemoveSp3Stereo(False)、SetRemoveBondStereo(False)，仍选择 canonical tautomer。该设置不保证选择到 source tautomer，不补回原本未定义的立体中心。
4. `frozen_old`：逐步复刻已冻结标准化，核对保存的 standardized SMILES/InChIKey 和两模型旧分数。v0.1 若训练身份映射覆盖则另留 `frozen_primary_effective`，不能与普通标准化混淆。
5. `source_stereo_stripped`：只删除 source 立体信息，作为模型是否能感知立体的阴性对照。

特征直接从每个 variant 的 RDKit Mol 构造，按 joblib 的 feature_type/scaler/estimator 计算；v0.1 原有 logit calibrator 按同一公式应用，v0.3 原有估计器分数保持原定义。绝不通过会再次 standardize 的旧预测入口计算 variant。

校验：原输入/模型/core SHA256 前后不变；13 个独立结构旧分数在 1e-8 精度内复现；同身份 5 对在各 variant 的两端结构与分数完全相同；纯去立体对照的 features / scores 不变；canonical 变体若解析失败、enumeration 不 Completed、元素组成或电荷异常，HOLD 不计为可用分数。历史冻结分数另列为 `historical_frozen_*_not_newly_admitted`，仅供复现与 source 诊断差值，不能绕过 HOLD。全模型 v0.4 的替代表示决策不在本轮计算范围。

输出：逐阶段结构/立体审计、原 source 原子编号的丢失立体中心与键审计、variant 特征与分数、配对 delta、身份对照、复现检查、带版本与输入哈希的 run manifest、中文结论。分数差为表示敏感性，不是 ER 变化、外排实验证据或替代表示校准结果。

机制参考：[RDKit 官方 TautomerEnumerator 文档](https://www.rdkit.org/docs/cppapi/classRDKit_1_1MolStandardize_1_1TautomerEnumerator.html)。文档说明默认去除参与互变异构的 sp3 / 双键立体信息；本地实际版本与运行参数以 manifest 为准。
