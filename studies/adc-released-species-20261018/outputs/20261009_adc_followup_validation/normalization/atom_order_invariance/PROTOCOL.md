# 固定物种的输入原子顺序敏感性实验

本协议在本轮运行前保存。输入为上一批 13 个已准入独立 source 结构，未准入/HOLD 身份不补猜。对每个 source 保存原始 SMILES，并以 `randomSeed=20261009+i`（i 为固定结构顺序的零基编号）调用 `MolToRandomSmilesVect(..., numSmiles=16, isomericSmiles=True)`；保留实际重复，额外统计唯一输入数。

每个随机 SMILES 在标准化前重新解析。其完整 raw InChIKey 必须与该 source 原始计算 key 完全相同；不相同者标失败、不标准化、不计为固定身份实验。成功输入使用冻结 core.standardize_smiles，caps 保持 MaxTautomers=1000 / MaxTransforms=1000。记录标准化状态、tautomer status、标准化 SMILES/完整 key、两模型 estimator feature hash、按原冻结公式的 v0.1/v0.3 score、0.5 方向及主模型 0.35/0.70 margin 类别。

本轮 cap-hit 数值专用于复现/表示不变性诊断，不改变上一实验 canonical variant HOLD 判定，不赋予部署或生物学验证资格。两模型权重和 scaler/calibrator 保持不变，直接构造 features，避免再经标准化入口。

按结构统计 key、feature、score 唯一数和 score range，并分别总结 cap-hit / Completed 对照。配对仅报告每端有限输入表示所得 score 范围构成的保守 delta 区间 `[release_min-parent_max, release_max-parent_min]`，不能解释为生物学 CI、ER 或真实化学互变分布；只审计观测范围内分数方向是否稳定。三终末/下游配对与两中间体配对阶段保持区分。

验收：13×17=221 行保留；raw identity 全相同或明确失败；原输入/core/model SHA256 前后不变；原始 SMILES 重算旧分数 1e-8 内复现；同身份控制构成区间须包含 0。有限随机化不能证明所有原子顺序不变、所有 tautomers 完备或替代标准化已验证。
