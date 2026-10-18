# 数据字典

## 数据关系

`compound_registry_v1.compound_id` 是 ADC 载荷数据的主键。

`adc_release_species_registry_v2.release_species_id` 是 ADC 实际释放物/代谢物的主键。它与游离载荷 `compound_id` 是不同层级，不得因为共享 payload core 而合并。v1 保留为第二批快照，第三批后续分析默认使用 v2。

- 一个 compound 可以出现在 Roth 27 载荷面板中；
- 一个 compound 可以有多个 transporter assay records；
- 每条 assay record 只对应一个 transporter；
- Daood 训练数据保留其原始 `PID`，在完成结构标准化和去重前不强行映射到 ADC `compound_id`。
- `*_standardized_v1.csv` 是派生层；源结构字段不覆盖，所有模型键均来自标准化连接性结构。

## 核心表

### `compound_registry_v1.csv`

| 字段 | 含义 |
|---|---|
| compound_id | 本项目唯一化学实体编号 |
| paper_compound_name | 合并后的标准论文名称 |
| paper_names_as_reported | 论文不同表格中的原始名称 |
| in_roth_27_panel | 是否出现在 Table 2 的 27 载荷面板 |
| in_roth_table3_confirmatory | 是否出现在 Table 3 确认实验 |
| pubchem_cid | PubChem 化合物编号 |
| isomeric_smiles | PubChem 返回的含立体化学 SMILES |
| connectivity_smiles | PubChem 返回的连接性 SMILES |
| structure_status | resolved、unresolved、manual_review 或 salt_form_review |
| structure_note | 结构使用限制和人工审核要求 |

### `roth_2026_table3_transporter_records_long_39.csv`

| 字段 | 含义 |
|---|---|
| assay_record_id | 一条化合物—转运体实验记录的唯一编号 |
| compound_id | 对应结构注册表编号 |
| transporter_gene | ABCB1、ABCG2 或 ABCC1 |
| control_cell_line | HEK293-pcDNA |
| transporter_cell_line | 对应稳定转染细胞 |
| control_ic50_mean/sd | 对照细胞 IC50 均值和标准差 |
| transporter_ic50_mean/sd | 转运体细胞 IC50 均值和标准差 |
| resistance_ratio | 论文报告的 RR |
| binary_substrate_label | 当前留空，防止擅自设置阈值 |
| evidence_tier | 当前为 B |

### `external_truth_set_v1_pgp_bcrp_continuous_26.csv`

从 Table 3 长表中仅保留 ABCB1 和 ABCG2。它是第一版 ADC 外部验证核心表，不含 ABCC1，也没有二分类标签。

### `daood_2025_tableS2_all_endpoints_8794.csv`

Daood 2025 补充数据 Table S2 的原样 CSV 转换。每行一个化合物，8 个转运端点分列；空值表示该端点没有标签，不表示阴性。

### `daood_2025_pgp_bcrp_substrate_training_long.csv`

| 字段 | 含义 |
|---|---|
| training_record_id | 本项目训练记录编号 |
| pid | Daood Table S2 的原始 PID |
| source_database | ChEMBL、PubChem、Metrabase 或 Shaikh et al. |
| smiles_as_published | Daood 提供的 SMILES，尚未二次标准化 |
| transporter_gene | ABCB1 或 ABCG2 |
| published_binary_label | Daood 数据中的 0/1 标签 |
| processing_status | 当前清洗和去重状态 |

### `adc_release_species_registry_v2.csv`

| 字段 | 含义 |
|---|---|
| release_species_id | 释放物/代谢物的唯一编号 |
| payload_core | 对应的游离载荷核，只用于分组，不代表化学实体相同 |
| residual_group | 释放后仍保留的氨基酸、linker 或代谢基团 |
| release_identity_status | 释放物身份的实验确认程度 |
| structure_status | 结构是否可用于建模，与身份证据分开记录 |
| existing_payload_compound_id | 与第一批游离载荷核的关联键，不是实体合并键 |
| modeling_status | 结构模型纳入、排除或审核状态 |

`regulatory_epimeric_mixture` 表示监管登记材料是 epimeric mixture，表内单一立体异构 SMILES 只能作为代表结构，不能在无处理规则时直接当作材料本身。

### `adc_release_observations_long_v2.csv`

每行是一条“ADC 格式—细胞/体系—释放物—时间/剂量—观测”。`measured_value` 只保留来源明确报告的数值；定性身份证据写在 `observation`，不伪造数值。

### `adc_release_species_transporter_evidence_long_v2.csv`

| 字段 | 含义 |
|---|---|
| transport_record_id | 释放物—转运体实验记录编号 |
| source_assigned_call | 来源作者/监管文件的原始解读 |
| source_binary_label | 来源明确给出的底物解读，不等于模型标签 |
| model_binary_label | 当前故意留空；需预注册证据映射后才生成 |
| evidence_tier | 直接双向转运/抑制剂逆转、定性直接实验或比较性细胞证据 |
| truth_set_candidate | 是否可以进入后续真值审核，不代表已经进入训练/测试集 |
| assay_endpoint | 原始实验端点类型，例如双向转运、细胞毒 IC50 或监管来源结论 |
| control/transporter/inhibitor condition/value | 保留功能性证据中的原始比较条件和值 |
| effect_ratio | 来源报告的功能性比值；仅在同一记录内部解释 |
| effect_ratio_definition | 比值分子、分母的明确定义 |

### `free_payload_release_pair_evidence_v1.csv`

| 字段 | 含义 |
|---|---|
| pair_id | 游离载荷—释放物成对证据编号 |
| free/release endpoint | 两侧各自的实验端点，不假设相同 |
| descriptor delta | 释放物减游离载荷的 MW、XlogP、TPSA 描述符差值 |
| directional_comparison | 仅用于生成方向性假设，不是模型标签 |
| endpoint_compatibility | 当前均为 `not_comparable_cross_assay` |
| causal_interpretation | 当前均为 `hypothesis_only_not_causal` |

v2 更新了 DM1—Lys-MCC-DM1 的直接双向转运近配对，并增加 DXd 的 ABCB1/ABCG2 独立方向证据；v1 保留第三批快照。

### 第五批直接与匹配证据表

| 表 | 主键/用途 | 关键边界 |
|---|---|---|
| `direct_bidirectional_transport_numeric_v1.csv` | `transport_record_id`；直接 Papp/ER 数值长表 | Figure 8 数值行必须按 `digitization_status` 排除真值集 |
| `foti_figure8_digitization_v1.csv` | `digitization_record_id`；像素级复现 | `pass1/pass2` 是算法校准，不是双人独立整理 |
| `matched_evidence_grade_definitions_v1.csv` | `matched_grade`；定义 M0-M4/M1H | 只有 M0 允许在预注册同端点内计算配对变化 |
| `matched_context_comparisons_v1.csv` | `comparison_id`；游离载荷—释放物可比性 | `cross_assay_fold_computed` 当前全部为 no |
| `paired_release_species_mechanistic_evidence_v1.csv` | `mechanistic_record_id`；同体系机制证据 | 不能转成 isolated transport label |
| `patent_mmaf_verapamil_functional_v1.csv` | `functional_record_id`；专利 Table 7 原始 IC50 | verapamil 10/40 uM 冲突，且 ADC 行不是释放物真值 |
| `adc_release_species_transporter_evidence_long_v3.csv` | `transport_record_id`；释放物证据 v3 | `source_binary_label` 仍不等于模型标签 |

`matched_grade` 定义：

- M0：同一对化学实体、同一转运体细胞背景、同一实验终点、同一浓度/扰动框架；
- M1：同一精确化学实体和转运体的独立直接实验，但不是同一实验批次；
- M1H：疑似同一实体，但材料身份或直接数值仍缺失，必须 hold；
- M2：同为直接双向转运终点，但细胞背景、浓度或实验室不同；
- M3：同细胞/同 ADC 体系的功能或细胞内保留对照，但并非两种孤立化学实体的直接转运；
- M4：跨终点方向性或同来源定性比较。

`direct_bidirectional_transport_numeric_v1.csv` 中 `truth_set_eligible=yes_continuous_no_binary_cutoff` 只表示可作为连续直接证据审核，不代表已有二分类标签；任何阈值仍需单独预注册。

### 第六批原始数值与身份审计表

| 表 | 主键/用途 | 关键边界 |
|---|---|---|
| `dxd_maaa1181a_caco2_transport_numeric_v1.csv` | `record_id`；4 个 1 uM DXd/MAAA-1181a Caco-2 抑制剂条件 | 只报告 ER；共表达体系不能拆成单基因 ER 真值 |
| `direct_bidirectional_transport_numeric_v2.csv` | `transport_record_id`；第五批 23 行加 DXd 4 行 | 不同实体、体系和端点仍不直接合并 |
| `adc_release_species_transporter_evidence_long_v4.csv` | `transport_record_id`；释放物证据 v4 | DXd 记录的 `target_attribution` 是多转运体，不生成模型二分类标签 |
| `roth_deruxtecan_identity_audit_v1.csv` | `audit_id`；Roth 试剂元数据逐项审核 | 语义称呼、数据库 exact-name 命中都不等于实验材料身份证明 |
| `external_truth_set_v2_pgp_bcrp_continuous_identity_audited.csv` | `assay_record_id`；Roth 连续 RR 身份审计版 | 原始连续记录保留，身份不明记录仅从结构指标排除，不删除来源证据 |
| `baseline_roth_external_predictions_v2_identity_audited.csv` | `assay_record_id` + `model_name`；预测有效性审计 | `prediction_validity=invalid_for_material_specific_interpretation` 时不得用于外部指标或材料解释 |
| `baseline_roth_external_metrics_v2_identity_safe.csv` | transporter + model；剔除身份不明记录后的 Spearman | 每基因 n=12、小样本、无阈值，只是方向性压力测试 |
| `roth_deruxtecan_structure_sensitivity_v1.csv` | transporter + model + structure hypothesis；C26/C52 双假设 | 所有行均要求身份 abstention，不能挑选更有利的结构 |
| `roth_deruxtecan_structure_sensitivity_delta_v1.csv` | transporter + model；C26-C52 概率/相似性差 | 衡量输入身份敏感性，不衡量预测准确性 |
| `lys_mcc_dm1_primary_table_audit_v1.csv` | `audit_id`；Rock/Foti 溯源 | 未取得 Supplemental Table 2 时不得升级为作者原始数值 |
| `pubchem_assay_availability_recheck_20260728.csv` | `aid`；CSV 端点状态 | HTTP 404 是可得性记录，不表示实验为阴性 |
| `baseline_model_gate_summary_v2_identity_safe.csv` | transporter + model；claim gate v2 | 当前 6 个组合全部仅限开发基线 |

`dxd_maaa1181a_caco2_transport_numeric_v1.csv` 关键字段：

| 字段 | 含义 |
|---|---|
| efflux_ratio | 来源报告的 Papp B-to-A / Papp A-to-B；不反推未给出的单向 Papp |
| baseline_to_condition_er_reduction_fold | 无抑制剂 ER / 当前条件 ER，仅作同面板描述 |
| target_attribution | 表明是 P-gp/BCRP 共存体系或抑制剂条件，不把 ER 唯一归因到单一基因 |
| truth_set_use | 当前固定为连续多转运体外部证据，不进入基因特异训练真值 |

身份审计规则：Roth 的 `Deruxtecan` 只有在材料元数据把测试对象唯一指向 C26 DXd 或 C52 linker-payload 后，才能恢复结构外部评估。PubChem CID 117888634 与 118305111 只定义两个候选结构，不解决论文试剂身份。

### 第七批释放载荷直接转运与同体系上下文表

| 表 | 主键/用途 | 关键边界 |
|---|---|---|
| `obi992_caco2_payload_transport_numeric_v1.csv` | `record_id`；Table S5 的 16 个 Papp/SD/ER 条件 | `transporter_arm` 是抑制剂归因，不是单基因细胞系；抑制剂浓度未报告 |
| `pmda_mmae_caco2_transport_numeric_v1.csv` | `record_id`；Study 8004 与 RPT-01709 的 15 个原始条件 | `uncertainty_type` 区分 SD 与 SEM；相同浓度的不同 run 不合并 |
| `mmae_mdck_mdr1_text_numeric_v1.csv` | `record_id`；正文精确 ER 与 ratio-of-ratios | Papp 柱形图未数字化；ER 是作者正文数值 |
| `bay716_caco2_transport_numeric_v1.csv` | `record_id`；BAY-716 原始 Papp/ER | 结构未闭环；nm/s 同时换算为 `10^-6 cm/s` |
| `zd06519_payload_transporter_calls_v1.csv` | `record_id`；来源 check/X/marginal 调用 | 无数值、无项目二分类标签 |
| `obi992_hcc827_payload_adc_resistance_v1.csv` | `record_id`；同背景载荷与完整 ADC IC50/RR | `entity_type=complete_ADC` 的行禁止当成载荷转运真值 |
| `obi992_adc_inhibitor_rescue_v1.csv` | `record_id`；完整 ADC 抑制剂逆转 | 只作 ADC 层功能证据 |
| `same_system_payload_adc_context_v1.csv` | `comparison_id`；载荷—完整 ADC 同体系桥 | `adc_to_payload_rr_ratio` 仅作上下文描述，不解释为转运动力学 |
| `s_methyl_dm1_cellular_retention_v1.csv` | `record_id`；细胞滞留/外排与再结合 | `transporter_specific=no`，进入膜分配机制层而非转运真值层 |
| `adc_release_species_registry_v3.csv` | `release_species_id`；12 个释放物实体 | 新增 exatecan 和 BAY-716；BAY-716 必须结构 abstain |
| `direct_bidirectional_transport_numeric_v3.csv` | `transport_record_id`；61 行统一直接数值表 | RPT-01709 的 SEM 保存在专表，不写入继承的 `*_sd` 字段 |
| `adc_release_species_transporter_evidence_long_v5.csv` | `transport_record_id`；34 条释放物证据 | `model_binary_label` 全空；来源二分类解读不自动迁移 |
| `matched_context_comparisons_v3.csv` | `comparison_id`；9 组可比性审计 | MMAE 首次达到 M1；M0 仍为 0 |
| `seventh_batch_validation.json` | 自动一致性检查 | 校验 ID 唯一、计数、模型标签为空及原始文件 SHA256 |

`obi992_caco2_payload_transport_numeric_v1.csv` 关键字段：

| 字段 | 含义 |
|---|---|
| transporter_gene / transporter_arm | 抑制剂设计所针对的 ABCB1/P-gp 或 ABCG2/BCRP；不是基因隔离体系 |
| papp_ab_mean_1e6_cm_s / papp_ba_mean_1e6_cm_s | 来源表列的双向平均表观通透系数，单位 `10^-6 cm/s` |
| papp_ab_sd / papp_ba_sd | 来源表列 SD |
| baseline_to_condition_er_reduction_fold | 同一化合物、同一 target arm 的 ER(no inhibitor)/ER(inhibitor) |
| truth_set_use | 固定为连续外部归因审计，不进入单基因训练真值 |

`pmda_mmae_caco2_transport_numeric_v1.csv` 中 `Study 8004` 为 mean(SD)，`RPT-01709` 为 mean ± SEM、n=3。统一 v3 表沿用旧版 `papp_*_sd` 字段，因此只将真正 SD 写入；SEM 只保存在本专表，避免把 SEM 伪装成 SD。

`matched_context_comparisons_v3.csv` 的 `MATCH_B7_001` 是 M1，不是 M0：两篇独立研究都测相同 MMAE 与 ABCB1 方向性，且一个含 MDCK-MDR1 单转运体对照，但它们并非同一实验 run/背景。M1 证明可复现方向，不证明游离载荷与带残基释放物之间的因果改变。

### 第八批 SG3199 与被动膜通透表

| 表 | 主键/用途 | 关键边界 |
|---|---|---|
| `khera_figure4b_pampa_digitization_v1.csv` | `digitization_record_id`；9 个 Figure 4B PAMPA 柱平均值 | 单图算法数字化；误差棒为空；`truth_set_eligible=no` |
| `adc_release_species_passive_permeability_v1.csv` | `release_species_id`；6 个已注册释放物的 PAMPA 映射 | 无转运体的被动膜通透层，不是 P-gp/ABCG2 真值 |
| `adc_release_species_registry_v4.csv` | `release_species_id`；13 个释放物实体 | 新增 SG3199，明确映射已有 `ADC_PAYLOAD_032` 结构 |
| `adc_release_species_transporter_evidence_long_v6.csv` | `transport_record_id`；39 条分层证据 | SG3199 的 FDA 来源标签与 Roth 功能 RR 分存；`model_binary_label` 全空 |
| `sg3199_cross_endpoint_conflict_v1.csv` | `conflict_record_id`；ABCB1/ABCG2 审计 | ABCG2 为明确上下文冲突，不做真值裁决 |
| `matched_context_comparisons_v4.csv` | `comparison_id`；11 组匹配审计 | SG3199 两行为 M4 跨终点比较；M0 仍为 0 |
| `eighth_batch_validation.json` | 自动一致性检查 | 检查计数、冲突保留、空模型标签、SHA256 和直接数值不虚增 |

`khera_figure4b_pampa_digitization_v1.csv` 关键字段：

| 字段 | 含义 |
|---|---|
| `x_interval_pixels` / `bar_top_*_pixel` | 220 dpi 全页 PNG 上每根柱的检测区间和柱顶坐标 |
| `pampa_pass1_1e6_cm_s` | 7 个 y 轴刻度最小二乘校准的数字化均值 |
| `pampa_pass2_1e6_cm_s` | 0 和 3.0 端点校准的数字化均值 |
| `pampa_consensus_1e6_cm_s` | 两种算法校准的均值；不代表双人独立整理 |
| `one_pixel_value_uncertainty_1e6_cm_s` | y 方向一个源像素对应的轴值，只是像素级尺度 |
| `error_bar_numeric_value` | 故意留空；未对误差棒做单人数字化 |
| `truth_set_eligible` | 固定为 `no`；PAMPA 不是基因特异转运真值 |

`sg3199_cross_endpoint_conflict_v1.csv` 中 `directional_relationship` 只记载两个来源的方向关系。`discordant_context_dependent` 不是对任一来源判错，而是说当前细胞背景、浓度和终点不同，无法合法生成唯一项目标签。

### 第十批训练标签来源审计表

| 表 | 主键/用途 | 关键边界 |
|---|---|---|
| `transporter_label_provenance_audit_v1.csv` | `audit_record_id`；1,114 条训练标签逐行审计 | `published_binary_label` 原样保留；`provenance_tier` 不能解释为来源作者给出的质量等级 |
| `transporter_label_source_matrix_v1.csv` | transporter + source；来源类别组成 | 基于 1,081 条泄漏清洁记录；单类来源不能计算平衡二分类外部指标 |
| `transporter_truthset_AB_candidate_v1.csv` | `audit_record_id`；512 条作者人工整理候选 | 名称中的 truthset 表示候选工作集；Tier A 仍为 0 |
| `transporter_auxiliary_labels_v1.csv` | `audit_record_id`；602 条辅助/PU/排除记录 | 不与 Tier B 等权进入主监督训练；外部泄漏排除继续有效 |
| `transporter_label_conflicts_v1.csv` | `issue_id`；7 个开放审计问题 | 包括来源单类、逐行 assay 缺失、无跨来源共识和结构 QC 限制 |
| `source_holdout_diagnostics_v1.csv` | transporter + held-out source；来源留出诊断 | 训练或测试单类时 MCC/BAcc/AUC 留空，不用 accuracy 冒充二分类泛化 |
| `tier_AB_scaffold_cv_diagnostics_v1.csv` | transporter；A/B 候选 5 折 scaffold CV | 固定 Morgan/Logistic/seed；项目门槛不是行业标准 |
| `targeted_label_data_gap_queue_v1.csv` | `gap_id`；P0/P1 定向补数任务 | 接收/拒绝标准在采集前固定，防止为改善模型分数而改标签 |
| `tenth_batch_label_audit_summary.json` | 批次机器可读摘要 | 明确 `claim_boundary`：不能从 Daood Table S2 重建 Tier A |

`transporter_label_provenance_audit_v1.csv` 关键字段：

| 字段 | 含义 |
|---|---|
| `published_binary_label` | Daood 发表的原始二分类调用；审计不修改 |
| `source_selection_role` | Daood 方法中该来源进入数据集的作用，例如主来源或类别补充 |
| `source_label_confounding_flag` | 同一 transporter 下该来源是否只贡献一个类别 |
| `row_level_primary_*_available` | Table S2 能否逐行闭环 assay ID、原始引用或连续数值；当前均为 `no` |
| `provenance_tier` | 本项目证据完整度等级；B 为作者人工整理但缺逐行实验，C 为数据库/二次汇编 |
| `strict_AB_candidate` | 是否允许进入 Tier A/B 候选监督基线；仍受 leakage action 约束 |
| `model_use_role` | 主监督候选、辅助/PU、外部泄漏排除或待映射 |
| `negative_evidence_basis` | 阴性标签来自何种来源调用；不把空值解释为阴性 |

第十批准入规则：只有逐行可追溯的直接转运数值、实验系统、条件/对照和原始引用才能进入 Tier A。作者人工整理但在补充表无法逐行复核的标签最多为 Tier B；数据库汇编和二次文献汇编为 Tier C。所有层级均继续与 ADC 释放物连续证据层分开。

### 第十一批 ChEMBL 逐行实验回填表

| 表 | 主键/用途 | 关键边界 |
|---|---|---|
| `chembl_candidate_activity_hits_v1.csv` | `backfill_record_id`；2,565 条活性命中 | API 行是原文定位线索，`tier_A_auto_upgrade` 全为 `no` |
| `chembl_candidate_assay_registry_v1.csv` | `assay_chembl_id`；1,342 个 assay | 保留 organism、confidence、BAO 和文献关联 |
| `chembl_candidate_document_registry_v1.csv` | `document_chembl_id`；520 篇 document | DOI 只解决可追溯性，不自动证明实验合格 |
| `chembl_candidate_provenance_triage_v1.csv` | `audit_record_id`；512 条候选 | 保留发表标签，只生成复核状态 |
| `chembl_reconstructed_binary_label_candidates_v1.csv` | `reconstruction_id`；181 条一致候选 | `strict_tier_A_training_use=no_pending_primary_fulltext_QC` |
| `chembl_reconstructed_label_conflicts_v1.csv` | `conflict_id`；36 条冲突 | 严格二分类训练排除，不做多数投票覆盖 |
| `eleventh_batch_manual_verified_binary_labels_v1.csv` | `manual_label_id`；7 条人工核验 | 2 条 A1 进严格核心；5 条 A2 仅降权敏感性 |
| `eleventh_batch_binary_training_supplement_v1.csv` | `manual_label_id`；建模接口 | 字段与人工核验表相同，保留证据权重 |
| `eleventh_batch_primary_source_queue_v1.csv` | `queue_id`；474 条复核队列 | P0 优先 ABCG2 阴性和 ABCB1 阳性 |
| `eleventh_batch_manual_primary_source_registry.csv` | `source_id`；4 个原始文件 | URL、本地路径、文件大小和 SHA256 同时保存 |
| `eleventh_batch_raw_download_log.csv` | `raw_file`；原始下载日志 | ChEMBL 分页 JSON 和人工 PDF/XML 共用溯源格式 |
| `eleventh_batch_chembl_backfill_summary.json` | 批次机器可读摘要 | 新增是证据可信度，不是新结构样本数 |

`chembl_candidate_activity_hits_v1.csv` 关键字段：

| 字段 | 含义 |
|---|---|
| `candidate_match_type` | 精确 molecule ChEMBL ID 或 parent molecule ID 匹配 |
| `evidence_stratum` | `direct_transport_candidate`、`inhibitor_or_probe_activity`、`functional_indirect_activity` 或 `other_target_activity` |
| `explicit_direct_substrate_call` | 仅在直接底物上下文中保留明确 0/1；数值或模糊结果留空 |
| `explicit_call_basis` | 得出明确调用的文本规则或留空原因 |
| `tier_A_auto_upgrade` | 固定为 `no`；防止把 API 重叙述当成原文真值 |

`eleventh_batch_manual_verified_binary_labels_v1.csv` 关键字段：

| 字段 | 含义 |
|---|---|
| `evidence_grade` | A1 为全文直接实验；A2 为原始摘要加可追溯 assay 行 |
| `transporter_comparison_control` | 亲本/空载体、mock、ATP/AMP 或选择性抑制剂等对照 |
| `chembl_primary_discrepancy` | ChEMBL 与原文的条件差异；Cefadroxil 的 6 h/6 min 差异在此保留 |
| `training_use` | `yes_strict_binary_core` 或 `yes_sensitivity_lower_weight` |
| `weight_recommendation` | 严格 A1 为 1.0，A2 建议为 0.75；是预注册建议，不是学习得到的最优权重 |

第十一批准入规则：严格核心需全文核对精确实体、人转运体、直接底物终点、转运体对照、实验 QC 和明确结果。A2 不与 A1 等权；36 个冲突对象不在本批裁决为唯一全局标签。

### 第十三、十四批原文复核与 evidence overlay

| 表 | 主键/用途 | 关键边界 |
|---|---|---|
| `targeted_primary_review_abcb1_top30_abcg2_negative_weak_v1.csv` | `targeted_review_id`；ABCB1 前30与ABCG2优先队列共42条 | A2 仅 `external_sensitivity_only`；`not_upgraded` 和身份受限记录不生成标签 |
| `targeted_primary_review_source_registry_v1.csv` | `source_id`；24个 PubMed/PMC/既有原始文件 | URL、本地文件、访问日期、大小和 SHA256 |
| `targeted_primary_review_interface_v1.csv` | `interface_id`；42条复核与标准化候选结构连接 | `review_binary_label` 全部为空；冻结 published label 不被覆盖 |
| `fourteenth_batch_primary_verified_binary_labels_v1.csv` | `manual_label_id`；3条 exact ABCG2 A2敏感性记录 | MW-150、COMT compound 28、exact itraconazole；本批不新增 strict core |
| `fourteenth_batch_primary_review_queue_resolution_v1.csv` | `audit_record_id`；10条优先 ABCG2 队列动作 | 数值 0.9、ChEMBL Not active、身份受限 series 继续保留空标签 |
| `binary_supplement_training_interface_v2_primary_overlay.csv` | `interface_id`；3条 T14 provenance overlay | 不增加重复结构；标签分布不变 |
| `binary_supplement_model_diagnostics_v2_primary_overlay.csv` | `transporter_gene × model × holdout_type` | 同一候选配方下重跑；仍不是 claim-gate 通过 |

T13/T14 的 A2 仅用于外部敏感性、索取全文/SI 和数据优先级；不把未报告 transport、抑制剂/ATPase、系列级“not transported”或溶解度不足改写为二分类阴性。

### ChEMBL direct-candidate locator 补充层

| 表 | 主键/用途 | 关键边界 |
|---|---|---|
| `chembl_direct_candidate_locator_v1.csv` | `locator_id`；987 条 ChEMBL direct-transport candidate locator | 仅保留定位字段、数值/文本结果和 assay 描述；`binary_label_proposal` 全空，`training_use=no_locator_only` |
| `chembl_direct_candidate_locator_summary_v1.json` | locator 行数、ABCB1/ABCG2 构成、端点兼容性和源表 SHA256 | ATPase、抑制/结合、transport 描述分层；不把数值自动变为 binary label |
| `fourteenth_batch_unresolved_numeric_locator_overlay_v1.csv` | 第十四批 5 个未决 ABCG2 对象与 ChEMBL activity/assay 行的连接 | ER=0.9 和 `Not active` 只作为待全文复核的 locator，`overlay_action=retain_blank_no_binary_label_transfer` |
| `fourteenth_batch_unresolved_primary_access_audit_v1.csv` | 5 个未决对象的 DOI、访问状态、ChEMBL 行和原始文件哈希 | ACS 403/Elsevier metadata-only 说明数据不可升级，不是阴性结论 |
| `fourteenth_batch_unresolved_primary_access_audit_summary_v1.json` | 访问状态计数和训练用途计数 | `binary_label_nonblank=0` 必须保持 |

### 第十六批独立 ABCB1/P-gp 来源与泄漏审计层

| 表 | 主键/用途 | 关键边界 |
|---|---|---|
| `public_independent_bikadi_2011_abcb1_table_s1_v1.csv` | `independent_record_id`；PLOS 2011 Table S1 的 195 条逐行 class、training/external set、assay 摘要、`assay_family`、reference 与身份字段 | 1=source substrate、0=source non-substrate；来源为异质文献整理，不能视为 Tier-A row-level direct transport truth；`strict_training_default=no` |
| `public_independent_bikadi_2011_abcb1_overlap_audit_v1.csv` | `independent_record_id`；与 Daood 标准化训练记录的 exact InChIKey、connectivity block 和标签一致性 | 102 条无 Daood 结构重叠可作来源留出；2 条跨来源冲突（Vindoline、Lovastatin）必须 hold，不覆盖冻结标签 |
| `public_independent_bikadi_2011_source_registry_v1.csv` | `source_id`；官方 PLOS DOCX URL、文件大小、SHA256、解析行数和允许用途 | 记录正文“197”与附件实际“195”计数差异；不凭空补行 |
| `public_independent_bikadi_2011_summary_v1.json` | 解析/身份/重叠/标签转移摘要 | `binary_label_transfer=none`；冻结接口前后行数均 1,137 |
| `public_independent_bikadi_model_diagnostics_v1.csv`、`public_independent_bikadi_model_predictions_v1.csv` | 102 条无重叠结构的 ABCB1 来源留出敏感性评分 | 与训练标签完全隔离；仅报告 source-heldout sensitivity，不能用来重调阈值或宣称 gate 通过 |
| `collect_plos_pgp_2011.py`、`run_bikadi_external_sensitivity.py` | DOCX 解析、PubChem CID/名称身份解析、重叠审计和冻结配方外部评分 | PubChem 仅用于结构身份与泄漏审计；名称解析标为 `ok_by_name`，不等同于原始 DrugBank ID 确认 |
| `beery_2012_chlorothiazide_abcg2_primary_observations_v1.csv` | `observation_id`；5 条 Beéry 2012 ABCG2 连续观察（MDCKII-BCRP ER、Caco-2 ER/Ko143 rescue、vesicle Km、ATPase EC50） | `model_binary_label` 全空；MDCKII-BCRP author-text ER 约 36 也只作为连续观察，Caco-2 与 ATPase 不进入 gene-specific binary target |
| `beery_2012_chlorothiazide_identity_overlap_audit_v1.csv` | chlorothiazide 标准化身份与 Daood/AB candidate 的 exact overlap | 命中 1 条既有 Metrabase 辅助记录、0 条 A/B candidate；`action=retain_continuous_external_overlay_no_binary_transfer` |
| `beery_2012_chlorothiazide_source_registry_v1.csv`、`beery_2012_chlorothiazide_summary_v1.json` | J-STAGE 官方 PDF 的 URL、大小、SHA256、端点计数和使用边界 | 只作 primary continuous evidence；不根据 ER/Km/ATPase 自动生成二分类标签 |
| `zhao_2009_bcrp_ortholog_primary_observations_v1.csv` | `observation_id`；4 个化合物的 MDCKII-Bcrp1 parent/chrysin Papp、SD、ER 与 mass recovery | `transporter_species=murine_Abcg2_Bcrp1`；正向/救援行只作 ortholog continuous evidence，LY2228820 结构保持 identity hold |
| `zhao_2009_bcrp_ortholog_identity_overlap_audit_v1.csv` | 4 个 Zhao 实体与 Daood/AB candidate 的 exact identity overlay | 3 个与既有 Daood auxiliary 命中、0 个 A/B candidate；不覆盖标签 |
| `zhao_2009_bcrp_ortholog_source_registry_v1.csv`、`zhao_2009_bcrp_ortholog_summary_v1.json` | 来源 PDF SHA256、4/8 行计数、species 和训练边界 | `binary_label_transfer=none`；不能当作人 ABCG2 strict truth |
| `shaikh_2009_bcrp_ortholog_primary_observations_v1.csv` | `observation_id`；4 个化合物的 MDCK WT/Bcrp1/Ko143 条件、Papp/SD/ER、TEER 和 mannitol QC | `transporter_species=murine_Abcg2_Bcrp1`；abacavir/zidovudine 为新结构但仅 ortholog external evidence |
| `shaikh_2009_bcrp_ortholog_identity_overlap_audit_v1.csv` | Shaikh/Giri 4 个实体与 Daood/AB candidate 的 exact identity overlay | imatinib/prazosin 命中现有 ChEMBL candidate；没有标签迁移 |
| `shaikh_2009_bcrp_ortholog_source_registry_v1.csv`、`shaikh_2009_bcrp_ortholog_summary_v1.json` | 来源 PDF SHA256、12 行计数、species、Ko143 rescue 和训练边界 | `binary_label_transfer=none`；不能当作人 ABCG2 strict truth |
| `haslam_2011_ciprofloxacin_human_abcg2_primary_observations_v1.csv` | `observation_id`；4 条 WT/human-BCRP/murine-bcrp1/Caco-2 ER 与 6 条三位人空肠供体 +/-Ko143 Papp/ER | human BCRP-MDCKII ER=0.8 为直接弱/非外排连续对照；人组织行是多转运体 context；`model_binary_label` 全空 |
| `haslam_2011_ciprofloxacin_human_abcg2_identity_overlap_audit_v1.csv` | ciprofloxacin 标准化身份与 Daood/AB candidate 的 exact overlap | 1 条既有 Daood ABCG2 auxiliary、0 条 A/B candidate；`retain_human_numeric_external_overlay_no_binary_transfer` |
| `haslam_2011_ciprofloxacin_human_abcg2_source_registry_v1.csv`、`haslam_2011_ciprofloxacin_human_abcg2_summary_v1.json` | Haslam 2011 官方开放 PDF URL、大小、SHA256、物种/端点计数和使用边界 | 人源弱对照可作 external sensitivity；鼠源/Caco-2/人组织不自动生成 gene-specific binary label |
| `xiao_2006_human_abcg2_abstract_observations_v1.csv` | `observation_id`；4 条 human MDCKII/ABCG2 摘要级 mean ER（prazosin、SN-38、CI1033、mitoxantrone） | `evidence_tier=A3_primary_abstract_numeric`；只保留连续阳性外部参照，full text/parental ER/SD/QC 不可见，`model_binary_label` 全空 |
| `xiao_2006_human_abcg2_abstract_identity_overlap_audit_v1.csv` | 4 个摘要实体与 Daood/AB candidate 的 exact identity overlay | 2 条 Daood exact、1 条 A/B candidate exact；不覆盖原标签 |
| `xiao_2006_human_abcg2_abstract_source_registry_v1.csv`、`xiao_2006_human_abcg2_abstract_summary_v1.json` | ACS 摘要 URL、DOI、4 行计数和 abstract-only 边界 | 不把摘要数值当成完整 Tier-A 表格；`binary_label_transfer=none` |
| `chen_2021_human_abcg2_cannabinoid_observations_v1.csv` | `observation_id`；9 条 human BCRP-MDCKII 直接 ER/AUC ratio（11-OH-THC、THC-COOH、THC、prazosin QC）及 Ko143 rescue | 11-OH-THC 为非底物对照，THC-COOH 为弱底物；THC 因胞内滞留没有 ER；TEER 和 parent/Ko143 条件保留；`model_binary_label` 全空 |
| `chen_2021_human_abcg2_cannabinoid_identity_overlap_audit_v1.csv` | 4 个 cannabinoid/QC 实体与 Daood/AB candidate 的 exact identity overlay | 1 条 Daood exact、1 条 A/B candidate exact；不覆盖标签 |
| `chen_2021_human_abcg2_cannabinoid_source_registry_v1.csv`、`chen_2021_human_abcg2_cannabinoid_summary_v1.json` | Chen 2021 官方开放 PDF URL、大小、SHA256、端点计数和使用边界 | 直接人源连续证据和 QC；禁止把 ER/AUC 自动迁移为二分类 truth |
| `blanco_2022_human_abcg2_albendazole_observations_v1.csv` | `observation_id`；12 条 parent/murine/human/Ko143 比较观察，ABZSO2 与 ABZSO2-NH2 的双向 Papp、SD、ER | 两个新结构的 human ABCG2 直接阳性及 Ko143 rescue；`model_binary_label` 全空；物种层不可混合 |
| `blanco_2022_human_abcg2_albendazole_identity_overlap_audit_v1.csv` | 两个 albendazole metabolite 与 Daood/AB candidate 的 exact identity overlay | Daood/A-B candidate 均 0 exact，适合 source/scaffold 外部连续核查 |
| `blanco_2022_human_abcg2_albendazole_source_registry_v1.csv`、`blanco_2022_human_abcg2_albendazole_summary_v1.json` | Blanco 2022 官方开放 PDF URL、大小、SHA256、物种/端点计数和使用边界 | comparative human/murine continuous evidence；不自动生成 binary label |
| `bleasby_2021_human_abcg2_islatravir_observations_v1.csv` | `observation_id`；5 条 human BCRP-MDCKII/parent/Ko143 双向 ER（islatravir 与 prazosin QC） | islatravir ER 4.3±1.1→0.9±0.1；`model_binary_label` 全空；直接连续外部层 |
| `bleasby_2021_human_abcg2_islatravir_identity_overlap_audit_v1.csv` | islatravir/prazosin 与 Daood/AB candidate 的 exact identity overlay | islatravir 为新结构；prazosin 重叠；不覆盖标签 |
| `bleasby_2021_human_abcg2_islatravir_source_registry_v1.csv`、`bleasby_2021_human_abcg2_islatravir_summary_v1.json` | Bleasby 2021 开放 PDF URL、大小、SHA256、5 行计数和使用边界 | human direct ER/Papp；不自动生成 binary label |
| `human_abcg2_external_source_manifest_v1.csv` | `source_id`；28 个人/鼠/组织/监管来源 ID 的 198 条观察统计，其中 155 条 human-direct | source-level external continuous/context benchmark；不把来源数或行数当作独立 binary n；`model_label_nonblank=0` |
| `human_abcg2_external_entity_overlay_v1.csv` | `standardized_inchikey`；41 个有结构实体的来源、物种、端点、重复/泄漏和 `human_external_comparability_class` 审计 | 同结构多来源/多端点只作重复证据；14 个 Daood exact、3 个 A/B candidate exact；所有 `binary_training_eligible=no` |
| `human_abcg2_external_overlay_summary_v1.json` | 外部 manifest 的行数、human-direct 计数、实体计数和标签空值摘要 | 155 条 human-direct 连续/定性记录、55 个结构实体、`model_label_nonblank=0`；不是二分类训练集 |
| `lentzas_2024_human_abcg2_accumulation_context_v1.csv` | `observation_id`；16 条 human BCRP accumulation-only 的 elacridar/tariquidar IC50/NA 观察，8 个实体 | IC50/NA 不是 ER；loperamide/verapamil/ibrutinib 的 NA 只表示无可拟合值；`model_binary_label` 全空 |
| `lentzas_2024_human_abcg2_accumulation_identity_overlap_audit_v1.csv` | 8 个 Lentzas 实体与 Daood/AB candidate 的 exact identity overlay | 2 条 Daood exact、0 条 A/B candidate；仅作 accumulation context，不覆盖标签 |
| `lentzas_2024_human_abcg2_accumulation_source_registry_v1.csv`、`lentzas_2024_human_abcg2_accumulation_summary_v1.json` | Lentzas 2024 开放 PDF URL、大小、SHA256、IC50/NA 计数和使用边界 | 单独端点层；不能替代双向 Papp/ER，也不能自动生成 binary label |
| `bora_2025_human_bbb_abcg2_context_observations_v1.csv` | `observation_id`；16 条 hCMEC/D3 人脑微血管 context 观察（15 条数值 Papp/ER、1 条无数值定性记录） | 内源人 BCRP、Ko143 1 µM；DCLB/OXC 救援、VPA/LTG 弱方向性、CBZ/TPM 无清楚信号；`model_binary_label` 全空 |
| `bora_2025_human_bbb_abcg2_context_identity_overlap_audit_v1.csv` | 9 个 Bora 实体的标准化结构与 Daood/A-B candidate exact overlap | 4 个 Daood ABCG2 exact、0 个 A/B candidate；human-BBB context 只作连续/边界证据，不覆盖标签 |
| `bora_2025_human_bbb_abcg2_context_source_registry_v1.csv`、`bora_2025_human_bbb_abcg2_context_summary_v1.json`、`collect_bora_2025_human_bbb_abcg2_context.py` | Bora 2025 PMCID/DOI、BioC URL/SHA256、16 行计数和端点边界 | 非转染同源体系；“未运输/未给数值”不能填 ER=1，不并入当前 101 行直接 human/ortholog manifest |
| `li_2016_human_abcg2_benzylpenicillin_observations_v1.csv` | `observation_id`；4 条 human BCRP benzylpenicillin uptake/flux-boundary 观察 | Ko143 uptake 1.6±0.2 倍、197226 uptake 3.8±0.3 倍；无数值 ER/Papp；`model_binary_label` 全空 |
| `li_2016_human_abcg2_benzylpenicillin_identity_overlap_audit_v1.csv` | benzylpenicillin 标准化身份与 Daood/A-B candidate exact overlap | Daood/A-B candidate 均 0 exact；只作 accumulation/flux-boundary 外部证据 |
| `li_2016_human_abcg2_benzylpenicillin_source_registry_v1.csv`、`li_2016_human_abcg2_benzylpenicillin_summary_v1.json`、`collect_li_2016_human_abcg2_benzylpenicillin.py` | Li 2016 PMCID/DOI、BioC SHA256、4 行计数和无 ER 边界 | uptake rescue 不能替代二分类标签；“通量无变化”不补 ER=1 |
| `langthaler_2024_human_abcb1_reference_panel_observations_v1.csv` | `observation_id`；16 个参考化合物在 human MDR1-transfected MDCKII 中的 0.5 µM、1 h 双向 Papp/ER（14 条数值，2 条 `<LLOQ/NA`） | 3 likely、3 possible、9 unlikely、1 NA 仅保留原文 source 分类；直接人源连续 source-heldout 层，`model_binary_label` 全空 |
| `langthaler_2024_human_abcb1_reference_panel_identity_overlap_audit_v1.csv` | 16 个结构与 Daood ABCB1/A-B candidate 的 exact identity overlay | Daood exact 6、A/B candidate exact 5；重叠仅用于泄漏审计，不覆盖冻结标签 |
| `langthaler_2024_human_abcb1_reference_panel_source_registry_v1.csv`、`langthaler_2024_human_abcb1_reference_panel_summary_v1.json`、`collect_langthaler_2024_human_abcb1_reference_panel.py` | Langthaler 2024 PMCID/DOI、BioC SHA256、转染体系、回收率和 16 行计数 | human MDR1 direct ER/Papp 连续核查；`<LLOQ` 不转成阴性，`binary_label_transfer=none` |
| `wegler_2021_human_abcg2_crispr_panel_observations_v1.csv` | `observation_id`；15 个化合物在 hBCRP-transfected/cMDR1-KO MDCKII 与匹配 knockout 对照中的 figure-level ER 分类 | 11 条 human ABCG2 positive reference、4 条 direct negative-boundary（3 个 MDR1-selective、rifampicin）；不伪造逐实体数值 ER，`model_binary_label` 全空 |
| `wegler_2021_human_abcg2_crispr_panel_identity_overlap_audit_v1.csv` | 15 个结构与 Daood ABCG2/A-B candidate 的 exact identity overlay | Daood exact 11、A/B candidate exact 3；只作 source-heldout 泄漏审计 |
| `wegler_2021_human_abcg2_crispr_panel_source_registry_v1.csv`、`wegler_2021_human_abcg2_crispr_panel_summary_v1.json`、`collect_wegler_2021_human_abcg2_crispr_panel.py` | Wegler 2021 DOI/PDF SHA256、hBCRP/cMDR1-KO 体系、15 行计数和正/负边界 | direct human ABCG2 qualitative source layer；figure-level range 不能替代 per-row ER，也不能生成二分类标签 |
| `pmda_2018_olaparib_human_abcg2_observations_v1.csv` | `observation_id`；olaparib 在 parental 与 human BCRP-MDCKII 中 1/3/10 µM 双向 Papp/ER（6 条）及 5 µM Ko143 单向 B→A 救援（3 条） | human BCRP corrected ER=0.57/0.95/0.79；n=3；PMDA ATPase 解释为 unlikely substrate；连续弱/阴性外部层，`model_binary_label` 全空 |
| `pmda_2018_olaparib_human_abcg2_identity_overlap_audit_v1.csv` | olaparib 标准化身份与 Daood ABCG2/A-B candidate exact identity overlay | 当前 exact 命中均为 0；新结构适合 source-heldout 连续核查 |
| `pmda_2018_olaparib_human_abcg2_source_registry_v1.csv`、`pmda_2018_olaparib_human_abcg2_summary_v1.json`、`collect_pmda_olaparib_human_abcg2.py` | PMDA PDF URL/SHA256、端点条件、9 行计数和 Ko143/ATPase 边界 | regulatory human ABCG2 direct numeric weak-negative evidence；不把 corrected ER<1 自动转成二分类阴性 |
| `data_gap_closure_status_v1.csv`、`build_data_gap_closure_status.py` | `gap_id`；T10_GAP_001–005 的已补证据、剩余漏洞、下一步和冻结二分类池许可 | 仅为缺口状态审计；T10_GAP_001 已达 `closed_for_external_validation`，T10_GAP_005 为 `substantially_closed_external_only`；5 项均 `frozen_binary_pool=no` |
| `pmda_regulatory_human_transporter_contrast_panel_observations_v1.csv` | `observation_id`；acalabrutinib、alectinib、selumetinib 的 42 条 PMDA ABCB1/ABCG2 Papp、ER、parent、抑制剂和回收率观察 | ABCB1 20、ABCG2 22；31 条人源转染直接记录；14 条抑制剂 context；`model_binary_label` 全空 |
| `pmda_regulatory_human_transporter_contrast_panel_similarity_leakage_audit_v1.csv` | `entity_name × transporter_gene`；6 行 exact InChIKey、Murcko scaffold 和 Morgan 最近邻审计 | 三种新结构对两基因均无 exact/scaffold 重叠，最大 Morgan Tanimoto 0.213–0.273；只控制外部评估泄漏，不生成标签 |
| `pmda_regulatory_human_transporter_contrast_panel_source_registry_v1.csv`、`..._summary_v1.json`、`collect_pmda_regulatory_human_transporter_contrast_panel.py` | 三份 PMDA PDF URL、原始 SHA256/大小、逐来源行数和哨兵文本检查 | 监管连续对照层；源表 ER 不从四舍五入后的 Papp 均值重算；`binary_label_transfer=none` |
| `human_abcb1_external_source_manifest_v1.csv` | `source_id`；9 个来源、76 条 ABCB1 观察 | 70 条 human-direct；alectinib Caco-2 两行单列为非基因隔离 context；不把来源调用转二分类 |
| `human_abcb1_external_entity_overlay_v1.csv`、`human_abcb1_external_overlay_summary_v1.json`、`build_human_abcb1_external_overlay_manifest.py` | 40 个结构实体的来源、端点、source call、exact overlap 和可比性分层 | Daood exact 10、candidate exact 8；1 个跨来源结构复现；`model_label_nonblank=0` |
| `independent_human_transporter_validation_panel_observations_v1.csv` | `observation_id`；Zhang 2021 与 Barrera 2025 的 26 条 human ABCB1/ABCG2 直接观察 | ABCB1 20、ABCG2 6；含 20 条数值抑制剂 context；二分类标签全空 |
| `independent_human_transporter_validation_panel_similarity_leakage_audit_v1.csv` | `entity_name × transporter_gene`；10 行 exact、Murcko 和 Morgan 泄漏审计 | Daood/candidate exact 都为 0；地塞米松/倍他米松的甾体同骨架重叠显式标记 |
| `independent_human_transporter_validation_panel_source_registry_v1.csv`、`..._summary_v1.json`、`collect_independent_human_transporter_validation_panels.py` | Europe PMC 原始 XML、SHA256、逐来源行数和结构属性审计 | NU-9 类似物的文献 SwissADME TPSA 与 IUPAC 二酮互变异构标准化值分开保留；不以 TPSA 反推身份 |
| `adc_m0_targeted_search_audit_v1.csv` | `search_record_id`；5 个 ADC 游离载荷—释放物定向检索单元 | 记录双侧已有证据、当前最高 M 等级和 M0 排除原因；M0 仍为 0，不计算跨来源倍数 |
| `wegler_2021_human_abcb1_crispr_panel_observations_v1.csv` | `observation_id`；14 个 human ABCB1/cMDR1-knockout 观察 | 12 个阳性参照、2 个弱阴性边界；fluvastatin ER 1.8、dantrolene ER 0.6；余下仅保留组范围 |
| `wegler_2021_human_abcb1_crispr_panel_similarity_leakage_audit_v1.csv` | `entity_name`；14 行 exact、既有外部重复、Murcko 和 Morgan 审计 | 5 个 Daood exact、4 个 candidate exact、1 个跨来源 exact；两个弱阴性结构 exact 均为 0 |
| `wegler_2021_human_abcb1_crispr_panel_source_registry_v1.csv`、`..._summary_v1.json`、`collect_wegler_2021_human_abcb1_crispr_panel.py` | 原文 PDF/纯文本、出版社 `mmc1.xlsx`、SHA256 和行数质控 | 补充表 S1 只有 LC–MS/MS 参数，不伪造浓度或逐化合物 ER |
| `wegler_2021_abcb1_external_predictions_v1.csv` | `entity_name × holdout_type × model`；Wegler 14 个外部实体在三种冻结基线上的概率、预测、适用域和泄漏标记 | 仅为 source-heldout 外部压力测试；来源调用不参加训练、不回填标签，也不用于重调阈值 |
| `wegler_2021_abcb1_external_diagnostics_v1.csv`、`..._external_stress_test_summary_v1.json`、`run_wegler_abcb1_source_scaffold_stress_test.py` | all n=14、exact-new n=9、strict source×scaffold-new n=5 的描述性 BAcc/MCC/AUC/AP 与混淆矩阵 | strict n=5 中 candidate-supervised 的两条弱阴性均误报；PU 排序虽 5/5 正确，但样本过小且不是校准分类器，不能重开 predictor claim gate |
| `independent_human_transporter_negative_contrast_panel_observations_v1.csv` | `observation_id`；CNO 与 secnidazole 的 12 条 human ABCB1/ABCG2 直接观察 | 6 条 ABCB1、6 条 ABCG2；CNO ABCB1 ER 20→valspodar≈1.2，ABCG2 ER≈3→Ko143≈3.2；secnidazole 三浓度 parent-corrected NFR 均<2；`model_binary_label` 全空 |
| `independent_human_transporter_negative_contrast_panel_similarity_leakage_audit_v1.csv` | `entity_name × transporter_gene`；4 个来源调用的 exact、Murcko 和 Morgan 审计 | 4 个组合均 Daood/candidate/既有外部 exact=0；CNO strict scaffold-new，secnidazole 各有 1 个低相似度相同 Murcko；泄漏审计不生成标签 |
| `independent_human_transporter_negative_contrast_panel_source_registry_v1.csv`、`..._summary_v1.json`、`collect_independent_human_transporter_negative_contrast_panel.py` | Raper CNO BioC、FDA secnidazole PDF、PubChem 身份文件的 URL、SHA256、大小和哨兵文本核验；12 行/2 来源/2 实体摘要 | direct human external source-call/continuous layer；3 阴性+1 阳性仅作外部诊断，`binary_label_transfer=none` |
| `negative_contrast_frozen_predictions_v1.csv` | `entity_name × transporter_gene × holdout_type × model`；4 个 exact-new 调用在三套冻结配方上的概率与预测 | 所有来源调用完全 held out；阈值不回调；只暴露失败模式，不回填二分类标签 |
| `negative_contrast_frozen_diagnostics_v1.csv`、`negative_contrast_frozen_probe_summary_v1.json`、`run_negative_contrast_source_scaffold_probe.py` | exact-new n=4、strict source×scaffold-new n=2 的描述性诊断和 1,137 行冻结接口写入审计 | 三种配方分别正确 1/4、2/4、2/4；每基因 n=2，禁止作正式性能声明；`model_binary_label_writes=0` |
| `adc_m0_targeted_search_audit_v2.csv` | `search_record_id`；5 个 free payload/released species 检索单元的二次复核 | M0=0；cys-mcMMAF GSK/PMDA 措辞冲突保留；DM4/S-methyl-DM4 联合定性调用不升级为同 run 配对 |
| `adc_m0_second_pass_source_registry_v1.csv`、`adc_m0_targeted_search_summary_v2.json`、`audit_adc_m0_second_pass.py` | GSK protocol 与复用 PMDA/FDA 原始文件的 SHA256、用途、冲突和 5 单元结果 | `strict_binary_labels_added=0`、`frozen_binary_pool_yes=0`；负向检索结果不伪装成 M0 |
| `payload_relevant_abcg2_negative_boundary_observations_v1.csv` | `observation_id`；Litman 2000 与 Miwa 2003 的 7 条 ABCG2 上下文观察、6 个结构实体 | vinblastine/paclitaxel 为 selected human overexpression 亚系未外排边界；doxorubicin 为 human wild-type BCRP 阴性而 R482 变体阳性；所有 `model_binary_label` 为空 |
| `payload_relevant_abcg2_negative_boundary_similarity_leakage_audit_v1.csv` | `entity_name`；vinblastine、paclitaxel、doxorubicin 三条阴性探针的 exact/Murcko/Morgan 与项目载荷身份双重审计 | 三者均 Daood/candidate exact-new；vinblastine 有同 Murcko 重叠；doxorubicin 通过 `source_inchikey` 解析项目标准化互变异构差异；审计不生成标签 |
| `payload_relevant_abcg2_negative_boundary_source_registry_v1.csv`、`..._summary_v1.json`、`collect_payload_relevant_abcg2_negative_boundaries.py` | 两份 PubMed XML、6 份 PubChem JSON 的 URL、SHA256、大小、摘要哨兵和身份硬校验 | 2 来源、7 观察、6 实体、2 个项目登记载荷精确阴性；定性/体系受限，`binary_label_transfer=none` |
| `assay_attribution_counterexamples_v1.csv` | `counterexample_id`；CNO 抑制剂不救援、olaparib 亲本校正、doxorubicin R482 序列上下文三个反例 | 用于阻止“ER 阈值或化合物唯一标签”的错误迁移；不构成训练行 |
| `payload_relevant_abcg2_negative_frozen_predictions_v1.csv` | 历史快照：`entity_name × holdout_type × model`；三个当时载荷邻近阴性在三种冻结配方上的概率、0.5 描述性误报与适用域 | 当时 candidate-supervised 为 3/3；第二十二批已撤回 doxorubicin `target=0` 的正误解释，当前评价必须使用 v2 |
| `payload_relevant_abcg2_negative_frozen_diagnostics_v1.csv`、`..._probe_summary_v1.json`、`run_payload_relevant_abcg2_negative_probe.py` | 历史 all n=3、registered payload n=2、strict source×scaffold-new n=2 单类探针 | 只为审计链保留；禁止报告 BAcc/MCC/AUC/AP 或据此选模；训练基础 1,121、完整冻结接口 1,137，`model_binary_label_writes=0` |
| `payload_abcg2_numeric_context_observations_v1.csv` | `observation_id`；Litman、Robey、Otter、Fung 的 23 条载荷邻近 ABCG2 观察，其中 22 条定量 | vinblastine/paclitaxel 保留上下文阴性边界；doxorubicin 为 WT/变体/直接转运冲突；`model_binary_label` 全空 |
| `payload_abcg2_system_conflict_matrix_v1.csv` | `entity_name`；3 个实体的阴性/低效背景、阳性/冲突背景、机制轴与二分类处理 | doxorubicin 为 `reject_compound_only_label`；vinblastine/paclitaxel 只用于外部负边界，不转训练真值 |
| `payload_abcg2_numeric_context_similarity_leakage_audit_v1.csv` | 结构 exact/Murcko/Morgan、载荷身份、旧来源与新增定量上下文的联合审计 | 明确 clean-negative eligibility 与 universal binary resolution；不靠多数投票消除体系冲突 |
| `payload_relevant_abcg2_negative_frozen_predictions_v2_context_audited.csv` | vinblastine/paclitaxel × 三个冻结配方 × 三个单类面板；共 12 条预测 | candidate-supervised 对 all panel 误报 2/2；两者均域外且应拒绝；不是性能验证 |
| `payload_abcg2_context_conflict_frozen_predictions_v1.csv` | doxorubicin 的 3 条历史冻结概率，无 target/correctness | `not_scored_context_conflict`；禁止重新解释成阳性或阴性性能 |
| `payload_relevant_abcg2_negative_frozen_probe_summary_v2_context_audited.json`、`run_payload_abcg2_context_audited_probe.py` | all n=2、registered payload n=1、strict source×scaffold-new n=1；训练基础 1,121、完整冻结接口 1,137 | v1 文件保留历史快照，但 doxorubicin `target=0` 的正误解释已撤回；`model_binary_label_writes=0` |
| `payload_abcg2_direct_knockout_observations_v1.csv` | `observation_id`；paclitaxel/vinblastine 的 38 条人细胞 direct-KO/transfectant/inhibitor 观察 | 37 条精确数值或下界；Volpe paclitaxel 27 条 Papp/ER、Seamon 6 条 IC50/RR、Kaliszczak 3 条 ER/抑制剂归因、Plise 2 条 KO ER；`model_binary_label` 全空 |
| `payload_abcg2_direct_knockout_similarity_leakage_audit_v1.csv` | 两个标准化实体的 exact/Murcko/Morgan、载荷身份、旧/新增来源和 evidence-upgrade 联合审计 | 结构与相似性结果沿用并硬核验：paclitaxel strict source×scaffold-new；vinblastine 与项目载荷 exact、与 Daood 有同 Murcko；重复 assay 条件不增加独立 compound n |
| `payload_abcg2_negative_evidence_upgrade_matrix_v1.csv` | `entity_name`；旧证据、新直接证据、解释、剩余限制与推荐模型角色 | 两者均升级为 direct multi-KO/transfectant/inhibitor negative-attribution boundary；`binary_training_eligible=no`、`strict_binary_label_added=0` |
| `payload_abcg2_direct_knockout_source_registry_v1.csv`、`..._summary_v1.json`、`collect_payload_abcg2_direct_knockout_boundaries.py` | 4 个来源的原始 HTML/XML/JPEG/PDF、访问记录、SHA256、逐来源行数和哨兵文本核验 | Seamon 出版社全文 CLI 403 时保留 web-verified extraction record；Plise PDF 已本地取得并通过全文哨兵；所有来源只进外部证据层 |
| `payload_relevant_abcg2_negative_frozen_predictions_v3_evidence_upgraded.csv`、`..._diagnostics_v3_evidence_upgraded.csv` | 在 v2 12 条预测/9 条诊断上附加新直接证据、限制和允许用途 | 概率、阈值和预测语义哈希与 v2 完全相同；candidate-supervised 误报 2/2 且拒绝 2/2；n=2 单类域外面板不能用于选模或性能声明 |
| `payload_relevant_abcg2_negative_frozen_probe_summary_v3_evidence_upgraded.json`、`build_payload_abcg2_evidence_upgraded_probe.py` | 证据升级、冻结预测不变与 1,137 行二分类接口不变性的机器摘要 | v2/v3 预测语义哈希相同，接口哈希 `0c098d...`、changed rows=0；只增强外部 target 可信度，不改变模型 |
| `payload_abcg2_direct_knockout_validation_v1.json`、`validate_payload_abcg2_direct_knockout_batch.py` | 38 行逐值、结构、来源文件哈希、manifest、标签空值、预测与接口不变量的自动验证 | `pass=true`；证明数据管线和冻结边界一致，不等于生物学二分类验证通过 |
| `payload_near_human_abcg2_reference_observations_v1.csv` | `observation_id`；Roth、Takagi、Yoshikawa 共 23 条人 ABCG2 载荷/近邻观察 | 21 条精确 RR/成对数值、2 条 SN-38/SN-38 glucuronide 囊泡定性直接转运；全部 `model_binary_label` 为空 |
| `payload_near_human_abcg2_reference_source_registry_v1.csv`、`payload_near_abcg2_balanced_reference_summary_v1.json`、`collect_payload_near_abcg2_balanced_reference.py` | 3 个来源的 PubMed/Crossref/原文抽取/PubChem/既有 Roth 原文哈希与 20 实体裁决摘要 | 出版社 CLI 403 时保留访问记录与 web-verified Table 2 抽取；允许外部压力测试，不转训练标签 |
| `takagi_cpt_keto_structure_resolution_v1.csv` | S38809/S39625 的机器结构解析审计 | S38809 未取得可核验结构；S39625 PubChem 名称查询返回立体歧义记录；两者均不进入模型探针 |
| `payload_near_abcg2_reference_similarity_leakage_audit_v1.csv` | 20 个裁决候选的 full InChIKey、connectivity、Murcko、Morgan 与冻结角色重叠 | 平衡面板对 candidate-supervised 连接性新颖 11/11；对 auxiliary/PU 有 5 个 connectivity overlap；禁止只凭完整 InChIKey 宣称无泄漏 |
| `payload_near_abcg2_reference_adjudication_matrix_v1.csv` | 20 个候选的纳入/排除、作者调用、来源、身份和冲突裁决 | 保留 11 个外部调用，排除 9 个冲突/弱响应/身份未闭环对象；`binary_training_eligible=no` |
| `payload_near_abcg2_balanced_external_reference_calls_v1.csv` | 11 个结构已解析外部参照：6 阳性、5 阴性/阴性边界 | 只为 source-call 外部压力测试；7 个同一 Roth assay，paclitaxel/vinblastine 仍为上下文边界；不等于通用真值 |
| `payload_near_abcg2_balanced_frozen_predictions_v1.csv` | 11 实体 × 3 个冻结模型，共 33 条概率 | 固定阈值 0.5；逐模型标出 connectivity overlap、适用域和正误；重叠行只作预测审计，不作性能诊断 |
| `payload_near_abcg2_balanced_frozen_diagnostics_v1.csv`、`..._probe_summary_v1.json`、`run_payload_near_abcg2_balanced_frozen_probe.py` | candidate connectivity-new n=11 平衡诊断，以及 auxiliary/PU 去重后的残余小面板 | candidate 在 0.5 下 sensitivity=1、specificity=0、BAcc=0.5、MCC=0；ROC-AUC=0.9 仅为小样本排序信号，禁止据此调阈值或宣称验证通过 |
| `payload_near_abcg2_balanced_validation_v1.json`、`validate_payload_near_abcg2_balanced_batch.py` | 逐值、来源哈希、结构/连接性泄漏、清单计数、冻结接口与探针结果的自动验证 | `pass=true`；训练新增行=0、模型标签写入=0、阈值回调=0，接口仍 1,137 行 |
| `adc_release_species_abcg2_structure_extension_v1.csv` | `entity_name`；5 个 ADC 释放物/代谢物的标准化结构扩展 | 补齐 N-Ac-γ-calicheamicin DMH 和 SN-38G；保留 cys-mcMMAF 名称/互变异构体与 calicheamicin 质量上下文警示，不覆盖原始注册表 |
| `adc_release_species_human_abcg2_calibration_observations_v1.csv` | `observation_id`；5 个来源、11 条真实释放物或活性代谢物 ABCG2 观察 | 7 条定量 ER、4 条定性来源调用；`model_binary_label` 全空，混合 Caco-2 和无公开数值表的监管/方案结论分层保存 |
| `adc_release_species_abcg2_calibration_source_registry_v1.csv` | `source_id`；5 份 PMDA/TGA/FDA/NCI 来源及人工抽取记录的本地路径、SHA256 与大小 | 以 source ID 而非条件行计算独立来源；同一人工抽取 JSON 的复用不增加独立证据数 |
| `adc_release_species_abcg2_calibration_similarity_leakage_audit_v1.csv` | 6 个结构已解析候选的 exact/connectivity/Murcko/Morgan 与冻结模型、第二十四批重叠 | 严格 4 实体候选对 candidate 和第二十四批均为 entity/source/connectivity-new；SN-38G 因来源重叠排除，DM4 因跨来源冲突排除 |
| `adc_release_species_abcg2_calibration_adjudication_matrix_v1.csv` | 8 个候选的身份、证据、冲突、来源独立性与纳排裁决 | FD1/FD5 仅保留未解析排除记录；`binary_training_eligible=no`、`threshold_tuning_eligible=no` |
| `adc_release_species_abcg2_independent_calibration_candidates_v1.csv` | 4 个 source/entity-heldout 真实 ADC 相关候选：1 阳性、3 阴性 | 只是锁定的校准候选池；证据异质且类别不足，不能直接拟合阈值或校准器 |
| `adc_release_species_abcg2_frozen_predictions_v1.csv` | 4 实体 × 3 个冻结模型，共 12 条预测 | candidate-supervised 固定 0.5 下 TP=1、TN=0、FP=3、FN=0；3 个阴性触发相似性拒判；无阈值/模型回调 |
| `adc_release_species_abcg2_fixed_threshold_diagnostics_v1.csv`、`..._probe_summary_v1.json` | 三个冻结配方的小样本描述性诊断与冻结接口哈希 | n=4、仅 1 阳性；BAcc/AUC/AP 不构成验证或校准性能，禁止据此选阈值/选模型 |
| `adc_release_species_abcg2_calibration_feasibility_v1.csv` | 4 个预声明项目门槛的 current/required/shortfall/pass | 候选类计数还缺 +4 阳性/+2 阴性；direct numeric gene-isolated 还缺 +3/+2；这是项目收数门槛，不是通用统计规则 |
| `adc_release_species_abcg2_calibration_validation_v1.json`、`validate_adc_release_species_abcg2_calibration_batch.py` | 结构、逐值、来源哈希、泄漏、诊断、manifest 与冻结接口的自动验证 | `pass=true`；模型标签写入、阈值拟合、校准拟合、训练接口新增均为 0 |

### 第十五批 assay-aware 端点与机制基础层

| 表 | 主键/用途 | 关键边界 |
|---|---|---|
| `assay_endpoint_ontology_v1.csv` | `endpoint_family`；7 类观察端点的定义、兼容性和二分类接口 | 不允许把 Papp/ER、功能耐受、PAMPA、滞留或完整 ADC 上下文池化成一个标签 |
| `assay_observations_unified_v1.csv` | `observation_id`；139 条统一观察 | `entity_layer`、`transporter_gene`、`gene_specific_context`、`truth_role` 和 `model_use` 明确观察能否进入训练；缺失数值不补 0 |
| `theory_entity_registry_v1.csv` | `chemical_entity_key`；游离载荷、释放物、完整 ADC 与无标签理论参照的实体注册 | `label_transfer_policy` 固定为不跨 exposure layer 转移；Topotecan 不带 ADC 标签 |
| `benchmark_split_manifest_v1.csv` | `manifest_id`；普通开发/PU/独立来源/ADC上下文/人工审计角色 | split 在接触 ADC 结局前固定；不按 ADC 结果调参 |
| `benchmark_preregistration_v1.json` | 预注册协议和输入 SHA256 | 普通 gate 当前关闭；禁止跨端点池化、阈值回调和标签迁移 |
| `physiological_microstates_v1.csv` | `microstate_id`；197 个规则枚举候选质子化状态 | 候选状态不是 pKa 概率、主导状态或生理分布真值 |
| `conformer_access_descriptors_v1.csv` | `microstate_id`；低成本 RDKit 真空构象/2D 描述符 | `geometry_status` 需与数值一起读取；安全跳过/排队不等于负标签，真空描述符不等于膜 PMF |
| `efflux_phase_diagram_v1.csv` | 无量纲 `access/recognition/retention` 情景 | 合成机制假设，用于敏感性和实验设计，不是生物学验证 |
| `efflux_identifiability_simulation_v1.csv`、`efflux_identifiability_summary_v1.csv` | 合成 endpoint 组合的参数可辨识性 | 只说明单一 ER 无法分离多个机制参数；不能估计真实 transporter 参数 |
| `computational_anchor_simulation_plan_v1.csv` | `anchor_system_id`；4 组 matched 体系的后续计算/实验计划 | 允许的 claim 为“可检验机制假设”；禁止称为残基因果效应已证明 |
| `external_bro5_dataset_manifest_v1.csv` | CREMP-CycPeptMPDB/NPMMP-DB 外部邻域数据注册 | 只作 bRo5 被动通透/构象校准，不作 ABCB1/ABCG2 底物真值 |
| `external_bro5_download_registry_v1.csv`、`external_bro5_download_summary_v1.json` | 实际下载文件的状态、大小、SHA256 和归档成员数 | `download_status` 不是实验质量等级；任何记录均禁止进入二分类标签 |
| `external_bro5_cremp_annotations_v1.csv` | `source_record_id`；CREMP 3,258 条结构/被动通透/构象注释 | `endpoint_unit_status=source_summary_unit_not_rewritten`；只能进入固定外部 bRo5 校准，不进入 transporter target |
| `external_bro5_npmmp_raw_v1.csv` | 5,638 条 NPMMP 原始观察；PAMPA/Caco-2/MDCK/RRCK/Others 与 Papp/ER 分类 | `transporter_gene=not_applicable`、`gene_specific_context=no`；不得把 ER 解释为 ABCB1/ABCG2 真值 |
| `external_bro5_npmmp_archive_audit_v1.json` | NPMMP 压缩包完整性、成员数和表格成员名审计 | `partial_or_incomplete_archive` 时不读取、不抽取、不用于模型 |
| `external_bro5_permeability_validation_v1.csv` | `external_record_id`；8,896 条 CREMP/NPMMP 外部通透性记录与固定 folds | `model_use` 仅为外部 permeability/conformer validation；不进入 ABCB1/ABCG2 binary target |
| `external_bro5_source_endpoint_counts_v1.csv`、`external_bro5_validation_summary_v1.json` | assay/endpoint/source group 计数、结构解析和删失审计 | 不把不同 assay family 合并成一个观测 |
| `external_bro5_oof_predictions_v1.csv` | 冻结 scaffold/source/source×scaffold 的逐行 OOF 预测 | 仅含 CREMP/NPMMP 被动通透目标；禁止解释为 ABCB1/ABCG2 预测 |
| `external_bro5_model_diagnostics_v1.csv`、`external_bro5_gate_summary_v1.json` | 逐折/聚合指标与外部低成本 Go/No-Go | 当前 CREMP 构象集合增量未过门；只否定当前特征模块，不否定膜机制 |
| `theory_runtime_capability_audit_v1.json` | 本地计算能力与缺失工具审计 | 无 gmx/xtb/crest/obabel/MD/Bayesian sampler 时，不声称高成本模拟已完成 |
| `theory_foundation_binary_invariance_audit_v1.json` | assay-aware 层接入后的冻结二分类接口不变性审计 | 目标标签、角色和转运体逐行比较；0 个改变才允许继续使用辅助层 |
| `theory_foundation_validation_v1.json` | 15 项数据/机制不变量和关键文件 SHA256 | 通过只说明实现、门控状态与合成理论自洽；不代表生物学或转运体预测验证 |

本层是机制/审计辅助层，不改变现有二分类训练表。任何后续新增直接实验必须先映射到本体中的兼容端点，再决定是否能进入 frozen binary core；跨 exposure layer、跨 transporter 或跨 assay family 的自动标签转移均禁止。

### ADC预测失效机制案例表

| 表 | 主键/用途 | 关键边界 |
|---|---|---|
| `adc_prediction_failure_case_matrix_v1.csv` | `case_id`；13个实体×2个转运体的26行案例 | `formal_failure_call` 当前全部待定；旧模型未通过gate，只能作历史候选定位 |
| `adc_residual_pair_failure_hypotheses_v1.csv` | `pair_id`；9组游离载荷—释放物比较 | `M0_same_assay_pair_available` 当前全为 `no`，不得计算残基因果倍数 |
| `adc_failure_mechanism_gap_queue_v1.csv` | `gap_id`；6类机制缺口 | 接收标准在冻结模型和新增数据前固定 |
| `adc_frozen_predictor_import_template_v1.csv` | `case_id`；预测器任务接口 | 只导入不可变模型版本/hash、盲测概率和域判定；禁止ADC调参 |
| `adc_failure_mechanism_summary_v1.json` | 案例库机器可读摘要 | 明确只组织失效假设，不证明普遍ADC特异失效 |

`adc_prediction_failure_case_matrix_v1.csv` 关键字段：

| 字段 | 含义 |
|---|---|
| `entity_layer` | 游离载荷/释放战斗部或带残基释放物；完整ADC不进入案例粒度 |
| `frozen_*` | 等待另一个任务交付的正式冻结预测及适用域判定；当前留空 |
| `legacy_*` | 未通过claim gate的历史基线输出，只能用于候选定位和版本对照 |
| `actual_outcome_evaluability` | 区分单转运体直接数值、抑制剂归因、多转运体、数字化、功能旁证和无结局 |
| `failure_test_readiness` | 冻结预测到位后能否做正式方向失效检验 |
| `pampa_mean_digitized_1e6_cm_s` | 被动膜通透数字化值；映射到两个基因案例但不是基因特异证据 |
| `pair_*` / `matched_*` | 游离—释放物配对、匹配等级和方向；无M0时只允许机制假设 |
| `domain_shift_hypothesis` | 适用域偏移假设；历史相似性拒绝不能替代冻结模型域规则 |
| `membrane_access_hypothesis` | PAMPA、胞内滞留或直接转运存在时的膜可达性数据状态 |
| `endpoint_mismatch_hypothesis` | 直接转运、功能RR、监管调用等是否跨终点不可比 |
| `source_context_conflict_state` | 显式来源上下文冲突，例如SG3199/ABCG2 |
| `identity_uncertainty_hypothesis` | 活动身份hold、注册状态不一致复核或当前无hold |
| `claim_boundary` | 正式ADC特异失效需要普通域有效、预先冻结且未用ADC调参的模型 |

正式失效判据：冻结普通预测器方向与兼容的直接转运结局不一致。适用域拒绝本身回答“模型无资格判断”，不是方向错误；PAMPA/胞内滞留回答膜可达性层，不是转运体底物标签；完整ADC仍只作上下文桥。

### 第四批结构标准化公共字段

这些字段出现在 `compound_registry_standardized_v1.csv`、`adc_release_species_registry_standardized_v1.csv` 和 `daood_training_records_standardized_v1.csv`；Daood 共识表保留建模需要的子集。

| 字段 | 含义 |
|---|---|
| source_smiles | 不可变输入 SMILES；不作为标准化结果覆盖 |
| standardization_status | standardized、unresolved_blank_structure、parse_failed 或 standardization_failed |
| standardized_isomeric_smiles | Cleanup、FragmentParent、Uncharger 和 canonical tautomer 后的 canonical isomeric SMILES |
| standardized_connectivity_smiles | 去掉立体标注后的连接性模型键 |
| standardized_inchikey / inchikey_connectivity_block | Standard InChIKey 及其第一段，用于精确连接性审计 |
| murcko_scaffold_smiles / generic_murcko_scaffold_smiles | exact 与 generic Bemis–Murcko scaffold |
| fragment_removed_flag | 是否在 FragmentParent 步骤移除了次要片段 |
| source_formal_charge / standardized_formal_charge | 标准化前后总形式电荷 |
| tautomer_changed_flag | canonical tautomer 是否改变了连接性 SMILES 表达 |
| tautomer_enumeration_status | Completed、MaxTautomersReached、MaxTransformsReached、Canceled 或 not_run |
| tautomer_enumerated_count | 实际枚举的 tautomer 数量 |
| tautomer_limit_flag | 非 Completed 时为 yes；表示 best-effort QC，不等于结构解析失败 |
| defined/unassigned_stereocenter_count | 已定义和未定义手性中心计数 |
| std_* | 标准化结构的 MW、XlogP、TPSA、HBD、HBA、可旋转键、形式电荷和重原子数 |

### `daood_training_consensus_v1.csv`

每行是一个“标准化连接性结构—transporter”对象。相同对象标签一致时合并来源；标签冲突时 `consensus_binary_label` 留空并进入冲突表。本次没有发现标签冲突。

### `daood_training_leakage_decisions_v1.csv`

| 字段 | 含义 |
|---|---|
| exact_external_reference_ids | 与 ADC 外部对象具有相同连接性结构的引用 ID |
| near_duplicate_external_reference_ids | Morgan Tanimoto ≥0.80 的引用 ID |
| complex_scaffold_external_reference_ids | 同 exact Murcko scaffold、scaffold ≥10 重原子且 Tanimoto ≥0.50 的引用 ID |
| max_external_morgan_tanimoto | 对全部可解析 ADC 参考结构的最高 Tanimoto |
| leakage_model_action | 三类排除决定或 `include_leakage_clean_training` |

排除优先级为 exact > near duplicate > complex scaffold family；简单/common scaffold 只标记，不自动排除。

### 基线与适用域表

| 表 | 含义 |
|---|---|
| baseline_scaffold_cv_metrics_v1.csv | 5 折 scaffold-grouped 每折及 mean/sd 指标 |
| baseline_scaffold_oof_predictions_v1.csv | 每条清洁记录的 OOF 概率和折号 |
| baseline_similarity_support_thresholds_v1.csv | 各 transporter 训练 self-nearest Tanimoto 的 q10/q25 与拒绝规则 |
| baseline_roth_external_predictions/metrics_v1.csv | Roth 连续 RR 的预测及 Spearman 方向性审计，无二分类标签 |
| baseline_release_evidence_predictions/metrics_v1.csv | 按 evidence stratum 聚合的释放物描述性外部审计 |
| baseline_pair_probability_shifts_v1.csv | 游离载荷与释放物的模型结构敏感性；不作因果解释 |
| baseline_model_gate_summary_v1.csv | 内部 gate、外部信号与最终 predictor claim 状态 |
| baseline_model_manifest_v1.csv | 6 个 joblib 文件的训练数量、路径、SHA256 和 claim 状态 |

## 禁止的处理

- 不用空值填 0；
- 不覆盖原始 SMILES；
- 不把 inhibitor 标签当成 substrate 标签；
- 不把同一 compound 在不同 transporter 上的记录合并；
- 不把模型预测结果写回真值列；
- 不在完成结构去重前把 Roth 记录并入训练集。
- 不把游离 MMAF 的转运标签直接复制给 cys-mcMMAF，也不把 DM1 标签复制给 Lys-MCC-DM1。
- 不把“细胞内保留更多”或“ADC 细胞毒更强”直接当成转运体底物标签。
- 不跨细胞系或端点计算游离载荷—释放物倍数；描述符变化也不能单独定义 P-gp/ABCG2 底物行为。
- 不把相似性支持域外的高概率解释为高置信；`abstain_recommendation=yes` 优先于概率值。
- 不以 ROC-AUC 单指标宣布成功；当前门槛为 scaffold CV mean MCC ≥0.50 且 balanced accuracy ≥0.75，并仍要求独立 matched released-species 验证。
- 不把 M1H 当作已确认同实体重复；不把 M2-M4 的数值相除得到“残基效应倍数”。
- 不把同一图像的两种算法校准称为双人独立数字化。
- 不把 P-gp/BCRP 共存 Caco-2 体系的同一 ER 同时复制成两条基因特异训练标签。
- 不在实验试剂身份未解析时以数据库 exact-name 返回结构代替材料结构；双结构敏感性输出必须保留 abstain。

## 人源 ABCG2 普通域桥接与前瞻确认表

### `human_abcg2_direct_bridge_balanced_panel_v1.csv`

| 字段 | 含义 |
|---|---|
| `bridge_candidate_call` | 原始来源支持的普通域外部调用，1=substrate、0=non-substrate；不是训练标签 |
| `source_id` | 每个面板实体唯一的证据来源标识 |
| `human_abcg2_transfectant_rows` | 该实体在人源 ABCG2 转染体系中的观察数 |
| `candidate_connectivity_match_count` | 与 candidate-supervised 模型输入的 connectivity 重叠数 |
| `daood_connectivity_match_count` | 与 Daood 辅助训练层的 connectivity 重叠数 |
| `candidate_supervised_leakage_status` | candidate 模型可否进入 connectivity-new 诊断 |
| `auxiliary_pu_leakage_status` | 辅助模型可否进入 connectivity-new 诊断 |
| `threshold_tuning_eligible` | 固定为 `no_until_separate_preregistration`；本表不允许调阈值 |
| `model_binary_label` | 必须为空；外部来源调用不得自动迁移到冻结标签 |

### `human_abcg2_prospective_evidence_observations_v1.csv`

| 字段 | 含义 |
|---|---|
| `prospective_source_call` | 锁定前由来源证据确定的外部调用；不是模型标签 |
| `transporter_system` | human BCRP 转染细胞、亲本或对照体系 |
| `substrate_concentration_*` | 被测物浓度值/范围及单位 |
| `inhibitor*` | Ko143、Ko134 等抑制剂及浓度；空值表示该观察不含抑制剂 |
| `papp_ab_*` / `papp_ba_*` | 双向 Papp 均值、误差和单位；不把缺失项反推为 0 |
| `efflux_ratio_*` | 来源报告的 ER 均值、误差或范围 |
| `relative_or_net_efflux_ratio*` | 亲本校正或净 ER；与原始 ER 分列保存 |
| `recovery_*` | 来源报告的回收率；用于保留 QC 限制 |
| `source_assigned_call` | 来源中的 substrate/non-substrate/对照语义 |
| `qc_flag` | 低回收、高亲本外排等解释限制 |
| `strict_binary_label_added` | 必须为 0 |

### `human_abcg2_prospective_structure_registry_v1.csv`

| 字段 | 含义 |
|---|---|
| `pubchem_cid` | 锁定的 PubChem 化合物 ID，防止同名结构漂移 |
| `source_*` | PubChem 原始结构、名称、式量和本地 JSON 哈希 |
| `standardized_*` | 项目清洗、去片段、去电荷和 canonical tautomer 后的派生结构 |
| `inchikey_connectivity_block` | 忽略立体/质子层的 connectivity 身份键 |
| `murcko_scaffold_smiles` | scaffold 分析所用的 Murcko 骨架 |
| `source_to_project_standardized_inchikey_match` | 原始 PubChem 身份与项目标准化身份是否一致；差异必须保留说明 |

### `human_abcg2_prospective_similarity_leakage_audit_v1.csv`

| 字段 | 含义 |
|---|---|
| `*_exact_match_count` | 与指定参考层的完整标准化 InChIKey 重叠数 |
| `*_connectivity_match_count` | 与指定参考层的 connectivity 重叠数 |
| `*_same_murcko_scaffold_count` | 同 Murcko scaffold 的参考实体数；不等同 exact 泄漏 |
| `*_max_morgan_tanimoto` | 到指定参考层的最大 Morgan 相似度 |
| `*_nearest_name` | 最大 Morgan 相似度对应实体 |
| `prior_external_panel_overlap_status` | 对第二十四、第二十五和第二十六批外部面板的隔离状态 |
| `prediction_status` | 锁定时为 `locked_not_run`，证明选择先于预测 |

### `human_abcg2_prospective_confirmation_queue_v1.csv`

该表把来源、结构、标准化、泄漏和选择字段连接为一次性预测队列。`selection_status=locked_before_prediction`、`threshold_tuning_eligible=no`、`sample_weight=0.0` 和空 `model_binary_label` 共同限定：该队列只作外部诊断，不参与训练、阈值拟合或样本加权。

### 冻结预测与诊断

| 表/字段 | 含义 |
|---|---|
| `human_abcg2_prospective_frozen_predictions_v1.csv` | 每实体 × 三模型的冻结概率、0.5 分类、正确性和适用域建议 |
| `similarity_abstain_recommendation` | `no_similarity_abstention`、`yes_caution` 或 `yes_abstain`；不是第四种标签 |
| `external_reference_not_model_label` | 必须为 `yes` |
| `threshold_tuned_on_panel` / `calibration_model_fit_on_panel` | 必须为 `no` |
| `human_abcg2_prospective_fixed_threshold_diagnostics_v1.csv` | 每模型的 TP/TN/FP/FN、敏感性、特异性、BAcc、MCC、AUC、AP、Brier |
| `human_abcg2_prospective_lock_manifest_v1.json` | 选择规则、主要/次要终点和输入文件哈希；任何输入变化都应使验证器失败 |

普通域桥接和前瞻表不能与 ADC 实际释放物表合并计算一个“ADC 准确率”。普通域用于诊断基础 ABCG2 二分类器；ADC 特异主张必须在真实释放物、对应结构和直接实验层独立满足门槛。
