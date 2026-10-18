#!/usr/bin/env python3
"""Frozen finite ADC case checks and one-source TOP1 order comparisons (stdlib)."""
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def read(path):
    return list(csv.DictReader(path.open(encoding='utf-8-sig')))


def write(name, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/name).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def ranks(xs):
    out = [0.0] * len(xs)
    ordered = sorted(range(len(xs)), key=lambda i: xs[i])
    start = 0
    while start < len(xs):
        end = start+1
        while end < len(xs) and xs[ordered[end]] == xs[ordered[start]]:
            end += 1
        rank = (start+1+end)/2
        for i in ordered[start:end]:
            out[i] = rank
        start = end
    return out


def spearman(xs, ys):
    x, y = ranks(xs), ranks(ys)
    mx, my = sum(x)/len(x), sum(y)/len(y)
    num = sum((a-mx)*(b-my) for a,b in zip(x,y))
    den = math.sqrt(sum((a-mx)**2 for a in x)*sum((b-my)**2 for b in y))
    return num/den if den else None


inputs = [
    HERE/'adc_source_context_positive_anchors.csv',
    ROOT/'数据/processed/obi992_caco2_payload_transport_numeric_v1.csv',
    ROOT/'outputs/20261009_adc_pair_experiments/pairs/results_curated_v2/entities_scored_or_hold.csv',
    ROOT/'数据/raw/adc13_supplementary_sources/mct-24-0588_supplementary_methods_and_tables_suppsmt.pdf',
    ROOT/'数据/raw/second_batch_release_species/PMDA_belantamab_mafodotin_English_report.pdf',
]
hashes = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
core = read(inputs[0])
strict = [r for r in core if r['strict_human_transfectant_benchmark_eligible']=='yes']
paired = []
for r in core:
    er0, eri = float(r['no_inhibitor_ER_same_run']), float(r['matched_inhibitor_ER'])
    paired.append({
        'entity_id':r['candidate_id'], 'entity_name':r['entity_name'],
        'baseline_ER_matched_run':er0, 'inhibitor_ER_matched_run':eri,
        'inhibitor_ER_fraction_of_baseline':eri/er0,
        'inhibitor_associated_ER_reduction_fraction':1-eri/er0,
        'baseline_to_inhibitor_ER_fold':er0/eri,
        'concentration':r['concentration'], 'inhibitor':r['inhibitor'],
        'source_experiment_cluster':r['source_experiment_cluster'], 'endpoint_context':r['endpoint_context'],
        'gene_isolated_transfected_assay':'no_endogenous_Caco2_plus_selective_inhibitor_and_competence' if r['entity_name']=='MMAE' else 'human_Pgp_transfectant',
        'v01_score':float(r['v01_score']), 'v03_score':float(r['v03_score']),
        'full_v04_accepted':r['full_v04_accepted'], 'full_v04_call':r['full_v04_call'],
        'analysis_status':'purposefully_selected_known_positive_cases;no_cross_source_ER_ranking_or_accuracy',
    })
write('direct_positive_matched_run_results.csv', paired)

metrics=[]
# Native full selector determines each pool's k. No tuning to cases.
for pool,poolrows in [('strict_human_transfectant',strict),('source_context_positive_anchors',core)]:
    k = sum(r['full_v04_accepted']=='yes' for r in poolrows)
    selectors = {
        'full_native': [r for r in poolrows if r['full_v04_accepted']=='yes'],
        'v01_own_margin_equal_k': sorted(poolrows,key=lambda r:(-abs(float(r['v01_score'])-.5),r['candidate_id']))[:k],
        'v03_own_margin_equal_k': sorted(poolrows,key=lambda r:(-abs(float(r['v03_score'])-.5),r['candidate_id']))[:k],
        'v01_forced_all':poolrows[:], 'v03_forced_all':poolrows[:],
    }
    for strategy,selected in selectors.items():
        scorekey='v03_score' if strategy.startswith('v03') else 'v01_score'
        tp=sum(float(r[scorekey])>=.5 for r in selected)
        fn=len(selected)-tp
        metrics.append({
            'analysis_population':pool,'native_full_k':k,
            'strategy':strategy, 'eligible_unique_compounds_n':len(poolrows),
            'eligible_positive_n':len(poolrows), 'eligible_negative_n':0,
            'accepted_n':len(selected), 'coverage':len(selected)/len(poolrows),
            'accepted_entity_ids':'|'.join(r['candidate_id'] for r in selected),
            'TP_count':tp, 'FN_count':fn, 'TN_count':0, 'FP_count':0,
            'BAcc':None, 'NPV':None, 'specificity':None, 'performance_CI':None,
            'metric_limitation':'No negative class; no low-risk predictions. Strict native accepts zero; source-context pool includes endogenous Caco2. Counts are not independent accuracy validation.',
        })
write('direct_positive_strategy_counts.csv', metrics)

# Only the same Table S5 P-gp arm is compared. No main-source or BCRP-arm pooling.
obs = [r for r in read(inputs[1]) if r['transporter_gene']=='ABCB1' and r['compound'] in ['SN-38','DXd','Exatecan']]
preds={r['entity_name']:r for r in read(inputs[2]) if r['score_status']=='scored_frozen_model' and r['entity_name'] in ['SN-38','DXd','Exatecan']}
top1=[]
for name in ['SN-38','DXd','Exatecan']:
    pair={r['inhibitor']:r for r in obs if r['compound']==name}
    b,i=pair['none'],pair['verapamil']
    assert b['cell_system']==i['cell_system']
    assert b['test_concentration_uM']==i['test_concentration_uM']=='10'
    assert b['incubation_hours']==i['incubation_hours']=='2'
    er0,eri=float(b['efflux_ratio']),float(i['efflux_ratio'])
    p=preds[name]
    top1.append({
        'entity_name':name,'entity_id':b['entity_id'], 'baseline_record_id':b['record_id'], 'inhibitor_record_id':i['record_id'],
        'source_id':b['source_id'],'source_table':'Supplementary Table S5 P-gp arm',
        'cell_system':b['cell_system'],'concentration_uM':10,'incubation_hours':2,
        'inhibitor':'verapamil','inhibitor_concentration_uM':None,
        'ER_no_inhibitor_author':er0,'ER_with_verapamil_author':eri,
        'ER_with_inhibitor_fraction_of_baseline':eri/er0,
        'inhibitor_associated_ER_reduction_fraction':1-eri/er0,
        'ER_reduction_fold':er0/eri,
        'baseline_Papp_AB_mean':float(b['papp_ab_mean_1e6_cm_s']),
        'baseline_Papp_AB_SD':float(b['papp_ab_sd']),
        'baseline_Papp_BA_mean':float(b['papp_ba_mean_1e6_cm_s']),
        'baseline_Papp_BA_SD':float(b['papp_ba_sd']),
        'v01_frozen_score':float(p['v0_1_risk_score']), 'v03_frozen_score':float(p['v0_3_risk_score']),
        'full_v04_accepted':p['accepted_prediction'],'full_v04_call':p['risk_call'],
        'identity_seen_label_status':'previously_curated;posthoc_not_blind',
        'endpoint_limitation':'Caco2 multitransporter ER; missing S5 inhibitor concentration; score is not predicted ER or Pgp-only transport.',
    })
write('top1_same_source_Pgp_arm_results.csv',top1)

comparisons=[]
correlations=[]
for endpoint in ['ER_no_inhibitor_author','inhibitor_associated_ER_reduction_fraction']:
    for score in ['v01_frozen_score','v03_frozen_score']:
        rho=spearman([r[endpoint] for r in top1],[r[score] for r in top1])
        correlations.append({'source_table':'OBI992 Supplementary Table S5 P-gp arm','n_unique':3,'endpoint':endpoint,'score_channel':score,'descriptive_spearman_rho':rho,'CI':None,'P_value':None,'accuracy':None,'interpretation':'Descriptive order check only; n=3 purposive known compounds, no performance/causal claim.'})
        for a,b in itertools.combinations(top1,2):
            de,ds=a[endpoint]-b[endpoint],a[score]-b[score]
            se=(de>0)-(de<0);ss=(ds>0)-(ds<0)
            comparisons.append({'entity_A':a['entity_name'],'entity_B':b['entity_name'],'endpoint':endpoint,'endpoint_A_minus_B':de,'score_channel':score,'score_A_minus_B':ds,'pair_order_status':'same_order' if se==ss else 'opposite_order','analysis_unit':'one_compound_pair_within_same_source;not_independent_accuracy_observation'})
write('top1_same_source_pairwise_order_checks.csv',comparisons)
write('top1_same_source_descriptive_spearman.csv',correlations)
summary={
    'status':'finite_retrospective_descriptive_checks_executed',
    'ADC_source_context_positive_unique_n':3,'ADC_source_context_eligible_negative_n':0,'source_context_full_native_k':1,
    'ADC_strict_human_transfectant_positive_n':2,'ADC_strict_human_transfectant_negative_n':0,'strict_human_transfectant_full_native_k':0,
    'direct_strategy_counts':metrics,'top1_same_source_results':top1,'top1_spearman':correlations,
    'same_system_parent_release_delta_available':False,
    'protected_inputs_unchanged':all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in hashes.items()),
    'input_sha256':hashes,'model_fit_or_threshold_change':False,
    'new_wet_lab_experiment':False,'new_blind_source_validation':False,
}
assert summary['protected_inputs_unchanged']
(HERE/'descriptive_validation_results.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
report='''# ADC组件运输验证：实际执行结果（2026-10-12）

本轮在冻结模型和已有来源上执行了有限回顾性计算：纠正来源实验配对、同覆盖率接受名单、同一Table S5内三种TOP1组分的运输排序。没有重新训练、改阈值或生成实验室数据。

## 1. 三个直接运输阳性病例的正确实验配对

|实体|同run ER：无抑制剂→抑制剂|ER降低比例|来源体系|
|---|---:|---:|---|
|DM1|5.90→1.10|81.36%|human MDR1 MDCKII，0.5 uM，PSC833 10 uM|
|cys-mcMMAF|34→1.1|96.76%|human P-gp MDCKII，100 nM，verapamil 60 uM|
|MMAE|28→1.3|95.36%|Caco2，10 uM，PSC833 10 uM；非gene-isolated|

MMAE旧锚49→1.3混用了100 uM浓度扫描与10 uM抑制剂实验。新分析只用来源同run的28→1.3；旧文件保持原样。三个来源的细胞体系、浓度和材料不同，不跨来源排列ER或估计相关性。

严格human P-gp转染池只有DM1和cys-mcMMAF两个阳性，完整选择器接受0/2。同coverage基线也必须保留0个，不额外选择1个。该池没有可评估的选择性准确率。

扩展的source-context阳性锚池另含内源Caco2的MMAE，完整选择器接受1/3（MMAE）；v0.1自身margin同1名额接受DM1；v0.3自身margin接受cys-mcMMAF。强制分类在各自池中将全部病例判高。两个池均没有负类，也没有低风险输出，BAcc、NPV和specificity均为null。1个接受病例不能包装为100%准确率。三锚池不是统一human-transfectant benchmark。

## 2. Table S5同来源TOP1组分比较

|实体|P-gp arm ER：无抑制剂→verapamil|ER降低比例|v0.1冻结分数|v0.3冻结分数|
|---|---:|---:|---:|---:|
|SN-38|20.0→10.5|47.50%|0.85405|0.52877|
|DXd|22.1→13.8|37.56%|0.82798|0.63007|
|Exatecan|2.7→1.7|37.04%|0.78437|0.56162|

同一Caco2 C2BBe1体系、10 uM、2 h。仅比较Table S5 P-gp arm，未混入Ko143 arm或其他论文。两通道与无抑制剂ER的Spearman均为0.5；与ER降低比例的Spearman分别为1.0和−0.5。以上只是3个既知化合物的描述性排序，未给CI、p值或准确率；同组三个pairwise comparisons也不是3个独立样本。

Caco2表达多个转运体，Table S5未报告抑制剂浓度，ER降低比例不等于ABCB1运输占比。DXd低A→B通透性0.12±0.10也限制ER解释。三者均被完整选择器接受高风险，但Exatecan的“相对低ER”不等于可赋予阴性底物标签。该实验比较不同TOP1组分，未证明它们具有母体→释放物化学关系，也未验证双payload协同。

## 3. 可复算与本轮边界

```bash
cd '/Users/hongyanliu/Desktop/论文/adc药物项目/预测器2'
python3 outputs/20261012_adc_completion/design_audit/build_eligibility_audit.py
python3 outputs/20261012_adc_completion/design_audit/run_descriptive_validation.py
python3 outputs/20261012_adc_completion/design_audit/verify_numeric_results.py
```

采用Python标准库，复用冻结13结构分数。CSV保留逐实体、逐pair、接受名单、null指标和来源实验族；JSON记录输入SHA256和无旧文件变更结果。另一不调用主helper的闭式Spearman及池/空值审核共12项通过，回执见independent_numeric_verification.json；它是数值复核，不是新生物验证。严格同体系母体—释放物运输配对仍为0。下一步缺口是材料身份明确且同体系直接测定的释放物负例及配对真值。

原始来源：DM1见[PMDA Kadcyla CTD](https://www.pmda.go.jp/drugs/2013/P201300137/450045000_22500AMX01816_I100_2.pdf)表2.6.5.15-1；cys-mcMMAF见[PMDA belantamab](https://www.pmda.go.jp/files/000280862.pdf)4.5.3；MMAE见本地PMDA Padcev CTD Study8004（PDF p76 / printed39）。TOP1见[OBI992原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11791482/)Supplementary Table S5，本地原补充PDF已重新抽取核对。
'''
(HERE/'本轮ADC有限验证结果.md').write_text(report)
print(json.dumps({'status':summary['status'],'source_context_n':3,'source_context_k':1,'strict_human_transfectant_n':2,'strict_human_transfectant_k':0,'negative_n':0,'top1_spearman':correlations,'protected_inputs_unchanged':summary['protected_inputs_unchanged']},ensure_ascii=False,indent=2))
