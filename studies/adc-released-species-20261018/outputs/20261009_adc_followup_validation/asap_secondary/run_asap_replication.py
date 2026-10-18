#!/usr/bin/env python3
"""Paired, scaffold-aware retrospective model comparison on ASAP 2025."""
from pathlib import Path
import csv
import hashlib
import json
from collections import defaultdict
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
SOURCE = ROOT / "数据/processed/abcb1_selective_v0_4_asap2025_posthoc_predictions_v1.csv"
SEED, N_BOOT = 20261009, 2000

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(p):
    with p.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def clean(value):
    if isinstance(value, dict): return {k: clean(v) for k,v in value.items()}
    if isinstance(value, list): return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)): return float(value) if np.isfinite(value) else None
    if isinstance(value, (np.integer,)): return int(value)
    return value

def write_csv(name, data):
    fields = list(dict.fromkeys(k for r in data for k in r))
    with (OUT/name).open("w", encoding="utf-8-sig", newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows([clean(r) for r in data])

def interval(values):
    v=np.array(values,dtype=float); v=v[np.isfinite(v)]
    return {"ci95_low": float(np.quantile(v,.025)) if len(v) else None,
            "ci95_high": float(np.quantile(v,.975)) if len(v) else None,
            "valid_draws": int(len(v))}

def direction_metrics(y,pred):
    tp=int(np.sum((y==1)&pred)); tn=int(np.sum((y==0)&~pred))
    fp=int(np.sum((y==0)&pred)); fn=int(np.sum((y==1)&~pred))
    sen=tp/(tp+fn) if tp+fn else None
    spec=tn/(tn+fp) if tn+fp else None
    return {"selected_n":len(y),"positive_n":int(np.sum(y==1)),"negative_n":int(np.sum(y==0)),
            "tp":tp,"tn":tn,"fp":fp,"fn":fn,
            "balanced_accuracy":(sen+spec)/2 if sen is not None and spec is not None else None,
            "npv":tn/(tn+fn) if tn+fn else None,
            "bacc_unavailable_reason":"" if sen is not None and spec is not None else "selected_subset_missing_true_class"}

def evaluate(name, cohort):
    ids=[r["id"] for r in cohort]
    clusters=sorted({r["cluster"] for r in cohort}); lookup={c:i for i,c in enumerate(clusters)}
    group=np.array([lookup[r["cluster"]] for r in cohort])
    index_groups=[np.flatnonzero(group==i) for i in range(len(clusters))]
    target=np.array([r["log2_er"] for r in cohort])
    y=(target>=1).astype(int)
    p=np.array([r["primary"] for r in cohort]); q=np.array([r["auxiliary"] for r in cohort])
    full=np.array([r["full"] for r in cohort])
    assert len(set(ids))==len(ids) and np.isfinite(target).all()
    draws=np.random.default_rng(SEED).multinomial(len(clusters),np.full(len(clusters),1/len(clusters)),size=N_BOOT)
    boot={"v01":[],"v03":[],"v03_minus_v01":[]}
    for counts in draws:
        ix=np.concatenate([np.tile(ix,int(n)) for ix,n in zip(index_groups,counts) if n])
        a=float(spearmanr(p[ix],target[ix]).statistic); b=float(spearmanr(q[ix],target[ix]).statistic)
        boot["v01"].append(a); boot["v03"].append(b); boot["v03_minus_v01"].append(b-a)
    output=[]
    for channel,scores in [("v01",p),("v03",q)]:
        rho=float(spearmanr(scores,target).statistic)
        output.append({"cohort":name,"role":"primary_unique_compound" if name=="unique_compound_56" else "historical_row_sensitivity",
                       "model":channel,"population_n":len(cohort),"scaffold_clusters":len(clusters),
                       "spearman_rho":rho,**{f"spearman_{k}":v for k,v in interval(boot[channel]).items()},
                       "auc_weak_er2":float(roc_auc_score(y,scores)),**direction_metrics(y,scores>=.5)})
    a=float(spearmanr(p,target).statistic); b=float(spearmanr(q,target).statistic)
    contrast={"cohort":name,"metric":"paired_delta_Spearman","lhs":"v03","rhs":"v01","point_delta":b-a,**interval(boot["v03_minus_v01"]),
              "bootstrap_unit":"whole_Murcko_scaffold","shared_draws":True}
    selected=[]; k=int(full.sum())
    for policy,scores,mask in [("v01_own_margin",p,None),("v03_own_margin",q,None),("full_v04_native",p,full)]:
        if mask is None:
            ix=sorted(range(len(cohort)),key=lambda i:(-abs(scores[i]-.5),ids[i]))[:k]
            mask=np.zeros(len(cohort),bool); mask[ix]=True
        selected.append({"cohort":name,"policy":policy,"population_n":len(cohort),"coverage":float(mask.mean()),
                         **direction_metrics(y[mask],scores[mask]>=.5)})
    return output,contrast,selected

def main():
    OUT.mkdir(exist_ok=True)
    before=sha(SOURCE)
    protocol={"source":str(SOURCE.relative_to(ROOT)),"sha256":before,"role":"retrospective_secondary_replication_not_ADC_validation",
              "already_observed_outcomes":True,"no_refit_or_recalibration":True,
              "primary_cohort":"unique exact standardized compound, median log2ER across existing observations; scores must agree",
              "secondary_cohort":"historical62 measurement rows, whole scaffold resampling retains all repeats together",
              "primary_test":"paired delta Spearman(v03 vs v01) using same2000 Murcko scaffold draws",
              "binary_ER2":"secondary assay-specific sensitivity endpoint without parental/inhibitor confirmation",
              "matched_count":"full native accepted count, own-margin deterministic UID ties; BAcc unavailable if selected class absent",
              "seed":SEED,"bootstrap_draws":N_BOOT,
              "no_new_source_collection":True,"scope":"One vendor MDR1-MDCK dataset; not completeADC, matched parent-release, or universal substrate status"}
    (OUT/"protocol.json").write_text(json.dumps(protocol,ensure_ascii=False,indent=2)+"\n")
    source=read(SOURCE)
    assert len(source)==62
    groups=defaultdict(list)
    for r in source:
        assert r["exact_reference_identity"]=="no" and r["v0_3_exact_development_identity_match"]=="no"
        assert abs(float(r["v0_1_calibrated_probability"])-float(r["primary_calibrated_probability"]))<1e-10
        assert abs(float(r["v0_3_candidate_probability"])-float(r["direct_panel_probability"]))<1e-10
        groups[r["standardized_inchikey"]].append(r)
    assert len(groups)==56
    unique=[]; observations=[]; audits=[]
    for key,rs in sorted(groups.items()):
        assert len({r["v0_1_calibrated_probability"] for r in rs})==1 and len({r["v0_3_candidate_probability"] for r in rs})==1
        assert len({r["murcko_scaffold"] for r in rs})==1 and len({r["accepted_prediction"] for r in rs})==1
        log=np.array([float(r["log2_efflux_ratio"]) for r in rs]); r=rs[0]
        common={"key":key,"cluster":r["murcko_scaffold"] or "acyclic::"+key,
                "primary":float(r["v0_1_calibrated_probability"]),"auxiliary":float(r["v0_3_candidate_probability"]),
                "full":r["accepted_prediction"]=="yes"}
        unique.append({**common,"id":key,"log2_er":float(np.median(log)),"n_observation_rows":len(rs)})
        observations.extend({**common,"id":x["observation_id"],"log2_er":float(x["log2_efflux_ratio"]),"n_observation_rows":1} for x in rs)
        if len(rs)>1:
            audits.append({"key":key,"observation_rows":len(rs),"er_min":float(2**log.min()),"er_max":float(2**log.max()),
                           "crosses_weak_ER2":bool(log.min()<1<=log.max()),"median_log2ER":float(np.median(log)),
                           "score_agreement":True,"independent_biological_replicate_status":"not_reported_not_assumed"})
    outputs=[]; contrasts=[]; selections=[]
    for name,cohort in [("unique_compound_56",unique),("measurement_rows_62",observations)]:
        o,c,s=evaluate(name,cohort); outputs.extend(o); contrasts.append(c); selections.extend(s)
    assert sha(SOURCE)==before
    write_csv("model_metrics.csv",outputs); write_csv("paired_rho_contrasts.csv",contrasts)
    write_csv("matched_count_metrics.csv",selections); write_csv("compound_aggregation_ledger.csv",unique)
    write_csv("duplicate_groups.csv",audits)
    (OUT/"completion.json").write_text(json.dumps(clean({"source_unchanged":True,"input_rows":62,"unique_compounds":56,"duplicate_groups":4,
                                                        "weak_ER2_crossing_duplicate_group_count":sum(r["crosses_weak_ER2"] for r in audits),
                                                        "metrics":outputs,"paired_deltas":contrasts,"matched_selection":selections}),ensure_ascii=False,indent=2)+"\n")
    report=["# ASAP来源的独立化合物与成对模型比较", "",
            "这是既有ASAP2025数据的回顾性再分析，未新增来源或ADC运输观测。62条可计算观察对应56个唯一化合物；主分析先按精确standardized身份聚合median log2ER，避免重复观察作为独立分子。62行分析保留作历史敏感性，bootstrap仍以完整Murcko scaffold为单位。",
            "", "|人群|模型|Spearman rho [95% CI]|弱ER≥2 BAcc（全部分子）|", "|---|---|---|---:|"]
    for r in outputs:
        report.append(f"|{r['cohort']}|{r['model']}|{r['spearman_rho']:.3f} [{r['spearman_ci95_low']:.3f}, {r['spearman_ci95_high']:.3f}]|{r['balanced_accuracy']:.3f}|")
    report.extend(["", "|人群|v0.3−v0.1 Δrho|共同scaffold抽样95% CI|", "|---|---:|---|"])
    for r in contrasts:
        report.append(f"|{r['cohort']}|{r['point_delta']:+.3f}|[{r['ci95_low']:+.3f}, {r['ci95_high']:+.3f}]|")
    report.extend(["", "full只接受5个，全部为真实弱ER阳性，无法估计负类性能或BAcc。相同5名额比较见matched_count_metrics.csv；缺失真类时保留指标为空，不能填成完美BAcc。",
                   "", "四个重复组中，一组两个观察的ER为1.944与2.641，跨ER=2边界。来源未报告这些行的独立培养日、样本编号或生物重复类型；不能把R2当成新增生物独立验证。",
                   "", "成对CI是固定两组模型分数、一个vendor来源内的scaffold抽样不确定性，不包含重新训练、来源转移或ADC领域偏移。没有parental或抑制剂臂，ER=2只作为测定特定的弱方向指标。即使相关性为正，也不构成ADC实际释放物预测准确率。"])
    (OUT/"ASAP复核.md").write_text("\n".join(report)+"\n")
    print(json.dumps(clean({"cohorts":[56,62],"paired_deltas":contrasts,"matched_selection":selections}),ensure_ascii=False))

if __name__=="__main__": main()
