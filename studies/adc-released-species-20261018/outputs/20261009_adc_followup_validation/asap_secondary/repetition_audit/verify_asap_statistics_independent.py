#!/usr/bin/env python3
"""Independent rank/statistic reproduction from legacy predictions, no root helpers."""
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import statistics

import numpy as np
from scipy.stats import rankdata

HERE=Path(__file__).resolve().parent
OUT=HERE.parent
ROOT=HERE.parents[3]
SOURCE=ROOT/'数据/processed/abcb1_asap2025_sourceheldout_predictions_v1.csv'
FULL=ROOT/'数据/processed/abcb1_selective_v0_4_asap2025_posthoc_predictions_v1.csv'

def read(p):
    with p.open(encoding='utf-8-sig',newline='') as h:return list(csv.DictReader(h))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rho(a,b):return float(np.corrcoef(rankdata(a),rankdata(b))[0,1])
def auc(y,s):
    pos=s[y==1];neg=s[y==0]
    return float(np.mean([(p>n)+.5*(p==n) for p in pos for n in neg]))
def binary(y,s):
    pred=s>=.5
    tp=int(np.count_nonzero((y==1)&pred));tn=int(np.count_nonzero((y==0)&~pred))
    fp=int(np.count_nonzero((y==0)&pred));fn=int(np.count_nonzero((y==1)&~pred))
    return dict(tp=tp,tn=tn,fp=fp,fn=fn,balanced_accuracy=.5*(tp/(tp+fn)+tn/(tn+fp)) if (tp+fn)*(tn+fp) else None,
        npv=tn/(tn+fn) if tn+fn else None)

def main():
    before=digest(SOURCE)
    rows=read(SOURCE)
    full={r['observation_id']:r for r in read(FULL)}
    protocol=json.loads((OUT/'protocol.json').read_text())
    bykey=defaultdict(list)
    for r in rows:bykey[r['standardized_inchikey']].append(r)
    compound,measurements=[],[]
    for key,rs in sorted(bykey.items()):
        r=rs[0]
        common=dict(uid=key,scaffold=r['murcko_scaffold'] or 'acyclic::'+key,
            p=float(r['v0_1_calibrated_probability']),q=float(r['v0_3_candidate_probability']))
        compound.append(common|dict(y=statistics.median(float(x['log2_efflux_ratio']) for x in rs)))
        measurements.extend(common|dict(uid=x['observation_id'],y=float(x['log2_efflux_ratio'])) for x in rs)
    metrics={ (r['cohort'],r['model']):r for r in read(OUT/'model_metrics.csv') }
    contrasts={r['cohort']:r for r in read(OUT/'paired_rho_contrasts.csv')}
    selected_metrics={(r['cohort'],r['policy']):r for r in read(OUT/'matched_count_metrics.csv')}
    ledger={r['key']:r for r in read(OUT/'compound_aggregation_ledger.csv')}
    checks=[]; own=[]; selections=[]
    for r in compound:
        c=ledger[r['uid']]
        for a,b in [('p','primary'),('q','auxiliary'),('y','log2_er')]:
            assert abs(r[a]-float(c[b]))<1e-12
        assert r['scaffold']==c['cluster']
    for name,cohort in [('unique_compound_56',compound),('measurement_rows_62',measurements)]:
        y=np.array([r['y'] for r in cohort]);p=np.array([r['p'] for r in cohort]);q=np.array([r['q'] for r in cohort])
        binary_y=(y>=1).astype(int)
        clusters=sorted({r['scaffold'] for r in cohort})
        positions={c:[i for i,r in enumerate(cohort) if r['scaffold']==c] for c in clusters}
        rng=np.random.default_rng(protocol['seed'])
        draw_counts=rng.multinomial(len(clusters),np.repeat(1/len(clusters),len(clusters)),size=protocol['bootstrap_draws'])
        bootp=[];bootq=[];delta=[]
        for weights in draw_counts:
            picked=[]
            for cluster,n in zip(clusters,weights):
                for _ in range(int(n)):picked.extend(positions[cluster])
            rp=rho(p[picked],y[picked]);rq=rho(q[picked],y[picked])
            bootp.append(rp);bootq.append(rq);delta.append(rq-rp)
        for model,scores,bootstrap in [('v01',p,bootp),('v03',q,bootq)]:
            lo,hi=np.percentile(bootstrap,[2.5,97.5])
            computed=dict(population_n=len(cohort),scaffold_clusters=len(clusters),
                spearman_rho=rho(scores,y),spearman_ci95_low=float(lo),spearman_ci95_high=float(hi),
                spearman_valid_draws=len(bootstrap),auc_weak_er2=auc(binary_y,scores),**binary(binary_y,scores))
            reported=metrics[(name,model)]
            for field,value in computed.items():
                err=abs(float(reported[field])-value) if value is not None else None
                checks.append(dict(cohort=name,model=model,field=field,computed=value,reported=reported[field],
                    absolute_difference=err,pass_check=err<1e-12 if err is not None else reported[field]==''))
            own.append(dict(cohort=name,model=model,**computed))
        c=contrasts[name];lo,hi=np.percentile(delta,[2.5,97.5])
        values=dict(point_delta=rho(q,y)-rho(p,y),ci95_low=float(lo),ci95_high=float(hi),valid_draws=len(delta))
        for field,value in values.items():
            err=abs(float(c[field])-value)
            checks.append(dict(cohort=name,model='v03_minus_v01',field=field,computed=value,
                reported=c[field],absolute_difference=err,pass_check=err<1e-12))
        by_source_key=defaultdict(list)
        for source_row in full.values():by_source_key[source_row['standardized_inchikey']].append(source_row)
        full_mask=[]
        for sample in cohort:
            if name=='unique_compound_56':
                accepted={x['accepted_prediction'] for x in by_source_key[sample['uid']]}
                assert len(accepted)==1
                full_mask.append(accepted.pop()=='yes')
            else:
                full_mask.append(full[sample['uid']]['accepted_prediction']=='yes')
        k=sum(full_mask)
        for policy,scores in [('v01_own_margin',p),('v03_own_margin',q),('full_v04_native',p)]:
            if policy=='full_v04_native':
                selected=np.flatnonzero(full_mask)
            else:
                selected=np.array(sorted(range(len(cohort)),key=lambda i:(-abs(float(scores[i])-.5),cohort[i]['uid']))[:k])
            y_selected=binary_y[selected]
            computed=dict(population_n=len(cohort),coverage=len(selected)/len(cohort),selected_n=len(selected),
                positive_n=int((y_selected==1).sum()),negative_n=int((y_selected==0).sum()),
                **binary(y_selected,scores[selected]))
            reported=selected_metrics[(name,policy)]
            for field,value in computed.items():
                error=abs(float(reported[field])-value) if value is not None else None
                checks.append(dict(cohort=name,model=policy,field=field,computed=value,reported=reported[field],
                    absolute_difference=error,pass_check=error<1e-12 if error is not None else reported[field]==''))
            assert computed['selected_n']==5 and computed['negative_n']==0 and computed['fn']==0
            assert computed['balanced_accuracy'] is None and computed['npv'] is None
            selections.append(dict(cohort=name,policy=policy,selected_ids=[cohort[i]['uid'] for i in selected],**computed))
        # The legacy and later joined tables agree exactly; no model refit is used.
        assert all(r['v0_1_calibrated_probability']==full[r['observation_id']]['v0_1_calibrated_probability']
            and r['v0_3_candidate_probability']==full[r['observation_id']]['v0_3_candidate_probability']
            and r['log2_efflux_ratio']==full[r['observation_id']]['log2_efflux_ratio'] for r in rows)
    assert all(c['pass_check'] for c in checks)
    assert digest(SOURCE)==before
    with (HERE/'independent_statistic_checks.csv').open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(checks[0]));w.writeheader();w.writerows(checks)
    result=dict(status='PASS_INDEPENDENT_STATISTIC_REPRODUCTION',total_checks=len(checks),
        all_checks_passed=True,legacy_input_sha256=before,legacy_input_unchanged=True,
        independent_implementation='scipy rankdata + numpy Pearson correlation, explicit pairwise AUC, explicit confusion counts; no run_asap_replication helper imports',
        primary_ledger_matched=True,source_vs_joined_score_agreement=True,
        seed=protocol['seed'],bootstrap_draws=protocol['bootstrap_draws'],
        bootstrap_design='same recorded multinomial whole-scaffold draws, preserving compound observations together and pairing models',
        metrics=own,matched_selection=selections)
    (HERE/'independent_statistics_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),'all_passed':True}))

if __name__=='__main__':main()
