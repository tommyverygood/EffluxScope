#!/usr/bin/env python3
"""Independent formula and population/null checks, without importing main helpers."""
import csv
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
read=lambda fn:list(csv.DictReader((HERE/fn).open(encoding='utf-8-sig')))
checks=[]

def check(name,passed):
    checks.append({'check':name,'passed':bool(passed)})
    assert passed,name

strict=read('adc_strict_human_transfectant_positive_anchors.csv')
extended=read('adc_source_context_positive_anchors.csv')
check('Strict human-transfectant population exactly DM1 and cys-mcMMAF',{r['entity_name'] for r in strict}=={'DM1','cys-mcMMAF'})
check('Source-context population adds only endogenous-Caco2 MMAE',{r['entity_name'] for r in extended}=={'DM1','cys-mcMMAF','MMAE'})
metrics=read('direct_positive_strategy_counts.csv')
check('Each population has five strategies',len(metrics)==10 and all(sum(r['analysis_population']==g for r in metrics)==5 for g in ['strict_human_transfectant','source_context_positive_anchors']))
check('All strict equal-coverage strategies select zero',all(r['accepted_n']=='0' and r['native_full_k']=='0' for r in metrics if r['analysis_population']=='strict_human_transfectant' and 'forced' not in r['strategy']))
check('All empty denominator metrics remain blank in CSV',all(all(r[x]=='' for x in ['BAcc','NPV','specificity','performance_CI']) for r in metrics))
expected={('SN-38','none'):20.0,('SN-38','verapamil'):10.5,('DXd','none'):22.1,('DXd','verapamil'):13.8,('Exatecan','none'):2.7,('Exatecan','verapamil'):1.7}
observations=list(csv.DictReader((ROOT/'数据/processed/obi992_caco2_payload_transport_numeric_v1.csv').open(encoding='utf-8-sig')))
check('Source S5 arm matches independent literal primary-table readings',all(float(next(r for r in observations if r['compound']==name and r['transporter_gene']=='ABCB1' and r['inhibitor']==inh)['efflux_ratio'])==value for (name,inh),value in expected.items()))
rows=read('top1_same_source_Pgp_arm_results.csv')
correlations=read('top1_same_source_descriptive_spearman.csv')
for endpoint in ['ER_no_inhibitor_author','inhibitor_associated_ER_reduction_fraction']:
    for score in ['v01_frozen_score','v03_frozen_score']:
        xs=[float(r[endpoint]) for r in rows]
        ys=[float(r[score]) for r in rows]
        # Independent no-ties closed-form Spearman, no main helper.
        rx=[sorted(xs).index(x)+1 for x in xs]
        ry=[sorted(ys).index(y)+1 for y in ys]
        rho=1-6*sum((x-y)**2 for x,y in zip(rx,ry))/(3*(3**2-1))
        actual=float(next(r for r in correlations if r['endpoint']==endpoint and r['score_channel']==score)['descriptive_spearman_rho'])
        check(f'Independent closed-form Spearman {endpoint} {score}',abs(rho-actual)<1e-12)
check('Exactly 12 pairwise checks share three compounds',len(read('top1_same_source_pairwise_order_checks.csv'))==12)
check('MMAE uses 28 baseline and not 49',float(next(r for r in read('direct_positive_matched_run_results.csv') if r['entity_name']=='MMAE')['baseline_ER_matched_run'])==28)
proof={}
for fn,pattern in [('obi992_primary_supplement.txt','Supplementary Table S5'),('belantamab_primary_text.txt','4.5.3 Transporters')]:
    source=(HERE/fn).read_text()
    proof[fn]={'page_1based_first_match':next(i+1 for i,t in enumerate(source.split('\f')) if pattern in t),'sha256':hashlib.sha256((HERE/fn).read_bytes()).hexdigest()}
output={'status':'PASS','checks_n':len(checks),'all_checks_passed':True,'independent_statistics_method':'no-ties Spearman closed-form versus main average-rank Pearson; literal primary S5 Pgp arm readings','checks':checks,'primary_text_location_proof':proof,'no_model_fit_or_rescoring':True,'scope':'numerical and protocol-layer verification, not independent biological validation'}
(HERE/'independent_numeric_verification.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks_n':len(checks)},ensure_ascii=False))
