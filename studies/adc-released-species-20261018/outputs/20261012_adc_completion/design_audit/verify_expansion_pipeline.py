#!/usr/bin/env python3
"""Independent DMH preprocessing/features/score reconstruction; no project helpers."""
import csv
import hashlib
import json
import math
import platform
import sys
import warnings
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
read=lambda p:list(csv.DictReader(p.open(encoding='utf-8-sig')))
print('Independent pipeline importing matching numerical packages',flush=True)
import joblib
import numpy as np
import rdkit
import sklearn
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import Crippen,Descriptors,Lipinski,rdFingerprintGenerator,rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog('rdApp.warning')

PATHS=[ROOT/'outputs/20261012_adc_completion/scoring/species_predictions.csv',
       ROOT/'outputs/20261012_adc_completion/scoring/historical_reproduction.csv',
       ROOT/'outputs/20261012_adc_completion/scoring/completion.json',
       ROOT/'outputs/20261009_adc_pair_experiments/pairs/results_curated_v2/entities_scored_or_hold.csv',
       ROOT/'模型/ABCB1_screener_v0_1/ABCB1_screener_v0_1.joblib',
       ROOT/'模型/ABCB1_assay_aware_candidate_v0_3_dev/ABCB1_assay_aware_candidate_v0_3_dev.joblib',
       ROOT/'模型/ABCB1_selective_screener_v0_4_dev/ABCB1_selective_uncertainty_v0_4_dev.joblib']
hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in PATHS}
root_predictions=read(PATHS[0])
r=next(x for x in root_predictions if x['entity_id']=='NEWF_CALICH_DMH')
checks=[]
def check(name, passed,details=None):
    checks.append({'check':name,'passed':bool(passed),'details':details})
    assert passed,name

source=Chem.MolFromSmiles(r['source_smiles'])
source_key=Chem.MolToInchiKey(source)
source_formula=rdMolDescriptors.CalcMolFormula(source)
check('RDKit independently confirms exact source key',source_key==r['source_inchikey'],source_key)
check('Source molecular formula',source_formula==r['source_formula'],source_formula)
mol=rdMolStandardize.Cleanup(source)
mol=rdMolStandardize.FragmentParent(mol,skipStandardize=True)
mol=rdMolStandardize.Uncharger(canonicalOrder=True).uncharge(mol)
enumerator=rdMolStandardize.TautomerEnumerator()
tautomers=enumerator.Enumerate(mol)
mol=enumerator.PickCanonical(tautomers)
Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
key=Chem.MolToInchiKey(mol)
status=getattr(tautomers.status,'name',str(tautomers.status))
iso=Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)
check('Independent frozen pipeline standardized key',key==r['standardized_inchikey'],key)
check('Independent frozen canonical SMILES',iso==r['standardized_isomeric_smiles'])
check('Independent enumeration Completed',status==r['tautomer_status']=='Completed',status)

generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
fp=generator.GetFingerprint(mol)
fp_array=np.zeros((2048,),dtype=np.float32)
DataStructs.ConvertToNumpyArray(fp,fp_array)
descriptors=np.asarray([Descriptors.MolWt(mol),Crippen.MolLogP(mol),rdMolDescriptors.CalcTPSA(mol),Lipinski.NumHDonors(mol),Lipinski.NumHAcceptors(mol),Lipinski.NumRotatableBonds(mol),sum(a.GetFormalCharge() for a in mol.GetAtoms()),mol.GetNumHeavyAtoms()],dtype=np.float64)
print('Independent pipeline loading fixed models',flush=True)
with warnings.catch_warnings(record=True) as captured:
    warnings.simplefilter('always')
    primary,direct,selective=[joblib.load(p) for p in PATHS[4:]]
    def features(package):
        kind=package['feature_type']
        fingerprint=fp_array.reshape(1,-1)
        if kind=='morgan':return fingerprint
        scaled=package['descriptor_scaler'].transform(descriptors.reshape(1,-1))
        if kind=='physchem':return scaled
        assert kind in {'morgan_physchem','morgan_plus_physchem'}
        return np.hstack([fingerprint,scaled])
    raw=float(primary['estimator'].predict_proba(features(primary))[0,1])
    if primary.get('calibrator') is not None:
        clipped=float(np.clip(raw,1e-6,1-1e-6))
        logit=math.log(clipped/(1-clipped))
        p01=float(primary['calibrator'].predict_proba(np.asarray([[logit]]))[0,1])
    else:p01=raw
    p03=float(direct['estimator'].predict_proba(features(direct))[0,1])
    check('Independent v01 calibrated score',round(p01,8)==float(r['primary_calibrated_probability']),p01)
    check('Independent v03 direct score',round(p03,8)==float(r['direct_panel_probability']),p03)
    pkeys=set(primary['training_structure_keys'])
    dkeys={x['standardized_inchikey'] for x in direct['training_rows']}
    qkeys={source_key,key}
    qconn={x.split('-')[0] for x in qkeys}
    overlaps={
       'primary_exact_identity':bool(qkeys&pkeys),
       'direct_exact_identity':bool(qkeys&dkeys),
       'primary_connectivity_overlap':bool(qconn&{x.split('-')[0] for x in pkeys}),
       'direct_connectivity_overlap':bool(qconn&{x.split('-')[0] for x in dkeys}),
    }
    for name,value in overlaps.items():check('Source plus frozen '+name,str(value)==r[name] and not value,value)
    sim01=max(DataStructs.BulkTanimotoSimilarity(fp,primary['training_fingerprints']))
    sim03=max(DataStructs.BulkTanimotoSimilarity(fp,direct['training_fingerprints']))
    z=float(np.max(np.abs((descriptors-np.asarray(selective['physchem_center']))/np.asarray(selective['physchem_scale']))))
    check('Primary max similarity',round(sim01,6)==float(r['primary_max_tanimoto']),sim01)
    check('Direct max similarity',round(sim03,6)==float(r['direct_max_tanimoto']),sim03)
    check('Physicochemical robust z',round(z,8)==float(r['physchem_max_abs_robust_z']),z)
    bootstrap=[]
    for item in selective['bootstrap_models']:
        x=np.hstack([fp_array.reshape(1,-1),item['descriptor_scaler'].transform(descriptors.reshape(1,-1))])
        prob=float(item['estimator'].predict_proba(x)[0,1])
        clipped=float(np.clip(prob,1e-6,1-1e-6))
        bootstrap.append(float(primary['calibrator'].predict_proba(np.asarray([[math.log(clipped/(1-clipped))]]))[0,1]) if primary.get('calibrator') else prob)
    p10,median,p90=np.quantile(bootstrap,selective['bootstrap_probability_quantiles'])
    check('Independent bootstrap p10',round(float(p10),8)==float(r['bootstrap_probability_p10']),float(p10))
    check('Both reference support gates fail',sim01<selective['primary_similarity_q10'] and sim03<selective['direct_similarity_q10'])
    check('Physicochemical gate fails',z>selective['physchem_score_threshold_q99'])
    check('High-score bootstrap direction stable',p01>=selective['primary_probability_high_threshold'] and p10>=.5)
    expected_reasons=['outside_both_reference_panels','outside_physicochemical_envelope']
    check('DMH abstention two applicable gate reasons',json.loads(r['abstention_reasons_json'])==expected_reasons and r['risk_call']=='abstain' and r['accepted_prediction']=='no')

# No old panel rescore; verify all 52 old/current ledger entries from raw files.
old={x['entity_id']:x for x in read(PATHS[3])}
new={x['entity_id']:x for x in root_predictions}
ledger=read(PATHS[1]); independent_ledger=[]
for entity in {x['entity_id'] for x in ledger}:
    for field in ['primary_calibrated_probability','direct_panel_probability','risk_call','standardized_inchikey']:
        value_old,value_new=old[entity][field],new[entity][field]
        passed=abs(float(value_old)-float(value_new))<=1e-8 if field.endswith('probability') else value_old==value_new
        recorded=next(x for x in ledger if x['entity_id']==entity and x['field']==field)
        independent_ledger.append({'entity_id':entity,'field':field,'old':value_old,'new':value_new,'passed':passed})
        assert recorded['old']==value_old and recorded['current']==value_new and recorded['passed']=='True' and passed
check('13 prior structures and 52 historical comparisons verified from raw files',len({x['entity_id'] for x in ledger})==13 and len(ledger)==52 and len(independent_ledger)==52)
completion=json.loads(PATHS[2].read_text())
check('Root input hashes current and immutable',all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in completion['input_sha256'].items()))
check('Root completion counts agree',completion['structures_scored']==14 and completion['new_structures_scored']==1 and completion['reproduction_checks']==52)
nulls=read(HERE/'direct_positive_strategy_counts.csv')
check('Two pools preserve null performance metrics',all(all(x[y]=='' for y in ['BAcc','NPV','specificity','performance_CI']) for x in nulls))
check('Strict pool equal-coverage still k0',all(x['native_full_k']=='0' and x['accepted_n']=='0' for x in nulls if x['analysis_population']=='strict_human_transfectant' and 'forced' not in x['strategy']))
check('All audit inputs unchanged',all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in hashes.items()))
output={
 'status':'PASS_independent_DMH_pipeline_and_ledger_audit','checks_n':len(checks),'checks':checks,
 'source_key':source_key,'source_formula':source_formula,'frozen_standardized_key':key,
 'source_connectivity_block':source_key.split('-')[0],'model_connectivity_block':key.split('-')[0],
 'source_or_model_development_overlap':overlaps,
 'independent_v01_score':round(p01,8),'independent_v03_score':round(p03,8),
 'independent_abstention_reasons':expected_reasons,'historical_ledger_checks':independent_ledger,
 'input_sha256':hashes,'protected_inputs_unchanged':True,
 'runtime':{'executable':sys.executable,'python':platform.python_version(),'numpy':np.__version__,'rdkit':rdkit.__version__,'sklearn':sklearn.__version__,'joblib':joblib.__version__},
 'warnings':sorted({f'{w.category.__name__}: {w.message}' for w in captured}),
 'review_scope':'No root runner or project preprocessing/feature/scoring helpers imported. Independent API pipeline on DMH only; 52 prior ledger comparisons verified without rerunning prior models. Numerical reproducibility is not new biological transport evidence.',
}
(HERE/'independent_expansion_verification.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':output['status'],'checks_n':len(checks),'v01':round(p01,8),'v03':round(p03,8),'key':key},ensure_ascii=False),flush=True)
