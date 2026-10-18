#!/usr/bin/env python3
"""Retrospective source-level transport checks; no model fitting or old-file writes."""
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OLD = ROOT / 'outputs/20261009_adc_pair_experiments'

def read_csv(path):
    return list(csv.DictReader(path.open(encoding='utf-8-sig')))

def write_csv(name, rows):
    with (HERE/name).open('w', encoding='utf-8-sig', newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def dump(name, data):
    (HERE/name).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')

INPUTS = [
    '数据/processed/adc_release_species_transporter_evidence_long_v6.csv',
    '数据/processed/direct_bidirectional_transport_numeric_v3.csv',
    '数据/processed/free_payload_release_pair_evidence_v3.csv',
    '数据/processed/foti_figure8_digitization_v1.csv',
    'outputs/20261009_adc_central_pairing/local_pair_audit.csv',
    'outputs/20261009_adc_pair_experiments/pairs/results_curated_v2/entities_scored_or_hold.csv',
    'outputs/20261009_adc_pair_experiments/pairs/results_curated_v2/source_vs_frozen_model_structure_audit.csv',
    'outputs/20261009_adc_pair_experiments/identity/accepted_structure.json',
    '数据/raw/second_batch_release_species/Foti_2015_Figure8_LysMCCDM1_Pgp.jpg',
    '数据/raw/second_batch_release_species/Foti_2015_PubMed_metadata.xml',
    '数据/raw/second_batch_release_species/PMDA_belantamab_mafodotin_English_report.pdf',
    '数据/raw/third_batch_pair_resolution/FDA_Elahere_2022_MultidisciplineReview.pdf',
    '数据/raw/fifth_batch_matched_pairs/PMDA_Kadcyla_2013_Nonclinical_PK.pdf',
    '数据/raw/seventh_batch_direct_transport/PMDA_Padcev_MMAE_transport_CTD_2021.pdf',
    '数据/raw/fifth_batch_matched_pairs/US7750116B1_Auristatin_Patent.pdf',
    '数据/raw/roth_2026/PMC12883344_fulltext.xml',
]

protocol={
    'created_at':datetime.now(timezone.utc).isoformat(),
    'timing':'Retrospective follow-up: prior scores and curated literature outcomes were already known. Not blind, preregistered, or prospective validation.',
    'models':'Existing frozen v0.1, v0.3, and v0.4 scores. No retraining, recalibration, selector changes, or source-label changes.',
    'score_direction':'Operational high/low is score >= 0.5 / score < 0.5; a raw-score direction is not an accepted v0.4 prediction or calibrated substrate probability.',
    'source_direction':'Retain author direction. Strong direct positive means the primary source reports high bidirectional ER plus matched inhibitor rescue; low source evidence is limited to its tested system and material.',
    'endpoints':'Direct transport, cytotoxicity/GI50 resistance, and regulatory call are separate evidence layers. No cross-layer numerical pooling.',
    'unit_of_count':'Unique molecular entity in one evidence layer. Concentrations, source runs, duplicated free/released MMAE identities and inhibitor arms are not independent samples.',
    'parent_release_rule':'No cross-study delta ER or fold change. DM1 versus Lys-MCC-DM1 remains unmatched because cells, concentrations, labs and source materials differ.',
    'figure_review':'Manual source-image coordinates with conservative +/- 2 source pixels. Previous numerical extraction was visible to this research agent, so this is not a blinded second curator or independent human curation. No gold promotion.',
    'material_mapping':'Source material name/structure is separate from frozen standardized model representation. Mixture/stereochemical and tautomerization limitations remain visible.',
    'input_sha256':{s:sha(ROOT/s) for s in INPUTS},
}
dump('protocol_and_input_sha256.json',protocol)

# Visually selected coordinates from the original 708 x 719 figure. The two bars
# are ER in parental / MDR1 LLC-PK1, not Papp A-B / B-A. We do not infer Papp.
points = [
    ('Mannitol', '', '', 169,399,195,409),
    ('Cyclosporin A', '', '', 238,343,263,165),
    ('Cyclosporin A + elacridar', '', 'elacridar', 308,419,333,438),
    ('Lys-MCC-DM1 0.1 uM','0.1','',377,402,403,401),
    ('Lys-MCC-DM1 1 uM','1','',447,404,472,403),
    ('Lys-MCC-DM1 5 uM','5','',517,408,542,417),
    ('Lys-MCC-DM1 10 uM','10','',586,411,611,422),
    ('Lys-MCC-DM1 5 uM + elacridar','5','elacridar',656,428,681,426),
]
image_path=ROOT/'数据/raw/second_batch_release_species/Foti_2015_Figure8_LysMCCDM1_Pgp.jpg'
assert Image.open(image_path).size==(708,719)
axis0,axis12=440.0,69.0
slope=12/(axis0-axis12)
manual=[]
for group,c,inh,xp,yp,xm,ym in points:
    p=slope*(axis0-yp);m=slope*(axis0-ym)
    u=2*slope
    for cell,x,y,val in [('parental',xp,yp,p),('human MDR1',xm,ym,m)]:
        manual.append({'review_id':f'FOTI_VIS_{len(manual)+1:02d}', 'group':group,'concentration_uM':c,'inhibitor':inh,'cell':cell,'x_center_source_px':x,'mean_bar_top_y_source_px':y,'axis_zero_y_source_px':axis0,'axis_12_y_source_px':axis12,'ER_visual_review':round(val,4),'ER_reading_lower':round(max(0,val-u),4),'ER_reading_upper':round(val+u,4),'reading_uncertainty':'plus/minus 2 source pixels; not biological SD or CI','Papp_AB':'not_reported','Papp_BA':'not_reported','source_sha256':sha(image_path),'source_locator':'Foti 2015 Figure 8','review_status':'research_agent_visual_review_not_blinded_or_independent_human_curator','gold_eligible':'no'})
write_csv('Foti_Figure8_visual_coordinate_review.csv',manual)
ratio=[]
for group,c,inh,xp,yp,xm,ym in points:
    p=slope*(axis0-yp);m=slope*(axis0-ym);u=2*slope
    ratio.append({'group':group,'concentration_uM':c,'inhibitor':inh,'ER_parental':round(p,4),'ER_human_MDR1':round(m,4),'net_ER_human_over_parental':round(m/p,4),'net_ER_reading_lower':round(max(0,m-u)/(p+u),4),'net_ER_reading_upper':round((m+u)/max(0.0001,p-u),4),'definition':'ER(MDR1) / ER(parental) within this source and concentration; not Papp ratio or parent-to-release delta','uncertainty':'pixel reading bounds only; no replicated biological error propagation','gold_eligible':'no'})
write_csv('Foti_within_source_net_ER_visual_review.csv',ratio)
old=read_csv(ROOT/'数据/processed/foti_figure8_digitization_v1.csv')
agreement=[]
for new,oldr in zip(manual,old):
    v=float(oldr['efflux_ratio_consensus'])
    agreement.append({'review_id':new['review_id'],'old_record_id':oldr['digitization_record_id'],'old_algorithm_ER':v,'new_visual_ER':new['ER_visual_review'],'absolute_difference':round(abs(v-new['ER_visual_review']),4),'within_new_reading_bounds':new['ER_reading_lower']<=v<=new['ER_reading_upper'],'independent_blind_second_curator':False,'gold_eligible':False})
write_csv('Foti_old_vs_visual_review.csv',agreement)

preds={r['entity_id']:r for r in read_csv(OLD/'pairs/results_curated_v2/entities_scored_or_hold.csv')}
audits={r['source_entity_id']:r for r in read_csv(OLD/'pairs/results_curated_v2/source_vs_frozen_model_structure_audit.csv')}

OBS=[]
def add(id,entity,layer,system,c,unit,inh,inhc,pa,pb,er,source,locator,primaryreview,direction,notes='',sd=''):
    OBS.append({'observation_id':id,'entity_id':entity,'entity_name':preds[entity]['entity_name'],'evidence_layer':layer,'transporter':'human ABCB1','system':system,'concentration':c,'concentration_unit':unit,'inhibitor':inh,'inhibitor_concentration':inhc,'Papp_AB_source':pa,'Papp_BA_source':pb,'ER_source_or_visual':er,'source_error':sd,'effect_other':'','effect_other_definition':'','source_direction':direction,'source_url':source,'source_locator':locator,'primary_review':primaryreview,'material_identity_status':preds[entity]['identity_status'],'model_source_inchikey':preds[entity]['source_inchikey'],'model_standardized_inchikey':preds[entity]['standardized_inchikey'],'model_identity_mapping_limit':preds[entity]['identity_limitation'],'notes':notes})
    OBS[-1]['observation_role']=('parental_control_arm' if 'parental' in system else 'inhibitor_arm' if inh and inh.lower()!='none' else 'transporter_assay_arm' if layer.startswith('direct') else 'source_context_or_comparison')

src=read_csv(ROOT/'数据/processed/direct_bidirectional_transport_numeric_v3.csv')
chosen=[r for r in src if r['transport_record_id'].startswith(('B5_DM1_','B5_CYSMC_','B7_PMDA_8004_')) and r['entity_name'] in ['DM1','cys-mcMMAF','MMAE']]
for r in chosen:
    entity=r['entity_id']
    loc={'DM1':'PDF pp119-120; printed CTD pp80-81, Table 2.6.5.15-1','cys-mcMMAF':'PDF p33; printed p28, Section 4.5.3','MMAE':'PDF p76; printed p39, Study 8004'}[r['entity_name']]
    add(r['transport_record_id'],entity,'direct_bidirectional_numeric',r['transporter_system'],r['concentration_value'],r['concentration_unit'],r['inhibitor'],r['inhibitor_concentration_value']+' '+r['inhibitor_concentration_unit'],r['papp_ab_mean'],r['papp_ba_mean'],r['efflux_ratio_mean'],r['source_url'],loc,'text plus source PDF pixels reviewed','source_study_positive_transport_context_not_individual_arm_label',r['notes'],r['efflux_ratio_sd'])

for r in manual:
    if r['group'].startswith('Lys-MCC'):
        add(r['review_id'],'ADC_REL_001','direct_transport_figure_exploratory','LLC-PK1 '+r['cell'],r['concentration_uM'],'uM',r['inhibitor'],'not available','','',r['ER_visual_review'],'https://doi.org/10.1124/dmd.115.064576','Figure 8; original 708x719 JPEG','source JPEG pixels reviewed; full text live access blocked','low_transport_signal_under_tested_conditions','No reciprocal Papp values; epimer mixture/material limitation; not gold. Pixel error bounds in coordinate CSV.')

for entity in ['ADC_PAYLOAD_009','ADC_REL_004']:
    add('FDA_CALL_'+entity,entity,'regulatory_call_only','not reproduced in review','','','','','','','','https://www.accessdata.fda.gov/drugsatfda_docs/nda/2022/761310Orig1s000MultidisciplineR.pdf','PDF p77; printed p76, Transporter Substrate row','text plus source PDF pixels reviewed','author_substrate_call_only','No ER, concentration or inhibitor arm reproduced. Correct S-methyl thioether identity from previous independently audited registry source; old S2 graph excluded.')

# Preserve Roth's actual endpoint: table footnote says GI50, not IC50 or ER.
xml=ET.parse(ROOT/'数据/raw/roth_2026/PMC12883344_fulltext.xml')
tab=next(t for t in xml.findall('.//table-wrap') if t.get('id')=='t3')
functional=[]
for tr in tab.findall('./table/tbody/tr'):
    cells=[''.join(c.itertext()).strip() for c in tr]
    nm=cells[0].split(' (')[0]
    mapping={'DM1':'ADC_PAYLOAD_008','DM4':'ADC_PAYLOAD_009','MMAE':'ADC_PAYLOAD_017','MMAF':'ADC_PAYLOAD_018'}
    if nm not in mapping:continue
    entity=mapping[nm]
    add('ROTH_T3_'+nm,entity,'functional_GI50_resistance','HEK293 MDR19 vs empty-vector pcDNA','','','','','','','','https://pmc.ncbi.nlm.nih.gov/articles/PMC12883344/','Table 3 and footnote','primary JATS table and footnote reviewed','functional_high_cross_resistance' if nm!='MMAF' else 'borderline_functional_only_no_direct_substrate_label','Three independent experiments; means +/- SD; RR is GI50(MDR19)/GI50(pcDNA). This table is not bidirectional transport.')
    OBS[-1]['effect_other']=cells[3]
    OBS[-1]['effect_other_definition']='GI50(MDR19) / GI50(pcDNA)'
    OBS[-1]['source_error']=f'pcDNA {cells[1]}; MDR19 {cells[2]}; units from {cells[0]}'
    functional.append({'entity':nm,'pcDNA_GI50_mean_SD':cells[1],'MDR19_GI50_mean_SD':cells[2],'RR_author':cells[3],'endpoint':'GI50; not IC50 or ER','source_table':'Roth 2026 Table 3','primary_sha256':sha(ROOT/'数据/raw/roth_2026/PMC12883344_fulltext.xml')})
write_csv('Roth_Table3_functional_recheck.csv',functional)
add('PATENT_T7_MMAF','ADC_PAYLOAD_018','functional_inhibitor_rescue','NCI/ADR-RES cytotoxicity resazurin assay','','','verapamil','10 uM table / 40 uM methods conflict','','','','https://patents.google.com/patent/US7750116B1/en','PDF p37; patent printed columns61-62, Table7','text plus source PDF pixels reviewed','poor_substrate_author_call_with_functional_rescue','Table free MMAF0.45 vs0.07; units source text mM. Preserve units/conflict rather than silently converting. Not a transport endpoint or negative substrate gold.')
OBS[-1]['effect_other']='6.4';OBS[-1]['effect_other_definition']='IC50(-verapamil) / IC50(+verapamil)';OBS[-1]['source_error']='not reported'
write_csv('primary_observation_evidence_rows.csv',OBS)

summ=[]
targets=[('ADC_PAYLOAD_008','direct_positive','DM1 high ER and matched inhibitor rescue'),('ADC_REL_001','direct_low_exploratory','Foti low ER; exact material mixture/low permeability limit'),('ADC_REL_006','direct_positive','cys-mcMMAF ER34 to1.1'),('ADC_PAYLOAD_017','direct_positive','MMAE high ER and matched inhibitor rescue'),('ADC_PAYLOAD_009','call_only_positive','DM4 regulatory source call'),('ADC_REL_004','call_only_positive','S-methyl-DM4 regulatory source call'),('ADC_PAYLOAD_018','functional_ambiguous','MMAF weak/borderline functional RR and patent rescue; no direct ER')]
for entity,direction,description in targets:
    p=preds[entity];a=audits[entity]
    b01=float(p['v0_1_risk_score'])>=0.5;b03=float(p['v0_3_risk_score'])>=0.5
    expected={'direct_positive':True,'direct_low_exploratory':False,'call_only_positive':True}.get(direction)
    summ.append({'entity_id':entity,'entity_name':p['entity_name'],'source_evidence_category':direction,'source_summary':description,'v0_1_raw_score':p['v0_1_risk_score'],'v0_1_operational_direction':'high' if b01 else 'low','v0_3_raw_score':p['v0_3_risk_score'],'v0_3_operational_direction':'high' if b03 else 'low','v0_1_vs_source':'unevaluable' if expected is None else ('aligned' if b01==expected else 'discordant'),'v0_3_vs_source':'unevaluable' if expected is None else ('aligned' if b03==expected else 'discordant'),'v0_4_accepted':p['accepted_prediction'],'v0_4_call':p['risk_call'],'v0_4_reasons':p['abstention_reasons_json'],'v0_4_vs_source':'abstained_not_evaluated' if p['accepted_prediction']!='yes' else ('aligned' if expected==b01 else 'discordant'),'source_identity_limitation':p['identity_limitation'],'source_inchikey':p['source_inchikey'],'model_standardized_inchikey':p['standardized_inchikey'],'frozen_preprocessing_connectivity_changed':a['connectivity_block_changed_by_frozen_standardizer'],'full_stereochemical_identity_validated':'no','blind_or_prospective_validation':False,'same_system_parent_release_delta_available':False})
write_csv('unique_entity_model_direction_crosscheck.csv',summ)

status={
    'completed':'retrospective primary-observation/source-direction crosschecks',
    'primary_observation_rows':len(OBS),
    'unique_entities_with_direction_or_context':len(summ),
    'direct_entity_comparisons_including_exploratory_material_class':4,
    'direct_table_numeric_positive_entities':3,
    'direct_low_exploratory_entities':1,
    'raw_v0_1_and_v0_3_direct_directions_aligned':3,
    'raw_v0_1_and_v0_3_direct_directions_discordant':1,
    'discordant_entity':'Lys-MCC-DM1; high raw scores versus low Foti Figure8 transport under tested conditions',
    'v0_4_direct_entities_accepted':1,
    'v0_4_direct_entities_abstained':3,
    'v0_4_accepted_direct_entity':'MMAE',
    'call_only_positive_entities':2,
    'call_only_v0_4_abstained':2,
    'functional_direction_held_entities':1,
    'direct_same_system_parent_release_pairs':0,
    'new_wet_lab_experiments':0,
    'Foti_manual_points':len(manual),
    'old_algorithm_within_visual_bounds':sum(r['within_new_reading_bounds'] for r in agreement),
    'Foti_gold_promoted':False,
    'no_accuracy_estimate':'Small, heterogeneous, purposively selected retrospective evidence is not an independent ADC benchmark. Neither 3/4 alignment nor 1/1 accepted alignment is a validated accuracy estimate.',
    'important_limits':['Foti original full text unavailable in live follow-up; stored primary figure re-read and DOI/title metadata checked.','No reciprocal Foti Papp values exist in this source figure; do not invent them.','Lys-MCC mixture and low permeability limit transporter interpretation.','Frozen preprocessing removes some stereochemistry and changes tautomer connectivity blocks in cys-mcMMAF and Lys-MCC-DM1.','Sources are independent observations but had already been curated before this follow-up; not new blind source-heldout validation.'],
    'missing_or_held':['Lys-PEG4Mal-DM1 exact machine identity','S-methyl-DM1 stereochemistry','isolated Lys-SPDB-DM4 transport','isolated cys-VC-MMAE transport','same-system parent versus released transport for all changed pairs'],
    'protected_inputs_unchanged':all(sha(ROOT/s)==v for s,v in protocol['input_sha256'].items()),
}
dump('validation_status.json',status)

fig,ax=plt.subplots(figsize=(8.2,4.2))
lmd=ratio[3:]
x=list(range(len(lmd)))
p=[r['ER_parental'] for r in lmd];m=[r['ER_human_MDR1'] for r in lmd]
ax.errorbar([i-.08 for i in x],p,yerr=[2*slope]*len(x),fmt='o',color='#777777',label='Parental LLC-PK1')
ax.errorbar([i+.08 for i in x],m,yerr=[2*slope]*len(x),fmt='s',color='#1b4f72',label='Human MDR1 LLC-PK1')
ax.axhline(2,color='black',ls=':',lw=1)
ax.set_xticks(x,['0.1 uM','1 uM','5 uM','10 uM','5 uM + elacridar'])
ax.set_ylim(0,2.25);ax.set_ylabel('ER = Papp B-A / Papp A-B')
ax.set_title('Lys-MCC-DM1: retrospective visual review of Foti Figure 8')
ax.legend(frameon=False,ncol=2)
fig.text(.5,.01,'Error bars are +/- 2 source pixels, not biological SD. Figure reading remains exploratory.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.06,1,1));fig.savefig(HERE/'figures/Foti_visual_ER_review.png',dpi=200);fig.savefig(HERE/'figures/Foti_visual_ER_review.pdf');plt.close(fig)

assert status['protected_inputs_unchanged']
print(json.dumps(status,ensure_ascii=False,indent=2))
