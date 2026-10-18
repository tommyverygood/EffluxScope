#!/usr/bin/env python3
"""Independent ASAP repeat audit. Reads legacy data; writes only this directory."""
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

from rdkit import Chem, rdBase

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PROCESSED = ROOT / '数据/processed'
RAW = ROOT / '数据/raw/abcb1_asap_mdck_mdr1_v1'

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as h:
        return list(csv.DictReader(h))

def write(path, rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields); w.writeheader(); w.writerows(rows)

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def packed(values): return json.dumps(values,ensure_ascii=False,separators=(',',':'))
def unique(rows, field): return sorted({str(r[field]) for r in rows})

def main():
    predpath=PROCESSED/'abcb1_asap2025_sourceheldout_predictions_v1.csv'
    obspath=PROCESSED/'abcb1_asap_mdck_mdr1_observations_v1.csv'
    lockpath=PROCESSED/'abcb1_asap2025_sourceheldout_protocol_lock_v1.json'
    actpath=RAW/'CHEMBL5674852_activities.json'
    assaypath=RAW/'CHEMBL5674852_assays.json'
    docpath=RAW/'CHEMBL5674852_document.json'
    inputs=[predpath,obspath,lockpath,actpath,assaypath,docpath]
    before={str(p.relative_to(ROOT)):digest(p) for p in inputs}
    pred=read(predpath)
    obsrows=read(obspath)
    obs={r['observation_id']:r for r in obsrows}
    activities={str(a['activity_id']):a for a in json.loads(actpath.read_text())['activities']}
    assay=json.loads(assaypath.read_text())['assays'][0]
    lock=json.loads(lockpath.read_text())
    bykey=defaultdict(list)
    for r in pred: bykey[r['standardized_inchikey']].append(r)
    assert len(pred)==62 and len(bykey)==56
    assert len({r['observation_id'] for r in pred})==len(pred)
    field_checks=['document_chembl_id','assay_chembl_id','molecule_chembl_id',
        'parent_molecule_chembl_id','assay_group','assay_cell_type','assay_description',
        'inhibitor','source_canonical_smiles','source_units','upstream_experiment_cluster']
    repeated, summary, compounds, receipts = [], [], [], []
    source_keys=set()
    all_activity_ids=[]
    max_ratio_error=0.0
    max_log_error=0.0
    for key, rs in sorted(bykey.items()):
        os=[obs[r['observation_id']] for r in rs]
        er=[float(r['efflux_ratio']) for r in rs]
        logs=[float(r['log2_efflux_ratio']) for r in rs]
        contextual={f:unique(os,f) for f in field_checks}
        mol=Chem.MolFromSmiles(os[0]['source_canonical_smiles'])
        assert mol is not None
        source_key=Chem.MolToInchiKey(mol)
        source_keys.add(source_key)
        score01=unique(rs,'v0_1_calibrated_probability')
        score03=unique(rs,'v0_3_candidate_probability')
        scaffold=unique(rs,'murcko_scaffold')
        assert len(score01)==len(score03)==len(scaffold)==1
        assert all(len(contextual[f])==1 for f in field_checks)
        for r,o in zip(rs,os):
            for f in ['standardized_inchikey','molecule_chembl_id','papp_ab','papp_ba','efflux_ratio']:
                assert r[f]==o[f],(r['observation_id'],f)
            expected=float(o['papp_ba'])/float(o['papp_ab'])
            max_ratio_error=max(max_ratio_error,abs(float(r['efflux_ratio'])-expected))
            max_log_error=max(max_log_error,abs(float(r['log2_efflux_ratio'])-math.log2(float(r['efflux_ratio']))))
            assert int(r['weak_binary_label_er_ge_2'])==int(float(r['efflux_ratio'])>=2)
            ids=[x.strip() for x in o['endpoint_activity_ids'].split('|')]
            all_activity_ids.extend(ids)
            if len(rs)>1:
                for activity_id in ids:
                    a=activities[activity_id]
                    assert str(a['record_id'])==o['record_id']
                    assert a['molecule_chembl_id']==o['molecule_chembl_id']
                    assert a['assay_chembl_id']==o['assay_chembl_id']
                    match={'Papp A to B (mean)':'papp_ab','Papp B to A (mean)':'papp_ba',
                        'Recovery A to B (mean)':'recovery_ab_pct','Recovery B to A (mean)':'recovery_ba_pct'}
                    # Raw endpoint naming is retained; match by direction/endpoint class below.
                    typ=a['standard_type']
                    target=next((v for k,v in match.items() if k==typ),None)
                    if target is None:
                        target=('recovery_ab_pct' if 'A to B' in typ and 'Recovery' in typ else
                            'recovery_ba_pct' if 'B to A' in typ and 'Recovery' in typ else None)
                    assert target is not None,typ
                    assert float(a['standard_value'])==float(o[target])
                    receipts.append(dict(observation_id=o['observation_id'],
                        standardized_inchikey=key, raw_activity_id=activity_id,
                        raw_record_id=a['record_id'], raw_standard_type=typ,
                        raw_standard_value=a['standard_value'], raw_standard_relation=a['standard_relation'],
                        raw_standard_units=a['standard_units'], mapped_observation_field=target,
                        mapped_observation_value=o[target], raw_potential_duplicate=a.get('potential_duplicate'),
                        raw_data_validity_comment=a.get('data_validity_comment'),
                        value_and_identity_check='PASS',
                        explicit_direction_pair_sample_id='NOT_AVAILABLE_IN_CHEMBL_EXPORT'))
                repeated.append({f:o[f] for f in ['observation_id','document_chembl_id','assay_chembl_id',
                    'molecule_chembl_id','record_id','replicate_index_within_record','source_canonical_smiles',
                    'papp_ab','papp_ba','efflux_ratio','recovery_ab_pct','recovery_ba_pct',
                    'endpoint_activity_ids','source_units','upstream_experiment_cluster']}
                    | dict(standardized_inchikey=key,source_exact_inchikey_computed=source_key,
                        v0_1_calibrated_probability=r['v0_1_calibrated_probability'],
                        v0_3_candidate_probability=r['v0_3_candidate_probability'],
                        same_source_graph_in_group='yes', same_model_scores_in_group='yes',
                        context_description_consistent_in_group='yes',
                        source_replicate_kind='UNKNOWN',
                        experiment_date_or_day='NOT_REPORTED', tested_concentration='NOT_REPORTED_PER_OBSERVATION',
                        incubation_time='NOT_REPORTED_PER_OBSERVATION', material_lot='NOT_REPORTED',
                        endpoint_pairing_basis='ordinal endpoint occurrence within same source record, sorted by activity_id',
                        no_independence_claim='Source rows, not proven independent experiments or molecules'))
        median_log=statistics.median(logs)
        nrecords=len(unique(os,'record_id'))
        compound=dict(standardized_inchikey=key,
            molecule_chembl_id=contextual['molecule_chembl_id'][0],
            source_exact_inchikey_computed=source_key,
            source_stereo_key_differs_from_model_standardized_key=source_key!=key,
            source_canonical_smiles=contextual['source_canonical_smiles'][0],
            observation_count=len(rs), source_record_count=nrecords,
            observation_ids=packed([r['observation_id'] for r in rs]),
            source_record_ids=packed(unique(os,'record_id')),
            efflux_ratios=packed(er), median_log2_efflux_ratio=median_log,
            median_log2_derived_er=2**median_log,
            weak_threshold_crossing_within_compound=min(er)<2<=max(er),
            aggregated_weak_er_ge2_label=int(median_log>=1),
            v0_1_calibrated_probability=score01[0],v0_3_candidate_probability=score03[0],
            murcko_scaffold=scaffold[0], score_constancy='PASS',
            source_graph_constancy='PASS', metadata_constancy='PASS',
            material_identity_scope='same public molecule ID and source graph; actual vial/lot not reported',
            primary_unit='unique chemical entity in one vendor protocol',
            replicate_independence='not established from export')
        compounds.append(compound)
        if len(rs)>1:
            summary.append(compound | dict(er_min=min(er),er_max=max(er),
                er_range_fold=max(er)/min(er),
                same_record_has_multiple_rows=any(n>1 for n in Counter(o['record_id'] for o in os).values()),
                distinct_records_have_same_compound=nrecords>1,
                shared_context_fields=packed(field_checks),
                unspecified_context_fields=packed(['experiment_day','test_concentration','incubation_time',
                    'biological_vs_technical_repeat','material_lot','explicit_A_B_direction_pair_id']),
                history_row_weight=len(rs)/62,unique_compound_weight=1/56))
    assert len(summary)==4 and len(repeated)==10 and len(receipts)==40
    assert len(source_keys)==56
    assert len(set(all_activity_ids))==len(all_activity_ids)
    assert max_ratio_error<1e-7 and max_log_error<1e-7
    write(HERE/'unique_compound_identity_and_aggregation.csv',compounds)
    write(HERE/'repeated_group_summary.csv',summary)
    write(HERE/'repeated_observation_details.csv',repeated)
    write(HERE/'raw_activity_receipts.csv',receipts)
    matches={
        'observations_hash_matches_original_lock':digest(obspath)==lock['processed_interface']['sha256'],
        'raw_document_hash_matches_lock':digest(docpath)==lock['source_document']['document_raw_sha256'],
        'raw_assay_hash_matches_lock':digest(assaypath)==lock['source_document']['assay_raw_sha256'],
        'raw_activity_hash_matches_lock':digest(actpath)==lock['source_document']['activity_raw_sha256'],
    }
    result=dict(status='PASS_IDENTITY_REPEAT_AUDIT_WITH_CONTEXT_LIMITATIONS',
        observation_rows=62,unique_standardized_compounds=56,unique_source_exact_graph_keys=56,
        repeated_compounds=4,rows_in_repeated_compounds=10,extra_rows_over_unique_compounds=6,
        repeated_group_sizes=sorted([len(rs) for rs in bykey.values() if len(rs)>1],reverse=True),
        weak_threshold_crossing_compounds=[r['molecule_chembl_id'] for r in summary
            if r['weak_threshold_crossing_within_compound']],
        source_graph_constant_in_each_group=True,model_scores_constant_in_each_group=True,
        no_repeated_endpoint_activity_ids=True,raw_endpoint_receipts_checked=40,
        same_vendor_assay_description_in_repeated_groups=True,
        raw_assay_parameters=assay.get('assay_parameters'),
        legacy_score_graph_boundary='All grouping is audited against source SMILES; standardized key is a model representation, not procurement identity.',
        raw_pairing_limit='R1/R2 were reconstructed by ordinal endpoint occurrence within record, ordered by activity_id; export has no explicit paired sample ID.',
        material_context_limit='same chemical identity/document/assay, not proof of same lot, concentration, time, or independent experimental day.',
        statistical_unit='56 unique compounds; historical 62-row analysis preserves observation weighting and is sensitivity only.',
        prospective_status='Retrospective reanalysis of already scored/seen outcomes; old protocol prohibited independent treatment of repeated observations.',
        original_protocol_primary='Spearman on 62 estimable rows with whole Murcko-scaffold bootstrap; unique-compound aggregation is a new retrospective sensitivity/primary reporting choice.',
        inference_boundary='one MDCK-MDR1 vendor dataset; no matched parental, inhibitor rescue, strict ABCB1 or ADC-specific validation',
        max_absolute_er_reconstruction_error=max_ratio_error,max_absolute_log2_reconstruction_error=max_log_error,
        rdkit_version=rdBase.rdkitVersion,lock_checks=matches,
        inputs_sha256_before=before,inputs_sha256_after={str(p.relative_to(ROOT)):digest(p) for p in inputs})
    result['inputs_unchanged']=result['inputs_sha256_before']==result['inputs_sha256_after']
    assert all(matches.values()) and result['inputs_unchanged']
    (HERE/'audit_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','observation_rows','unique_standardized_compounds',
        'repeated_group_sizes','weak_threshold_crossing_compounds','inputs_unchanged']}))

if __name__=='__main__':main()
