#!/usr/bin/env python3
"""Read-only input audit; writes only into protocol_review. No score metrics."""
import csv
import json
import hashlib
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import Fragments

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BIOGEN = ROOT / 'outputs/20260810_soskuti_identity_double_review_v1/reanalysis/external/analysis/abcb1_quantity_biogen_analysis_rows_identity_audited_v2.csv'
ASAP = ROOT / '数据/processed/abcb1_selective_v0_4_asap2025_posthoc_predictions_v1.csv'
ADC = ROOT / 'outputs/20261012_adc_completion/scoring/species_predictions.csv'
ADC_FAMILY = ROOT / 'outputs/20261012_adc_completion/scored_family_inventory.csv'
AMIDE_QUERY = Chem.MolFromSmarts('[CX3](=[OX1])[NX3]')
REGIONS = {
    'MW_ge_600': lambda d: d['molecular_weight'] >= 600,
    'TPSA_ge_140': lambda d: d['tpsa'] >= 140,
    'amide_ge_4_and_MW_ge_500': lambda d: d['amide_count_frozen'] >= 4 and d['molecular_weight'] >= 500,
}

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

def descriptor(row):
    d = json.loads(row['physicochemical_descriptors_json'])
    mol = Chem.MolFromSmiles(row['standardized_isomeric_smiles'])
    if mol is None:
        raise ValueError('unparseable frozen standardized input')
    d['amide_count_frozen'] = len(mol.GetSubstructMatches(AMIDE_QUERY))
    d['amide_count_frozen_alternative_fr_amide'] = int(Fragments.fr_amide(mol))
    return d

def rows_biogen():
    rows = [r for r in read(BIOGEN) if r['identity_clean'] == 'yes' and r['model_usable'] == 'yes']
    assert len(rows) == 2624
    out = []
    for r in rows:
        assert r['exact_connectivity_overlap_primary_443'] == 'no'
        assert r['exact_connectivity_overlap_direct_84'] == 'no'
        out.append(dict(uid=r['external_row_id'], key=r['standardized_inchikey'],
                        cluster=r['bootstrap_cluster'], y=int(r['frozen_er2_binary_label']),
                        full=r['accepted_prediction'] == 'yes',
                        p=float(r['primary_calibrated_probability']),
                        q=float(r['direct_panel_probability']), d=descriptor(r)))
    assert len({r['uid'] for r in out}) == 2624
    return out

def rows_asap():
    groups = defaultdict(list)
    for r in read(ASAP):
        groups[r['standardized_inchikey']].append(r)
    assert len(groups) == 56
    out = []
    for key, rs in sorted(groups.items()):
        assert len({r['physicochemical_descriptors_json'] for r in rs}) == 1
        assert len({r['standardized_isomeric_smiles'] for r in rs}) == 1
        assert len({r['accepted_prediction'] for r in rs}) == 1
        assert len({r['primary_calibrated_probability'] for r in rs}) == 1
        logs = sorted(float(r['log2_efflux_ratio']) for r in rs)
        n = len(logs)
        med = logs[n//2] if n % 2 else (logs[n//2-1]+logs[n//2])/2
        r = rs[0]
        out.append(dict(uid=key, key=key, cluster=r['murcko_scaffold'] or 'acyclic::'+key,
                        y=int(med>=1), full=r['accepted_prediction']=='yes',
                        p=float(r['primary_calibrated_probability']),
                        q=float(r['direct_panel_probability']), d=descriptor(r)))
    return out

def counts(source, region, rows, policy):
    all_n = len(rows)
    k = sum(r['full'] for r in rows)
    if policy == 'population':
        selected = rows
    elif policy == 'full':
        selected = [r for r in rows if r['full']]
    else:
        score = 'p' if policy == 'primary_margin' else 'q'
        selected = sorted(rows, key=lambda r: (-abs(r[score]-.5), r['uid']))[:k]
    return dict(source=source, region=region, policy=policy,
                population_n=all_n, selected_n=len(selected),
                positive_n=sum(r['y'] for r in selected),
                negative_n=sum(1-r['y'] for r in selected),
                cluster_n=len({r['cluster'] for r in selected}),
                coverage=len(selected)/all_n if all_n else None)

def write_csv(name, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/name).open('w',encoding='utf-8-sig',newline='') as handle:
        w = csv.DictWriter(handle,fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def main():
    HERE.mkdir(parents=True,exist_ok=True)
    before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (BIOGEN,ASAP,ADC,ADC_FAMILY)}
    output=[]
    for source, rows in [('Biogen2624',rows_biogen()),('ASAP56',rows_asap())]:
        for region, pred in REGIONS.items():
            subset=[r for r in rows if pred(r['d'])]
            for policy in ('population','full','primary_margin','direct_margin'):
                output.append(counts(source,region,subset,policy))
    mapping=[]
    families={r['entity_id']:r['family'] for r in read(ADC_FAMILY)}
    for r in read(ADC):
        d=descriptor(r)
        sm=Chem.MolFromSmiles(r['source_smiles'])
        source_amides=len(sm.GetSubstructMatches(AMIDE_QUERY)) if sm is not None else None
        source_fr_amide=int(Fragments.fr_amide(sm)) if sm is not None else None
        mapping.append(dict(entity_id=r['entity_id'],entity_name=r['entity_name'],
                            family=families[r['entity_id']],role=r['role'],
                            molecular_weight=d['molecular_weight'],tpsa=d['tpsa'],
                            amide_count_frozen=d['amide_count_frozen'],amide_count_source=source_amides,
                            amide_count_frozen_alternative_fr_amide=d['amide_count_frozen_alternative_fr_amide'],
                            amide_count_source_alternative_fr_amide=source_fr_amide,
                            accepted_prediction=r['accepted_prediction'],
                            **{key:pred(d) for key,pred in REGIONS.items()}))
    write_csv('independent_region_feasibility.csv',output)
    write_csv('ADC14_fixed_descriptor_region_mapping.csv',mapping)
    after={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (BIOGEN,ASAP,ADC,ADC_FAMILY)}
    assert before == after
    (HERE/'feasibility_audit_manifest.json').write_text(json.dumps(dict(
        input_sha256=before,inputs_unchanged=True,score_metrics_computed=False,
        rdkit_version=rdBase.rdkitVersion,
        descriptor_policy='MW/TPSA use frozen table; [CX3](=[OX1])[NX3] unique matches computed on frozen standardized model input; alternative fr_amide and source counts separate diagnostic',
        outputs=['independent_region_feasibility.csv','ADC14_fixed_descriptor_region_mapping.csv'],
        analysis_scope='Post hoc chemical stress regions on previously observed sources; not ADC validation'),ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(output,ensure_ascii=False,indent=2))
    print(json.dumps(mapping,ensure_ascii=False,indent=2))

if __name__ == '__main__':
    main()
