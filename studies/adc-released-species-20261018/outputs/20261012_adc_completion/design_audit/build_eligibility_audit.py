#!/usr/bin/env python3
"""Finite retrospective ADC eligibility audit; stdlib, no fit or source edits."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def read(rel):
    return list(csv.DictReader((ROOT / rel).open(encoding='utf-8-sig')))


def write(name, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (HERE / name).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


INPUTS = [
    'outputs/20261009_adc_pair_experiments/pairs/results_curated_v2/entities_scored_or_hold.csv',
    'outputs/20261009_adc_followup_validation/transport_evidence/primary_observation_evidence_rows.csv',
    '数据/processed/adc_direct_numeric_binary_anchors_v1.csv',
    '数据/processed/direct_bidirectional_transport_numeric_v3.csv',
    '数据/processed/adc_release_species_transporter_evidence_long_v6.csv',
    '数据/processed/dxd_maaa1181a_caco2_transport_numeric_v1.csv',
    '数据/processed/bay716_caco2_transport_numeric_v1.csv',
    '数据/processed/free_payload_release_pair_evidence_v3.csv',
    '数据/processed/abcb1_screener_v0_1_training_audit.csv',
    '数据/processed/soskuti_2024_mdr1_species_training63_identity_audited_v2.csv',
    '数据/processed/abcb1_regulatory_expansion_batch4_negative_entities_v1.csv',
    '数据/processed/independent_human_transporter_negative_contrast_panel_observations_v1.csv',
    '模型/ABCB1_screener_v0_1/manifest.json',
    '模型/ABCB1_assay_aware_candidate_v0_3_dev/manifest.json',
    '模型/ABCB1_selective_screener_v0_4_dev/manifest.json',
    'outputs/20261009_adc_followup_validation/transport_evidence/source_extracts/PMDA_Padcev_local.txt',
    'outputs/20261009_adc_followup_validation/transport_evidence/source_extracts/PMDA_Kadcyla_local.txt',
]
before = {r: hashlib.sha256((ROOT / r).read_bytes()).hexdigest() for r in INPUTS}
train01 = read(INPUTS[8])
train03 = read(INPUTS[9])
keys01 = {r['standardized_inchikey'] for r in train01}
keys03 = {r['standardized_inchikey'] for r in train03}
conn01 = {k.split('-')[0] for k in keys01 if k}
conn03 = {k.split('-')[0] for k in keys03 if k}


def leak(source, normalized):
    candidates = {k for k in (source, normalized) if k}
    return {
        'v01_standardized_exact_development_overlap': 'yes' if normalized and normalized in keys01 else 'no',
        'v03_standardized_exact_development_overlap': 'yes' if normalized and normalized in keys03 else 'no',
        'v01_source_or_model_connectivity_overlap': 'yes' if any(k.split('-')[0] in conn01 for k in candidates) else 'no',
        'v03_source_or_model_connectivity_overlap': 'yes' if any(k.split('-')[0] in conn03 for k in candidates) else 'no',
        'source_or_model_key_exact_any_development_overlap': 'yes' if any(k in keys01 or k in keys03 for k in candidates) else 'no',
        'identity_leakage_audit_scope': '443 v01 and 63 identity-audited v03 development rows; source and frozen model keys; no proof of upstream literature novelty',
    }


# Labels refer to finite tested assay context. Old training labels remain untouched.
QUAL = {
    'DM1': ('core_direct_positive', '1', '5.90', '1.10', '0.5 uM', 'PSC833 10 uM', 'SRC_B5_001', 'source_table_MDCKII_MDR1_plus_parental_and_matched_inhibitor_run', 'yes', '', 'Both arms are the independent inhibitor run; do not combine 6.65/7.47 with the rescue arm.'),
    'cys-mcMMAF': ('core_direct_positive', '1', '34', '1.1', '100 nM', 'verapamil 60 uM', 'SRC_B5_009', 'regulatory_summary_human_Pgp_MDCKII_numeric_rescue', 'yes', 'ADC_PAYLOAD_018', 'Correct source catabolite graph retained; exact experimental lot/succinimide adduct stereochemical composition unreported. Frozen tautomer enumeration capped.'),
    'MMAE': ('core_direct_positive', '1', '28', '1.3', '10 uM', 'PSC833 10 uM', 'SRC_B7_003', 'source_table_Caco2_selective_inhibitor_run_and_competence_control', 'yes', '', '28 is the matched inhibitor-run baseline; 49 belongs to the separate 100 uM concentration series. Supplementary human-MDR1 transport supports attribution.'),
    'Lys-MCC-DM1': ('exploratory_direct_low_HOLD', '', '', '', '0.1 | 1 | 5 | 10 uM', 'elacridar arm at 5 uM', 'SRC_B5_002', 'source_figure_LLC_PK1_human_Pgp_and_parental', 'no', 'ADC_PAYLOAD_008', 'Low ER signal is retained as a failure probe; low passive permeability, epimer mixture, no reciprocal Papp or assay recovery, already-visible image digitization; no negative gold.'),
    'DXd': ('secondary_multitransporter_continuous', '', '10.7;22.1', '4.37;13.8', '1 uM PMDA;10 uM OBI', 'verapamil and other transporter inhibitors', 'SRC_B6_001|SRC_B7_001', 'Caco2_multitransporter_numeric', 'no', '', 'Do not pool sources; PMDA multidrug inhibitor panel and OBI partial rescue do not create gene-isolated binary transport truth.'),
    'SN-38': ('secondary_multitransporter_continuous', '', '20.0', '10.5', '10 uM', 'verapamil', 'SRC_B7_001', 'Caco2_multitransporter_numeric', 'no', '', '47.5% reduction; preserve Pgp contribution as source context, not universal binary label.'),
    'Exatecan': ('secondary_multitransporter_continuous', '', '2.7', '1.7', '10 uM', 'verapamil', 'SRC_B7_001', 'Caco2_multitransporter_numeric', 'no', '', 'Lower than DXd/SN38 does not mean negative; 37.0% reduction and later qualitative non-substrate call do not make endpoint-matched negative gold.'),
    'DM4': ('regulatory_call_only', '', '', '', '', '', 'FDA_Elahere', 'source_substrate_call_without_transport_arms', 'no', '', 'Call-only evidence excluded from direct numeric benchmark.'),
    'S-methyl-DM4': ('regulatory_call_only', '', '', '', '', '', 'FDA_Elahere', 'source_substrate_call_without_transport_arms', 'no', 'ADC_PAYLOAD_009', 'Corrected single-sulfur structure; call-only evidence does not provide continuous transport truth.'),
    'MMAF': ('functional_only', '', '', '', '', '', 'Roth_Table3|US7750116B1', 'GI50_resistance_or_cytotoxicity_rescue', 'no', '', 'GI50 and functional rescue are not direct ER; no binary negative assignment.'),
    'SG3199': ('regulatory_call_and_functional_only', '', '', '', '', '', 'FDA_Zynlonta', 'source_call_and_functional_resistance', 'no', '', 'No direct bidirectional numeric transport arms in local evidence.'),
    'Lys-SPDB-DM4': ('structure_only_no_direct_transport', '', '', '', '', '', '', 'no_direct_transport_truth', 'no', 'ADC_PAYLOAD_009', 'Processing intermediate, not a default terminal released species.'),
    'cys-VC-MMAE': ('structure_only_no_direct_transport', '', '', '', '', '', '', 'no_direct_transport_truth', 'no', 'ADC_PAYLOAD_017', 'Cys-MC-VC-PABC-MMAE processing intermediate, not free terminal MMAE.'),
}

all_entities = read(INPUTS[0])
unique = {}
holds = []
for r in all_entities:
    if r['score_status'] == 'scored_frozen_model':
        unique.setdefault(r['standardized_inchikey'], r)
    elif not r['entity_id'].startswith('DUAL_MATERIAL'):
        holds.append(r)
rows = []
for r in list(unique.values()) + holds:
    q = QUAL.get(r['entity_name'], ('identity_HOLD', '', '', '', '', '', '', 'not_admitted_identity', 'no', '', r['identity_hold']))
    category, label, er0, eri, concentration, inhibitor, source, endpoint, eligible, parent, reason = q
    row = {
        'candidate_id': r['entity_id'], 'entity_name': r['entity_name'],
        'analysis_population': 'ADC_component_candidate',
        'source_input_smiles': r['source_input_smiles'], 'source_inchikey': r['source_inchikey'],
        'model_standardized_inchikey': r['standardized_inchikey'],
        'full_source_material_identity_verified': 'no_exact_lot_or_stereomixture_validation',
        'model_score_status': r['score_status'], 'unique_structure_id': r['unique_structure_id'],
        'evidence_category': category, 'retrospective_assay_context_label': label,
        'core_direct_benchmark_eligible': 'yes' if r['entity_name'] in ['DM1','cys-mcMMAF'] and eligible=='yes' else 'no', 'no_inhibitor_ER_same_run': er0,
        'source_context_positive_anchor_eligible': eligible,
        'strict_human_transfectant_benchmark_eligible': 'yes' if r['entity_name'] in ['DM1','cys-mcMMAF'] and eligible=='yes' else 'no',
        'matched_inhibitor_ER': eri, 'concentration': concentration, 'inhibitor': inhibitor,
        'source_experiment_cluster': source, 'endpoint_context': endpoint,
        'seen_label_status': 'already_curated_before_this_audit;posthoc_not_blind',
        'parent_proxy_entity_id': parent,
        'parent_proxy_evaluation_scope': 'structural_hypothesis_only_no_matched_transport_pair' if parent else 'none_or_identity_preserving_not_new_pair',
        'same_system_parent_release_delta_eligible': 'no',
        'v01_score': r['v0_1_risk_score'], 'v03_score': r['v0_3_risk_score'],
        'full_v04_accepted': r['accepted_prediction'], 'full_v04_call': r['risk_call'],
        'full_v04_reason': r['abstention_reasons_json'], 'qualification_reason': reason,
        **leak(r['source_inchikey'], r['standardized_inchikey']),
    }
    if eligible == 'yes':
        assert row['source_or_model_key_exact_any_development_overlap'] == 'no'
        assert row['v01_source_or_model_connectivity_overlap'] == 'no'
        assert row['v03_source_or_model_connectivity_overlap'] == 'no'
    rows.append(row)
rows.append({
    'candidate_id': 'ADC_REL_012', 'entity_name': 'BAY-716', 'analysis_population': 'ADC_component_candidate',
    'evidence_category': 'identity_HOLD_low_Caco2', 'retrospective_assay_context_label': '',
    'core_direct_benchmark_eligible': 'no', 'no_inhibitor_ER_same_run': '1.6',
    'concentration': '2 uM', 'source_experiment_cluster': 'SRC_B7_005',
    'endpoint_context': 'Caco2_endogenous_multitransporter', 'seen_label_status': 'already_curated;posthoc',
    'same_system_parent_release_delta_eligible': 'no',
    'qualification_reason': 'Machine structure unresolved; low passive permeability and no gene-specific controls. No forced identity or negative label.',
})
write('adc_component_candidate_eligibility.csv', rows)
core = [r for r in rows if r.get('source_context_positive_anchor_eligible') == 'yes']
write('adc_source_context_positive_anchors.csv', core)
strict = [r for r in core if r['strict_human_transfectant_benchmark_eligible']=='yes']
write('adc_core_direct_context_anchors.csv', strict)
write('adc_strict_human_transfectant_positive_anchors.csv', strict)

# Auxiliary negative-control stress set remains separate from ADC main population.
negative = []
for r in read(INPUTS[10]):
    lim = r['tautomer_status'] != 'Completed'
    overlap = leak(r['pubchem_inchikey'], r['standardized_inchikey'])
    admitted = r['strict_external_candidate_eligible'] == 'yes' and overlap['source_or_model_key_exact_any_development_overlap'] == 'no' and overlap['v01_source_or_model_connectivity_overlap'] == 'no' and overlap['v03_source_or_model_connectivity_overlap'] == 'no'
    negative.append({
        'candidate_id': r['batch4_entity_id'], 'entity_name': r['entity_name'],
        'analysis_population': 'non_ADC_direct_negative_stress_set',
        'source_input_smiles': r['pubchem_smiles'], 'source_inchikey': r['pubchem_inchikey'],
        'model_standardized_inchikey': r['standardized_inchikey'],
        'assay_context_label': '0', 'direct_negative_evidence_eligible': 'yes' if admitted else 'no',
        'legacy_frozen_scoring_eligible': 'yes' if admitted else 'no',
        'new_enumeration_gate_diagnostic_eligible': 'no' if lim else ('yes' if admitted else 'no'),
        'tautomer_status': r['tautomer_status'],
        'source_experiment_cluster': r['upstream_experiment_cluster'],
        'source_url': r['source_url'], 'source_locator': r['source_location'],
        'local_source_file': r['local_source_file'], 'local_source_sha256': r['local_source_sha256'],
        'net_ER_source': r['net_er'], 'positive_control': r['positive_control'],
        'qc_note': r['qc_note'], 'seen_label_status': 'already_curated_before_scoring;retrospective_not_blind',
        'claim_boundary': 'Separate non-ADC negative stress panel; cannot establish ADC release-species specificity or be pooled with ADC-only positives for headline BAcc.',
        **overlap,
    })
write('non_adc_direct_negative_stress_eligibility.csv', negative)

metrics = []
for pool, poolrows in [('strict_human_transfectant',strict),('source_context_positive_anchors',core)]:
    for mode in ['full_native', 'v01_margin_equal_k', 'v03_margin_equal_k', 'v01_forced_all', 'v03_forced_all']:
        k = sum(r['full_v04_accepted'] == 'yes' for r in poolrows)
        if mode == 'full_native':
            selected = [r for r in poolrows if r['full_v04_accepted'] == 'yes']
        elif mode.endswith('equal_k'):
            score = 'v01_score' if mode.startswith('v01') else 'v03_score'
            selected = sorted(poolrows, key=lambda r: (-abs(float(r[score]) - .5), r['candidate_id']))[:k]
        else:
            selected = poolrows[:]
        metrics.append({
            'analysis_population': pool, 'strategy': mode, 'eligible_n': len(poolrows), 'eligible_positive_n': len(poolrows), 'eligible_negative_n': 0,
            'accepted_n': len(selected), 'accepted_positive_n': len(selected), 'accepted_negative_n': 0,
            'selection_ids': '|'.join(r['candidate_id'] for r in selected),
            'BAcc': '', 'specificity': '', 'NPV': '', 'accuracy_claim': '',
            'metric_status': 'not_estimable_both_classes_missing;purposefully_selected_known_positive_cases',
            'comparison_boundary': 'Forced-all coverage differs from selective k; strict pool k=0; source-context pool k=1.',
        })
write('core_comparison_feasibility_and_null_metrics.csv', metrics)
summary = {
    'created_on_user_date': '2026-10-12', 'status': 'eligibility_audit_complete_not_new_biological_validation',
    'ADC_candidate_entities_n': len(rows), 'ADC_unique_scored_structures_n': len(unique),
    'ADC_source_context_positive_anchors_n': len(core), 'ADC_source_context_negative_n': 0,
    'ADC_source_context_full_accepted_n': sum(r['full_v04_accepted'] == 'yes' for r in core),
    'ADC_strict_human_transfectant_n': len(strict), 'ADC_strict_human_transfectant_negative_n': 0,
    'ADC_strict_human_transfectant_full_accepted_n': sum(r['full_v04_accepted'] == 'yes' for r in strict),
    'same_system_parent_release_transport_pairs': 0,
    'core_candidate_categories': dict(Counter(r['evidence_category'] for r in rows)),
    'non_ADC_negative_candidates_n': len(negative),
    'non_ADC_direct_negative_legacy_scoring_eligible_n': sum(r['legacy_frozen_scoring_eligible'] == 'yes' for r in negative),
    'non_ADC_direct_negative_new_gate_diagnostic_eligible_n': sum(r['new_enumeration_gate_diagnostic_eligible'] == 'yes' for r in negative),
    'legacy_anchor_correction': 'MMAE same-run baseline 28 replaces invalid 49 concentration-series pairing only in new artifact; no old file edits.',
    'input_sha256': before,
    'protected_inputs_unchanged': all(hashlib.sha256((ROOT/r).read_bytes()).hexdigest() == h for r, h in before.items()),
}
assert summary['protected_inputs_unchanged']
(HERE/'audit_manifest.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k != 'input_sha256'},ensure_ascii=False,indent=2))
