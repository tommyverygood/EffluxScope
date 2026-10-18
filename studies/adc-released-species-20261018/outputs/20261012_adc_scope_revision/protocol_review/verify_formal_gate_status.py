#!/usr/bin/env python3
"""Independently verify root gate JSON against protocol and own counts."""
import csv
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
RUN=HERE.parent
REGION={'MW_ge600':'MW_ge_600','TPSA_ge140':'TPSA_ge_140','amide_ge4_MW_ge500':'amide_ge_4_and_MW_ge_500'}

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

protocol=json.loads((RUN/'protocol.json').read_text())
status=json.loads((RUN/'regional_evaluation_status.json').read_text())
with (HERE/'independent_region_feasibility.csv').open(encoding='utf-8-sig',newline='') as handle:
    counts={(r['source'],r['region'],r['policy']):r for r in csv.DictReader(handle)}

checks=[]
for row in status['source_region_rows']:
    source,region=row['source'],REGION[row['region']]
    pop=counts[(source,region,'population')]
    full=counts[(source,region,'full')]
    expected={
        'population_n':int(pop['selected_n']),
        'positive_n':int(pop['positive_n']),
        'negative_n':int(pop['negative_n']),
        'scaffold_n':int(pop['cluster_n']),
        'full_selected_n':int(full['selected_n']),
        'full_selected_positive_n':int(full['positive_n']),
        'full_selected_negative_n':int(full['negative_n']),
        'full_selected_scaffold_n':int(full['cluster_n']),
    }
    assert all(row[k]==v for k,v in expected.items()),(source,region,'count mismatch')
    coverage=float(full['coverage']) if full['coverage'] else None
    assert row['full_coverage']==coverage,(source,region,'coverage mismatch')
    pg,mg=protocol['population_gate'],protocol['matched_count_gate']
    declared={
        'population_n':('population_n',pg['minimum_n']),
        'population_positive_n':('positive_n',pg['minimum_each_measured_class']),
        'population_negative_n':('negative_n',pg['minimum_each_measured_class']),
        'population_scaffold_n':('scaffold_n',pg['minimum_scaffolds']),
        'full_selected_n':('full_selected_n',mg['minimum_k']),
        'full_selected_positive_n':('full_selected_positive_n',mg['minimum_each_selected_measured_class']),
        'full_selected_negative_n':('full_selected_negative_n',mg['minimum_each_selected_measured_class']),
        'full_selected_scaffold_n':('full_selected_scaffold_n',mg['minimum_selected_scaffolds']),
    }
    failures={key:dict(observed=row[col],required=minimum)
              for key,(col,minimum) in declared.items() if row[col]<minimum}
    assert failures==row['failed_report_gates'],(source,region,'gate mismatch')
    assert failures and row['regional_comparison_status']=='not_evaluable'
    for metric in ('BAcc','MCC','paired_delta_BAcc','paired_CI95'):
        assert row[metric] is None,(source,region,metric,'not null')
    assert row['bootstrap_executed'] is False
    assert row['null_reason']
    checks.append(dict(source=source,region=row['region'],count_fields=len(expected),
                       gate_failure_fields=len(failures),all_null=True,no_bootstrap=True,status='PASS'))

assert len(checks)==6
assert status['qualifying_regions']==0
assert status['no_threshold_relaxation'] is True
assert status['no_regional_performance_or_bootstrap_computed'] is True
assert status['no_new_ADC_accuracy_validation'] is True
for name,sha in status['input_sha256'].items():
    assert digest(RUN/name)==sha,('stale input receipt',name)

out=dict(status='PASS',root_protocol_sha256=digest(RUN/'protocol.json'),
         gate_status_sha256=digest(RUN/'regional_evaluation_status.json'),
         independent_count_sha256=digest(HERE/'independent_region_feasibility.csv'),
         source_region_rows=6,verified_numeric_fields=48,
         verified_coverage_fields=6,null_metric_fields=24,bootstrap_executed_regions=0,
         qualifying_regions=0,checks=checks,
         scope='posthoc chemical support feasibility, not ADC performance validation')
(HERE/'formal_gate_and_status_independent_verification.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False,indent=2))
