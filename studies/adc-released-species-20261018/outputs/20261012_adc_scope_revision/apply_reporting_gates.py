#!/usr/bin/env python3
"""Apply the frozen reporting protocol without computing underqualified metrics."""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
protocol_path = HERE / "protocol.json"
counts_path = HERE / "feasibility/count_feasibility.json"
protocol = json.loads(protocol_path.read_text())
data = json.loads(counts_path.read_text())
population = protocol["population_gate"]
selected = protocol["matched_count_gate"]
rows = []
for row in data["counts"]:
    if row["region"] == "all_identity_clean":
        continue
    checks = {
        "population_n": (row["unique_compound_n"], population["minimum_n"]),
        "population_positive_n": (row["positive_n"], population["minimum_each_measured_class"]),
        "population_negative_n": (row["negative_n"], population["minimum_each_measured_class"]),
        "population_scaffold_n": (row["scaffold_n"], population["minimum_scaffolds"]),
        "full_selected_n": (row["full_accepted_n"], selected["minimum_k"]),
        "full_selected_positive_n": (row["full_accepted_positive_n"], selected["minimum_each_selected_measured_class"]),
        "full_selected_negative_n": (row["full_accepted_negative_n"], selected["minimum_each_selected_measured_class"]),
        "full_selected_scaffold_n": (row["full_accepted_scaffold_n"], selected["minimum_selected_scaffolds"]),
    }
    failures = {key: {"observed": observed, "required": required}
                for key, (observed, required) in checks.items() if observed < required}
    n = row["unique_compound_n"]
    rows.append({"source": row["source"], "region": row["region"],
                 "population_n": n, "positive_n": row["positive_n"], "negative_n": row["negative_n"],
                 "scaffold_n": row["scaffold_n"], "full_selected_n": row["full_accepted_n"],
                 "full_selected_positive_n": row["full_accepted_positive_n"],
                 "full_selected_negative_n": row["full_accepted_negative_n"],
                 "full_selected_scaffold_n": row["full_accepted_scaffold_n"],
                 "full_coverage": row["full_accepted_n"] / n if n else None,
                 "adc14_region_n": row["adc14_region_n"],
                 "regional_comparison_status": "not_evaluable" if failures else "eligible",
                 "failed_report_gates": failures, "BAcc": None, "MCC": None,
                 "paired_delta_BAcc": None, "paired_CI95": None,
                 "bootstrap_executed": False,
                 "null_reason": "prespecified small-denominator safeguards failed" if failures else "not_run_yet"})
assert len(rows) == 6
assert all(row["regional_comparison_status"] == "not_evaluable" for row in rows)
assert all(row["BAcc"] is None and row["paired_CI95"] is None for row in rows)
result = {"status": "PASS_frozen_gate_application", "qualifying_regions": 0, "source_region_rows": rows,
          "no_threshold_relaxation": True, "no_regional_performance_or_bootstrap_computed": True,
          "no_new_ADC_accuracy_validation": True,
          "input_sha256": {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in [protocol_path, counts_path, Path(__file__).resolve()]}}
(HERE / "regional_evaluation_status.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
with (HERE / "regional_evaluation_status.csv").open("w", encoding="utf-8-sig", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows([{**r, "failed_report_gates": json.dumps(r["failed_report_gates"], separators=(",", ":"))}
                     for r in rows])
print(json.dumps({"qualifying_regions": 0, "underqualified_source_regions": len(rows),
                  "performance_metrics": "null", "no_threshold_relaxation": True}))
