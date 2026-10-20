#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reproduce the retrospective identity-audited Roth association without a model."""
from pathlib import Path
import csv
import hashlib
import json
import platform
import scipy
from scipy.stats import spearmanr


def main():
    source = Path(__file__).resolve().with_name("roth_identity_safe_12.csv")
    with source.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 12
    assert len({row["material_name"] for row in rows}) == 12
    assert all("included_after_identity_audit" in row["identity_include_flag"] for row in rows)
    assert "Deruxtecan" not in {row["material_name"] for row in rows}
    rr = [float(row["abcb1_GI50_resistance_ratio"]) for row in rows]
    scores = [float(row["in_house_Morgan_baseline_score"]) for row in rows]
    result = spearmanr(scores, rr, alternative="two-sided")
    rho, pvalue = float(result.statistic), float(result.pvalue)
    expected_rho = -0.5709290716957522
    expected_p = 0.05252355548220434
    assert abs(rho - expected_rho) < 1e-12
    assert abs(pvalue - expected_p) < 1e-12
    print(json.dumps({
        "status": "PASS",
        "n": len(rows),
        "spearman_rho": rho,
        "two_sided_asymptotic_p_value": pvalue,
        "csv_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "python_version": platform.python_version(),
        "scipy_version": scipy.__version__,
        "method": "scipy.stats.spearmanr; average ranks for ties; two-sided asymptotic p value",
        "exclusion": "Unresolved source-specific Deruxtecan identity; no score reassignment",
        "endpoint": "GI50 resistance ratio, not direct monolayer transport ER",
        "model_rerun": False,
        "new_model_fit": False,
        "ADC_classification_validation": False,
        "scope": "Numerical reproduction of an existing retrospective diagnostic; no bootstrap interval or new biological validation"
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
