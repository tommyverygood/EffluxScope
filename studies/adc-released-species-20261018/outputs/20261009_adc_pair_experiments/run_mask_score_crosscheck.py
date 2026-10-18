#!/usr/bin/env python3
"""Exploratory crossed mask x classifier diagnostic after the main benchmark."""
from pathlib import Path
import csv
import importlib.util
import itertools
import json
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "mask_score_diagnostic"
OUT.mkdir(exist_ok=True)
spec = importlib.util.spec_from_file_location("benchmark", HERE / "benchmark/run_benchmark.py")
b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b)

def write_csv(name, rows):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def main():
    source = HERE / "benchmark/selection_ledger.csv"
    protocol = {
        "role": "exploratory_posthoc_diagnostic_not_confirmatory",
        "trigger": "Primary matched-coverage experiment found auxiliary own-margin superior to full",
        "question": "Does the gap remain when classifier is held fixed and only mask changes?",
        "source_sha256": b.digest(source),
        "all_masks": "preexisting frozen 453-member masks from the benchmark ledger",
        "all_classifiers": "v0.1 and v0.3, each at fixed 0.5 threshold",
        "bootstrap": "2000 shared whole-cluster draws; same seed as main experiment",
        "refitting": False,
        "new_mask_tuning": False,
        "cross_population_version_join": False
    }
    (OUT / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    rows = sorted([r for r in b.read_rows(source) if r["version"] == "identity_audited_v2"],
                  key=lambda r: r["external_row_id"])
    assert len(rows) == 2624
    y = np.array([int(r["observed_er2_label"]) for r in rows])
    scores = {"primary_v01": np.array([float(r["primary_score"]) for r in rows]),
              "auxiliary_v03": np.array([float(r["auxiliary_score"]) for r in rows])}
    policies = ["v01_margin", "v03_margin", "full_v04_pool_primary_margin", "no_local_pool_primary_margin"]
    masks = {p: np.array([r[f"selected__{p}__fixed_full_coverage"] == "True" for r in rows]) for p in policies}
    assert all(m.sum() == 453 for m in masks.values())
    clusters = sorted({r["bootstrap_cluster"] for r in rows})
    lookup = {c: i for i, c in enumerate(clusters)}
    group = np.array([lookup[r["bootstrap_cluster"]] for r in rows])
    k = len(clusters)
    weights = np.random.default_rng(b.SEED).multinomial(k, np.full(k, 1 / k), size=b.N_BOOTSTRAP).astype(float)
    population = np.zeros((k, 2))
    for j, label in enumerate([1, 0]):
        np.add.at(population[:, j], group[y == label], 1)
    results, boots, point_by = [], {}, {}
    for policy, selected in masks.items():
        for channel, score in scores.items():
            pred = score >= .5
            grouped = np.zeros((k, 4))
            for j, mask in enumerate([selected & (y == 1) & pred, selected & (y == 0) & ~pred,
                                      selected & (y == 0) & pred, selected & (y == 1) & ~pred]):
                np.add.at(grouped[:, j], group[mask], 1)
            point = b.metrics(grouped.sum(axis=0)[None, :], population.sum(axis=0)[None, :])
            boot = b.metrics(weights @ grouped, weights @ population)
            name = policy + "::" + channel
            result = {"mask": policy, "classifier": channel, "selected_n": 453}
            for metric in b.METRICS:
                result[metric] = float(point[metric][0])
                for field, value in b.interval(boot[metric]).items():
                    result[metric + "_" + field] = value
            results.append(result)
            boots[name], point_by[name] = boot, result
    contrasts = []
    for lhs, rhs in itertools.combinations(boots, 2):
        same_mask = lhs.split("::")[0] == rhs.split("::")[0]
        same_channel = lhs.split("::")[1] == rhs.split("::")[1]
        if not (same_mask or same_channel):
            continue
        for metric in ["balanced_accuracy", "npv", "specificity", "sensitivity"]:
            contrasts.append({"lhs": lhs, "rhs": rhs,
                              "held_fixed": "mask" if same_mask else "classifier",
                              "metric": metric,
                              "point_delta": point_by[lhs][metric] - point_by[rhs][metric],
                              **b.interval(boots[lhs][metric] - boots[rhs][metric])})
    write_csv("crossed_metrics.csv", results)
    write_csv("paired_contrasts.csv", contrasts)
    report = ["# 选择名单与分类通道的交叉诊断", "",
              "这是看到主比较负结果后追加的探索性诊断。固定已有的453名额名单，分别用v0.1和v0.3在0.5阈值分类；不重新拟合，也不按标签选择名单。",
              "", "|453个分子的选择名单|v0.1 BAcc|v0.3 BAcc|", "|---|---:|---:|"]
    for policy in policies:
        report.append(f"|{policy}|{point_by[policy+'::primary_v01']['balanced_accuracy']:.3f}|{point_by[policy+'::auxiliary_v03']['balanced_accuracy']:.3f}|")
    (OUT / "交叉诊断.md").write_text("\n".join(report) + "\n")
    print(json.dumps({"rows": len(results), "paired_contrasts": len(contrasts),
                      "metrics": [{k: r[k] for k in ['mask', 'classifier', 'balanced_accuracy']} for r in results]}))

if __name__ == "__main__":
    main()
