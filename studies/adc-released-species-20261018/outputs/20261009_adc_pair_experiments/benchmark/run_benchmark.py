#!/usr/bin/env python3
"""Retrospective frozen-score comparison; writes only beside this script.

No models are loaded or fitted. All rules, thresholds, and target counts are
declared before the benchmark. Whole-cluster draws are shared across policies.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SEED = 20261009
N_BOOTSTRAP = 2000
LOCAL_REASONS = {
    "primary_local_neighborhood_strongly_opposite",
    "direct_local_neighborhood_strongly_opposite",
}
INPUTS = {
    "identity_audited_v2": ROOT / "outputs/20260810_soskuti_identity_double_review_v1/reanalysis/external/analysis/abcb1_quantity_biogen_analysis_rows_identity_audited_v2.csv",
    "historical_v1_sensitivity": ROOT / "数据/processed/abcb1_quantity_biogen_analysis_rows_v1.csv",
}
MODELS = {
    "v01_margin": "primary_calibrated_probability",
    "v03_margin": "direct_panel_probability",
    "full_v04_pool_primary_margin": "primary_calibrated_probability",
    "no_local_pool_primary_margin": "primary_calibrated_probability",
}
METRICS = [
    "coverage", "positive_class_coverage", "negative_class_coverage",
    "accuracy", "balanced_accuracy", "mcc", "sensitivity", "specificity",
    "ppv", "npv",
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def write_json(name: str, data) -> None:
    (HERE / name).write_text(json.dumps(clean(data), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(name: str, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Empty output {name}")
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (HERE / name).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{k: clean(v) for k, v in row.items()} for row in rows])


def divide(numerator, denominator):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.asarray(denominator) > 0, np.asarray(numerator) / np.asarray(denominator), np.nan)


def metrics(confusion: np.ndarray, population: np.ndarray) -> dict[str, np.ndarray]:
    """Columns of confusion: TP, TN, FP, FN; population: positive, negative."""
    tp, tn, fp, fn = confusion.T
    positive, negative = population.T
    selected = tp + tn + fp + fn
    sensitivity = divide(tp, tp + fn)
    specificity = divide(tn, tn + fp)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {
        "coverage": divide(selected, positive + negative),
        "positive_class_coverage": divide(tp + fn, positive),
        "negative_class_coverage": divide(tn + fp, negative),
        "accuracy": divide(tp + tn, selected),
        "balanced_accuracy": (sensitivity + specificity) / 2,
        "mcc": divide(tp * tn - fp * fn, denominator),
        "sensitivity": sensitivity,
        "specificity": specificity,
        "ppv": divide(tp, tp + fp),
        "npv": divide(tn, tn + fn),
    }


def interval(values: np.ndarray) -> dict:
    finite = values[np.isfinite(values)]
    return {
        "ci95_low": float(np.quantile(finite, 0.025)) if len(finite) else np.nan,
        "ci95_high": float(np.quantile(finite, 0.975)) if len(finite) else np.nan,
        "valid_bootstrap_replicates": int(len(finite)),
    }


def target_counts(n: int, full_n: int) -> list[tuple[str, int, float]]:
    return [
        ("coverage_05pct", int(math.floor(n * .05 + .5)), .05),
        ("coverage_10pct", int(math.floor(n * .10 + .5)), .10),
        ("fixed_full_coverage", full_n, full_n / n),
        ("coverage_17_3pct", int(math.floor(n * .173 + .5)), .173),
        ("coverage_20pct", int(math.floor(n * .20 + .5)), .20),
        ("coverage_25pct", int(math.floor(n * .25 + .5)), .25),
    ]


def evaluate(version: str, all_rows: list[dict[str, str]]):
    rows = sorted(
        [r for r in all_rows if r["identity_clean"] == "yes" and r["model_usable"] == "yes"],
        key=lambda r: r["external_row_id"],
    )
    n = len(rows)
    uids = [r["external_row_id"] for r in rows]
    assert n == 2624 and len(set(uids)) == n
    assert len({r["connectivity_block"] for r in rows}) == n
    assert all(r["exact_connectivity_overlap_primary_443"] == "no" for r in rows)
    assert all(r["exact_connectivity_overlap_direct_84"] == "no" for r in rows)
    y = np.array([int(r["frozen_er2_binary_label"]) for r in rows])
    p = np.array([float(r["primary_calibrated_probability"]) for r in rows])
    q = np.array([float(r["direct_panel_probability"]) for r in rows])
    assert np.isfinite(p).all() and np.isfinite(q).all()
    full = np.array([r["accepted_prediction"] == "yes" for r in rows])
    reasons = [json.loads(r["abstention_reasons_json"]) for r in rows]
    assert np.array_equal(full, np.array([len(z) == 0 for z in reasons]))
    no_local = np.array([not (set(z) - LOCAL_REASONS) for z in reasons])
    margin_pool = (p <= .35) | (p >= .70)
    assert int(full.sum()) == 453 and int(no_local.sum()) == 551 and int(margin_pool.sum()) == 935
    assert all((r["risk_call"] == "high_ABCB1_efflux_risk") == bool(score >= .5)
               for r, score, accepted in zip(rows, p, full) if accepted)
    # Each version is evaluated from its own complete source table. No score join.
    pools = {
        "v01_margin": np.ones(n, bool),
        "v03_margin": np.ones(n, bool),
        "full_v04_pool_primary_margin": full,
        "no_local_pool_primary_margin": no_local,
    }
    scores = {name: (q if name == "v03_margin" else p) for name in MODELS}
    specs = []
    unavailable = []
    for target_name, take, rate in target_counts(n, int(full.sum())):
        for policy in MODELS:
            candidate = np.flatnonzero(pools[policy])
            if take > len(candidate):
                unavailable.append({"version": version, "target": target_name, "policy": policy,
                                    "requested_n": take, "requested_coverage": rate,
                                    "available_pool_n": int(len(candidate)),
                                    "reason": "target_exceeds_frozen_accepted_pool; abstentions_not_force_filled"})
                continue
            ordered = sorted(candidate, key=lambda i: (-abs(scores[policy][i] - .5), uids[i]))
            selected = np.zeros(n, bool)
            selected[ordered[:take]] = True
            specs.append({"id": f"{policy}__{target_name}", "policy": policy, "target": target_name,
                          "target_coverage": rate, "selected": selected, "pred": scores[policy] >= .5,
                          "analysis_role": "primary_matched_count" if target_name == "fixed_full_coverage" else "coverage_sensitivity"})
    for name, selected, score in [
        ("v01_all", np.ones(n, bool), p), ("v03_all", np.ones(n, bool), q),
        ("full_v04_fixed", full, p), ("no_local_fixed", no_local, p),
        ("primary_threshold_only", margin_pool, p),
    ]:
        specs.append({"id": name, "policy": name, "target": "native_working_point",
                      "target_coverage": float(selected.sum() / n), "selected": selected,
                      "pred": score >= .5, "analysis_role": "frozen_reference" if name != "no_local_fixed" else "posthoc_simplification"})

    cluster_ids = sorted({r["bootstrap_cluster"] for r in rows})
    assert all(cluster_ids)
    cluster_lookup = {x: i for i, x in enumerate(cluster_ids)}
    group = np.array([cluster_lookup[r["bootstrap_cluster"]] for r in rows])
    k = len(cluster_ids)
    # One whole-cluster multinomial draw bank is reused for every strategy.
    # Every drawn cluster contributes every member, with identical multiplicity.
    weights = np.random.default_rng(SEED).multinomial(k, np.full(k, 1 / k), size=N_BOOTSTRAP).astype(float)
    population = np.zeros((k, 2), float)
    for j, label in enumerate([1, 0]):
        np.add.at(population[:, j], group[y == label], 1)
    boot_population = weights @ population
    point_population = population.sum(axis=0, keepdims=True)
    results, bootstrap, ledger = [], {}, []
    for spec in specs:
        selected, pred = spec["selected"], spec["pred"]
        grouped = np.zeros((k, 4), float)
        for j, mask in enumerate([
            selected & (y == 1) & pred, selected & (y == 0) & ~pred,
            selected & (y == 0) & pred, selected & (y == 1) & ~pred,
        ]):
            np.add.at(grouped[:, j], group[mask], 1)
        counts = grouped.sum(axis=0)
        boot = metrics(weights @ grouped, boot_population)
        point = metrics(counts[None, :], point_population)
        result = {"version": version, "strategy_id": spec["id"], "policy": spec["policy"],
                  "target": spec["target"], "analysis_role": spec["analysis_role"],
                  "target_coverage": spec["target_coverage"], "population_n": n,
                  "n_bootstrap_clusters": k, "selected_n": int(selected.sum()),
                  "selected_positive_n": int(y[selected].sum()),
                  "selected_negative_n": int(selected.sum() - y[selected].sum()),
                  "predicted_high_n": int((selected & pred).sum()),
                  "predicted_low_n": int((selected & ~pred).sum()),
                  "tp": int(counts[0]), "tn": int(counts[1]), "fp": int(counts[2]), "fn": int(counts[3])}
        for name in METRICS:
            result[name] = float(point[name][0])
            for field, value in interval(boot[name]).items():
                result[f"{name}_{field}"] = value
        bootstrap[spec["id"]] = boot
        results.append(result)
    fixed = next(r for r in results if r["strategy_id"] == "full_v04_fixed")
    assert (fixed["tp"], fixed["tn"], fixed["fp"], fixed["fn"]) == (221, 119, 84, 29)
    assert abs(fixed["balanced_accuracy"] - .735103448275862) < 1e-12
    primary = next(s for s in specs if s["id"] == "full_v04_pool_primary_margin__fixed_full_coverage")
    assert np.array_equal(primary["selected"], full)
    contrasts = []
    for target in [x[0] for x in target_counts(n, int(full.sum()))]:
        comparable = [s for s in specs if s["target"] == target]
        for lhs, rhs in itertools.combinations(comparable, 2):
            # Order principal contrasts as full/no-local minus margin when possible.
            if lhs["policy"] in {"v01_margin", "v03_margin"} and rhs["policy"] not in {"v01_margin", "v03_margin"}:
                lhs, rhs = rhs, lhs
            assert lhs["selected"].sum() == rhs["selected"].sum()
            lpoint = next(z for z in results if z["strategy_id"] == lhs["id"])
            rpoint = next(z for z in results if z["strategy_id"] == rhs["id"])
            for name in METRICS:
                # A correlated difference in EACH SHARED draw, never two independent CIs.
                values = bootstrap[lhs["id"]][name] - bootstrap[rhs["id"]][name]
                finite = values[np.isfinite(values)]
                contrasts.append({"version": version, "target": target,
                                  "lhs_strategy": lhs["id"], "rhs_strategy": rhs["id"],
                                  "metric": name, "selected_n_each": int(lhs["selected"].sum()),
                                  "point_difference": lpoint[name] - rpoint[name],
                                  **interval(values),
                                  "bootstrap_fraction_difference_gt_zero": float(np.mean(finite > 0)) if len(finite) else np.nan,
                                  "shared_selected_n": int((lhs["selected"] & rhs["selected"]).sum()),
                                  "analysis_role": "primary_matched_count" if target == "fixed_full_coverage" else "coverage_sensitivity"})
    for i, row in enumerate(rows):
        ledger.append({"version": version, "external_row_id": uids[i],
                       "connectivity_block": row["connectivity_block"], "standardized_inchikey": row["standardized_inchikey"],
                       "bootstrap_cluster": row["bootstrap_cluster"], "observed_er2_label": int(y[i]),
                       "primary_score": float(p[i]), "auxiliary_score": float(q[i]),
                       "accepted_full": bool(full[i]), "accepted_no_local_posthoc": bool(no_local[i]),
                       "abstention_reasons_json": row["abstention_reasons_json"],
                       **{f"selected__{s['id']}": bool(s["selected"][i]) for s in specs}})
    return results, contrasts, unavailable, ledger, {
        "version": version, "operational_n": len(all_rows), "identity_clean_n": n,
        "n_positive": int(y.sum()), "n_negative": int(n-y.sum()), "clusters": k,
        "full_selected_n": int(full.sum()), "no_local_selected_n": int(no_local.sum()),
        "primary_margin_pool_n": int(margin_pool.sum()),
        "model_version_values": sorted({r["model_version"] for r in all_rows}),
        "selected_scores_and_labels_from_same_input_table": True,
        "old_full_point_reproduced": True,
    }


def write_report(results: list[dict], contrasts: list[dict], audits: list[dict]) -> None:
    main = "identity_audited_v2"
    points = [r for r in results if r["version"] == main and (r["target"] == "fixed_full_coverage" or r["strategy_id"] in {"v01_all", "v03_all", "no_local_fixed"})]
    lines = ["# 冻结 score 与选择性拒判：同覆盖率回顾性比较", "",
             "这是 Biogen/Fang 小分子数据上的回顾性计算实验。没有训练模型、改变类别阈值或新增湿实验；不能称 ADC 配对验证、前瞻验证或新来源验证。", "",
             "## 共同分析人群与时间边界", "",
             "主分析使用 identity-audited v2 的整张 analysis table；原始 v1 整表独立复算作为敏感性。各自保留 identity-clean、model-usable 的 2624 个唯一 connectivity 分子，删除与 weak443/direct84 重叠者。标签固定为 measured ER ≥ 2。两个 score 与接受理由始终来自同版本同一行，未跨版本拼接。", "",
             "v2 是既有 identity correction 后的事后重分析，不是未观察结局的确认集。原 v1 与 v2 的 scores、标签、结构、cluster 和接受决策逐行相同；因此复算一致是 identity correction 不改变该实验的敏感性证据，不能视为两次独立验证。", "",
             "## 同 coverage 工作点", "",
             "固定 full v0.4 接受 453/2624 = 17.2637%。各 comparator 恰好接受 453 个，以各自 score 到 0.5 的绝对距离降序、UID 升序打破 ties。所有类别均按 score ≥ 0.5 判为高风险。full 和 no-local 的曲线仅在其 frozen accepted pool 内取前 k；pool 不足时标 unavailable，不强制填入拒判分子。literal 17.3% 目标四舍五入为454，与固定453工作点分别报告。", "",
             "|策略|接受n|coverage|BAcc [95% CI]|MCC|敏感度|特异度|PPV|NPV|正类coverage|负类coverage|", "|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in points:
        lines.append(f"|{r['policy']} ({r['target']})|{r['selected_n']}|{r['coverage']:.3%}|{r['balanced_accuracy']:.3f} [{r['balanced_accuracy_ci95_low']:.3f}, {r['balanced_accuracy_ci95_high']:.3f}]|{r['mcc']:.3f}|{r['sensitivity']:.3f}|{r['specificity']:.3f}|{r['ppv']:.3f}|{r['npv']:.3f}|{r['positive_class_coverage']:.3%}|{r['negative_class_coverage']:.3%}|")
    lines += ["", "## 成对差值而非独立置信区间比较", "",
              "2000 次 whole-Murcko-cluster bootstrap；每次从共同 population 的2312个 cluster 有放回抽2312个完整 cluster。所有策略使用完全相同的 cluster multiplicity。每次 draw 先分别算指标，再作 lhs−rhs；百分位 CI 由这些相关的成对差值给出。pool、排名及接受名单固定，不在 resample 中重选或调阈值。因此 CI 是在当前 frozen score 与选择名单条件下的 test-population sampling uncertainty，不包括重新训练或新数据来源偏移。", "",
              "|lhs−rhs (均453个)|ΔBAcc|成对95% CI|", "|---|---:|---|"]
    for r in contrasts:
        if r["version"] == main and r["target"] == "fixed_full_coverage" and r["metric"] == "balanced_accuracy":
            lines.append(f"|{r['lhs_strategy'].split('__')[0]} − {r['rhs_strategy'].split('__')[0]}|{r['point_difference']:+.3f}|[{r['ci95_low']:+.3f}, {r['ci95_high']:+.3f}]|")
    lines += ["", "完整 no-local 固定策略接受551个（21.0%），BAcc=0.763；同coverage比较的 no-local 则是在这个551人pool内按primary margin重新只取453个，BAcc约0.746。这是两个不同接受名单。历史 full−V08 在453对551的原始工作点比较中可能显著为负，不能直接等同本实验453对453的差值；本实验同coverage的full−no-local 95% CI跨0。", "", "## 如何使用这些结果", "",
              "保留 full 相对 auxiliary 的负结果。辅助模型单独加自身 margin 已经是必须面对的实质 comparator；不能只用弱 v0.1 的提升证明复杂组合有方法优势。no-local 是在既有 external outcomes 观察后产生的简化假设，其成绩不能用于在同一数据上选出新 deployment model。", "",
              "此实验解决了同 score、同类别阈值、同 coverage 与已拒判名单的比较定义；没有解决所有模型的来源独立性、ADC released-species 外部准确率、潜在 assay-specific transport 因果确认或低风险类别的独立验证。", "",
              "## 文件", "", "- `protocol.json`：预先写入的分析规则、输入SHA256、模型版本与回顾性边界。",
              "- `input_audit.json`：原v1与v2逐行可比性及所有验收检查。",
              "- `metrics.csv`：fixed points与coverage grid，包括类别coverage、PPV/NPV及各指标CI。",
              "- `paired_contrasts.csv`：相同接受数的相关cluster-bootstrap差值。",
              "- `unattainable_targets.csv`：固定pool不能实现的coverage目标。",
              "- `selection_ledger.csv`：每个分子、score、规则理由与各策略接受名单。",
              "- `run_benchmark.py`：独立可执行复算脚本；只在本目录写文件。",
              "- `plot_benchmark.py`：从metrics生成BAcc、类别coverage、NPV三panel图与coverage曲线；没有重新选择策略或阈值。", ""]
    (HERE / "实验报告.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    HERE.mkdir(parents=True, exist_ok=True)
    protocol = {
        "analysis_id": "20261009_frozen_efflux_selective_benchmark",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "protocol_before_this_run_metrics; retrospective_outcomes_already_seen",
        "is_preregistered": False, "is_prospective": False, "is_adc_paired_validation": False,
        "inputs": {name: {"path": str(path.relative_to(ROOT)), "sha256": digest(path)} for name, path in INPUTS.items()},
        "script_sha256": digest(Path(__file__)),
        "runtime": {"python_executable": sys.executable, "python_version": sys.version, "numpy_version": np.__version__},
        "population": "same-version identity_clean=yes and model_usable=yes; unique connectivity; no overlap weak443 or direct84",
        "endpoint": "unchanged measured MDR1-MDCK ER>=2",
        "classification_threshold": .5,
        "frozen_full_selection_thresholds": {"low": .35, "high": .70},
        "ranking": "absolute distance from own 0.5 threshold, descending; external_row_id ascending ties",
        "grid": [.05, .10, "fixed_full453/2624", .173, .20, .25],
        "rounding": "floor(population_n*target+0.5); fixed full uses exact existing selected count",
        "full_and_no_local_grid": "rank primary margin only inside corresponding accepted pool; above pool targets unavailable",
        "no_local_definition": "only remove the two named local-neighborhood-opposition reasons; all other frozen reasons retained",
        "removed_reasons": sorted(LOCAL_REASONS),
        "simplification_role": "posthoc diagnostic, cannot select/redeploy a best rule on this source",
        "bootstrap": {"replicates": N_BOOTSTRAP, "seed": SEED, "unit": "whole existing Murcko scaffold cluster with existing acyclic fallback",
                      "paired": "identical cluster multiplicity bank for all masks in one common population",
                      "selection_fixed": True, "CI": "percentile95 of per-draw correlated metric differences",
                      "uncertainty_scope": "test-population sampling conditional on frozen scores and fixed masks; no refit uncertainty"},
        "prohibitions": ["no training", "no external-outcome threshold optimization", "no pooled v1/v2 score join", "no independent-CI subtraction", "no filling full-rule abstentions to reach a larger coverage", "no ADC accuracy claim"],
    }
    write_json("protocol.json", protocol)
    inputs = {name: read_rows(path) for name, path in INPUTS.items()}
    by_uid = {name: {r["external_row_id"]: r for r in rows} for name, rows in inputs.items()}
    names = list(INPUTS)
    assert all(len(rows) == len(by_uid[name]) for name, rows in inputs.items())
    assert set(by_uid[names[0]]) == set(by_uid[names[1]])
    compare_cols = ["connectivity_block", "standardized_inchikey", "frozen_er2_binary_label", "bootstrap_cluster",
                    "identity_clean", "model_usable", "primary_calibrated_probability", "direct_panel_probability",
                    "accepted_prediction", "risk_call", "abstention_reasons_json"]
    changed = {c: sum(by_uid[names[0]][u][c] != by_uid[names[1]][u][c] for u in by_uid[names[0]]) for c in compare_cols}
    all_results, all_contrasts, all_unavailable, all_ledger, audits = [], [], [], [], []
    for name, rows in inputs.items():
        r, c, u, ledger, audit = evaluate(name, rows)
        all_results.extend(r); all_contrasts.extend(c); all_unavailable.extend(u); all_ledger.extend(ledger); audits.append(audit)
    metric_signature = lambda name: [(r["strategy_id"], r["selected_n"], r["tp"], r["tn"], r["fp"], r["fn"], r["balanced_accuracy"]) for r in all_results if r["version"] == name]
    assert metric_signature(names[0]) == metric_signature(names[1])
    write_csv("metrics.csv", all_results)
    write_csv("paired_contrasts.csv", all_contrasts)
    write_csv("unattainable_targets.csv", all_unavailable)
    write_csv("selection_ledger.csv", all_ledger)
    audit_output = {"passed": True, "no_score_join_across_versions": True,
                    "v1_v2_common_operational_UID_n": len(by_uid[names[0]]), "v1_v2_changed_row_counts": changed,
                    "all_checked_columns_identical": all(v == 0 for v in changed.values()),
                    "v1_v2_strategy_point_metrics_identical": True, "populations": audits,
                    "checks": ["UID and connectivity uniqueness", "excluded all documented weak443/direct84 connectivity overlaps", "full mask reconstructs no remaining reasons", "accepted risk calls agree with same 0.5 class threshold", "full 453 point TP221 TN119 FP84 FN29 and BAcc0.735103 exactly reproduced", "whole cluster shared resampling across policies", "all matched contrasts selected_n equal", "larger-than-full coverage unavailable rather than filling abstentions"]}
    write_json("input_audit.json", audit_output)
    write_report(all_results, all_contrasts, audits)
    write_json("completion.json", {"completed_utc": datetime.now(timezone.utc).isoformat(), "status": "completed_retrospective_frozen_score_benchmark", "bootstrap_replicates": N_BOOTSTRAP,
                                    "versions": names, "output_sha256": {p.name: digest(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.name != "completion.json"}})
    print(json.dumps({"status": "completed", "path": str(HERE), "metric_rows": len(all_results), "paired_contrast_rows": len(all_contrasts), "audit_passed": True}, ensure_ascii=False))


if __name__ == "__main__":
    main()
