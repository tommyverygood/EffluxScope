#!/usr/bin/env python3
"""ASAP 2025 (CHEMBL5674852) source-held-out one-shot scoring for frozen ABCB1 models.

Executes `abcb1_asap2025_sourceheldout_protocol_lock_v1.json` exactly as locked on
2026-08-03. No refit, no recalibration, no threshold selection, no model choice is
performed with any CHEMBL5674852 outcome.

Primary endpoint : Spearman(model probability, log2(efflux_ratio)) on the estimable rows.
Directional gate : point estimate > 0 AND 95% Murcko-scaffold bootstrap lower bound > 0.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import joblib
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

PROJECT = Path(__file__).resolve().parents[2]
SCRIPTS = PROJECT / "数据" / "scripts"
PROCESSED = PROJECT / "数据" / "processed"
MODELS = PROJECT / "模型"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(MODELS / "ABCB1_assay_aware_candidate_v0_3_dev"))

from abcb1_screener_core_v0_1 import load_package, predict_smiles  # noqa: E402
from predict_abcb1_assay_aware import predict_one as predict_assay_aware  # noqa: E402

LOCK_PATH = PROCESSED / "abcb1_asap2025_sourceheldout_protocol_lock_v1.json"
OBS_PATH = PROCESSED / "abcb1_asap_mdck_mdr1_observations_v1.csv"
V01_PATH = MODELS / "ABCB1_screener_v0_1" / "ABCB1_screener_v0_1.joblib"
V03_PATH = MODELS / "ABCB1_assay_aware_candidate_v0_3_dev" / "ABCB1_assay_aware_candidate_v0_3_dev.joblib"

BOOTSTRAP_ITERATIONS = 2000
RANDOM_SEED = 20260803
WEAK_BINARY_CUT = 2.0

RDLogger.DisableLog("rdApp.*")


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: str):
    text = (value or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        ranks[order[i : j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def spearman(x: np.ndarray, y: np.ndarray):
    if len(x) < 3:
        return None
    rx, ry = rankdata(np.asarray(x, float)), rankdata(np.asarray(y, float))
    if np.std(rx) == 0 or np.std(ry) == 0:
        return None
    return float(np.corrcoef(rx, ry)[0, 1])


def roc_auc(labels: np.ndarray, scores: np.ndarray):
    pos, neg = scores[labels == 1], scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    ranks = rankdata(scores)
    return float((ranks[labels == 1].sum() - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg)))


def balanced_accuracy(labels: np.ndarray, predictions: np.ndarray):
    pos, neg = labels == 1, labels == 0
    if pos.sum() == 0 or neg.sum() == 0:
        return None
    return float(0.5 * (predictions[pos].mean() + (1 - predictions[neg]).mean()))


def brier(labels: np.ndarray, probabilities: np.ndarray) -> float:
    return float(np.mean((probabilities - labels) ** 2))


def murcko_scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return "__unparsed__"
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or "__acyclic__"
    except Exception:
        return "__unparsed__"


def scaffold_bootstrap_ci(groups: dict, statistic, iterations: int, seed: int):
    """Resample whole Murcko-scaffold clusters, per protocol prohibition #4."""
    keys = sorted(groups)
    if len(keys) < 2:
        return None, None, 0
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(iterations):
        picked = rng.choice(len(keys), size=len(keys), replace=True)
        rows: list[dict] = []
        for idx in picked:
            rows.extend(groups[keys[idx]])
        value = statistic(rows)
        if value is not None and not math.isnan(value):
            samples.append(value)
    if len(samples) < iterations * 0.5:
        return None, None, len(samples)
    arr = np.asarray(samples, float)
    return float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5)), len(samples)


# --------------------------------------------------------------------------------------
# step 1 - pre-scoring validation
# --------------------------------------------------------------------------------------
def validate(lock: dict) -> tuple[list[dict], list[dict], dict]:
    checks: dict = {}

    obs_hash = sha256(OBS_PATH)
    checks["observations_sha256"] = obs_hash
    checks["observations_sha256_matches_lock"] = obs_hash == lock["processed_interface"]["sha256"]

    for label, path, expected in (
        ("v0_1", V01_PATH, lock["frozen_models_allowed_for_first_scoring"][0]["sha256"]),
        ("v0_3_dev", V03_PATH, lock["frozen_models_allowed_for_first_scoring"][1]["sha256"]),
    ):
        digest = sha256(path)
        checks[f"model_{label}_sha256"] = digest
        checks[f"model_{label}_sha256_matches_lock"] = digest == expected

    rows = read_csv(OBS_PATH)
    checks["observation_rows_total"] = len(rows)

    doc = lock["source_document"]["document_chembl_id"]
    assay = lock["source_document"]["assay_chembl_id"]
    scoped = [r for r in rows if r["document_chembl_id"] == doc and r["assay_chembl_id"] == assay]
    checks["scoped_rows"] = len(scoped)

    eligible = [r for r in scoped if as_float(r["papp_ab"]) is not None and as_float(r["papp_ba"]) is not None]
    checks["eligible_bidirectional_rows"] = len(eligible)
    checks["eligible_matches_lock"] = len(eligible) == lock["source_document"]["bidirectional_observation_rows"]

    estimable = [r for r in eligible if as_float(r["efflux_ratio"]) is not None and as_float(r["efflux_ratio"]) > 0]
    checks["estimable_efflux_ratio_rows"] = len(estimable)
    checks["estimable_matches_lock"] = len(estimable) == lock["source_document"]["estimable_efflux_ratio_rows"]

    checks["identity_status_counts"] = dict(Counter(r["identity_status"] for r in estimable))
    checks["overlap_exact_prior_sets_nonblank"] = sum(1 for r in estimable if (r["overlap_exact_prior_sets"] or "").strip())
    checks["overlap_connectivity_prior_sets_nonblank"] = sum(
        1 for r in estimable if (r["overlap_connectivity_prior_sets"] or "").strip()
    )
    checks["model_scoring_performed_flags"] = dict(Counter(r["model_scoring_performed"] for r in estimable))
    checks["unique_molecules_estimable"] = len({r["standardized_inchikey"] for r in estimable})
    checks["weak_binary_label_counts"] = dict(Counter(r["weak_binary_label"] for r in estimable))

    blocking = [
        "observations_sha256_matches_lock",
        "model_v0_1_sha256_matches_lock",
        "model_v0_3_dev_sha256_matches_lock",
        "eligible_matches_lock",
        "estimable_matches_lock",
    ]
    checks["all_blocking_checks_passed"] = all(bool(checks[k]) for k in blocking)
    checks["identity_audit_clean"] = (
        checks["overlap_exact_prior_sets_nonblank"] == 0 and checks["overlap_connectivity_prior_sets_nonblank"] == 0
    )
    return eligible, estimable, checks


# --------------------------------------------------------------------------------------
# step 2 - one-shot scoring
# --------------------------------------------------------------------------------------
def score(estimable: list[dict]) -> list[dict]:
    v01 = load_package(V01_PATH)
    v03 = joblib.load(V03_PATH)

    scored: list[dict] = []
    for row in estimable:
        smiles = row["standardized_isomeric_smiles"] or row["input_smiles"]
        er = as_float(row["efflux_ratio"])
        p1 = predict_smiles(smiles, v01)
        p3 = predict_assay_aware(smiles, v03)
        scored.append(
            {
                "observation_id": row["observation_id"],
                "molecule_chembl_id": row["molecule_chembl_id"],
                "molecule_pref_name": row["molecule_pref_name"],
                "standardized_inchikey": row["standardized_inchikey"],
                "murcko_scaffold": murcko_scaffold(smiles),
                "papp_ab": row["papp_ab"],
                "papp_ba": row["papp_ba"],
                "efflux_ratio": er,
                "log2_efflux_ratio": math.log2(er),
                "weak_binary_label_er_ge_2": int(er >= WEAK_BINARY_CUT),
                "reconstructed_class": row["reconstructed_class"],
                "evidence_tier": row["evidence_tier"],
                "v0_1_calibrated_probability": p1.get("calibrated_probability", ""),
                "v0_1_raw_probability": p1.get("raw_probability", ""),
                "v0_1_binary_prediction": p1.get("binary_prediction", ""),
                "v0_1_binary_threshold": p1.get("binary_threshold", ""),
                "v0_1_max_training_tanimoto": p1.get("max_training_morgan_tanimoto", ""),
                "v0_1_applicability_band": p1.get("applicability_band", ""),
                "v0_1_screening_call": p1.get("screening_call", ""),
                "v0_3_candidate_probability": p3.get("candidate_probability", ""),
                "v0_3_binary_prediction": p3.get("candidate_binary_prediction", ""),
                "v0_3_binary_threshold": p3.get("binary_threshold", ""),
                "v0_3_max_direct_panel_tanimoto": p3.get("max_direct_panel_tanimoto", ""),
                "v0_3_applicability_band": p3.get("applicability_band", ""),
                "v0_3_triage_call": p3.get("triage_call", ""),
                "v0_3_exact_development_identity_match": p3.get("exact_development_identity_match", ""),
                "claim_boundary": "source-held-out continuous evaluation; not a universal ABCB1 substrate label",
            }
        )
    return scored


# --------------------------------------------------------------------------------------
# step 3 - metrics + predeclared gate
# --------------------------------------------------------------------------------------
def build_baselines(scored: list[dict]) -> dict:
    """Protocol-listed comparators that need no CHEMBL5674852 outcome."""
    v01 = load_package(V01_PATH)
    prior = float(np.mean([r["target_label"] for r in v01["training_rows"]])) if "training_rows" in v01 else None
    return {"training_class_prior": prior}


def evaluate(scored: list[dict], model_key: str, probability_field: str) -> list[dict]:
    usable = [r for r in scored if r[probability_field] not in ("", None)]
    if len(usable) < 3:
        return []

    prob = np.asarray([float(r[probability_field]) for r in usable], float)
    log2er = np.asarray([float(r["log2_efflux_ratio"]) for r in usable], float)
    weak = np.asarray([int(r["weak_binary_label_er_ge_2"]) for r in usable], int)

    groups: dict[str, list[dict]] = defaultdict(list)
    for r in usable:
        groups[r["murcko_scaffold"]].append(r)

    out: list[dict] = []

    def add(metric, value, lo=None, hi=None, n=None, extra=None):
        row = {
            "model": model_key,
            "probability_field": probability_field,
            "subset": "estimable_efflux_ratio_rows",
            "metric": metric,
            "value": "" if value is None else round(float(value), 8),
            "ci95_lower": "" if lo is None else round(float(lo), 8),
            "ci95_upper": "" if hi is None else round(float(hi), 8),
            "n": n if n is not None else len(usable),
            "n_scaffold_clusters": len(groups),
            "bootstrap_iterations": BOOTSTRAP_ITERATIONS,
            "random_seed": RANDOM_SEED,
        }
        if extra:
            row.update(extra)
        out.append(row)

    # primary endpoint
    rho = spearman(prob, log2er)
    lo, hi, n_ok = scaffold_bootstrap_ci(
        groups,
        lambda rows: spearman(
            np.asarray([float(r[probability_field]) for r in rows], float),
            np.asarray([float(r["log2_efflux_ratio"]) for r in rows], float),
        ),
        BOOTSTRAP_ITERATIONS,
        RANDOM_SEED,
    )
    add("spearman_rho_vs_log2_efflux_ratio", rho, lo, hi, extra={"bootstrap_samples_used": n_ok, "endpoint_role": "primary"})

    # secondary endpoints
    auc = roc_auc(weak, prob)
    lo_a, hi_a, n_a = scaffold_bootstrap_ci(
        groups,
        lambda rows: roc_auc(
            np.asarray([int(r["weak_binary_label_er_ge_2"]) for r in rows], int),
            np.asarray([float(r[probability_field]) for r in rows], float),
        ),
        BOOTSTRAP_ITERATIONS,
        RANDOM_SEED,
    )
    add("roc_auc_weak_er_ge_2", auc, lo_a, hi_a, extra={"bootstrap_samples_used": n_a, "endpoint_role": "secondary"})

    threshold_field = "v0_1_binary_threshold" if model_key.startswith("frozen_ABCB1_v0_1") else "v0_3_binary_threshold"
    threshold = as_float(str(usable[0].get(threshold_field, "")))
    if threshold is not None:
        preds = (prob >= threshold).astype(int)
        add(
            "balanced_accuracy_at_frozen_threshold",
            balanced_accuracy(weak, preds),
            extra={"frozen_threshold": threshold, "endpoint_role": "secondary"},
        )
        add(
            "sensitivity_at_frozen_threshold",
            float(preds[weak == 1].mean()) if (weak == 1).any() else None,
            extra={"frozen_threshold": threshold, "endpoint_role": "secondary"},
        )
        add(
            "specificity_at_frozen_threshold",
            float(1 - preds[weak == 0].mean()) if (weak == 0).any() else None,
            extra={"frozen_threshold": threshold, "endpoint_role": "secondary"},
        )
    add("brier_score_weak_er_ge_2", brier(weak, prob), extra={"endpoint_role": "secondary"})

    # coverage / metrics by frozen applicability band
    band_field = "v0_1_applicability_band" if model_key.startswith("frozen_ABCB1_v0_1") else "v0_3_applicability_band"
    band_counts = Counter(r[band_field] for r in usable)
    abstain_tokens = ("outside_empirical_support", "outside_direct_panel_empirical_support")
    n_abstain = sum(v for k, v in band_counts.items() if k in abstain_tokens)
    add(
        "nonabstain_coverage",
        (len(usable) - n_abstain) / len(usable),
        extra={"band_counts": json.dumps(dict(band_counts), ensure_ascii=False), "endpoint_role": "secondary"},
    )
    for band, count in sorted(band_counts.items()):
        sub = [r for r in usable if r[band_field] == band]
        if len(sub) >= 3:
            sub_rho = spearman(
                np.asarray([float(r[probability_field]) for r in sub], float),
                np.asarray([float(r["log2_efflux_ratio"]) for r in sub], float),
            )
        else:
            sub_rho = None
        add(f"spearman_rho_within_band::{band}", sub_rho, n=count, extra={"endpoint_role": "secondary"})

    return out


def main() -> None:
    lock = json.loads(LOCK_PATH.read_text(encoding="utf-8"))
    eligible, estimable, checks = validate(lock)

    print("=== step 1  评分前校验 ===")
    for key in (
        "observations_sha256_matches_lock",
        "model_v0_1_sha256_matches_lock",
        "model_v0_3_dev_sha256_matches_lock",
        "scoped_rows",
        "eligible_bidirectional_rows",
        "eligible_matches_lock",
        "estimable_efflux_ratio_rows",
        "estimable_matches_lock",
        "unique_molecules_estimable",
        "overlap_exact_prior_sets_nonblank",
        "overlap_connectivity_prior_sets_nonblank",
        "identity_audit_clean",
        "all_blocking_checks_passed",
    ):
        print(f"  {key:42s} {checks[key]}")
    print(f"  weak_binary_label_counts                   {checks['weak_binary_label_counts']}")

    if not checks["all_blocking_checks_passed"]:
        print("\n阻断性校验未通过，按协议停止，不进行评分。")
        (PROCESSED / "abcb1_asap2025_sourceheldout_summary_v1.json").write_text(
            json.dumps({"status": "halted_precheck_failed", "checks": checks}, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        sys.exit(1)

    print("\n=== step 2  一次性评分 ===")
    scored = score(estimable)
    write_csv(PROCESSED / "abcb1_asap2025_sourceheldout_predictions_v1.csv", scored)
    print(f"  已写出逐行预测 {len(scored)} 行")

    print("\n=== step 3  指标与预声明方向门 ===")
    metrics: list[dict] = []
    metrics += evaluate(scored, "frozen_ABCB1_v0_1", "v0_1_calibrated_probability")
    metrics += evaluate(scored, "frozen_ABCB1_assay_aware_v0_3_dev", "v0_3_candidate_probability")
    write_csv(PROCESSED / "abcb1_asap2025_sourceheldout_metrics_v1.csv", metrics)

    gate_results = {}
    for model_key in ("frozen_ABCB1_v0_1", "frozen_ABCB1_assay_aware_v0_3_dev"):
        primary = next(
            (m for m in metrics if m["model"] == model_key and m["metric"] == "spearman_rho_vs_log2_efflux_ratio"),
            None,
        )
        if primary is None or primary["value"] == "":
            gate_results[model_key] = {"status": "not_estimable"}
            continue
        rho = float(primary["value"])
        lower = primary["ci95_lower"]
        lower_val = float(lower) if lower != "" else None
        passed = rho > 0 and lower_val is not None and lower_val > 0
        gate_results[model_key] = {
            "spearman_rho": rho,
            "ci95_lower": lower_val,
            "ci95_upper": float(primary["ci95_upper"]) if primary["ci95_upper"] != "" else None,
            "n": primary["n"],
            "n_scaffold_clusters": primary["n_scaffold_clusters"],
            "status": "pass" if passed else "fail",
            "rule": lock["predeclared_directional_gate"]["rule"],
            "interpretation": lock["predeclared_directional_gate"][
                "interpretation_if_pass" if passed else "interpretation_if_fail"
            ],
        }

    for model_key, res in gate_results.items():
        print(f"\n  【{model_key}】")
        for k, v in res.items():
            print(f"    {k}: {v}")

    summary = {
        "protocol_version": lock["protocol_version"],
        "scoring_date": "2026-08-03",
        "status": "scored_once_no_tuning",
        "script": str(Path(__file__).relative_to(PROJECT)),
        "random_seed": RANDOM_SEED,
        "bootstrap_iterations": BOOTSTRAP_ITERATIONS,
        "precheck": checks,
        "primary_endpoint": lock["primary_endpoint"],
        "predeclared_directional_gate": lock["predeclared_directional_gate"],
        "gate_results": gate_results,
        "cumulative_stopping_rule": (
            "连续端点方向此前已有两次不支持：Roth 2026 面板 ABCB1 Morgan Spearman = -0.520；"
            "Langthaler 2024 n=11 rho = +0.2096 (p = 0.536)。本次 ASAP 2025 为第三次独立检验。"
            "若第三次仍不通过方向门，宣布『普通二分类训练知识不能迁移到连续外排端点』，"
            "停止继续堆积弱标签，转向针对 ADC 化学空间的少量实验数据或证据分层/拒判系统，不再对该命题复测。"
        ),
        "prohibitions_observed": lock["prohibitions"],
        "outputs": {
            "predictions": "数据/processed/abcb1_asap2025_sourceheldout_predictions_v1.csv",
            "metrics": "数据/processed/abcb1_asap2025_sourceheldout_metrics_v1.csv",
        },
        "claim_boundary": lock["claim_boundary"],
    }
    (PROCESSED / "abcb1_asap2025_sourceheldout_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\n已写出 summary JSON。")


if __name__ == "__main__":
    main()
