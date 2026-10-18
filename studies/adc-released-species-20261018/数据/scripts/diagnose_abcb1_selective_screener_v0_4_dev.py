#!/usr/bin/env python3
"""Run frozen-rule, explicitly post-hoc diagnostics for v0.4-dev.

The regulatory v4 and ASAP 2025 outcomes were observed before the v0.4 rule
was locked. They are therefore diagnostics, not external validation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr


DATA = Path(__file__).resolve().parents[1]
PROJECT = DATA.parent
PROCESSED = DATA / "processed"
MODEL_DIR = PROJECT / "模型" / "ABCB1_selective_screener_v0_4_dev"
sys.path.insert(0, str(MODEL_DIR))

from predict_abcb1_selective import flatten_batch_result, load_packages, predict_one  # noqa: E402


SELECTIVE_MODEL = MODEL_DIR / "ABCB1_selective_uncertainty_v0_4_dev.joblib"
PRIMARY_MODEL = PROJECT / "模型" / "ABCB1_screener_v0_1" / "ABCB1_screener_v0_1.joblib"
DIRECT_MODEL = PROJECT / "模型" / "ABCB1_assay_aware_candidate_v0_3_dev" / "ABCB1_assay_aware_candidate_v0_3_dev.joblib"
REGULATORY = PROCESSED / "abcb1_regulatory_v4_blind_predictions_v1.csv"
ASAP_PRED = PROCESSED / "abcb1_asap2025_sourceheldout_predictions_v1.csv"
ASAP_OBS = PROCESSED / "abcb1_asap_mdck_mdr1_observations_v1.csv"
CONRAD_ENTITIES = PROCESSED / "abcb1_conrad2024_matched_pair_entities_v1.csv"
CONRAD_PAIRS = PROCESSED / "abcb1_conrad2024_figure5_pairs_v1.csv"

OUT_REG = PROCESSED / "abcb1_selective_v0_4_regulatory_v4_posthoc_predictions_v1.csv"
OUT_ASAP = PROCESSED / "abcb1_selective_v0_4_asap2025_posthoc_predictions_v1.csv"
OUT_CONRAD = PROCESSED / "abcb1_selective_v0_4_conrad_posthoc_predictions_v1.csv"
OUT_PAIR = PROCESSED / "abcb1_selective_v0_4_conrad_pairwise_posthoc_v1.csv"
OUT_METRICS = PROCESSED / "abcb1_selective_v0_4_posthoc_diagnostic_metrics_v1.csv"
OUT_SUMMARY = PROCESSED / "abcb1_selective_v0_4_posthoc_diagnostic_summary_v1.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def dump_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def safe_spearman(x: list[float], y: list[float]) -> float | None:
    if len(x) < 3 or len(set(x)) < 2 or len(set(y)) < 2:
        return None
    value = float(spearmanr(x, y).statistic)
    return value if math.isfinite(value) else None


def directional_metrics(rows: list[dict], label_field: str) -> dict:
    accepted = [row for row in rows if row["risk_call"] != "abstain"]
    high = [row for row in accepted if row["risk_call"] == "high_ABCB1_efflux_risk"]
    low = [row for row in accepted if row["risk_call"] == "low_ABCB1_efflux_risk"]
    correct = [
        row
        for row in accepted
        if (row["risk_call"] == "high_ABCB1_efflux_risk") == (int(float(row[label_field])) == 1)
    ]
    high_correct = [row for row in high if int(float(row[label_field])) == 1]
    low_correct = [row for row in low if int(float(row[label_field])) == 0]
    return {
        "n": len(rows),
        "accepted_n": len(accepted),
        "abstain_n": len(rows) - len(accepted),
        "coverage": len(accepted) / len(rows) if rows else None,
        "accepted_accuracy": len(correct) / len(accepted) if accepted else None,
        "high_n": len(high),
        "high_precision": len(high_correct) / len(high) if high else None,
        "low_n": len(low),
        "low_negative_predictive_value": len(low_correct) / len(low) if low else None,
    }


def metric_rows(dataset: str, metrics: dict, endpoint: str) -> list[dict]:
    return [
        {
            "dataset": dataset,
            "endpoint": endpoint,
            "metric": key,
            "value": value if value is not None else "",
            "status": "historical_posthoc_diagnostic_not_external_validation",
        }
        for key, value in metrics.items()
    ]


def main() -> None:
    selective, primary, direct = load_packages(SELECTIVE_MODEL, PRIMARY_MODEL, DIRECT_MODEL)

    regulatory_rows = read_csv(REGULATORY)
    scored_regulatory: list[dict] = []
    for row in regulatory_rows:
        result = flatten_batch_result(predict_one(row["standardized_isomeric_smiles"], selective, primary, direct))
        scored_regulatory.append({**row, **result})
    write_csv(OUT_REG, scored_regulatory)

    observations = {row["observation_id"]: row for row in read_csv(ASAP_OBS)}
    scored_asap: list[dict] = []
    for row in read_csv(ASAP_PRED):
        observation = observations[row["observation_id"]]
        result = flatten_batch_result(
            predict_one(observation["standardized_isomeric_smiles"], selective, primary, direct)
        )
        scored_asap.append({**row, "standardized_isomeric_smiles": observation["standardized_isomeric_smiles"], **result})
    write_csv(OUT_ASAP, scored_asap)

    conrad_rows = read_csv(CONRAD_ENTITIES)
    scored_conrad: list[dict] = []
    for row in conrad_rows:
        result = flatten_batch_result(predict_one(row["standardized_isomeric_smiles"], selective, primary, direct))
        scored_conrad.append({**row, **result})
    write_csv(OUT_CONRAD, scored_conrad)

    conrad_by_key = {row["standardized_inchikey"]: row for row in scored_conrad}
    pair_rows: list[dict] = []
    for pair in read_csv(CONRAD_PAIRS):
        left = conrad_by_key[pair["left_standardized_inchikey"]]
        right = conrad_by_key[pair["right_standardized_inchikey"]]
        left_probability = float(left["primary_calibrated_probability"])
        right_probability = float(right["primary_calibrated_probability"])
        predicted_higher = pair["left_compound_id"] if left_probability > right_probability else pair["right_compound_id"]
        expected_higher = pair["higher_efflux_compound"]
        pair_rows.append(
            {
                **pair,
                "left_primary_probability": round(left_probability, 8),
                "right_primary_probability": round(right_probability, 8),
                "predicted_higher_efflux_compound": predicted_higher,
                "pair_direction_correct": "yes" if predicted_higher == expected_higher else "no",
                "left_risk_call": left["risk_call"],
                "right_risk_call": right["risk_call"],
                "posthoc_status": "historical_posthoc_diagnostic_not_external_validation",
            }
        )
    write_csv(OUT_PAIR, pair_rows)

    reg_metrics = directional_metrics(scored_regulatory, "directional_label")
    asap_metrics = directional_metrics(scored_asap, "weak_binary_label_er_ge_2")
    asap_all_rho = safe_spearman(
        [float(row["primary_calibrated_probability"]) for row in scored_asap],
        [float(row["log2_efflux_ratio"]) for row in scored_asap],
    )
    asap_accepted = [row for row in scored_asap if row["risk_call"] != "abstain"]
    asap_accepted_rho = safe_spearman(
        [float(row["primary_calibrated_probability"]) for row in asap_accepted],
        [float(row["log2_efflux_ratio"]) for row in asap_accepted],
    )
    clean_pairs = [row for row in pair_rows if row["identity_clean_pair"] == "yes"]
    clean_pair_correct = sum(row["pair_direction_correct"] == "yes" for row in clean_pairs)
    pair_metrics = {
        "identity_clean_pair_n": len(clean_pairs),
        "pair_direction_correct_n": clean_pair_correct,
        "pair_direction_accuracy": clean_pair_correct / len(clean_pairs) if clean_pairs else None,
    }
    metrics = [
        *metric_rows("regulatory_v4", reg_metrics, "strict directional source call"),
        *metric_rows("ASAP_2025", asap_metrics, "weak ER >= 2 sensitivity call"),
        *metric_rows(
            "ASAP_2025",
            {
                "spearman_all_rows": asap_all_rho,
                "spearman_accepted_rows": asap_accepted_rho,
                "accepted_rows_for_spearman": len(asap_accepted),
            },
            "continuous log2 efflux ratio",
        ),
        *metric_rows("Conrad_2024_Figure5", pair_metrics, "within-series pair direction"),
    ]
    write_csv(OUT_METRICS, metrics)

    summary = {
        "model_version": selective["model_version"],
        "diagnostic_status": "historical_posthoc_diagnostic_not_external_validation",
        "rule_tuned_on_historical_external_outcomes": False,
        "package_sha256": sha256_file(SELECTIVE_MODEL),
        "input_hashes": {
            str(REGULATORY.relative_to(PROJECT)): sha256_file(REGULATORY),
            str(ASAP_PRED.relative_to(PROJECT)): sha256_file(ASAP_PRED),
            str(ASAP_OBS.relative_to(PROJECT)): sha256_file(ASAP_OBS),
            str(CONRAD_ENTITIES.relative_to(PROJECT)): sha256_file(CONRAD_ENTITIES),
            str(CONRAD_PAIRS.relative_to(PROJECT)): sha256_file(CONRAD_PAIRS),
        },
        "regulatory_v4": reg_metrics,
        "asap_2025_weak_binary": asap_metrics,
        "asap_2025_continuous": {
            "spearman_all_rows": asap_all_rho,
            "spearman_accepted_rows": asap_accepted_rho,
            "accepted_rows_for_spearman": len(asap_accepted),
        },
        "conrad_figure5": pair_metrics,
        "abstention_reason_counts": {
            "regulatory_v4": dict(Counter(reason for row in scored_regulatory for reason in json.loads(row["abstention_reasons_json"]))),
            "asap_2025": dict(Counter(reason for row in scored_asap for reason in json.loads(row["abstention_reasons_json"]))),
        },
        "outputs": {
            "regulatory_predictions": str(OUT_REG.relative_to(PROJECT)),
            "asap_predictions": str(OUT_ASAP.relative_to(PROJECT)),
            "conrad_predictions": str(OUT_CONRAD.relative_to(PROJECT)),
            "conrad_pairs": str(OUT_PAIR.relative_to(PROJECT)),
            "metrics": str(OUT_METRICS.relative_to(PROJECT)),
        },
        "formal_interpretation": (
            "The three-output implementation is operational and conservative. Historical diagnostics quantify behavior only; "
            "a future outcome-unseen source is still required before claiming external selective-classification validity."
        ),
        "claim_boundary": selective["claim_boundary"],
    }
    dump_json(OUT_SUMMARY, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
