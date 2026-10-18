#!/usr/bin/env python3
"""Run frozen evidence-channel ablations for the ABCB1 v0.4 rejector."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_abcb1_external_risk_experiments_v1 import (
    ABSTAIN,
    CLUSTER_COL,
    HIGH,
    LOW,
    N_BOOT,
    PROJECT,
    PROCESSED,
    SEED,
    Y_COL,
    clean_frame,
    counts_from_arrays,
    metrics_for_frame,
    result_row,
    round_frame,
    vector_metrics,
)


INPUT = PROCESSED / "abcb1_quantity_biogen_analysis_rows_v1.csv"
MODEL = PROJECT / "模型" / "ABCB1_selective_screener_v0_4_dev" / "ABCB1_selective_uncertainty_v0_4_dev.joblib"
PRIOR_PROTOCOL = PROCESSED / "abcb1_external_risk_experiment_protocol_v1.json"
PROTOCOL = PROCESSED / "abcb1_external_evidence_channel_ablation_protocol_v1.json"
METRICS = PROCESSED / "abcb1_external_evidence_channel_ablation_metrics_v1.csv"
SET_CONTRASTS = PROCESSED / "abcb1_external_evidence_channel_set_contrasts_v1.csv"
PAIRWISE_BOOTSTRAP = PROCESSED / "abcb1_external_evidence_channel_pairwise_bootstrap_v1.csv"
SUMMARY = PROCESSED / "abcb1_external_evidence_channel_ablation_summary_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def component_masks(work: pd.DataFrame, protocol: dict[str, object]) -> dict[str, np.ndarray]:
    thresholds = protocol["frozen_thresholds"]
    reasons = work["reason_set"]
    primary_supported = pd.to_numeric(work["primary_max_tanimoto"], errors="coerce").ge(float(thresholds["primary_similarity_support"]))
    direct_supported = pd.to_numeric(work["direct_max_tanimoto"], errors="coerce").ge(float(thresholds["direct_similarity_support"]))
    masks = {
        "margin": work["pass_probability_margin"],
        "bootstrap": work["pass_bootstrap_direction_stability"],
        "primary_support": primary_supported,
        "direct_support": direct_supported,
        "dual_support": primary_supported | direct_supported,
        "physchem": reasons.map(lambda value: "outside_physicochemical_envelope" not in value),
        "cross_model": reasons.map(lambda value: "supported_model_direction_disagreement" not in value),
        "primary_local": reasons.map(lambda value: "primary_local_neighborhood_strongly_opposite" not in value),
        "direct_local": reasons.map(lambda value: "direct_local_neighborhood_strongly_opposite" not in value),
        "stereo": reasons.map(lambda value: "stereo_blind_fingerprint_collision" not in value),
    }
    return {key: value.to_numpy(dtype=bool) for key, value in masks.items()}


def variant_masks(work: pd.DataFrame, protocol: dict[str, object]) -> dict[str, np.ndarray]:
    components = component_masks(work, protocol)
    variants: dict[str, np.ndarray] = {}
    for variant in protocol["locked_variants"]:
        mask = np.ones(len(work), dtype=bool)
        for component in variant["components"]:
            mask &= components[component]
        variants[variant["variant_id"]] = mask
    return variants


def variant_metrics(work: pd.DataFrame, population: str, protocol: dict[str, object]) -> tuple[list[dict[str, object]], dict[str, np.ndarray]]:
    masks = variant_masks(work, protocol)
    rows = []
    for index, variant in enumerate(protocol["locked_variants"]):
        variant_id = variant["variant_id"]
        selected = work[masks[variant_id]].copy()
        row = result_row(
            selected,
            len(work),
            population=population,
            experiment="evidence_channel_ablation",
            label=variant_id,
            seed_offset=12000 + 100 * (population == "scaffold_clean") + index,
        )
        row["variant_id"] = variant_id
        row["components"] = "|".join(variant["components"])
        row["accepted_high_n"] = int((selected["forced_prediction"] == 1).sum())
        row["accepted_low_n"] = int((selected["forced_prediction"] == 0).sum())
        rows.append(row)
    return rows, masks


def set_contrast_rows(
    work: pd.DataFrame,
    population: str,
    masks: dict[str, np.ndarray],
    protocol: dict[str, object],
) -> list[dict[str, object]]:
    rows = []
    for full_id, comparator_id in protocol["locked_set_contrasts"]:
        full = masks[full_id]
        comparator = masks[comparator_id]
        regions = {
            "intersection": full & comparator,
            "full_only": full & ~comparator,
            "comparator_only": ~full & comparator,
            "neither": ~full & ~comparator,
        }
        union_n = int((full | comparator).sum())
        for region, mask in regions.items():
            frame = work[mask].copy()
            row = {
                "population": population,
                "full_variant": full_id,
                "comparator_variant": comparator_id,
                "region": region,
                "n": int(len(frame)),
                "fraction_of_population": float(len(frame) / len(work)),
                "fraction_of_union": float(len(frame) / union_n) if union_n else np.nan,
                "full_accepted_n": int(full.sum()),
                "comparator_accepted_n": int(comparator.sum()),
                "accepted_set_jaccard": float((full & comparator).sum() / union_n) if union_n else np.nan,
            }
            if len(frame):
                row.update(metrics_for_frame(frame))
            rows.append(row)
    return rows


def pairwise_bootstrap_rows(
    work: pd.DataFrame,
    population: str,
    masks: dict[str, np.ndarray],
    protocol: dict[str, object],
) -> list[dict[str, object]]:
    rows = []
    cluster_groups = list(work.groupby(CLUSTER_COL, sort=False).groups.values())
    for contrast_index, (full_id, comparator_id) in enumerate(protocol["locked_set_contrasts"]):
        full_mask = masks[full_id]
        comparator_mask = masks[comparator_id]
        grouped = []
        for indices in cluster_groups:
            idx = np.asarray(list(indices), dtype=int)
            full_frame = work.iloc[idx[full_mask[idx]]]
            comparator_frame = work.iloc[idx[comparator_mask[idx]]]
            grouped.append(np.concatenate([
                counts_from_arrays(full_frame[Y_COL].to_numpy(dtype=int), full_frame["forced_prediction"].to_numpy(dtype=int)),
                counts_from_arrays(comparator_frame[Y_COL].to_numpy(dtype=int), comparator_frame["forced_prediction"].to_numpy(dtype=int)),
            ]))
        grouped_array = np.asarray(grouped, dtype=float)
        rng = np.random.default_rng(SEED + 13000 + 100 * (population == "scaffold_clean") + contrast_index)
        bootstrap_values = {"accuracy": [], "balanced_accuracy": [], "mcc": []}
        remaining = N_BOOT
        while remaining:
            batch = min(100, remaining)
            sampled = rng.integers(0, len(grouped_array), size=(batch, len(grouped_array)))
            totals = grouped_array[sampled].sum(axis=1)
            full_metrics = vector_metrics(totals[:, :4])
            comparator_metrics = vector_metrics(totals[:, 4:])
            for metric in bootstrap_values:
                delta = full_metrics[metric] - comparator_metrics[metric]
                bootstrap_values[metric].extend(delta[np.isfinite(delta)].tolist())
            remaining -= batch

        full_point = metrics_for_frame(work[full_mask])
        comparator_point = metrics_for_frame(work[comparator_mask])
        for metric, values in bootstrap_values.items():
            array = np.asarray(values, dtype=float)
            point = full_point[metric] - comparator_point[metric]
            rows.append({
                "population": population,
                "full_variant": full_id,
                "comparator_variant": comparator_id,
                "metric": metric,
                "full_value": full_point[metric],
                "comparator_value": comparator_point[metric],
                "full_minus_comparator": point,
                "ci95_low": float(np.quantile(array, 0.025)),
                "ci95_high": float(np.quantile(array, 0.975)),
                "bootstrap_probability_full_gt_comparator": float(np.mean(array > 0)),
                "valid_bootstrap_replicates": int(len(array)),
                "n_bootstrap_clusters": int(len(grouped_array)),
            })
    return rows


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text())
    assert protocol["status"].startswith("locked_posthoc_component_ablation")
    for entry in protocol["frozen_inputs"].values():
        path = PROJECT / entry["path"]
        assert path.exists(), path
        assert sha256(path) == entry["sha256"], f"Hash mismatch: {path}"
    assert protocol["bootstrap"]["replicates"] == N_BOOT
    assert protocol["bootstrap"]["seed"] == SEED

    raw = pd.read_csv(INPUT, low_memory=False)
    metric_rows = []
    contrast_rows = []
    pairwise_rows = []
    population_masks: dict[str, dict[str, np.ndarray]] = {}
    populations: dict[str, pd.DataFrame] = {}
    for population in protocol["analysis_populations"]:
        work = clean_frame(raw, population)
        populations[population] = work
        rows, masks = variant_metrics(work, population, protocol)
        metric_rows.extend(rows)
        contrast_rows.extend(set_contrast_rows(work, population, masks, protocol))
        pairwise_rows.extend(pairwise_bootstrap_rows(work, population, masks, protocol))
        population_masks[population] = masks

    metrics = round_frame(pd.DataFrame(metric_rows))
    contrasts = round_frame(pd.DataFrame(contrast_rows))
    pairwise = round_frame(pd.DataFrame(pairwise_rows))
    metrics.to_csv(METRICS, index=False)
    contrasts.to_csv(SET_CONTRASTS, index=False)
    pairwise.to_csv(PAIRWISE_BOOTSTRAP, index=False)

    identity_metrics = metrics[metrics["population"].eq("identity_clean")].set_index("variant_id")
    full = identity_metrics.loc["V11_full_v0_4"]
    frozen_ids = set(populations["identity_clean"].loc[
        populations["identity_clean"]["risk_call"].isin([HIGH, LOW]), "external_row_id"
    ].astype(str))
    reconstructed_ids = set(populations["identity_clean"].loc[
        population_masks["identity_clean"]["V11_full_v0_4"], "external_row_id"
    ].astype(str))
    descriptive_best = identity_metrics.sort_values("balanced_accuracy", ascending=False).iloc[0]

    def selected_metrics(variant_id: str) -> dict[str, object]:
        row = identity_metrics.loc[variant_id]
        keys = [
            "accepted_n", "coverage", "accuracy", "balanced_accuracy", "mcc", "sensitivity", "specificity",
            "high_risk_precision", "low_risk_npv", "accepted_high_n", "accepted_low_n",
            "balanced_accuracy_ci95_low", "balanced_accuracy_ci95_high", "mcc_ci95_low", "mcc_ci95_high",
        ]
        return {key: row[key].item() if isinstance(row[key], np.generic) else row[key] for key in keys}

    summary = {
        "analysis_version": "ABCB1_EXTERNAL_EVIDENCE_CHANNEL_ABLATION_V1",
        "analysis_date": "2026-08-05",
        "status": "completed_frozen_model_posthoc_evidence_channel_ablation",
        "protocol_path": str(PROTOCOL.relative_to(PROJECT)),
        "protocol_sha256": sha256(PROTOCOL),
        "full_v0_4_reconstruction": {
            "exactly_reproduced": frozen_ids == reconstructed_ids,
            "frozen_accepted_n": len(frozen_ids),
            "reconstructed_accepted_n": len(reconstructed_ids),
            "symmetric_difference_n": len(frozen_ids ^ reconstructed_ids),
            "balanced_accuracy": full["balanced_accuracy"],
            "mcc": full["mcc"],
        },
        "identity_clean_key_variants": {
            variant_id: selected_metrics(variant_id)
            for variant_id in [
                "V00_forced_binary_all",
                "V02_prior_stability_only",
                "V04_primary_channel_complete",
                "V06_direct_channel_complete",
                "V07_dual_support_only",
                "V08_dual_cross_model_no_local",
                "V09_dual_cross_primary_local",
                "V10_dual_cross_direct_local",
                "V11_full_v0_4",
            ]
        },
        "descriptive_maximum_not_for_model_selection": {
            "variant_id": descriptive_best.name,
            "coverage": descriptive_best["coverage"],
            "balanced_accuracy": descriptive_best["balanced_accuracy"],
            "mcc": descriptive_best["mcc"],
            "warning": "Outcome-observed post-hoc maximum; cannot be selected as a validated replacement model.",
        },
        "interpretation": protocol["interpretation_rules"],
        "output_hashes": {
            str(path.relative_to(PROJECT)): sha256(path)
            for path in [METRICS, SET_CONTRASTS, PAIRWISE_BOOTSTRAP]
        },
    }
    assert summary["full_v0_4_reconstruction"]["exactly_reproduced"]
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=lambda value: value.item() if isinstance(value, np.generic) else value) + "\n")

    print(json.dumps({
        "full_v0_4_reconstruction": summary["full_v0_4_reconstruction"],
        "identity_clean_key_variants": summary["identity_clean_key_variants"],
        "descriptive_maximum_not_for_model_selection": summary["descriptive_maximum_not_for_model_selection"],
        "outputs": [str(path.relative_to(PROJECT)) for path in [METRICS, SET_CONTRASTS, PAIRWISE_BOOTSTRAP, SUMMARY]],
    }, ensure_ascii=False, indent=2, default=lambda value: value.item() if isinstance(value, np.generic) else value))


if __name__ == "__main__":
    main()
