#!/usr/bin/env python3
"""Build source-traceable tables for the ABCB1 selective-screener manuscript."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"
RESULTS_OUT = PROCESSED / "abcb1_manuscript_primary_results_v1.csv"
CLAIMS_OUT = PROCESSED / "abcb1_manuscript_claim_ledger_v1.csv"
FIGURES_OUT = PROCESSED / "abcb1_manuscript_figure_manifest_v1.csv"
TABLES_OUT = PROCESSED / "abcb1_manuscript_table_manifest_v1.csv"
READINESS_OUT = PROCESSED / "abcb1_manuscript_readiness_summary_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def select(table: pd.DataFrame, **conditions: str) -> pd.Series:
    mask = pd.Series(True, index=table.index)
    for column, value in conditions.items():
        mask &= table[column].eq(value)
    rows = table[mask]
    if len(rows) != 1:
        raise RuntimeError(f"expected one row for {conditions}, found {len(rows)}")
    return rows.iloc[0]


def main() -> None:
    quantity = pd.read_csv(PROCESSED / "abcb1_quantity_external_metrics_v1.csv")
    binary = pd.read_csv(PROCESSED / "abcb1_weak_binary_v0_5_external_metrics_v1.csv")
    ledger = pd.read_csv(PROCESSED / "abcb1_quantity_external_evidence_ledger_v1.csv")
    qsummary = json.loads((PROCESSED / "abcb1_quantity_external_validation_summary_v1.json").read_text())
    risk_operating = pd.read_csv(PROCESSED / "abcb1_external_risk_operating_points_v1.csv")
    risk_enrichment = pd.read_csv(PROCESSED / "abcb1_external_selective_enrichment_bootstrap_v1.csv")
    risk_calibration = pd.read_csv(PROCESSED / "abcb1_external_calibration_summary_v1.csv")
    risk_loo = pd.read_csv(PROCESSED / "abcb1_external_rejector_leave_one_out_v1.csv")
    channel_metrics = pd.read_csv(PROCESSED / "abcb1_external_evidence_channel_ablation_metrics_v1.csv")
    channel_pairwise = pd.read_csv(PROCESSED / "abcb1_external_evidence_channel_pairwise_bootstrap_v1.csv")
    calibration_v2 = pd.read_csv(PROCESSED / "abcb1_external_calibration_metrics_v2.csv")
    calibration_bootstrap_v2 = pd.read_csv(PROCESSED / "abcb1_external_calibration_bootstrap_v2.csv")
    calibration_pairwise_v2 = pd.read_csv(PROCESSED / "abcb1_external_calibration_pairwise_bootstrap_v2.csv")
    error_counts_v2 = pd.read_csv(PROCESSED / "abcb1_external_error_counts_v2.csv")
    error_strata_v2 = pd.read_csv(PROCESSED / "abcb1_external_error_strata_v2.csv")

    binary_identity = select(
        binary,
        source="Biogen_Fang_2023",
        subset="identity_clean",
        support_layer="binary_all_in_scope",
    )
    selective_identity = select(
        binary,
        source="Biogen_Fang_2023",
        subset="identity_clean",
        support_layer="v0_4_supported_only",
    )
    forced_identity = select(
        binary,
        source="Biogen_Fang_2023",
        subset="identity_clean",
        support_layer="forced_low_confidence_only",
    )
    binary_scaffold = select(
        binary,
        source="Biogen_Fang_2023",
        subset="scaffold_clean",
        support_layer="binary_all_in_scope",
    )
    selective_scaffold = select(
        binary,
        source="Biogen_Fang_2023",
        subset="scaffold_clean",
        support_layer="v0_4_supported_only",
    )
    continuous_identity = select(
        quantity,
        source="Biogen_Fang_2023",
        subset="identity_clean",
        stratum="all",
    )
    continuous_scaffold = select(
        quantity,
        source="Biogen_Fang_2023",
        subset="scaffold_clean",
        stratum="all",
    )
    daood_identity = select(
        binary,
        source="Daood_2025",
        subset="identity_clean",
        support_layer="binary_all_in_scope",
    )
    asap = ledger[ledger["source"].eq("ASAP 2025 MDCK-MDR1")].iloc[0]
    regulatory = ledger[ledger["source"].eq("Regulatory strict pool v4")].iloc[0]
    margin_identity = select(risk_operating, population="identity_clean", label="frozen_probability_margin_only")
    enrichment_bacc = select(risk_enrichment, population="identity_clean", metric="balanced_accuracy")
    enrichment_mcc = select(risk_enrichment, population="identity_clean", metric="mcc")
    calibration_identity = select(risk_calibration, population="identity_clean_all")
    without_local = select(risk_loo, population="identity_clean", label="without_local_neighborhood_agreement")
    without_model_agreement = select(risk_loo, population="identity_clean", label="without_supported_model_direction_agreement")
    primary_channel = select(channel_metrics, population="identity_clean", variant_id="V04_primary_channel_complete")
    direct_channel = select(channel_metrics, population="identity_clean", variant_id="V06_direct_channel_complete")
    direct_no_local = select(channel_metrics, population="identity_clean", variant_id="V05_direct_support_no_local")
    full_vs_primary_channel_bacc = select(
        channel_pairwise,
        population="identity_clean",
        comparator_variant="V04_primary_channel_complete",
        metric="balanced_accuracy",
    )
    full_vs_primary_channel_mcc = select(
        channel_pairwise,
        population="identity_clean",
        comparator_variant="V04_primary_channel_complete",
        metric="mcc",
    )
    raw_calibration_v2 = select(calibration_v2, population="identity_clean", method="raw_frozen_score")
    platt_calibration_v2 = select(calibration_v2, population="identity_clean", method="cross_fitted_platt")
    raw_slope_bootstrap_v2 = select(
        calibration_bootstrap_v2,
        population="identity_clean",
        method="raw_frozen_score",
        metric="calibration_slope",
    )
    platt_brier_pairwise_v2 = select(
        calibration_pairwise_v2,
        population="identity_clean",
        calibrated_method="cross_fitted_platt",
        metric="brier_score",
    )
    false_positive_count_v2 = select(error_counts_v2, population="identity_clean", error_outcome="false_positive")
    mid_high_band_v2 = select(
        error_strata_v2,
        population="identity_clean",
        stratum_type="probability_band",
        stratum_level="[0.50,0.70)",
    )

    results = pd.DataFrame([
        {
            "result_id": "R1_BINARY_BASELINE_IDENTITY_CLEAN",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak binary sensitivity label",
            "analysis_subset": "identity-clean",
            "method_layer": "forced full-coverage weak binary baseline",
            "n_total": int(binary_identity["n_subset_total"]),
            "n_evaluated": int(binary_identity["n_evaluated"]),
            "coverage": binary_identity["coverage_within_subset"],
            "primary_metric": "balanced_accuracy",
            "primary_value": binary_identity["balanced_accuracy"],
            "ci95_low": binary_identity["balanced_accuracy_ci95_low"],
            "ci95_high": binary_identity["balanced_accuracy_ci95_high"],
            "mcc": binary_identity["mcc"],
            "mcc_ci95_low": binary_identity["mcc_ci95_low"],
            "mcc_ci95_high": binary_identity["mcc_ci95_high"],
            "sensitivity": binary_identity["sensitivity"],
            "specificity": binary_identity["specificity"],
            "ppv": binary_identity["positive_predictive_value"],
            "npv": binary_identity["negative_predictive_value"],
            "manuscript_role": "baseline",
            "evidence_timing": "threshold pre-existing; interface evaluation post-outcome",
            "allowed_interpretation": "full coverage has a stable but low-discrimination signal",
        },
        {
            "result_id": "R2_SELECTIVE_PRIMARY_IDENTITY_CLEAN",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak binary sensitivity label",
            "analysis_subset": "identity-clean",
            "method_layer": "v0.4 high/low accepted; abstentions excluded from accepted-performance denominator",
            "n_total": int(selective_identity["n_subset_total"]),
            "n_evaluated": int(selective_identity["n_evaluated"]),
            "coverage": selective_identity["coverage_within_subset"],
            "primary_metric": "balanced_accuracy",
            "primary_value": selective_identity["balanced_accuracy"],
            "ci95_low": selective_identity["balanced_accuracy_ci95_low"],
            "ci95_high": selective_identity["balanced_accuracy_ci95_high"],
            "mcc": selective_identity["mcc"],
            "mcc_ci95_low": selective_identity["mcc_ci95_low"],
            "mcc_ci95_high": selective_identity["mcc_ci95_high"],
            "sensitivity": selective_identity["sensitivity"],
            "specificity": selective_identity["specificity"],
            "ppv": selective_identity["positive_predictive_value"],
            "npv": selective_identity["negative_predictive_value"],
            "manuscript_role": "primary selective-classification result",
            "evidence_timing": "v0.4 rules and thresholds locked before Biogen outcomes were scored",
            "allowed_interpretation": "abstention identifies a smaller subset with materially improved reliability",
        },
        {
            "result_id": "R3_FORCED_REJECTED_LAYER_IDENTITY_CLEAN",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak binary sensitivity label",
            "analysis_subset": "identity-clean",
            "method_layer": "rows rejected by v0.4 but forced through the 0.5 binary interface",
            "n_total": int(forced_identity["n_subset_total"]),
            "n_evaluated": int(forced_identity["n_evaluated"]),
            "coverage": forced_identity["coverage_within_subset"],
            "primary_metric": "balanced_accuracy",
            "primary_value": forced_identity["balanced_accuracy"],
            "ci95_low": forced_identity["balanced_accuracy_ci95_low"],
            "ci95_high": forced_identity["balanced_accuracy_ci95_high"],
            "mcc": forced_identity["mcc"],
            "mcc_ci95_low": forced_identity["mcc_ci95_low"],
            "mcc_ci95_high": forced_identity["mcc_ci95_high"],
            "sensitivity": forced_identity["sensitivity"],
            "specificity": forced_identity["specificity"],
            "ppv": forced_identity["positive_predictive_value"],
            "npv": forced_identity["negative_predictive_value"],
            "manuscript_role": "ablation and falsification control",
            "evidence_timing": "post-outcome binary overlay; no threshold tuning",
            "allowed_interpretation": "most coverage gained by forcing abstentions is low-confidence",
        },
        {
            "result_id": "R4_BINARY_BASELINE_SCAFFOLD_CLEAN",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak binary sensitivity label",
            "analysis_subset": "identity- and reference-Murcko-scaffold-clean",
            "method_layer": "forced full-coverage weak binary baseline",
            "n_total": int(binary_scaffold["n_subset_total"]),
            "n_evaluated": int(binary_scaffold["n_evaluated"]),
            "coverage": binary_scaffold["coverage_within_subset"],
            "primary_metric": "balanced_accuracy",
            "primary_value": binary_scaffold["balanced_accuracy"],
            "ci95_low": binary_scaffold["balanced_accuracy_ci95_low"],
            "ci95_high": binary_scaffold["balanced_accuracy_ci95_high"],
            "mcc": binary_scaffold["mcc"],
            "mcc_ci95_low": binary_scaffold["mcc_ci95_low"],
            "mcc_ci95_high": binary_scaffold["mcc_ci95_high"],
            "sensitivity": binary_scaffold["sensitivity"],
            "specificity": binary_scaffold["specificity"],
            "ppv": binary_scaffold["positive_predictive_value"],
            "npv": binary_scaffold["negative_predictive_value"],
            "manuscript_role": "scaffold extrapolation sensitivity",
            "evidence_timing": "same as R1",
            "allowed_interpretation": "weak baseline direction is not explained by known reference scaffolds alone",
        },
        {
            "result_id": "R5_SELECTIVE_PRIMARY_SCAFFOLD_CLEAN",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak binary sensitivity label",
            "analysis_subset": "identity- and reference-Murcko-scaffold-clean",
            "method_layer": "v0.4 high/low accepted",
            "n_total": int(selective_scaffold["n_subset_total"]),
            "n_evaluated": int(selective_scaffold["n_evaluated"]),
            "coverage": selective_scaffold["coverage_within_subset"],
            "primary_metric": "balanced_accuracy",
            "primary_value": selective_scaffold["balanced_accuracy"],
            "ci95_low": selective_scaffold["balanced_accuracy_ci95_low"],
            "ci95_high": selective_scaffold["balanced_accuracy_ci95_high"],
            "mcc": selective_scaffold["mcc"],
            "mcc_ci95_low": selective_scaffold["mcc_ci95_low"],
            "mcc_ci95_high": selective_scaffold["mcc_ci95_high"],
            "sensitivity": selective_scaffold["sensitivity"],
            "specificity": selective_scaffold["specificity"],
            "ppv": selective_scaffold["positive_predictive_value"],
            "npv": selective_scaffold["negative_predictive_value"],
            "manuscript_role": "primary robustness result",
            "evidence_timing": "same as R2",
            "allowed_interpretation": "selective improvement persists beyond reference Murcko scaffolds",
        },
        {
            "result_id": "R6_CONTINUOUS_DIRECTION_IDENTITY_CLEAN",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "continuous log10 efflux ratio",
            "analysis_subset": "identity-clean",
            "method_layer": "frozen primary probability ranking",
            "n_total": int(continuous_identity["n_unique_structures"]),
            "n_evaluated": int(continuous_identity["n_unique_structures"]),
            "coverage": 1.0,
            "primary_metric": "Spearman_rho",
            "primary_value": continuous_identity["spearman_probability_vs_endpoint"],
            "ci95_low": continuous_identity["spearman_ci95_low"],
            "ci95_high": continuous_identity["spearman_ci95_high"],
            "mcc": "",
            "mcc_ci95_low": "",
            "mcc_ci95_high": "",
            "sensitivity": "",
            "specificity": "",
            "ppv": "",
            "npv": "",
            "manuscript_role": "primary continuous endpoint result",
            "evidence_timing": "protocol and model locked before new outcomes were scored",
            "allowed_interpretation": "large external dataset supports a positive ABCB1 efflux-risk direction",
        },
        {
            "result_id": "R7_CONTINUOUS_DIRECTION_SCAFFOLD_CLEAN",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "continuous log10 efflux ratio",
            "analysis_subset": "identity- and reference-Murcko-scaffold-clean",
            "method_layer": "frozen primary probability ranking",
            "n_total": int(continuous_scaffold["n_unique_structures"]),
            "n_evaluated": int(continuous_scaffold["n_unique_structures"]),
            "coverage": 1.0,
            "primary_metric": "Spearman_rho",
            "primary_value": continuous_scaffold["spearman_probability_vs_endpoint"],
            "ci95_low": continuous_scaffold["spearman_ci95_low"],
            "ci95_high": continuous_scaffold["spearman_ci95_high"],
            "mcc": "",
            "mcc_ci95_low": "",
            "mcc_ci95_high": "",
            "sensitivity": "",
            "specificity": "",
            "ppv": "",
            "npv": "",
            "manuscript_role": "continuous scaffold-robustness result",
            "evidence_timing": "same as R6",
            "allowed_interpretation": "continuous direction persists beyond known reference scaffolds",
        },
        {
            "result_id": "R8_DAOOD_LIMITATION",
            "dataset": "Daood 2025 curated P-gp labels",
            "endpoint": "heterogeneous author-curated substrate class",
            "analysis_subset": "identity-clean; only 4 positive labels",
            "method_layer": "forced full-coverage weak binary baseline",
            "n_total": int(daood_identity["n_subset_total"]),
            "n_evaluated": int(daood_identity["n_evaluated"]),
            "coverage": daood_identity["coverage_within_subset"],
            "primary_metric": "balanced_accuracy",
            "primary_value": daood_identity["balanced_accuracy"],
            "ci95_low": daood_identity["balanced_accuracy_ci95_low"],
            "ci95_high": daood_identity["balanced_accuracy_ci95_high"],
            "mcc": daood_identity["mcc"],
            "mcc_ci95_low": daood_identity["mcc_ci95_low"],
            "mcc_ci95_high": daood_identity["mcc_ci95_high"],
            "sensitivity": daood_identity["sensitivity"],
            "specificity": daood_identity["specificity"],
            "ppv": daood_identity["positive_predictive_value"],
            "npv": daood_identity["negative_predictive_value"],
            "manuscript_role": "limitation and source-leakage audit",
            "evidence_timing": "retrospective heterogeneous-source corroboration",
            "allowed_interpretation": "large nominal datasets can collapse after identity and provenance cleaning",
        },
        {
            "result_id": "R9_ASAP_SOURCE_HELDOUT",
            "dataset": "ASAP 2025 MDCK-MDR1",
            "endpoint": "continuous log2 efflux ratio",
            "analysis_subset": "62 estimable exact/connectivity-clean observations",
            "method_layer": "frozen primary probability ranking",
            "n_total": int(asap["n_total_or_estimable"]),
            "n_evaluated": int(asap["n_total_or_estimable"]),
            "coverage": 1.0,
            "primary_metric": "Spearman_rho",
            "primary_value": 0.55015301,
            "ci95_low": 0.20837074,
            "ci95_high": 0.74409192,
            "mcc": "",
            "mcc_ci95_low": "",
            "mcc_ci95_high": "",
            "sensitivity": "",
            "specificity": "",
            "ppv": "",
            "npv": "",
            "manuscript_role": "independent source-held-out continuous support",
            "evidence_timing": "pre-registered source-held-out scoring",
            "allowed_interpretation": "positive continuous direction replicated in a second MDCK-MDR1 source",
        },
        {
            "result_id": "R10_REGULATORY_POSTHOC",
            "dataset": "Regulatory strict pool v4",
            "endpoint": "evidence-tiered direct ABCB1 directional calls",
            "analysis_subset": "33 ordinary-drug entities",
            "method_layer": "v0.4 accepted calls",
            "n_total": int(regulatory["n_total_or_estimable"]),
            "n_evaluated": int(regulatory["n_accepted_v0_4"]),
            "coverage": regulatory["coverage_v0_4"],
            "primary_metric": "accepted_accuracy",
            "primary_value": 0.84615385,
            "ci95_low": "",
            "ci95_high": "",
            "mcc": 0.63960215,
            "mcc_ci95_low": "",
            "mcc_ci95_high": "",
            "sensitivity": "",
            "specificity": "",
            "ppv": "",
            "npv": "",
            "manuscript_role": "small direct-assay post-hoc stress test",
            "evidence_timing": "v0.4 post-hoc; outcomes observed before v0.4 rule freeze",
            "allowed_interpretation": "direct evidence is directionally compatible but not prospective confirmation",
        },
        {
            "result_id": "R11_PROBABILITY_MARGIN_ONLY_POSTHOC",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak external risk label",
            "analysis_subset": "identity-clean",
            "method_layer": "frozen probability margin only; p<=0.35 or p>=0.70",
            "n_total": int(margin_identity["denominator_n"]),
            "n_evaluated": int(margin_identity["accepted_n"]),
            "coverage": margin_identity["coverage"],
            "primary_metric": "balanced_accuracy",
            "primary_value": margin_identity["balanced_accuracy"],
            "ci95_low": margin_identity["balanced_accuracy_ci95_low"],
            "ci95_high": margin_identity["balanced_accuracy_ci95_high"],
            "mcc": margin_identity["mcc"],
            "mcc_ci95_low": margin_identity["mcc_ci95_low"],
            "mcc_ci95_high": margin_identity["mcc_ci95_high"],
            "sensitivity": margin_identity["sensitivity"],
            "specificity": margin_identity["specificity"],
            "ppv": margin_identity["high_risk_precision"],
            "npv": margin_identity["low_risk_npv"],
            "manuscript_role": "post-hoc confidence-only comparator",
            "evidence_timing": "secondary analysis locked after primary external outcome observation; no threshold selection",
            "allowed_interpretation": "probability margin alone does not explain the full selective gain",
        },
        {
            "result_id": "R12_SELECTIVE_ENRICHMENT_DIFFERENCE_POSTHOC",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak external risk label",
            "analysis_subset": "identity-clean accepted versus rejected-forced layers",
            "method_layer": "joint Murcko-scaffold bootstrap contrast",
            "n_total": int(binary_identity["n_subset_total"]),
            "n_evaluated": int(binary_identity["n_evaluated"]),
            "coverage": 1.0,
            "primary_metric": "accepted_minus_rejected_balanced_accuracy",
            "primary_value": enrichment_bacc["favorable_difference"],
            "ci95_low": enrichment_bacc["ci95_low"],
            "ci95_high": enrichment_bacc["ci95_high"],
            "mcc": enrichment_mcc["favorable_difference"],
            "mcc_ci95_low": enrichment_mcc["ci95_low"],
            "mcc_ci95_high": enrichment_mcc["ci95_high"],
            "sensitivity": "",
            "specificity": "",
            "ppv": "",
            "npv": "",
            "manuscript_role": "post-hoc selective-enrichment contrast",
            "evidence_timing": "secondary analysis locked after primary external outcome observation",
            "allowed_interpretation": "the frozen rejector enriches lower-reliability rows under the weak external endpoint",
        },
        {
            "result_id": "R13_EXTERNAL_CALIBRATION_POSTHOC",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak external risk label",
            "analysis_subset": "identity-clean",
            "method_layer": "frozen primary score calibration diagnostic",
            "n_total": int(calibration_identity["n"]),
            "n_evaluated": int(calibration_identity["n"]),
            "coverage": 1.0,
            "primary_metric": "calibration_slope",
            "primary_value": calibration_identity["calibration_slope"],
            "ci95_low": "",
            "ci95_high": "",
            "mcc": "",
            "mcc_ci95_low": "",
            "mcc_ci95_high": "",
            "sensitivity": "",
            "specificity": "",
            "ppv": "",
            "npv": "",
            "manuscript_role": "post-hoc calibration limitation",
            "evidence_timing": "secondary analysis locked after primary external outcome observation",
            "allowed_interpretation": f"score is not an absolute probability; ECE={calibration_identity['ece_equal_width_10']:.3f}",
        },
        {
            "result_id": "R14_LOCAL_RULE_ABLATION_HYPOTHESIS",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak external risk label",
            "analysis_subset": "identity-clean",
            "method_layer": "full frozen rejector without local-neighborhood agreement rule",
            "n_total": int(without_local["denominator_n"]),
            "n_evaluated": int(without_local["accepted_n"]),
            "coverage": without_local["coverage"],
            "primary_metric": "balanced_accuracy",
            "primary_value": without_local["balanced_accuracy"],
            "ci95_low": without_local["balanced_accuracy_ci95_low"],
            "ci95_high": without_local["balanced_accuracy_ci95_high"],
            "mcc": without_local["mcc"],
            "mcc_ci95_low": without_local["mcc_ci95_low"],
            "mcc_ci95_high": without_local["mcc_ci95_high"],
            "sensitivity": without_local["sensitivity"],
            "specificity": without_local["specificity"],
            "ppv": without_local["high_risk_precision"],
            "npv": without_local["low_risk_npv"],
            "manuscript_role": "post-hoc rule-simplification hypothesis",
            "evidence_timing": "same external outcome used to identify this hypothesis; requires new frozen validation",
            "allowed_interpretation": "local-neighborhood guard was not supported here; do not promote this ablation as a validated replacement model",
        },
        {
            "result_id": "R15_DIRECT_EVIDENCE_CHANNEL_CONTRIBUTION",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak external risk label",
            "analysis_subset": "identity-clean full v0.4 versus primary-reference-only complete channel",
            "method_layer": "post-hoc evidence-channel ablation with joint Murcko-scaffold bootstrap",
            "n_total": int(primary_channel["denominator_n"]),
            "n_evaluated": int(primary_channel["accepted_n"]),
            "coverage": primary_channel["coverage"],
            "primary_metric": "full_minus_primary_channel_balanced_accuracy",
            "primary_value": full_vs_primary_channel_bacc["full_minus_comparator"],
            "ci95_low": full_vs_primary_channel_bacc["ci95_low"],
            "ci95_high": full_vs_primary_channel_bacc["ci95_high"],
            "mcc": full_vs_primary_channel_mcc["full_minus_comparator"],
            "mcc_ci95_low": full_vs_primary_channel_mcc["ci95_low"],
            "mcc_ci95_high": full_vs_primary_channel_mcc["ci95_high"],
            "sensitivity": "",
            "specificity": "",
            "ppv": "",
            "npv": "",
            "manuscript_role": "post-hoc evidence-channel decomposition",
            "evidence_timing": "component variants locked after primary and first rule-ablation outcomes were observed",
            "allowed_interpretation": "direct-experiment reference evidence contributes non-redundant selective information; not a causal or prospective channel comparison",
        },
        {
            "result_id": "R16_RAW_EXTERNAL_CALIBRATION_LIMITATION",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak external risk label",
            "analysis_subset": "identity-clean",
            "method_layer": "frozen raw score external calibration",
            "n_total": int(raw_calibration_v2["n"]),
            "n_evaluated": int(raw_calibration_v2["n"]),
            "coverage": 1.0,
            "primary_metric": "calibration_slope",
            "primary_value": raw_calibration_v2["calibration_slope"],
            "ci95_low": raw_slope_bootstrap_v2["ci95_low"],
            "ci95_high": raw_slope_bootstrap_v2["ci95_high"],
            "mcc": "",
            "mcc_ci95_low": "",
            "mcc_ci95_high": "",
            "sensitivity": "",
            "specificity": "",
            "ppv": "",
            "npv": "",
            "manuscript_role": "post-hoc external calibration limitation",
            "evidence_timing": "deep calibration protocol locked after prior external diagnostics were observed",
            "allowed_interpretation": f"raw score is over-scaled and overpredicts high risk; ECE={raw_calibration_v2['ece_equal_width_10']:.3f}; do not interpret as universal probability",
        },
        {
            "result_id": "R17_WITHIN_SOURCE_PLATT_DIAGNOSTIC",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak external risk label",
            "analysis_subset": "identity-clean 5-fold Murcko-group cross-fit",
            "method_layer": "same-source Platt recalibration diagnostic",
            "n_total": int(platt_calibration_v2["n"]),
            "n_evaluated": int(platt_calibration_v2["n"]),
            "coverage": 1.0,
            "primary_metric": "favorable_brier_improvement",
            "primary_value": platt_brier_pairwise_v2["favorable_improvement"],
            "ci95_low": platt_brier_pairwise_v2["ci95_low"],
            "ci95_high": platt_brier_pairwise_v2["ci95_high"],
            "mcc": "",
            "mcc_ci95_low": "",
            "mcc_ci95_high": "",
            "sensitivity": "",
            "specificity": "",
            "ppv": "",
            "npv": "",
            "manuscript_role": "post-hoc same-source calibratability diagnostic",
            "evidence_timing": "all folds come from the same outcome-observed Biogen source",
            "allowed_interpretation": f"assay-specific scale is correctable within source; Platt slope={platt_calibration_v2['calibration_slope']:.3f}; not independent validation",
        },
        {
            "result_id": "R18_FALSE_POSITIVE_ERROR_CONCENTRATION",
            "dataset": "Biogen/Fang 2023 MDR1-MDCK",
            "endpoint": "ER>=2 weak external risk label",
            "analysis_subset": "identity-clean forced binary",
            "method_layer": "post-hoc error taxonomy and frozen score bands",
            "n_total": int(raw_calibration_v2["n"]),
            "n_evaluated": int(mid_high_band_v2["n"]),
            "coverage": float(mid_high_band_v2["n"] / raw_calibration_v2["n"]),
            "primary_metric": "error_rate_in_score_0_50_to_0_70",
            "primary_value": mid_high_band_v2["error_rate"],
            "ci95_low": mid_high_band_v2["error_rate_wilson95_low"],
            "ci95_high": mid_high_band_v2["error_rate_wilson95_high"],
            "mcc": "",
            "mcc_ci95_low": "",
            "mcc_ci95_high": "",
            "sensitivity": "",
            "specificity": "",
            "ppv": "",
            "npv": "",
            "manuscript_role": "post-hoc error localization",
            "evidence_timing": "deep error protocol locked after prior external diagnostics were observed",
            "allowed_interpretation": f"forced-binary errors are concentrated just above 0.5 and are dominated by {int(false_positive_count_v2['n'])} false positives",
        },
    ])
    results.to_csv(RESULTS_OUT, index=False)

    claims = pd.DataFrame([
        ["C1", "primary", "A frozen ordinary-drug prior shows a reproducible positive association with continuous MDR1-MDCK efflux risk.", "Biogen identity-clean rho=0.247 CI 0.207-0.284; ASAP rho=0.550 CI 0.208-0.744", "Do not call this accurate substrate truth or causal transport attribution.", "allowed"],
        ["C2", "primary", "Uncertainty-aware abstention improves reliability at the cost of coverage.", "Biogen binary baseline BAcc=0.582/MCC=0.168 at 100% coverage versus selective BAcc=0.735/MCC=0.499 at 17.3% coverage", "Always report coverage and class-specific accepted counts.", "allowed"],
        ["C3", "primary", "The improvement persists after exact identity and reference-Murcko-scaffold exclusion.", "Scaffold-clean selective n=439; BAcc=0.723; MCC=0.479; continuous rho=0.242 CI 0.205-0.281", "Do not equate scaffold-clean with prospective source independence.", "allowed"],
        ["C4", "primary", "Forced decisions among rows rejected by the selective system are substantially less reliable.", "Forced layer n=2171; BAcc=0.547; MCC=0.096", "Do not hide this layer when presenting a full-coverage binary interface.", "allowed"],
        ["C5", "secondary", "Abstention is explainable through probability gray zones, domain limits, bootstrap instability, and model/neighborhood conflicts.", "abcb1_quantity_external_abstention_reasons_v1.csv", "Abstention is a decision policy, not a biological third class.", "allowed"],
        ["C6", "secondary", "Nominally large public binary data can become unsuitable for external validation after provenance and identity cleaning.", "Daood 607 rows -> 156 identity-clean rows with only 4 positives", "Do not report Daood clean ROC-AUC without the four-positive denominator.", "allowed"],
        ["C11", "secondary", "The selective gain is not explained by probability extremeness alone and is driven mainly by supported cross-model directional agreement.", f"Margin-only BAcc={margin_identity['balanced_accuracy']:.3f}; full v0.4 BAcc=0.735; without model agreement BAcc={without_model_agreement['balanced_accuracy']:.3f}", "This is a post-hoc external rule decomposition, not a newly tuned operating policy.", "allowed_with_boundary"],
        ["C12", "secondary", "The frozen output is an external risk-ranking score rather than an absolute substrate probability.", f"Calibration slope={calibration_identity['calibration_slope']:.3f}; ECE={calibration_identity['ece_equal_width_10']:.3f}", "Do not translate a score of 0.70 into a 70% probability of substrate truth.", "allowed"],
        ["C13", "secondary", "The local-neighborhood guard did not improve this external benchmark and should be re-tested in a future frozen model revision.", f"Full BAcc=0.735/MCC=0.499 versus without-local BAcc={without_local['balanced_accuracy']:.3f}/MCC={without_local['mcc']:.3f}", "Do not call the post-hoc without-local variant externally validated.", "allowed_with_boundary"],
        ["C14", "secondary", "The direct-experiment reference channel contributes non-redundant selective information beyond the 443-row weak-label reference channel.", f"Full versus primary-channel-only BAcc difference={full_vs_primary_channel_bacc['full_minus_comparator']:.3f}, 95% CI {full_vs_primary_channel_bacc['ci95_low']:.3f}-{full_vs_primary_channel_bacc['ci95_high']:.3f}; direct-complete BAcc={direct_channel['balanced_accuracy']:.3f}", "Channel variants have different coverage and accepted sets; describe as post-hoc decomposition, not a causal or prospectively validated model comparison.", "allowed_with_boundary"],
        ["C15", "secondary", "The direct-channel local-neighborhood guard accounts for most of the local-rule performance loss in this benchmark.", f"Direct-support no-local BAcc={direct_no_local['balanced_accuracy']:.3f}; direct-channel complete BAcc={direct_channel['balanced_accuracy']:.3f}; full v0.4 BAcc=0.735", "This identifies a future rule-revision hypothesis only; do not replace v0.4 using the same external outcome.", "allowed_with_boundary"],
        ["C16", "secondary", "The frozen raw output is significantly miscalibrated in the Biogen MDR1-MDCK system.", f"Calibration intercept={raw_calibration_v2['calibration_intercept']:.3f}; slope={raw_calibration_v2['calibration_slope']:.3f}, 95% CI {raw_slope_bootstrap_v2['ci95_low']:.3f}-{raw_slope_bootstrap_v2['ci95_high']:.3f}; ECE={raw_calibration_v2['ece_equal_width_10']:.3f}", "Call the output a risk score, not an absolute universal substrate probability.", "allowed"],
        ["C17", "secondary", "A simple Platt mapping can correct probability scale within the same Biogen source without improving discrimination.", f"Cross-fitted Brier {raw_calibration_v2['brier_score']:.3f}->{platt_calibration_v2['brier_score']:.3f}; ECE {raw_calibration_v2['ece_equal_width_10']:.3f}->{platt_calibration_v2['ece_equal_width_10']:.3f}; ROC-AUC {raw_calibration_v2['roc_auc']:.3f}->{platt_calibration_v2['roc_auc']:.3f}", "This is same-source group-cross-fitting, not an independently validated calibrator and not a replacement for v0.4.", "allowed_with_boundary"],
        ["C18", "secondary", "Forced-binary error is dominated by false positives and concentrated in the score interval immediately above 0.5.", f"FP={int(false_positive_count_v2['n'])}; score [0.50,0.70) n={int(mid_high_band_v2['n'])}, error rate={mid_high_band_v2['error_rate']:.3f}", "Do not tune a new threshold on this same outcome; use the finding to justify abstention and future frozen validation.", "allowed"],
        ["C7", "application", "The framework can triage standardized free or ADC-released small-molecule structures.", "v0.4 ADC release-species demonstration and context guards", "Application demonstration is not ADC-specific accuracy validation.", "allowed_with_boundary"],
        ["C8", "forbidden", "This is a universal assay-independent ABCB1 substrate classifier.", "No qualifying matched-parental/inhibitor-rescue future source with balanced positives and negatives", "Replace with selective ABCB1 efflux-risk screening framework.", "forbidden"],
        ["C9", "forbidden", "The study validates intact-ADC efflux, intracellular unbound exposure, resistance, efficacy, or clinical DDI.", "No such endpoint is modeled", "Keep ADC as an implications/application section only.", "forbidden"],
        ["C10", "forbidden", "Abstain is a third biological class.", "v0.4 protocol defines abstention as insufficient/conflicting/out-of-domain evidence", "Use high risk, low risk, or abstain decision outputs; two biological classes only when discussing the weak baseline.", "forbidden"],
    ], columns=["claim_id", "role", "manuscript_wording", "supporting_evidence", "mandatory_boundary_or_replacement", "status"])
    claims.to_csv(CLAIMS_OUT, index=False)

    figures = pd.DataFrame([
        ["Figure 1", "Study design, evidence layers, and leakage controls", "Flow diagram: primary 443 and direct 84 references; Biogen/Daood/ASAP/regulatory external roles; identity and scaffold filters", "protocol locks; source registry; analysis-row files", "required", "ready_for_drawing"],
        ["Figure 2", "Coverage-reliability trade-off", "Confidence-only retention curves plus frozen full-binary, margin-only, selective-accepted, and forced-rejected operating points", "outputs/20260805_abcb1_external_risk_experiments_v1/figures/Figure2_external_coverage_reliability.*", "primary", "figure_ready"],
        ["Figure 3", "External rejector ablation", "Cumulative rejection-rule path, leave-one-rule-out comparison, and overlapping rule burden", "outputs/20260805_abcb1_external_risk_experiments_v1/figures/Figure3_external_rejector_ablation.*", "primary", "figure_ready"],
        ["Figure 4", "External calibration and error-layer diagnostics", "Calibration curves, probability versus continuous log10 ER, and endpoint distributions by selective layer", "outputs/20260805_abcb1_external_risk_experiments_v1/figures/Figure4_external_calibration_error.*", "primary", "figure_ready"],
        ["Figure 5", "Cross-source evidence and application boundary", "Source-separated evidence forest/table for Biogen, ASAP, regulatory, and Daood; ADC release species shown as application examples only", "abcb1_quantity_external_evidence_ledger_v1.csv; ADC registry predictions", "secondary/application", "data_ready"],
        ["Figure S1", "Exact identity and Murcko-scaffold overlap audit", "Counts before/after duplicate, identity, and scaffold filtering", "quantity validation summary and analysis rows", "supplementary", "data_ready"],
        ["Figure S2", "Daood class collapse after leakage cleaning", "607 nominal rows, 443 primary overlaps, 156 clean with four positives", "Daood analysis rows", "supplementary", "data_ready"],
        ["Figure S3", "Probability and endpoint distributions by decision output", "Violin/ECDF and high-versus-low Cliff delta", "Biogen analysis rows and quantity metrics", "supplementary", "data_ready"],
        ["Figure S4", "External evidence-channel ablation", "Coverage-reliability, channel-level BAcc, sensitivity/specificity, and primary/direct local-guard decomposition", "outputs/20260805_abcb1_external_evidence_channel_ablation_v1/figures/FigureS4_external_evidence_channel_ablation.*", "supplementary", "figure_ready"],
        ["Figure S5", "Deep external calibration and forced-decision errors", "Raw and cross-fitted calibration curves, calibration metrics, TP/TN/FP/FN counts, and error by frozen score band", "outputs/20260805_abcb1_external_calibration_error_v2/figures/FigureS5_external_calibration_error_v2.*", "supplementary", "figure_ready"],
        ["Figure S6", "External error phenotypes and scaffold hotspots", "FP/FN feature effect sizes, error by selection/support layer, and shrunken Murcko-scaffold hotspots", "outputs/20260805_abcb1_external_calibration_error_v2/figures/FigureS6_external_error_phenotypes_v2.*", "supplementary", "figure_ready"],
    ], columns=["figure_id", "title", "content", "source_artifacts", "role", "status"])
    figures.to_csv(FIGURES_OUT, index=False)

    tables = pd.DataFrame([
        ["Table 1", "Dataset roles, endpoints, overlap, and allowed use", "source registry; quantity summary", "main", "data_ready"],
        ["Table 2", "Binary baseline versus selective screening", "R1-R5", "main", "data_ready"],
        ["Table 3", "Source-separated continuous and direct evidence", "R6-R10", "main", "data_ready"],
        ["Table S1", "Frozen model and protocol hashes", "protocol/config/manifest files", "supplementary", "data_ready"],
        ["Table S2", "Per-structure Biogen predictions and overlap audit", "abcb1_quantity_biogen_analysis_rows_v1.csv", "supplementary", "data_ready"],
        ["Table S3", "Per-structure Daood audit", "abcb1_quantity_daood_analysis_rows_v1.csv", "supplementary", "data_ready"],
        ["Table S4", "Abstention reason counts", "abcb1_quantity_external_abstention_reasons_v1.csv", "supplementary", "data_ready"],
        ["Table S5", "Claim ledger and prohibited interpretations", "abcb1_manuscript_claim_ledger_v1.csv", "supplementary/internal governance", "data_ready"],
        ["Table S6", "External coverage-reliability curve and rejector ablation", "abcb1_external_risk_coverage_curve_v1.csv; abcb1_external_rejector_cumulative_ablation_v1.csv; abcb1_external_rejector_leave_one_out_v1.csv", "supplementary", "data_ready"],
        ["Table S7", "External calibration and error strata", "abcb1_external_calibration_summary_v1.csv; abcb1_external_calibration_bins_v1.csv; abcb1_external_error_strata_v1.csv", "supplementary", "data_ready"],
        ["Table S8", "External evidence-channel ablation and accepted-set contrasts", "abcb1_external_evidence_channel_ablation_metrics_v1.csv; abcb1_external_evidence_channel_set_contrasts_v1.csv; abcb1_external_evidence_channel_pairwise_bootstrap_v1.csv", "supplementary", "data_ready"],
        ["Table S9", "Deep external calibration and error diagnostics", "abcb1_external_calibration_metrics_v2.csv; abcb1_external_calibration_bootstrap_v2.csv; abcb1_external_error_counts_v2.csv; abcb1_external_error_feature_comparisons_v2.csv; abcb1_external_error_strata_v2.csv", "supplementary", "data_ready"],
    ], columns=["table_id", "title", "source", "role", "status"])
    tables.to_csv(TABLES_OUT, index=False)

    readiness = {
        "package_version": "ABCB1_SELECTIVE_MANUSCRIPT_PACKAGE_V1",
        "build_date": "2026-08-05",
        "recommended_article_type": "methods_and_data_resource_article_with_large_external_evaluation",
        "recommended_primary_title": "An uncertainty-aware selective screening framework for ABCB1 efflux risk: large-scale external evaluation and implications for ADC-released small molecules",
        "central_thesis": "weak-label full-coverage binary classification has reproducible but low discrimination; an uncertainty-aware abstention policy identifies a smaller, more reliable subset and exposes where evidence is insufficient or conflicting",
        "main_method": "v0.4 selective high-risk/low-risk/abstain screening",
        "baseline": "v0.5 weak-label forced binary interface at frozen 0.5 threshold",
        "gates": {
            "large_external_continuous_direction": "pass",
            "identity_clean_analysis": "pass",
            "scaffold_clean_sensitivity": "pass",
            "selective_reliability_improvement": "pass",
            "posthoc_external_risk_experiments": "pass",
            "posthoc_evidence_channel_ablation": "pass",
            "posthoc_deep_calibration_error_analysis": "pass",
            "binary_full_coverage_high_discrimination": "fail",
            "future_strict_matched_parental_or_inhibitor_validation": "pending",
            "adc_specific_external_validation": "fail",
        },
        "manuscript_readiness": "go_for_draft_under_bounded_methods_claim",
        "not_ready_for": [
            "universal ABCB1 substrate classifier claim",
            "ADC-specific predictive accuracy claim",
            "intact-ADC transport, intracellular exposure, resistance, efficacy, or clinical DDI claim",
        ],
        "primary_numeric_anchors": {
            "biogen_identity_clean_n": int(continuous_identity["n_unique_structures"]),
            "continuous_spearman": continuous_identity["spearman_probability_vs_endpoint"],
            "continuous_spearman_ci95": [continuous_identity["spearman_ci95_low"], continuous_identity["spearman_ci95_high"]],
            "binary_baseline_bacc": binary_identity["balanced_accuracy"],
            "binary_baseline_mcc": binary_identity["mcc"],
            "selective_coverage": selective_identity["coverage_within_subset"],
            "selective_bacc": selective_identity["balanced_accuracy"],
            "selective_mcc": selective_identity["mcc"],
            "forced_rejected_bacc": forced_identity["balanced_accuracy"],
            "forced_rejected_mcc": forced_identity["mcc"],
            "accepted_minus_rejected_bacc": enrichment_bacc["favorable_difference"],
            "accepted_minus_rejected_bacc_ci95": [enrichment_bacc["ci95_low"], enrichment_bacc["ci95_high"]],
            "external_calibration_slope": calibration_identity["calibration_slope"],
            "external_ece_10bin": calibration_identity["ece_equal_width_10"],
            "without_local_posthoc_bacc": without_local["balanced_accuracy"],
            "full_minus_primary_channel_bacc": full_vs_primary_channel_bacc["full_minus_comparator"],
            "full_minus_primary_channel_bacc_ci95": [full_vs_primary_channel_bacc["ci95_low"], full_vs_primary_channel_bacc["ci95_high"]],
            "direct_channel_complete_bacc": direct_channel["balanced_accuracy"],
            "direct_support_no_local_bacc": direct_no_local["balanced_accuracy"],
            "raw_external_calibration_slope": raw_calibration_v2["calibration_slope"],
            "raw_external_calibration_slope_ci95": [raw_slope_bootstrap_v2["ci95_low"], raw_slope_bootstrap_v2["ci95_high"]],
            "crossfitted_platt_brier": platt_calibration_v2["brier_score"],
            "false_positive_n": int(false_positive_count_v2["n"]),
            "score_0_50_to_0_70_error_rate": mid_high_band_v2["error_rate"],
            "quantity_gate": qsummary["predeclared_quantity_gate"]["quantity_support_status"],
        },
        "outputs": {},
    }
    for path in [RESULTS_OUT, CLAIMS_OUT, FIGURES_OUT, TABLES_OUT]:
        readiness["outputs"][str(path.relative_to(PROJECT))] = sha256(path)
    READINESS_OUT.write_text(json.dumps(readiness, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({
        "manuscript_readiness": readiness["manuscript_readiness"],
        "gates": readiness["gates"],
        "outputs": [str(path.relative_to(PROJECT)) for path in [RESULTS_OUT, CLAIMS_OUT, FIGURES_OUT, TABLES_OUT, READINESS_OUT]],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
