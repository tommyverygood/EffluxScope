#!/usr/bin/env python3
"""Analyze the locked ABCB1 v0.4 quantity-oriented external validation.

The script never refits or changes the frozen model. It collapses repeated
connectivity structures within each source, audits reference identity and
Murcko-scaffold overlap, and reports source-separated metrics. The Biogen
continuous MDR1-MDCK endpoint is the primary quantity gate; the heterogeneous
Daood labels are a retrospective corroboration layer only.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"

BIOGEN_PRED = PROCESSED / "abcb1_quantity_biogen_v0_4_predictions_v1.csv"
DAOOD_PRED = PROCESSED / "abcb1_quantity_daood_v0_4_predictions_v1.csv"
PRIMARY_REF = PROCESSED / "abcb1_screener_v0_1_training_audit.csv"
DIRECT_REF = PROCESSED / "soskuti_2024_mdr1_species_panel_v1.csv"
PROTOCOL = PROCESSED / "abcb1_quantity_external_validation_protocol_lock_v1.json"

BIOGEN_ROWS = PROCESSED / "abcb1_quantity_biogen_analysis_rows_v1.csv"
DAOOD_ROWS = PROCESSED / "abcb1_quantity_daood_analysis_rows_v1.csv"
METRICS = PROCESSED / "abcb1_quantity_external_metrics_v1.csv"
ABSTENTION = PROCESSED / "abcb1_quantity_external_abstention_reasons_v1.csv"
LEDGER = PROCESSED / "abcb1_quantity_external_evidence_ledger_v1.csv"
SUMMARY = PROCESSED / "abcb1_quantity_external_validation_summary_v1.json"

HIGH_CALL = "high_ABCB1_efflux_risk"
LOW_CALL = "low_ABCB1_efflux_risk"
ABSTAIN_CALL = "abstain"
SEED = 20260805
BOOTSTRAP_REPLICATES = 2000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value).strip()


def connectivity_block(inchikey: object) -> str:
    value = text(inchikey)
    if not value:
        return ""
    return value.split("-", 1)[0]


def murcko(smiles: object) -> str:
    value = text(smiles)
    if not value:
        return ""
    try:
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            return ""
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold.GetNumAtoms() == 0:
            return ""
        return Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=False)
    except Exception:
        return ""


def reference_sets() -> dict[str, set[str]]:
    primary = pd.read_csv(PRIMARY_REF, low_memory=False)
    direct = pd.read_csv(DIRECT_REF, low_memory=False)

    primary_keys = primary["standardized_inchikey"].fillna("").map(connectivity_block)
    derived_keys = primary["derived_standardized_inchikey"].fillna("").map(connectivity_block)
    primary_blocks = {x for x in pd.concat([primary_keys, derived_keys]) if x}
    direct_blocks = {
        x for x in direct["standardized_inchikey"].fillna("").map(connectivity_block) if x
    }

    primary_scaffolds = {
        text(x) for x in primary["murcko_scaffold_smiles"].fillna("") if text(x)
    }
    for value in primary["standardized_isomeric_smiles"].fillna(""):
        scaffold = murcko(value)
        if scaffold:
            primary_scaffolds.add(scaffold)
    direct_scaffolds = {
        scaffold
        for scaffold in direct["standardized_isomeric_smiles"].fillna("").map(murcko)
        if scaffold
    }
    return {
        "primary_blocks": primary_blocks,
        "direct_blocks": direct_blocks,
        "primary_scaffolds": primary_scaffolds,
        "direct_scaffolds": direct_scaffolds,
    }


def unique_join(values: pd.Series) -> str:
    return "|".join(sorted({text(value) for value in values if text(value)}))


def collapse_source(df: pd.DataFrame, source: str, refs: dict[str, set[str]]) -> pd.DataFrame:
    work = df.copy()
    work["connectivity_block"] = work["standardized_inchikey"].map(connectivity_block)
    work["murcko_scaffold_smiles"] = work["standardized_isomeric_smiles"].map(murcko)
    work["collapse_key"] = work["connectivity_block"]
    missing = work["collapse_key"].eq("")
    work.loc[missing, "collapse_key"] = "UNRESOLVED::" + work.loc[missing, "external_row_id"].astype(str)

    rows: list[pd.Series] = []
    endpoint_col = "log10_mdr1_mdck_er" if source == "Biogen_Fang_2023" else None
    label_col = "frozen_er2_binary_label" if source == "Biogen_Fang_2023" else "author_pgp_substrate_label"
    median_cols = [
        "primary_calibrated_probability",
        "bootstrap_probability_p10",
        "bootstrap_probability_median",
        "bootstrap_probability_p90",
        "bootstrap_interval_width",
        "direct_panel_probability",
        "primary_max_tanimoto",
        "direct_max_tanimoto",
        "physchem_max_abs_robust_z",
    ]
    if endpoint_col:
        median_cols += [endpoint_col, "mdr1_mdck_er"]

    for _, group in work.groupby("collapse_key", sort=False, dropna=False):
        base = group.iloc[0].copy()
        base["within_source_row_count"] = int(len(group))
        base["within_source_external_row_ids"] = unique_join(group["external_row_id"])
        base["within_source_record_ids"] = unique_join(group["source_record_id"])
        if "source_database" in group:
            base["source_database"] = unique_join(group["source_database"])

        calls = {text(value) for value in group["risk_call"] if text(value)}
        base["decision_conflict_within_connectivity"] = "yes" if len(calls) > 1 else "no"
        if len(calls) == 1:
            base["risk_call"] = next(iter(calls))
        elif len(calls) > 1:
            base["risk_call"] = ABSTAIN_CALL
            base["decision_zh"] = "拒判（连接体内决策冲突）"

        for column in median_cols:
            if column in group:
                numeric = pd.to_numeric(group[column], errors="coerce")
                base[column] = numeric.median() if numeric.notna().any() else np.nan

        labels = {
            int(value)
            for value in pd.to_numeric(group[label_col], errors="coerce").dropna().tolist()
        }
        base["endpoint_label_conflict_within_connectivity"] = "yes" if len(labels) > 1 else "no"
        if len(labels) == 1:
            base[label_col] = next(iter(labels))
        else:
            base[label_col] = np.nan

        if endpoint_col and pd.notna(base[endpoint_col]):
            base["frozen_er2_binary_label"] = int(float(base[endpoint_col]) >= math.log10(2.0))

        block = text(base["connectivity_block"])
        scaffold = text(base["murcko_scaffold_smiles"])
        usable = bool(block) and not text(base["standardization_status"]).lower().startswith("failed")
        overlap_primary = usable and block in refs["primary_blocks"]
        overlap_direct = usable and block in refs["direct_blocks"]
        scaffold_primary = usable and bool(scaffold) and scaffold in refs["primary_scaffolds"]
        scaffold_direct = usable and bool(scaffold) and scaffold in refs["direct_scaffolds"]

        base["model_usable"] = "yes" if usable else "no"
        base["exact_connectivity_overlap_primary_443"] = "yes" if overlap_primary else "no"
        base["exact_connectivity_overlap_direct_84"] = "yes" if overlap_direct else "no"
        base["identity_clean"] = "yes" if usable and not overlap_primary and not overlap_direct else "no"
        base["murcko_overlap_primary_443"] = "yes" if scaffold_primary else "no"
        base["murcko_overlap_direct_84"] = "yes" if scaffold_direct else "no"
        base["scaffold_clean"] = (
            "yes"
            if usable and not overlap_primary and not overlap_direct and not scaffold_primary and not scaffold_direct
            else "no"
        )
        base["bootstrap_cluster"] = scaffold if scaffold else ("ACYCLIC::" + block if block else text(base["external_row_id"]))
        base["analysis_source"] = source
        rows.append(base)

    result = pd.DataFrame(rows).drop(columns=["collapse_key"], errors="ignore")
    return result.reset_index(drop=True)


def safe_metric(function: Callable[[], float]) -> float:
    try:
        value = float(function())
        return value if np.isfinite(value) else np.nan
    except Exception:
        return np.nan


def accepted_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    if len(y_true) == 0:
        return {key: np.nan for key in [
            "accuracy", "balanced_accuracy", "mcc", "sensitivity", "specificity",
            "high_risk_precision", "low_risk_npv",
        ]}
    return {
        "accuracy": safe_metric(lambda: accuracy_score(y_true, y_pred)),
        "balanced_accuracy": safe_metric(lambda: balanced_accuracy_score(y_true, y_pred)),
        "mcc": safe_metric(lambda: matthews_corrcoef(y_true, y_pred)),
        "sensitivity": safe_metric(lambda: recall_score(y_true, y_pred, pos_label=1, zero_division=0)),
        "specificity": safe_metric(lambda: recall_score(y_true, y_pred, pos_label=0, zero_division=0)),
        "high_risk_precision": safe_metric(lambda: precision_score(y_true, y_pred, pos_label=1, zero_division=0)),
        "low_risk_npv": safe_metric(lambda: precision_score(y_true, y_pred, pos_label=0, zero_division=0)),
    }


def cluster_bootstrap(
    df: pd.DataFrame,
    metric: Callable[[pd.DataFrame], float],
    *,
    seed_offset: int,
) -> tuple[float, float, int]:
    if df.empty:
        return np.nan, np.nan, 0
    grouped = [indices.to_numpy() for _, indices in df.groupby("bootstrap_cluster").groups.items()]
    if len(grouped) < 2:
        return np.nan, np.nan, 0
    rng = np.random.default_rng(SEED + seed_offset)
    values: list[float] = []
    for _ in range(BOOTSTRAP_REPLICATES):
        selected = rng.integers(0, len(grouped), size=len(grouped))
        indices = np.concatenate([grouped[index] for index in selected])
        sample = df.loc[indices]
        value = safe_metric(lambda: metric(sample))
        if np.isfinite(value):
            values.append(value)
    if not values:
        return np.nan, np.nan, 0
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)), len(values)


def subset_frame(df: pd.DataFrame, subset: str) -> pd.DataFrame:
    if subset == "operational_all":
        return df.copy()
    if subset == "identity_clean":
        return df[df["identity_clean"].eq("yes")].copy()
    if subset == "scaffold_clean":
        return df[df["scaffold_clean"].eq("yes")].copy()
    raise ValueError(subset)


def metric_row(df: pd.DataFrame, source: str, subset: str, stratum: str = "all") -> dict[str, object]:
    work = subset_frame(df, subset)
    if stratum != "all" and "source_database" in work:
        work = work[work["source_database"].eq(stratum)].copy()
    label_col = "frozen_er2_binary_label" if source == "Biogen_Fang_2023" else "author_pgp_substrate_label"

    counts = work["risk_call"].value_counts().to_dict()
    accepted = work[work["risk_call"].isin([HIGH_CALL, LOW_CALL])].copy()
    accepted[label_col] = pd.to_numeric(accepted[label_col], errors="coerce")
    accepted = accepted[accepted[label_col].isin([0, 1])]
    accepted["predicted_binary"] = accepted["risk_call"].map({HIGH_CALL: 1, LOW_CALL: 0})
    y_true = accepted[label_col].astype(int).to_numpy()
    y_pred = accepted["predicted_binary"].astype(int).to_numpy()
    binary = accepted_metrics(y_true, y_pred)

    result: dict[str, object] = {
        "source": source,
        "subset": subset,
        "stratum": stratum,
        "n_unique_structures": int(len(work)),
        "n_endpoint_positive": int((pd.to_numeric(work[label_col], errors="coerce") == 1).sum()),
        "n_endpoint_negative": int((pd.to_numeric(work[label_col], errors="coerce") == 0).sum()),
        "n_high_call": int(counts.get(HIGH_CALL, 0)),
        "n_low_call": int(counts.get(LOW_CALL, 0)),
        "n_abstain": int(counts.get(ABSTAIN_CALL, 0)),
        "coverage": float(len(accepted) / len(work)) if len(work) else np.nan,
        "accepted_n": int(len(accepted)),
        "accepted_positive_n": int((accepted[label_col] == 1).sum()),
        "accepted_negative_n": int((accepted[label_col] == 0).sum()),
        **binary,
        "roc_auc_all_probabilities": np.nan,
        "spearman_probability_vs_endpoint": np.nan,
        "spearman_ci95_low": np.nan,
        "spearman_ci95_high": np.nan,
        "spearman_bootstrap_valid_replicates": 0,
        "balanced_accuracy_ci95_low": np.nan,
        "balanced_accuracy_ci95_high": np.nan,
        "balanced_accuracy_bootstrap_valid_replicates": 0,
        "mcc_ci95_low": np.nan,
        "mcc_ci95_high": np.nan,
        "mcc_bootstrap_valid_replicates": 0,
        "accepted_high_endpoint_median": np.nan,
        "accepted_low_endpoint_median": np.nan,
        "accepted_high_minus_low_median": np.nan,
        "accepted_high_vs_low_cliffs_delta": np.nan,
        "accepted_high_vs_low_mannwhitney_p": np.nan,
        "n_murcko_bootstrap_clusters": int(work["bootstrap_cluster"].nunique()) if len(work) else 0,
    }

    probability_rows = work.copy()
    probability_rows["primary_calibrated_probability"] = pd.to_numeric(
        probability_rows["primary_calibrated_probability"], errors="coerce"
    )
    probability_rows[label_col] = pd.to_numeric(probability_rows[label_col], errors="coerce")
    probability_rows = probability_rows.dropna(subset=["primary_calibrated_probability", label_col])
    if probability_rows[label_col].nunique() == 2:
        result["roc_auc_all_probabilities"] = safe_metric(
            lambda: roc_auc_score(probability_rows[label_col], probability_rows["primary_calibrated_probability"])
        )

    if source == "Biogen_Fang_2023":
        endpoint = "log10_mdr1_mdck_er"
        continuous = work.copy()
        continuous[endpoint] = pd.to_numeric(continuous[endpoint], errors="coerce")
        continuous["primary_calibrated_probability"] = pd.to_numeric(
            continuous["primary_calibrated_probability"], errors="coerce"
        )
        continuous = continuous.dropna(subset=[endpoint, "primary_calibrated_probability"])
        result["spearman_probability_vs_endpoint"] = safe_metric(
            lambda: spearmanr(continuous["primary_calibrated_probability"], continuous[endpoint]).statistic
        )
        low, high, valid = cluster_bootstrap(
            continuous,
            lambda sample: spearmanr(sample["primary_calibrated_probability"], sample[endpoint]).statistic,
            seed_offset=1000 + len(work),
        )
        result["spearman_ci95_low"] = low
        result["spearman_ci95_high"] = high
        result["spearman_bootstrap_valid_replicates"] = valid

        high_values = pd.to_numeric(
            work.loc[work["risk_call"].eq(HIGH_CALL), endpoint], errors="coerce"
        ).dropna().to_numpy()
        low_values = pd.to_numeric(
            work.loc[work["risk_call"].eq(LOW_CALL), endpoint], errors="coerce"
        ).dropna().to_numpy()
        if len(high_values) and len(low_values):
            test = mannwhitneyu(high_values, low_values, alternative="two-sided")
            result["accepted_high_endpoint_median"] = float(np.median(high_values))
            result["accepted_low_endpoint_median"] = float(np.median(low_values))
            result["accepted_high_minus_low_median"] = float(np.median(high_values) - np.median(low_values))
            result["accepted_high_vs_low_cliffs_delta"] = float(
                2.0 * test.statistic / (len(high_values) * len(low_values)) - 1.0
            )
            result["accepted_high_vs_low_mannwhitney_p"] = float(test.pvalue)

    if len(accepted) and accepted[label_col].nunique() == 2:
        low, high, valid = cluster_bootstrap(
            accepted,
            lambda sample: balanced_accuracy_score(sample[label_col].astype(int), sample["predicted_binary"].astype(int))
            if sample[label_col].nunique() == 2 else np.nan,
            seed_offset=2000 + len(work),
        )
        result["balanced_accuracy_ci95_low"] = low
        result["balanced_accuracy_ci95_high"] = high
        result["balanced_accuracy_bootstrap_valid_replicates"] = valid
        low, high, valid = cluster_bootstrap(
            accepted,
            lambda sample: matthews_corrcoef(sample[label_col].astype(int), sample["predicted_binary"].astype(int))
            if sample[label_col].nunique() == 2 else np.nan,
            seed_offset=3000 + len(work),
        )
        result["mcc_ci95_low"] = low
        result["mcc_ci95_high"] = high
        result["mcc_bootstrap_valid_replicates"] = valid
    return result


def abstention_rows(df: pd.DataFrame, source: str) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for subset in ["operational_all", "identity_clean", "scaffold_clean"]:
        work = subset_frame(df, subset)
        abstained = work[work["risk_call"].eq(ABSTAIN_CALL)]
        counts: Counter[str] = Counter()
        zh: dict[str, str] = {}
        for _, row in abstained.iterrows():
            try:
                reasons = json.loads(text(row.get("abstention_reasons_json")) or "[]")
            except json.JSONDecodeError:
                reasons = ["unparseable_reason_json"]
            try:
                reasons_zh = json.loads(text(row.get("abstention_reasons_zh_json")) or "[]")
            except json.JSONDecodeError:
                reasons_zh = []
            if not reasons:
                reasons = ["unspecified_or_collapsed_conflict"]
            for index, reason in enumerate(reasons):
                counts[reason] += 1
                if index < len(reasons_zh):
                    zh[reason] = reasons_zh[index]
        for reason, count in counts.most_common():
            output.append({
                "source": source,
                "subset": subset,
                "reason": reason,
                "reason_zh": zh.get(reason, ""),
                "n_structures": count,
                "fraction_of_all": count / len(work) if len(work) else np.nan,
                "fraction_of_abstentions": count / len(abstained) if len(abstained) else np.nan,
            })
    return output


def historical_ledger(new_metrics: pd.DataFrame) -> pd.DataFrame:
    def select(source: str, subset: str) -> pd.Series:
        return new_metrics[(new_metrics["source"] == source) & (new_metrics["subset"] == subset) & (new_metrics["stratum"] == "all")].iloc[0]

    biogen = select("Biogen_Fang_2023", "identity_clean")
    daood = select("Daood_2025", "identity_clean")

    asap = json.loads((PROCESSED / "abcb1_asap2025_sourceheldout_summary_v1.json").read_text())
    asap_gate = asap["gate_results"]["frozen_ABCB1_v0_1"]
    reg_posthoc = pd.read_csv(PROCESSED / "abcb1_selective_v0_4_regulatory_v4_posthoc_predictions_v1.csv")
    reg_accepted = reg_posthoc[reg_posthoc["risk_call"].isin([HIGH_CALL, LOW_CALL])].copy()
    reg_pred = reg_accepted["risk_call"].map({HIGH_CALL: 1, LOW_CALL: 0}).astype(int)
    reg_true = reg_accepted["directional_label"].astype(int)
    reg_binary = accepted_metrics(reg_true.to_numpy(), reg_pred.to_numpy())

    asap_posthoc = pd.read_csv(PROCESSED / "abcb1_selective_v0_4_asap2025_posthoc_predictions_v1.csv")
    rows = [
        {
            "source": "Biogen Fang 2023 public MDR1-MDCK",
            "evidence_role": "new locked primary quantity benchmark",
            "endpoint": "continuous log10 MDR1-MDCK efflux ratio",
            "n_total_or_estimable": int(biogen["n_unique_structures"]),
            "n_accepted_v0_4": int(biogen["accepted_n"]),
            "coverage_v0_4": biogen["coverage"],
            "primary_result": f"Spearman rho={biogen['spearman_probability_vs_endpoint']:.4f}; 95% scaffold CI [{biogen['spearman_ci95_low']:.4f}, {biogen['spearman_ci95_high']:.4f}]",
            "binary_or_selective_result": f"accepted BAcc={biogen['balanced_accuracy']:.4f}; MCC={biogen['mcc']:.4f}",
            "independence": "identity-clean against primary 443 and direct 84; endpoint has no matched parental or inhibitor rescue",
            "allowed_claim": "external ABCB1 efflux-risk ranking and selective triage",
        },
        {
            "source": "Daood 2025 curated P-gp classes",
            "evidence_role": "new secondary heterogeneous corroboration",
            "endpoint": "author-curated P-gp substrate class",
            "n_total_or_estimable": int(daood["n_unique_structures"]),
            "n_accepted_v0_4": int(daood["accepted_n"]),
            "coverage_v0_4": daood["coverage"],
            "primary_result": f"all-probability ROC-AUC={daood['roc_auc_all_probabilities']:.4f}",
            "binary_or_selective_result": f"accepted BAcc={daood['balanced_accuracy']:.4f}; MCC={daood['mcc']:.4f}",
            "independence": "identity-clean subset only; original source is heterogeneous and overlaps model-training provenance",
            "allowed_claim": "retrospective corroboration only",
        },
        {
            "source": "ASAP 2025 MDCK-MDR1",
            "evidence_role": "historical source-held-out continuous diagnostic",
            "endpoint": "log2 efflux ratio",
            "n_total_or_estimable": int(asap_gate["n"]),
            "n_accepted_v0_4": int(asap_posthoc["risk_call"].isin([HIGH_CALL, LOW_CALL]).sum()),
            "coverage_v0_4": float(asap_posthoc["risk_call"].isin([HIGH_CALL, LOW_CALL]).mean()),
            "primary_result": f"frozen primary-component rho={asap_gate['spearman_rho']:.4f}; 95% scaffold CI [{asap_gate['ci95_lower']:.4f}, {asap_gate['ci95_upper']:.4f}]",
            "binary_or_selective_result": "v0.4 accepted calls are coverage diagnostics, not a new tuned gate",
            "independence": "62 exact/connectivity-clean observations; no parental or inhibitor arm",
            "allowed_claim": "source-held-out continuous direction",
        },
        {
            "source": "Regulatory strict pool v4",
            "evidence_role": "historical direct-assay post-hoc v0.4 diagnostic",
            "endpoint": "evidence-tiered directional ABCB1 calls",
            "n_total_or_estimable": int(len(reg_posthoc)),
            "n_accepted_v0_4": int(len(reg_accepted)),
            "coverage_v0_4": float(len(reg_accepted) / len(reg_posthoc)),
            "primary_result": f"v0.4 accepted accuracy={reg_binary['accuracy']:.4f}",
            "binary_or_selective_result": f"accepted BAcc={reg_binary['balanced_accuracy']:.4f}; MCC={reg_binary['mcc']:.4f}",
            "independence": "direct/regulatory ordinary-drug evidence; v0.4 evaluation is post-hoc",
            "allowed_claim": "small direct-assay stress test only",
        },
    ]
    return pd.DataFrame(rows)


def round_frame(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    for column in result.select_dtypes(include=["number"]).columns:
        result[column] = result[column].round(8)
    return result


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text())
    assert protocol["status"] == "locked_before_new_external_outcome_scoring"
    assert protocol["bootstrap"]["replicates"] == BOOTSTRAP_REPLICATES

    refs = reference_sets()
    biogen_raw = pd.read_csv(BIOGEN_PRED, low_memory=False)
    daood_raw = pd.read_csv(DAOOD_PRED, low_memory=False)
    biogen = collapse_source(biogen_raw, "Biogen_Fang_2023", refs)
    daood = collapse_source(daood_raw, "Daood_2025", refs)

    metric_rows: list[dict[str, object]] = []
    for source, frame in [("Biogen_Fang_2023", biogen), ("Daood_2025", daood)]:
        for subset in ["operational_all", "identity_clean", "scaffold_clean"]:
            metric_rows.append(metric_row(frame, source, subset))
    for database in sorted({x for x in daood["source_database"].dropna().astype(str) if "|" not in x}):
        metric_rows.append(metric_row(daood, "Daood_2025", "identity_clean", database))
    metrics = round_frame(pd.DataFrame(metric_rows))

    abstentions = round_frame(pd.DataFrame(
        abstention_rows(biogen, "Biogen_Fang_2023") + abstention_rows(daood, "Daood_2025")
    ))
    ledger = historical_ledger(metrics)

    biogen_identity = metrics[
        (metrics["source"] == "Biogen_Fang_2023")
        & (metrics["subset"] == "identity_clean")
        & (metrics["stratum"] == "all")
    ].iloc[0]
    gate_direction = (
        float(biogen_identity["spearman_probability_vs_endpoint"]) > 0
        and float(biogen_identity["spearman_ci95_low"]) > 0
    )
    gate_selective = (
        float(biogen_identity["accepted_high_endpoint_median"])
        > float(biogen_identity["accepted_low_endpoint_median"])
    )
    quantity_gate_pass = bool(gate_direction and gate_selective)

    round_frame(biogen).to_csv(BIOGEN_ROWS, index=False)
    round_frame(daood).to_csv(DAOOD_ROWS, index=False)
    metrics.to_csv(METRICS, index=False)
    abstentions.to_csv(ABSTENTION, index=False)
    round_frame(ledger).to_csv(LEDGER, index=False)

    summary = {
        "analysis_version": "ABCB1_QUANTITY_EXTERNAL_VALIDATION_V1",
        "analysis_date": "2026-08-05",
        "status": "completed_frozen_model_no_tuning",
        "protocol_path": str(PROTOCOL.relative_to(PROJECT)),
        "protocol_sha256": sha256(PROTOCOL),
        "model_sha256": protocol["frozen_model"]["sha256"],
        "reference_panel_counts": {
            "primary_rows": 443,
            "primary_unique_connectivity_blocks": len(refs["primary_blocks"]),
            "direct_rows": 84,
            "direct_unique_connectivity_blocks": len(refs["direct_blocks"]),
        },
        "new_data_quantity": {
            "raw_scored_rows": int(len(biogen_raw) + len(daood_raw)),
            "biogen_raw_rows": int(len(biogen_raw)),
            "daood_raw_rows": int(len(daood_raw)),
            "collapsed_unique_connectivity_or_unresolved_rows": int(len(biogen) + len(daood)),
            "identity_clean_unique_rows": int(biogen["identity_clean"].eq("yes").sum() + daood["identity_clean"].eq("yes").sum()),
            "scaffold_clean_unique_rows": int(biogen["scaffold_clean"].eq("yes").sum() + daood["scaffold_clean"].eq("yes").sum()),
        },
        "overlap_audit": {
            "biogen_exact_overlap_primary": int(biogen["exact_connectivity_overlap_primary_443"].eq("yes").sum()),
            "biogen_exact_overlap_direct": int(biogen["exact_connectivity_overlap_direct_84"].eq("yes").sum()),
            "biogen_identity_clean": int(biogen["identity_clean"].eq("yes").sum()),
            "biogen_scaffold_clean": int(biogen["scaffold_clean"].eq("yes").sum()),
            "daood_exact_overlap_primary": int(daood["exact_connectivity_overlap_primary_443"].eq("yes").sum()),
            "daood_exact_overlap_direct": int(daood["exact_connectivity_overlap_direct_84"].eq("yes").sum()),
            "daood_identity_clean": int(daood["identity_clean"].eq("yes").sum()),
            "daood_scaffold_clean": int(daood["scaffold_clean"].eq("yes").sum()),
            "daood_label_conflict_connectivities": int(daood["endpoint_label_conflict_within_connectivity"].eq("yes").sum()),
        },
        "primary_biogen_identity_clean_result": {
            key: (None if pd.isna(value) else value.item() if isinstance(value, np.generic) else value)
            for key, value in biogen_identity.to_dict().items()
        },
        "predeclared_quantity_gate": {
            "continuous_direction_positive_ci_excludes_zero": gate_direction,
            "accepted_high_median_exceeds_low_median": gate_selective,
            "quantity_support_status": "pass" if quantity_gate_pass else "fail",
            "interpretation_if_pass": "supports a bounded external ABCB1 efflux-risk ranking/selective-triage claim",
            "interpretation_if_fail": "does not support an external quantity-based risk-ranking claim; retain only failure-boundary evidence",
        },
        "strict_predictor_upgrade": {
            "status": "not_permitted",
            "reason": "no new outcome-unseen matched-parental or inhibitor-rescue source with both positive and negative calls",
        },
        "claim_boundary": "free or released small-molecule ABCB1 efflux-risk ranking and selective triage only; not intact-ADC transport, universal substrate truth, intracellular unbound exposure, resistance, efficacy, or clinical DDI",
        "output_hashes": {},
    }
    for path in [BIOGEN_ROWS, DAOOD_ROWS, METRICS, ABSTENTION, LEDGER]:
        summary["output_hashes"][str(path.relative_to(PROJECT))] = sha256(path)
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")

    print(json.dumps({
        "new_data_quantity": summary["new_data_quantity"],
        "overlap_audit": summary["overlap_audit"],
        "primary_biogen_identity_clean_result": summary["primary_biogen_identity_clean_result"],
        "predeclared_quantity_gate": summary["predeclared_quantity_gate"],
        "outputs": [str(path.relative_to(PROJECT)) for path in [BIOGEN_ROWS, DAOOD_ROWS, METRICS, ABSTENTION, LEDGER, SUMMARY]],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
