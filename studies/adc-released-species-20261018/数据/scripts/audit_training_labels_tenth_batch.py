#!/usr/bin/env python3
"""Audit provenance and source-label confounding in the Daood substrate set.

This script never changes Daood's published binary labels.  It assigns an
evidence/provenance tier and a model-use role, then measures whether source-
held-out evaluation is identifiable.  Tier names are project audit terms, not
claims made by the source authors or regulatory standards.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold


DATA = Path(__file__).resolve().parents[1]
PROCESSED = DATA / "processed"
RAW_LONG = PROCESSED / "daood_2025_pgp_bcrp_substrate_training_long.csv"
DECISIONS = PROCESSED / "daood_training_leakage_decisions_v1.csv"

RUN_DATE = "2026-07-28"
SEED = 20260728
DAOOD_DOI = "10.1021/acs.molpharmaceut.5c01065"
DAOOD_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC12587445/"
MORGAN_RADIUS = 2
MORGAN_BITS = 2048


SOURCE_RULES = {
    "ChEMBL": {
        "source_origin": "Daood author-curated ChEMBL records",
        "source_selection_role": "major source; authors manually reviewed assay context and primary literature",
        "provenance_tier": "B_author_curated_no_row_primary_assay",
        "tier_rank": 2,
        "default_model_role": "AB_candidate_primary_supervised",
        "row_level_primary_assay_id_available": "no",
        "row_level_primary_reference_available": "no",
        "row_level_numeric_transport_available": "no",
        "database_url": "https://www.ebi.ac.uk/chembl/",
    },
    "PubChem": {
        "source_origin": "Daood author-curated PubChem records",
        "source_selection_role": "supplementary public bioassay source reviewed by authors",
        "provenance_tier": "B_author_curated_no_row_primary_assay",
        "tier_rank": 2,
        "default_model_role": "AB_candidate_primary_supervised",
        "row_level_primary_assay_id_available": "no",
        "row_level_primary_reference_available": "no",
        "row_level_numeric_transport_available": "no",
        "database_url": "https://pubchem.ncbi.nlm.nih.gov/",
    },
    "Metrabase": {
        "source_origin": "Metrabase labels imported by Daood",
        "source_selection_role": "added to supplement transporter substrate/non-substrate class balance",
        "provenance_tier": "C_database_compilation_no_row_primary_assay",
        "tier_rank": 3,
        "default_model_role": "auxiliary_or_PU_only",
        "row_level_primary_assay_id_available": "no",
        "row_level_primary_reference_available": "no",
        "row_level_numeric_transport_available": "no",
        "database_url": "https://www-metrabase.ch.cam.ac.uk/metrabaseui/",
    },
    "Shaikh et al": {
        "source_origin": "Shaikh et al. 2017 compiled data imported by Daood",
        "source_selection_role": "added specifically to increase BCRP nonsubstrates",
        "provenance_tier": "C_secondary_compilation_no_row_primary_assay",
        "tier_rank": 3,
        "default_model_role": "auxiliary_or_PU_only",
        "row_level_primary_assay_id_available": "no",
        "row_level_primary_reference_available": "no",
        "row_level_numeric_transport_available": "no",
        "database_url": "https://doi.org/10.1021/acs.jcim.6b00508",
    },
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if fields is None:
        if not rows:
            raise ValueError(f"No rows and no explicit fields for {path}")
        fields = list(rows[0])
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def source_label_matrix(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    counts: dict[tuple[str, str], Counter] = defaultdict(Counter)
    for row in rows:
        counts[(row["transporter_gene"], row["source_databases"])][row["consensus_binary_label"]] += 1

    result = []
    for (gene, source), labels in sorted(counts.items()):
        total = labels["0"] + labels["1"]
        result.append(
            {
                "transporter_gene": gene,
                "source_database": source,
                "n_total": total,
                "n_negative": labels["0"],
                "n_positive": labels["1"],
                "positive_fraction": round(labels["1"] / total, 6) if total else "",
                "source_class_coverage": "single_class" if not labels["0"] or not labels["1"] else "both_classes",
                "source_label_confounding_flag": "yes" if not labels["0"] or not labels["1"] else "no",
                "daood_source_selection_note": SOURCE_RULES.get(source, {}).get("source_selection_role", "unmapped source"),
            }
        )
    return result


def cramers_v(rows: list[dict[str, str]]) -> float:
    sources = sorted({row["source_databases"] for row in rows})
    table = np.array(
        [
            [sum(row["source_databases"] == source and row["consensus_binary_label"] == label for row in rows) for label in ("0", "1")]
            for source in sources
        ],
        dtype=float,
    )
    n = table.sum()
    expected = table.sum(axis=1, keepdims=True) @ table.sum(axis=0, keepdims=True) / n
    valid = expected > 0
    chi2 = float((((table - expected) ** 2 / np.where(valid, expected, 1))[valid]).sum())
    denom = n * min(table.shape[0] - 1, table.shape[1] - 1)
    return math.sqrt(chi2 / denom) if denom else math.nan


def audit_rows(decisions: list[dict[str, str]], raw_by_id: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    matrix_lookup = {
        (row["transporter_gene"], row["source_database"]): row
        for row in source_label_matrix(decisions)
    }
    audited = []
    for row in decisions:
        source = row["source_databases"]
        rule = SOURCE_RULES.get(
            source,
            {
                "source_origin": "unmapped",
                "source_selection_role": "unmapped",
                "provenance_tier": "U_unmapped_source",
                "tier_rank": 9,
                "default_model_role": "exclude_until_mapped",
                "row_level_primary_assay_id_available": "no",
                "row_level_primary_reference_available": "no",
                "row_level_numeric_transport_available": "no",
                "database_url": "",
            },
        )
        raw_ids = row["source_training_record_ids"].split(" | ")
        raw = raw_by_id.get(raw_ids[0], {})
        matrix = matrix_lookup[(row["transporter_gene"], source)]
        qc_flags = []
        if row["tautomer_limit_flag"] == "yes":
            qc_flags.append("tautomer_enumeration_limit")
        if row["source_record_count"] == "1":
            qc_flags.append("no_cross_source_replication")
        if matrix["source_label_confounding_flag"] == "yes":
            qc_flags.append("single_class_source_for_transporter")
        if row["leakage_model_action"] != "include_leakage_clean_training":
            qc_flags.append(row["leakage_model_action"])

        role = rule["default_model_role"]
        strict = "no"
        if row["leakage_model_action"] != "include_leakage_clean_training":
            role = "exclude_external_evaluation_leakage"
        elif rule["tier_rank"] <= 2:
            strict = "yes"
        elif rule["tier_rank"] == 3:
            role = "auxiliary_or_PU_only"

        negative_basis = "not_applicable_positive_label"
        if row["consensus_binary_label"] == "0":
            if source == "ChEMBL":
                negative_basis = "Daood_author_curated_nonsubstrate_call"
            elif source == "Metrabase":
                negative_basis = "Metrabase_database_nonsubstrate_compilation"
            elif source == "Shaikh et al":
                negative_basis = "Shaikh_secondary_nonsubstrate_compilation"
            else:
                negative_basis = "source_negative_without_row_numeric_provenance"

        audited.append(
            {
                "audit_record_id": f"T10_AUD_{len(audited) + 1:05d}",
                "consensus_record_id": row["consensus_record_id"],
                "training_structure_id": row["training_structure_id"],
                "transporter_gene": row["transporter_gene"],
                "transporter_protein": row["transporter_protein"],
                "published_binary_label": row["consensus_binary_label"],
                "preferred_names": row["preferred_names"],
                "source_pids": row["source_pids"],
                "source_database": source,
                "source_record_count": int(row["source_record_count"]),
                "source_origin": rule["source_origin"],
                "source_selection_role": rule["source_selection_role"],
                "source_has_both_labels_for_transporter": "no" if matrix["source_label_confounding_flag"] == "yes" else "yes",
                "source_n_negative": matrix["n_negative"],
                "source_n_positive": matrix["n_positive"],
                "source_label_confounding_flag": matrix["source_label_confounding_flag"],
                "row_level_primary_assay_id_available": rule["row_level_primary_assay_id_available"],
                "row_level_primary_reference_available": rule["row_level_primary_reference_available"],
                "row_level_numeric_transport_available": rule["row_level_numeric_transport_available"],
                "direct_bidirectional_transport_confirmed": "no",
                "inhibitor_rescue_confirmed": "no",
                "negative_evidence_basis": negative_basis,
                "provenance_tier": rule["provenance_tier"],
                "strict_AB_candidate": strict,
                "model_use_role": role,
                "manual_audit_status": "source_method_audited_row_primary_assay_unavailable",
                "chemical_qc_flags": " | ".join(qc_flags),
                "leakage_model_action": row["leakage_model_action"],
                "standardized_inchikey": row["standardized_inchikey"],
                "standardized_isomeric_smiles": row["standardized_isomeric_smiles"],
                "murcko_scaffold_smiles": row["murcko_scaffold_smiles"],
                "std_molecular_weight": row["std_molecular_weight"],
                "std_xlogp": row["std_xlogp"],
                "std_tpsa": row["std_tpsa"],
                "original_endpoint_field": raw.get("endpoint", ""),
                "original_label_source": raw.get("label_source", ""),
                "daood_source_doi": DAOOD_DOI,
                "daood_source_url": DAOOD_URL,
                "upstream_database_url": rule["database_url"],
                "audit_note": "Published label retained unchanged; tier describes project provenance completeness, not label correctness.",
            }
        )
    return audited


def fp_matrix(rows: list[dict[str, object]]) -> np.ndarray:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=MORGAN_RADIUS, fpSize=MORGAN_BITS, includeChirality=False)
    output = np.zeros((len(rows), MORGAN_BITS), dtype=np.uint8)
    for index, row in enumerate(rows):
        mol = Chem.MolFromSmiles(str(row["standardized_isomeric_smiles"]))
        if mol is None:
            raise ValueError(f"Could not parse audited structure {row['consensus_record_id']}")
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(mol), output[index])
    return output


def evaluate_predictions(y_true: np.ndarray, scores: np.ndarray) -> tuple[float, float, float]:
    predictions = (scores >= 0.5).astype(int)
    return (
        float(matthews_corrcoef(y_true, predictions)),
        float(balanced_accuracy_score(y_true, predictions)),
        float(roc_auc_score(y_true, scores)),
    )


def source_holdout_diagnostics(audited: list[dict[str, object]]) -> list[dict[str, object]]:
    results = []
    for gene in ("ABCB1", "ABCG2"):
        rows = [
            row for row in audited
            if row["transporter_gene"] == gene and row["leakage_model_action"] == "include_leakage_clean_training"
        ]
        x = fp_matrix(rows)
        y = np.array([int(row["published_binary_label"]) for row in rows], dtype=int)
        sources = np.array([row["source_database"] for row in rows])
        for held_source in sorted(set(sources.tolist())):
            train = sources != held_source
            test = sources == held_source
            train_classes = sorted(set(y[train].tolist()))
            test_classes = sorted(set(y[test].tolist()))
            base = {
                "transporter_gene": gene,
                "held_out_source": held_source,
                "n_train": int(train.sum()),
                "train_n_negative": int((y[train] == 0).sum()),
                "train_n_positive": int((y[train] == 1).sum()),
                "n_test": int(test.sum()),
                "test_n_negative": int((y[test] == 0).sum()),
                "test_n_positive": int((y[test] == 1).sum()),
                "model": "Morgan_r2_2048_logistic_C0.1_balanced",
                "mcc": "",
                "balanced_accuracy": "",
                "roc_auc": "",
            }
            if len(train_classes) < 2:
                base["metric_status"] = "not_identifiable_training_single_class"
            elif len(test_classes) < 2:
                model = LogisticRegression(C=0.1, class_weight="balanced", max_iter=4000, solver="liblinear", random_state=SEED)
                model.fit(x[train], y[train])
                scores = model.predict_proba(x[test])[:, 1]
                predictions = (scores >= 0.5).astype(int)
                if test_classes == [0]:
                    base["single_class_accuracy"] = round(float((predictions == 0).mean()), 6)
                    base["single_class_metric"] = "specificity_on_all_negative_holdout"
                else:
                    base["single_class_accuracy"] = round(float((predictions == 1).mean()), 6)
                    base["single_class_metric"] = "sensitivity_on_all_positive_holdout"
                base["metric_status"] = "test_single_class_no_MCC_BAcc_AUC"
            else:
                model = LogisticRegression(C=0.1, class_weight="balanced", max_iter=4000, solver="liblinear", random_state=SEED)
                model.fit(x[train], y[train])
                scores = model.predict_proba(x[test])[:, 1]
                mcc, bacc, auc = evaluate_predictions(y[test], scores)
                base.update(
                    {
                        "metric_status": "valid_binary_source_holdout",
                        "mcc": round(mcc, 6),
                        "balanced_accuracy": round(bacc, 6),
                        "roc_auc": round(auc, 6),
                        "single_class_accuracy": "",
                        "single_class_metric": "",
                    }
                )
            base.setdefault("single_class_accuracy", "")
            base.setdefault("single_class_metric", "")
            results.append(base)
    return results


def tier_b_scaffold_cv(audited: list[dict[str, object]]) -> list[dict[str, object]]:
    results = []
    for gene in ("ABCB1", "ABCG2"):
        rows = [
            row for row in audited
            if row["transporter_gene"] == gene
            and row["strict_AB_candidate"] == "yes"
            and row["leakage_model_action"] == "include_leakage_clean_training"
        ]
        x = fp_matrix(rows)
        y = np.array([int(row["published_binary_label"]) for row in rows], dtype=int)
        groups = np.array(
            [row["murcko_scaffold_smiles"] or f"ACYCLIC::{row['training_structure_id']}" for row in rows]
        )
        splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED)
        fold_metrics = []
        invalid = 0
        for train, test in splitter.split(x, y, groups):
            if len(set(y[train].tolist())) < 2 or len(set(y[test].tolist())) < 2:
                invalid += 1
                continue
            model = LogisticRegression(C=0.1, class_weight="balanced", max_iter=4000, solver="liblinear", random_state=SEED)
            model.fit(x[train], y[train])
            fold_metrics.append(evaluate_predictions(y[test], model.predict_proba(x[test])[:, 1]))
        if fold_metrics:
            values = np.array(fold_metrics)
            mean_values = values.mean(axis=0)
            sd_values = values.std(axis=0)
        else:
            mean_values = [math.nan, math.nan, math.nan]
            sd_values = [math.nan, math.nan, math.nan]
        results.append(
            {
                "transporter_gene": gene,
                "training_set": "tier_AB_candidate_author_curated_rows",
                "n_records": len(rows),
                "n_negative": int((y == 0).sum()),
                "n_positive": int((y == 1).sum()),
                "source_databases": " | ".join(sorted({str(row["source_database"]) for row in rows})),
                "model": "Morgan_r2_2048_logistic_C0.1_balanced",
                "split": "5-fold_StratifiedGroupKFold_exact_Murcko",
                "valid_folds": len(fold_metrics),
                "invalid_single_class_folds": invalid,
                "mean_mcc": round(float(mean_values[0]), 6),
                "sd_mcc": round(float(sd_values[0]), 6),
                "mean_balanced_accuracy": round(float(mean_values[1]), 6),
                "sd_balanced_accuracy": round(float(sd_values[1]), 6),
                "mean_roc_auc": round(float(mean_values[2]), 6),
                "sd_roc_auc": round(float(sd_values[2]), 6),
                "gate_rule": "MCC >= 0.50 AND balanced_accuracy >= 0.75",
                "gate_status": "pass" if mean_values[0] >= 0.5 and mean_values[1] >= 0.75 and invalid == 0 else "fail",
            }
        )
    return results


def issue_registry(audited: list[dict[str, object]], matrix: list[dict[str, object]]) -> list[dict[str, object]]:
    issues = []
    def add(issue_type: str, scope: str, evidence: str, impact: str, action: str, severity: str) -> None:
        issues.append(
            {
                "issue_id": f"T10_ISSUE_{len(issues) + 1:03d}",
                "issue_type": issue_type,
                "scope": scope,
                "evidence": evidence,
                "model_impact": impact,
                "required_action": action,
                "severity": severity,
                "status": "open",
            }
        )

    for row in matrix:
        if row["source_label_confounding_flag"] == "yes":
            add(
                "single_class_source",
                f"{row['transporter_gene']} / {row['source_database']}",
                f"negative={row['n_negative']}; positive={row['n_positive']}",
                "source-held-out binary validation is undefined or label-source selection can be encoded through chemical-space differences",
                "obtain both classes from an independent source or move this source to auxiliary/PU use",
                "critical",
            )
    add(
        "row_level_assay_provenance_missing",
        "all 1114 Daood substrate rows",
        "Table S2 exposes PID, source database, structure and binary endpoints but not per-row assay ID, primary citation, numeric result or conditions",
        "cannot independently distinguish direct transport from other author-curated evidence at row level",
        "recover assay-level ChEMBL/Metrabase/Shaikh provenance for prioritized rows",
        "critical",
    )
    cross_source = sum(int(row["source_record_count"]) > 1 for row in audited)
    add(
        "no_cross_source_consensus_records",
        "ABCB1 and ABCG2 substrate training set",
        f"records with source_record_count > 1: {cross_source}",
        "source-independent label agreement cannot be estimated",
        "structure-crosswalk independent databases and retain agreements/conflicts",
        "high",
    )
    tier_a = sum(str(row["provenance_tier"]).startswith("A_") for row in audited)
    add(
        "no_tier_A_direct_rows",
        "ordinary substrate training truth set",
        f"tier A rows: {tier_a}",
        "the strict set is an author-curated Tier B candidate rather than a direct-assay truth set",
        "collect bidirectional Papp/net ER and inhibitor-rescue rows with primary provenance",
        "high",
    )
    leakage = Counter(row["leakage_model_action"] for row in audited)
    add(
        "external_reference_leakage_exclusions",
        "ADC/Roth external evaluation",
        "; ".join(f"{key}={value}" for key, value in sorted(leakage.items()) if key != "include_leakage_clean_training"),
        "removing related training examples reduces optimistic ADC external performance",
        "retain existing exclusions and never tune them using external outcomes",
        "high",
    )
    tautomer = sum("tautomer_enumeration_limit" in str(row["chemical_qc_flags"]) for row in audited)
    add(
        "tautomer_enumeration_limit",
        "chemical standardization",
        f"flagged records: {tautomer}",
        "large/flexible molecules may have unstable canonical representations",
        "run sensitivity analysis with source and connectivity representations",
        "medium",
    )
    return issues


def gap_queue() -> list[dict[str, object]]:
    return [
        {
            "priority": "P0",
            "gap_id": "T10_GAP_001",
            "transporter_gene": "ABCB1",
            "missing_data": "independent source containing confirmed substrates and nonsubstrates with row-level assay provenance",
            "why_needed": "holding out ChEMBL leaves a training set with no positives",
            "target_sources": "ChEMBL assay records | TransPortal | ISTransbase | primary literature | regulatory reviews",
            "acceptance_criteria": "exact structure; human ABCB1; direct substrate endpoint; experimental system and numeric result or explicit author call; primary citation",
            "reject_if": "inhibitor-only; unspecified species; database label without primary trace; unlabeled treated as negative",
            "model_use": "enable identifiable source-held-out validation",
        },
        {
            "priority": "P0",
            "gap_id": "T10_GAP_002",
            "transporter_gene": "ABCG2",
            "missing_data": "additional independently sourced direct nonsubstrates and weak substrates",
            "why_needed": "Shaikh contributes only negatives and ChEMBL has only 15 clean negatives",
            "target_sources": "ChEMBL assay records | TransPortal | ISTransbase | primary literature | FDA/EMA/PMDA",
            "acceptance_criteria": "human ABCG2/BCRP direct transport or vesicle endpoint with controls and primary citation",
            "reject_if": "ABCB1-only inference; cytotoxicity without transporter-isogenic control; missing exact identity",
            "model_use": "reduce source-label confounding and improve specificity",
        },
        {
            "priority": "P0",
            "gap_id": "T10_GAP_003",
            "transporter_gene": "ABCB1 | ABCG2",
            "missing_data": "cross-source exact-structure overlaps",
            "why_needed": "all current consensus records have source_record_count=1",
            "target_sources": "Metrabase | TransPortal | ISTransbase | ChEMBL | DrugBank source literature",
            "acceptance_criteria": "same standardized connectivity and transporter; independent source; agreement or conflict retained",
            "reject_if": "same paper re-exported by two databases counted as independent",
            "model_use": "estimate label reliability and build consensus/conflict priors",
        },
        {
            "priority": "P1",
            "gap_id": "T10_GAP_004",
            "transporter_gene": "ABCB1 | ABCG2",
            "missing_data": "continuous direct endpoints for compounds near ADC payload/release chemical space",
            "why_needed": "binary training compounds are structurally and mechanistically remote from many bRo5 release species",
            "target_sources": "primary ADME papers | regulatory transport studies | patents with methods",
            "acceptance_criteria": "Papp A-B/B-A, ER/net ER, concentration, cell background, inhibitor, recovery or source caveat",
            "reject_if": "whole ADC result substituted for released species; figure-only value without reproducible digitization",
            "model_use": "continuous auxiliary head and ADC-domain transfer",
        },
        {
            "priority": "P1",
            "gap_id": "T10_GAP_005",
            "transporter_gene": "ABCB1 | ABCG2",
            "missing_data": "true negative controls within each assay family",
            "why_needed": "absence of a report is not evidence of nonsubstrate status",
            "target_sources": "bidirectional transport panels and regulatory study tables",
            "acceptance_criteria": "tested negative with valid monolayer/recovery and positive-control performance",
            "reject_if": "no result; low cytotoxicity alone; negative database annotation without assay provenance",
            "model_use": "high-confidence negative class",
        },
    ]


def main() -> None:
    raw = read_csv(RAW_LONG)
    decisions = read_csv(DECISIONS)
    raw_by_id = {row["training_record_id"]: row for row in raw}
    if len(raw) != 1114 or len(decisions) != 1114:
        raise ValueError(f"Unexpected input sizes: raw={len(raw)}, decisions={len(decisions)}")

    audited = audit_rows(decisions, raw_by_id)
    matrix = source_label_matrix(
        [row for row in decisions if row["leakage_model_action"] == "include_leakage_clean_training"]
    )
    truthset = [row for row in audited if row["strict_AB_candidate"] == "yes"]
    auxiliary = [row for row in audited if row["strict_AB_candidate"] != "yes"]
    issues = issue_registry(audited, matrix)
    holdouts = source_holdout_diagnostics(audited)
    tier_cv = tier_b_scaffold_cv(audited)
    queue = gap_queue()

    write_csv(PROCESSED / "transporter_label_provenance_audit_v1.csv", audited)
    write_csv(PROCESSED / "transporter_label_source_matrix_v1.csv", matrix)
    write_csv(PROCESSED / "transporter_truthset_AB_candidate_v1.csv", truthset)
    write_csv(PROCESSED / "transporter_auxiliary_labels_v1.csv", auxiliary)
    write_csv(PROCESSED / "transporter_label_conflicts_v1.csv", issues)
    write_csv(PROCESSED / "source_holdout_diagnostics_v1.csv", holdouts)
    write_csv(PROCESSED / "tier_AB_scaffold_cv_diagnostics_v1.csv", tier_cv)
    write_csv(PROCESSED / "targeted_label_data_gap_queue_v1.csv", queue)

    clean = [row for row in audited if row["leakage_model_action"] == "include_leakage_clean_training"]
    summary = {
        "batch": 10,
        "date": RUN_DATE,
        "purpose": "training_label_provenance_audit_and_source_confounding_repair",
        "input_records": len(audited),
        "leakage_clean_records": len(clean),
        "strict_tier_A_records": sum(str(row["provenance_tier"]).startswith("A_") for row in audited),
        "tier_AB_candidate_records": len(truthset),
        "tier_AB_by_gene_and_label": {
            gene: {
                "negative": sum(row["transporter_gene"] == gene and row["published_binary_label"] == "0" for row in truthset),
                "positive": sum(row["transporter_gene"] == gene and row["published_binary_label"] == "1" for row in truthset),
            }
            for gene in ("ABCB1", "ABCG2")
        },
        "auxiliary_or_excluded_records": len(auxiliary),
        "row_level_primary_assay_provenance_records": sum(row["row_level_primary_assay_id_available"] == "yes" for row in audited),
        "cross_source_consensus_records": sum(int(row["source_record_count"]) > 1 for row in audited),
        "source_label_cramers_v_clean": {
            gene: round(cramers_v([row for row in decisions if row["transporter_gene"] == gene and row["leakage_model_action"] == "include_leakage_clean_training"]), 6)
            for gene in ("ABCB1", "ABCG2")
        },
        "single_class_source_rows": [
            {key: row[key] for key in ("transporter_gene", "source_database", "n_negative", "n_positive")}
            for row in matrix if row["source_label_confounding_flag"] == "yes"
        ],
        "open_audit_issues": len(issues),
        "targeted_gap_queue_rows": len(queue),
        "tier_AB_scaffold_cv": tier_cv,
        "claim_boundary": (
            "No Tier A direct-assay training rows can be reconstructed from Table S2 alone. "
            "Tier AB is a candidate author-curated set, not a final transporter truth set."
        ),
    }
    (PROCESSED / "tenth_batch_label_audit_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
