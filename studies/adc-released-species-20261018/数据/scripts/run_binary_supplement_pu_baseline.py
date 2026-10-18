#!/usr/bin/env python3
"""Build an evidence-aware binary supplement and a PU sensitivity baseline.

This batch is intentionally a *data-method* experiment, not a rescue of the
predictor claim.  ChEMBL/author-curated A/B candidates remain the only rows
used as an ordinary supervised target.  Metrabase/Shaikh-style auxiliary rows
are retained as traceable unlabeled data for positive-unlabeled (PU) bagging;
they are never silently converted into negatives.  ADC direct-transport rows
are aggregated only when the same assay context contains a no-inhibitor and a
selective-inhibitor measurement satisfying the preregistered rule
net-efflux >= 2 and >=50% inhibitor rescue.  Those rows are held out as
external ADC anchors, not used to train the ordinary small-molecule model.

Outputs are CSV/JSON artifacts under data/processed.  No source label is
overwritten and no model_binary_label is written into the existing evidence
tables.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, matthews_corrcoef, roc_auc_score


DATA = Path(__file__).resolve().parents[1]
PROCESSED = DATA / "processed"
RUN_DATE = "2026-07-29"
SEED = 20260729
MORGAN_BITS = 2048
MORGAN_RADIUS = 2
PU_ITERATIONS = 50
FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(
    radius=MORGAN_RADIUS, fpSize=MORGAN_BITS, includeChirality=False
)
RDLogger.DisableLog("rdApp.warning")


def read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(PROCESSED / name, dtype=str, keep_default_na=False)


def write_csv(name: str, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path = PROCESSED / name
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def safe_float(value: object) -> float | None:
    try:
        if value is None or str(value).strip() == "":
            return None
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def fp_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return FP_GENERATOR.GetFingerprint(mol) if mol is not None else None


def fp_matrix(rows: list[dict]) -> tuple[np.ndarray, list[dict]]:
    valid: list[dict] = []
    vectors: list[np.ndarray] = []
    for row in rows:
        fp = fp_from_smiles(row.get("standardized_isomeric_smiles", ""))
        if fp is None:
            continue
        arr = np.zeros((MORGAN_BITS,), dtype=np.float64)
        DataStructs.ConvertToNumpyArray(fp, arr)
        vectors.append(arr)
        valid.append(row)
    return np.vstack(vectors) if vectors else np.empty((0, MORGAN_BITS)), valid


def similarity_support(rows: list[dict]) -> dict:
    """Return empirical self-nearest thresholds for a training subset."""
    fps = [fp_from_smiles(row.get("standardized_isomeric_smiles", "")) for row in rows]
    fps = [fp for fp in fps if fp is not None]
    maxima: list[float] = []
    for idx, fp in enumerate(fps):
        others = fps[:idx] + fps[idx + 1 :]
        maxima.append(max(DataStructs.BulkTanimotoSimilarity(fp, others)) if others else 0.0)
    if not maxima:
        return {"q10": 0.0, "q25": 0.0}
    return {"q10": float(np.quantile(maxima, 0.10)), "q25": float(np.quantile(maxima, 0.25))}


def support_fields(row: dict, training_rows: list[dict], thresholds: dict) -> dict:
    fp = fp_from_smiles(row.get("standardized_isomeric_smiles", ""))
    train_fps = [fp_from_smiles(item.get("standardized_isomeric_smiles", "")) for item in training_rows]
    train_fps = [item for item in train_fps if item is not None]
    maximum = max(DataStructs.BulkTanimotoSimilarity(fp, train_fps)) if fp is not None and train_fps else 0.0
    if maximum < thresholds["q10"]:
        band, abstain = "outside_empirical_candidate_support", "yes_abstain"
    elif maximum < thresholds["q25"]:
        band, abstain = "low_empirical_candidate_support", "yes_caution"
    else:
        band, abstain = "within_empirical_candidate_support", "no_similarity_abstention"
    return {
        "max_training_morgan_tanimoto": round(float(maximum), 6),
        "similarity_support_band": band,
        "similarity_abstain_recommendation": abstain,
    }


def support_fields_with_fps(row: dict, train_fps: list, thresholds: dict) -> dict:
    """Same as support_fields, with precomputed training fingerprints."""
    fp = fp_from_smiles(row.get("standardized_isomeric_smiles", ""))
    maximum = max(DataStructs.BulkTanimotoSimilarity(fp, train_fps)) if fp is not None and train_fps else 0.0
    if maximum < thresholds["q10"]:
        band, abstain = "outside_empirical_candidate_support", "yes_abstain"
    elif maximum < thresholds["q25"]:
        band, abstain = "low_empirical_candidate_support", "yes_caution"
    else:
        band, abstain = "within_empirical_candidate_support", "no_similarity_abstention"
    return {
        "max_training_morgan_tanimoto": round(float(maximum), 6),
        "similarity_support_band": band,
        "similarity_abstain_recommendation": abstain,
    }


def structure_key(row: dict) -> str:
    return row.get("standardized_connectivity_smiles", "") or row.get("standardized_inchikey", "")


def build_interface() -> tuple[list[dict], dict[str, dict], dict[str, dict]]:
    candidates = read_csv("transporter_truthset_AB_candidate_v1.csv")
    auxiliary = read_csv("transporter_auxiliary_labels_v1.csv")
    manual = read_csv("eleventh_batch_manual_verified_binary_labels_v1.csv")
    interface: list[dict] = []
    holdouts: dict[str, dict] = {}
    manual_by_audit = {row["audit_record_id"]: row for row in candidates.to_dict("records")}

    common_fields = [
        "standardized_inchikey",
        "standardized_isomeric_smiles",
        "standardized_connectivity_smiles",
        "murcko_scaffold_smiles",
        "std_molecular_weight",
        "std_xlogp",
        "std_tpsa",
        "std_hbond_donor_count",
        "std_hbond_acceptor_count",
        "std_rotatable_bond_count",
        "std_formal_charge",
        "std_heavy_atom_count",
    ]

    for row in candidates.to_dict("records"):
        if row.get("published_binary_label", "") not in {"0", "1"}:
            continue
        # ChEMBL is the only current source that contains both classes for
        # both genes; other A/B rows are retained but down-weighted.
        source = row.get("source_database", "")
        weight = 1.0 if source == "ChEMBL" else 0.35
        interface.append(
            {
                "interface_id": f"CAND::{row['audit_record_id']}",
                "transporter_gene": row["transporter_gene"],
                "preferred_name": row.get("preferred_names", ""),
                "target_label": row["published_binary_label"],
                "published_label_for_trace": row["published_binary_label"],
                "target_role": "candidate_supervised",
                "evidence_tier": row.get("provenance_tier", ""),
                "source_database": source,
                "source_origin": row.get("source_origin", ""),
                "sample_weight": weight,
                "standardization_key": structure_key(row),
                **{field: row.get(field, "") for field in common_fields},
                "source_record_id": row["audit_record_id"],
                "source_url": row.get("daood_source_url", ""),
                "claim_boundary": "author-curated candidate; not Tier A direct assay truth",
            }
        )

    for row in auxiliary.to_dict("records"):
        if row.get("published_binary_label", "") not in {"0", "1"}:
            continue
        # The published database label is preserved only for provenance.  The
        # PU target is deliberately blank so a missing/weak call is not turned
        # into a negative class.
        interface.append(
            {
                "interface_id": f"AUX::{row['audit_record_id']}",
                "transporter_gene": row["transporter_gene"],
                "preferred_name": row.get("preferred_names", ""),
                "target_label": "",
                "published_label_for_trace": row["published_binary_label"],
                "target_role": "pu_unlabeled",
                "evidence_tier": row.get("provenance_tier", ""),
                "source_database": row.get("source_database", ""),
                "source_origin": row.get("source_origin", ""),
                "sample_weight": 0.0,
                "standardization_key": structure_key(row),
                **{field: row.get(field, "") for field in common_fields},
                "source_record_id": row["audit_record_id"],
                "source_url": row.get("daood_source_url", ""),
                "claim_boundary": "auxiliary/PU only; published label is not a binary target",
            }
        )

    for row in manual.to_dict("records"):
        base = manual_by_audit.get(row["audit_record_id"], {})
        holdout = {
            "interface_id": f"MANUAL::{row['manual_label_id']}",
            "transporter_gene": row["transporter_gene"],
            "preferred_name": row["preferred_name"],
            "target_label": row["binary_label"],
            "published_label_for_trace": row["binary_label"],
            "target_role": "audit_holdout",
            "evidence_tier": row["evidence_grade"],
            "source_database": "primary_manual_review",
            "source_origin": row.get("primary_doi", ""),
            "sample_weight": row.get("weight_recommendation", ""),
            "standardization_key": structure_key(base),
            **{field: base.get(field, "") for field in common_fields},
            "source_record_id": row["manual_label_id"],
            "source_url": row.get("primary_source_url", ""),
            "claim_boundary": row.get("claim_boundary", ""),
        }
        holdouts[row["manual_label_id"]] = holdout
        interface.append(holdout)

    return interface, holdouts, {row["audit_record_id"]: row for row in candidates.to_dict("records")}


def derive_adc_anchors() -> list[dict]:
    direct = read_csv("direct_bidirectional_transport_numeric_v3.csv")
    # Combine the standardized free-payload and released-species registries.
    # v4 is the raw release registry and does not carry standardized fields;
    # it is intentionally not used as a modeling structure source here.
    registry_rows: list[dict] = []
    registry_rows.extend(read_csv("compound_registry_standardized_v1.csv").to_dict("records"))
    registry_rows.extend(read_csv("adc_release_species_registry_standardized_v1.csv").to_dict("records"))
    reg_by_name = {row.get("reference_name", "").strip().lower(): row for row in registry_rows}
    groups: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in direct.to_dict("records"):
        gene = row.get("transporter_gene", "")
        if gene not in {"ABCB1", "ABCG2"}:
            continue
        name = row.get("entity_name", "").strip()
        if not name or row.get("entity_type") not in {"free_payload", "released_species"}:
            continue
        group_id = row.get("comparability_group_id", "") or row.get("source_id", "")
        groups[(gene, name, group_id)].append(row)

    anchors: list[dict] = []
    for (gene, name, group_id), rows in sorted(groups.items()):
        no_inhib = [r for r in rows if not r.get("inhibitor", "").strip() or r.get("inhibitor", "").strip().lower() in {"none", "no inhibitor"}]
        with_inhib = [r for r in rows if r.get("inhibitor", "").strip() and r.get("inhibitor", "").strip().lower() not in {"none", "no inhibitor"}]

        def er(row: dict) -> float | None:
            return safe_float(row.get("background_corrected_efflux_ratio")) or safe_float(row.get("efflux_ratio_mean"))

        er0_values = [er(r) for r in no_inhib if er(r) is not None]
        eri_values = [er(r) for r in with_inhib if er(r) is not None]
        er0 = max(er0_values) if er0_values else None
        eri = min(eri_values) if eri_values else None
        rescue = 1.0 - eri / er0 if er0 is not None and eri is not None and er0 > 0 else None
        strict_label = ""
        derivation = "no_same_context_inhibitor_pair"
        if er0 is not None and eri is not None:
            derivation = "net_efflux_ge_2_and_inhibitor_rescue_ge_50pct"
            if er0 >= 2.0 and rescue is not None and rescue >= 0.50:
                strict_label = "1"
            elif er0 < 2.0:
                # A low ER can support a negative only in a gene-isolated,
                # QC-passing assay.  Failure of a selective inhibitor to rescue
                # a high ER is not a nonsubstrate call; another transporter or
                # assay artifact may explain it, so retain the label blank.
                strict_label = "0" if gene_specific else ""
            else:
                strict_label = ""

        registry_row = reg_by_name.get(name.lower(), {})
        # Names such as “DXd / MAAA-1181a” are intentionally not fuzzy-mapped.
        if not registry_row:
            continue
        source_ids = sorted({r.get("source_id", "") for r in rows if r.get("source_id", "")})
        source_urls = sorted({r.get("source_url", "") for r in rows if r.get("source_url", "")})
        gene_specific = all("Caco-2" not in r.get("transporter_system", "") and "multi-transporter" not in r.get("transporter_system", "") for r in rows)
        anchors.append(
            {
                "anchor_id": f"ADC_ANCHOR_{len(anchors)+1:03d}",
                "entity_name": name,
                "entity_id": rows[0].get("entity_id", ""),
                "entity_type": rows[0].get("entity_type", ""),
                "transporter_gene": gene,
                "comparability_group_id": group_id,
                "no_inhibitor_efflux_ratio": "" if er0 is None else round(er0, 6),
                "inhibitor_efflux_ratio": "" if eri is None else round(eri, 6),
                "inhibitor_rescue_fraction": "" if rescue is None else round(rescue, 6),
                "derived_binary_label": strict_label,
                "derivation_rule": derivation,
                "gene_specific_assay_context": "yes" if gene_specific else "no",
                "external_anchor_role": "external_ADC_anchor_only",
                "truth_set_eligible": "yes_same_context_numeric_pair" if strict_label and gene_specific else "no_context_limited_or_unpaired",
                "standardized_inchikey": registry_row.get("standardized_inchikey", ""),
                "standardized_isomeric_smiles": registry_row.get("standardized_isomeric_smiles", ""),
                "standardized_connectivity_smiles": registry_row.get("standardized_connectivity_smiles", ""),
                "std_molecular_weight": registry_row.get("std_molecular_weight", ""),
                "std_xlogp": registry_row.get("std_xlogp", ""),
                "std_tpsa": registry_row.get("std_tpsa", ""),
                "source_ids": " | ".join(source_ids),
                "source_urls": " | ".join(source_urls),
                "source_record_count": len(rows),
                "claim_boundary": "binary companion is an assay-context anchor; not a universal ADC truth label",
            }
        )
    return anchors


def build_public_external_holdout() -> list[dict]:
    """Build a structure-disjoint, source-held-out sensitivity set.

    The two public collections are B2 curated literature datasets, not Tier-A
    direct assays.  We therefore use only records with no Daood overlap and
    keep them as an external source diagnostic.  Duplicate connectivity keys
    inside a source are collapsed to the first traceable row to avoid counting
    repeated literature entries as independent molecules.
    """
    # The strict candidate view explicitly excludes Lavecchia inhibitors and
    # non-active-only broad negatives, and is already restricted to parseable
    # structures with no Daood exact/connectivity overlap.
    labels = read_csv("public_independent_strict_binary_candidate_view_v1.csv")
    merged = labels.copy()
    merged["binary_label_numeric"] = pd.to_numeric(merged["binary_label"], errors="coerce")
    merged = merged[
        merged["binary_label_numeric"].isin([0, 1])
        & (merged["overlap_type"] == "no_daood_overlap")
        & (merged["source_structure_status"] != "failed")
    ].copy()
    rows: list[dict] = []
    seen: set[tuple[str, str]] = set()
    for row in merged.to_dict("records"):
        structure = row.get("standardized_smiles", "") or row.get("source_smiles", "")
        key = (row.get("transporter_gene", ""), row.get("standardized_inchikey", "") or structure)
        if not structure or key in seen:
            continue
        seen.add(key)
        rows.append(
            {
                "interface_id": f"EXT::{row['independent_record_id']}",
                "transporter_gene": row["transporter_gene"],
                "preferred_name": row.get("compound_name", ""),
                "target_label": str(int(row["binary_label_numeric"])),
                "published_label_for_trace": str(int(row["binary_label_numeric"])),
                "target_role": "public_source_holdout",
                "evidence_tier": row.get("evidence_grade", "B2"),
                "source_database": row.get("source_doi", ""),
                "source_origin": row.get("source_class", ""),
                "sample_weight": 0.0,
                "standardization_key": key[1],
                "standardized_inchikey": row.get("standardized_inchikey", ""),
                "standardized_isomeric_smiles": structure,
                "standardized_connectivity_smiles": "",
                "murcko_scaffold_smiles": "",
                "std_molecular_weight": "",
                "std_xlogp": "",
                "std_tpsa": "",
                "std_hbond_donor_count": "",
                "std_hbond_acceptor_count": "",
                "std_rotatable_bond_count": "",
                "std_formal_charge": "",
                "std_heavy_atom_count": "",
                "source_record_id": row["independent_record_id"],
                "source_url": row.get("source_url", ""),
                "claim_boundary": "B2 independent curated source; source-heldout sensitivity only",
            }
        )
    return rows


def metrics(y_true: np.ndarray, prob: np.ndarray) -> dict:
    pred = (prob >= 0.5).astype(int)
    both = len(set(y_true.tolist())) == 2
    out = {
        "n": len(y_true),
        "n_positive": int(y_true.sum()),
        "mean_probability": float(np.mean(prob)) if len(prob) else math.nan,
        "balanced_accuracy": float(balanced_accuracy_score(y_true, pred)) if both else math.nan,
        "mcc": float(matthews_corrcoef(y_true, pred)) if both else math.nan,
        "roc_auc": float(roc_auc_score(y_true, prob)) if both else math.nan,
        "average_precision": float(average_precision_score(y_true, prob)) if both else math.nan,
        "metric_status": "binary_holdout_metrics" if both else "single_class_holdout_descriptive_only",
    }
    if not both and len(y_true):
        out["single_class_correct_fraction"] = float(np.mean(pred == y_true))
    return out


def train_models(interface: list[dict], holdouts: dict[str, dict], anchors: list[dict], external_holdout: list[dict]) -> tuple[list[dict], list[dict]]:
    rng = np.random.default_rng(SEED)
    diagnostics: list[dict] = []
    prediction_rows: list[dict] = []
    manual_rows = list(holdouts.values())
    holdout_keys = {row.get("standardization_key", "") for row in manual_rows if row.get("standardization_key", "")}

    for gene in ["ABCB1", "ABCG2"]:
        candidates = [
            row for row in interface
            if row["target_role"] == "candidate_supervised"
            and row["transporter_gene"] == gene
            and row["target_label"] in {"0", "1"}
            and row.get("standardization_key", "") not in holdout_keys
        ]
        aux = [row for row in interface if row["target_role"] == "pu_unlabeled" and row["transporter_gene"] == gene and row.get("standardization_key", "") not in holdout_keys]
        holdout_gene = [row for row in manual_rows if row["transporter_gene"] == gene and row.get("standardization_key", "")]
        anchor_gene = [row for row in anchors if row["transporter_gene"] == gene and row.get("standardized_isomeric_smiles", "")]
        external_gene = [row for row in external_holdout if row["transporter_gene"] == gene and row.get("standardized_isomeric_smiles", "")]

        x_cand, cand_valid = fp_matrix(candidates)
        y_cand = np.array([int(row["target_label"]) for row in cand_valid], dtype=int)
        if len(cand_valid) < 10 or len(set(y_cand.tolist())) < 2:
            continue
        support_thresholds = similarity_support(cand_valid)
        train_support_fps = [fp_from_smiles(row.get("standardized_isomeric_smiles", "")) for row in cand_valid]
        train_support_fps = [fp for fp in train_support_fps if fp is not None]
        estimators: dict[str, tuple[object, list[dict], np.ndarray]] = {}
        est = LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000, solver="liblinear", random_state=SEED)
        est.fit(x_cand, y_cand)
        estimators["candidate_supervised"] = (est, cand_valid, y_cand)

        # Weak sensitivity model: auxiliary labels are admitted only with
        # explicit source weights, and this result can never replace PU output.
        weighted_rows = candidates + aux
        x_weak, weak_valid = fp_matrix(weighted_rows)
        y_weak = np.array([int(row["published_label_for_trace"]) for row in weak_valid], dtype=int)
        w_weak = np.array([float(row.get("sample_weight") or 0.25) if row["target_role"] == "candidate_supervised" else 0.25 for row in weak_valid])
        if len(weak_valid) and len(set(y_weak.tolist())) == 2:
            weak = LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000, solver="liblinear", random_state=SEED)
            weak.fit(x_weak, y_weak, sample_weight=w_weak)
            estimators["weighted_auxiliary_sensitivity"] = (weak, weak_valid, y_weak)

        # Positive-unlabeled bagging: sample temporary negatives from the
        # auxiliary pool, repeat, and report mean + SD instead of a hard label.
        pu_models = []
        x_aux, aux_valid = fp_matrix(aux)
        if len(aux_valid) and len(set(y_cand.tolist())) == 2:
            n_temp_neg = min(len(aux_valid), max(10, int(y_cand.sum())))
            for _ in range(PU_ITERATIONS):
                neg_idx = rng.choice(len(aux_valid), size=n_temp_neg, replace=len(aux_valid) < n_temp_neg)
                train_rows = cand_valid + [aux_valid[idx] for idx in neg_idx]
                x_train, train_valid = fp_matrix(train_rows)
                y_train = np.array([1] * len(cand_valid) + [0] * (len(train_valid) - len(cand_valid)), dtype=int)
                if len(set(y_train.tolist())) < 2:
                    continue
                model = LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000, solver="liblinear", random_state=SEED)
                model.fit(x_train, y_train)
                pu_models.append(model)

        # Holdout diagnostics are audit checks only; the seven labels overlap
        # the historical candidate set and are not an independent validation.
        for model_name, (model, _, _) in estimators.items():
            x_hold, hold_valid = fp_matrix(holdout_gene)
            if len(hold_valid):
                prob = model.predict_proba(x_hold)[:, 1]
                y = np.array([int(row["target_label"]) for row in hold_valid], dtype=int)
                diagnostics.append({"run_date": RUN_DATE, "transporter_gene": gene, "model": model_name, "holdout_type": "manual_primary_audit_holdout", **metrics(y, prob), "claim_boundary": "audit diagnostic only; no independent validation"})
                for row, p in zip(hold_valid, prob):
                    prediction_rows.append({"prediction_set": "manual_primary_audit_holdout", "model": model_name, "transporter_gene": gene, "entity_name": row["preferred_name"], "source_record_id": row["source_record_id"], "target_label": row["target_label"], "predicted_probability": round(float(p), 8), "support_status": "audit_holdout_overlap_excluded_from_training", **support_fields_with_fps(row, train_support_fps, support_thresholds), "claim_boundary": "audit diagnostic only"})

            x_anchor, anchor_valid = fp_matrix(anchor_gene)
            if len(anchor_valid):
                prob = model.predict_proba(x_anchor)[:, 1]
                for row, p in zip(anchor_valid, prob):
                    prediction_rows.append({"prediction_set": "adc_external_anchor", "model": model_name, "transporter_gene": gene, "entity_name": row["entity_name"], "source_record_id": row["anchor_id"], "target_label": row.get("derived_binary_label", ""), "predicted_probability": round(float(p), 8), "support_status": row.get("truth_set_eligible", ""), **support_fields_with_fps(row, train_support_fps, support_thresholds), "claim_boundary": "external ADC anchor; not training data"})

            x_ext, ext_valid = fp_matrix(external_gene)
            if len(ext_valid):
                ext_prob = model.predict_proba(x_ext)[:, 1]
                ext_y = np.array([int(row["target_label"]) for row in ext_valid], dtype=int)
                diagnostics.append({"run_date": RUN_DATE, "transporter_gene": gene, "model": model_name, "holdout_type": "public_independent_source_no_overlap", **metrics(ext_y, ext_prob), "claim_boundary": "B2 curated literature source; source-heldout sensitivity only"})
                for row, p in zip(ext_valid, ext_prob):
                    prediction_rows.append({"prediction_set": "public_independent_source_holdout", "model": model_name, "transporter_gene": gene, "entity_name": row["preferred_name"], "source_record_id": row["source_record_id"], "target_label": row["target_label"], "predicted_probability": round(float(p), 8), "support_status": "no_daood_overlap", **support_fields_with_fps(row, train_support_fps, support_thresholds), "claim_boundary": "B2 curated literature source; sensitivity only"})

        if pu_models:
            # PU model is evaluated as a ranking/probability distribution only.
            def pu_predict(rows: list[dict]) -> tuple[np.ndarray, np.ndarray]:
                x, valid = fp_matrix(rows)
                probs = np.vstack([m.predict_proba(x)[:, 1] for m in pu_models])
                return probs.mean(axis=0), probs.std(axis=0)

            x_hold, hold_valid = fp_matrix(holdout_gene)
            if len(hold_valid):
                mean_p, sd_p = pu_predict(holdout_gene)
                y = np.array([int(row["target_label"]) for row in hold_valid], dtype=int)
                diagnostics.append({"run_date": RUN_DATE, "transporter_gene": gene, "model": "pu_bagging_rank", "holdout_type": "manual_primary_audit_holdout", **metrics(y, mean_p), "mean_probability_sd": float(np.mean(sd_p)), "claim_boundary": "PU ranking diagnostic only; not a calibrated binary truth"})
                for row, p, s in zip(hold_valid, mean_p, sd_p):
                    prediction_rows.append({"prediction_set": "manual_primary_audit_holdout", "model": "pu_bagging_rank", "transporter_gene": gene, "entity_name": row["preferred_name"], "source_record_id": row["source_record_id"], "target_label": row["target_label"], "predicted_probability": round(float(p), 8), "pu_probability_sd": round(float(s), 8), "support_status": "audit_holdout_overlap_excluded_from_training", **support_fields_with_fps(row, train_support_fps, support_thresholds), "claim_boundary": "PU ranking diagnostic only"})
            if len(aux_valid):
                mean_p, sd_p = pu_predict(aux_valid)
                for row, p, s in zip(aux_valid, mean_p, sd_p):
                    prediction_rows.append({"prediction_set": "auxiliary_pu_unlabeled", "model": "pu_bagging_rank", "transporter_gene": gene, "entity_name": row["preferred_name"], "source_record_id": row["source_record_id"], "target_label": "", "published_label_for_trace": row.get("published_label_for_trace", ""), "predicted_probability": round(float(p), 8), "pu_probability_sd": round(float(s), 8), "support_status": "unlabeled_auxiliary", **support_fields_with_fps(row, train_support_fps, support_thresholds), "claim_boundary": "ranking only; published auxiliary label not used as target"})
            if len(anchor_gene):
                mean_p, sd_p = pu_predict(anchor_gene)
                for row, p, s in zip(anchor_gene, mean_p, sd_p):
                    prediction_rows.append({"prediction_set": "adc_external_anchor", "model": "pu_bagging_rank", "transporter_gene": gene, "entity_name": row["entity_name"], "source_record_id": row["anchor_id"], "target_label": row.get("derived_binary_label", ""), "predicted_probability": round(float(p), 8), "pu_probability_sd": round(float(s), 8), "support_status": row.get("truth_set_eligible", ""), **support_fields_with_fps(row, train_support_fps, support_thresholds), "claim_boundary": "PU ranking only; external ADC anchor not training data"})
            if len(external_gene):
                mean_p, sd_p = pu_predict(external_gene)
                ext_y = np.array([int(row["target_label"]) for row in external_gene], dtype=int)
                diagnostics.append({"run_date": RUN_DATE, "transporter_gene": gene, "model": "pu_bagging_rank", "holdout_type": "public_independent_source_no_overlap", **metrics(ext_y, mean_p), "mean_probability_sd": float(np.mean(sd_p)), "claim_boundary": "PU ranking on B2 source-heldout set; not calibrated"})
                for row, p, s in zip(external_gene, mean_p, sd_p):
                    prediction_rows.append({"prediction_set": "public_independent_source_holdout", "model": "pu_bagging_rank", "transporter_gene": gene, "entity_name": row["preferred_name"], "source_record_id": row["source_record_id"], "target_label": row["target_label"], "predicted_probability": round(float(p), 8), "pu_probability_sd": round(float(s), 8), "support_status": "no_daood_overlap", **support_fields_with_fps(row, train_support_fps, support_thresholds), "claim_boundary": "PU ranking on B2 source-heldout set; sensitivity only"})

    return diagnostics, prediction_rows


def main() -> None:
    interface, holdouts, _ = build_interface()
    anchors = derive_adc_anchors()
    external_holdout = build_public_external_holdout()
    # Add the external-anchor rows to a separate interface table so the role
    # and non-training boundary are machine-readable.
    for row in anchors:
        interface.append({
            "interface_id": f"ADC::{row['anchor_id']}",
            "transporter_gene": row["transporter_gene"],
            "preferred_name": row["entity_name"],
            "target_label": row.get("derived_binary_label", ""),
            "published_label_for_trace": "",
            "target_role": "external_adc_anchor",
            "evidence_tier": "A_same_context_numeric_pair" if row.get("derived_binary_label") else "A_numeric_unpaired_or_context_limited",
            "source_database": "project_direct_transport",
            "source_origin": row.get("source_ids", ""),
            "sample_weight": 0.0,
            "standardization_key": row.get("standardized_connectivity_smiles", ""),
            "standardized_inchikey": row.get("standardized_inchikey", ""),
            "standardized_isomeric_smiles": row.get("standardized_isomeric_smiles", ""),
            "standardized_connectivity_smiles": row.get("standardized_connectivity_smiles", ""),
            "murcko_scaffold_smiles": "",
            "std_molecular_weight": row.get("std_molecular_weight", ""),
            "std_xlogp": row.get("std_xlogp", ""),
            "std_tpsa": row.get("std_tpsa", ""),
            "std_hbond_donor_count": "",
            "std_hbond_acceptor_count": "",
            "std_rotatable_bond_count": "",
            "std_formal_charge": "",
            "std_heavy_atom_count": "",
            "source_record_id": row["anchor_id"],
            "source_url": row.get("source_urls", ""),
            "claim_boundary": row["claim_boundary"],
        })

    interface_fields = list(interface[0]) if interface else []
    write_csv("binary_supplement_training_interface_v1.csv", interface, interface_fields)
    write_csv("adc_direct_numeric_binary_anchors_v1.csv", anchors)
    diagnostics, predictions = train_models(interface, holdouts, anchors, external_holdout)
    write_csv("binary_supplement_model_diagnostics_v1.csv", diagnostics)
    write_csv("binary_supplement_predictions_v1.csv", predictions)

    summary = {
        "run_date": RUN_DATE,
        "seed": SEED,
        "interface_counts": {
            f"{role}|{gene}": int(count)
            for (role, gene), count in pd.DataFrame(interface).groupby(["target_role", "transporter_gene"]).size().items()
        } if interface else {},
        "adc_anchor_counts": {
            f"{gene}|{label}": int(count)
            for (gene, label), count in pd.DataFrame(anchors).groupby(["transporter_gene", "derived_binary_label"]).size().items()
        } if anchors else {},
        "manual_holdout_count": len(holdouts),
        "public_independent_source_holdout_count": len(external_holdout),
        "pu_iterations": PU_ITERATIONS,
        "pu_rule": "labeled positives from A/B candidates; auxiliary rows sampled as temporary negatives across repeated fits; output is ranking/probability distribution only",
        "derived_binary_rule": "same comparability group, no-inhibitor net efflux ratio >= 2 and selective-inhibitor rescue >= 50%; gene-specific context required for truth_set_eligible",
        "model_claim": "development sensitivity baseline only; no universal ordinary-substrate or ADC predictor claim",
        "artifacts": [
            "binary_supplement_training_interface_v1.csv",
            "adc_direct_numeric_binary_anchors_v1.csv",
            "binary_supplement_model_diagnostics_v1.csv",
            "binary_supplement_predictions_v1.csv",
            "public_independent_transporter_labels_long_v1.csv",
            "public_independent_cross_source_overlap_audit_v1.csv",
        ],
    }
    (PROCESSED / "binary_supplement_summary_v1.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
