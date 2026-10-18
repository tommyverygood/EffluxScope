#!/usr/bin/env python3
"""Run the frozen post-hoc public-predictor rejection benchmark."""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef,
    roc_auc_score,
)
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "数据" / "processed"
VENDOR = ROOT / "数据" / "vendor" / "public_predictor_benchmark_v1"
INPUT = PROCESSED / "abcb1_public_predictor_benchmark_input_v1.csv"
PROTOCOL = PROCESSED / "abcb1_public_predictor_benchmark_protocol_v1.json"
AZ_RAW = PROCESSED / "abcb1_public_predictor_astrazeneca_gnn_mtl_raw_v1.csv"
ADMET_AI_RAW = PROCESSED / "abcb1_public_predictor_admet_ai_v2_raw_v1.csv"
DEEP_PK_RESPONSES = [
    PROCESSED / "abcb1_public_predictor_deeppk_response_batch1_v1.json",
    PROCESSED / "abcb1_public_predictor_deeppk_response_batch2_v1.json",
]
DEEP_PK_DATASETS = {
    "train": VENDOR / "deeppk" / "pgp_substrate_train.csv",
    "validation": VENDOR / "deeppk" / "pgp_substrate_val.csv",
    "test": VENDOR / "deeppk" / "pgp_substrate_test.csv",
}

PREDICTIONS = PROCESSED / "abcb1_public_predictor_predictions_long_v1.csv"
DEEP_PK_PARSED = PROCESSED / "abcb1_public_predictor_deeppk_parsed_v1.csv"
OVERLAP = PROCESSED / "abcb1_public_predictor_deeppk_overlap_audit_v1.csv"
METRICS = PROCESSED / "abcb1_public_predictor_metrics_v1.csv"
COVERAGE = PROCESSED / "abcb1_public_predictor_coverage_curve_v1.csv"
BOOTSTRAP = PROCESSED / "abcb1_public_predictor_cluster_bootstrap_v1.csv"
SUMMARY = PROCESSED / "abcb1_public_predictor_benchmark_summary_v1.json"

SEED = 20260805
N_BOOT = 2000
MATCHED_N = 453
MATCHED_COVERAGE = MATCHED_N / 2624
COVERAGE_GRID = [0.05, 0.10, MATCHED_COVERAGE, 0.25, 0.50, 0.75, 1.00]
DEEP_PK_PROBABILITY_KEY = "[Absorption/P-Glycoprotein Substrate] Probability"
UNCHARGER = rdMolStandardize.Uncharger(canonicalOrder=True)
TAUTOMER_ENUMERATOR = rdMolStandardize.TautomerEnumerator()


COMPARATORS = [
    {
        "comparator_id": "internal_v0_4_primary_score",
        "display_name": "Internal weak-label score",
        "role": "internal_reference",
        "endpoint": "heterogeneous weak ABCB1 substrate risk",
        "threshold": 0.5,
        "source": "input",
        "column": "v0_4_primary_score",
        "eligible_primary": True,
        "training_overlap": "known reference identities already removed",
    },
    {
        "comparator_id": "astrazeneca_gnn_mtl_mdck_er_2025",
        "display_name": "AZ GNN-MTL MDCK ER",
        "role": "same_assay_family_continuous",
        "endpoint": "MDR1-MDCK LogER",
        "threshold": math.log10(2),
        "source": "az",
        "column": "pred_2",
        "eligible_primary": True,
        "training_overlap": "not auditable: proprietary AZ structures",
    },
    {
        "comparator_id": "astrazeneca_gnn_mtl_nih_mdck_er_2025",
        "display_name": "AZ GNN-MTL NIH-MDCK ER",
        "role": "same_transporter_assay_continuous",
        "endpoint": "NIH MDR1-MDCK LogER",
        "threshold": math.log10(2),
        "source": "az",
        "column": "pred_3",
        "eligible_primary": True,
        "training_overlap": "not auditable: proprietary AZ structures",
    },
    {
        "comparator_id": "deeppk_pgp_substrate_2026_08_05",
        "display_name": "Deep-PK P-gp substrate",
        "role": "public_weak_substrate_probability",
        "endpoint": "P-gp substrate probability",
        "threshold": 0.5,
        "source": "deeppk",
        "column": "deeppk_pgp_substrate_probability",
        "eligible_primary": True,
        "training_overlap": "audited from downloaded train/validation/test structures",
    },
    {
        "comparator_id": "admet_ai_v2_pgp_broccatelli_inhibition",
        "display_name": "ADMET-AI P-gp inhibition",
        "role": "endpoint_mismatch_negative_control",
        "endpoint": "P-gp inhibition probability, not substrate",
        "threshold": 0.5,
        "source": "admet_ai",
        "column": "Pgp_Broccatelli",
        "eligible_primary": False,
        "training_overlap": "not evaluated because endpoint is ineligible",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def connectivity_block(smiles: str) -> str:
    try:
        mol = Chem.MolFromSmiles(str(smiles), sanitize=True)
        if mol is None:
            return ""
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol, skipStandardize=True)
        mol = UNCHARGER.uncharge(mol)
        result = TAUTOMER_ENUMERATOR.Enumerate(mol)
        mol = TAUTOMER_ENUMERATOR.PickCanonical(result)
        return Chem.MolToInchiKey(mol).split("-")[0]
    except Exception:
        return ""


def parse_deeppk() -> pd.DataFrame:
    rows: list[dict] = []
    for batch, path in enumerate(DEEP_PK_RESPONSES, start=1):
        outer = json.loads(path.read_text())
        payload = json.loads(outer) if isinstance(outer, str) else outer
        ordered = sorted(payload.items(), key=lambda item: int(item[0]))
        for within_batch_index, (_, record) in enumerate(ordered):
            rows.append(
                {
                    "batch": batch,
                    "within_batch_index": within_batch_index,
                    "smiles_returned": record["SMILES"],
                    "deeppk_pgp_substrate_probability": float(record[DEEP_PK_PROBABILITY_KEY]),
                    "deeppk_pgp_substrate_call": record[
                        "[Absorption/P-Glycoprotein Substrate] Predictions"
                    ],
                    "deeppk_pgp_substrate_interpretation": record[
                        "[Absorption/P-Glycoprotein Substrate] Interpretation"
                    ],
                }
            )
    parsed = pd.DataFrame(rows)
    parsed.to_csv(DEEP_PK_PARSED, index=False)
    return parsed


def audit_deeppk_overlap(benchmark: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    mapping: dict[str, set[str]] = defaultdict(set)
    vendor_rows: list[dict] = []
    for split, path in DEEP_PK_DATASETS.items():
        data = pd.read_csv(path)
        for index, row in data.iterrows():
            block = connectivity_block(row["smiles_standarized"])
            if block:
                mapping[block].add(split)
            vendor_rows.append(
                {
                    "deeppk_split": split,
                    "deeppk_row_index": int(index),
                    "deeppk_smiles": row["smiles_standarized"],
                    "deeppk_label": row["label"],
                    "connectivity_block": block,
                    "standardization_status": "ok" if block else "failed",
                }
            )

    audit = benchmark[
        ["external_row_id", "standardized_inchikey", "smiles", "observed_er2_label"]
    ].copy()
    audit["connectivity_block"] = audit["standardized_inchikey"].str.split("-").str[0]
    audit["deeppk_overlap_splits"] = audit["connectivity_block"].map(
        lambda key: "|".join(sorted(mapping.get(key, set())))
    )
    for split in DEEP_PK_DATASETS:
        audit[f"deeppk_exact_{split}_overlap"] = audit["deeppk_overlap_splits"].str.split("|").map(
            lambda values: split in values
        )
    audit["deeppk_development_overlap"] = (
        audit["deeppk_exact_train_overlap"] | audit["deeppk_exact_validation_overlap"]
    )
    audit["deeppk_development_identity_clean"] = ~audit["deeppk_development_overlap"]
    audit.to_csv(OVERLAP, index=False)

    vendor = pd.DataFrame(vendor_rows)
    overlap_summary = {
        "vendor_rows": int(len(vendor)),
        "vendor_standardization_failures": int(vendor["connectivity_block"].eq("").sum()),
        "benchmark_exact_train_overlap": int(audit["deeppk_exact_train_overlap"].sum()),
        "benchmark_exact_validation_overlap": int(audit["deeppk_exact_validation_overlap"].sum()),
        "benchmark_exact_test_overlap": int(audit["deeppk_exact_test_overlap"].sum()),
        "benchmark_development_overlap": int(audit["deeppk_development_overlap"].sum()),
        "benchmark_development_identity_clean": int(audit["deeppk_development_identity_clean"].sum()),
    }
    return audit, overlap_summary


def safe_div(num: float, den: float) -> float:
    return float(num / den) if den else math.nan


def calculate_metrics(frame: pd.DataFrame) -> dict:
    y = frame["observed_er2_label"].astype(int).to_numpy()
    score = frame["score"].astype(float).to_numpy()
    pred = frame["predicted_label"].astype(int).to_numpy()
    tp = int(((y == 1) & (pred == 1)).sum())
    tn = int(((y == 0) & (pred == 0)).sum())
    fp = int(((y == 0) & (pred == 1)).sum())
    fn = int(((y == 1) & (pred == 0)).sum())
    both = len(np.unique(y)) == 2
    rho = spearmanr(score, frame["observed_log10_mdr1_mdck_er"].to_numpy()).statistic
    return {
        "n": int(len(frame)),
        "positive_prevalence": float(y.mean()) if len(y) else math.nan,
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)) if both else math.nan,
        "mcc": float(matthews_corrcoef(y, pred)) if both else math.nan,
        "sensitivity": safe_div(tp, tp + fn),
        "specificity": safe_div(tn, tn + fp),
        "ppv": safe_div(tp, tp + fp),
        "npv": safe_div(tn, tn + fn),
        "roc_auc": float(roc_auc_score(y, score)) if both else math.nan,
        "average_precision": float(average_precision_score(y, score)) if both else math.nan,
        "spearman_rho": float(rho),
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
    }


def build_long_predictions() -> tuple[pd.DataFrame, dict]:
    benchmark = pd.read_csv(INPUT)
    az = pd.read_csv(AZ_RAW)
    admet_ai = pd.read_csv(ADMET_AI_RAW)
    deeppk = parse_deeppk()

    if len(benchmark) != 2624 or len(az) != 2624 or len(admet_ai) != 2624 or len(deeppk) != 2624:
        raise ValueError("Every available comparator must return exactly 2,624 rows")
    if not benchmark["smiles"].equals(az["smiles"]):
        raise ValueError("AstraZeneca output order differs from the frozen input")
    if not benchmark["smiles"].equals(admet_ai["smiles"]):
        raise ValueError("ADMET-AI output order differs from the frozen input")
    if not benchmark["smiles"].equals(deeppk["smiles_returned"]):
        mismatch = int((benchmark["smiles"] != deeppk["smiles_returned"]).sum())
        raise ValueError(f"Deep-PK returned {mismatch} SMILES in a different order")

    overlap, overlap_summary = audit_deeppk_overlap(benchmark)
    score_sources = {
        "input": benchmark,
        "az": az,
        "deeppk": deeppk,
        "admet_ai": admet_ai,
    }
    base_columns = [
        "external_row_id",
        "smiles",
        "standardized_inchikey",
        "observed_log10_mdr1_mdck_er",
        "observed_mdr1_mdck_er",
        "observed_er2_label",
        "bootstrap_cluster",
        "v0_4_accepted",
        "v0_4_risk_call",
    ]
    frames: list[pd.DataFrame] = []
    for spec in COMPARATORS:
        frame = benchmark[base_columns].copy()
        frame["comparator_id"] = spec["comparator_id"]
        frame["display_name"] = spec["display_name"]
        frame["role"] = spec["role"]
        frame["endpoint"] = spec["endpoint"]
        frame["eligible_primary_comparison"] = spec["eligible_primary"]
        frame["training_overlap_status"] = spec["training_overlap"]
        frame["score"] = score_sources[spec["source"]][spec["column"]].astype(float).to_numpy()
        frame["binary_threshold"] = spec["threshold"]
        frame["predicted_label"] = (frame["score"] >= spec["threshold"]).astype(int)
        frame["native_margin"] = (frame["score"] - spec["threshold"]).abs()
        frame["native_margin_rank"] = frame["native_margin"].rank(
            method="first", ascending=False
        ).astype(int)
        frame["native_margin_matched_453"] = frame["native_margin_rank"] <= MATCHED_N
        if spec["source"] == "deeppk":
            frame = frame.merge(
                overlap[
                    [
                        "external_row_id",
                        "deeppk_overlap_splits",
                        "deeppk_development_overlap",
                        "deeppk_development_identity_clean",
                    ]
                ],
                on="external_row_id",
                how="left",
                validate="one_to_one",
            )
        else:
            frame["deeppk_overlap_splits"] = "not_applicable"
            frame["deeppk_development_overlap"] = False
            frame["deeppk_development_identity_clean"] = True
        frames.append(frame)

    long = pd.concat(frames, ignore_index=True)
    long.to_csv(PREDICTIONS, index=False)
    return long, overlap_summary


def analyze_metrics(long: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for comparator_id, group in long.groupby("comparator_id", sort=False):
        populations = {
            "all_identity_clean": group,
            "v0_4_accepted_mask": group.loc[group["v0_4_accepted"]],
            "v0_4_rejected_forced": group.loc[~group["v0_4_accepted"]],
            "native_margin_matched_453": group.loc[group["native_margin_matched_453"]],
        }
        if comparator_id == "deeppk_pgp_substrate_2026_08_05":
            populations["deeppk_development_identity_clean"] = group.loc[
                group["deeppk_development_identity_clean"]
            ]
            populations["deeppk_development_clean_v0_4_accepted"] = group.loc[
                group["deeppk_development_identity_clean"] & group["v0_4_accepted"]
            ]
        for population, frame in populations.items():
            metrics = calculate_metrics(frame)
            rows.append(
                {
                    "comparator_id": comparator_id,
                    "display_name": group["display_name"].iloc[0],
                    "role": group["role"].iloc[0],
                    "eligible_primary_comparison": bool(group["eligible_primary_comparison"].iloc[0]),
                    "population": population,
                    "coverage_of_2624": len(frame) / 2624,
                    **metrics,
                }
            )
    metrics = pd.DataFrame(rows)
    metrics.to_csv(METRICS, index=False)
    return metrics


def analyze_coverage(long: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for comparator_id, group in long.groupby("comparator_id", sort=False):
        ordered = group.sort_values(["native_margin", "external_row_id"], ascending=[False, True])
        for target in COVERAGE_GRID:
            n_keep = len(group) if target == 1.0 else max(1, int(round(target * len(group))))
            frame = ordered.iloc[:n_keep]
            rows.append(
                {
                    "comparator_id": comparator_id,
                    "display_name": group["display_name"].iloc[0],
                    "role": group["role"].iloc[0],
                    "target_coverage": target,
                    "n": n_keep,
                    "coverage": n_keep / len(group),
                    **calculate_metrics(frame),
                }
            )
    coverage = pd.DataFrame(rows)
    coverage.to_csv(COVERAGE, index=False)
    return coverage


def analyze_bootstrap(long: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows: list[dict] = []
    for comparator_id, group in long.groupby("comparator_id", sort=False):
        group = group.reset_index(drop=True)
        cluster_to_indices = {
            cluster: indices.to_numpy()
            for cluster, indices in group.groupby("bootstrap_cluster").groups.items()
        }
        clusters = np.array(list(cluster_to_indices), dtype=object)
        for replicate in range(N_BOOT):
            sampled_clusters = rng.choice(clusters, size=len(clusters), replace=True)
            sampled_indices = np.concatenate([cluster_to_indices[c] for c in sampled_clusters])
            sample = group.iloc[sampled_indices].copy()
            native_n = max(1, int(round(MATCHED_COVERAGE * len(sample))))
            populations = {
                "all_identity_clean": sample,
                "v0_4_accepted_mask": sample.loc[sample["v0_4_accepted"]],
                "v0_4_rejected_forced": sample.loc[~sample["v0_4_accepted"]],
                "native_margin_matched": sample.nlargest(native_n, "native_margin"),
            }
            for population, frame in populations.items():
                metrics = calculate_metrics(frame)
                rows.append(
                    {
                        "comparator_id": comparator_id,
                        "replicate": replicate,
                        "population": population,
                        "n": metrics["n"],
                        "balanced_accuracy": metrics["balanced_accuracy"],
                        "mcc": metrics["mcc"],
                        "roc_auc": metrics["roc_auc"],
                        "spearman_rho": metrics["spearman_rho"],
                    }
                )
    bootstrap = pd.DataFrame(rows)
    bootstrap.to_csv(BOOTSTRAP, index=False)
    return bootstrap


def summarize(metrics: pd.DataFrame, coverage: pd.DataFrame, bootstrap: pd.DataFrame, overlap: dict) -> dict:
    key_metrics = metrics.loc[
        metrics["population"].isin(
            ["all_identity_clean", "v0_4_accepted_mask", "native_margin_matched_453"]
        )
    ].copy()
    intervals: list[dict] = []
    for (comparator_id, population), group in bootstrap.groupby(["comparator_id", "population"]):
        if population not in {"all_identity_clean", "v0_4_accepted_mask", "native_margin_matched"}:
            continue
        row = {"comparator_id": comparator_id, "population": population}
        for metric in ["balanced_accuracy", "mcc", "roc_auc", "spearman_rho"]:
            values = group[metric].dropna().to_numpy()
            row[f"{metric}_ci95_low"] = float(np.quantile(values, 0.025))
            row[f"{metric}_ci95_high"] = float(np.quantile(values, 0.975))
        intervals.append(row)

    pivot = metrics.pivot_table(
        index="comparator_id", columns="population", values=["balanced_accuracy", "mcc", "spearman_rho"]
    )
    transfer = {}
    for comparator_id in pivot.index:
        transfer[comparator_id] = {
            "accepted_minus_all_balanced_accuracy": float(
                pivot.loc[comparator_id, ("balanced_accuracy", "v0_4_accepted_mask")]
                - pivot.loc[comparator_id, ("balanced_accuracy", "all_identity_clean")]
            ),
            "accepted_minus_native_margin_matched_balanced_accuracy": float(
                pivot.loc[comparator_id, ("balanced_accuracy", "v0_4_accepted_mask")]
                - pivot.loc[comparator_id, ("balanced_accuracy", "native_margin_matched_453")]
            ),
        }

    summary = {
        "version": "ABCB1_PUBLIC_PREDICTOR_BENCHMARK_SUMMARY_V1",
        "status": "completed_available_comparators_admetlab_pending_server_failure",
        "protocol": str(PROTOCOL.relative_to(ROOT)),
        "protocol_sha256": sha256(PROTOCOL),
        "n_benchmark": 2624,
        "n_bootstrap": N_BOOT,
        "deeppk_overlap_audit": overlap,
        "key_metrics": key_metrics.replace({np.nan: None}).to_dict(orient="records"),
        "bootstrap_intervals": intervals,
        "v0_4_mask_transfer_differences": transfer,
        "admetlab_status": {
            "status": "pending_official_server_failure",
            "observed_error": "HTTP 500 KeyError ['BSEP'] not in index at /api/single/admet; documented /api/admet returned 404",
        },
        "claim_boundary": (
            "Post-hoc public-predictor comparison. AZ training identities are proprietary and not auditable; "
            "ADMET-AI is an endpoint-mismatch inhibition negative control; ADMETlab is not imputed."
        ),
        "outputs": {},
    }
    for path in [PREDICTIONS, DEEP_PK_PARSED, OVERLAP, METRICS, COVERAGE, BOOTSTRAP]:
        summary["outputs"][str(path.relative_to(ROOT))] = sha256(path)
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    return summary


def main() -> None:
    long, overlap = build_long_predictions()
    metrics = analyze_metrics(long)
    coverage = analyze_coverage(long)
    bootstrap = analyze_bootstrap(long)
    summary = summarize(metrics, coverage, bootstrap, overlap)
    print(
        json.dumps(
            {
                "status": summary["status"],
                "deeppk_overlap_audit": overlap,
                "outputs": summary["outputs"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
