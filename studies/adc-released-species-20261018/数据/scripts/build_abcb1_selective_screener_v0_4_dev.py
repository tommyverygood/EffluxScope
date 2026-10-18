#!/usr/bin/env python3
"""Build the frozen-rule ABCB1 selective research screener v0.4-dev.

The biological output remains a two-direction risk signal. ``abstain`` is a
decision state produced by uncertainty, conflict, identity, and applicability
checks; it is never fitted as a third biological class.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
from pathlib import Path

import joblib
import numpy as np
import rdkit
import scipy
import sklearn
from rdkit import Chem
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from abcb1_screener_core_v0_1 import descriptor_array, morgan_array


DATA = Path(__file__).resolve().parents[1]
PROJECT = DATA.parent
PROCESSED = DATA / "processed"
MODEL_DIR = PROJECT / "模型" / "ABCB1_selective_screener_v0_4_dev"
MODEL_DIR.mkdir(parents=True, exist_ok=True)

PROTOCOL = PROCESSED / "abcb1_selective_screener_v0_4_protocol_lock_v1.json"
PRIMARY_MODEL = PROJECT / "模型" / "ABCB1_screener_v0_1" / "ABCB1_screener_v0_1.joblib"
DIRECT_MODEL = PROJECT / "模型" / "ABCB1_assay_aware_candidate_v0_3_dev" / "ABCB1_assay_aware_candidate_v0_3_dev.joblib"
PRIMARY_AUDIT = PROCESSED / "abcb1_screener_v0_1_training_audit.csv"
PRIMARY_OOF = PROCESSED / "abcb1_screener_v0_1_oof_predictions.csv"
DIRECT_PANEL = PROCESSED / "soskuti_2024_mdr1_species_panel_v1.csv"

PACKAGE_PATH = MODEL_DIR / "ABCB1_selective_uncertainty_v0_4_dev.joblib"
INTERNAL_METRICS = PROCESSED / "abcb1_selective_internal_metrics_v1.csv"
THRESHOLD_SWEEP = PROCESSED / "abcb1_selective_probability_threshold_sweep_v1.csv"
SUMMARY = PROCESSED / "abcb1_selective_screener_v0_4_build_summary_v1.json"

MODEL_VERSION = "ABCB1_selective_screener_v0.4_dev"
RANDOM_SEED = 20260803
N_BOOTSTRAP_ESTIMATORS = 64
LOW_THRESHOLD = 0.35
HIGH_THRESHOLD = 0.70
LOCAL_K = 5
LOCAL_LOW = 0.35
LOCAL_HIGH = 0.65

EXPECTED_HASHES = {
    PRIMARY_MODEL: "fbc92dfa99583386ba124141fed5e77c6ecf8a59850301a6533db3883b298ac3",
    DIRECT_MODEL: "ebe6d1a8d5f06afc5c2d0b306dfd5266af56fd28fefd130147d20c5e015b907a",
    PRIMARY_AUDIT: "54979f1037817d6a218d31301c59a1056b0123747f1aacbbf2e12ddc85c8d747",
    PRIMARY_OOF: "330d1124886229e4f6cd22988e1e438a3250a208339a3336b49357edbef5e7b4",
    DIRECT_PANEL: "9892cd2bb4be3ff0c834113c5ec7b8d9b7752bc69667782cd5980f19f2ac3bf7",
}

DESCRIPTOR_NAMES = [
    "molecular_weight",
    "xlogp",
    "tpsa",
    "hbond_donors",
    "hbond_acceptors",
    "rotatable_bonds",
    "formal_charge",
    "heavy_atoms",
]


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
        raise ValueError(f"refusing to write an empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
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
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def verify_dependencies(protocol: dict) -> dict[str, str]:
    actual: dict[str, str] = {}
    for path, expected in EXPECTED_HASHES.items():
        observed = sha256_file(path)
        if observed != expected:
            raise RuntimeError(f"frozen dependency hash mismatch: {path}: {observed} != {expected}")
        actual[str(path.relative_to(PROJECT))] = observed
    if protocol.get("status") != "locked_before_selective_screener_build":
        raise RuntimeError("selective screener protocol is not in the required pre-build locked state")
    rule = protocol["frozen_decision_rule"]
    if float(rule["primary_low_probability_threshold"]) != LOW_THRESHOLD:
        raise RuntimeError("low threshold differs from lock")
    if float(rule["primary_high_probability_threshold"]) != HIGH_THRESHOLD:
        raise RuntimeError("high threshold differs from lock")
    if int(rule["bootstrap_ensemble"]["n_estimators"]) != N_BOOTSTRAP_ESTIMATORS:
        raise RuntimeError("bootstrap estimator count differs from lock")
    return actual


def mols_from_smiles(smiles_values: list[str], context: str) -> list[Chem.Mol]:
    mols: list[Chem.Mol] = []
    for index, smiles in enumerate(smiles_values):
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            raise RuntimeError(f"unparseable standardized SMILES in {context} row {index + 1}")
        mols.append(mol)
    return mols


def wilson_lower(correct: int, n: int, z: float = 1.96) -> float:
    if n == 0:
        return math.nan
    proportion = correct / n
    denominator = 1.0 + z * z / n
    center = proportion + z * z / (2.0 * n)
    spread = z * math.sqrt(proportion * (1.0 - proportion) / n + z * z / (4.0 * n * n))
    return (center - spread) / denominator


def selective_probability_metrics(y: np.ndarray, probability: np.ndarray, low: float, high: float) -> dict:
    high_mask = probability >= high
    low_mask = probability <= low
    accepted = high_mask | low_mask
    high_correct = int(y[high_mask].sum())
    low_correct = int((1 - y[low_mask]).sum())
    accepted_correct = high_correct + low_correct
    accepted_n = int(accepted.sum())
    high_n = int(high_mask.sum())
    low_n = int(low_mask.sum())
    return {
        "n": int(len(y)),
        "low_probability_threshold": low,
        "high_probability_threshold": high,
        "accepted_n": accepted_n,
        "abstain_n": int(len(y) - accepted_n),
        "coverage": accepted_n / len(y),
        "accepted_accuracy": accepted_correct / accepted_n if accepted_n else math.nan,
        "high_n": high_n,
        "high_precision": high_correct / high_n if high_n else math.nan,
        "high_precision_wilson95_lower": wilson_lower(high_correct, high_n),
        "low_n": low_n,
        "low_negative_predictive_value": low_correct / low_n if low_n else math.nan,
        "low_npv_wilson95_lower": wilson_lower(low_correct, low_n),
    }


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    dependency_hashes = verify_dependencies(protocol)
    primary_package = joblib.load(PRIMARY_MODEL)
    direct_package = joblib.load(DIRECT_MODEL)
    primary_rows = read_csv(PRIMARY_AUDIT)
    oof_rows = read_csv(PRIMARY_OOF)

    if len(primary_rows) != 443 or len(oof_rows) != 443:
        raise RuntimeError("expected 443 aligned primary rows and OOF rows")
    primary_keys = [row["derived_standardized_inchikey"] for row in primary_rows]
    if primary_keys != list(primary_package["training_structure_keys"]):
        raise RuntimeError("primary audit order does not match frozen v0.1 fingerprints")
    if len(direct_package["training_rows"]) != 63:
        raise RuntimeError("expected 63 direct-panel training rows")

    primary_mols = mols_from_smiles([row["derived_isomeric_smiles"] for row in primary_rows], "primary audit")
    direct_mols = mols_from_smiles(
        [row["standardized_isomeric_smiles"] for row in direct_package["training_rows"]],
        "direct package",
    )
    x_fp = np.stack([morgan_array(mol) for mol in primary_mols])
    x_desc = np.stack([descriptor_array(mol) for mol in primary_mols])
    y = np.asarray([int(row["target_label"]) for row in primary_rows], dtype=int)
    groups = np.asarray([row["exact_scaffold_group"] for row in primary_rows], dtype=object)
    unique_groups = np.unique(groups)

    rng = np.random.default_rng(RANDOM_SEED)
    bootstrap_models: list[dict] = []
    for bootstrap_index in range(N_BOOTSTRAP_ESTIMATORS):
        sampled_groups = rng.choice(unique_groups, size=len(unique_groups), replace=True)
        sampled_indices = np.concatenate([np.flatnonzero(groups == group) for group in sampled_groups])
        if len(np.unique(y[sampled_indices])) != 2:
            raise RuntimeError(f"bootstrap replicate {bootstrap_index} contains one class")
        scaler = StandardScaler().fit(x_desc[sampled_indices])
        x_train = np.hstack([x_fp[sampled_indices], scaler.transform(x_desc[sampled_indices])])
        estimator = LogisticRegression(
            C=0.1,
            class_weight="balanced",
            max_iter=5000,
            solver="liblinear",
            random_state=RANDOM_SEED + bootstrap_index,
        )
        estimator.fit(x_train, y[sampled_indices])
        bootstrap_models.append(
            {
                "replicate": bootstrap_index + 1,
                "random_seed": RANDOM_SEED + bootstrap_index,
                "descriptor_scaler": scaler,
                "estimator": estimator,
            }
        )

    combined_descriptors = np.vstack(
        [x_desc, np.stack([descriptor_array(mol) for mol in direct_mols])]
    )
    physchem_center = np.median(combined_descriptors, axis=0)
    q25 = np.quantile(combined_descriptors, 0.25, axis=0)
    q75 = np.quantile(combined_descriptors, 0.75, axis=0)
    physchem_scale = q75 - q25
    mad = np.median(np.abs(combined_descriptors - physchem_center), axis=0) * 1.4826
    physchem_scale = np.where(physchem_scale > 1e-12, physchem_scale, np.where(mad > 1e-12, mad, 1.0))
    physchem_scores = np.max(np.abs((combined_descriptors - physchem_center) / physchem_scale), axis=1)
    physchem_threshold = float(np.quantile(physchem_scores, 0.99))

    oof_y = np.asarray([int(row["target_label"]) for row in oof_rows], dtype=int)
    oof_probability = np.asarray(
        [float(row["repeated_exact_scaffold_oof_calibrated_probability"]) for row in oof_rows],
        dtype=float,
    )
    selected_metrics = selective_probability_metrics(oof_y, oof_probability, LOW_THRESHOLD, HIGH_THRESHOLD)
    internal_rows = [
        {
            "analysis_set": "v0_1_repeated_exact_Murcko_grouped_OOF",
            "rule_scope": "probability_margin_only_before_operational_rejectors",
            "metric": key,
            "value": round(float(value), 8) if isinstance(value, (float, np.floating)) and math.isfinite(value) else value,
            "claim_boundary": "single author-curated source interface; internal development only",
        }
        for key, value in selected_metrics.items()
    ]
    write_csv(INTERNAL_METRICS, internal_rows)

    sweep_rows: list[dict] = []
    for low in [0.20, 0.25, 0.30, 0.35, 0.40]:
        for high in [0.60, 0.65, 0.70, 0.75, 0.80, 0.85]:
            row = selective_probability_metrics(oof_y, oof_probability, low, high)
            sweep_rows.append(
                {
                    **{key: round(float(value), 8) if isinstance(value, (float, np.floating)) and math.isfinite(value) else value for key, value in row.items()},
                    "selected_frozen_rule": "yes" if low == LOW_THRESHOLD and high == HIGH_THRESHOLD else "no",
                    "selection_dataset": "v0_1 grouped OOF only",
                }
            )
    write_csv(THRESHOLD_SWEEP, sweep_rows)

    package = {
        "model_version": MODEL_VERSION,
        "build_date": "2026-08-03",
        "status": "development_selective_risk_screener_not_future_source_validated",
        "protocol_path": str(PROTOCOL.relative_to(PROJECT)),
        "protocol_sha256": sha256_file(PROTOCOL),
        "dependency_hashes": dependency_hashes,
        "primary_probability_low_threshold": LOW_THRESHOLD,
        "primary_probability_high_threshold": HIGH_THRESHOLD,
        "decision_threshold": 0.5,
        "bootstrap_probability_quantiles": [0.10, 0.50, 0.90],
        "bootstrap_resampling_unit": "exact_Murcko_scaffold_group",
        "bootstrap_random_seed": RANDOM_SEED,
        "bootstrap_models": bootstrap_models,
        "primary_training_labels": y.tolist(),
        "primary_reference_rows": [
            {
                "preferred_name": row["preferred_name"],
                "target_label": int(row["target_label"]),
                "standardized_inchikey": row["derived_standardized_inchikey"],
                "evidence_tier": row["evidence_tier"],
                "source_url": row["source_url"],
            }
            for row in primary_rows
        ],
        "direct_training_labels": [int(row["target_label"]) for row in direct_package["training_rows"]],
        "local_neighbor_k": LOCAL_K,
        "local_low_boundary": LOCAL_LOW,
        "local_high_boundary": LOCAL_HIGH,
        "primary_similarity_q10": float(primary_package["applicability_thresholds"]["self_nearest_q10"]),
        "direct_similarity_q10": float(direct_package["applicability_thresholds"]["self_nearest_q10"]),
        "physchem_descriptor_names": DESCRIPTOR_NAMES,
        "physchem_center": physchem_center,
        "physchem_scale": physchem_scale,
        "physchem_score_threshold_q99": physchem_threshold,
        "physchem_reference_n": int(combined_descriptors.shape[0]),
        "selected_internal_probability_margin_metrics": selected_metrics,
        "claim_boundary": (
            "ABCB1 efflux-risk triage for a standardized free/released small-molecule structure; "
            "not intact-ADC transport, intracellular unbound exposure, resistance, efficacy, clinical DDI, "
            "or universal assay-independent substrate truth"
        ),
    }
    joblib.dump(package, PACKAGE_PATH, compress=3)
    package_hash = sha256_file(PACKAGE_PATH)

    summary = {
        "model_version": MODEL_VERSION,
        "status": package["status"],
        "protocol_path": str(PROTOCOL.relative_to(PROJECT)),
        "protocol_sha256": package["protocol_sha256"],
        "package_path": str(PACKAGE_PATH.relative_to(PROJECT)),
        "package_sha256": package_hash,
        "frozen_dependency_hashes_verified": True,
        "bootstrap_estimators": N_BOOTSTRAP_ESTIMATORS,
        "primary_reference_n": len(primary_rows),
        "direct_reference_n": len(direct_package["training_rows"]),
        "physchem_reference_n": int(combined_descriptors.shape[0]),
        "physchem_score_threshold_q99": round(physchem_threshold, 8),
        "probability_margin_internal_metrics": selected_metrics,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scikit_learn": sklearn.__version__,
            "scipy": scipy.__version__,
            "rdkit": rdkit.__version__,
            "joblib": joblib.__version__,
        },
        "historical_external_outcomes_used_for_rule_selection": False,
        "future_source_confirmation_required": True,
        "claim_boundary": package["claim_boundary"],
    }
    dump_json(SUMMARY, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
