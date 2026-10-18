#!/usr/bin/env python3
"""Freeze label-free structural-shift assignments for the ABCB1 OOD audit.

This stage reads structures and previously frozen model artifacts, but it never
uses MDR1-MDCK outcomes or prediction correctness.  Its output is hashed and
locked before the separate metric-generation stage is allowed to run.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import rdkit
import sklearn
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from sklearn.cluster import MiniBatchKMeans


PROJECT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT / "数据" / "processed"

PROTOCOL = PROCESSED / "abcb1_ood_structural_shift_protocol_v1.json"
PRIMARY_REFERENCE = PROCESSED / "abcb1_screener_v0_1_training_audit.csv"
DIRECT_EVIDENCE_PANEL = PROCESSED / "soskuti_2024_mdr1_species_panel_v1.csv"
DIRECT_MODEL = PROJECT / "模型" / "ABCB1_assay_aware_candidate_v0_3_dev" / "ABCB1_assay_aware_candidate_v0_3_dev.joblib"
BIOGEN_ANALYSIS = PROCESSED / "abcb1_quantity_biogen_analysis_rows_v1.csv"
SELECTIVE_MODEL = PROJECT / "模型" / "ABCB1_selective_screener_v0_4_dev" / "ABCB1_selective_uncertainty_v0_4_dev.joblib"
OOD_PAPER = PROJECT / "论文参考" / "ci5c00475.pdf"

STRUCTURE_ROWS = PROCESSED / "abcb1_ood_structural_assignments_v1.csv"
REFERENCE_LOO = PROCESSED / "abcb1_ood_reference_loo_distances_v1.csv"
STRUCTURE_LOCK = PROCESSED / "abcb1_ood_structure_assignment_lock_v1.json"

FP_RADIUS = 2
FP_BITS = 2048
TOP_K = 5
LOHI_SIMILARITY_CUTOFF = 0.40
BUTINA_DISTANCE_CUTOFF = 0.60
KMEANS_CLUSTERS = 20
KMEANS_SEED = 20260807

FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(
    radius=FP_RADIUS,
    fpSize=FP_BITS,
    includeChirality=False,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_connectivity(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f"unparseable SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def connectivity_block(inchikey: object, smiles: str) -> str:
    value = "" if pd.isna(inchikey) else str(inchikey).strip()
    if value:
        return value.split("-", 1)[0]
    return "SMILES::" + canonical_connectivity(smiles)


def make_reference_frame(rows: list[dict[str, Any]], pool: str) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    required = {"reference_id", "reference_name", "standardized_smiles", "standardized_inchikey"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"missing reference fields for {pool}: {sorted(missing)}")
    frame["reference_pool"] = pool
    frame["connectivity_block"] = [
        connectivity_block(key, smiles)
        for key, smiles in zip(frame["standardized_inchikey"], frame["standardized_smiles"])
    ]
    frame["mol"] = frame["standardized_smiles"].map(Chem.MolFromSmiles)
    if frame["mol"].isna().any():
        raise RuntimeError(f"unparseable structure in {pool}")
    frame["fp"] = frame["mol"].map(FP_GENERATOR.GetFingerprint)
    return frame


def deduplicate_connectivity(frame: pd.DataFrame, pool: str) -> pd.DataFrame:
    work = frame.drop_duplicates("connectivity_block", keep="first").copy().reset_index(drop=True)
    work["reference_pool"] = pool
    return work


def top_similarity(query_fp: Any, reference_fps: list[Any], k: int = TOP_K) -> tuple[float, float, list[float]]:
    similarities = np.asarray(DataStructs.BulkTanimotoSimilarity(query_fp, reference_fps), dtype=float)
    if similarities.size < k:
        raise RuntimeError(f"reference pool has fewer than {k} structures")
    top = np.sort(similarities)[::-1][:k]
    return float(top[0]), float(np.mean(top)), [float(x) for x in top]


def leave_one_connectivity_out(frame: pd.DataFrame, pool: str) -> list[dict[str, Any]]:
    fps = frame["fp"].tolist()
    blocks = frame["connectivity_block"].tolist()
    output: list[dict[str, Any]] = []
    for index, query_fp in enumerate(fps):
        keep = [j for j, block in enumerate(blocks) if block != blocks[index]]
        if len(keep) < TOP_K:
            raise RuntimeError(f"too few leave-one-connectivity-out references in {pool}")
        s1, s5, top = top_similarity(query_fp, [fps[j] for j in keep])
        output.append(
            {
                "reference_pool": pool,
                "reference_id": frame.iloc[index]["reference_id"],
                "connectivity_block": blocks[index],
                "loo_s1": s1,
                "loo_s5": s5,
                "loo_d1": 1.0 - s1,
                "loo_d5": 1.0 - s5,
                "loo_top5_similarities_json": json.dumps(top, separators=(",", ":")),
            }
        )
    return output


def d5_thresholds(loo_rows: list[dict[str, Any]], pool: str) -> dict[str, float]:
    values = np.asarray([row["loo_d5"] for row in loo_rows if row["reference_pool"] == pool], dtype=float)
    if values.size < 20:
        raise RuntimeError(f"too few LOO rows for {pool}")
    return {
        "q50": float(np.quantile(values, 0.50)),
        "q80": float(np.quantile(values, 0.80)),
        "q95": float(np.quantile(values, 0.95)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
    }


def distance_stratum(d5: float, thresholds: dict[str, float]) -> str:
    if d5 <= thresholds["q50"]:
        return "ID_like"
    if d5 <= thresholds["q80"]:
        return "transition"
    if d5 <= thresholds["q95"]:
        return "OOD"
    return "extreme_OOD"


def fp_array(fp: Any) -> np.ndarray:
    array = np.zeros((FP_BITS,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, array)
    return array


def butina_assignments(fps: list[Any]) -> tuple[np.ndarray, list[tuple[int, ...]]]:
    distances: list[float] = []
    for index in range(1, len(fps)):
        similarities = DataStructs.BulkTanimotoSimilarity(fps[index], fps[:index])
        distances.extend(1.0 - float(value) for value in similarities)
    clusters = list(
        Butina.ClusterData(
            distances,
            len(fps),
            BUTINA_DISTANCE_CUTOFF,
            isDistData=True,
            reordering=True,
        )
    )
    assignments = np.full(len(fps), -1, dtype=int)
    for cluster_id, members in enumerate(clusters, start=1):
        assignments[list(members)] = cluster_id
    if np.any(assignments < 0):
        raise RuntimeError("Butina assignment left unassigned rows")
    return assignments, clusters


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if protocol.get("status") != "locked_before_structure_only_assignment_and_metric_generation":
        raise RuntimeError("OOD protocol is not in the required locked state")

    path_map = {
        "数据/processed/abcb1_screener_v0_1_training_audit.csv": PRIMARY_REFERENCE,
        "数据/processed/soskuti_2024_mdr1_species_panel_v1.csv": DIRECT_EVIDENCE_PANEL,
        "模型/ABCB1_assay_aware_candidate_v0_3_dev/ABCB1_assay_aware_candidate_v0_3_dev.joblib": DIRECT_MODEL,
        "数据/processed/abcb1_quantity_biogen_analysis_rows_v1.csv": BIOGEN_ANALYSIS,
        "模型/ABCB1_selective_screener_v0_4_dev/ABCB1_selective_uncertainty_v0_4_dev.joblib": SELECTIVE_MODEL,
        "论文参考/ci5c00475.pdf": OOD_PAPER,
    }
    for relative, path in path_map.items():
        observed = sha256(path)
        expected = protocol["input_hashes_sha256"][relative]
        if observed != expected:
            raise RuntimeError(f"input hash mismatch for {relative}: {observed} != {expected}")

    primary_raw = pd.read_csv(PRIMARY_REFERENCE, low_memory=False)
    if len(primary_raw) != 443:
        raise RuntimeError(f"expected 443 primary rows, observed {len(primary_raw)}")
    primary_rows = [
        {
            "reference_id": row["interface_id"],
            "reference_name": row.get("preferred_name", ""),
            "standardized_smiles": row["derived_isomeric_smiles"],
            "standardized_inchikey": row["derived_standardized_inchikey"],
        }
        for _, row in primary_raw.iterrows()
    ]
    primary = make_reference_frame(primary_rows, "primary_model_443")

    direct_package = joblib.load(DIRECT_MODEL)
    direct_training_rows = direct_package.get("training_rows", [])
    if len(direct_training_rows) != 63:
        raise RuntimeError(f"expected 63 direct-model rows, observed {len(direct_training_rows)}")
    direct_model_rows = [
        {
            "reference_id": row["panel_entity_id"],
            "reference_name": row.get("entity_name", ""),
            "standardized_smiles": row["standardized_isomeric_smiles"],
            "standardized_inchikey": row["standardized_inchikey"],
        }
        for row in direct_training_rows
    ]
    direct_model = make_reference_frame(direct_model_rows, "direct_model_63")

    direct_evidence_raw = pd.read_csv(DIRECT_EVIDENCE_PANEL, low_memory=False)
    if len(direct_evidence_raw) != 84:
        raise RuntimeError(f"expected 84 direct-evidence rows, observed {len(direct_evidence_raw)}")
    direct_evidence_rows = [
        {
            "reference_id": row["panel_entity_id"],
            "reference_name": row.get("entity_name_source", ""),
            "standardized_smiles": row["standardized_isomeric_smiles"],
            "standardized_inchikey": row["standardized_inchikey"],
        }
        for _, row in direct_evidence_raw.iterrows()
    ]
    direct_evidence = make_reference_frame(direct_evidence_rows, "direct_evidence_84")

    joint_model = deduplicate_connectivity(
        pd.concat([primary, direct_model], ignore_index=True),
        "joint_model_connectivity_unique",
    )
    joint_evidence = deduplicate_connectivity(
        pd.concat([primary, direct_evidence], ignore_index=True),
        "joint_evidence_connectivity_unique",
    )

    pool_frames = {
        "primary_model_443": primary,
        "direct_model_63": direct_model,
        "joint_model_connectivity_unique": joint_model,
        "joint_evidence_connectivity_unique": joint_evidence,
    }
    loo_rows: list[dict[str, Any]] = []
    for pool_name, frame in pool_frames.items():
        loo_rows.extend(leave_one_connectivity_out(frame, pool_name))
    loo_frame = pd.DataFrame(loo_rows)
    loo_frame.to_csv(REFERENCE_LOO, index=False)
    thresholds = {pool: d5_thresholds(loo_rows, pool) for pool in pool_frames}

    biogen = pd.read_csv(BIOGEN_ANALYSIS, low_memory=False)
    structure_columns = [
        "external_row_id",
        "source_record_id",
        "standardized_isomeric_smiles",
        "standardized_inchikey",
        "connectivity_block",
        "murcko_scaffold_smiles",
        "bootstrap_cluster",
        "identity_clean",
        "scaffold_clean",
    ]
    missing_columns = set(structure_columns) - set(biogen.columns)
    if missing_columns:
        raise RuntimeError(f"Biogen structural columns missing: {sorted(missing_columns)}")
    work = biogen.loc[biogen["identity_clean"].eq("yes"), structure_columns].copy().reset_index(drop=True)
    if len(work) != 2624:
        raise RuntimeError(f"expected 2624 identity-clean Biogen structures, observed {len(work)}")
    if work["external_row_id"].duplicated().any():
        raise RuntimeError("duplicate external_row_id in identity-clean Biogen structures")

    work["mol"] = work["standardized_isomeric_smiles"].map(Chem.MolFromSmiles)
    if work["mol"].isna().any():
        raise RuntimeError("unparseable identity-clean Biogen structure")
    work["fp"] = work["mol"].map(FP_GENERATOR.GetFingerprint)

    pool_fps = {name: frame["fp"].tolist() for name, frame in pool_frames.items()}
    prefix_map = {
        "primary_model_443": "primary443",
        "direct_model_63": "direct63",
        "joint_model_connectivity_unique": "joint_model",
        "joint_evidence_connectivity_unique": "joint_evidence",
    }
    similarity_columns: dict[str, list[Any]] = {}
    for pool_name, fps in pool_fps.items():
        prefix = prefix_map[pool_name]
        s1_values: list[float] = []
        s5_values: list[float] = []
        top5_json: list[str] = []
        for query_fp in work["fp"]:
            s1, s5, top = top_similarity(query_fp, fps)
            s1_values.append(s1)
            s5_values.append(s5)
            top5_json.append(json.dumps(top, separators=(",", ":")))
        similarity_columns[f"{prefix}_s1"] = s1_values
        similarity_columns[f"{prefix}_s5"] = s5_values
        similarity_columns[f"{prefix}_d1"] = [1.0 - value for value in s1_values]
        similarity_columns[f"{prefix}_d5"] = [1.0 - value for value in s5_values]
        similarity_columns[f"{prefix}_top5_similarities_json"] = top5_json
    for column, values in similarity_columns.items():
        work[column] = values

    joint_thresholds = thresholds["joint_model_connectivity_unique"]
    primary_thresholds = thresholds["primary_model_443"]
    direct_thresholds = thresholds["direct_model_63"]
    work["reference_adaptive_stratum"] = work["joint_model_d5"].map(
        lambda value: distance_stratum(float(value), joint_thresholds)
    )
    work["strict_lohi_ood"] = np.where(work["joint_model_s1"] <= LOHI_SIMILARITY_CUTOFF, "yes", "no")
    work["dual_ood"] = np.where(
        (work["primary443_d5"] > primary_thresholds["q80"])
        & (work["direct63_d5"] > direct_thresholds["q80"]),
        "yes",
        "no",
    )

    ranked = work["joint_model_d5"].rank(method="first")
    work["external_d5_quintile"] = pd.qcut(
        ranked,
        5,
        labels=["Q1_nearest", "Q2", "Q3", "Q4", "Q5_farthest"],
    ).astype(str)

    fp_matrix = np.vstack([fp_array(fp) for fp in work["fp"]]).astype(np.float32)
    kmeans = MiniBatchKMeans(
        n_clusters=KMEANS_CLUSTERS,
        random_state=KMEANS_SEED,
        batch_size=512,
        n_init=20,
        max_iter=300,
    )
    work["kmeans_cluster_id"] = kmeans.fit_predict(fp_matrix).astype(int)
    cluster_distance = (
        work.groupby("kmeans_cluster_id", sort=True)["joint_model_d5"]
        .agg(["median", "size"])
        .sort_values(["median", "size"], ascending=[False, False])
    )
    selected_clusters: list[int] = []
    selected_n = 0
    target_n = int(math.ceil(0.20 * len(work)))
    for cluster_id, row in cluster_distance.iterrows():
        selected_clusters.append(int(cluster_id))
        selected_n += int(row["size"])
        if selected_n >= target_n:
            break
    work["kmeans_ood20"] = np.where(work["kmeans_cluster_id"].isin(selected_clusters), "yes", "no")

    butina_ids, butina_clusters = butina_assignments(work["fp"].tolist())
    work["butina_cluster_id"] = butina_ids

    output_columns = [column for column in work.columns if column not in {"mol", "fp"}]
    output = work[output_columns].copy()
    output.to_csv(STRUCTURE_ROWS, index=False)

    selective_package = joblib.load(SELECTIVE_MODEL)
    if selective_package.get("model_version") != "ABCB1_selective_screener_v0.4_dev":
        raise RuntimeError("unexpected selective model version")

    lock = {
        "lock_id": "ABCB1_OOD_STRUCTURE_ASSIGNMENT_LOCK_V1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "structure_only_assignments_frozen_before_new_metric_generation",
        "protocol_path": str(PROTOCOL.relative_to(PROJECT)),
        "protocol_sha256": sha256(PROTOCOL),
        "outcomes_or_correctness_used_for_assignment": False,
        "runtime": {
            "rdkit": rdkit.__version__,
            "scikit_learn": sklearn.__version__,
            "fingerprint": {"radius": FP_RADIUS, "n_bits": FP_BITS, "include_chirality": False},
        },
        "reference_counts": {
            "primary_model_rows": int(len(primary)),
            "direct_model_rows": int(len(direct_model)),
            "direct_evidence_rows": int(len(direct_evidence)),
            "joint_model_connectivity_unique": int(len(joint_model)),
            "joint_evidence_connectivity_unique": int(len(joint_evidence)),
        },
        "reference_d5_thresholds": thresholds,
        "external_counts": {
            "identity_clean": int(len(output)),
            "scaffold_clean": int(output["scaffold_clean"].eq("yes").sum()),
            "reference_adaptive_strata": output["reference_adaptive_stratum"].value_counts().sort_index().to_dict(),
            "strict_lohi_ood": int(output["strict_lohi_ood"].eq("yes").sum()),
            "dual_ood": int(output["dual_ood"].eq("yes").sum()),
            "kmeans_ood20": int(output["kmeans_ood20"].eq("yes").sum()),
            "kmeans_ood20_clusters": selected_clusters,
            "kmeans_total_clusters": int(output["kmeans_cluster_id"].nunique()),
            "butina_total_clusters": int(len(butina_clusters)),
            "butina_singletons": int(sum(len(cluster) == 1 for cluster in butina_clusters)),
        },
        "outputs": {
            str(STRUCTURE_ROWS.relative_to(PROJECT)): sha256(STRUCTURE_ROWS),
            str(REFERENCE_LOO.relative_to(PROJECT)): sha256(REFERENCE_LOO),
        },
        "interpretation_boundary": (
            "Label-free structural assignments for a retrospective post-outcome audit; "
            "not a prospective OOD validation and not permission to retune v0.4."
        ),
    }
    STRUCTURE_LOCK.write_text(json.dumps(lock, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(lock, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
