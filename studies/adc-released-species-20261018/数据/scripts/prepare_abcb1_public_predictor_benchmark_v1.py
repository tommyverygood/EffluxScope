#!/usr/bin/env python3
"""Prepare the frozen identity-clean Biogen interface for public predictors."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "数据" / "processed"
SOURCE = PROCESSED / "abcb1_quantity_biogen_analysis_rows_v1.csv"
OUTPUT = PROCESSED / "abcb1_public_predictor_benchmark_input_v1.csv"
MANIFEST = PROCESSED / "abcb1_public_predictor_benchmark_input_manifest_v1.json"
DEEP_PK_REQUESTS = [
    PROCESSED / "abcb1_public_predictor_deeppk_request_batch1_v1.csv",
    PROCESSED / "abcb1_public_predictor_deeppk_request_batch2_v1.csv",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    data = pd.read_csv(SOURCE)
    data = data.loc[data["identity_clean"].eq("yes") & data["model_usable"].eq("yes")].copy()
    data = data.sort_values("external_row_id", kind="stable").reset_index(drop=True)

    columns = {
        "external_row_id": "external_row_id",
        "standardized_isomeric_smiles": "smiles",
        "standardized_inchikey": "standardized_inchikey",
        "log10_mdr1_mdck_er": "observed_log10_mdr1_mdck_er",
        "mdr1_mdck_er": "observed_mdr1_mdck_er",
        "frozen_er2_binary_label": "observed_er2_label",
        "bootstrap_cluster": "bootstrap_cluster",
        "murcko_scaffold_smiles": "murcko_scaffold_smiles",
        "accepted_prediction": "v0_4_accepted",
        "risk_call": "v0_4_risk_call",
        "primary_calibrated_probability": "v0_4_primary_score",
    }
    output = data[list(columns)].rename(columns=columns)
    output["v0_4_accepted"] = output["v0_4_accepted"].eq("yes")

    if len(output) != 2624:
        raise ValueError(f"Expected 2,624 identity-clean rows, found {len(output):,}")
    if output["external_row_id"].duplicated().any():
        raise ValueError("external_row_id must be unique")
    if output["smiles"].isna().any():
        raise ValueError("All public-predictor rows require a standardized SMILES")
    if output["observed_er2_label"].isna().any():
        raise ValueError("All benchmark rows require the frozen ER>=2 label")

    output.to_csv(OUTPUT, index=False)
    midpoint = len(output) // 2
    output.loc[: midpoint - 1, ["smiles"]].to_csv(DEEP_PK_REQUESTS[0], index=False)
    output.loc[midpoint:, ["smiles"]].to_csv(DEEP_PK_REQUESTS[1], index=False)
    manifest = {
        "version": "ABCB1_PUBLIC_PREDICTOR_BENCHMARK_INPUT_V1",
        "status": "frozen_posthoc_comparator_interface",
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": sha256(SOURCE),
        "output": str(OUTPUT.relative_to(ROOT)),
        "output_sha256": sha256(OUTPUT),
        "deeppk_request_batches": [
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "n_rows": int(sum(1 for _ in path.open()) - 1),
            }
            for path in DEEP_PK_REQUESTS
        ],
        "n_rows": int(len(output)),
        "n_v0_4_accepted": int(output["v0_4_accepted"].sum()),
        "n_v0_4_rejected": int((~output["v0_4_accepted"]).sum()),
        "observed_positive_prevalence": float(output["observed_er2_label"].mean()),
        "claim_boundary": (
            "Biogen outcomes were already observed before this public-predictor comparison; "
            "the interface is frozen for post-hoc comparator analysis, not prospective validation."
        ),
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
