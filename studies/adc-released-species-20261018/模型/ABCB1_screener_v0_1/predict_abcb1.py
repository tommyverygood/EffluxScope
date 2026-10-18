#!/usr/bin/env python3
"""Command-line interface for the ABCB1 research screener v0.1."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
SCRIPT_DIR = PROJECT / "数据" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from abcb1_screener_core_v0_1 import load_package, predict_smiles  # noqa: E402


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Predict ABCB1 structural substrate propensity with explicit applicability-domain flags."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--smiles", help="one SMILES string")
    source.add_argument("--input-csv", type=Path, help="batch CSV input")
    parser.add_argument("--smiles-column", default="smiles", help="SMILES column for --input-csv")
    parser.add_argument("--output", type=Path, help="batch CSV output path")
    parser.add_argument(
        "--model",
        type=Path,
        default=HERE / "ABCB1_screener_v0_1.joblib",
        help="model package path",
    )
    args = parser.parse_args()
    package = load_package(args.model)

    if args.smiles is not None:
        print(json.dumps(predict_smiles(args.smiles, package), ensure_ascii=False, indent=2))
        return

    if args.output is None:
        parser.error("--output is required with --input-csv")
    rows = read_csv(args.input_csv)
    if rows and args.smiles_column not in rows[0]:
        parser.error(f"SMILES column not found: {args.smiles_column}")
    output_rows = []
    for row in rows:
        prediction = predict_smiles(row.get(args.smiles_column, ""), package)
        output_rows.append({**row, **prediction})
    write_csv(args.output, output_rows)
    print(json.dumps({"input_rows": len(rows), "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
