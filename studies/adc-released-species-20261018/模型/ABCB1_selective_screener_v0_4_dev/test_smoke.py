#!/usr/bin/env python3
"""Operational and boundary smoke tests for v0.4-dev."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(HERE))

from predict_abcb1_selective import load_packages, predict_one  # noqa: E402


SELECTIVE = HERE / "ABCB1_selective_uncertainty_v0_4_dev.joblib"
PRIMARY = PROJECT / "模型" / "ABCB1_screener_v0_1" / "ABCB1_screener_v0_1.joblib"
DIRECT = PROJECT / "模型" / "ABCB1_assay_aware_candidate_v0_3_dev" / "ABCB1_assay_aware_candidate_v0_3_dev.joblib"
SUMMARY = PROJECT / "数据" / "processed" / "abcb1_selective_screener_v0_4_build_summary_v1.json"
DIAGNOSTIC = PROJECT / "数据" / "processed" / "abcb1_selective_v0_4_posthoc_diagnostic_summary_v1.json"

LOW_EXAMPLE = "CC(=O)Nc1ccc(O)cc1"  # acetaminophen; historical diagnostic, not a tuning label
HIGH_EXAMPLE = "O=C(NCc1c(F)cc(F)cc1F)c1cn2c(c(O)c1=O)C(=O)N1[C@H]3CC[C@H](C3)O[C@@H]1C2"  # bictegravir
ABSTAIN_EXAMPLE = "CCO"


class SelectiveScreenerSmokeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.selective, cls.primary, cls.direct = load_packages(SELECTIVE, PRIMARY, DIRECT)

    def score(self, smiles: str) -> dict:
        return predict_one(smiles, self.selective, self.primary, self.direct)

    def test_all_three_operational_outputs_are_reachable(self) -> None:
        self.assertEqual(self.score(HIGH_EXAMPLE)["risk_call"], "high_ABCB1_efflux_risk")
        self.assertEqual(self.score(LOW_EXAMPLE)["risk_call"], "low_ABCB1_efflux_risk")
        self.assertEqual(self.score(ABSTAIN_EXAMPLE)["risk_call"], "abstain")

    def test_accepted_calls_have_no_abstention_reason(self) -> None:
        for smiles in [HIGH_EXAMPLE, LOW_EXAMPLE]:
            result = self.score(smiles)
            self.assertEqual(result["accepted_prediction"], "yes")
            self.assertEqual(result["abstention_reasons"], [])

    def test_abstention_is_explained(self) -> None:
        result = self.score(ABSTAIN_EXAMPLE)
        self.assertEqual(result["accepted_prediction"], "no")
        self.assertIn("probability_gray_zone", result["abstention_reasons"])
        self.assertEqual(len(result["abstention_reasons"]), len(result["abstention_reasons_zh"]))

    def test_invalid_structure_abstains(self) -> None:
        result = self.score("not_a_smiles")
        self.assertEqual(result["risk_call"], "abstain")
        self.assertEqual(result["abstention_reasons"], ["unresolved_structure"])

    def test_incompatible_object_and_context_abstain(self) -> None:
        intact = predict_one(
            HIGH_EXAMPLE,
            self.selective,
            self.primary,
            self.direct,
            chemical_object_type="intact_ADC",
        )
        ddi = predict_one(
            HIGH_EXAMPLE,
            self.selective,
            self.primary,
            self.direct,
            requested_context="clinical_DDI",
        )
        self.assertEqual(intact["risk_call"], "abstain")
        self.assertIn("object_scope_incompatible", intact["abstention_reasons"])
        self.assertEqual(ddi["risk_call"], "abstain")
        self.assertIn("requested_context_incompatible", ddi["abstention_reasons"])

    def test_training_identity_is_not_a_prospective_call(self) -> None:
        norethindrone = "C#C[C@]1(O)CC[C@H]2[C@@H]3CC=C4CC(=O)CCC4[C@H]3CC[C@@]21C"
        result = self.score(norethindrone)
        self.assertEqual(result["risk_call"], "abstain")
        self.assertIn("reference_identity_not_prospective", result["abstention_reasons"])

    def test_frozen_thresholds(self) -> None:
        self.assertEqual(self.selective["primary_probability_low_threshold"], 0.35)
        self.assertEqual(self.selective["primary_probability_high_threshold"], 0.70)
        self.assertEqual(len(self.selective["bootstrap_models"]), 64)

    def test_prediction_is_deterministic(self) -> None:
        first = json.dumps(self.score(HIGH_EXAMPLE), ensure_ascii=False, sort_keys=True)
        second = json.dumps(self.score(HIGH_EXAMPLE), ensure_ascii=False, sort_keys=True)
        self.assertEqual(first, second)

    def test_batch_cli_preserves_rows_and_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            input_path = Path(temporary) / "queries.csv"
            output_path = Path(temporary) / "predictions.csv"
            with input_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["name", "smiles"])
                writer.writeheader()
                writer.writerows(
                    [
                        {"name": "high", "smiles": HIGH_EXAMPLE},
                        {"name": "low", "smiles": LOW_EXAMPLE},
                        {"name": "abstain", "smiles": ABSTAIN_EXAMPLE},
                    ]
                )
            subprocess.run(
                [
                    sys.executable,
                    str(HERE / "predict_abcb1_selective.py"),
                    "--input-csv",
                    str(input_path),
                    "--output",
                    str(output_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            with output_path.open("r", encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 3)
            self.assertEqual(
                [row["risk_call"] for row in rows],
                ["high_ABCB1_efflux_risk", "low_ABCB1_efflux_risk", "abstain"],
            )

    def test_build_and_diagnostic_boundaries_are_explicit(self) -> None:
        build = json.loads(SUMMARY.read_text(encoding="utf-8"))
        diagnostic = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
        self.assertFalse(build["historical_external_outcomes_used_for_rule_selection"])
        self.assertTrue(build["future_source_confirmation_required"])
        self.assertEqual(
            diagnostic["diagnostic_status"],
            "historical_posthoc_diagnostic_not_external_validation",
        )
        self.assertFalse(diagnostic["rule_tuned_on_historical_external_outcomes"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
