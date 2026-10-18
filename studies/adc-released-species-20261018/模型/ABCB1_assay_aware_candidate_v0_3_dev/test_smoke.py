#!/usr/bin/env python3
"""Minimal behavioral tests for the assay-aware candidate CLI core."""

from pathlib import Path

import joblib

from predict_abcb1_assay_aware import predict_one


HERE = Path(__file__).resolve().parent
package = joblib.load(HERE / "ABCB1_assay_aware_candidate_v0_3_dev.joblib")

blank = predict_one("", package)
assert blank["triage_call"] == "unresolved_structure"

methane = predict_one("C", package)
assert methane["triage_call"] == "abstain_outside_direct_panel_domain"

acrivastine = predict_one("CC1=CC=C(C=C1)/C(=C\\CN2CCCC2)/C3=CC=CC(=N3)/C=C/C(=O)O", package)
assert acrivastine["exact_development_identity_match"] == "yes"
assert acrivastine["triage_call"] == "training_identity_reference_only_not_validation"

assert package["status"] == "single_source_internal_development_candidate_not_externally_validated"
print("PASS: assay-aware v0.3 development candidate smoke tests")
