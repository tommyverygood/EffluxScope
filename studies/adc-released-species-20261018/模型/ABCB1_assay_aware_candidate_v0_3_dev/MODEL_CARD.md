# ABCB1 assay-aware candidate v0.3-dev

## Intended use

This package ranks standardized free/released small-molecule structures for a narrowly defined human MDR1 matched-parental screening context. It returns a probability, a binary development call, empirical-domain status, and the nearest direct-panel reference.

It is not validated for intact ADCs, intracellular free exposure, resistance, efficacy, clinical drug interactions, or a universal ABCB1 substrate label.

## Development evidence

- Source: Sóskuti et al., *Pharmaceutics* 2024, Supplementary Table S1.
- Context: hMDR1 overexpressing Abcb1KO-MDCKII versus matched Mock, 1 µM, 120 min.
- Source-QC exclusions: 12 low-recovery compounds and the famotidine permeability control.
- Primary label: approximate 95% interval for net ER entirely above or below 2; 63 entities (41 positive-like, 22 negative-like).
- Model: Morgan fingerprint plus eight physicochemical descriptors, class-balanced logistic regression.
- Internal evaluation: 10×5 repeated scaffold-grouped outer CV with nested regularization selection.

Primary internal OOF results: MCC 0.511 (scaffold-cluster bootstrap 95% CI 0.297–0.710), balanced accuracy 0.756, ROC AUC 0.849, and Brier score 0.162.

These are single-source internal-development estimates and are not external validation.

## Interpretation

- `training_identity_reference_only_not_validation`: exact identity is in the 63-entity development panel; the score is not a validation result.
- `candidate_positive_like` / `candidate_negative_like`: query is within the direct-panel empirical domain.
- `manual_review_low_direct_panel_support`: similarity lies between the frozen q10 and q25 domain thresholds.
- `abstain_outside_direct_panel_domain`: similarity is below the frozen q10 threshold.

## Required next gate

Freeze and score a new untouched human ABCB1 matched-parental source or a prospective assay panel. Do not tune this package on the previously used 2026-08-03 development external benchmark.
