# Figure 3 inputs and locally computed comparator scores

The scripts preserve the historical `Figure1_endpoint_shift` filename; the current manuscript calls this Figure 3. Run `python run.py` in this component, or use the release-level verifier. No model is fitted by the plotting program.

`data/author_computed_az_predictions.csv` has 5,248 rows: two endpoint predictions for each of the same 2,624 external molecules. These rows are not 5,248 independent molecules. The corresponding checkpoint is identified in `provenance.json`; model weights are obtained separately from its official provider.

| Field | Meaning |
|---|---|
| external_row_id | Frozen Biogen row identifier, used for joining analysis records. |
| smiles, standardized_inchikey | Standardized structure and identity key; original source and MIT notice are retained. |
| observed_log10_mdr1_mdck_er | Measured MDR1-MDCK efflux ratio on the base-10 logarithmic scale. |
| observed_mdr1_mdck_er | Measured ratio on its original scale. |
| observed_er2_label | Binary measured endpoint at ER >= 2. |
| bootstrap_cluster | Frozen scaffold grouping used by the original analysis. |
| comparator_id, display_name, role, endpoint | Frozen model identity and assay endpoint description. |
| score | Locally computed comparator output on the log10 ER scale; not a calibrated probability. |
| binary_threshold | log10(2) = 0.3010299956639812 for both exported comparators. |
| predicted_label | Frozen binary call obtained at the stated threshold. |
| native_margin | Absolute distance of the score from its own binary threshold. |
| native_margin_rank | Frozen ordering used by the archived margin-based selection. |
| native_margin_matched_453 | Membership of the 453-molecule own-margin list. |
| v0_4_accepted, v0_4_risk_call | Acceptance and call from the authors' frozen composite policy, distinct from the comparator's own classification. |
| eligible_primary_comparison, training_overlap_status | Archived comparison eligibility and disclosed overlap-status categories. |

The `handoff/fig1ab.json` six scatter arrays preserve the approved own-score and locally computed MDCK-ER points in Q1/Q3/Q5. `fig1ab_grid.json` stores frozen group metrics; `fig1new_stats.json` and `fig1_permband.json` store historical confidence limits and permutation bands. The current program redraws these archived quantities and validates full-set metrics against the supplied frozen CSV. It does not rerun all historical Figure 3 confidence intervals or regenerate the GNN model outputs.

No Deep-PK individual service scores are included in this component. Fields in the original long comparison table that were specific to that service were excluded from the AZ export. All retained numerical cells are copied without rounding or recomputation.
