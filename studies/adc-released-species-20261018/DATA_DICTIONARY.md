# Study data dictionary

Paths are relative to this study directory. CSV files are UTF-8; some retain
a UTF-8 byte-order mark, so readers should use `encoding="utf-8-sig"`.
JSON fields embedded in CSV are serialized strings, not extra table rows.
The earlier [`数据/DATA_DICTIONARY.md`](数据/DATA_DICTIONARY.md) describes legacy
inputs. This document covers the ADC follow-up and chemical-support tables.

## Shared conventions

| Field or concept | Meaning |
|---|---|
| `entity_id`, `uid`, `record_id` | Stable within the named table/source; an observation ID is not necessarily a distinct molecule |
| `source_id`, `source_doi`, `source_url`, `primary_attachment_url`, `source_location` | Citation and exact table/figure/page or attachment locator; URLs are recovery pointers, not a license grant |
| Source vs standardized SMILES/InChIKey | Curated source graph versus frozen model-input graph; preserve both and their transformation status |
| Connectivity block | First InChIKey block used for connectivity matching; this deliberately loses full stereochemical/protonation detail |
| Booleans | Depending on historical table: `yes`/`no`, `True`/`False`, or `1`/`0`; interpret using that table's schema |
| Missing values | Blank, JSON `null`, or a stated status denotes missing/not applicable/not reportable; never automatically replace with zero |
| `HOLD` | Chemical identity/release-role curation is unresolved; not a model class or measured negative |
| `abstain` | Model withholds a class for a specified structure; not a biological third class |
| `*_sha256`, manifest hashes | Byte-level provenance. Historical manifests may refer to pre-export files; use `EXPORT_PROVENANCE.json` for original/export mapping |

## Current inventory and model outputs

Files:

- `outputs/20261012_adc_completion/scoring/species_predictions.csv`
- `outputs/20261012_adc_completion/scored_family_inventory.csv`
- `outputs/20261012_adc_completion/family_expansion/new_machine_entities.csv`

The current inventory has 14 scored structures, not 14 independently measured
transport outcomes. Family counts are 4/5/3/1/1 for auristatin/maytansinoid/
TOP1/PBD/calicheamicin.

| Fields | Meaning / unit |
|---|---|
| `entity_id`, `entity_name`, `family`, `role` | Chemical identity, payload family, and free/terminal/downstream/intermediate role |
| `source_smiles`, `source_inchikey`, `source_formula`, `input_source` | Curated source identity and locator |
| `input_smiles`, `standardization_status`, `standardized_isomeric_smiles`, `standardized_inchikey` | Actual model input and resulting frozen standardized representation |
| `chemical_object_type`, `requested_context`, `context_compatibility` | Declared small-molecule role, requested use, and scope check |
| `primary_calibrated_probability`, `direct_panel_probability` | Frozen primary and auxiliary unitless scores in [0,1]; the historical probability field name is not external ADC calibration |
| `primary_low_threshold`, `primary_high_threshold` | Fixed primary selection cutoffs, 0.35 and 0.70 |
| `bootstrap_probability_p10`, `bootstrap_probability_median`, `bootstrap_probability_p90`, `bootstrap_interval_width` | Frozen bootstrap-model score distribution summaries; not a clinical prediction interval or external ADC calibration interval |
| `risk_call`, `decision_zh`, `accepted_prediction` | High/low/abstain decision, Chinese label, and acceptance flag; 5 high and 9 abstain in this inventory |
| `primary_max_tanimoto`, `direct_max_tanimoto` | Maximum fingerprint similarity to each reference panel, unitless [0,1] |
| `*_local_positive_fraction`, `*_local_conflict_index` | Local reference-label summaries used by frozen selection checks |
| `physchem_max_abs_robust_z`, `physchem_threshold_q99` | Largest absolute robust descriptor deviation and frozen envelope threshold |
| `applicability_status`, `uncertainty_status`, `abstention_reasons_json`, `abstention_reasons_zh_json` | Machine-readable support/stability status and ordered reason lists |
| `nearest_primary_references_json`, `nearest_direct_references_json`, `physicochemical_descriptors_json` | Embedded reference evidence and descriptor objects; source-specific licenses apply to reference records |
| `exact_reference_identity`, `primary_exact_identity`, `direct_exact_identity`, `*_connectivity_overlap` | Development/reference identity checks; an overlap is not validation |
| `stereo_blind_collision_warning` | Warning that the model representation may collapse distinct stereochemical identities |
| `tautomer_status`, `enumeration_complete`, `diagnostic_enumeration_hold` | Diagnostic enumeration status; limits do not vanish because a numeric score exists |
| `model_version`, `validation_status`, `claim_boundary`, `identity_is_not_biological_validation`, `identity_limitation` | Version, evidence limits, and interpretation restrictions |
| `chemically_changed_pair_family` | Family-level membership in changed comparisons; not the number of measured pairs |

## Identity, pairs, and dual-component comparisons

Key files under `outputs/20261009_adc_pair_experiments/`:

- `identity/approved_structures_v2.csv` and `identity/identity_sensitivity_inputs.csv`
- `pairs/results_curated_v2/entities_scored_or_hold.csv`
- `pairs/results_curated_v2/parent_release_score_contrasts.csv`
- `pairs/results_curated_v2/pair_release_stage_classification.csv`
- `pairs/results_curated_v2/source_vs_frozen_model_structure_audit.csv`
- `pairs/results_curated_v2/dual_component_scores.csv`

Entity and pair identifiers link back to curated structures. Parent/species
columns specify the actual compared chemical objects. Score differences are
unitless differences on frozen model scales; they are not changes in measured
ER. Release-stage/eligibility columns distinguish a terminal/downstream product,
an intermediate, a same-chemical control, and a HOLD. Five changed comparisons
span two families; same-identity controls must not be counted as additional
changed families. Dual rows contain individual component predictions, not
combined ER, synergy, or an optimized mixture ratio.

## Transport evidence and source recovery

Files under `outputs/20261012_adc_completion/`:

- `public_search/source_inventory.csv` and `public_search/new_source_observations.csv`
- `public_search/parent_release_pair_eligibility.csv`
- `family_expansion/family_candidate_evidence.csv`
- `family_expansion/direct_transport_observations.csv`
- `family_expansion/Weng2023_FigS1E_transport.csv`
- `design_audit/adc_component_candidate_eligibility.csv` and related anchor/eligibility tables

| Fields | Meaning / unit |
|---|---|
| `compound`, `entity_id`, `entity_type`, `object_role` | Observed identity and role; a publication name alone may not resolve an exact released species |
| `transporter_gene`, `transporter_arm`, `cell_system` | Gene/context and assay arm; separate P-gp, BCRP, matched-parental, and inhibitor conditions |
| `test_concentration_uM`, `inhibitor_concentration_uM` | Micromolar; a missing inhibitor concentration remains missing |
| `incubation_hours` | Hours |
| `papp_ab_mean_1e6_cm_s`, `papp_ba_mean_1e6_cm_s` | Numeric Papp in units of 10^-6 cm/s; associated `*_sd` uses the same scale |
| `efflux_ratio` | Dimensionless ER as reported/defined by that source, normally Papp BA / Papp AB; do not assume all sources use identical normalization |
| `baseline_to_condition_er_reduction_fold` | Dimensionless paired baseline-to-condition ratio where the source supports the pairing; not the model score delta |
| `uncertainty_type`, `numeric_precision` | Source SD/SEM/interval or extraction precision; these are not interchangeable |
| `source_assigned_call`, `transport_evidence_tier`, `truth_set_use`, `permitted_analysis` | Original qualitative call, evidence type, and allowable use in this study |
| `identity_review_status`, `identity_limitation`, `source_review_status` | Identity/evidence review status and unresolved limitations |
| `parent_release_same_chemical`, `changed_pair_eligible`, `parent_release_pair_role` | Whether the source describes a chemical change and a qualifying comparison |
| `original_curation_seen`, `prior_curation_seen`, `development_identity_block_seen_v01/v03` | Prior exposure or development overlap; these prevent a retrospective check being labeled blind |
| `qualifying_new_direct_ABCB1_negative` | Eligibility for the strict negative pool, not merely ER < 2 in an unrelated general dataset |

Repeated concentrations, inhibitor arms, or measurements can create several
rows for one entity. Source rows can also be repeated across extraction,
eligibility, and reporting tables. Do not sum them into an experiment count.
Blank matched baselines cannot be borrowed from another arm. Functional
cytotoxicity/resistance measurements remain distinct from direct transport.

## Benchmark and compound aggregation tables

`outputs/20261009_adc_pair_experiments/benchmark/` contains `metrics.csv`,
`selection_ledger.csv`, `paired_contrasts.csv`, `unattainable_targets.csv`,
and `protocol.json`. `version`, `policy`, and `target` identify a comparison;
`selected_n`/class counts/confusion counts give its denominator. Coverage,
accuracy, BAcc, sensitivity, specificity, PPV, NPV, and MCC are unitless.
Margin policies rank on their own channel; identical counts do not imply
identical selected molecules. Historical bootstrap contrasts retain their
recorded interval and resampling unit; the public verifier checks point
estimates without rerunning those intervals.

`outputs/20261009_adc_followup_validation/asap_secondary/` and its
`repetition_audit/` record the 62-measurement to 56-compound aggregation.
`compound_aggregation_ledger.csv` and `unique_compound_identity_and_aggregation.csv`
preserve identity and observation membership. Aggregate ER is derived through
the recorded median-log2ER rule, not by counting repeats as new compounds.
`paired_rho_contrasts.csv` is a retrospective rank-correlation sensitivity
analysis. Missing class-dependent metrics are not zero performance.

## Representation audit tables

`outputs/20261009_adc_followup_validation/normalization/` contains variant
scores, stage/source-structure audits, and same-identity controls. Its
`atom_order_invariance/` tables retain the equivalent-SMILES inputs, identity
checks, score ranges, and enumeration status. A variant/representation row is
not an independent chemical observation. The original 221 rows represent 13
structures. Equal scores for equivalent SMILES test computational stability,
not stereo-specific prediction or biological correctness.

## Chemical support features and reporting gates

Files under `outputs/20261012_adc_scope_revision/`:

- `feasibility/features.csv`: 2,694 records = Biogen 2,624 + ASAP 56 + ADC 14.
- `feasibility/count_feasibility.csv` / `.json`: within-source denominators and selections.
- `regional_evaluation_status.csv` / `.json`: six source-region statuses.
- `protocol.json`: final thresholds and gates; `feasibility/protocol.json` is an earlier, superseded proposal.

| Fields | Meaning / unit |
|---|---|
| `source`, `uid`, `entity_name` | Source and within-source entity key; 2,694 rows are not a globally deduplicated molecule count |
| `source_smiles`, `source_exact_inchikey_computed`, `standardized_isomeric_smiles`, `standardized_inchikey`, `connectivity_block` | Source and frozen input identity, with stereo-aware and connectivity matching kept separate |
| `molecular_weight` | RDKit molecular weight, Da |
| `tpsa` | RDKit topological polar surface area, square angstroms |
| `amide_count` | Unique matches of `[CX3](=[OX1])[NX3]`; operational C(=O)-N motif count, not a validated peptide-bond count |
| `murcko_scaffold` | Frozen scaffold representation used for distinct-scaffold counts |
| `measured_er`, `er_ge2_weak_label`, `observation_count`, `posthoc_aggregation` | Source endpoint/weak direction and aggregation metadata; ADC rows do not acquire labels from neighboring general-data compounds |
| `full_accepted`, `primary_score`, `auxiliary_score`, `risk_call`, `abstention_reasons_json` | Frozen selection and score outputs |
| `MW_ge600`, `TPSA_ge140`, `amide_ge4_MW_ge500` | Inclusive fixed region memberships on the standardized graph |
| `development_*_overlap`, `validation_role` | Development-overlap checks and allowed analytic role |
| `population_n`, `positive_n`, `negative_n`, `scaffold_n` | Region population, measured class counts, and distinct scaffolds within one source |
| `full_selected_n`, `full_selected_positive_n`, `full_selected_negative_n`, `full_selected_scaffold_n`, `full_coverage` | Native full-selector accepted denominator, classes/scaffolds, and fraction |
| `adc14_region_n` | ADC coverage description for that region; not an ADC accuracy denominator |
| `failed_report_gates`, `regional_comparison_status`, `null_reason` | Exact failure reasons and `not_evaluable` status under the final safeguards |
| `BAcc`, `MCC`, `paired_delta_BAcc`, `paired_CI95`, `bootstrap_executed` | All six statuses retain null comparison metrics and `false` bootstrap execution; no zeros are imputed |

Population gates require n >= 30, each class >= 5, and scaffolds >= 10;
selection gates require k >= 20, each selected class >= 5, and selected
scaffolds >= 10. A null reflects the reporting safeguards, not a claim that
the quantity is mathematically undefined for every possible small sample.
Sources and regions can overlap; region counts cannot be added to infer
globally independent compounds.

