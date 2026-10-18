# EffluxScope ADC payload and released-species reproducibility materials

Release identifier: `adc-released-species-20261018`. Prepared on 2026-10-18.
This study package contains the data and methods supporting the ADC application assessment. The offline checks recorded in [VALIDATION.json](VALIDATION.json) passed in a newly installed Python 3.12.14 environment. No archival DOI is assigned to this GitHub study release.

Versioned repository location:
<https://github.com/tommyverygood/EffluxScope/tree/adc-released-species-20261018/studies/adc-released-species-20261018>.

## Research scope

EffluxScope evaluates structure-based ABCB1/P-glycoprotein efflux-risk scores
for chemically specified ADC payloads, released small molecules, and selected
release intermediates. It combines frozen scores, chemical-support checks,
uncertainty, and explicit abstention reasons. Chemical identity is curated
from sources before a molecule is scored. The tool does not infer every
released species from an intact ADC name or construct.

The application inventory contains **14 structures in five payload families**:
auristatin (4), maytansinoid (5), TOP1/camptothecin (3), PBD (1), and
calicheamicin (1). The full frozen strategy produces **5 accepted high-risk
calls, 9 abstentions, and no accepted low-risk calls**. Five chemically
changed parent--species comparisons span **two families** and include
terminal/downstream products and intermediates. These are computational
comparisons, not five measured transport pairs. Dual-payload examples are
component-wise analyses only.

This is a development-stage application assessment. The available source
review did not yield a qualified same-system, chemically changed
parent--release transport truth set or a qualified ADC negative class.
Consequently the materials do not establish ADC-specific accuracy,
superiority of released-species inputs over parent proxies, joint dual-payload
transport, synergy, or intact-ADC efficacy. General transport benchmarks
evaluate the underlying scores and selection rules; their performance does
not establish performance on the ADC inventory. Null and unfavorable findings
are retained.

## File navigation

The original Chinese relative directory names are preserved. Run commands
from this study directory, without flattening its directory tree.

| Entry | Contents |
|---|---|
| [METHODS.md](METHODS.md) | Frozen workflow, analysis chronology, chemical identity, and reporting safeguards |
| [DATA_DICTIONARY.md](DATA_DICTIONARY.md) | New ADC, transport, pair, benchmark, representation, and chemical-region fields |
| [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) | Source-specific rights, attribution, modifications, and exclusions |
| [CITATION.cff](CITATION.cff) | Citation metadata for these reproducibility materials |
| [requirements-public.txt](requirements-public.txt) | Pinned dependencies for the public verification entry point |
| [scripts/verify_public_bundle.py](scripts/verify_public_bundle.py) | Verification against frozen scores, features, counts, and point estimates |
| [EXPORT_PROVENANCE.json](EXPORT_PROVENANCE.json) | Exported paths, original/export hashes, and provenance-only transformations |
| [模型/](模型/) | Frozen v0.1, v0.3, and v0.4 model bundles and original inference code |
| [数据/](数据/) | Selected processed inputs, historical scripts, and the earlier data dictionary |
| [outputs/20261009_adc_pair_experiments/](outputs/20261009_adc_pair_experiments/) | Curated identities, parent--species and dual-component scores, matched-count comparisons |
| [outputs/20261009_adc_followup_validation/](outputs/20261009_adc_followup_validation/) | Representation checks, transport-evidence review, and compound-level ASAP sensitivity analysis |
| [outputs/20261012_adc_completion/](outputs/20261012_adc_completion/) | Fourteen-species inventory, family expansion, source locators, and eligibility audit |
| [outputs/20261012_adc_scope_revision/](outputs/20261012_adc_scope_revision/) | Fixed chemical-region features, counts, final protocol, and reportability status |
| [handoff/](handoff/) | Historical figure input arrays; these are not a new experimental dataset |

Start with
[species_predictions.csv](outputs/20261012_adc_completion/scoring/species_predictions.csv),
[scored_family_inventory.csv](outputs/20261012_adc_completion/scored_family_inventory.csv),
[family_candidate_evidence.csv](outputs/20261012_adc_completion/family_expansion/family_candidate_evidence.csv),
and the final
[regional_evaluation_status.json](outputs/20261012_adc_scope_revision/regional_evaluation_status.json).
The final scope protocol supersedes the earlier feasibility-stage proposal.

## Verification in a fresh environment

Use **Python 3.12** and the pinned public requirements. Model bundles depend
on the recorded scikit-learn/RDKit runtime. A virtual environment copied from
another computer is not a substitute for installing these dependencies.

```sh
python3.12 -m venv .venv-public
. .venv-public/bin/activate
python -m pip install -r requirements-public-lock.txt
python scripts/verify_public_bundle.py --root . --output-dir /absolute/fresh-check
```

Replace `/absolute/fresh-check` with an absolute path to a **new, nonexistent
directory**, preferably outside the study tree. Every generated verification
file is written there. Read the resulting `verification.json`: a successful
receipt, rather than this instruction, establishes which checks completed.

The verification entry point is designed to recompute fourteen frozen ADC
scores, molecular descriptors and fixed-region memberships, reporting-gate
failures, and matched-count point estimates. It performs no model fitting,
recalibration, network requests, new transport experiments, or bootstrap
reruns. Historical confidence intervals are supplied as recorded results.
`--skip-inference` explicitly omits the fourteen model reruns and gives a
different receipt status; it is not equivalent to full verification.

## Source recovery and historical reconstruction

Third-party article/assessment full texts, supplemental documents, external
model weights, private reviewer links, and Deep-PK per-molecule outputs are
excluded. Retained observations give source attribution and page/table/figure
locators. Recover original documents from
[source_inventory.csv](outputs/20261012_adc_completion/public_search/source_inventory.csv)
and the source URL/DOI/attachment fields in
[direct_transport_observations.csv](outputs/20261012_adc_completion/family_expansion/direct_transport_observations.csv)
and
[family_candidate_evidence.csv](outputs/20261012_adc_completion/family_expansion/family_candidate_evidence.csv),
subject to the provider's terms. The manifest identifies historical inputs;
it does not grant permission to redistribute a source document.

Historical collection and reconstruction scripts are retained for methods
transparency. Some require excluded original documents, upstream downloads,
optional dependencies, or original historical paths. **The entire historical
pipeline is not promised to run offline from this package.** The supported
offline entry point is the verification command above. Rebuilding historical
analyses should be done in a separate working copy after recovering their
documented inputs.

Historical run manifests describe the original runs. Where host paths or
instruction-only metadata were normalized for this export,
`EXPORT_PROVENANCE.json` records original and exported hashes separately;
scientific values and frozen models are retained. Neither successful file
verification nor redistribution closes the missing ADC biological validation.


The historical `handoff/fig1ab.json` export retains the internal-model scatter and aggregate grid only. Three external-model individual-prediction scatter arrays are omitted; this transformation is recorded in `EXPORT_PROVENANCE.json`.

## Recorded computational verification

The extracted archive reproduced all 14 species (eight numerical and seven categorical fields each), the five changed comparisons and five same-identity controls, the ASAP 62-observation/56-compound aggregation, eight matched-count point-estimate rows, 2694 descriptor rows, eight source/region counts and six non-evaluable regional rows with 24 explicit nulls. No model was fitted and no bootstrap confidence interval was rerun. The 490 scientific export files and the verifier are byte-identical to the checked extraction.
