# Methods and interpretation

These materials accompany the ADC application analyses conducted with frozen
EffluxScope models on 2026-10-09 to 2026-10-12. This document describes those
methods and the public verification entry point. It does not report a new
training run, new wet-lab experiment, independent source collection, or
prospective validation on 2026-10-18.

## Chemical objects and identity

Curated records distinguish free payloads, actual terminal/downstream release
products, and release intermediates. Exact chemical structures and release
roles are checked against original source locators. A named payload is not
automatically the structure liberated by every linker. The same-structure
controls are distinct from comparisons that change chemical composition.

`source_smiles` and source InChIKeys identify the curated source representation.
The frozen standardization pipeline generates the graph used for model
features and scoring. Both are retained because charge/tautomer handling can
change a standardized graph or key. A source-to-standardized score change is
a representation diagnostic; it does not establish a chemical reaction or
an experimentally changed transport rate. Exact and connectivity-only
development overlap checks are recorded separately.

An unresolved identity is **HOLD**, a curation state that prevents the
corresponding biological interpretation. **Abstention** is a model decision
for a specified structure when support, stability, agreement, or context is
insufficient. The two states must not be pooled as measured negatives.

## Frozen scores and selection

The v0.1 channel is a structure-based model developed from 443 author-curated
ABCB1 label records. The auxiliary v0.3 channel uses a 63-compound direct
matched-parental MDR1 development panel. The full v0.4 strategy combines
frozen score thresholds with bootstrap-model uncertainty, two structural
reference panels, a physicochemical envelope, channel agreement, local
neighborhood evidence, and context restrictions. The three original model
bundles and their manifests define the computational specification.

The primary thresholds are 0.70 for a high-risk candidate and 0.35 for a
low-risk candidate, followed by the support checks. Failure of a required
check produces explicit abstention reasons. No threshold was changed for
this release. The historical field `primary_calibrated_probability` is a
development-model score; its name does not imply external calibration on
ADC species. It is not an experimental ER or a universal assay-independent
substrate probability. An accepted score remains a risk-screening output.

## General benchmarks and equal-count comparisons

The general external transport comparison uses 2,624 identity-audited
Biogen/Fang compounds. Endpoint labels are study-specific ER >= 2 weak
directions. These are not ADC negative controls or gene-isolated universal
substrate truth. Direct transport and weak-label channels can differ in
endpoint alignment.

Selection comparisons hold the accepted count fixed, including k = 453 in
the principal Biogen comparison. Each simple margin policy ranks by absolute
distance from 0.5 on its own score, with deterministic identifiers breaking
ties. The full selector's abstentions are not force-filled. Per-policy class
coverage and confusion counts accompany accepted-subset metrics. The data
retain cases in which a simple margin rule outperforms the full selector;
the full rule has no universal superiority claim.

The ASAP sensitivity analysis aggregates 62 measurements to 56 unique
compounds using the recorded median-log2ER rule and preserves the observation
ledger. It is retrospective, using previously seen outcomes. Five accepted
compounds are all measured positives under that endpoint; specificity,
balanced accuracy, or negative predictive performance cannot be inferred
from that selected subset. A paired confidence interval crossing zero does
not demonstrate equivalence or a ranking advantage.

## Representation and transport-evidence audits

Normalization variants and atom-order/equivalent-SMILES checks distinguish
input representation effects from parent--species differences. The original
13-structure check contains 221 equivalent-SMILES rows; these are repeated
representations, not 221 independent molecules. Enumeration limits for some
species are recorded and must remain visible. A nonchiral fingerprint does
not establish discrimination between stereoisomers.

Transport records preserve cell system, transporter arm, concentration,
duration, inhibitor and its concentration, Papp/ER units, uncertainty type,
source precision, and eligibility. Assay conditions, transporter arms, and
baseline/inhibitor rows are not exchanged to fabricate matched controls.
Functional resistance or cytotoxicity is not relabeled as direct transport.
Conflicting source calls and incomplete controls remain explicit. The five
changed computational pairs do not close the missing qualified same-system
parent--release truth set. Dual-component scores do not define a combined ER.

## Fixed chemical-region support audit

The final
[`protocol.json`](outputs/20261012_adc_scope_revision/protocol.json)
defines a retrospective secondary diagnostic whose source outcomes were
already known. The region rules were fixed before the new regional
calculation. Regions are calculated on the frozen standardized graph:

| Region | Operational rule |
|---|---|
| High molecular weight | Molecular weight >= 600 Da |
| High polarity | TPSA >= 140 square angstroms |
| Multiple C(=O)-N motifs | At least 4 unique SMARTS `[CX3](=[OX1])[NX3]` matches and molecular weight >= 500 Da |

The motif count is an operational carbonyl-nitrogen pattern count, not a
validated peptide-bond count. These are chemical stress regions, not a
validated ADC applicability domain. Counts are separate within each source;
concentrations and repeats do not create additional molecules. Development
identity overlaps are audited under the recorded source/model exact and
connectivity rules.

The population reporting gate requires at least 30 compounds, 5 in each
measured class, and 10 Murcko scaffolds. The matched-count gate requires
at least 20 selected compounds, 5 in each selected measured class, and 10
selected scaffolds. These are reporting safeguards, not power guarantees.
The final protocol supersedes the earlier feasibility proposal.

Biogen has 6, 18, and 1 compounds in the three regions; ASAP has none in
all three. All six source-region combinations fail the gates. Performance
comparisons and paired intervals are therefore `null` with reasons, and no
regional bootstrap was executed. A conditional 2,000-draw scaffold-bootstrap
plan in the protocol does not mean those draws were run. Thresholds were not
relaxed. The 14 ADC structures supply region coverage descriptions only.

## Public verification

The supported verifier reads frozen inputs and writes to a new output
directory. It checks recorded matched-count point estimates, 2,694 feature
records, six reporting statuses, and fourteen model predictions. It does not
recompute historical bootstrap intervals or reconstruct excluded original
documents. Verification status is established by its generated receipt.

