# Third-party notices for the ADC study materials

This notice supplements the repository-level
[`THIRD_PARTY_NOTICES.md`](REPOSITORY_THIRD_PARTY_NOTICES.md) and
[`LICENSE`](LICENSE). The MIT license applies to author-owned software
and documentation. It does not replace a third-party data license, convert
source records to MIT, or grant rights to excluded documents or models.
Attribution follows the underlying record even when it appears in a derived
table, fingerprint, nearest-neighbor record, or frozen model bundle.

## 1. ChEMBL-derived reference structures: CC BY-SA 3.0

The 443-record primary reference pool includes ChEMBL-derived structures and
fingerprints. The corresponding training audit and embedded reference
records in `ABCB1_screener_v0_1.joblib` and
`ABCB1_selective_uncertainty_v0_4_dev.joblib` retain the obligations described
in the repository notice.

- Source: ChEMBL, EMBL-EBI, <https://www.ebi.ac.uk/chembl/>.
- License: CC BY-SA 3.0 Unported,
  <https://creativecommons.org/licenses/by-sa/3.0/>.
- Source-publication access information: April 2025; a numbered ChEMBL
  release was not documented and is not invented here.
- Modifications include structure standardization, fingerprint generation,
  reference selection, and structural-support calculations.

Redistributed adaptations of those structures/fingerprints retain CC BY-SA
3.0 and ChEMBL attribution. Their presence inside a mixed model bundle does
not remove these obligations. The authors' MIT software license does not
extend to those records.

## 2. Daood et al. curated labels: CC BY 4.0

The primary reference labels derive from Daood et al., *Machine Learning
Modeling for ABC Transporter Efflux and Inhibition: Data Curation, Model
Development, and New Compound Interaction Predictions*, *Molecular
Pharmaceutics* (2025), DOI
<https://doi.org/10.1021/acs.molpharmaceut.5c01065>, PMCID PMC12587445.
The recorded license is CC BY 4.0,
<https://creativecommons.org/licenses/by/4.0/>.

Labels are author-curated literature evidence, not new primary measurements
by this project. Filtering, identity resolution, standardization, aggregation,
and training-interface construction are modifications. Source structures
that originated in ChEMBL retain the separate structure license above.

## 3. Sóskuti et al. direct panel: CC BY 4.0

The direct model and full selector embed records/labels derived from the
63-compound retained subset of Sóskuti et al., *Applicability of MDR1
Overexpressing Abcb1KO-MDCKII Cell Lines for Investigating In Vitro Species
Differences and Brain Penetration Prediction*, *Pharmaceutics* 2024;16:736,
DOI <https://doi.org/10.3390/pharmaceutics16060736>, PMCID PMC11207571.

- Source: supplementary Table S1, original measurements at 1 micromolar,
  120 minutes, <https://www.mdpi.com/article/10.3390/pharmaceutics16060736/s1>.
- License: CC BY 4.0, <https://creativecommons.org/licenses/by/4.0/>.
- Copyright: 2024 by the authors; licensee MDPI, Basel, Switzerland.
- Attribution basis: article-level CC BY notice and MDPI's supplementary
  material policy, as detailed in the repository notice.
- Modifications: parsed values, identity correction and standardization,
  source-QC and interval exclusions, hMDR1/Mock net-ratio calculation,
  log2 transformation, and study-specific thresholding. The retained panel
  contains 63 entities; the source table is not redistributed as a PDF.

The derived screening-evidence label is not compound-level inhibitor-rescue
confirmation. The original license's warranty/disclaimer provisions apply.

## 4. General external data and study-derived calculations

Biogen/Fang structures and measured MDR1-MDCK endpoints originate from the
public Computational-ADME dataset, distributed under MIT:
<https://github.com/molecularinformatics/Computational-ADME>.
Its attribution and license remain applicable to copied source fields.
Scores, feature calculations, matched-count selections, and statistical
summaries are derived analyses by this project. Neither model performance
nor a derived label changes the provenance of a source measurement.

ASAP and other literature-derived records retain their publication DOI/URL,
source identifier, exact locator, and the modifications recorded in their
curation tables. This package does not assign all mixed-source records a
blanket MIT or Creative Commons license. Readers must preserve the applicable
source attribution and any source-specific license when redistributing
those records. Independently calculated model scores do not license an
underlying third-party structure or database collection.

## 5. PubChem identity facts and source-specific ADC evidence

Selected PubChem property-response JSON files contain factual identity
crosschecks for the ADC case study. Cite PubChem, the indicated CID, and
the recorded retrieval/source locator. For example, CID 86567977 is used in
the calicheamicin DMH identity crosscheck. PubChem aggregates multiple
contributors; a factual identity check is not a claim of ownership over all
upstream contributor records or a license grant for an entire source dump.

Other ADC records contain limited transport/identity facts extracted by the
authors from articles, regulatory assessments, database records, and patents.
Sources, page/table/figure locators, concentrations, assay arms, extraction
precision, and limitations remain in the observation and source-inventory
tables. These are modified, source-attributed research records, not original
transport experiments performed by this project. Short factual fields and
URLs do not grant permission to republish the original document. Original
third-party prose, figures, full articles, and supplemental documents are
not included in this study export.

## 6. External predictors and excluded content

No third-party model weights or third-party per-molecule predictor outputs
are included as part of this ADC export. Aggregate comparison metrics may
name upstream methods for attribution. Obtain external models from their
original distributions under their applicable terms:

- AstraZeneca GNN-MTL: upstream checkpoint described as Apache-2.0 in the
  existing project record; checkpoints are excluded here.
- ADMET-AI: upstream software described as MIT; its code/weights and
  per-molecule outputs are excluded here.
- Deep-PK: the historical record did not identify a license/terms statement
  at access. Its downloaded training/validation/test files, archives, and
  per-molecule outputs are excluded from this extension; no permission is
  inferred from availability on a web server.

Also excluded are third-party PDF/HTML/DOCX full texts and full-text
extractions, private reviewer access URLs, credentials, personal/internal
materials, and original manuscript backups. Retained source hashes and
historical paths document provenance only; they do not mean the excluded
file is supplied or licensed for redistribution.

## 7. Export and reuse

`EXPORT_PROVENANCE.json` records the source and exported byte hashes and any
provenance-only normalization. Scientific thresholds, arrays, identities,
labels, and frozen model bundles are retained under their recorded source
terms. Historical manifests refer to the original analysis files and can
therefore differ from an export with normalized host paths.

When reusing these materials, cite this release metadata, the EffluxScope
repository, and the specific upstream publications/datasets actually used.
The release citation does not substitute for source attribution. No
published-article status or archival DOI is asserted by this notice.

