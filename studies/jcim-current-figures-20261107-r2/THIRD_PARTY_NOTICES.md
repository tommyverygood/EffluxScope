# Source-specific attribution and scope

This update retains the original study's source-specific attribution in `licenses/CORE_STUDY_NOTICES.md` and component provenance files. The MIT code licence does not replace source-data terms.

## Biogen measured endpoints

Fang et al., *Prospective validation of machine learning algorithms for absorption, distribution, metabolism, and excretion prediction: an industrial perspective*, J. Chem. Inf. Model. 2023, 63, 3263–3274, DOI: 10.1021/acs.jcim.3c00160. Public repository: https://github.com/molecularinformatics/Computational-ADME. Its MIT notice is preserved in `licenses/Biogen_MIT_LICENSE.txt`. Structures/outcomes are standardized, filtered and joined to author-generated model scores; these transformations are documented in the frozen study.

## Locally generated GNN-MTL predictions

Ivers Ohlsson et al., *Prediction of permeability and efflux using multi-task learning*, https://doi.org/10.5281/zenodo.16948542 (model v1; Chemprop 2.1.0). The official record labels the model Apache-2.0 and open access. `licenses/Apache-2.0.txt` preserves the model licence for attribution; the checkpoint itself is not included. Model outputs 2 and 3 are MDCK ER and NIH-MDCK ER. The exported table preserves the authors' local scores, measured endpoint, IDs, native margins and fixed selection mask; it excludes unrelated external-service fields. Predicted values are not new experimental measurements. The predictor remains distinct from the authors' weak-label/assay-aware models.

## Reference-derived and source-specific material

ChEMBL structures and derivatives retain CC BY-SA 3.0; Daood et al. curated labels and the Sóskuti direct panel retain their source CC BY 4.0 terms. Frozen model coefficients, trace/example records and source-summary inputs retain provenance rather than being represented as new transport measurements. Exact selected fields and original file hashes are listed in the overview export record. See the retained core-study notices for full citations and modification notices.

## Comparator aggregates and Figure 6 sufficient statistics

Deep-PK and ADMET-AI remain attributed to their publications/providers as in the manuscript and core study. ADMET-AI source licence is preserved in `licenses/ADMET_AI_MIT_LICENSE.txt`. This archive includes the authors' evaluation metrics, curve summaries and fixed-scaffold confusion counts. It does not redistribute Deep-PK per-molecule probabilities, service job identifiers or a third-party service database. The count matrix is an author-derived numerical representation for repeating the stated bootstrap; its components are explained in `components/uncertainty/README.md`. Singleton groups exist, and no privacy/non-inference guarantee is asserted.

## Conceptual imagery and fonts

The approved background asset is an AI-assisted conceptual illustration supplied as an illustration, not scientific raw data. No third-party paper figures or microscopy images were copied into this release. Arial and other system font files are not supplied and must be acquired from their respective providers.
