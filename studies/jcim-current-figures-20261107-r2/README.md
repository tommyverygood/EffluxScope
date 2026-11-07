# Current JCIM figures and reproducible plotting inputs

Version: `jcim-figures-20261107-r2`. Companion to *ABCB1 efflux prediction at matched coverage and across experimental sources*.

This self-contained archive reproduces the six current main figures and the TOC graphic from frozen inputs, without applying earlier overlays. It also recomputes the Figure 6 balanced-accuracy and retained-fraction percentile intervals from frozen scaffold sufficient statistics. No predictor, embedding, decision threshold or scientific result was refitted or replaced.

This revision enlarges four mathematical subscripts in Figure 2 for legibility at manuscript size. All scientific inputs, numerical results, other figures and AI-use disclosures are unchanged from the [first current-figure release](https://github.com/tommyverygood/EffluxScope/releases/tag/jcim-figures-20261107).

## Reproduce

Use Python 3.12.14, install Ghostscript (`gs`) and Poppler (`pdftoppm`), and make Arial available to Matplotlib. Font files are not redistributed. The reference artwork uses Arial; the verifier stops instead of silently substituting a different font. The recorded validation is on macOS arm64. Other platforms may vary in font/rasterizer rendering even when numerical results agree.

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python verify_reproduce.py --output-dir verification
```

Run from an extracted writable copy. The script verifies all supplied-file checksums, executes the plotting programs and numerical checks, runs 2,000 paired scaffold resamples, and compares all seven output PDFs with the approved references by text, page geometry and pixels at 150 dpi. It exits nonzero on a failed check. `VALIDATION.json` records the release's clean-directory verification; a new run writes its own report. EPS, PNG and SVG variants are produced by the original scripts where supported. PDF byte hashes can vary with creation metadata; text and rendered-pixel equality are checked separately.

## Figure and input map

| Current figure | Component and executable | Main frozen inputs |
|---|---|---|
| 1: task background | `overview/scripts/make_figure_background.py` | approved conceptual background image |
| 2: evaluation overview | `overview/scripts/make_figure_overview.py` | frozen weights, 2,624 score/membership rows, exact source-data excerpts and aggregate evaluations |
| 3: endpoint comparison | `endpoint/run.py` | own-score and author-computed MDCK-ER scatter arrays; frozen comparison metrics, intervals and permutation bands |
| 4: source diagnostics | `results/run_results_figures.py` | source-specific scores, audited Roth summary and regulatory counts |
| 5: structural splits | same results runner | frozen density coordinates, split memberships and performance summaries |
| 6: threshold errors | same results runner | frozen rolling curves and confusion-count summaries |
| TOC | same results runner | frozen matched-coverage results |

Paths in this table are relative to `components/`. `figure_map.json` maps the current numbers, historical filenames and exact output paths. Historical generator names are preserved to avoid confusing them with recomputed scientific analyses.

## What is newly supplied

- Current Figures 1–2 scripts, source asset, frozen model coefficients and complete plotting dependencies.
- Current Figures 3–6 and TOC source versions, shared palette, approved references and all required plotting inputs.
- A 5,248-row table of the authors' locally computed predictions for two released AstraZeneca GNN-MTL endpoints on 2,624 Biogen molecules, with provenance and source/model attribution. This is not a download of an external web-service response.
- Figure 6's 2,312-row scaffold count matrix and resampling program. It reproduces all 100 checked fields, including ten BAcc intervals and the retained-fraction contrast, using seed 20261018 and 2,000 shared draws. The groups retain the frozen ordering; no molecule identifiers, structures or service probability scores are supplied in this matrix. Some groups are singletons; it is analysis input, not a claim of absolute anonymization.

## Scope and limits

Figure 1 re-embeds the approved AI-assisted conceptual illustration; it is not microscopy, an experimental observation or a new image-generation run. Figure 2 checks its traced molecule's frozen linear-model calculation and every displayed numerical quantity. Figures 3–6 read frozen numerical inputs; this release does not rerun model training, every historical bootstrap, calibration fitting or t-SNE.

The Figure 6 interval calculation starts from frozen per-scaffold sufficient statistics. It does not regenerate or independently validate the historical third-party model/service predictions. Historical per-molecule Deep-PK service probabilities and third-party model checkpoints are not redistributed here. Their providers and acquisition code remain documented in the core study release. The new locally computed AZ predictions remove the former missing-scatter-input limitation for current Figure 3. The cluster-count route removes the former missing-input limitation for the stated Figure 6 interval calculation, without claiming full end-to-end regeneration of all comparator outputs.

This figure archive is independent of the manuscript's LaTeX compilation package. It contains no manuscript text, cover letter, reviewer information, credentials or third-party article full texts. Reproduction and file-integrity checks are not additional biological validation.

## Earlier models and analyses

The original frozen v0.1/v0.3/v0.4 models, analysis code, protocols, data dictionaries and offline numerical verifier remain in the unchanged [core study release](https://github.com/tommyverygood/EffluxScope/releases/tag/adc-released-species-20261018). Earlier corrections are historical records; they need not be overlaid to run this package. The current paper evaluates matched coverage and transfer across sources; historical ADC directory/tag names are retained only for provenance.

## Attribution and licences

Author-owned code uses MIT (`LICENSE`); author-owned aggregate evaluation tables and figures are provided under CC BY 4.0 to the extent of the authors' rights. Original source terms and attributions are retained; no blanket relicensing of source-derived material is intended. See `THIRD_PARTY_NOTICES.md`, `licenses/` and component provenance records. Fonts and model checkpoints are obtained separately from their providers.
