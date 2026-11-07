# Current manuscript Figures 1 and 2

This component rebuilds the approved background illustration (Figure 1) and the method/results overview (Figure 2). It is self-contained: it does not search a research checkout, refit a model, access a remote service, or require the manuscript source.

Run with the release's pinned Python 3.12 environment:

```sh
python scripts/make_figure_background.py
python scripts/make_figure_overview.py
python scripts/verify_overview.py
```

The scripts can also be invoked by absolute or relative path from another working directory. They write only to this component's `outputs/` folder.

## What is reproduced

- Figure 1 re-embeds the previously approved AI-generated conceptual raster and reapplies its editable panel headings. It does not rerun image generation, and the illustration is not experimental microscopy or measured molecular evidence.
- Figure 2 uses the same drawing geometry, colours, font sizes and numerical inputs as the current paper. It recalculates the traced compound's 2,048 Morgan bits and eight descriptors, frozen logistic weighted sum, sigmoid and original calibration. The resulting primary score must agree with the stored value within 1e-6. It also checks 2,624 molecules, three retained lists of 453, 56-compound ASAP aggregation, 11 Langthaler compounds, 12 identity-admitted Roth compounds, 33 regulatory cases, and the reported correlation/accuracy values.
- Bootstrap intervals shown in the overview are frozen recorded results read from the supplied tables. This drawing script does not recompute those intervals or refit the model.

## Fonts and rendering

The maintained artwork uses Arial. Install Arial through your platform's normal licensed font installation if exact typography is required. No font files are redistributed. Figure 1 uses the standard macOS Arial Bold location when present and otherwise the ReportLab Helvetica-Bold fallback; Figure 2 uses Matplotlib's font lookup. Missing/different fonts or different rasterizers can change pixels without changing numerical values. The release validation records the exact tested Python packages and renderer.

The source code uses Pillow, ReportLab, NumPy, Matplotlib, SciPy, RDKit and joblib. The copied chemistry helper is the original project helper; it is included so feature calculations do not depend on another checkout.

## Inputs and provenance

`data/EXPORT_PROVENANCE.json` identifies every original project-relative input, its SHA-256, the exported path/hash, and whether the file was copied unchanged or reduced to exactly the rows/columns read by this figure. Row order and all selected numeric strings are preserved. The only change to frozen model parameter metadata removes a personal absolute machine path; numerical values are untouched.

`data/source_inputs/` contains:

- one Biogen molecule's identifier and standardized structure for the deterministic traced example;
- 62 ASAP rows, grouped by structure identifier to 56 compounds with median log2(ER);
- the 11 eligible Langthaler rows;
- the 12 identity-admitted Roth rows for the historical Morgan-logistic baseline, with GI50 resistance ratios (not transport ER);
- the frozen regulatory 2-by-2 counts and balanced accuracy.

These inputs contain the authors' predictions and published source facts, not third-party external-service predictions, article full texts, credentials, or model checkpoints. The Biogen/Fang source is https://github.com/molecularinformatics/Computational-ADME (MIT). ASAP is ChEMBL assay CHEMBL5674853 / document CHEMBL5674852; maintain its source-data attribution and licence. Langthaler and Roth source locators are retained in the exported rows. The release-level `THIRD_PARTY_NOTICES.md` governs per-source attribution; the code's MIT licence does not override source-data terms.

`data/expected_numeric_checks.json` records the approved figure's numerical/provenance values (with location-only fields removed). `scripts/verify_overview.py` checks those values against the freshly generated provenance. Rendering comparison is recorded separately; passing numerical checks alone is not a claim of pixel identity.
