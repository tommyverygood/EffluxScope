# EffluxScope figure-revision addendum, 18 October 2026

This addendum supplies the maintained plotting sources and corrected Figure 4
plotting data for the ADC-focused study. It overlays a **new working copy** of
the original public study archive; it does not modify that frozen release,
retrain the models, or replace the underlying frozen predictions.

Base release: [adc-released-species-20261018](https://github.com/tommyverygood/EffluxScope/releases/tag/adc-released-species-20261018).
Required asset: `EffluxScope_ADC_reproducibility_20261018.zip`.
SHA256: `19fd294d893345af41edb08a1a8b86b9b9a4bfe7a859e50d52c44998be26938f`.

## Apply and redraw

Unzip the base archive into a fresh directory. Its study root is named
`adc-released-species-20261018`. Copy the **contents** of this addendum's
`overlay/` into that study root, preserving relative paths and replacing the
listed files. Keep the original archive unchanged.

Use Python with NumPy, Pandas, Matplotlib, Pillow and RDKit. The recorded
environment is Python 3.12.14, NumPy 2.5.1, Pandas 2.2.3, Matplotlib 3.11.2,
Pillow 11.3.0 and RDKit 2026.03.4. Matching versions are listed in
`requirements-plotting.txt`. Ghostscript (`gs`) must be available for the EPS
exports of Figures 3 and 4. Fonts are **not** redistributed: Figure 5 uses the
original Songti SC where available and the palette module's documented font
fallback otherwise. Font substitution can change text metrics.

From the overlaid study root:

```sh
python -c "from pathlib import Path; [Path(p).mkdir(parents=True, exist_ok=True) for p in ['论文主稿/latex_jcim/figures', '论文主稿/latex_jcim/tables', '论文主稿/figures_src/fig5_structures']]"
python 论文主稿/latex_jcim/tools/make_figure3_jcim.py
python 论文主稿/latex_jcim/tools/make_figure4_jcim.py
python 论文主稿/figures_src/make_figure5.py
```

Figures are written under `论文主稿/latex_jcim/figures/`. Figure 5 additionally
writes its ADC data table and verification outputs. The commands consume the
base archive's existing data: `handoff/fig3_data.json`, public metric summaries,
regulatory/ASAP summary files, and
`outputs/20261012_adc_completion/scoring/species_predictions.csv`.

## What changed

- Figure 3 corrects the evaluation-timing terminology to a pre-score-locked
  audit and retains the frozen numerical results.
- Figure 4 retains the layout, palette and numerical point estimates. Its
  upper score-axis panel now uses signed `p−0.5`, so the asymmetric score
  abstention interval `(0.35, 0.70)` is represented correctly. The old absolute
  distance display incorrectly suggested a symmetric interval.
- Figure 4's untraceable historical confidence intervals were replaced by a
  documented calculation: 2000 full-set Murcko scaffold-cluster bootstrap
  draws, seed 20261018, 2312 clusters, including one acyclic connectivity-block
  fallback. All models and both fixed endpoint bands share each draw. Intervals
  use the 2.5th and 97.5th percentiles. The outcome-aware band remains fixed at
  `abs(log10(ER/2)) <= 0.49`; the intervals are conditional on this retrospective
  band. It is neither a prospective selection rule nor a guarantee of error.
- The retained-fraction difference remains −0.04407098940391574. Its replacement
  95% interval is `[−0.0742511782727669, −0.013685523854519778]`, displayed as
  **−4.4 percentage points [−7.4, −1.4]**.
- Figure 4 removes the incorrect assertion that pooled balanced accuracy must
  lie between subgroup balanced accuracies. Confusion counts now supply the
  arithmetic checks.
- Figure 5 uses the existing corrected 14-species predictions. It removes
  the historical incorrect two-sulfur S-methyl-DM4 input and unresolved
  S-methyl-DM1 illustration, identifies the cys-VC-MMAE intermediate by its full
  name, and marks the plotted Lys-MCC-DM1 stereoisomer as representative, distinct
  from the reported epimeric-mixture material. Its frozen
  decisions remain 5 accepted high-risk and 9 abstained, with no accepted
  low-risk species. The original palette and figure style are retained.

## Reproducibility boundary

**The public addendum supports redrawing Figures 3, 4 and 5 from the included
plot data and the original public archive. It does not support a standalone
recomputation of all Figure 4 confidence intervals.**

`handoff/fig4_data.json` includes plotted aggregates, confusion counts,
bootstrap parameters and input SHA256 hashes. `verification/figure4_plot_points.csv`
contains the 482 plotted summary points. `make_figure4.py` checks aggregate
arithmetic against the public frozen metric table when the per-molecule input
is absent; it explicitly reports that per-molecule validation and bootstrap
were **not** rerun. When that input is present, it performs the stricter
per-molecule checks.

The optional `rebuild_figure4_data.py` documents the full reconstruction but
requires these upstream inputs, which are **not redistributed in this addendum
or the base study archive**:

- `数据/processed/abcb1_public_predictor_predictions_long_v1.csv`
- `数据/processed/abcb1_public_predictor_benchmark_input_v1.csv`

Their identities are recorded in the handoff provenance and
`verification/input_provenance.json`. The local reconstruction verified all
unchanged historical rolling-error, measured-coverage, absolute-margin and
balanced-accuracy points against frozen per-molecule outputs at absolute
tolerance `1e-14`. An independent local reconstruction without the old handoff
produced identical new handoff bytes. This local input verification is distinct
from the public-archive redraw receipt.

The smooth score curve uses windows of width 0.05 score units; the measured
curves use width 0.15 log10 units. A window can straddle a decision boundary.
The statement that none of the 1689 molecules inside `(0.35,0.70)` were accepted
is an individual-level count, not an inference from smoothing.

## Contents, provenance and scope

`MANIFEST.json` records every packaged payload's SHA256; `SHA256SUMS.txt` lists
the same files for command-line verification. `verification/redraw_receipt.json`
records the actual overlay-on-public-archive redraw checks. Optional previews
contain standalone data figures, not manuscript pages.

No manuscript text, full manuscript PDF, author information file, proprietary
font, new per-molecule third-party predictor CSV, or model bundle is added.
The historical `handoff/fig5_data.json` in the base archive is not consumed by
the current Figure 5 generator; do not use it to regenerate the corrected ADC
figure. Existing data/structure attribution and restrictions remain those in
the base archive's `THIRD_PARTY_NOTICES.md` and
`REPOSITORY_THIRD_PARTY_NOTICES.md`. The included MIT license applies to
author-owned code and documentation, not to relabeling third-party records.
