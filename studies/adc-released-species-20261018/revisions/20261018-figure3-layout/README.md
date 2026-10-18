# Figure 3 panel-layout addendum, 18 October 2026

This second figure addendum updates the placement of Figure 3's panel a/b
annotations and comparison labels. It retains the scientific inputs, frozen
numerical results, gate definitions and palette. It contains the current
Figure 3 plotting tool, a standalone figure PDF, a numerical-guard log and
verification hashes. No manuscript text or full manuscript PDF is included.

Apply the releases in this order to a **fresh working copy**:

1. [Original ADC study archive](https://github.com/tommyverygood/EffluxScope/releases/tag/adc-released-species-20261018):
   `EffluxScope_ADC_reproducibility_20261018.zip`, SHA256
   `19fd294d893345af41edb08a1a8b86b9b9a4bfe7a859e50d52c44998be26938f`.
2. [First figure-correction addendum](https://github.com/tommyverygood/EffluxScope/releases/tag/adc-figure-corrections-20261018):
   `EffluxScope_figures_revision_20261018.zip`, SHA256
   `4f53ee9d1e6e3b17d43c48b822761fe99e2da5042df01925fe4987fc3042eac2`.
   Copy its `overlay/` contents into the extracted study root.
3. Copy this addendum's `overlay/` contents into that same study root.

From the study root, using the plotting dependencies and fonts described in
the first addendum:

```sh
python -c "from pathlib import Path; Path('论文主稿/latex_jcim/figures').mkdir(parents=True, exist_ok=True)"
python 论文主稿/latex_jcim/tools/make_figure3_jcim.py
```

The script writes the updated Figure 3 PDF/PNG/EPS. NumPy, Matplotlib and
Ghostscript are required. `palette.py`, `palette.json` and
`jcim_fig34_style.py` come from the first addendum. Fonts are not redistributed.

`verification/numerical_guard.log` records the actual guard check and redraw
on the original public archive plus both overlays. `verification/receipt.json`
records unchanged public-input SHA256 hashes, scientific-content equivalence
to the local summaries, source equivalence to the maintained local tool,
and text/pixel comparison with the current standalone
figure. This is a plotting and layout verification, not new model training,
scoring, bootstrap resampling or validation data.

The earlier releases and their assets remain unchanged. Figures 4 and 5 and
the limitations documented in the first addendum are unaffected. The software
license and third-party data attributions remain those of the original study
archive and its first addendum. `MANIFEST.json` and `SHA256SUMS.txt` describe
the payload files of this second addendum.
