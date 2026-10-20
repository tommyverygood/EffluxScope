# Figure 1/2 plotting-source completion, 20 October 2026

This component supplies the maintained Figure 1/2 plotting code and the frozen Figure 2 reference coordinates that were absent from the earlier public study export. It does not change the models, observations, selected sets, estimates, uncertainty intervals, or manuscript figures. The files belong in a fresh working copy of the original public study archive. This component does not contain the manuscript, author-contact tables, third-party model weights, or third-party per-molecule predictor outputs.

## Base releases and overlay order

Download these public assets anonymously from their versioned release pages:

1. [Core study](https://github.com/tommyverygood/EffluxScope/releases/tag/adc-released-species-20261018): `EffluxScope_ADC_reproducibility_20261018.zip`, SHA256 `19fd294d893345af41edb08a1a8b86b9b9a4bfe7a859e50d52c44998be26938f`.
2. [First figure addendum](https://github.com/tommyverygood/EffluxScope/releases/tag/adc-figure-corrections-20261018): `EffluxScope_figures_revision_20261018.zip`, SHA256 `4f53ee9d1e6e3b17d43c48b822761fe99e2da5042df01925fe4987fc3042eac2`.
3. [Figure 3 layout addendum](https://github.com/tommyverygood/EffluxScope/releases/tag/adc-figure3-layout-20261018): `EffluxScope_Figure3_layout_20261018.zip`, SHA256 `788c042c2b2bd7fb84a001a418139a033deb44fcb42f65dd08eb299f00d11887`.
4. This component's `overlay/` contents, preserving relative paths.

The extracted core study root is `adc-released-species-20261018`. Copy each addendum's **overlay contents** into that root, in the order above. Do not place an extra `overlay` directory inside the study root. Keep the original archives unchanged. If the enclosing addendum also includes a corrected Figure 3 component, its own instructions supersede the historical Figure 3 layout while leaving this Figure 2 command unchanged.

## Redraw Figure 2 from the public inputs

Use Python 3.12 with the versions in `requirements-figures12.txt`, and install Ghostscript so that `gs` is on PATH. Arial is the preferred Figure 2 font; DejaVu Sans is its documented fallback. Proprietary fonts are not included. Font substitution and library differences can change text metrics even with identical coordinates and estimates.

From the overlaid study root:

```sh
python -c "from pathlib import Path; Path('论文主稿/latex_jcim/figures').mkdir(parents=True, exist_ok=True)"
python 论文主稿/latex_jcim/tools/make_figure2_jcim.py
```

Outputs are `论文主稿/latex_jcim/figures/Figure2_split_separation.pdf`, `.png`, and `.eps`. The script checks the five models' full-set balanced accuracy and AUROC against the public frozen metric table before drawing. It reads these inputs:

| Relative path | Source |
|---|---|
| `handoff/fig2_ref_xy.npy` | This component; 443 × 2 float32 reference t-SNE coordinates |
| `handoff/fig2_emb.csv` | Core study; external coordinates and split membership |
| `handoff/fig2_perf.csv` | Core study; aggregate performance values |
| `handoff/fig2_leak.json` | Core study; aggregate reference-similarity summaries |
| `handoff/fig2_hist.json` | Core study; histogram and structural-threshold information |
| `handoff/fig2_ci.json` | Core study; recorded confidence intervals |
| `数据/processed/abcb1_public_predictor_metrics_v1.csv` | Core study; frozen aggregate-metric crosscheck |
| `论文主稿/figures_src/palette.py` and `palette.json` | Included unchanged for style consistency |

The new coordinate array contains no molecular prediction values, outcome labels, or third-party service outputs. These are the authors' frozen visualization coordinates derived from the reference structures. They support the current Figure 2 drawing; this command does **not** rerun t-SNE, regenerate bootstrap intervals, fit a model, or establish independent biological validation. Preserve the underlying source-data attribution described in `THIRD_PARTY_NOTICES_FIGURES12.md`.

## Figure 1: source provided, complete public redraw still limited

The maintained scripts are supplied at:

- `论文主稿/latex_jcim/tools/make_figure1_jcim.py`
- `论文主稿/figures_src/make_figure1_v2.py`

The wrapper retains the current figure's typography and stroke widths. It reads `handoff/fig1ab.json`, `fig1ab_grid.json`, `fig1new_stats.json`, `fig1_permband.json`, and the frozen aggregate metric table. Those public files support the internal-model scatter and aggregate panels. However, the published core export deliberately omits the AstraZeneca GNN-MTL MDCK-ER per-molecule scatter arrays `astrazeneca_gnn_mtl_mdck_er_2025|Q1`, `|Q3`, and `|Q5` required by the complete current Figure 1. **Adding these scripts does not make the full Figure 1 independently redrawable from this public package alone.** No excluded per-molecule arrays have been copied into this component.

To reconstruct those panels, obtain the exact external predictor/checkpoint from its cited provider under its terms, recover the recorded external query structures and preprocessing/stratum assignments, and regenerate the historical comparator outputs. The existing historical analysis scripts, model/status summaries, aggregate values, and source locators document the comparison; they are not a guaranteed one-command recreation of the external service or its historical runtime. Do not substitute a different model/version and call the result the original figure. After lawful recovery and verification of the required inputs, the full drawing command is:

```sh
python 论文主稿/latex_jcim/tools/make_figure1_jcim.py
```

Without the missing external arrays, that command is expected to fail at the missing scatter keys. This explicit limit is distinct from the supported public Figure 2 redraw.

## Figure 4 and historical-analysis limits

The first figure addendum supplies public aggregate data sufficient to redraw Figure 4. Recomputing all of its scaffold-bootstrap intervals still requires upstream per-molecule comparator inputs that are not redistributed. This component neither supplies those inputs nor changes that limitation. Figure 2 likewise uses its recorded intervals rather than recalculating them. See the first addendum's README and input provenance for the exact bootstrap/data boundary.

The core offline verification command remains a separate workflow for frozen ADC scores, point estimates, descriptors and reporting gates. It does not reproduce every historical collection, training or external-model comparison run. Successful plotting is not proof of ADC-specific accuracy, matched parent–released-species transport truth, or complete-ADC efficacy.

## Verification and provenance

`verification/figures12_source_provenance.json` records the exact copied-source hashes and coordinate shape. `verification/figure2_public_redraw_receipt.json` records the actual redraw against anonymously downloaded public archives, input hashes, frozen-metric checks and comparison with the maintained local Figure 2. A successful receipt establishes the checks performed; the instructions alone do not.

This component belongs to the 20 October source-audit addendum. Publication of the addendum is separate from journal submission; no journal submission is recorded here.
