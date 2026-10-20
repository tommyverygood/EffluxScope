# ADC study source-audit addendum — 20 October 2026

This addendum corrects one historical source diagnostic and completes plotting-source distribution. It supplements the ADC payload/released-species study; it is not a new trained predictor or new experimental validation. No full manuscript, author-contact records, third-party prediction arrays, model checkpoints, article full texts, or fonts are included.

## Changes and verified scope

- **Figure 3b:** the Roth panel uses a cytotoxicity-derived GI50 resistance ratio, not a measured monolayer transport efflux ratio. The current display excludes the source material named Deruxtecan because its chemical identity was unresolved. The recorded in-house Morgan baseline is evaluated on 12 admitted compounds: Spearman rho = -0.570929, p = 0.052524. This retrospective identity audit supersedes the old 13-row display; it does not turn the diagnostic into a new scoring-locked transport validation. The negative result is retained and no confidence interval is invented.
- **Plot guard:** Figure 3 now checks its Roth n/rho/p and endpoint/exclusion against the identity-audited record, and checks its ASAP interval against panel a. The old 13-row/transport-label combination is rejected. Panels a, c and d, their inputs, the model decisions and the existing palette are unchanged.
- **Figures 1/2:** the maintained plotting scripts and the missing 443-by-2 Figure 2 reference-coordinate array are provided. Figure 2 was redrawn using only anonymously downloaded public archives plus this supplement and matched the maintained local figure pixel for pixel under the recorded environment.
- **Remaining limits:** complete Figure 1 redraw still needs excluded upstream per-molecule MDCK-ER outputs. Figure 4's public aggregate data support redraw but not complete independent recomputation of its bootstrap intervals. Historical bootstrap/training runs are not newly reproduced here. See `README_FIGURES12.md` for exact inputs and limits.

## Apply in order

Work in a **fresh extracted copy** of the [core study archive](https://github.com/tommyverygood/EffluxScope/releases/tag/adc-released-species-20261018), whose root is `adc-released-species-20261018`.

1. Copy the **contents** of `overlay/` from the [first figure addendum](https://github.com/tommyverygood/EffluxScope/releases/tag/adc-figure-corrections-20261018) into that root.
2. Apply the [Figure 3 layout addendum](https://github.com/tommyverygood/EffluxScope/releases/tag/adc-figure3-layout-20261018) in the same way.
3. Apply this addendum's `overlay/` contents last. Do not nest an extra `overlay` folder inside the study root.

The earlier immutable releases remain available for provenance; their Figure 3b n=13 display is superseded by this correction. From the overlaid study root, using the documented Python plotting environment and Ghostscript:

```sh
python -c "from pathlib import Path; Path('论文主稿/latex_jcim/figures').mkdir(parents=True, exist_ok=True)"
python 论文主稿/latex_jcim/tools/make_figure2_jcim.py
python 论文主稿/latex_jcim/tools/make_figure3_jcim.py
python outputs/20261020_roth_endpoint_audit/reproduce_roth_identity_safe.py
```

The last command recomputes n/rho/p from 12 published functional observations and previously recorded **in-house** baseline scores. It does not load a model. All revised scores, outcomes and exclusion semantics are described alongside that CSV. Figure 4/5 commands and dependencies remain those in the first addendum; they use unchanged plotting inputs. `requirements-figures12.txt` records the successfully tested plotting environment; SciPy is required for the Roth association check.

## Provenance and licensing

Author-owned plotting/reproduction code uses the project MIT license. ChEMBL-derived visualization coordinates retain their source attribution and applicable CC BY-SA 3.0 terms. The small Roth Table 3 extract retains CC BY 4.0 attribution. Those terms are scoped to their respective materials, not replaced by the software license. See `THIRD_PARTY_NOTICES_FIGURES12.md` and `overlay/outputs/20261020_roth_endpoint_audit/README.md`.

`MANIFEST.json` and `SHA256SUMS.txt` identify all payload files. Verification receipts record actual redraw and numerical checks, including their limits. Successful verification does not establish ADC-specific accuracy, matched parent–released-species transport truth, intact-ADC efficacy, or journal submission.
