# Figure 6: scaffold bootstrap reproducibility

This component recomputes the point estimates and confidence intervals in the current Figure 6 from the authors' frozen scaffold-level sufficient statistics. It requires no network access, trained model, service response, molecular structure, or private input.

Run from this directory:

```bash
python -m pip install numpy==2.5.1
python recompute_intervals.py --output verification.json
```

The released verification used Python 3.12.14 and NumPy 2.5.1. All 100 checked numeric fields match the approved `handoff/fig4_data.json` values exactly (maximum absolute difference 0.0). The historical source name `fig4_data` predates renumbering to Figure 6. `expected_figure6_intervals.json` is used only for verification; estimates and intervals are computed from the count matrix.

## Files and scope

- `figure6_scaffold_sufficient_statistics.csv`: 2,312 ordered scaffold groups and 28 integer count columns, plus an ordinal `cluster_index`.
- `provenance.json`: input hashes, group ordering, frozen band, model identifiers and resampling settings.
- `recompute_intervals.py`: NumPy-only recomputation and input-integrity checks.
- `expected_figure6_intervals.json`: approved point estimates and intervals used as a comparison target.
- `verification.json`: released verification result, numerical outputs and runtime versions.

For each side of the fixed measured-endpoint band, the matrix supplies total (`n`), positive (`pos`), negative (`neg`), accepted (`accepted`), and model-specific true-positive (`tp`) and true-negative (`tn`) counts. `in` means `abs(log10(measured ER) - log10(2)) <= 0.49`; `out` means its complement. There are 2,624 compounds, 1,349 inside the band and 453 accepted compounds in total. The band is outcome-aware and held fixed during resampling; it is not a prospective rejection rule.

Each of 2,000 draws samples 2,312 scaffold groups with replacement using `numpy.random.default_rng(20261018)`. Group order is the original sorted `numpy.unique(bootstrap_cluster)` order, preserved as ordinals. The same draw is used for every model and both bands. Balanced accuracy is `(TP / positive + TN / negative) / 2`. Undefined single-class draws are excluded. Limits are the 2.5th and 97.5th percentiles. The accepted-fraction contrast is the inside fraction minus the outside fraction.

This reproduces the statistical resampling stage and the corresponding point estimates. It does not regenerate or independently audit the upstream predictions, refit models, reconstruct molecular identities, or reproduce the rolling curves in panels a and c. Those curves are supplied by the main figure package. No per-molecule probability scores, SMILES, InChIKeys, original scaffold strings or service job identifiers are included here. Some groups are singletons: the count representation is not a guarantee that group-level binary outcomes cannot be inferred.

## Comparator provenance

The AstraZeneca MDCK-ER and NIH-MDCK-ER scores were computed locally using the [GNN-MTL Chemprop 2.1.0 checkpoint](https://zenodo.org/records/16948542). The [record's machine-readable license](https://zenodo.org/api/records/16948542) is Apache-2.0. The checkpoint is not redistributed in this component.

The ADMET-AI Pgp_Broccatelli inhibition comparator was computed locally with ADMET-AI 2.0.1. The [official repository license](https://github.com/swansonk14/admet_ai/blob/main/LICENSE.txt) is MIT. This is an endpoint-mismatched negative control, not a substrate classifier.

Deep-PK substrate predictions were obtained from the public service on 2026-08-05. The [official help](https://biosig.lab.uq.edu.au/deeppk/help) and [API documentation](https://biosig.lab.uq.edu.au/deeppk/api_docs) describe prediction and result download. No explicit redistribution license for individual service outputs was identified in the examined official pages. This component supplies the authors' evaluation counts, not the service's probability-score table. Regenerating contemporary service outputs is a separate operation and cannot be claimed to reproduce the historical deployment exactly.

Measured endpoints derive from the [Computational-ADME repository](https://github.com/molecularinformatics/Computational-ADME), whose [repository license is MIT](https://github.com/molecularinformatics/Computational-ADME/blob/main/LICENSE). The internal score and acceptance mask are the authors' frozen outputs. The upstream benchmark was post hoc; publishing its analysis inputs does not change its validation status.

## Licenses

The authors' recomputation code is covered by the included MIT `LICENSE`. The authors' count matrix and schema are offered under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); cite this manuscript and release when reusing them. This grant does not relicense upstream models, source datasets or service outputs. Their attribution and applicable conditions remain as documented above and in the release's third-party notices.
