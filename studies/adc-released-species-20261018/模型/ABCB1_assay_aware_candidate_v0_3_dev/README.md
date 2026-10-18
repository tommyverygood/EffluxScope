# Run

```bash
uv run --python 3.12 \
  --with scikit-learn==1.9.0 \
  --with rdkit==2026.3.4 \
  --with scipy==1.18.0 \
  --with joblib==1.5.3 \
  --with numpy==2.5.1 \
  predict_abcb1_assay_aware.py --smiles 'CCO'
```

Batch input requires a `smiles` column unless `--smiles-column` is supplied:

```bash
uv run --python 3.12 \
  --with scikit-learn==1.9.0 \
  --with rdkit==2026.3.4 \
  --with scipy==1.18.0 \
  --with joblib==1.5.3 \
  --with numpy==2.5.1 \
  predict_abcb1_assay_aware.py --input-csv queries.csv --output predictions.csv
```

This is a single-source development candidate, not an externally validated universal ABCB1 substrate predictor.
