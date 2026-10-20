# Roth 2026 endpoint and identity correction

This small addendum reproduces the corrected Figure 3b association for 12 identity-admitted compounds. It uses previously recorded **in-house Morgan-fingerprint logistic baseline scores**, not EffluxScope v0.1 scores and not third-party predictions. No model is loaded, trained or scored by this reproduction.

## Run

With Python and SciPy installed (the study reproduction environment already includes SciPy):

```sh
python reproduce_roth_identity_safe.py
```

Expected n = 12, Spearman rho = -0.5709290716957522, two-sided asymptotic p = 0.05252355548220434. Average ranks handle the tied RR values. The p value is the historical asymptotic Spearman calculation; no exact permutation test or bootstrap interval is added.

## Endpoint and exclusion

The numeric endpoint is the **GI50 resistance ratio** between ABCB1-expressing MDR19 and control pcDNA cells, from Roth et al. Table 3. It is a functional cytotoxicity endpoint, **not a monolayer transport efflux ratio (ER)**. These values are not assigned as direct-transport ADC classification truth. A negative association here supplies no positive directional support and cannot validate transport prediction.

The source-specific material named **Deruxtecan** is excluded because its exact experimental identity was unresolved (C26 DXd versus C52 linker-payload interpretation). Its biological observation is not erased or reassigned; it is omitted from structure-based correlation. This is an identity-audited retrospective sensitivity analysis, not a new scoring-locked or outcome-blind validation.

Only the 12 required material names, RR values, our scores and source/identity flags are included. This is a selected and column-reduced data extract; manuscript outcomes, predictions and assay values were not recomputed or altered. The diagnostic does not establish ADC-specific accuracy or intact-ADC effects.

## Attribution and licences

Source: Roth JS, Guo H, Chen L, Shen M, Gbadegesin O, Robey RW, Gottesman MM, Hall MD. Identification of antibody-drug conjugate payloads that are substrates of ATP-binding cassette drug efflux transporters. Cancer Drug Resist. 2026;9:2. DOI: [10.20517/cdr.2025.151](https://doi.org/10.20517/cdr.2025.151). [Publisher article and Table 3](https://www.oaepublish.com/articles/cdr.2025.151).

The publisher's article-specific Copyright section was checked on 2026-10-20 and identifies copyright as The Author(s) 2026, licensed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/). The Table 3 material-name/RR extract retains CC BY 4.0 attribution and terms. Changes: selected 12 identity-admitted rows and the ABCB1 RR column, renamed columns, and appended our frozen scores and audit flags. The underlying article is not redistributed.

Our original reproduction script, diagnostic documentation and in-house derived score contributions are offered under the project's MIT terms (see `LICENSE-MIT.txt`); this does not relicense the source-derived table values. No third-party model, checkpoint, service output, or article full text is included.
