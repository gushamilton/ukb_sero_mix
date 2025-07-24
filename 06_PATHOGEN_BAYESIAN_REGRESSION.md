# Step 6: Pathogen-Level Bayesian Regression

This step replaces the fragile multivariate‐mixture attempts with a robust, fully Bayesian logistic-regression model that learns data-driven weights for each antigen.

## Rationale

1. Mixture models on raw MFI are powerful but can fail to converge for multi-antigen pathogens.
2. We already have reliable single-antigen posterior probabilities (`*_sero_soft`) from Step 2.
3. Treating those probabilities as predictors lets us combine them in a principled way: a Bayesian logistic regression learns how much each antigen contributes to the latent *pathogen* infection state.
4. Using the noisy UK Biobank “hard rule” (≥2 antigens positive, or Pgp3 positive for CT) as the outcome gives the model a weakly-informative target while still allowing it to smooth over rule mis-classifications.

## Model specification (per pathogen)

```
z_j   ~  Bernoulli(π_j)
logit(π_j) = β₀  + Σ_k β_k · p_{jk}
β_k  ~ Normal(0, 1)
```

*   `p_{jk}`  … per-antigen posterior probability from Step 2.
*   `z_j`     … noisy hard status defined by the UKB rule.
*   Weak Normal(0,1) priors are sufficient; the predictors are already on [0,1].

## Implementation Highlights (`06_pathogen_level_phenotypes_bayesReg.R`)

1. Reads `phenotypes_master_postqc.tsv` for all `_sero_soft` and `_sero_hard` columns.
2. Builds the UKB hard label per pathogen.
3. Fits the above model with **brms** (Stan backend).
4. Extracts the posterior mean \(\hat π_j\) as the new
   * `{pathogen}_sero_soft`
   * `{pathogen}_sero_hard` = 1 if \(\hat π_j ≥ 0.5\).
5. Writes
   * `phenotypes_pathogen_bayesReg.tsv`
   * weight files `{pathogen}_p_soft.weights` and `{pathogen}_sero_hard.weights`.

## Development vs. Production

* While prototyping the script subsamples **200** complete cases per pathogen for speed. Comment-out the `slice_sample()` line to run on the full cohort.
* Typical run-time on the subsample: ~1–2 minutes per pathogen on a laptop; full data will scale linearly.

## Advantages

* Fully Bayesian: uncertainty in β propagates into per-sample probabilities.
* Automatically down-weights “weak” antigens (small β) and amplifies informative ones.
* No convergence issues observed; logistic regression is easy for HMC.

## Next Steps

* Inspect posterior summaries (`summary(fit)`) to see which antigens carry the most weight.
* Run posterior predictive checks: `pp_check(fit)`.
* Remove subsampling and run on full UKB cohort before feeding phenotypes into Quickdraws GWAS. 