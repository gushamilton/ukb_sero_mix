# Step 6: Pathogen-Level Bayesian Regression

This step replaces the fragile multivariate‐mixture attempts with a robust, fully Bayesian logistic-regression model that learns data-driven weights for each antigen.

## Rationale

1. Mixture models on raw MFI are powerful but can fail to converge for multi-antigen pathogens.
2. We already have reliable single-antigen posterior probabilities (`*_sero_soft`) from Step 2.
3. Treating those probabilities as predictors lets us combine them in a principled way: a Bayesian logistic regression learns how much each antigen contributes to the latent *pathogen* infection state.
4. Using the noisy UK Biobank "hard rule" (≥2 antigens positive, or Pgp3 positive for CT) as the outcome gives the model a weakly-informative target while still allowing it to smooth over rule mis-classifications.

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
2. Builds the UKB hard label per pathogen using exact Butler-Laporte et al. 2023 rules.
3. Fits the above model with **brms** (Stan backend) on the full dataset.
4. Extracts the posterior mean \(\hat π_j\) as the new
   * `{pathogen}_sero_soft`
   * `{pathogen}_sero_hard` = 1 if \(\hat π_j ≥ 0.5\).
5. Generates comprehensive outputs for GWAS comparison.

## Enhanced Features

### Credible Intervals
* **Individual-level CIs**: Each soft probability comes with 95% credible intervals (`*_sero_soft_lower`, `*_sero_soft_upper`).
* **CI width**: Uncertainty measure (`*_sero_soft_ci_width`) for quality assessment.
* **Population-level CIs**: Posterior-based credible intervals for seroprevalence estimates using `posterior_epred()`.

### Multiple Serostatus Thresholds
The script generates multiple binary classifications:
* `*_sero_hard_30`: ≥30% probability threshold (sensitive)
* `*_sero_hard_50`: ≥50% probability threshold (balanced)
* `*_sero_hard_70`: ≥70% probability threshold (specific)
* `*_sero_hard_90`: ≥90% probability threshold (very specific)

### Butler-Laporte Comparison
* **Exact MFI thresholds**: Uses precise thresholds from Butler-Laporte et al. 2023 for direct paper comparison.
* **Exact seropositivity rules**: Implements the exact rules from their Table 1.
* **Comparison outputs**: Generates `{pathogen}_bl_hard.pheno` files for direct comparison.

### Comprehensive Weights
* **Soft phenotype weights**: `{pathogen}_bayes_soft.weights` using soft probability as weight.
* **Hard phenotype weights**: `{pathogen}_bayes_hard.weights` using confidence (2 × |soft_prob - 0.5|) as weight.
* **No weights for hard cutoffs**: Butler-Laporte and UKB hard cutoffs are deterministic.

### Model Diagnostics
* **Coefficient summaries**: Posterior means and credible intervals for each antigen weight.
* **Seroprevalence estimates**: With posterior-based credible intervals.
* **Mean CI width**: Average uncertainty across individuals.

## Outputs

### Phenotype Files
* `phenotypes_pathogen_bayesReg.tsv`: Main phenotype file with soft probabilities, CIs, and multiple thresholds.
* `phenotypes_gwas_comparison.tsv`: Combined file with all phenotype types for comparison.

### Individual Phenotype Files
* `{pathogen}_ukb_hard.pheno`: UKB rule hard cutoffs (noisy target).
* `{pathogen}_bayes_hard.pheno`: Bayesian hard cutoffs (≥0.5).
* `{pathogen}_bayes_soft.pheno`: Bayesian soft probabilities.

### Weight Files
* `{pathogen}_bayes_soft.weights`: Weights for soft phenotype GWAS.
* `{pathogen}_bayes_hard.weights`: Weights for hard phenotype GWAS.

### Diagnostic Files
* `results/bayesian_regression_seroprevalence_summary.tsv`: Summary statistics.
* `results/gwas_comparison_summary.tsv`: Comparison table with all seroprevalence estimates.
* `results/{pathogen}_antigen_vs_pathogen_scatters.pdf`: Scatter plots comparing individual antigens to pathogen-level probabilities.

## Butler-Laporte et al. 2023 Implementation

### MFI Thresholds (from Table 1)
* **EBV**: VCA p18 (250), EBNA-1 (250), ZEBRA (100), EA-D (100)
* **CMV**: pp150 (100), pp52 (150), pp28 (200)
* **H. pylori**: CagA (400), VacA (100), OMP (170), GroEL (80), Catalase (180), UreA (130)
* **C. trachomatis**: momp A (100), momp D (100), tarp-D F1 (100), tarp-D F2 (100), PorB (80), pGP3 (200)

### Seropositivity Rules
* **EBV**: Positive for 2 or more antigens
* **CMV**: Positive for 2 or more antigens
* **H. pylori**: Positive for 2 or more antigens, **except CagA**
* **C. trachomatis**: Positive for pGP3 **OR** negative for pGP3 but positive for 2 out of 5 remaining antigens

## Advantages

* **Fully Bayesian**: Uncertainty in β propagates into per-sample probabilities.
* **Multiple thresholds**: Flexible serostatus classification for different use cases.
* **Comprehensive diagnostics**: Model validation and uncertainty quantification.
* **Direct paper comparison**: Exact implementation of Butler-Laporte methodology.
* **Confidence-based weights**: Maximizes GWAS power by weighting uncertain calls appropriately.
* **No convergence issues observed**; logistic regression is easy for HMC.

## Interpretation

### Coefficient Interpretation
* Positive coefficients indicate the antigen increases pathogen seropositivity probability.
* Larger coefficients indicate stronger predictive power.
* Credible intervals show uncertainty in antigen importance.

### Threshold Selection
* **30% threshold**: High sensitivity, includes borderline cases.
* **50% threshold**: Balanced, traditional cutoff.
* **70% threshold**: High specificity, confident seropositive calls.
* **90% threshold**: Very high specificity, excludes uncertain cases.

### Weight Interpretation
* **Soft weights**: Direct probability weighting for soft phenotype GWAS.
* **Hard weights**: Confidence-based weighting (distance from 0.5) for hard phenotype GWAS.
* Higher weights indicate more confident predictions.

### CI Width Interpretation
* Narrow CIs indicate high confidence in the prediction.
* Wide CIs suggest the model is uncertain, possibly due to conflicting antigen signals.

## Next Steps

* Inspect posterior summaries (`summary(fit)`) to see which antigens carry the most weight.
* Compare seroprevalence estimates with Butler-Laporte paper results.
* Use CI width to identify individuals with uncertain classifications.
* Run GWAS with both weighted and unweighted approaches to assess power gains.
* Compare performance across different serostatus thresholds. 