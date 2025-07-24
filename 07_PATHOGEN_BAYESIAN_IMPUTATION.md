# Step 7: Pathogen-Level Bayesian Regression with Missing Data Imputation

This step extends the Bayesian regression approach to handle missing antigen values directly in the model, allowing us to use ALL individuals rather than dropping those with incomplete data.

## Why Bayesian Imputation?

### Current Limitation (Script 06)
- Drops individuals with any missing antigens for a pathogen
- Loses potentially valuable data
- May introduce bias if missingness is not random

### Bayesian Solution (Script 07)
- Models missing values as parameters in the Bayesian model
- Uses correlations between antigens to inform imputation
- Preserves uncertainty in imputed values
- Uses ALL available data

## How It Works

### 1. Missing Data as Parameters
Instead of treating missing values as "unknown", we treat them as parameters to be estimated:

```
z_j ~ Bernoulli(π_j)
logit(π_j) = β₀ + Σ_k β_k · p_{jk}
β_k ~ Normal(0, 1)

# For missing p_{jk}:
p_{jk} ~ Normal(μ_k, σ_k)  # Prior for missing values
```

### 2. Imputation During MCMC
- Missing values are sampled during MCMC iterations
- Each iteration provides one possible imputation
- Posterior distribution includes uncertainty in imputed values
- Predictions naturally account for imputation uncertainty

### 3. Key Advantages
- **Principled**: Missing data mechanism is part of the model
- **Uncertainty**: Imputed values have proper uncertainty quantification
- **Efficiency**: Uses all available information
- **Robust**: Less sensitive to missing data patterns

## Implementation Details

### Missing Data Handling
```r
# Keep individuals with at least 1 non-missing antigen
filter(rowSums(!is.na(select(., all_of(soft_cols)))) >= 1)

# Handle missing hard calls in target variable
hard_matrix[is.na(hard_matrix)] <- 0
```

### brms Configuration
```r
fit <- brm(formula, data = df_model,
           missing = "mi",  # Enable missing data imputation
           chains = 2, iter = 2000, warmup = 500)
```

### Prediction with Imputation
```r
# Predict on ALL individuals (including those with missing data)
pred_full <- fitted(fit, newdata = df, scale = "response", probs = c(0.025, 0.25, 0.75, 0.975))
```

## Expected Benefits

### 1. Increased Sample Size
- Use individuals with partial antigen data
- Potentially significant power gains for rare pathogens

### 2. Better Uncertainty Quantification
- Credible intervals reflect imputation uncertainty
- More realistic uncertainty estimates

### 3. Reduced Bias
- Avoids bias from complete-case analysis
- Uses all available information

## Outputs

### Phenotype Files
- `phenotypes_pathogen_bayesReg_imputed.tsv`: Main phenotype file with imputed values
- `{pathogen}_p_soft_imputed.weights`: Weight files for soft phenotypes
- `{pathogen}_sero_hard_imputed.weights`: Weight files for hard phenotypes

### Diagnostic Information
- Missing data summaries per pathogen
- Comparison of sample sizes (with vs. without imputation)
- Uncertainty quantification for imputed values

## Comparison with Script 06

| Aspect | Script 06 (Complete Cases) | Script 07 (With Imputation) |
|--------|---------------------------|----------------------------|
| Sample Size | Reduced (complete cases only) | Full (all individuals) |
| Missing Data | Dropped | Imputed with uncertainty |
| Uncertainty | Underestimated | Properly quantified |
| Bias Risk | Higher | Lower |
| Computational Cost | Lower | Higher |

## When to Use Each Approach

### Use Script 06 (Complete Cases) When:
- Missing data is minimal (<5%)
- Computational resources are limited
- You want conservative estimates

### Use Script 07 (With Imputation) When:
- Missing data is substantial (>10%)
- You want maximum power
- You have computational resources
- Missingness is likely informative

## Next Steps

1. **Compare Results**: Run both scripts and compare seroprevalence estimates
2. **Assess Power Gain**: Compare sample sizes and expected power
3. **Validate Imputation**: Check if imputed values make biological sense
4. **GWAS Analysis**: Use imputed phenotypes in Quickdraws GWAS 