# Step 3: Pathogen-Level Phenotype Generation

This step builds upon the single-antigen probabilities generated in Step 2. For pathogens with multiple measured antigens (like EBV, CMV, HP), it combines the evidence from all of them to produce a single, more robust pathogen-level serostatus call.

## Motivation

The DKFZ/UKB validation paper provides hard cut-offs and logical rules (e.g. “positive for ≥2 of 4 antigens means EBV-seropositive”). These rules, while useful, ignore:

*   Continuous antibody titre information
*   Uncertainty near the cut-offs
*   Correlation between antigens for the same pathogen

This script replaces those hard rules with a unified probabilistic model.

## Statistical Model

### 2.1 Latent-Class Mixture

For a pathogen with \(m\) antigens and participant \(j\), let:

*   \(\mathbf y_j = (y_{j1},\ldots,y_{jm})\) be the vector of log-MFI values (the `*_IgG_raw` columns).
*   \(z_j \in \{0,1\}\) denote the unobserved true infection state for the pathogen.

We assume a **two-class multivariate Gaussian mixture model** with a full (unstructured) covariance matrix, fitted using the `mclust` R package:

\[
\mathbf y_j \mid z_j = 0 \;\sim\; \mathcal N(\boldsymbol\mu_0, \mathbf\Sigma_0),\qquad
\mathbf y_j \mid z_j = 1 \;\sim\; \mathcal N(\boldsymbol\mu_1, \mathbf\Sigma_1).
\]

**Key Points:**

1.  **Full Covariance \(\mathbf\Sigma_z\)**: The model learns the antigen-antigen correlations directly from the data, avoiding any assumption of independence.
2.  **Seeded from Step 2**: The EM algorithm is initialized using the single-antigen probabilities (`p_soft`) to ensure faster and more stable convergence. A participant's initial class is determined by whether their average `p_soft` across the pathogen's antigens is > 0.5.
3.  **Robustness**: The implementation includes safeguards to handle cases where an antigen has zero variance (by dropping it) or where the `mclust` algorithm fails for a specific pathogen (by catching the error and skipping that pathogen).

### 2.2 Pathogens to Model

The script applies this multivariate model to the following pathogens, using the specified antigen sets. Note that antigens that failed in the upstream single-antigen modeling will be silently skipped.

| Pathogen | Antigen set |
| :--- | :--- |
| EBV | `ebv_vca`, `ebv_ebna1`, `ebv_zebra`, `ebv_ead` |
| CMV | `cmv_pp150`, `cmv_pp52`, `cmv_pp28` |
| HHV-6 | `hhv6_ie1a`, `hhv6_ie1b`, `hhv6_p101k` |
| C. trachomatis | `ct_pgp3`, `ct_mompa`, `ct_mompd`, `ct_tarpf1`, `ct_tarpf2` |
| H. pylori | `hp_caga`, `hp_vaca`, `hp_omp`, `hp_groel`, `hp_catalase`, `hp_urea` |

## Implementation

The script `ukb_potential_power/02_creating_pheno/04_pathogen_level_phenotypes.R` performs this analysis.

## Outputs

For each pathogen modeled, the script generates a new set of phenotype and weight columns, analogous to the single-antigen ones:

*   **`{pathogen}_sero_soft`**: The final posterior probability of being infected with the pathogen.
*   **`{pathogen}_sero_hard`**: The binary (0/1) call based on whether `{pathogen}_sero_soft` is ≥ 0.5.
*   **`{pathogen}_w_soft`**: Weight file for use with quantitative traits.
*   **`{pathogen}_w_hard`**: Weight file for use with binary traits.

These are saved to a new file, `quickdraws_input/phenotypes_pathogen.tsv`, and the corresponding weight files are also written to the `quickdraws_input` directory. A file containing the fitted `mclust` objects (`pathogen_mixture_model_fits.rds`) is also saved for diagnostic purposes. 