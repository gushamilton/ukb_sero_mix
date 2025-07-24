# Step 2: Single-Antigen Phenotype Generation

This step is the core of the power-enhancement strategy. It processes the raw MFI data from Step 1 to generate a suite of "soft" (probabilistic) and "hard" (binary) phenotypes for each individual antigen.

## Objective

-   **Model Distributions**: Fit a two-component, skew-t mixture model to the log-MFI distribution of each selected antigen to handle the bimodal and overlapping nature of serology data.
-   **Generate Probabilistic Phenotypes**: Use the model fits to calculate the posterior probability of seropositivity (`p_soft`) for each individual and antigen.
-   **Create a Comprehensive Phenotype Suite**: For each antigen, generate multiple phenotype definitions to support different types of GWAS analyses (e.g., case/control vs. quantitative trait, weighted vs. unweighted).
-   **Derive a Novel Covariate**: Create a `latent_factor_1` covariate to control for shared technical noise across assays.

## Implementation

The script `ukb_potential_power/02_creating_pheno/02_create_serology_phenotypes.R` performs this entire process.

### 1. Core Antigen Selection

The selection of antigens for analysis is based on a combination of factors: evidence of bimodality in their distributions, previously reported genetic signals (primarily from Butler-Laporte et al., 2020), and the need for diverse pathogen representation.

The final list of **core antigens** is:

```
"cmv_pp150", "cmv_pp28", "cmv_pp52",
"hsv1", "hsv2",
"ebv_vca", "ebv_ebna1", "ebv_zebra", "ebv_ead",
"hp_omp", "hp_urea", "hp_caga", "hp_vaca", "hp_groel", "hp_catalase", # All HP
"toxo_p22",
"bkv_vp1", "jcv_vp1", "mcv_vp1",
"hhv6_ie1a", "hhv6_ie1b", "hhv7_u14",
"ct_pgp3", "ct_mompd", "ct_tarpf2", "ct_mompa", "ct_tarpf1", # All CT except porb
"kshv_lana"
```

### 2. Phenotype Definitions

For each core antigen, the script generates four primary phenotype definitions:

1.  **`{antigen}_IgG_raw`** (Quantitative): The natural log of the raw MFI value. This tests for variants affecting antibody levels across the whole population.
2.  **`{antigen}_sero_hard`** (Binary): A standard 0/1 case/control phenotype based on the official hard MFI cutoff. Used for direct replication of previous studies.
3.  **`{antigen}_sero_soft`** (Probabilistic): The posterior probability of being seropositive (`p_soft`), derived from the mixture model. This is our most powerful phenotype for detecting serostatus associations and is treated as a continuous trait.
4.  **`{antigen}_IgG_seropos_only`** (Quantitative, Conditional): The log-MFI value, but only for individuals classified as seropositive by the hard cutoff (it is `NA` for seronegatives). This replicates the quantitative analysis strategy of Butler-Laporte et al.

### 3. Weight Generation

The script also generates corresponding weights files to modulate the contribution of individuals in the GWAS:

-   **`{antigen}_w_igg`**: Weight is equal to `p_soft`. Used with the `_IgG_raw` phenotype to focus the analysis on likely seropositive individuals.
-   **`{antigen}_w_hard`**: Weight is `2 * |p_soft - 0.5|`. Used with the `_sero_hard` phenotype to up-weight confidently classified cases and controls and down-weight ambiguous ones.

### 4. Novel Covariate: `latent_factor_1` (formerly `PC1_assay_noise`)

A key innovation is the creation of a covariate to model shared technical noise.

-   **Rationale**: To control for non-biological factors (e.g., background antibody "stickiness", batch effects) that might affect an individual's baseline MFI readings across all antigens.
-   **Method**:
    1.  Identify **"confidently seronegative"** individuals (`p_soft < 0.1`) for each antigen.
    2.  Create a matrix of their `log(MFI)` values, masking out values for any antigen where they are not confidently negative.
    3.  Impute the matrix and run a Principal Component Analysis (PCA).
-   **Result**: The first principal component, `latent_factor_1`, represents the major axis of shared variance among seronegative individuals and is used as a covariate to increase the precision of all association tests.

## Outputs

This script generates the core input files for the downstream GWAS analyses, all placed in the `quickdraws_input/` directory:

-   **`phenotypes_master_postqc.tsv`**: The master table containing all generated phenotype columns for all individuals.
-   **`covariates_*.tsv`**: Covariate files, including a base set (age, sex, genetic PCs) and the novel `latent_factor_1`.
-   **`analysis_manifest.tsv`**: A machine-readable file that defines every GWAS to be run, specifying the phenotype, analysis type (linear/logistic), and any associated weights file.
-   **`*.weights`**: A series of weight files, one for each weighted analysis defined in the manifest.
-   **`mixture_model_fits.rds`**: An R object containing the fitted `mixsmsn` model for each antigen, for diagnostic purposes.
-   **QC Reports**: Various plots and tables for quality control, such as phenotype correlations and lists of excluded traits. 