# Comprehensive Phenotype Generation for Serology GWAS

This document describes the methodology implemented in the `02_create_serology_phenotypes.R` script. The script's primary purpose is to process raw UK Biobank serology data and generate a suite of phenotype and covariate files suitable for a comprehensive genome-wide association study (GWAS) using Quickdraws.

## 1. Core Objective

The goal is to test for genetic associations with pathogen exposure using multiple, complementary phenotype definitions. This allows us to:
1.  **Compare directly** with previous results (e.g., Butler-Laporte et al.) by using standard hard-cutoff case/control definitions.
2.  **Increase statistical power** by leveraging a two-component mixture model to create "soft" or probabilistic phenotypes, which better handle the ambiguity of individuals with intermediate antibody levels.
3.  **Dissect genetic effects** on both the binary outcome (serostatus) and the quantitative trait (IgG levels).

## 2. Phenotype Definitions

For each of the **24 core antigens**, the script generates four distinct phenotype definitions.

### A. IgG-level Phenotypes (Quantitative Traits)

These phenotypes test for genetic variants associated with the continuous antibody level (a γ-effect).

1.  **`{antigen}_IgG_raw`**
    *   **Definition**: The natural log of the raw Mean Fluorescence Intensity (MFI) value for all individuals.
    *   **Purpose**: Serves as the baseline quantitative trait. A GWAS on this phenotype identifies variants affecting the antibody level across the entire population, without distinguishing between seropositive and seronegative states.

2.  **`{antigen}_IgG_wgt`** (Weighted IgG)
    *   **Definition**: The same `log(MFI)` value as above, but analyzed in conjunction with a **weights file**.
    *   **Weight**: The posterior probability of being seropositive (`p_soft`) from the mixture model.
    *   **Purpose**: To focus the quantitative trait analysis on individuals who are likely seropositive. By down-weighting the influence of seronegative individuals (whose `p_soft` is close to 0), this analysis preferentially detects variants that modulate antibody levels *within the seropositive group*.

### B. Serostatus Phenotypes (Binary Traits)

These phenotypes test for genetic variants associated with the risk of being seropositive (a β-effect).

3.  **`{antigen}_sero_hard`**
    *   **Definition**: A binary (0/1) phenotype based on the official hard cutoff threshold for that antigen. This is the traditional case/control definition.
    *   **Purpose**: To replicate the standard analysis method used in Butler-Laporte et al., providing a direct comparison. It is analyzed using logistic regression in Quickdraws.
    *   **Associated Weight**: This phenotype is paired with a weights file where the weight is `2 * |p_soft - 0.5|`. This down-weights individuals with ambiguous serostatus (where `p_soft` is near 0.5), effectively focusing the logistic regression on the most confidently classified cases and controls.

4.  **`{antigen}_sero_soft`**
    *   **Definition**: The posterior probability of being seropositive (`p_soft`) itself, treated as a continuous quantitative trait ranging from 0 to 1.
    *   **Purpose**: This is our most powerful method for detecting serostatus associations. It avoids information loss by not forcing individuals into a hard 0/1 classification. For GWAS-scale effects, analyzing this quantitative trait with a linear model yields p-values that are nearly identical to a more complex quasi-binomial logistic regression, making it both powerful and practical.

## 3. Novel Covariate: `PC1_assay_noise`

A key innovation in this script is the generation of a novel covariate to account for shared technical noise.

*   **Rationale**: We hypothesized that there might be systematic, non-biological factors that affect an individual's baseline MFI readings across all antigens (e.g., background antibody stickiness, assay batch effects). Controlling for this "assay noise" should increase the precision of all association tests.
*   **Methodology**:
    1.  For each antigen, we first identify the set of **"confidently seronegative"** individuals (defined as having a posterior probability `p_soft < 0.1`).
    2.  We create a matrix containing the `log(MFI)` values for these seronegative individuals only. For a given person, if they are seropositive for an antigen, their value for that antigen is set to `NA`.
    3.  This matrix is then imputed (filling `NA`s with the column mean for that antigen) to create a complete dataset.
    4.  Principal Component Analysis (PCA) is run on this matrix.
*   **Result**: The first principal component, **`PC1_assay_noise`**, represents the major axis of shared variance in IgG levels *among the negative controls*. It is included as an additional covariate in all GWAS models to absorb and control for this shared background noise.

## 4. Within-Antigen Correlation Analysis

As a final quality control step, the script calculates and plots the raw Pearson correlations between these different phenotype definitions for each antigen. This provides a direct, unweighted comparison of the different approaches.

*   **`cor_igg_vs_soft`**: The standard Pearson correlation between the raw `log(MFI)` and the soft probability (`p_soft`).
*   **`cor_hard_vs_soft`**: The correlation between the hard 0/1 classification and the soft probability. This shows how well the continuous probability aligns with the traditional binary classification.

These correlations provide valuable insight into how the different analytical approaches relate to one another and help validate the mixture model's output. 