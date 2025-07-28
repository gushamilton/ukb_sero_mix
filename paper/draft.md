


# A probabilistic phenotyping framework boosts power and reveals novel sources of variance in serological GWAS

**Authors:** [Your Name], et al.

**Journal:** Nature Communications (Target)

## Abstract

Genome-wide association studies (GWAS) of antibody responses are critical for understanding the genetic basis of immunity, but their power is often limited by crude phenotyping. Typically, continuous antibody measurements are converted into a binary serostatus using a fixed threshold, a practice that causes misclassification, reducing statistical power and biasing effect size estimates. Here, we present a framework that overcomes this limitation. First, we employ two-component mixture models to transform raw antibody titres into a probabilistic measure of serostatus (`p_soft`). Through extensive simulation mirroring real data, we demonstrate that this probabilistic approach substantially increases statistical power and reduces bias compared to traditional hard-cutoff methods. Second, as an exploratory analysis, we introduce a novel method to dissect sources of shared variance across antigens by applying Principal Component Analysis (PCA) to the antibody titre matrix, stratified by serostatus. This allows us to derive covariates that capture both technical assay noise and shared biological variance in immune response magnitude. Applying this framework to UK Biobank serology data, we provide a robust pipeline for discovering the genetic architecture of the human immune response, validating our primary phenotyping method and exploring novel sources of heritable variance.

## Introduction

The human immune response to pathogens is highly variable between individuals, and a significant portion of this variability is heritable [ref]. Identifying the specific genetic variants that shape this response is a central goal of modern immunogenetics, offering pathways to new vaccines, therapies, and diagnostics. Genome-wide association studies (GWAS) using serological data—the measurement of antibody levels in blood—have emerged as a powerful, high-throughput tool for this purpose [ref].

However, the statistical power of serology-based GWAS has been persistently hampered by a fundamental challenge in phenotype definition. Serological assays produce continuous data (e.g., Mean Fluorescence Intensity, MFI), but for analysis, this data is often dichotomized into a binary "seropositive" or "seronegative" status based on a pre-defined hard threshold. This approach is problematic for two reasons. First, the bimodal distribution of antibody titres in a population often has significant overlap between the seronegative and seropositive components, meaning any single threshold will inevitably misclassify a substantial number of individuals [ref]. This phenotypic misclassification is known to dramatically reduce the statistical power of GWAS and lead to an underestimation of genetic effect sizes [ref]. Second, it discards all the quantitative information within the seropositive group, precluding the discovery of variants that modulate the *magnitude* of an antibody response without affecting susceptibility.

To address these limitations, we propose and validate a new analytical framework with two main components. Our primary advance is the replacement of hard-threshold classification with a probabilistic approach. We fit two-component statistical mixture models directly to the continuous antibody titre distributions. This allows us to calculate, for each individual and each antigen, the posterior probability of being seropositive—a "soft" phenotype (`p_soft`) that retains quantitative uncertainty. We validate through simulation that this approach, by avoiding definitive misclassification, substantially boosts GWAS power and reduces estimation bias.

As a secondary, exploratory advance, we address the shared variance that exists across measurements for different antigens. This variance can be technical (e.g., batch effects) or biological (e.g., an individual's general propensity to mount a strong antibody response). We introduce a data-driven method to capture this shared variance by performing Principal Component Analysis (PCA) directly on the matrix of antibody titres, stratified by serostatus. This yields novel quantitative traits representing shared "background noise" and shared "response magnitude," respectively. We then use these traits in a GWAS to discover loci with pleiotropic effects on the immune system.

This paper first presents the validation of our probabilistic phenotyping method via simulation, then describes the application of the full framework to the UK Biobank serology dataset, demonstrating its utility in discovering novel, robust genetic associations.

## Methods

### Study Cohort and Serological Data

All data was from the UK Biobank, a large-scale prospective cohort study [ref]. We utilized serological data generated using a custom Olink Proteomics panel measuring IgG antibody reactivity against a wide range of pathogen-derived antigens. Raw data were provided as Mean Fluorescence Intensity (MFI) units.

### Probabilistic Phenotype Generation (`p_soft`)

For each antigen, raw MFI values were log-transformed to better approximate normality within the bimodal distribution. We then fitted a two-component skew-t mixture model to the log(MFI) distribution for all individuals using the R package `mixsmsn` [ref]. The skew-t distribution was chosen for its robustness to outliers and its ability to flexibly model the often-asymmetric shapes of the seronegative and seropositive components.

From the fitted model, we calculated the posterior probability for each individual belonging to the component with the higher mean. This value, termed `p_soft` (ranging from 0 to 1), served as our primary probabilistic phenotype for downstream analyses. For comparison, we also generated traditional binary phenotypes based on manufacturer-provided MFI thresholds (`sero_hard_BL`).

### Pathogen-Level Phenotyping

For pathogens represented by multiple antigens (e.g., EBV, CMV, HP), we aggregated evidence to create a single, robust pathogen-level phenotype. We used a Bayesian logistic regression model where the outcome was the hard-call serostatus defined by clinical rules (e.g., 2 of 3 antigens positive for CMV [ref]), and the predictors were the `p_soft` values for the constituent antigens. The fitted probability from this model for each individual was used as the pathogen-level `p_soft`.

### Derivation of Exploratory Latent Factor Phenotypes

To explore shared sources of variance across antigens, we performed three separate Principal Component Analyses (PCA) on the `log(MFI)` matrix:

1.  **`latent_factor_1` (Seronegative Noise Component):** We created a sub-matrix containing `log(MFI)` data only for individuals deemed confidently seronegative for a given antigen (`p_soft` < 0.1). The first principal component (PC1) from this matrix was extracted. This factor is hypothesized to capture systematic, non-biological variance, such as assay background noise.
2.  **`latent_factor_2` (Seropositive Response Component):** We created a sub-matrix for confidently seropositive individuals (`p_soft` > 0.9). PC1 from this matrix was extracted. This factor is hypothesized to capture shared *biological* variance related to the magnitude of the immune response.
3.  **`latent_factor_igg` (Global IgG Component):** We performed PCA on the full `log(MFI)` matrix for all individuals.

These latent factors were treated as novel quantitative traits for exploratory GWAS.

### GWAS and Statistical Analysis

GWAS was performed using a manifest-driven approach with Quickdraws [ref]. We ran linear regression models for quantitative traits (e.g., standardized `p_soft` values, latent factors) and logistic regression for binary traits. All models were adjusted for age, sex, genotyping batch, and the first 20 genetic principal components to control for population stratification.

### Power and Bias Simulation

We conducted a comprehensive simulation study to quantify the performance gains of our probabilistic phenotyping approach. Using parameters derived from fits to real antigen data (for `cmv_pp150`, `hsv1`, and a high-overlap `hsv1_overlap` scenario), we simulated datasets with known genetic effects. We simulated a SNP affecting both susceptibility (`beta`, a logistic effect on serostatus) and antibody magnitude (`gamma`, a linear effect on log(MFI) in seropositives). We then ran mock GWAS on these simulated datasets using three methods: (1) using the true (known) serostatus as a gold standard, (2) using a traditional hard-cutoff method, and (3) using our `p_soft` method. We compared the statistical power (rate of detecting the known effect) and estimation bias across thousands of simulations under a wide grid of `beta` and `gamma` effect sizes.


