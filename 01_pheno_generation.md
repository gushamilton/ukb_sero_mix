Serology GWAS: Phenotype Generation and Analysis Plan
1. Overall Strategy
This project aims to conduct a comprehensive and powerful genome-wide association study (GWAS) on serology data from the UK Biobank. The core of the strategy is to move beyond simple binary (positive/negative) classifications based on hard thresholds. Instead, we generate a suite of robust, nuanced phenotypes using probabilistic models.

This approach involves three main stages:

Single-Antigen Phenotype Generation: Modeling the distribution of each antigen individually to derive a foundational probability of seropositivity (p_soft).

Pathogen-Level Phenotype Aggregation: Combining the evidence from multiple antigens for a given pathogen (e.g., EBV, CMV) into a single, refined pathogen-level serostatus probability using a Bayesian model.

Systematic GWAS Execution: Defining a manifest of distinct analysis types to interrogate different biological questions and systematically running them for every antigen, pathogen, and novel latent factor.

2. Single-Antigen Phenotype Generation
This is the foundational step where we process the raw Mean Fluorescence Intensity (MFI) data for each antigen.

2.1 Methodology: Skew-T Mixture Models
To account for the bimodal and often overlapping distributions of serology data, we model the log-MFI distribution for each antigen using a two-component, skew-t mixture model.

The primary output from this model is p_soft, a posterior probability (from 0 to 1) that an individual is seropositive for that specific antigen. This probabilistic measure is central to most of our downstream analyses.

2.2 Novel Phenotypes: Latent Factors
A key innovation of this pipeline is the creation of three novel quantitative phenotypes derived from Principal Component Analysis (PCA) to capture shared variance across antigens. These are GWAS'd directly.

latent_factor_1: Derived from a PCA on the IgG levels of confidently seronegative individuals (p_soft < 0.1). This factor aims to capture shared technical noise or non-specific background binding.

latent_factor_2: Derived from a PCA on the IgG levels of confidently seropositive individuals (p_soft > 0.9). This factor aims to capture shared biological variance among true responders, such as common components of the humoral response.

latent_factor_igg: Derived from a PCA on the IgG levels of all individuals. This represents the main axis of antibody-level variation across the entire dataset.

3. Pathogen-Level Phenotype Aggregation
For pathogens with multiple antigens (e.g., EBV, CMV, HP), we combine the single-antigen evidence to produce a more robust pathogen-level serostatus.

3.1 Rationale
Simple rules (e.g., "positive if ≥2 antigens are positive") ignore the continuous nature of antibody titres, the uncertainty near cutoffs, and the relative importance of each antigen. To overcome this, we replace these hard rules with a principled, probabilistic model.

3.2 Methodology: Bayesian Logistic Regression
We use a Bayesian logistic regression model to learn data-driven weights for each antigen's contribution to the overall pathogen status.

Model Specification (per pathogen):

z_j   ~  Bernoulli(π_j)
logit(π_j) = β₀  + Σ_k β_k · p_{jk}
β_k  ~ Normal(0, 1)
p_{jk}: The single-antigen posterior probability (p_soft) from Step 2, used as a predictor.

z_j: The "noisy" hard status defined by the established UK Biobank/Butler-Laporte rule (e.g., ≥2 antigens positive), used as the target outcome for the model.

This model learns the optimal weight (β_k) for each antigen and outputs a new, refined pathogen-level probability (π_j), which forms the basis for the pathogen-level phenotypes.

4. The Five-Analysis GWAS Plan
Our final script generates phenotypes to support five distinct analysis strategies for each antigen and pathogen. This is all defined in the analysis_manifest.tsv file, which systematically controls the downstream GWAS execution.

4.1 Covariate Set
All analyses use a standard set of covariates from the covariates_base.tsv file:

Age

Sex

Genetic Principal Components 1-20

4.2 Analysis Manifest
ID	Phenotype Column	Model	Weights?	Rationale
BL_SER0	{trait}_sero_hard_BL	linear	None	Butler-Laporte Replication (Binary): An exact replication of the original hard-cutoff binary analysis. Uses a linear model for direct comparability with their fastGWA approach.
BL_IGG_SERPOS	{trait}_IgG_seropos_BL	linear	None	Butler-Laporte Replication (Quantitative): IgG levels in hard seropositives only.
MIX_SER0	{trait}_sero_hard_mix	logistic	None	Optimized Hard Call: A new binary phenotype based on an optimal threshold derived from the mixture/Bayesian model to improve classification accuracy.
P_SOFT	{trait}_psoft_std	linear	None	Primary Discovery (Probabilistic): The standardized p_soft treated as a quantitative trait. This is our most powerful approach for discovering serostatus-associated loci.
IGG_WGT	{trait}_IgG_raw	linear	p_soft	Primary Discovery (Quantitative): Raw IgG levels for all individuals, weighted by p_soft. This focuses the quantitative analysis on likely true positives without discarding ambiguous individuals.
LATENT_FACTOR	latent_factor_1, _2, _igg	linear	None	Novel Quantitative Phenotypes: Direct GWAS of the three latent factors representing shared variance across antigens.

Export to Sheets
5. Key Outputs
The unified pipeline produces two sets of files: one containing all participants for visualization and a second, filtered set in quickdraws_input/gwas_ready_files/ containing only individuals with genotype data, ready for GWAS.

phenotypes_quantitative.tsv & phenotypes_binary.tsv: The phenotype files containing all generated traits.

covariates_base.tsv: The covariate file for use in the GWAS.

analysis_manifest.tsv: The machine-readable file defining every GWAS to be run.

*.weights files: A series of weight files, one for each weighted analysis.

mixture_model_fits.rds: The saved R model objects for each antigen, for diagnostic plotting.
