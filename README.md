# UKB Serology GWAS: A Power-Enhanced Workflow

This project implements a comprehensive workflow for performing a power-enhanced genome-wide association study (GWAS) on UK Biobank serology data. It moves beyond traditional hard-cutoff serostatus calls to a more nuanced, statistically powerful mixture-model-based approach.

## Core Objective

The central goal is to increase statistical power and discover novel genetic associations with pathogen exposure and antibody response. This is achieved by:

1.  **Probabilistic Phenotyping**: Using two-component mixture models to convert continuous antibody titre data (MFI) into a posterior probability of seropositivity (`p_soft`). This retains information lost in binary classification.
2.  **Bayesian Pathogen-Level Calls**: Combining evidence from multiple antigens for a single pathogen (e.g., EBV, CMV) into a unified, more robust serostatus call using Bayesian logistic regression.
3.  **Comprehensive Analysis Strategy**: Testing for genetic effects on both binary serostatus and quantitative antibody levels, using a variety of weighted and unweighted models.
4.  **Direct Paper Comparison**: Implementing exact Butler-Laporte et al. 2023 methodology for direct comparison with published results.

## Project Workflow

The analysis is structured as a clear, iterative pipeline. Each step is documented in a corresponding markdown file in this directory.

1.  **[Data Extraction](./01_DATA_EXTRACTION.md)**
    *   **Goal**: Extract raw serology and genetic data for the UK Biobank cohort.
    *   **Script**: `ukb_potential_power/02_creating_pheno/extract_data.sh`

2.  **[Single-Antigen Phenotype Generation](./02_PHENOTYPE_GENERATION.md)**
    *   **Goal**: Process raw MFI data to generate a suite of "soft" (probabilistic) and "hard" (binary) phenotypes for each individual antigen. This is the core of the power-enhancement strategy.
    *   **Script**: `ukb_potential_power/02_creating_pheno/02_create_serology_phenotypes.R`

3.  **[Pathogen-Level Phenotype Generation](./03_PATHOGEN_LEVEL_PHENOTYPES.md)**
    *   **Goal**: Combine single-antigen data for multi-antigen pathogens (like EBV, CMV, HP) into a single, more reliable serostatus probability using multivariate mixture models.
    *   **Script**: `ukb_potential_power/02_creating_pheno/04_pathogen_level_phenotypes.R`

4.  **[Bayesian Pathogen-Level Regression](./06_PATHOGEN_BAYESIAN_REGRESSION.md)**
    *   **Goal**: Implement robust Bayesian logistic regression to combine per-antigen probabilities into pathogen-level serostatus, with comprehensive uncertainty quantification and direct comparison with Butler-Laporte et al. 2023.
    *   **Script**: `ukb_potential_power/02_creating_pheno/06_pathogen_level_phenotypes_bayesReg.R`

5.  **[GWAS Analysis with Quickdraws](./04_GWAS_ANALYSIS.md)**
    *   **Goal**: Execute a comprehensive GWAS using the generated phenotypes and a manifest-driven approach with the Quickdraws software.
    *   **Scripts**: `03_quickdraw_gwas/` contains the run plan and supporting scripts.

## Key Outputs

-   **Phenotype & Covariate Files**: Located in `quickdraws_input/`, these files are formatted for direct use with Quickdraws and include master phenotype tables, covariate files, and per-analysis weight files.
-   **GWAS Summary Statistics**: The ultimate output of the pipeline, ready for downstream analysis and interpretation.
-   **Diagnostic Plots & Reports**: Various intermediate outputs used for quality control, such as mixture model fit diagnostics and seroprevalence estimates.
-   **Comparison Files**: Direct comparison with Butler-Laporte et al. 2023 methodology for validation and benchmarking.

## Usage

Follow the numbered markdown documents in order to understand the project's logic and execute the analysis scripts. Each document details the purpose, methods, and implementation of its corresponding step. 