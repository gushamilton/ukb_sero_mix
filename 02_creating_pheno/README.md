# UKB Serology Phenotype Creation Workflow

This directory contains the scripts for creating comprehensive serology phenotypes and diagnostics for UK Biobank data.

## Workflow

1. **02_serology_setup.sh** - Setup script for environment and dependencies
2. **02_create_serology_phenotypes.R** - Main script that:
   - Processes raw serology data
   - Fits 2-component skew-t mixture models
   - Generates multiple phenotype definitions (hard/soft IgG/serostatus)
   - Creates PC1_assay_noise and PC1_seropositive covariates
   - Saves mixture model fits as RDS file
   - Outputs phenotype and covariate files for GWAS

3. **03_comprehensive_antigen_diagnostics.R** - Diagnostic script that:
   - Loads saved mixture model fits
   - Creates comprehensive 8-panel diagnostic plots per antigen
   - Generates single PDF with all antigens for sharing

4. **serology_power_analysis.R** - Additional power analysis and distribution plots

## Output Files

All outputs are saved to the `../results/` directory:
- `phenotypes_master.tsv` - All phenotype definitions
- `covariates_master.tsv` - Covariates including PC1 components
- `mixture_model_fits.rds` - Saved mixture model objects
- `comprehensive_antigen_diagnostics.pdf` - Diagnostic plots for all antigens

## Usage

```bash
# Run the main phenotype generation
Rscript 02_create_serology_phenotypes.R

# Generate comprehensive diagnostics
Rscript 03_comprehensive_antigen_diagnostics.R
```

## Key Features

- **Multiple phenotype definitions**: Hard cutoffs, soft probabilities, weighted analyses
- **Mixture model validation**: Actual fitted densities overlaid on distributions
- **Technical noise control**: PC1_assay_noise and PC1_seropositive covariates
- **Comprehensive diagnostics**: 8-panel plots per antigen for quality assessment 