# Power-Enhanced Serology GWAS using Mixture Models

This project aims to improve the statistical power of genome-wide association studies (GWAS) on serological data by replacing traditional "hard-cutoff" serostatus calls with a more nuanced mixture model-based approach.

## Project Overview

The standard approach to analyzing serology data in GWAS involves classifying individuals as seropositive or seronegative based on a fixed threshold. However, serological measurements often show bimodal distributions with significant overlap between the two components, leading to misclassification and reduced statistical power.

This project demonstrates that using mixture model-derived probabilities as quantitative phenotypes provides substantial gains in statistical power and reduces estimation bias compared to hard-cutoff methods.

## Project Structure

```
sim_power_gain/
├── sims/                           # GWAS power simulation scripts and outputs
│   ├── gwas_power_simulation.R     # Basic simulation framework
│   ├── gwas_power_sim_comprehensive.R  # Comprehensive power/bias analysis
│   └── comprehensive_power_bias_plots.pdf  # Simulation results
├── ukb_potential_power/            # UK Biobank analysis scripts and data
│   ├── extract_data.sh             # Data extraction script (DNAnexus)
│   ├── serology_power_analysis.R   # Mixture model fitting and power analysis
│   ├── UKB_sero_test.R             # Distribution visualization
│   ├── UKB_sero_simple.R           # Simple analysis scripts
│   ├── fields/                     # Field definitions
│   │   └── serology_fields.tsv     # UKB field IDs for serology data
│   ├── serology_diagnostic_plots.pdf  # Mixture model diagnostic plots
│   └── serology_power_plots.pdf    # Power analysis results
└── draft/                          # Development and exploratory scripts
    ├── draft_*.R                   # Iterative development versions
    ├── debug_*.R                   # Debugging scripts
    ├── em_mixture_model.R          # Custom EM algorithm implementation
    └── test*.R                     # Testing scripts
```

## Workflow

### Stage 1: Theoretical Justification via Simulation

**Goal:** Demonstrate the theoretical power gain of the mixture model approach.

**Implementation:**
1. **Run Simulation:** Execute `sims/gwas_power_sim_comprehensive.R` to compare three methods:
   - **Gold Standard:** Uses true, latent serostatus
   - **Hard Cutoff:** Standard approach using fixed MFI threshold
   - **Mixture Model:** Realistic approach using posterior probabilities

2. **Analyze Results:** Review `sims/comprehensive_power_bias_plots.pdf` to confirm mixture model superiority.

**Intermediate Step 1.5 (Recommended):** Refine simulation parameters using empirical data from Stage 2 to ensure realism.

### Stage 2: UK Biobank Data Analysis & Phenotype Generation

**Goal:** Characterize UKB serology data and generate probabilistic phenotypes for GWAS.

**Implementation:**
1. **Data Extraction:** Use `ukb_potential_power/extract_data.sh` to extract serology data from UKB
2. **Visualize Distributions:** Run `ukb_potential_power/UKB_sero_test.R` to generate distribution plots
3. **Fit Mixture Models:** Execute `ukb_potential_power/serology_power_analysis.R` to fit skew-t mixture models

**Critical Intermediate Step 2.5:** Phenotype QC and Power Check
- Plot histogram of posterior probabilities (should be bimodal)
- Correlate with hard-cutoff phenotype and covariates
- Confirm inflation factor (1/λ²) indicates meaningful power gain

### Stage 3: Power-Enhanced GWAS

**Goal:** Perform GWAS using mixture model-derived probabilities to discover genetic associations.

**Implementation:**
1. **Prepare Inputs:** Assemble genotype data, probabilistic phenotype file, and covariates
2. **Pilot GWAS:** Run on chromosome 22 to validate pipeline
3. **Full GWAS:** Execute complete association analysis

## Selected Antigens

Based on diagnostic plots, power analysis, and [Butler-Laporte et al. (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7641500/), the following antigens are prioritized:

| Short Code | Full Description | Threshold | Rationale |
|------------|------------------|-----------|-----------|
| cmv_pp150 | CMV pp150 | 100 | Strong bimodality, high power gain |
| cmv_pp28 | CMV pp28 | 200 | Clear genetic signal in Butler-Laporte |
| hsv1 | HSV-1 gG | 150 | High prevalence, good separation |
| hsv2 | HSV-2 mgG-2 | 150 | Significant genetic associations |
| ebv_vca | EBV VCA p18 | 250 | Multiple MHC associations |
| ebv_ebna1 | EBV EBNA-1 | 250 | Strong genetic signal |
| ebv_zebra | EBV ZEBRA | 100 | Clear bimodal distribution |
| ebv_ead | EBV EA-D | 100 | Good component separation |
| hp_omp | H. pylori OMP | 170 | Bacterial pathogen representation |
| hp_urea | H. pylori UreA | 130 | Additional H. pylori antigen |
| toxo_sag1 | T. gondii sag1 | 160 | Protozoal pathogen representation |

## Key Findings

1. **Simulation Results:** Mixture models show substantial power gains, especially for detecting effects on antibody levels (gamma effects)
2. **Distribution Analysis:** UKB serology data shows clear bimodality with significant overlap, justifying mixture model approach
3. **Power Loss Quantification:** Hard cutoffs result in 2-5x power loss compared to mixture model approach

## Requirements

- R with packages: `mixsmsn`, `sn`, `dplyr`, `ggplot2`, `patchwork`
- DNAnexus account for UKB data access
- UKB application approval for serology data

## Usage

1. **Data Extraction:**
   ```bash
   cd ukb_potential_power
   ./extract_data.sh
   ```

2. **Run Simulations:**
   ```r
   source("sims/gwas_power_sim_comprehensive.R")
   ```

3. **Analyze UKB Data:**
   ```r
   source("ukb_potential_power/serology_power_analysis.R")
   ```

## References

- Butler-Laporte, G., et al. (2020). Genetic determinants of antibody-mediated immune responses to infectious diseases agents: A genome-wide and HLA association study. *Open Forum Infectious Diseases*, 7(11), ofaa450. [DOI: 10.1093/ofid/ofaa450](https://pmc.ncbi.nlm.nih.gov/articles/PMC7641500/)

## Next Steps

- [ ] Generate phenotype file for selected antigens
- [ ] Implement GWAS pipeline
- [ ] Validate results with external cohorts
- [ ] Extend to additional pathogens

## Contributing

This project is under active development. Please refer to the `draft/` directory for development history and experimental approaches. 