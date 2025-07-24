# Step 1: Data Extraction

This step covers the initial extraction of serology and genetic data from the UK Biobank RAP (Research Analysis Platform).

## Objective

To obtain the necessary raw data files for downstream analysis, including:
- Serology MFI (Mean Fluorescence Intensity) values for all antigens.
- Covariate information (age, sex, genetic principal components).
- Genotype data in PLINK format.

## Implementation

The primary script for this stage is `ukb_potential_power/02_creating_pheno/extract_data.sh`. This script uses the `dx extract_assay` tool, which is a powerful utility for exporting phenotypic data from the UKB RAP database into a tabular format.

### Key Files and Fields

-   **`serology_fields.tsv`**: This file, located in `ukb_potential_power/02_creating_pheno/fields/`, contains the specific UK Biobank field IDs for the serology data. This ensures that the correct data columns are extracted.
-   **Covariates**: The script also extracts standard covariates necessary for any GWAS:
    -   `21003`: Age at recruitment
    -   `31`: Sex
    -   `22009`: Genetic principal components (PC1-PC40)

### Execution

The extraction is performed by running the shell script on a DNAnexus worker:

```bash
# Navigate to the correct directory
cd ukb_potential_power/02_creating_pheno/

# Run the extraction script
./extract_data.sh
```

## Outputs

This script produces several key text files that serve as inputs for the next stage of the pipeline:

-   **`serology_export_title.tsv`**: The main phenotype file containing the raw MFI values for all individuals and antigens, with descriptive column headers.
-   **`age_sex.txt`**: A file containing participant IDs, age, and sex.
-   **`covariates_clean.txt`**: A formatted file containing the genetic principal components for each participant.

These files are the foundational data upon which all subsequent phenotype generation and analysis will be built. 