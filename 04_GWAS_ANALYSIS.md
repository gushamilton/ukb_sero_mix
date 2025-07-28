# Step 4: GWAS Analysis with Quickdraws

This document lays out the exact workflow for running the serology GWAS on the RAP using the **Quickdraws** software. The workflow is designed around the manifest-driven output from the phenotype generation steps.

## Core Design

1.  **Step-1 (Null Model Fitting) is run once** across *all* phenotypes. This is efficient and captures the shared null model (kinship, covariates) for every trait.
2.  **Step-2 (Association Testing) is run separately** for every phenotype/weight combination defined in the `analysis_manifest.tsv`. This is essential because weights are applied only in Step-2, and analysis types (linear vs. logistic) differ.

---

## 1. Analysis Catalogue

The `analysis_manifest.tsv` file defines every analysis to be performed. For each antigen, it specifies multiple analysis strategies designed to answer different biological questions.

### Planned Analysis Grid

| ID | Phenotype Column | Model in QD | Weights File? | Rationale |
| :--- | :--- | :--- | :--- | :--- |
| A | `{ab}_sero_hard` | `logistic` | None | Direct replication of classic binary GWAS. |
| B | `{ab}_IgG_seropos_only`| `linear` | None | IgG levels in hard seropositives only; replicates Butler-Laporte's quantitative scan. |
| C | `{ab}_IgG_raw` | `linear` | `{ab}_IgG_wgt.weights` | IgG levels weighted by `p_soft`; focuses quantitative analysis on likely positives. |
| D | `{ab}_sero_hard` | `logistic` | `{ab}_w_hard.weights` | Hard serostatus weighted by confidence (`2*|p-0.5|`); focuses on confident cases/controls. |
| E | `{ab}_sero_soft` | **`linear`** | None | **Primary probabilistic analysis:** `p_soft` as a quantitative trait. Simplest, fastest, and most powerful approach for serostatus genetics. |

*(Note: A previous version included a methodological cross-check analysis. This has been streamlined to the core set above for clarity).*

### Covariate Set

All analyses will use a master covariate file (`quickdraws_input/covariates_base.tsv`) containing:
-   Age
-   Sex
-   Genetic Principal Components 1-20
-   **`latent_factor_1`** (the shared assay noise covariate)

---

## 2. Quickdraws Input Structure

The `02_...R` and `04_...R` scripts produce the following file structure in `quickdraws_input/`, which is ready for upload to the RAP:

```
quickdraws_input/
 ├─ covariates_base.tsv
 ├─ phenotypes_quantitative.tsv
 ├─ phenotypes_binary.tsv
 ├─ analysis_manifest.tsv
 ├─ *.weights                     # All weight files
 └─ kinship.txt                   # (Pre-computed)
```

---

## 3. Execution on RAP

### 3.1 Data Upload

-   Upload the entire `quickdraws_input/` directory to a project folder on the RAP.
-   Upload the genotype data (e.g., PLINK `*.bed/bim/fam` files).
-   Upload the Quickdraws Docker image (`quickdraws.tar`).

### 3.2 Step-1 (Model Fitting)

Run a single `quickdraws-step-1` command. As this step can be computationally intensive, a GPU-enabled instance is recommended.

```bash
# Example command inside dx run app-swiss-army-knife
cmd="docker load -i quickdraws.tar && \
     docker run --security-opt seccomp=unconfined --gpus all --rm \
       -v ./:/mnt/ quickdraws \
       quickdraws-step-1 \
         --out /mnt/step1_output \
         --bed /mnt/genotypes  \
         --phenoFile /mnt/phenotypes_quantitative.tsv \
         --phenoFile /mnt/phenotypes_binary.tsv \
         --covarFile /mnt/covariates_base.tsv \
         --kinship /mnt/kinship.txt"
```

This produces a single `step1_output` directory containing the fitted null models for all traits.

### 3.3 Step-2 (Association Testing)

This stage is parallelized, with one job per row in `analysis_manifest.tsv`. The following pseudo-script illustrates the logic.

```bash
# For each line in analysis_manifest.tsv...
while IFS=$'\t' read -r phenotype analysis_type weights_file; do

  # Skip header
  [[ $phenotype == "phenotype_name" ]] && continue

  out_prefix="step2_results/${phenotype}"
  cmd="quickdraws-step-2 --out /mnt/${out_prefix} \
       --bed /mnt/genotypes \
       --out_step1 /mnt/step1_output \
       --phenoName ${phenotype} \
       --covarFile /mnt/covariates_base.tsv"

  # Set analysis type
  if [[ $analysis_type == "logistic" ]]; then
    cmd+=" --logistic"
    # Specify which phenotype file to use
    cmd+=" --phenoFile /mnt/phenotypes_binary.tsv"
  else
    cmd+=" --phenoFile /mnt/phenotypes_quantitative.tsv"
  fi

  # Add weights if specified
  if [[ $weights_file != "NA" ]]; then
    cmd+=" --weightsFile /mnt/${weights_file}"
  fi

  # ... wrap in Docker and launch dx run job ...
done < analysis_manifest.tsv
```

This manifest-driven approach ensures every planned analysis is executed systematically, producing a separate summary statistics file for each one. 