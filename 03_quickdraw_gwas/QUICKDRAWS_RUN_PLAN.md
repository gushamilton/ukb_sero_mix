# Quickdraws GWAS Workflow

This document lays out **exactly** how we will run Quickdraws on the RAP for the serology GWAS.  The workflow reflects two key design choices:

1. **Step-1 (model fitting) is run **once** across *all* phenotypes.**  This captures the shared null model for every trait.
2. **Step-2 (test-statistics) is run *separately* for every phenotype/weight combination.**  This is necessary because weights are applied *only* in step-2.

---
## 0  Phenotype catalogue & planned analyses  
*(summary extracted from `PHENOTYPE_GENERATION_README.md`)*

For each of the 24 antigens the phenotype-generation script outputs **four base traits** plus two weight columns.  They map to the γ- and β-effects as follows:

| Label                    | Type        | Definition (recap) | Intended GWAS question |
|--------------------------|-------------|--------------------|------------------------|
| `{ab}_IgG_raw`           | quantitative| log-MFI for everyone | Does a variant modulate antibody *titre* in the population? *(γ-effect)* |
| `{ab}_IgG_wgt.weights`   | weight file | weight = `p_soft`   | Focuses IgG analysis on seropositives (down-weights negatives) |
| `{ab}_sero_hard`         | binary      | 0/1 by official cutoff | Classic case/control replication of Butler-Laporte et al. *(β-effect)* |
| `{ab}_w_hard.weights`    | weight file | weight = `2·|p_soft−0.5|` | Emphasises confident cases/controls |
| `{ab}_sero_soft`         | quantitative| posterior `p_soft∈[0,1]` | Soft serostatus, maximises power *(β)* |
| `{ab}_IgG_seropos_only`  | quantitative| IgG levels only in seropositives (NA for seronegatives) | **NEW** — replicates Butler-Laporte IgG analysis restricted to hard-positives |

### Analysis grid we will run

| ID | Phenotype column | Model in QD | Weights? | Rationale |
|----|------------------|-------------|----------|-----------|
| A  | `{ab}_sero_hard` | logistic    | none     | Direct replication of previous binary GWAS |
| B  | `{ab}_IgG_seropos_only` | linear      | none     | IgG analysed **only in hard seropositives** (NA for seronegatives); Butler-Laporte quantitative scan |
| C  | `{ab}_IgG_raw`   | linear      | `{ab}_IgG_wgt.weights` | IgG weighted by p; extends B into probabilistic space |
| D1 | `{ab}_sero_hard` | logistic    | `{ab}_p_soft.weights` | **Validation for E:** see discussion below. |
| D2 | `{ab}_sero_hard` | logistic    | `{ab}_w_hard.weights` | Weight = 2·|p−0.5|; focuses on confident cases & controls |
| E  | `{ab}_sero_soft` | **linear**  | none                   | **Primary probabilistic analysis:** Analyse p directly; simplest & fastest |

We **do not run** an unweighted whole-population IgG analysis.

Row **B** is conceptually identical to the log-MFI analyses performed by Butler-Laporte et&nbsp;al. [Open Forum Infect Dis 2020;7:ofaa450](https://pubmed.ncbi.nlm.nih.gov/33204752/).

### Rationale for Analysis D1 vs E
Analysis E, which performs a simple linear regression on the soft-probability phenotype (`sero_soft`), is our primary and most powerful method for detecting serostatus associations.

Analysis D1, which performs a *weighted logistic regression* on the hard-binary phenotype (`sero_hard`), may seem counter-intuitive. As the weights are simply the `p_soft` values, seronegative individuals (`sero_hard=0`) have their contribution almost entirely removed from the analysis.

This is intentional. The purpose of Analysis D1 is **not** to be a standalone, improved case/control GWAS. Rather, it serves as a crucial **methodological cross-check** for Analysis E. According to score-test theory, the p-values from a weighted logistic regression (D1) should be statistically equivalent to those from a linear regression on the probability itself (E).

If the results from D1 and E are highly concordant, it provides strong evidence that the simpler, faster, and more practical linear model (Analysis E) is a valid and robust approach.

### Logistic vs linear discussion
* `sero_hard`: must be logistic to respect Bernoulli variance.
* `sero_soft`: could be quasi-binomial. However, as confirmed by the D1 vs E comparison, using a linear model is a valid and efficient approximation that avoids potential issues with model convergence or over-dispersion.

### Covariate sets
`covariates_master.tsv` contains **age, sex and the first 20 genetic PCs only**.  No assay-noise PCs are included in this first run, matching the latest RAP plan.  If we wish to add or drop covariates later we can regenerate this file with the script’s flags.

The script now also writes minimal per-analysis phenotype files in `quickdraws_input/step2/`, one `<trait>.phen.tsv` for every row of the manifest, so that each `quickdraws-step-2` invocation can mount only the columns it needs.

---

The file structure produced by `02_create_serology_phenotypes.R` and consumed by Quickdraws looks like this:

```
quickdraws_input/
 ├─ covariates_master.tsv             # Base + PCs + assay-noise PCs
 ├─ phenotypes_master.tsv             # All traits (IgG_raw, sero_soft, sero_hard)
 ├─ analysis_manifest.tsv             # One row per analysis (phenotype, type, weights)
 ├─ <phenotype>_*.weights             # Weight files referenced in manifest (optional)
 └─ kinship.txt                       # KING-robust kinship matrix (pre-computed)
```

---
## 1  Data upload to RAP
* Upload all files inside `quickdraws_input/` to a project folder on RAP, e.g. `QuickdrawsGWAS/`.
* Upload genotype PLINK files (`*.bed/bim/fam`).
* Upload the Quickdraws docker image `quickdraws.tar` (see RAP instructions).

---
## 2  Step-1 – Model fitting (GPU recommended)
We run **one** `quickdraws-step-1` command that references the **combined** phenotype table.

```bash
# inside dx run app-swiss-army-knife
cmd="docker load -i quickdraws.tar && \
     docker run --security-opt seccomp=unconfined --gpus all --rm --shm-size=16g \
       -v ./:/mnt/ quickdraws \
       quickdraws-step-1 \
         --out /mnt/step1 \
         --bed /mnt/genotypes  \
         --phenoFile /mnt/phenotypes_master.tsv \
         --covarFile /mnt/covariates_master.tsv \
         --kinship /mnt/kinship.txt"  # optional: other flags (threads, etc.)
```

The output directory `step1/` contains the fitted null models shared by every downstream test.

---
## 3  Step-2 – Test statistics
### 3.1  Overview
`analysis_manifest.tsv` drives this stage.  Each row has:

* `phenotype_name`  – column in `phenotypes_master.tsv`
* `analysis_type`   – `linear` or `logistic`
* `weights_file`    – filename **or** `NA`

| phenotype_name          | analysis_type | weights_file              |
|-------------------------|---------------|---------------------------|
| cmv_pp150_sero_hard     | logistic      | NA                        |
| cmv_pp150_IgG_seropos_only | linear    | NA                        |
| cmv_pp150_IgG_raw       | linear        | cmv_pp150_IgG_wgt.weights |
| cmv_pp150_sero_hard     | logistic      | cmv_pp150_p_soft.weights  |
| cmv_pp150_sero_hard     | logistic      | cmv_pp150_sero_hard.weights |
| cmv_pp150_sero_soft     | linear        | NA                        |
| ⋯                       | ⋯             | ⋯                         |

Any row where `weights_file = NA` is an *unweighted* analysis; otherwise the weights are passed with `--weightsFile`.

### 3.2  Batch pseudo-script
```bash
# step1 output already untarred into ./step1/
while IFS=$'\t' read -r phenotype analysis_type weights_file; do
  #  skip header
  [[ $phenotype == "phenotype_name" ]] && continue

  out_prefix="step2_${phenotype}"
  cmd="quickdraws-step-2 --out /mnt/${out_prefix} \
       --bed /mnt/genotypes \
       --out_step1 /mnt/step1 \
       --phenoName ${phenotype} \
       --covarFile /mnt/covariates_master.tsv"

  # analysis type
  if [[ $analysis_type == "logistic" ]]; then
    cmd+=" --logistic"
  fi
  # optional weights
  if [[ $weights_file != "NA" ]]; then
    cmd+=" --weightsFile /mnt/${weights_file}"
  fi

  # wrap in docker
  docker_cmd="docker load -i quickdraws.tar && docker run --rm -v ./:/mnt/ quickdraws ${cmd}"

  # launch on RAP (CPU instance; GPUs not required)
  dx run app-swiss-army-knife \
     --folder "QuickdrawsGWAS/" \
     -iin="docker/quickdraws.tar" \
     -iin="QuickdrawsGWAS/genotypes.*" \   # bed/bim/fam wildcard (simplified)
     -iin="QuickdrawsGWAS/${weights_file}" \ # only if not NA
     -iin="QuickdrawsGWAS/step1_output.tar" \
     -icmd="${docker_cmd}" \
     --instance-type mem1_ssd1_v2_x8 \
     --name qd_step2_${phenotype} \
     --priority high -y

done < analysis_manifest.tsv
```
Notes
* `step1_output.tar` is untarred **once** per job.
* Unweighted phenotypes omit `--weightsFile`.
* Logistic flag is added only for binary traits.

---
## 4  Imputed data (optional)
After obtaining `step2_<phenotype>.calibration` for each trait on common variants, reuse the same loop but point `--bgen` / `--sample` to the UKB imputation files **and** supply `--calibrationFile`.

---
## 5  Key design checks
1. **Step-1 once** – quick; saves GPU cost.
2. **Per-phenotype Step-2** – required for weight handling & differing test families.
3. **Manifest-driven** – guarantees every trait/weight combo is tested.
4. **Docker-wrapped** – single environment regardless of RAP node type.

If the above matches expectations, we are ready to script & launch the jobs. 