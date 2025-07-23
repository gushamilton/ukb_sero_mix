# Quickdraws GWAS Workflow

This document lays out **exactly** how we will run Quickdraws on the RAP for the serology GWAS.  The workflow reflects two key design choices:

1. **Step-1 (model fitting) is run **once** across *all* phenotypes.**  This captures the shared null model for every trait.
2. **Step-2 (test-statistics) is run *separately* for every phenotype/weight combination.**  This is necessary because weights are applied *only* in step-2.

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
| cmv_pp150_IgG_raw       | linear        | NA                        |
| cmv_pp150_IgG_raw       | linear        | cmv_pp150_IgG_wgt.weights |
| cmv_pp150_sero_hard     | logistic      | cmv_pp150_sero_hard.weights |
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