#!/usr/bin/env bash
# -------------------------------------------------
#  Step 2 for *unweighted quantitative* traits
#  44 jobs: 22 chromosomes × 2 phenotype batches
# -------------------------------------------------
set -euo pipefail

# ---------- batch p1 ----------------------------------------------------------
for CHR in {1..20}; do
  dx run app-swiss-army-knife \
    --folder="QuickdrawsOutput/" \
    -iin="QuickdrawsGWAS/quickdraws.tar" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.bed" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.bim" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.fam" \
    -iin="QuickdrawsGWAS/phenotypes_quantitative1.tsv" \
    -iin="QuickdrawsGWAS/covariates_master.tsv" \
    -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
    -iin="QuickdrawsOutput/step1_output_p1.tar" \
    -iin="QuickdrawsOutput/step2_p1.calibration" \
    -icmd="docker load -i quickdraws.tar && \
           tar -xf step1_output_p1.tar && \
           docker run --security-opt seccomp=unconfined --rm -v \$PWD:/mnt quickdraws \
             quickdraws-step-2 \
               --out /mnt/step2_p1_chr${CHR}_impute \
               --bed /mnt/chr${CHR}_imputed_antibody \
               --out_step1 /mnt/step1_p1 \
               --covarFile /mnt/covariates_master.tsv \
               --unrel_sample_list /mnt/unrel_FID_IID.txt \
               --calibrationFile /mnt/step2_p1.calibration && \
           rm -rf step1_p1*" \
    --instance-type=mem3_ssd1_v2_x16 \
    --name="qd_step2_quant_p1_chr${CHR}" \
    --priority=high \
    -y
    sleep 1
done

# ---------- batch p2 ----------------------------------------------------------
for CHR in {1..22}; do
  dx run app-swiss-army-knife \
    --folder="QuickdrawsOutput/" \
    -iin="QuickdrawsGWAS/quickdraws.tar" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.bed" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.bim" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.fam" \
    -iin="QuickdrawsGWAS/phenotypes_quantitative2.tsv" \
    -iin="QuickdrawsGWAS/covariates_master.tsv" \
    -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
    -iin="QuickdrawsOutput/step1_output_p2.tar" \
    -iin="QuickdrawsOutput/step2_p2.calibration" \
    -icmd="docker load -i quickdraws.tar && \
           tar -xf step1_output_p2.tar && \
           docker run --security-opt seccomp=unconfined --rm -v \$PWD:/mnt quickdraws \
             quickdraws-step-2 \
               --out /mnt/step2_p2_chr${CHR}_impute \
               --bed /mnt/chr${CHR}_imputed_antibody \
               --out_step1 /mnt/step1_p2 \
               --covarFile /mnt/covariates_master.tsv \
               --unrel_sample_list /mnt/unrel_FID_IID.txt \
               --calibrationFile /mnt/step2_p2.calibration && \
           rm -rf step1_p2*" \
    --instance-type=mem2_ssd1_v2_x16 \
    --name="qd_step2_quant_p2_chr${CHR}" \
    --priority=high \
    -y
sleep 1 
done
