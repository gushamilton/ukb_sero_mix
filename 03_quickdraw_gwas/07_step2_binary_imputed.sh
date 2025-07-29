#!/usr/bin/env bash
# -------------------------------------------------
#  Step 2 for *binary* traits, imputed data
#  22 jobs (one per chromosome)
# -------------------------------------------------
set -euo pipefail

for CHR in {22..23}; do
  dx run app-swiss-army-knife \
    --folder="QuickdrawsOutputBinary/" \
    -iin="QuickdrawsGWAS/quickdraws.tar" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.bed" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.bim" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.fam" \
    -iin="QuickdrawsGWAS/phenotypes_binary.tsv" \
    -iin="QuickdrawsGWAS/covariates_master.tsv" \
    -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
    -iin="QuickdrawsOutputBinary/step2.calibration" \
    -iin="QuickdrawsOutputBinary/step1_output.tar" \
    -icmd="docker load -i quickdraws.tar && \
           tar -xf step1_output.tar && \
           docker run --security-opt seccomp=unconfined --rm -v \$PWD:/mnt quickdraws \
             quickdraws-step-2 \
               --out /mnt/step2_bin_chr${CHR}_impute \
               --bed /mnt/chr${CHR}_imputed_antibody \
               --out_step1 /mnt/step1 \
               --covarFile /mnt/covariates_master.tsv \
               --unrel_sample_list /mnt/unrel_FID_IID.txt \
               --calibrationFile /mnt/step2.calibration \
               --binary --firth && \
           rm -rf step1*" \
    --instance-type=mem3_ssd1_v2_x16 \
    --name="qd_step2_bin_chr${CHR}" \
    --priority=high \
    -y
    sleep 1
done
