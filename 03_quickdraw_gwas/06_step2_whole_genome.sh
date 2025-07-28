BINARY


dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
    -iin="QuickdrawsGWAS/chr22_imputed_antibody.bed" \
    -iin="QuickdrawsGWAS/chr22_imputed_antibody.bim" \
    -iin="QuickdrawsGWAS/chr22_imputed_antibody.fam" \
  -iin="QuickdrawsGWAS/covariates_master.tsv" \
  -iin="QuickdrawsGWAS/step2_binary.calibration" \
  -iin="QuickdrawsGWAS/step1_output_binary.tar" \
  -iin="QuickdrawsGWAS/step2_binary22.firth_null" \
  -icmd="docker load -i quickdraws.tar && tar -xvf step1_output_binary.tar && docker run  --security-opt seccomp=unconfined --rm -v ./:/mnt/ quickdraws quickdraws-step-2 --out /mnt/step2_22_impute --bed /mnt/chr22_imputed_antibody --out_step1 /mnt/binary_step1 --covarFile /mnt/covariates_master.tsv --calibrationFile /mnt/step2_binary.calibration --binary --firth && rm -r binary_step1* step1_output_binary.tar " \
  --instance-type mem2_ssd1_v2_x8 \
  --name qd_step2_imputed \
  --priority high \
  -y