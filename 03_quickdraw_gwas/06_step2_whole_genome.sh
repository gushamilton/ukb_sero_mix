BINARY


dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
    -iin="QuickdrawsGWAS/chr22_imputed_antibody.bed" \
    -iin="QuickdrawsGWAS/chr22_imputed_antibody.bim" \
    -iin="QuickdrawsGWAS/chr22_imputed_antibody.fam" \
  -iin="QuickdrawsGWAS/phenotypes_binary.tsv" \
  -iin="QuickdrawsGWAS/covariates_master.tsv" \
  -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
  -iin="QuickdrawsOutputBinary/step2.calibration" \
  -iin="QuickdrawsOutputBinary/step1_output.tar" \
  -icmd="docker load -i quickdraws.tar && tar -xvf step1_output.tar && docker run  --security-opt seccomp=unconfined --rm -v ./:/mnt/ quickdraws quickdraws-step-2 --out /mnt/step2_22_impute --bed /mnt/chr22_imputed_antibody --out_step1 /mnt/step1 --covarFile /mnt/covariates_master.tsv --unrel_sample_list /mnt/unrel_FID_IID.txt --calibrationFile /mnt/step2.calibration --binary --firth && rm -r step1* step1_output.tar " \
  --instance-type mem2_ssd1_v2_x8 \
  --name qd_step2_imputed \
  --priority high \
  -y