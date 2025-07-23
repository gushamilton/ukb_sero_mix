
#  now to run step 2

dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="docker/quickdraws.tar" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bed" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bim" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.fam" \
  -iin="QuickdrawsGWAS/phenotypes_quantitative.tsv" \
  -iin="QuickdrawsGWAS/covariates_no_lf.txt" \
  -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
  -iin="QuickdrawsGWAS/step1_output.tar" \
  -icmd="docker load -i quickdraws.tar && tar -xvf step1_output.tar && docker run  --security-opt seccomp=unconfined --rm -v ./:/mnt/ quickdraws quickdraws-step-2 --out /mnt/step2_22 --bed /mnt/genotype_antibody_only --out_step1 /mnt/step1 --covarFile /mnt/covariates_no_lf.txt --unrel_sample_list /mnt/unrel_FID_IID.txt && rm -r step1*" \
  --instance-type mem1_ssd1_v2_x8 \
  --name qd_step2 \
  --priority high \
  -y