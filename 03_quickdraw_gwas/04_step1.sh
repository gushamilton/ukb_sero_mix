# RUN for quantitative and binary


dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bed" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bim" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.fam" \
  -iin="QuickdrawsGWAS/phenotypes_quantitative.tsv" \
  -iin="QuickdrawsGWAS/covariates_no_lf.txt" \
  -iin="QuickdrawsGWAS/kinship.txt" \
  -icmd="docker load -i quickdraws.tar && docker run  --security-opt seccomp=unconfined --gpus all --rm --shm-size=16g -v ./:/mnt/ quickdraws quickdraws-step-1 --out /mnt/step1 --bed /mnt/genotype_antibody_only --phenoFile /mnt/phenotypes_quantitative.tsv --covarFile /mnt/covariates_no_lf.txt --kinship /mnt/kinship.txt && tar -cvf step1_output_quant.tar ./step1*" \
  --instance-type mem2_ssd2_gpu1_x8 \
  --name qd_step1_quant \
  --priority high \
  -y


dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bed" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bim" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.fam" \
  -iin="QuickdrawsGWAS/phenotypes_binary.tsv" \
  -iin="QuickdrawsGWAS/covariates_no_lf.txt" \
  -iin="QuickdrawsGWAS/kinship.txt" \
  -icmd="docker load -i quickdraws.tar && docker run  --security-opt seccomp=unconfined --gpus all --rm --shm-size=16g -v ./:/mnt/ quickdraws quickdraws-step-1 --out /mnt/step1 --bed /mnt/genotype_antibody_only --phenoFile /mnt/phenotypes_binary.txt --covarFile /mnt/covariates_no_lf.txt --kinship /mnt/kinship.txt && tar -cvf step1_output_binary.tar ./step1*" \
  --instance-type mem2_ssd2_gpu1_x8 \
  --name qd_step1_quant \
  --priority high \
  -y


  



