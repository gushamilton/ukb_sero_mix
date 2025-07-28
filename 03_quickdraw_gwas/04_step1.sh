# RUN for quantitative and binary


dx run app-swiss-army-knife \
  --folder="QuickdrawsOutput/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bed" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bim" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.fam" \
  -iin="QuickdrawsGWAS/phenotypes_quantitative.tsv" \
  -iin="QuickdrawsGWAS/covariates_master.tsv" \
  -iin="QuickdrawsGWAS/kinship.txt" \
  -icmd="docker load -i quickdraws.tar && docker run  --security-opt seccomp=unconfined --gpus all --rm --shm-size=16g -v ./:/mnt/ quickdraws quickdraws-step-1 --out /mnt/step1 --bed /mnt/genotype_antibody_only --lowmem --phenoFile /mnt/phenotypes_quantitative.tsv --covarFile /mnt/covariates_master.tsv --kinship /mnt/kinship.txt --chunksize 1024 --rhe_random_vectors 100 && tar -cvf step1_output_quant.tar ./step1*" \
  --instance-type mem2_ssd2_gpu1_x8 \
  --name qd_step1_quant \
  --priority high \
  -y



BINARY



dx run app-swiss-army-knife \
  --folder="QuickdrawsOutputBinary/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bed" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bim" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.fam" \
  -iin="QuickdrawsGWAS/phenotypes_binary.tsv" \
  -iin="QuickdrawsGWAS/covariates_master.tsv" \
  -iin="QuickdrawsGWAS/kinship.txt" \
  -icmd="docker load -i quickdraws.tar && docker run  --security-opt seccomp=unconfined --gpus all --rm --shm-size=16g -v ./:/mnt/ quickdraws quickdraws-step-1 --out /mnt/step1 --bed /mnt/genotype_antibody_only --lowmem --phenoFile /mnt/phenotypes_binary.tsv --covarFile /mnt/covariates_master.tsv --kinship /mnt/kinship.txt --binary --h2_grid && tar -cvf step1_output.tar ./step1*" \
  --instance-type mem2_ssd2_gpu1_x8 \
  --name qd_step1_binary \
  --priority high \
  -y


