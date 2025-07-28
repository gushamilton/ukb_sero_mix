






#  now to run step 2

dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bed" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bim" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.fam" \
  -iin="QuickdrawsGWAS/phenotypes_quantitative.tsv" \
  -iin="QuickdrawsGWAS/covariates_master.tsv" \
  -iin="QuikdrawsGWAS/unrel_FID_IID.txt" \
  -iin="QuickdrawsGWAS/step1_output_quant.tar" \
  -icmd="docker load -i quickdraws.tar && tar -xvf step1_output_quant.tar && docker run  --security-opt seccomp=unconfined --rm -v ./:/mnt/ quickdraws quickdraws-step-2 --out /mnt/step2_quant --bed /mnt/genotype_antibody_only --out_step1 /mnt/step1 --covarFile /mnt/covariates_master.tsv --unrel_sample_list /mnt/unrel_FID_IID.txt && rm -r step1*" \
  --instance-type mem1_ssd1_v2_x8 \
  --name qd_step2_quant \
  --priority high \
  -yc


BINARY

dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bed" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bim" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.fam" \
  -iin="QuickdrawsGWAS/phenotypes_binary.tsv" \
  -iin="QuickdrawsGWAS/covariates_master.tsv" \
  -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
  -iin="QuickdrawsGWAS/step1_output_binary.tar" \
  -icmd="docker load -i quickdraws.tar && tar -xvf step1_output_binary.tar && docker run  --security-opt seccomp=unconfined --rm -v ./:/mnt/ quickdraws quickdraws-step-2 --out /mnt/step2_binary --bed /mnt/genotype_antibody_only --out_step1 /mnt/binary_step1 --covarFile /mnt/covariates_master.tsv --unrel_sample_list /mnt/unrel_FID_IID.txt  --binary --firth && rm -rf step1_output_binary.tar" \
  --instance-type mem3_ssd1_v2_x8 \
  --name qd_step2_binary \
  --priority high \
  -y