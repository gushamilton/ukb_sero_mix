






#  now to run step 2 part 1 and part 2

dx run app-swiss-army-knife \
  --folder="QuickdrawsOutput/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bed" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bim" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.fam" \
  -iin="QuickdrawsGWAS/phenotypes_quantitative1.tsv" \
  -iin="QuickdrawsGWAS/covariates_master.tsv" \
  -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
  -iin="QuickdrawsOutput/step1_output_p1.tar" \
  -icmd="docker load -i quickdraws.tar && tar -xvf step1_output_p1.tar && docker run  --security-opt seccomp=unconfined --rm -v ./:/mnt/ quickdraws quickdraws-step-2 --out /mnt/step2_p1 --bed /mnt/genotype_antibody_only --out_step1 /mnt/step1_p1 --covarFile /mnt/covariates_master.tsv --unrel_sample_list /mnt/unrel_FID_IID.txt && rm -r step1_p1*" \
  --instance-type mem1_ssd1_v2_x8 \
  --name qd_step2_quant \
  --priority high \
  -y

dx run app-swiss-army-knife \
  --folder="QuickdrawsOutput/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bed" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bim" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.fam" \
  -iin="QuickdrawsGWAS/phenotypes_quantitative2.tsv" \
  -iin="QuickdrawsGWAS/covariates_master.tsv" \
  -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
  -iin="QuickdrawsOutput/step1_output_p2.tar" \
  -icmd="docker load -i quickdraws.tar && tar -xvf step1_output_p2.tar && docker run  --security-opt seccomp=unconfined --rm -v ./:/mnt/ quickdraws quickdraws-step-2 --out /mnt/step2_p2 --bed /mnt/genotype_antibody_only --out_step1 /mnt/step1_p2 --covarFile /mnt/covariates_master.tsv --unrel_sample_list /mnt/unrel_FID_IID.txt && rm -r step1_p2*" \
  --instance-type mem1_ssd1_v2_x8 \
  --name qd_step2_quant \
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
  -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
  -iin="QuickdrawsOutputBinary/step1_output.tar" \
  -icmd="docker load -i quickdraws.tar && tar -xvf step1_output.tar && docker run  --security-opt seccomp=unconfined --rm -v ./:/mnt/ quickdraws quickdraws-step-2 --out /mnt/step2 --bed /mnt/genotype_antibody_only --out_step1 /mnt/step1 --covarFile /mnt/covariates_master.tsv --unrel_sample_list /mnt/unrel_FID_IID.txt  --binary --firth && rm -rf step1*" \
  --instance-type mem3_ssd1_v2_x8 \
  --name qd_step2_binary \
  --priority high \
  -y