# RUN for quantitative and binary

Run the first one: split the quantitative phenotypes into two files

awk 'BEGIN{FS=OFS="\t"} NR==1{n=NF;h=int((n-2)/2)+2;f1="phenotypes_quantitative1.tsv";f2="phenotypes_quantitative2.tsv"; for(i=1;i<=h;i++)printf "%s%s",$i,(i<h?OFS:ORS)>f1; for(i=1;i<=2;i++)printf "%s%s",$i,OFS>f2; for(i=h+1;i<=n;i++)printf "%s%s",$i,(i<n?OFS:ORS)>f2; next} { for(i=1;i<=h;i++)printf "%s%s",$i,(i<h?OFS:ORS)>>f1; for(i=1;i<=2;i++)printf "%s%s",$i,OFS>>f2; for(i=h+1;i<=NF;i++)printf "%s%s",$i,(i<NF?OFS:ORS)>>f2 }' phenotypes_quantitative.tsv


dx run app-swiss-army-knife \
  --folder="QuickdrawsOutput/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bed" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.bim" \
  -iin="QuickdrawsGWAS/genotype_antibody_only.fam" \
  -iin="QuickdrawsGWAS/phenotypes_quantitative1.tsv" \
  -iin="QuickdrawsGWAS/covariates_master.tsv" \
  -iin="QuickdrawsGWAS/kinship.txt" \
  -icmd="docker load -i quickdraws.tar && docker run  --security-opt seccomp=unconfined --gpus all --rm   --shm-size=16g -v ./:/mnt/ quickdraws quickdraws-step-1 --out /mnt/step1_p1 --bed /mnt/genotype_antibody_only --batch_size 64 --lowmem --phenoFile /mnt/phenotypes_quantitative1.tsv --covarFile /mnt/covariates_master.tsv --kinship /mnt/kinship.txt --chunksize 1024 --rhe_random_vectors 150 && tar -cvf step1_output_p1.tar ./step1_p1*" \
  --instance-type mem2_ssd2_gpu1_x8 \
  --name qd_step1_quant \
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
  -iin="QuickdrawsGWAS/kinship.txt" \
  -icmd="docker load -i quickdraws.tar && docker run  --security-opt seccomp=unconfined --gpus all --rm   --shm-size=16g -v ./:/mnt/ quickdraws quickdraws-step-1 --out /mnt/step1_p2 --bed /mnt/genotype_antibody_only --batch_size 64 --lowmem --phenoFile /mnt/phenotypes_quantitative2.tsv --covarFile /mnt/covariates_master.tsv --kinship /mnt/kinship.txt --chunksize 1024 --rhe_random_vectors 150 && tar -cvf step1_output_p2.tar ./step1_p2*" \
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
  -icmd="docker load -i quickdraws.tar && docker run  --security-opt seccomp=unconfined --gpus all --rm --shm-size=16g -v ./:/mnt/ quickdraws quickdraws-step-1 --out /mnt/step1 --bed /mnt/genotype_antibody_only --lowmem --phenoFile /mnt/phenotypes_binary.tsv --covarFile /mnt/covariates_master.tsv --chunksize 2048 --kinship /mnt/kinship.txt --binary --h2_grid && tar -cvf step1_output.tar ./step1*" \
  --instance-type mem2_ssd2_gpu1_x8 \
  --name qd_step1_binary \
  --priority high \
  -y


