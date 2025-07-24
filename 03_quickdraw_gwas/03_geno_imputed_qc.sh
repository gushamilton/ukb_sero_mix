dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="Bulk/Imputation/UKB imputation from genotype/ukb22828_c22_b0_v3.bgen" \
  -iin="Bulk/Imputation/UKB imputation from genotype/ukb22828_c22_b0_v3.bgen.bgi" \
  -iin="Bulk/Imputation/UKB imputation from genotype/ukb22828_c22_b0_v3.sample" \
  -iin="QuickdrawsGWAS/keep.ids" \
  -icmd="plink2 \
           --bgen ukb22828_c22_b0_v3.bgen ref-first \
           --sample ukb22828_c22_b0_v3.sample \
           --keep keep.ids \
           --make-bed \
           --out chr22_imputed_antibody" \
  --instance-type mem1_ssd1_v2_x8 \
  --name make_chr22_imputed_subset \
  -y


  for chr in {1..21}; do
  sleep 1
  dx run app-swiss-army-knife \
    --folder="QuickdrawsGWAS/" \
    -iin="Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3.bgen" \
    -iin="Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3.bgen.bgi" \
    -iin="Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3.sample" \
  -iin="QuickdrawsGWAS/keep.ids" \
    -icmd="plink2 \
             --bgen ukb22828_c${chr}_b0_v3.bgen ref-first \
             --sample ukb22828_c${chr}_b0_v3.sample \
           --keep keep.ids \
             --make-bed \
             --out chr${chr}_imputed_antibody" \
    --instance-type mem1_ssd2_v2_x8 \
    --name make_chr${chr}_imputed_subset \
    -y &
done

wait



for chr in {21..22}; do
  sleep 1

  dx run app-swiss-army-knife \
    --folder="QuickdrawsGWAS/" \
    -iin="QuickdrawsGWAS/chr${chr}_imputed_antibody.bed" \
    -iin="QuickdrawsGWAS/chr${chr}_imputed_antibody.bim" \
    -iin="QuickdrawsGWAS/chr${chr}_imputed_antibody.fam" \
    -iin="QuickdrawsGWAS/genotype_antibody_only.fam" \
    -icmd="plink2 \
             --bfile chr${chr}_imputed_antibody \
             --keep genotype_antibody_only.fam \
             --make-bed \
             --out chr${chr}_antibody_white_british" \
    --instance-type mem1_ssd2_v2_x8 \
    --name subset_chr${chr}_white_british \
    -y &
done

wait
echo "All chromosomes subsetted."