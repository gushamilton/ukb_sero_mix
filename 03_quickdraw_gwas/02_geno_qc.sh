remove bad variants



dx run app-swiss-army-knife \
  -iin="Bulk/Genotype Results/Genotype calls/ukb_snp_qc.txt" \
  --folder="QuickdrawsGWAS/" \
  -icmd="cat ukb_snp_qc.txt | awk '{ print \$1, \$159; }' > SNPlist.unfiltered.txt && \
         cat SNPlist.unfiltered.txt | sed '1d' | awk '{ if (\$2 == 1) print \$1; }' > SNPlist.filtered.QC.txt" \
  --instance-type mem1_ssd1_v2_x2 \
  --name qc_snp \
  --priority normal \
  -y

dx run app-swiss-army-knife \
  --folder="/QuickdrawsGWAS/" \
  -iin="Bulk/Genotype Results/Genotype calls/ukb_sqc_v2.txt" \
  -iin="Bulk/Genotype Results/Genotype calls/ukb22418_c1_b0_v2.fam" \
  -icmd="
    cat ukb_sqc_v2.txt | cut -d ' ' -f 24 > ukb_white_british.txt && \
    cat ukb22418_c1_b0_v2.fam | cut -d ' ' -f 1 > samples.tmp.txt && \
    paste samples.tmp.txt ukb_white_british.txt > INDlist.unfiltered.txt && \
    cat INDlist.unfiltered.txt | sed '1d' | awk '{ if (\$2 == 1) print \$1, \$1; }' > INDlist.filtered.QC.txt && \
    rm samples.tmp.txt && rm ukb_white_british.txt && rm INDlist.unfiltered.txt" \
  --instance-type mem1_ssd1_v2_x2 \
  --name subset_genotype \
  --priority normal \
  -y


  dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="QuickdrawsGWAS/genotype_500k.bed" \
  -iin="QuickdrawsGWAS/genotype_500k.bim" \
  -iin="QuickdrawsGWAS/genotype_500k.fam" \
  -iin="QuickdrawsGWAS/INDlist.filtered.QC.txt" \
  -iin="QuickdrawsGWAS/SNPlist.filtered.QC.txt" \
  -icmd="
    plink2 --bfile genotype_500k --keep INDlist.filtered.QC.txt --extract SNPlist.filtered.QC.txt --make-bed --out genotype_400k && \
    rm genotype_500k.bed && rm genotype_500k.bim && rm genotype_500k.fam && \
    plink2 --bfile genotype_400k --maf 0.01 --make-bed --out genotype_400k_common && \
    rm genotype_400k.bed && rm genotype_400k.bim && rm genotype_400k.fam" \
  --instance-type mem1_ssd1_v2_x8 \
  --name genotype_maf_filter \
  --priority normal \
  -y

dx upload make_unrel_homogenous.py --path QuickdrawsGWAS/
dx run app-dxjupyterlab \
  --folder="QuickdrawsGWAS/" \
  -iin="Bulk/Genotype\ Results/Genotype\ calls/ukb_rel.dat" \
  -iin="QuickdrawsGWAS/INDlist.filtered.QC.txt" \
  -iin="QuickdrawsGWAS/make_unrel_homogenous.py" \
  -icmd="
    python make_unrel_homogenous.py INDlist.filtered.QC.txt ukb_rel.dat unrel_FID_IID.txt && \
    cp ukb_rel.dat kinship.txt" \
  --instance-type mem1_ssd1_v2_x2 \
  --name make_unrel_homo_file \
  --priority normal \
  -y

#MAKE THE FINAL FILE WITH ONLY THE INCLUDED SAMPLES

  dx run app-swiss-army-knife \
  --folder="QuickdrawsGWAS/" \
  -iin="QuickdrawsGWAS/genotype_400k_common.bed" \
  -iin="QuickdrawsGWAS/genotype_400k_common.bim" \
  -iin="QuickdrawsGWAS/genotype_400k_common.fam" \
  -iin="QuickdrawsGWAS/included_samples.txt" \
  -icmd="plink2 --bfile genotype_400k_common \
                --keep included_samples.txt \
                --make-bed --out genotype_antibody_only" \
  --instance-type mem1_ssd1_v2_x8 \
  --name genotype_filtered \
  --priority normal \
  -y