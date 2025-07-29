dx run app-swiss-army-knife \
  --folder="QuickdrawsOutput/" \
  -iin="QuickdrawsGWAS/quickdraws.tar" \
  -iin="QuickdrawsGWAS/chr22_imputed_antibody.bed" \
  -iin="QuickdrawsGWAS/chr22_imputed_antibody.bim" \
  -iin="QuickdrawsGWAS/chr22_imputed_antibody.fam" \
  -iin="QuickdrawsGWAS/phenotypes_quantitative1.tsv" \
  -iin="QuickdrawsGWAS/phenotypes_quantitative2.tsv" \
  -iin="QuickdrawsGWAS/covariates_master.tsv" \
  -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
  -iin="QuickdrawsOutput/step1_output_p1.tar" \
  -iin="QuickdrawsOutput/step1_output_p2.tar" \
  -iin="QuickdrawsOutput/step2_p1.calibration" \
  -iin="QuickdrawsOutput/step2_p2.calibration" \
    -iin="chr22_weighted.sh" \
  -iin="QuickdrawsGWAS/weights/bkv_vp1_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/cmv_pp150_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/cmv_pp28_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/cmv_pp52_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/ct_mompa_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/ct_mompd_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/ct_tarpf1_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/ct_tarpf2_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/ct_pgp3_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/ebv_ead_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/ebv_ebna1_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/ebv_vca_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/ebv_zebra_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/hhv6_ie1a_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/hhv6_ie1b_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/hhv7_u14_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/hp_caga_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/hp_catalase_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/hp_groel_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/hp_omp_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/hp_urea_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/hp_vaca_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/hsv1_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/hsv2_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/jcv_vp1_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/kshv_lana_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/mcv_vp1_p_soft.weights" \
  -iin="QuickdrawsGWAS/weights/toxo_p22_p_soft.weights" \
  -icmd="bash chr22_weighted.sh" \
  --instance-type=mem1_ssd1_v2_x16 \
  --name="qd_step2_weighted_chr22" \
  --priority=high \
  -y


  # Run weighted step across autosomes 1‑22
for CHR in {1..19}; do
  dx run app-swiss-army-knife \
    --folder="QuickdrawsOutput/" \
    -iin="QuickdrawsGWAS/quickdraws.tar" \
    -iin="QuickdrawsGWAS/phenotypes_quantitative1.tsv" \
    -iin="QuickdrawsGWAS/phenotypes_quantitative2.tsv" \
    -iin="QuickdrawsGWAS/covariates_master.tsv" \
    -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
    -iin="QuickdrawsOutput/step1_output_p1.tar" \
    -iin="QuickdrawsOutput/step1_output_p2.tar" \
    -iin="QuickdrawsOutput/step2_p1.calibration" \
    -iin="QuickdrawsOutput/step2_p2.calibration" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.bed" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.bim" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.fam" \
    -iin="chr${CHR}_weighted.sh" \
    -iin="QuickdrawsGWAS/weights/bkv_vp1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/cmv_pp150_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/cmv_pp28_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/cmv_pp52_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ct_mompa_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ct_mompd_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ct_tarpf1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ct_tarpf2_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ct_pgp3_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ebv_ead_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ebv_ebna1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ebv_vca_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ebv_zebra_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hhv6_ie1a_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hhv6_ie1b_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hhv7_u14_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_caga_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_catalase_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_groel_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_omp_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_urea_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_vaca_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hsv1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hsv2_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/jcv_vp1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/kshv_lana_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/mcv_vp1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/toxo_p22_p_soft.weights" \
    -icmd="bash chr${CHR}_weighted.sh" \
    --instance-type=mem2_ssd1_v2_x8 \
    --name="qd_step2_weighted_chr${CHR}" \
    --priority=high \
    -y
done


#!/usr/bin/env bash
# Rerun qd_step2_weighted only for chromosomes that failed on 2025‑07‑29

FAILED=(1 3 12)

for CHR in "${FAILED[@]}"; do
  dx run app-swiss-army-knife \
    --folder="QuickdrawsOutput/" \
    -iin="QuickdrawsGWAS/quickdraws.tar" \
    -iin="QuickdrawsGWAS/phenotypes_quantitative1.tsv" \
    -iin="QuickdrawsGWAS/phenotypes_quantitative2.tsv" \
    -iin="QuickdrawsGWAS/covariates_master.tsv" \
    -iin="QuickdrawsGWAS/unrel_FID_IID.txt" \
    -iin="QuickdrawsOutput/step1_output_p1.tar" \
    -iin="QuickdrawsOutput/step1_output_p2.tar" \
    -iin="QuickdrawsOutput/step2_p1.calibration" \
    -iin="QuickdrawsOutput/step2_p2.calibration" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.bed" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.bim" \
    -iin="QuickdrawsGWAS/chr${CHR}_imputed_antibody.fam" \
    -iin="chr${CHR}_weighted.sh" \
    -iin="QuickdrawsGWAS/weights/bkv_vp1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/cmv_pp150_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/cmv_pp28_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/cmv_pp52_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ct_mompa_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ct_mompd_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ct_tarpf1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ct_tarpf2_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ct_pgp3_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ebv_ead_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ebv_ebna1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ebv_vca_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/ebv_zebra_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hhv6_ie1a_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hhv6_ie1b_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hhv7_u14_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_caga_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_catalase_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_groel_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_omp_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_urea_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hp_vaca_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hsv1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/hsv2_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/jcv_vp1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/kshv_lana_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/mcv_vp1_p_soft.weights" \
    -iin="QuickdrawsGWAS/weights/toxo_p22_p_soft.weights" \
    -icmd="bash chr${CHR}_weighted.sh" \
    --instance-type=mem2_ssd1_v2_x16 \
    --name="qd_step2_weighted_chr${CHR}_retry" \
    --priority=high \
    -y
done
