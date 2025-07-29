#!/usr/bin/env bash
# -------------------------------------------------
#  Step 2 for *weighted* IgG_raw traits
#  Chromosome 01 — single SAK job
# -------------------------------------------------
set -euo pipefail

# 1. Load container once, unpack Step 1 outputs
docker load -i quickdraws.tar
tar -xf step1_output_p1.tar
tar -xf step1_output_p2.tar

# 2. Trait lists (derived from p1 / p2 batches)
TRAITS_P1=(
  hsv1_IgG_raw hsv2_IgG_raw ebv_vca_IgG_raw ebv_ebna1_IgG_raw
  ebv_zebra_IgG_raw ebv_ead_IgG_raw cmv_pp150_IgG_raw cmv_pp52_IgG_raw
  cmv_pp28_IgG_raw hhv6_ie1a_IgG_raw hhv6_ie1b_IgG_raw hhv7_u14_IgG_raw
  bkv_vp1_IgG_raw jcv_vp1_IgG_raw
)

TRAITS_P2=(
  mcv_vp1_IgG_raw ct_mompa_IgG_raw ct_mompd_IgG_raw ct_tarpf1_IgG_raw
  ct_tarpf2_IgG_raw ct_pgp3_IgG_raw hp_caga_IgG_raw hp_vaca_IgG_raw
  hp_omp_IgG_raw hp_groel_IgG_raw hp_catalase_IgG_raw hp_urea_IgG_raw
  toxo_p22_IgG_raw kshv_lana_IgG_raw
)

BED=/mnt/chr19_imputed_antibody
COVAR=/mnt/covariates_master.tsv
UNREL=/mnt/unrel_FID_IID.txt

run_trait () {
  local TRAIT=$1
  local STEP1_DIR=$2       # /mnt/step1_p{1,2}
  local CAL=$3             # /mnt/step2_p{1,2}.calibration
  local WFILE=${TRAIT%%_IgG_raw}_p_soft.weights

  docker run --security-opt seccomp=unconfined --rm -v "$PWD":/mnt quickdraws \
  python /mnt/quickdraws_step2_patched.py \
    quickdraws-step-2 \
      --out /mnt/step2_chr19_${TRAIT} \
      --bed ${BED} \
      --out_step1 ${STEP1_DIR} \
      --covarFile ${COVAR} \
      --unrel_sample_list ${UNREL} \
      --calibrationFile ${CAL} \
      --sample_weights /mnt/${WFILE}
}

for T in "${TRAITS_P1[@]}"; do
  run_trait "$T" /mnt/step1_p1 /mnt/step2_p1.calibration
done

for T in "${TRAITS_P2[@]}"; do
  run_trait "$T" /mnt/step1_p2 /mnt/step2_p2.calibration
done

rm -rf step1_p{1,2}