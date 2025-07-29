#!/usr/bin/env bash
# -----------------------------------------------------------
#  Quickdraws Step 2 – run every IgG_raw trait individually
#  Chromosome 01 — bash + awk only, zero Python
# -----------------------------------------------------------
set -euo pipefail

# 0.  Environment set-up
docker load -i quickdraws.tar
tar -xf step1_output_p1.tar    # gives step1_p1_* outputs
tar -xf step1_output_p2.tar    # gives step1_p2_* outputs

# 1.  Lists of traits you really want
TRAITS_P1=( hsv1_IgG_raw hsv2_IgG_raw ebv_vca_IgG_raw ebv_ebna1_IgG_raw
            ebv_zebra_IgG_raw ebv_ead_IgG_raw cmv_pp150_IgG_raw cmv_pp52_IgG_raw
            cmv_pp28_IgG_raw hhv6_ie1a_IgG_raw hhv6_ie1b_IgG_raw hhv7_u14_IgG_raw
            bkv_vp1_IgG_raw jcv_vp1_IgG_raw )

TRAITS_P2=( mcv_vp1_IgG_raw ct_mompa_IgG_raw ct_mompd_IgG_raw ct_tarpf1_IgG_raw
            ct_tarpf2_IgG_raw ct_pgp3_IgG_raw hp_caga_IgG_raw hp_vaca_IgG_raw
            hp_omp_IgG_raw hp_groel_IgG_raw hp_catalase_IgG_raw hp_urea_IgG_raw
            toxo_p22_IgG_raw kshv_lana_IgG_raw )

# 2.  Common inputs
BED=/mnt/chr5_imputed_antibody
COVAR=/mnt/covariates_master.tsv
UNREL=/mnt/unrel_FID_IID.txt

##############################################################
# helper: build a mini-prefix for <trait> and launch step-2  #
##############################################################
run_trait () {
    local TRAIT=$1          # hsv1_IgG_raw, …
    local BATCH=$2          # p1 | p2
    local STEP1=step1_${BATCH}
    local PREFIX=${TRAIT}   # mini-prefix will be the trait name

    # ---------------- 2a. work out line number in original _pred.list
    local IDX
    IDX=$(awk -v t="$TRAIT" '$1==t{print NR;exit}' ${STEP1}_pred.list)

    # ---------------- 2b. mini-traits file (FID IID Trait)
    if [[ ! -f ${PREFIX}.traits ]]; then
        awk -v col="$TRAIT" '
            BEGIN{FS=OFS="\t"}
            NR==1{
                for(i=1;i<=NF;i++) if($i==col) c=i;
                print $1,$2,col;
                next
            }
            {print $1,$2,$c}
        ' ${STEP1}.traits > ${PREFIX}.traits
    fi

    # ---------------- 2c. symlinks to other Step-1 artefacts
    ln -sf ${STEP1}.covar_effects  ${PREFIX}.covar_effects
    ln -sf ${STEP1}.neff           ${PREFIX}.neff

    # LOCO: rename the right file to  “<trait>_1.loco”
    ln -sf ${STEP1}_${IDX}.loco    ${PREFIX}_1.loco

    # ---------------- 2d. weight file for participation-bias regression
    local WFILE=${TRAIT%%_IgG_raw}_p_soft.weights   # hsv1_p_soft.weights

    # ---------------- 2e. fire up Step 2
    docker run --security-opt seccomp=unconfined --rm -v "$PWD":/mnt quickdraws \
      quickdraws-step-2 \
        --out /mnt/step2_chr5_${TRAIT} \
        --bed ${BED} \
        --out_step1 /mnt/${PREFIX} \
        --covarFile ${COVAR} \
        --unrel_sample_list ${UNREL} \
        --sample_weights /mnt/${WFILE}
}

# 3.  Loop over every requested trait
for T in "${TRAITS_P1[@]}"; do
    run_trait "$T" p1
done

for T in "${TRAITS_P2[@]}"; do
    run_trait "$T" p2
done

# 4.  (Optional) tidy-up
# rm -f *_1.loco *.traits *.covar_effects *.neff
