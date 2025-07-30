#!/usr/bin/env bash
# ---------------------------------------------------------------------------
#  Quickdraws Step 2 – run every IgG_raw trait individually (Parameterized)
#  Accepts chromosome and trait group as arguments.
# ---------------------------------------------------------------------------
set -euo pipefail

# --- Script Parameterization ---
if [[ "$#" -ne 2 ]]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: $0 <chromosome_number> <trait_group>"
    echo "  - chromosome_number: 1-22"
    echo "  - trait_group: p1a, p1b, p2a, or p2b"
    echo "Example: $0 1 p1a"
    exit 1
fi

CHR=$1
GROUP=$2

echo "--- Starting Quickdraws Step 2 ---"
echo "Chromosome: ${CHR}"
echo "Trait Group: ${GROUP}"
echo "------------------------------------"

# 0.  Environment set-up
docker load -i quickdraws.tar
tar -xf step1_output_p1.tar    # gives step1_p1_* outputs
tar -xf step1_output_p2.tar    # gives step1_p2_* outputs

# 1.  Lists of traits you really want
# MODIFIED: Split p1 and p2 into smaller groups for better parallelization
TRAITS_P1A=( hsv1_IgG_raw hsv2_IgG_raw ebv_vca_IgG_raw ebv_ebna1_IgG_raw
             ebv_zebra_IgG_raw ebv_ead_IgG_raw cmv_pp150_IgG_raw )

TRAITS_P1B=( cmv_pp52_IgG_raw cmv_pp28_IgG_raw hhv6_ie1a_IgG_raw hhv6_ie1b_IgG_raw
             hhv7_u14_IgG_raw bkv_vp1_IgG_raw jcv_vp1_IgG_raw )

TRAITS_P2A=( mcv_vp1_IgG_raw ct_mompa_IgG_raw ct_mompd_IgG_raw ct_tarpf1_IgG_raw
             ct_tarpf2_IgG_raw ct_pgp3_IgG_raw hp_vaca_IgG_raw )

TRAITS_P2B=( hp_omp_IgG_raw hp_groel_IgG_raw hp_catalase_IgG_raw hp_urea_IgG_raw
             toxo_p22_IgG_raw kshv_lana_IgG_raw )


# 2.  Common inputs
# MODIFIED: Corrected filenames to remove '_antibody' suffix.
BGEN=/mnt/chr${CHR}_imputed.bgen
# The BGI file is used by the environment but not passed as an argument to quickdraws
SAMPLE=/mnt/chr${CHR}_imputed.sample
COVAR=/mnt/covariates_master.tsv
UNREL=/mnt/unrel_FID_IID.txt

##############################################################
# helper: build a mini-prefix for <trait> and launch step-2  #
##############################################################
run_trait () {
    local TRAIT=$1
    local CURRENT_GROUP=$2 # This will be 'p1a', 'p1b', 'p2a', or 'p2b'

    # Determine the base batch ('p1' or 'p2') for finding step1/calibration files
    local BASE_BATCH
    if [[ "$CURRENT_GROUP" == "p1a" || "$CURRENT_GROUP" == "p1b" ]]; then
        BASE_BATCH="p1"
    elif [[ "$CURRENT_GROUP" == "p2a" || "$CURRENT_GROUP" == "p2b" ]]; then
        BASE_BATCH="p2"
    fi

    local STEP1=step1_${BASE_BATCH}
    local PREFIX=${TRAIT}

    local IDX
    IDX=$(awk -v t="$TRAIT" '$1==t{print NR;exit}' ${STEP1}_pred.list)

    if [[ ! -f ${PREFIX}.traits ]]; then
        awk -v col="$TRAIT" '
            BEGIN{FS=OFS="\t"}
            NR==1{for(i=1;i<=NF;i++) if($i==col) c=i; print $1,$2,col; next}
            {print $1,$2,$c}
        ' ${STEP1}.traits > ${PREFIX}.traits
    fi

    ln -sf ${STEP1}.covar_effects  ${PREFIX}.covar_effects
    ln -sf ${STEP1}.neff           ${PREFIX}.neff
    # CORRECTED: Reverted LOCO filename to use a static '_1.loco' suffix,
    # as expected by the quickdraws python script.
    ln -sf "${STEP1}_${IDX}.loco"  "${PREFIX}_1.loco"

    local WFILE=${TRAIT%%_IgG_raw}_p_soft.weights
    local CAL_FILE="/mnt/step2_${BASE_BATCH}.calibration"

    docker run --security-opt seccomp=unconfined --rm -v "$PWD":/mnt quickdraws \
      quickdraws-step-2 \
        --out /mnt/step2_chr${CHR}_${TRAIT} \
        --bgen ${BGEN} \
        --sample ${SAMPLE} \
        --out_step1 /mnt/${PREFIX} \
        --covarFile ${COVAR} \
        --calibrationFile ${CAL_FILE} \
        --unrel_sample_list ${UNREL} \
        --sample_weights /mnt/${WFILE}
}

# 3.  Loop over traits based on the GROUP argument
if [[ "$GROUP" == "p1a" ]]; then
    for T in "${TRAITS_P1A[@]}"; do run_trait "$T" "p1a"; done
elif [[ "$GROUP" == "p1b" ]]; then
    for T in "${TRAITS_P1B[@]}"; do run_trait "$T" "p1b"; done
elif [[ "$GROUP" == "p2a" ]]; then
    for T in "${TRAITS_P2A[@]}"; do run_trait "$T" "p2a"; done
elif [[ "$GROUP" == "p2b" ]]; then
    for T in "${TRAITS_P2B[@]}"; do run_trait "$T" "p2b"; done
else
    echo "Error: Invalid trait group '${GROUP}'. Must be 'p1a', 'p1b', 'p2a', or 'p2b'."
    exit 1
fi

# 4.  (Optional) tidy-up
rm -f *_1.loco *.traits *.covar_effects *.neff
echo "--- Quickdraws Step 2 for Chr ${CHR}, Group ${GROUP} finished successfully. ---"
