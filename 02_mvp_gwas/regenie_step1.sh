#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail
set -o xtrace

# Submit with: bsub -env all < 02_mvp_gwas/regenie_step1.sh
# Export the required variables described below before submission.

#BSUB -J regenie_step1
#BSUB -G mvp001
#BSUB -M 10000
#BSUB -a "multithread(8)"
#BSUB -q regenie_step1
#BSUB -o logs/regenie_step1.%J.stdout
#BSUB -e logs/regenie_step1.%J.stderr

: "${REGENIE_STEP1_BIN:?Set REGENIE_STEP1_BIN to the REGENIE 3.1.3 executable}"
: "${BED_PREFIX:?Set BED_PREFIX to the PLINK 1 binary-file prefix}"
: "${EXTRACT_FILE:?Set EXTRACT_FILE to the one-variant-ID-per-line file}"
: "${KEEP_FILE:?Set KEEP_FILE to the participant keep file}"
: "${PHENO_FILE:?Set PHENO_FILE to the REGENIE phenotype file}"
: "${COVAR_FILE:?Set COVAR_FILE to the REGENIE covariate file}"
: "${RESULTS_DIR:?Set RESULTS_DIR to the output directory}"

ANCESTRY=${ANCESTRY:-AMR}
OUT_PREFIX="${RESULTS_DIR}/RA.${ANCESTRY}.step1"

mkdir -p "${RESULTS_DIR}" logs

# Step 1 fits the null logistic mixed model using common genotyped variants.
"${REGENIE_STEP1_BIN}" \
  --step 1 \
  --bed "${BED_PREFIX}" \
  --extract "${EXTRACT_FILE}" \
  --phenoFile "${PHENO_FILE}" \
  --phenoCol binary_ppv_90_rm \
  --covarFile "${COVAR_FILE}" \
  --covarColList age,sex,pc1,pc2,pc3,pc4,pc5 \
  --keep "${KEEP_FILE}" \
  --bsize 1000 \
  --bt \
  --lowmem \
  --loocv \
  --threads 4 \
  --verbose \
  --out "${OUT_PREFIX}"
