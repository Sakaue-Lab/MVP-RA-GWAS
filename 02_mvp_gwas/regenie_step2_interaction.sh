#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail
set -o xtrace

# Example submission for a 289-row seed file:
# bsub -env all -J "RA_step2_interaction.eur[1-289]" \
#   < 02_mvp_gwas/regenie_step2_interaction.sh

#BSUB -J regenie_step2_interaction
#BSUB -G mvp001
#BSUB -M 15000
#BSUB -a "multithread(8)"
#BSUB -q regenie_step2
#BSUB -o logs/regenie_step2_interaction.%J.%I.stdout
#BSUB -e logs/regenie_step2_interaction.%J.%I.stderr

: "${REGENIE_INTERACTION_BIN:?Set REGENIE_INTERACTION_BIN to the REGENIE 3.4 executable}"
: "${GENOTYPE_DIR:?Set GENOTYPE_DIR to the directory containing chromosome subdirectories}"
: "${SEED_FILE:?Set SEED_FILE to the chromosome/chunk seed file}"
: "${EXTRACT_FILE:?Set EXTRACT_FILE to the targeted variant list}"
: "${KEEP_FILE:?Set KEEP_FILE to the participant keep file}"
: "${PHENO_FILE:?Set PHENO_FILE to the REGENIE phenotype file}"
: "${COVAR_FILE:?Set COVAR_FILE to the REGENIE covariate file}"
: "${PRED_FILE:?Set PRED_FILE to the Step 1 prediction-list file}"
: "${RESULTS_DIR:?Set RESULTS_DIR to the output directory}"
: "${LSB_JOBINDEX:?LSB_JOBINDEX is set automatically by an LSF job array}"

ANCESTRY=${ANCESTRY:-EUR}
PHENO_COLUMN=${PHENO_COLUMN:-binary_ppv_90_rm}
INTERACTION_COVARIATE=${INTERACTION_COVARIATE:-sex}
INTERACTION_LEVEL=${INTERACTION_LEVEL:-0}
INTERACTION_TERM="${INTERACTION_COVARIATE}[${INTERACTION_LEVEL}]"

seed_row=$(awk -v row="${LSB_JOBINDEX}" 'NR == row {print; exit}' "${SEED_FILE}")
if [[ -z "${seed_row}" ]]; then
  echo "No seed-file row found for LSF array index ${LSB_JOBINDEX}." >&2
  exit 1
fi

read -r chromosome chunk extra <<< "${seed_row}"
if [[ -z "${chromosome:-}" || -z "${chunk:-}" || -n "${extra:-}" ]]; then
  echo "Expected exactly two whitespace-delimited fields in seed row: chromosome chunk" >&2
  exit 1
fi

BGEN_FILE="${GENOTYPE_DIR}/${chromosome}/${chunk}.bgen"
SAMPLE_FILE="${GENOTYPE_DIR}/${chromosome}/${chunk}.sample"
OUT_PREFIX="${RESULTS_DIR}/${ANCESTRY}.${chunk}.sex_interaction.step2"

for required_file in \
  "${BGEN_FILE}" "${SAMPLE_FILE}" "${EXTRACT_FILE}" \
  "${KEEP_FILE}" "${PHENO_FILE}" "${COVAR_FILE}" "${PRED_FILE}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Required input does not exist: ${required_file}" >&2
    exit 1
  fi
done

mkdir -p "${RESULTS_DIR}" logs

# Test the SNP-by-sex interaction for targeted variants in one BGEN chunk.
"${REGENIE_INTERACTION_BIN}" \
  --step 2 \
  --bgen "${BGEN_FILE}" \
  --sample "${SAMPLE_FILE}" \
  --phenoFile "${PHENO_FILE}" \
  --phenoCol "${PHENO_COLUMN}" \
  --covarFile "${COVAR_FILE}" \
  --covarColList age,pc1,pc2,pc3,pc4,pc5 \
  --catCovarList "${INTERACTION_COVARIATE}" \
  --interaction "${INTERACTION_TERM}" \
  --ref-first \
  --pred "${PRED_FILE}" \
  --keep "${KEEP_FILE}" \
  --extract "${EXTRACT_FILE}" \
  --gz \
  --bsize 400 \
  --bt \
  --firth \
  --approx \
  --firth-se \
  --pThresh 0.05 \
  --threads 6 \
  --out "${OUT_PREFIX}"
