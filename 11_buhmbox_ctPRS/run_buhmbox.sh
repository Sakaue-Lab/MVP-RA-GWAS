#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

if [[ $# -lt 6 || $# -gt 8 ]]; then
  echo "Usage: bash run_buhmbox.sh <cell_type> <pgen_prefix> <phenotype.txt> <snp_bed> <covariates.txt|-> <output_dir> [plink2] [buhmbox.R]" >&2
  exit 1
fi

CELL_TYPE=$1
PGEN_PREFIX=$2
PHENO_FILE=$3
SNP_BED=$4
COVAR_FILE=$5
OUTPUT_ROOT=$6
PLINK2_BIN=${7:-plink2}
BUHMBOX_R=${8:-$(cd "$(dirname "$0")" && pwd)/buhmbox.R}

OUT_DIR="${OUTPUT_ROOT}/${CELL_TYPE}"
TMP_DIR="${OUT_DIR}/tmp"
OUT_PREFIX="${OUT_DIR}/${CELL_TYPE}_buhmbox"

for required_file in \
  "${PGEN_PREFIX}.pgen" "${PGEN_PREFIX}.pvar" "${PGEN_PREFIX}.psam" \
  "${PHENO_FILE}" "${SNP_BED}" "${BUHMBOX_R}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Required input does not exist: ${required_file}" >&2
    exit 1
  fi
done

if [[ "${COVAR_FILE}" != "-" && ! -f "${COVAR_FILE}" ]]; then
  echo "Covariate file does not exist: ${COVAR_FILE}" >&2
  exit 1
fi

if ! command -v "${PLINK2_BIN}" >/dev/null 2>&1 && [[ ! -x "${PLINK2_BIN}" ]]; then
  echo "PLINK 2 executable was not found: ${PLINK2_BIN}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

# BED columns: CHROM, START, END, ID, EFFECT_ALLELE, BETA, OR.
# BUHMBOX SNP columns: ID, risk/effect allele, frequency, OR.
awk -v OFS="\t" 'NF >= 7 {print $4, $5, "NA", $7}' "${SNP_BED}" \
  > "${TMP_DIR}/buhmbox_snps.txt"
awk '{print $1}' "${TMP_DIR}/buhmbox_snps.txt" > "${TMP_DIR}/snp_ids.txt"

# Phenotype columns are FID, IID, STATUS with 1=control and 2=case.
awk '$3 == 1 {print $1, $2}' "${PHENO_FILE}" > "${TMP_DIR}/controls.txt"
awk '$3 == 2 {print $1, $2}' "${PHENO_FILE}" > "${TMP_DIR}/cases.txt"

if [[ ! -s "${TMP_DIR}/controls.txt" || ! -s "${TMP_DIR}/cases.txt" ]]; then
  echo "The phenotype file must contain at least one control (1) and one case (2)." >&2
  exit 1
fi

# Estimate LD in controls and retain variants pruned at r2 < 0.1.
"${PLINK2_BIN}" \
  --pfile "${PGEN_PREFIX}" \
  --keep "${TMP_DIR}/controls.txt" \
  --extract "${TMP_DIR}/snp_ids.txt" \
  --indep-pairwise 1000kb 0.1 \
  --out "${TMP_DIR}/ld_prune" \
  --threads 4

if [[ ! -s "${TMP_DIR}/ld_prune.prune.in" ]]; then
  echo "LD pruning retained no variants." >&2
  exit 1
fi

awk 'NR == FNR {keep[$1] = 1; next} ($1 in keep)' \
  "${TMP_DIR}/ld_prune.prune.in" "${TMP_DIR}/buhmbox_snps.txt" \
  > "${TMP_DIR}/buhmbox_snps_pruned.txt"
cp "${TMP_DIR}/ld_prune.prune.in" "${TMP_DIR}/snp_ids_pruned.txt"
awk '{print $1, $2}' "${TMP_DIR}/buhmbox_snps_pruned.txt" \
  > "${TMP_DIR}/risk_alleles.txt"

for group in cases controls; do
  "${PLINK2_BIN}" \
    --pfile "${PGEN_PREFIX}" \
    --keep "${TMP_DIR}/${group}.txt" \
    --extract "${TMP_DIR}/snp_ids_pruned.txt" \
    --export A \
    --export-allele "${TMP_DIR}/risk_alleles.txt" \
    --out "${TMP_DIR}/${group}" \
    --threads 4
done

buhmbox_command=(
  Rscript "${BUHMBOX_R}"
  "${TMP_DIR}/buhmbox_snps_pruned.txt"
  "${TMP_DIR}/cases.raw"
  "${TMP_DIR}/controls.raw"
  YY Y Y
  "${OUT_PREFIX}"
)

if [[ "${COVAR_FILE}" != "-" ]]; then
  # Expected columns: FID IID sex age PC1 ... PC10.
  awk '{print $1, $2, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' \
    "${COVAR_FILE}" > "${TMP_DIR}/pcs.txt"
  buhmbox_command+=("${TMP_DIR}/pcs.txt")
fi

"${buhmbox_command[@]}"

if [[ ! -f "${OUT_PREFIX}.BBrst" ]]; then
  echo "BUHMBOX did not create the expected result: ${OUT_PREFIX}.BBrst" >&2
  exit 1
fi

echo "BUHMBOX result: ${OUT_PREFIX}.BBrst"
echo "BUHMBOX log: ${OUT_PREFIX}.log"
