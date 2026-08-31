#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: $0 <chromosome> <full_chunked_pgen_dir> <chromosome_weights.csv.gz> <output_dir> [plink2]" >&2
  exit 1
fi

chromosome=$1
pgen_dir=$2
weights_file=$3
output_dir=$4
plink_bin=${5:-plink2}

mkdir -p "${output_dir}/chr${chromosome}"
if [[ ! -s "${weights_file}" ]]; then
  echo "Skipping chromosome ${chromosome}: no clumped weight file."
  exit 0
fi
score_file=$(mktemp "${TMPDIR:-/tmp}/prs-score.XXXXXX")
trap 'rm -f "${score_file}"' EXIT

# prepare_clumped_weights.R writes ID,Beta,P,A1. PLINK expects ID,A1,Beta.
gzip -cd -- "${weights_file}" | awk -F',' -v OFS='\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) column[$i] = i
    print "ID", "A1", "BETA"
    next
  }
  { print $(column["ID"]), $(column["A1"]), $(column["Beta"]) }
' > "${score_file}"

shopt -s nullglob
chunks=("${pgen_dir}/chr${chromosome}"/*.pgen)
if [[ ${#chunks[@]} -eq 0 ]]; then
  echo "No full-cohort PGEN chunks found for chromosome ${chromosome}." >&2
  exit 1
fi

for chunk in "${chunks[@]}"; do
  prefix=${chunk%.pgen}
  chunk_base=$(basename "${chunk}" .pgen)
  chunk_variants=$(mktemp "${TMPDIR:-/tmp}/chunk-variants.XXXXXX")
  chunk_scores=$(mktemp "${TMPDIR:-/tmp}/chunk-scores.XXXXXX")

  awk 'NR > 1 {print $3}' "${prefix}.pvar" > "${chunk_variants}"
  awk -v OFS='\t' 'NR == FNR {ids[$1]; next} FNR == 1 || ($1 in ids)' \
    "${chunk_variants}" "${score_file}" > "${chunk_scores}"

  if [[ $(wc -l < "${chunk_scores}") -gt 1 ]]; then
    "${plink_bin}" \
      --pfile "${prefix}" \
      --score "${chunk_scores}" header cols=scoresums ignore-dup-ids \
      --threads 1 \
      --memory 8000 \
      --silent \
      --out "${output_dir}/chr${chromosome}/${chunk_base}"
    gzip -f "${output_dir}/chr${chromosome}/${chunk_base}.sscore"
  else
    echo "Skipping ${chunk_base}: no selected variants in this chunk."
  fi

  rm -f "${chunk_variants}" "${chunk_scores}"
done
