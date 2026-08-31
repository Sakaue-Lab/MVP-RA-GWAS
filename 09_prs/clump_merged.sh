#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: $0 <merged_pgen_dir> <gwas_weights.csv> <p_threshold> <output_dir> [plink2]" >&2
  exit 1
fi

merged_dir=$1
gwas_file=$2
p_threshold=$3
output_dir=$4
plink_bin=${5:-plink2}
r_squared=0.2

for chromosome in $(seq 1 22); do
  pgen_prefix="${merged_dir}/chr${chromosome}/merged_chr${chromosome}"
  if [[ ! -f "${pgen_prefix}.pgen" ]]; then
    echo "Skipping chromosome ${chromosome}: merged PGEN not found."
    continue
  fi

  chromosome_dir="${output_dir}/${p_threshold}/chr${chromosome}"
  mkdir -p "${chromosome_dir}"
  clump_input=$(mktemp "${TMPDIR:-/tmp}/clump-input.XXXXXX")
  trap 'rm -f "${clump_input}"' EXIT

  printf 'ID\tP\n' > "${clump_input}"
  awk -F',' -v chromosome="chr${chromosome}" -v OFS='\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) column[$i] = i
      if (!("ID" in column) || !("P" in column)) exit 2
      next
    }
    $(column["ID"]) ~ ("^" chromosome ":") {
      print $(column["ID"]), $(column["P"])
    }
  ' "${gwas_file}" >> "${clump_input}"

  if [[ $(wc -l < "${clump_input}") -le 1 ]]; then
    echo "No GWAS variants found for chromosome ${chromosome}."
    rm -f "${clump_input}"
    trap - EXIT
    continue
  fi

  "${plink_bin}" \
    --pfile "${pgen_prefix}" \
    --clump "${clump_input}" \
    --clump-p1 "${p_threshold}" \
    --clump-r2 "${r_squared}" \
    --clump-unphased \
    --chr "${chromosome}" \
    --threads 1 \
    --memory 8000 \
    --out "${chromosome_dir}/merged_chr${chromosome}"

  rm -f "${clump_input}"
  trap - EXIT
done
