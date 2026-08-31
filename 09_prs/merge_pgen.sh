#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 <chromosome> <downsampled_dir> <output_dir> [plink2]" >&2
  exit 1
fi

chromosome=$1
input_dir=$2
output_dir=$3
plink_bin=${4:-plink2}

mkdir -p "${output_dir}/chr${chromosome}"
prefix_list=$(mktemp "${TMPDIR:-/tmp}/pgen-prefixes.XXXXXX")
trap 'rm -f "${prefix_list}"' EXIT

shopt -s nullglob
chunks=("${input_dir}/chr${chromosome}"/*.pgen)
if [[ ${#chunks[@]} -eq 0 ]]; then
  echo "No downsampled PGEN chunks found for chromosome ${chromosome}." >&2
  exit 1
fi

for chunk in "${chunks[@]}"; do
  printf '%s\n' "${chunk%.pgen}" >> "${prefix_list}"
done

"${plink_bin}" \
  --pmerge-list "${prefix_list}" \
  --make-pgen \
  --threads 8 \
  --memory 8000 \
  --out "${output_dir}/chr${chromosome}/merged_chr${chromosome}"
