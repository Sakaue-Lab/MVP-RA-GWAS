#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: $0 <chromosome> <chunked_pgen_dir> <sample_keep.txt> <output_dir> [plink2]" >&2
  exit 1
fi

chromosome=$1
input_dir=$2
keep_file=$3
output_dir=$4
plink_bin=${5:-plink2}

mkdir -p "${output_dir}/chr${chromosome}"
shopt -s nullglob
chunks=("${input_dir}/chr${chromosome}"/*.pgen)
if [[ ${#chunks[@]} -eq 0 ]]; then
  echo "No PGEN chunks found for chromosome ${chromosome}." >&2
  exit 1
fi

for chunk in "${chunks[@]}"; do
  prefix=${chunk%.pgen}
  chunk_base=$(basename "${chunk}" .pgen)
  echo "Downsampling participants in ${chunk_base}"
  "${plink_bin}" \
    --pfile "${prefix}" \
    --keep "${keep_file}" \
    --make-pgen \
    --memory 32000 \
    --out "${output_dir}/chr${chromosome}/${chunk_base}_down"
done
