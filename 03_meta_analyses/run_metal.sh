#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

usage() {
  echo "Usage: METAL_BIN=/path/to/metal OUTPUT_PREFIX=results/name bash $0 <regenie.csv[.gz]> [more cohorts ...]" >&2
}

if [[ $# -lt 1 ]]; then
  usage
  exit 1
fi

METAL_BIN=${METAL_BIN:-generic-metal/metal}
: "${OUTPUT_PREFIX:?Set OUTPUT_PREFIX to the desired METAL output prefix}"

if [[ ! -x "${METAL_BIN}" ]]; then
  echo "METAL executable not found or not executable: ${METAL_BIN}" >&2
  exit 1
fi

for cohort_file in "$@"; do
  if [[ ! -f "${cohort_file}" ]]; then
    echo "Input cohort file does not exist: ${cohort_file}" >&2
    exit 1
  fi
done

output_dir=$(dirname "${OUTPUT_PREFIX}")
mkdir -p "${output_dir}"

temporary_dir=$(mktemp -d "${TMPDIR:-/tmp}/mvp-ra-metal.XXXXXX")
trap 'rm -rf "${temporary_dir}"' EXIT
metal_script="${temporary_dir}/metal_script.txt"

read_cohort_file() {
  local cohort_file=$1
  if [[ "${cohort_file}" == *.gz ]]; then
    gzip -cd -- "${cohort_file}"
  else
    cat -- "${cohort_file}"
  fi
}

prepare_regenie_input() {
  local input_file=$1
  local output_file=$2

  read_cohort_file "${input_file}" | awk -F',' -v OFS='\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        gsub(/^[[:space:]\"]+|[[:space:]\"\r]+$/, "", $i)
        column[$i] = i
      }

      required_count = split("ID ALLELE0 ALLELE1 BETA SE PVAL N TEST", required, " ")
      for (i = 1; i <= required_count; i++) {
        if (!(required[i] in column)) {
          print "Missing required REGENIE column " required[i] " in " FILENAME > "/dev/stderr"
          exit 2
        }
      }

      print "SNP", "EFFECT", "OTHER", "BETA", "PVALUE", "SE", "N"
      next
    }

    {
      for (i = 1; i <= NF; i++) {
        gsub(/^[[:space:]\"]+|[[:space:]\"\r]+$/, "", $i)
      }

      # REGENIE can emit multiple tests; this analysis meta-analyzes ADD only.
      if ($(column["TEST"]) != "ADD") next

      beta = $(column["BETA"])
      se = $(column["SE"])
      pvalue = $(column["PVAL"])
      n = $(column["N"])
      if (beta == "" || beta == "NA" || se == "" || se == "NA" ||
          pvalue == "" || pvalue == "NA" || n == "" || n == "NA") next

      # In REGENIE output, ALLELE1 is the tested/effect allele and ALLELE0 is
      # the other/reference allele.
      print $(column["ID"]), $(column["ALLELE1"]), $(column["ALLELE0"]),
            beta, pvalue, se, n
    }
  ' > "${output_file}"
}

printf 'SEPARATOR TAB\n' > "${metal_script}"
printf 'SCHEME STDERR\n' >> "${metal_script}"
printf 'MARKER SNP\n' >> "${metal_script}"
printf 'ALLELE EFFECT OTHER\n' >> "${metal_script}"
printf 'EFFECT BETA\n' >> "${metal_script}"
printf 'PVALUE PVALUE\n' >> "${metal_script}"
printf 'STDERR SE\n' >> "${metal_script}"
printf 'CUSTOMVARIABLE TOTALN\n' >> "${metal_script}"
printf 'LABEL TOTALN AS N\n' >> "${metal_script}"

cohort_number=0
for cohort_file in "$@"; do
  cohort_number=$((cohort_number + 1))
  prepared_file="${temporary_dir}/cohort_${cohort_number}.tsv"
  prepare_regenie_input "${cohort_file}" "${prepared_file}"
  printf 'PROCESS %s\n' "${prepared_file}" >> "${metal_script}"
done

# METAL requires a space between the output prefix and suffix.
printf 'OUTFILE %s_out .txt\n' "${OUTPUT_PREFIX}" >> "${metal_script}"
printf 'ANALYZE\n' >> "${metal_script}"
printf 'QUIT\n' >> "${metal_script}"

echo "Running METAL on ${cohort_number} cohort file(s)."
echo "Output prefix: ${OUTPUT_PREFIX}_out"
"${METAL_BIN}" < "${metal_script}"
