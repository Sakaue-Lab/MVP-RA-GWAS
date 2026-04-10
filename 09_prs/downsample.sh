#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o xtrace
set +x

#BSUB -G mvp001
#BSUB -M 93000
#BSUB -a "multithread(8)"
#BSUB -q short
#BSUB -o [PATH]/%J.stdout
#BSUB -e [PATH]/%J.stderr

# Paths and parameters
chr=22
dat=[PATH_PGEN]
odir=[PATH]
mkdir -p $odir
keep=[PATH]/[FILE].txt # list of SNPs to keep from clumping

for chunk in "$dat/chr$chr"/*.pgen; do
	echo "Processing chunk: $chunk"
	prefix="${chunk%.pgen}"
	chunk_base=$(basename "$chunk" .pgen)

	plink2a \
	--pfile "$prefix" \
	--keep "$keep" \
	--make-pgen \
	--memory 32000 \
	--out "$odir/${chunk_base}_down"
done


