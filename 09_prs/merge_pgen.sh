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

chr=21

mkdir -p [PATH]

ls [PATH]/*.pgen |
	sed 's/\.pgen$//' > tmp.chunks.prefixes.chr"$chr".txt

plink2a \
--pmerge-list tmp.chunks.prefixes.chr"$chr".txt \
--make-pgen \
--threads 8 \
--memory 8000 \
--out [PATH]/merged/chr$chr/merged_chr$chr
