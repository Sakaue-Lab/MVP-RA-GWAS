#usage: bsub < regenie_step1.afr.sh

#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o xtrace

#BSUB -J step1
#BSUB -G mvp001
#BSUB -M 10000
#BSUB -a "multithread(8)"
#BSUB -q regenie_step1
#BSUB -o [PATH]/%J.stdout
#BSUB -e [PATH]/%J.stderr

datadir=[PATH]
phedir=[PATH]
resdir=[PATH]

HARE=AMR # ancestry

BED=$datadir/[FILE] # BED file prefix
EXTRACT=$datadir/[FILE]
KEEP=$phedir/[FILE].txt # SNPs to keep
PHENO=$phedir/[FILE].txt # phenotype file
COVAR=$phedir/[FILE].txt # covariate file

OUT=$resdir/RA.$HARE.step1

# Single variant association tests 
# Step 1: fitting the null logistic/linear mixed model
# input plink file

[PATH]/regenie/regenie-3.1.3/regenie --step 1 --bed $BED --extract $EXTRACT --phenoFile $PHENO --covarFile $COVAR --keep $KEEP --bsize 1000 --bt --lowmem --loocv --threads 4 --verbose --out $OUT
