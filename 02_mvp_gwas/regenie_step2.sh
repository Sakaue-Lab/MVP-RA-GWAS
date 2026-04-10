##  usage: bsub -J "RA_s2.amr[1-289]" < regenie_step2.afr.sh

#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o xtrace

#BSUB -J step2
#BSUB -G mvp001
#BSUB -M 15000
#BSUB -a "multithread(8)"
#BSUB -q regenie_step2
#BSUB -o [PATH]/%J.stdout
#BSUB -e [PATH]/%J.stderr

genodir=[PATH] # BGEN file directory
datadir=[PATH]
phedir=[PATH] # phenotype file directory
resdir=[PATH]

SEEDFILE=$datadir/[FILE].txt
SEED=$(cat $SEEDFILE | head -n $LSB_JOBINDEX | tail -n 1)
COLS=($SEED)
CHR=${COLS[0]}
CHUNK=${COLS[1]}

HARE=AMR # ancestry

BGEN=$genodir/$CHR/$CHUNK.bgen
SAM=$genodir/$CHR/$CHUNK.sample

EXTRACT=$datadir/[FILE]snplist.$CHR.txt
KEEP=$phedir/[FILE].$HARE.txt

PHENO=$phedir/[FILE].txt
COVAR=$phedir/[FILE].txt
PRED=$resdir/RA.$HARE.step1_pred.list

OUT=$resdir/$HARE.$CHUNK.step2

# Single variant association tests 
# Step 2: running association test

[PATH]/regenie/regenie-2.2.4/regenie --step 2 --bgen $BGEN --sample $SAM --phenoFile $PHENO --covarFile $COVAR --ref-first --pred $PRED --keep $KEEP --extract $EXTRACT --gz --bsize 400 --bt --firth --approx --firth-se --pThresh 0.05 --approx --threads 6 --out $OUT
