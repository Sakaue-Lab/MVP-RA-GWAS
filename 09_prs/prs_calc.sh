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

top=100
chr=1

prs_out=[PATH]
dat=[PATH]
odir=[PATH]
mkdir -p "$odir"/$chr

# Change depending on p-val cutoff
for pth in 1e-01;do
	for chunk in "$dat/chr${chr}"/*.pgen; do
	echo "Processing chunk: $chunk"
		# filter beta info file by p values
		echo -e "SNP\tA1\tBETA" > $odir/$chr/tmp."$chr".beta.info

		zcat $prs_out/"$chr".gwasres.csv.gz |
		sed -e "1d" |
		awk -F, -v pth=$pth '{
			if( $3 < pth ){
				print $1"\t"$4"\t"$2
			}
		}' >> $odir/$chr/tmp."$chr".beta.info

		# Derive prefix for plink
		prefix="${chunk%.pgen}"
		chunk_base=$(basename "$chunk" .pgen)

		# Check whether there are target SNPs in this chunk
		cut -f3 "$prefix".pvar > $odir/$chr/tmp.chunk"$chunk_base".variants
		echo -e "SNP\tA1\tBETA" > $odir/$chr/tmp.chr"$chr".chunk"$chunk_base".beta.info
		grep -F -w -f $odir/$chr/tmp.chunk"$chunk_base".variants $odir/$chr/tmp."$chr".beta.info >> $odir/$chr/tmp.chr"$chr".chunk"$chunk_base".beta.info || true
		
		N=$( cat $odir/$chr/tmp.chr"$chr".chunk"$chunk_base".beta.info | wc -l )
		echo "calculate N is $N"
		if [ $N -gt 1 ] ;then
			echo "PLINK calculation"
			# Calculate PRS
			plink2 \
			--pfile $prefix \
			--score $odir/$chr/tmp.chr"$chr".chunk"$chunk_base".beta.info header cols=scoresums ignore-dup-ids \
			--threads 1 \
			--memory 8000 \
			--silent \
			--out $odir/$chr/$chunk_base

			gzip -f $odir/$chr/"$chunk_base".sscore
		else
			echo "skip PLINK"
			
		fi

		rm -f $odir/$chr/tmp.*
	done
done
		
