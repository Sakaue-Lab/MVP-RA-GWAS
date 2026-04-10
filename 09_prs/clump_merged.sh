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

p=0.00000005

for chr in $(seq 1 22);do
# Paths and parameters
gwasres=[PATH]/[FILE].csv # csv with columns: (1) SNP name, (2) beta, (3) p-value
datdir=[PATH]
odir=[PATH]/$p/$chr
mkdir -p "$odir"
dat=[PATH]
rsq=0.2

# Prepare SNP list
input_snps_csv="$datdir/[FILE].csv" # list of SNPs
tail -n +2 "$input_snps_csv" | cut -d, -f1 > "$odir/tmp01"

# Filter GWAS results to common SNPs
head -1 "$gwasres" > "$odir/tmp.gwasres.filt_common.csv"
tail -n +2 "$gwasres" |
  awk -F, 'NR==FNR { snps[$1]; next } $1 in snps' "$odir/tmp01" - >> "$odir/tmp.gwasres.filt_common.csv"

# Loop over top-percent thresholds
top=100
	total=$(( $(wc -l < "$odir/tmp.gwasres.filt_common.csv") - 1 ))
	keep=$(( total * top / 100 ))
	echo "Keeping top $keep SNPs ($top%)"

	# Create top-% GWAS file
	awkgwas="$odir/tmp.${top}.gwasres.csv"
	head -1 "$odir/tmp.gwasres.filt_common.csv" > "$awkgwas"
	tail -n +2 "$odir/tmp.gwasres.filt_common.csv" |
	  sort -t, -k3,3g |
	  head -n "$keep" >> "$awkgwas"

	  # Process each .pgen chunk
	  for chunk in $dat/merged/chr$chr/*.pgen; do
	    echo "Processing chunk: $chunk"

	    # Derive prefix for plink and sanitized base name for files
	    prefix="${chunk%.pgen}"
	    chunk_base=$(basename "$chunk" .pgen)

	    # Prepare per-chr GWAS input for PLINK
	    plink_input="$odir/tmp.${top}.${chr}.gwasres"
	    echo "SNP P" > "$plink_input"
	    awk -F, -v chr="$chr" '
		{ id = $1 
		  gsub(/^"|"$/, "", id) }
		id ~ ("^chr" chr ":") { print $1" "$3 }' \
	      "$awkgwas" >> "$plink_input"
	
	    if [[ ! -s "$plink_input" ]]; then
		echo "WARNING: $plink_input is empty; skipping chr$chr plink clump"
		continue
	    fi	

	    # Run PLINK2 clumping
	    plink2a \
		--pfile "$prefix" \
		--clump "$plink_input" \
		--clump-p1 "$p" \
		--clump-r2 "$rsq" \
		--clump-unphased \
		--chr "$chr" \
		--threads 1 \
		--memory 8000 \
		--out "$odir/${top}.${chr}.${chunk_base}"

	    # Extract kept SNPs
	    clumped_file="$odir/${top}.${chr}.${chunk_base}.clumps"
	    
	    # Guard: skip if no clumps found
	    if [[ ! -f "$clumped_file" ]]; then
		echo "No .clumps file found -- likely no significant clumps for chr$chr chunk $chunk_base"
		rm -f "$odir/${top}.${chr}.${chunk_base}.*"
		continue
	    fi
	    # Cleanup chunk temps
	    rm -f "$odir/tmp.${top}.${chr}.${chunk_base}.*"
	  done
	rm $odir/tmp*
done

