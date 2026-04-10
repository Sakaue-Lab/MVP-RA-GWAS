#!/bin/bash

script="metal_script_trans.txt"

# Define the cohort file paths
cohort_1="[FILE].tsv"
cohort_2="[FILE].tsv"
cohort_3="[FILE].tsv"
cohort_4="[FILE].tsv"

# Define the output prefix
meta_stats_prefix="[PATH]"

# Start writing the METAL script
echo -e "SEPARATOR TAB" > $script
echo -e "SCHEME STDERR" >> $script
echo -e "MARKER SNP" >> $script
echo -e "ALLELE EFFECT OTHER" >> $script
echo -e "EFFECT BETA" >> $script
echo -e "PVALUE PVALUE" >> $script
echo -e "STDERR SE" >> $script
echo -e "CUSTOMVARIABLE TOTALN" >> $script
echo -e "LABEL TOTALN AS N" >> $script

# Add cohort files
for cohort in "$cohort_1" "$cohort_2" "$cohort_3" "$cohort_4"; do
  echo -e "PROCESS $cohort" >> $script
done

# Specify output (leave the random space!)
echo -e "OUTFILE ${meta_stats_prefix}_out .txt" >> $script
echo -e "ANALYZE" >> $script
echo -e "QUIT" >> $script

# Debugging: Display the script content
echo "Generated METAL script:"
cat $script

# Debugging: Confirm the METAL command
echo "Running METAL with: generic-metal/metal < $script"

