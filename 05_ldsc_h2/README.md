[← Repository overview](../README.md)

# 05. Heritability and inflation (LDSC)

Estimate SNP-based heritability and quantify confounding/inflation of the RA
meta-analyses with LD score regression (LDSC).

- `ldsc.sh` : munge the meta-analysis summary statistics (HapMap3 SNPs, MHC
  removed), run LDSC against the 1000G Phase 3 baselineLD model with EUR
  frequency and weight files, and collect the total observed-scale h2, mean
  chi-square, LDSC intercept, and ratio into a summary table.
- `liability_scale.R` : convert the observed-scale h2 (and its SE) to the
  liability scale given the population prevalence (K) and the case fraction in
  the sample.

Replace the `${munge}`/`${ldsc}` tool paths and the `LDSCORE` reference paths
with your local installation.
