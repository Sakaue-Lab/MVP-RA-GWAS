#!/usr/bin/env Rscript
# BUHMBOX heterogeneity test on cell-type-specific PRS variant sets.
# For each set of RA risk variants, BUHMBOX (Han et al., Nat Genet 2016) tests
# for a latent case subgroup from excess positive correlation among risk-increasing
# alleles beyond that expected under homogeneity and LD.
#
# For each variant set we run BUHMBOX on RA cases (with controls to estimate the
# null allele-correlation structure) and write:
#   <set>_buhmbox.weights : LD-adjusted expected pairwise correlation matrix
#   <set>_buhmbox.BBrst   : test result (PVALUE LOG10P N NCASE NCONT ZSCORE NLOCUS ...)
#
# Usage: Rscript run_buhmbox.R
suppressPackageStartupMessages({library(data.table)})
source("/path_to/BUHMBOX/BUHMBOX.R")   # Han et al. 2016 implementation

SETS <- c("Bcell", "DC", "endothelial", "fibroblast", "monocyte", "Tcell",
          "unioncells", "all1e-05")     # per-cell-type, union, and full PRS set

# RA risk alleles + weights (risk-increasing allele and its GWAS effect)
risk <- fread("stats/prs_best_variants.risk_alleles.txt")   # varID risk_allele beta

# individual-level genotype dosage of RA cases and controls (variants x samples)
geno_cases    <- fread("stats/geno_cases.dosage.txt.gz")     # rows: varID ; cols: samples
geno_controls <- fread("stats/geno_controls.dosage.txt.gz")

run_one <- function(setname){
  vars <- readLines(sprintf("stats/celltype_variants/%s.variants.txt", setname))
  vars <- intersect(vars, risk$varID)
  rr   <- risk[match(vars, varID)]
  # orient dosages to the risk-increasing allele
  Xcase <- as.matrix(geno_cases[match(vars, varID), -1, with = FALSE])
  Xctrl <- as.matrix(geno_controls[match(vars, varID), -1, with = FALSE])

  res <- BUHMBOX(genotype_case   = t(Xcase),
                 genotype_control = t(Xctrl),
                 risk_allele_freq = rr$af,
                 gamma            = exp(rr$beta),      # per-allele OR as risk weights
                 out_prefix       = sprintf("stats/buhmbox/%s_buhmbox", setname))
  cat(sprintf("%-12s NLOCUS=%d  Z=%.3f  P=%.3g\n", setname, res$NLOCUS, res$ZSCORE, res$PVALUE))
  data.table(set = setname, NLOCUS = res$NLOCUS, ZSCORE = res$ZSCORE,
             PVALUE = res$PVALUE, N = res$N, NCASE = res$NCASE, NCONT = res$NCONT)
}

dir.create("stats/buhmbox", showWarnings = FALSE, recursive = TRUE)
summ <- rbindlist(lapply(SETS, run_one))
fwrite(summ, "stats/buhmbox_summary.tsv", sep = "\t")
cat("wrote stats/buhmbox_summary.tsv\n")
