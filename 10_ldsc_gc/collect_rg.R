#!/usr/bin/env Rscript
# Parse the LDSC --rg log into a tidy rg_summary.csv (one row per trait).
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
logf <- args[1]; out <- if (length(args) >= 2) args[2] else "stats/rg_summary.csv"

L <- readLines(logf)
# LDSC prints a whitespace-delimited "Summary of Genetic Correlation Results" table
i <- grep("Summary of Genetic Correlation Results", L)
stopifnot(length(i) == 1)
tab <- L[(i + 1):length(L)]
tab <- tab[tab != "" & !grepl("^Analysis finished", tab)]
dt <- fread(text = paste(tab, collapse = "\n"))
# columns: p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se
dt[, trait := gsub("\\.sumstats\\.gz$", "", basename(p2))]
fwrite(dt[, .(trait, rg, se, z, p, h2_obs, h2_obs_se, h2_int, h2_int_se,
              gcov_int, gcov_int_se)], out)
cat(sprintf("wrote %s (%d traits)\n", out, nrow(dt)))
