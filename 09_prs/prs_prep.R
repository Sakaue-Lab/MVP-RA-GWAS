#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    paste(
      "Usage: Rscript 09_prs/prs_prep.R <metal_output.txt>",
      "<snp_list.csv> <gwas_weights.csv>"
    ),
    call. = FALSE
  )
}

metal_file <- args[1]
snp_output <- args[2]
gwas_output <- args[3]

dat <- fread(metal_file)
required <- c("MarkerName", "Effect", "StdErr", "P-value", "Allele1", "Allele2")
missing_columns <- setdiff(required, names(dat))
if (length(missing_columns) > 0) {
  stop("Missing METAL column(s): ", paste(missing_columns, collapse = ", "))
}

marker_fields <- tstrsplit(dat$MarkerName, ":", fixed = TRUE)
if (length(marker_fields) != 4) {
  stop("MarkerName must use chr:position:ref:alt format.")
}

gwas <- dat %>%
  transmute(
    ID = as.character(MarkerName),
    Chr = marker_fields[[1]],
    Pos = as.numeric(marker_fields[[2]]),
    Beta = as.numeric(Effect),
    P = as.numeric(`P-value`),
    SE = as.numeric(StdErr),
    A1 = as.character(Allele1),
    A2 = as.character(Allele2)
  ) %>%
  filter(
    !is.na(Chr),
    nzchar(ID),
    is.finite(Pos),
    is.finite(Beta),
    is.finite(P),
    !(Chr == "chr6" & Pos >= 25e6 & Pos <= 35e6)
  )

dir.create(dirname(snp_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(gwas_output), recursive = TRUE, showWarnings = FALSE)
write.csv(select(gwas, ID), snp_output, row.names = FALSE, quote = FALSE)
write.csv(
  select(gwas, ID, Beta, P, SE, A1, A2),
  gwas_output,
  row.names = FALSE,
  quote = FALSE
)
