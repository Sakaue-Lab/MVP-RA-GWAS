#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    paste(
      "Usage: Rscript 09_prs/prepare_clumped_weights.R",
      "<gwas_weights.csv> <clump_threshold_dir> <output_dir>"
    ),
    call. = FALSE
  )
}

gwas_file <- args[1]
clump_dir <- args[2]
output_dir <- args[3]

gwas <- fread(gwas_file)
required <- c("ID", "Beta", "P", "A1")
missing_columns <- setdiff(required, names(gwas))
if (length(missing_columns) > 0) {
  stop("Missing GWAS column(s): ", paste(missing_columns, collapse = ", "))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
n_clumps <- integer(22)

for (chromosome in seq_len(22)) {
  chromosome_dir <- file.path(clump_dir, paste0("chr", chromosome))
  clump_files <- list.files(
    chromosome_dir,
    pattern = "\\.clumps$",
    full.names = TRUE
  )

  if (length(clump_files) == 0) {
    message("No clump output for chromosome ", chromosome)
    next
  }

  clumps <- bind_rows(lapply(clump_files, fread))
  if (!("ID" %in% names(clumps))) {
    stop("PLINK clump output is missing the ID column: ", clump_files[1])
  }

  selected_ids <- unique(clumps$ID)
  chromosome_weights <- gwas %>%
    filter(ID %in% selected_ids) %>%
    select(ID, Beta, P, A1)

  n_clumps[chromosome] <- nrow(chromosome_weights)
  fwrite(
    chromosome_weights,
    file.path(output_dir, paste0(chromosome, ".gwasres.csv.gz"))
  )
}

message("Total clumped variants: ", sum(n_clumps))
