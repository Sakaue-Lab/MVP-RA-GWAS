#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    paste(
      "Usage: Rscript 04_locus_define_finemap/define_loci.R",
      "<trans.txt> <eur.txt> <trans_seropositive.txt>",
      "<eur_seropositive.txt> <output.csv>"
    ),
    call. = FALSE
  )
}

trans_file <- args[1]
eur_file <- args[2]
trans_seropositive_file <- args[3]
eur_seropositive_file <- args[4]
output_file <- args[5]

genome_wide_threshold <- 5e-8
locus_window <- 1e6
mhc_chromosome <- "chr6"
mhc_start <- 25e6
mhc_end <- 35e6

read_significant_variants <- function(path, analysis_name) {
  dat <- fread(path)
  required <- c(
    "MarkerName", "Effect", "StdErr", "P-value",
    "Direction", "Allele1", "Allele2"
  )
  missing_columns <- setdiff(required, names(dat))
  if (length(missing_columns) > 0) {
    stop(
      sprintf(
        "%s is missing required METAL column(s): %s",
        path,
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  marker_fields <- tstrsplit(dat$MarkerName, ":", fixed = TRUE)
  if (length(marker_fields) != 4) {
    stop(
      sprintf("MarkerName in %s must use chr:position:ref:alt format.", path),
      call. = FALSE
    )
  }

  dat <- dat %>%
    transmute(
      Chr = marker_fields[[1]],
      Pos = as.numeric(marker_fields[[2]]),
      RefAllele = marker_fields[[3]],
      AltAllele = marker_fields[[4]],
      Effect = as.numeric(Effect),
      StdErr = as.numeric(StdErr),
      PValue = as.numeric(`P-value`),
      Direction,
      Allele1,
      Allele2,
      ancestry = analysis_name
    ) %>%
    filter(
      is.finite(Pos),
      is.finite(Effect),
      is.finite(StdErr),
      is.finite(PValue),
      PValue <= genome_wide_threshold,
      !(Chr == mhc_chromosome & Pos >= mhc_start & Pos <= mhc_end),
      Chr != "chrX"
    ) %>%
    mutate(chisq_stat = (Effect / StdErr)^2) %>%
    arrange(desc(chisq_stat))

  dat
}

# Greedy physical-distance pruning. Input order defines priority: the first
# remaining variant is retained and same-chromosome variants within +/- 1 Mb
# are removed before the next iteration.
select_lead_variants <- function(dat, window_size = locus_window) {
  selected <- dat[0, , drop = FALSE]
  remaining <- dat

  while (nrow(remaining) > 0) {
    lead <- remaining[1, , drop = FALSE]
    selected <- bind_rows(selected, lead)

    remaining <- remaining %>%
      filter(
        Chr != lead$Chr |
          Pos < lead$Pos - window_size |
          Pos > lead$Pos + window_size
      )
  }

  selected
}

trans <- read_significant_variants(trans_file, "trans_combined") %>%
  select_lead_variants()

secondary_analyses <- bind_rows(
  read_significant_variants(eur_file, "eur_combined"),
  read_significant_variants(trans_seropositive_file, "trans_seropositive"),
  read_significant_variants(eur_seropositive_file, "eur_seropositive")
) %>%
  arrange(desc(chisq_stat))

# The overall trans-ancestry leads are placed first deliberately. Therefore,
# they define their loci before signals found only in secondary analyses.
selected_variants <- bind_rows(trans, secondary_analyses) %>%
  select_lead_variants()

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(selected_variants, output_file, row.names = FALSE)
