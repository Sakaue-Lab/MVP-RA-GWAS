#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    paste(
      "Usage: Rscript 04_locus_define_finemap/fine_map.R",
      "<metal_output.txt> <lead_variants.csv> <output.csv>"
    ),
    call. = FALSE
  )
}

metal_file <- args[1]
lead_file <- args[2]
output_file <- args[3]

window_size <- 5e5
prior_variance <- 0.04
credible_level <- 0.95

require_columns <- function(dat, required, file_name) {
  missing_columns <- setdiff(required, names(dat))
  if (length(missing_columns) > 0) {
    stop(
      sprintf(
        "%s is missing required column(s): %s",
        file_name,
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

normalize_chromosome <- function(x) {
  sub("^chr", "", as.character(x), ignore.case = TRUE)
}

metal <- fread(metal_file)
require_columns(
  metal,
  c("MarkerName", "Effect", "StdErr", "P-value", "Allele1", "Allele2"),
  metal_file
)

marker_fields <- tstrsplit(metal$MarkerName, ":", fixed = TRUE)
if (length(marker_fields) != 4) {
  stop("METAL MarkerName must use chr:position:ref:alt format.", call. = FALSE)
}

metal <- metal %>%
  transmute(
    Chr = normalize_chromosome(marker_fields[[1]]),
    Pos = as.numeric(marker_fields[[2]]),
    RefAllele = marker_fields[[3]],
    AltAllele = marker_fields[[4]],
    SNP = paste(Chr, Pos, RefAllele, AltAllele, sep = ":"),
    Beta = as.numeric(Effect),
    SE = as.numeric(StdErr),
    P = as.numeric(`P-value`),
    Allele1,
    Allele2
  ) %>%
  filter(
    is.finite(Pos),
    is.finite(Beta),
    is.finite(SE),
    SE > 0
  )

lead_variants <- fread(lead_file)
require_columns(
  lead_variants,
  c("Chr", "Pos", "RefAllele", "AltAllele"),
  lead_file
)
lead_variants <- lead_variants %>%
  mutate(
    Chr = normalize_chromosome(Chr),
    Pos = as.numeric(Pos),
    locus_id = paste(Chr, Pos, RefAllele, AltAllele, sep = ":")
  )

fine_map_locus <- function(lead) {
  locus <- metal %>%
    filter(
      Chr == lead$Chr,
      Pos >= lead$Pos - window_size,
      Pos <= lead$Pos + window_size
    )

  if (nrow(locus) == 0) {
    warning("No METAL variants found within 500 kb of lead variant ", lead$locus_id)
    return(NULL)
  }

  # Wakefield approximate Bayes factor with prior variance w = 0.04.
  z_squared <- (locus$Beta / locus$SE)^2
  locus$log_abf <-
    0.5 * log(locus$SE^2 / (locus$SE^2 + prior_variance)) +
    prior_variance * z_squared / (2 * (locus$SE^2 + prior_variance))

  # Subtracting the maximum log-ABF prevents numerical overflow without
  # changing the normalized posterior inclusion probabilities.
  scaled_abf <- exp(locus$log_abf - max(locus$log_abf))
  locus$pip <- scaled_abf / sum(scaled_abf)
  locus <- locus %>%
    arrange(desc(pip)) %>%
    mutate(cumulative_pip = cumsum(pip))

  last_index <- which(locus$cumulative_pip >= credible_level)[1]
  locus <- locus[seq_len(last_index), , drop = FALSE]
  locus$credible_set_rank <- seq_len(nrow(locus))
  locus$credible_set_size <- nrow(locus)
  locus$lead_chr <- lead$Chr
  locus$lead_pos <- lead$Pos
  locus$lead_ref_allele <- lead$RefAllele
  locus$lead_alt_allele <- lead$AltAllele
  locus$locus_id <- lead$locus_id
  locus
}

credible_sets <- bind_rows(lapply(
  seq_len(nrow(lead_variants)),
  function(i) fine_map_locus(lead_variants[i, , drop = FALSE])
))

if (nrow(credible_sets) == 0) {
  stop("No credible sets could be constructed.", call. = FALSE)
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(credible_sets, output_file)
