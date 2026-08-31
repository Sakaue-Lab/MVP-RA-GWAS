#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    paste(
      "Usage: Rscript 04_locus_define_finemap/annotate_previous_study.R",
      "<lead_variants.csv> <previous_study_loci.csv> <output.csv>"
    ),
    call. = FALSE
  )
}

lead_file <- args[1]
previous_study_file <- args[2]
output_file <- args[3]
novelty_window <- 1e6

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

find_variant_id_column <- function(dat, file_name) {
  candidates <- c("VariantID", "Variant.ID", "ID", "MarkerName")
  matches <- candidates[candidates %in% names(dat)]
  if (length(matches) == 0) {
    stop(
      sprintf(
        "%s must contain one variant-ID column named one of: %s",
        file_name,
        paste(candidates, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  matches[1]
}

parse_variant_ids <- function(ids, file_name) {
  fields <- tstrsplit(as.character(ids), ":", fixed = TRUE)
  if (length(fields) != 4) {
    stop(
      sprintf("Variant IDs in %s must use chr:position:ref:alt format.", file_name),
      call. = FALSE
    )
  }

  data.frame(
    previous_chr = normalize_chromosome(fields[[1]]),
    previous_pos = as.numeric(fields[[2]]),
    previous_ref = fields[[3]],
    previous_alt = fields[[4]],
    stringsAsFactors = FALSE
  )
}

lead_variants <- fread(lead_file)
require_columns(
  lead_variants,
  c("Chr", "Pos", "RefAllele", "AltAllele"),
  lead_file
)
lead_variants <- lead_variants %>%
  mutate(
    Chr = normalize_chromosome(Chr),
    Pos = as.numeric(Pos)
  )

previous_study <- fread(previous_study_file)
variant_id_column <- find_variant_id_column(previous_study, previous_study_file)
parsed_previous_ids <- parse_variant_ids(
  previous_study[[variant_id_column]],
  previous_study_file
)
if (any(!is.finite(parsed_previous_ids$previous_pos)) ||
    any(is.na(parsed_previous_ids$previous_chr)) ||
    any(is.na(parsed_previous_ids$previous_ref)) ||
    any(is.na(parsed_previous_ids$previous_alt))) {
  stop(
    sprintf("At least one variant ID in %s is malformed.", previous_study_file),
    call. = FALSE
  )
}

# Retain all supplied previous-study metadata, adding a generic prefix so its
# origin remains clear after it is joined to the new lead-variant table.
previous_metadata <- previous_study
names(previous_metadata) <- paste0("previous_nearest_", names(previous_metadata))
previous_study <- bind_cols(parsed_previous_ids, previous_metadata)

nearest_rows <- lapply(seq_len(nrow(lead_variants)), function(i) {
  candidates <- previous_study %>%
    filter(previous_chr == lead_variants$Chr[i])

  if (nrow(candidates) == 0) {
    empty <- previous_study[NA_integer_, , drop = FALSE]
    empty$distance_to_previous_locus <- NA_real_
    return(empty)
  }

  distances <- abs(candidates$previous_pos - lead_variants$Pos[i])
  nearest <- candidates[which.min(distances), , drop = FALSE]
  nearest$distance_to_previous_locus <- min(distances)
  nearest
})

nearest_previous <- bind_rows(nearest_rows)
annotated <- bind_cols(lead_variants, nearest_previous) %>%
  mutate(
    novel_vs_previous_study = is.na(distance_to_previous_locus) |
      distance_to_previous_locus > novelty_window
  )

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(annotated, output_file)
