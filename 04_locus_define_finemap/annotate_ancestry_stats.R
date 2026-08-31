#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    paste(
      "Usage: Rscript 04_locus_define_finemap/annotate_ancestry_stats.R",
      "<sig_alleles.csv> <mvp_eur.csv[.gz]> <mvp_afr.csv[.gz]>",
      "<mvp_amr.csv[.gz]> <output.csv>"
    ),
    call. = FALSE
  )
}

lead_file <- args[1]
eur_file <- args[2]
afr_file <- args[3]
amr_file <- args[4]
output_file <- args[5]

required_columns <- function(dat, required, file_name) {
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

read_ancestry_statistics <- function(path, ancestry) {
  dat <- fread(path)
  required_columns(
    dat,
    c("CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N"),
    path
  )

  # No TEST filtering is applied: the annotation uses all rows supplied in the
  # post-liftover REGENIE file.
  dat <- dat %>%
    transmute(
      Chr = normalize_chromosome(CHROM),
      Pos = as.numeric(GENPOS),
      RefAllele = as.character(ALLELE0),
      AltAllele = as.character(ALLELE1),
      AltFreq = as.numeric(A1FREQ),
      INFO = as.numeric(INFO),
      N = as.numeric(N)
    )

  duplicate_keys <- dat %>%
    count(Chr, Pos, RefAllele, AltAllele) %>%
    filter(n > 1)
  if (nrow(duplicate_keys) > 0) {
    stop(
      sprintf(
        "%s has more than one row for at least one variant. Supply one row per variant.",
        path
      ),
      call. = FALSE
    )
  }

  names(dat)[names(dat) == "AltFreq"] <- paste0(ancestry, "_AltFreq")
  names(dat)[names(dat) == "INFO"] <- paste0(ancestry, "_INFO")
  names(dat)[names(dat) == "N"] <- paste0(ancestry, "_N")
  dat[[paste0(ancestry, "_MAF")]] <- pmin(
    dat[[paste0(ancestry, "_AltFreq")]],
    1 - dat[[paste0(ancestry, "_AltFreq")]]
  )

  dat
}

lead_variants <- fread(lead_file)
required_columns(
  lead_variants,
  c("Chr", "Pos", "RefAllele", "AltAllele"),
  lead_file
)
lead_variants <- lead_variants %>%
  mutate(
    Chr = normalize_chromosome(Chr),
    Pos = as.numeric(Pos),
    RefAllele = as.character(RefAllele),
    AltAllele = as.character(AltAllele)
  )

join_columns <- c("Chr", "Pos", "RefAllele", "AltAllele")
annotated <- lead_variants %>%
  left_join(read_ancestry_statistics(eur_file, "EUR"), by = join_columns) %>%
  left_join(read_ancestry_statistics(afr_file, "AFR"), by = join_columns) %>%
  left_join(read_ancestry_statistics(amr_file, "AMR"), by = join_columns)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(annotated, output_file)
