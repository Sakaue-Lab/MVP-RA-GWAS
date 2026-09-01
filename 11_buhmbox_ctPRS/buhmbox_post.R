#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(
    "Usage: Rscript buhmbox_post.R <cell_types.txt> <cell_score_dir> ",
    "<full_prs.csv> <training_ids.txt> <covariates.txt> <serology.txt> ",
    "<output.csv>"
  )
}

cell_file <- args[1]
cell_score_dir <- args[2]
full_prs_file <- args[3]
training_file <- args[4]
covariate_file <- args[5]
serology_file <- args[6]
output_file <- args[7]

require_columns <- function(dat, columns, label) {
  missing_columns <- setdiff(columns, names(dat))
  if (length(missing_columns) > 0) {
    stop(label, " is missing column(s): ", paste(missing_columns, collapse = ", "))
  }
}

cells <- fread(cell_file, header = FALSE)[[1]]
cells <- cells[nzchar(cells)]
if (length(cells) == 0) stop("The cell-type file is empty.")

full_prs <- fread(full_prs_file)
require_columns(full_prs, c("id", "raw_score"), "Full PRS file")
full_prs <- full_prs %>%
  transmute(id = as.character(id), raw_score_full = raw_score)

training_ids <- fread(training_file)
require_columns(training_ids, "IID", "Training-ID file")
training_ids$IID <- as.character(training_ids$IID)

covariates <- fread(covariate_file)
require_columns(
  covariates,
  c("id", "sex", "age", paste0("pc", 1:5)),
  "Covariate file"
)
covariates$id <- as.character(covariates$id)

serology <- fread(serology_file)
require_columns(
  serology,
  c("id", "ccflab", "rflab", "seropositive"),
  "Serology file"
)
serology$id <- as.character(serology$id)

results <- lapply(cells, function(cell_type) {
  score_file <- file.path(
    cell_score_dir,
    cell_type,
    paste0("scores_", cell_type, ".csv")
  )
  scores <- fread(score_file)
  require_columns(scores, c("id", "raw_score"), score_file)

  analysis_data <- scores %>%
    mutate(id = as.character(id)) %>%
    filter(!(id %in% training_ids$IID)) %>%
    inner_join(covariates, by = "id") %>%
    inner_join(serology, by = "id") %>%
    filter(!is.na(ccflab) | !is.na(rflab)) %>%
    mutate(seropositive = if_else(is.na(seropositive), 0, seropositive)) %>%
    inner_join(full_prs, by = "id")

  model <- glm(
    seropositive ~ raw_score + raw_score_full + sex + age +
      pc1 + pc2 + pc3 + pc4 + pc5,
    data = analysis_data,
    family = binomial()
  )

  coefficient_table <- summary(model)$coefficients
  tested_variables <- c("raw_score", "raw_score_full")
  if (!all(tested_variables %in% rownames(coefficient_table))) {
    stop("PRS coefficient was not estimable for cell type: ", cell_type)
  }

  data.frame(
    cell_type = cell_type,
    variable = c("prs_ct", "prs_full"),
    beta = coefficient_table[tested_variables, "Estimate"],
    se = coefficient_table[tested_variables, "Std. Error"],
    z = coefficient_table[tested_variables, "z value"],
    p = coefficient_table[tested_variables, "Pr(>|z|)"],
    n = nobs(model),
    row.names = NULL
  )
})

results <- bind_rows(results)
fwrite(results, output_file)
