#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(boot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5 || length(args) > 6) {
  stop(
    paste(
      "Usage: Rscript 09_prs/prs_calc.R <score_dir> <training_ids.txt>",
      "<covariates.txt> <ancestry> <output.csv> [bootstrap_replicates]"
    ),
    call. = FALSE
  )
}

score_dir <- args[1]
training_file <- args[2]
covariate_file <- args[3]
analysis_ancestry <- args[4]
output_file <- args[5]
bootstrap_replicates <- if (length(args) == 6) as.integer(args[6]) else 1000L

score_files <- list.files(
  score_dir,
  pattern = "\\.sscore\\.gz$",
  full.names = TRUE,
  recursive = TRUE
)
if (length(score_files) == 0) stop("No .sscore.gz files found in ", score_dir)

read_score <- function(path) {
  dat <- fread(path)
  required <- c("#IID", "SCORE1_SUM")
  missing_columns <- setdiff(required, names(dat))
  if (length(missing_columns) > 0) {
    stop(path, " is missing column(s): ", paste(missing_columns, collapse = ", "))
  }
  dat %>% transmute(id = as.character(`#IID`), chunk_score = SCORE1_SUM)
}

# Sum chromosome/chunk score contributions by participant using IDs, without
# assuming that every score file has the same row order.
scores <- bind_rows(lapply(score_files, read_score)) %>%
  group_by(id) %>%
  summarize(raw_score = sum(chunk_score), .groups = "drop")

training_ids <- fread(training_file)
if (!("IID" %in% names(training_ids))) stop("Training-ID file must contain IID.")
training_ids$IID <- as.character(training_ids$IID)

covariates <- fread(covariate_file)
required_covariates <- c(
  "id", "ancestry", "age", "sex", paste0("pc", 1:5),
  "binary_ppv_90_rm_trans"
)
missing_covariates <- setdiff(required_covariates, names(covariates))
if (length(missing_covariates) > 0) {
  stop("Covariate file is missing: ", paste(missing_covariates, collapse = ", "))
}

evaluation_data <- covariates %>%
  mutate(id = as.character(id)) %>%
  filter(
    ancestry == analysis_ancestry,
    !(id %in% training_ids$IID)
  ) %>%
  inner_join(scores, by = "id")

null_formula <- binary_ppv_90_rm_trans ~ sex + age + pc1 + pc2 + pc3 + pc4 + pc5
full_formula <- update(null_formula, . ~ . + raw_score)
null_model <- glm(null_formula, data = evaluation_data, family = binomial())
full_model <- glm(full_formula, data = evaluation_data, family = binomial())

nagelkerke_r2 <- function(full_model, null_model) {
  n <- nobs(full_model)
  ll_full <- as.numeric(logLik(full_model))
  ll_null <- as.numeric(logLik(null_model))
  r2_cs <- 1 - exp((2 / n) * (ll_null - ll_full))
  r2_cs / (1 - exp((2 / n) * ll_null))
}

point_r2 <- nagelkerke_r2(full_model, null_model)
bootstrap_statistic <- function(dat, indices) {
  sampled <- dat[indices, , drop = FALSE]
  tryCatch({
    sampled_null <- glm(null_formula, data = sampled, family = binomial())
    sampled_full <- glm(full_formula, data = sampled, family = binomial())
    nagelkerke_r2(sampled_full, sampled_null)
  }, error = function(e) NA_real_)
}

set.seed(123)
bootstrap <- boot(
  data = evaluation_data,
  statistic = bootstrap_statistic,
  R = bootstrap_replicates
)
valid_bootstrap <- bootstrap$t[is.finite(bootstrap$t)]
confidence_interval <- if (length(valid_bootstrap) > 0) {
  quantile(valid_bootstrap, c(0.025, 0.975), na.rm = TRUE)
} else {
  c(NA_real_, NA_real_)
}

prs_coefficient <- summary(full_model)$coefficients["raw_score", ]
result <- data.frame(
  ancestry = analysis_ancestry,
  n = nobs(full_model),
  n_cases = sum(model.response(model.frame(full_model)) == 1),
  prs_beta = prs_coefficient["Estimate"],
  prs_se = prs_coefficient["Std. Error"],
  prs_pvalue = prs_coefficient["Pr(>|z|)"],
  nagelkerke_r2 = point_r2,
  nagelkerke_r2_lower = confidence_interval[1],
  nagelkerke_r2_upper = confidence_interval[2],
  bootstrap_replicates = bootstrap_replicates
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(result, output_file)
