#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    paste(
      "Usage: Rscript 01_phenotype/komap_calib.R",
      "<chart_labels.csv> <komap_scores.csv> <covariates.txt>",
      "<phecode_counts.csv> <output.txt>"
    ),
    call. = FALSE
  )
}

chart_labels_file <- args[1]
komap_scores_file <- args[2]
covariates_file <- args[3]
phecode_counts_file <- args[4]
output_file <- args[5]

# Load helpers relative to this script, not relative to the working directory.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
  "01_phenotype"
}
source(file.path(script_dir, "library_roc.R"))

require_columns <- function(dat, required, input_name) {
  missing_columns <- setdiff(required, names(dat))
  if (length(missing_columns) > 0) {
    stop(
      sprintf(
        "%s is missing required column(s): %s",
        input_name,
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

# Because the smoothed ROC curve generally has no row with PPV exactly 0.90,
# select the finite row whose PPV is closest to the requested target.
select_ppv_cutoff <- function(observed, predicted, target_ppv = 0.90) {
  roc_result <- ROC.Est.FUN(
    observed,
    predicted,
    yy0 = 0.5,
    fpr0 = seq(0, 1, 0.0001),
    yes.smooth = TRUE
  )

  roc_table <- data.frame(matrix(roc_result[-1], ncol = 6))
  names(roc_table) <- c("cut", "p.pos", "fpr", "tpr", "ppv", "npv")
  eligible <- which(is.finite(roc_table$cut) & is.finite(roc_table$ppv))
  if (length(eligible) == 0) {
    stop("No finite ROC thresholds and PPV values were produced.", call. = FALSE)
  }

  selected_index <- eligible[which.min(abs(roc_table$ppv[eligible] - target_ppv))]
  list(
    auc = unname(roc_result[1]),
    cutoff = roc_table$cut[selected_index],
    achieved_ppv = roc_table$ppv[selected_index]
  )
}

chart_labels <- read.csv(chart_labels_file, sep = "|")
komap_scores <- read.csv(komap_scores_file, sep = "|")
covariates <- read.table(covariates_file, sep = "\t", header = TRUE)
phecode_counts <- fread(phecode_counts_file)

require_columns(chart_labels, c("id", "chartlabel"), "chart-label file")
require_columns(komap_scores, c("id", "komap_score"), "KOMAP-score file")
require_columns(covariates, c("id", "ancestry"), "covariate file")
require_columns(phecode_counts, c("id", "phe_714_1"), "PheCode-count file")

# Prevent joins from failing when an identifier is inferred as numeric in one
# file and character in another.
chart_labels$id <- as.character(chart_labels$id)
komap_scores$id <- as.character(komap_scores$id)
covariates$id <- as.character(covariates$id)
phecode_counts$id <- as.character(phecode_counts$id)

allowed_labels <- c("Yes", "No", "Probable", "Possible")
unexpected_labels <- setdiff(unique(na.omit(chart_labels$chartlabel)), allowed_labels)
if (length(unexpected_labels) > 0) {
  warning("Unexpected chartlabel value(s): ", paste(unexpected_labels, collapse = ", "))
}

calibration_data <- chart_labels %>%
  left_join(select(covariates, id, ancestry), by = "id") %>%
  left_join(select(komap_scores, id, komap_score), by = "id") %>%
  filter(!is.na(komap_score), !is.na(ancestry), !is.na(chartlabel)) %>%
  mutate(
    # Ancestry-specific calibration treats only chart-confirmed "Yes" as RA.
    chartlabel_strict = as.integer(chartlabel == "Yes"),
    # Trans-ancestry calibration also treats Probable/Possible as RA.
    chartlabel_trans = as.integer(chartlabel != "No")
  )

calibration_afr <- filter(calibration_data, ancestry == "AFR")
calibration_eur <- filter(calibration_data, ancestry == "EUR")
calibration_trans <- filter(calibration_data, ancestry %in% c("AFR", "EUR", "AMR"))

fit_afr <- glm(chartlabel_strict ~ komap_score, family = "binomial", data = calibration_afr)
fit_eur <- glm(chartlabel_strict ~ komap_score, family = "binomial", data = calibration_eur)
fit_trans <- glm(chartlabel_trans ~ komap_score, family = "binomial", data = calibration_trans)

calibration_afr$komap_probability <- predict(fit_afr, type = "response")
calibration_eur$komap_probability <- predict(fit_eur, type = "response")
calibration_trans$komap_probability <- predict(fit_trans, type = "response")

roc_afr <- select_ppv_cutoff(calibration_afr$chartlabel_strict, calibration_afr$komap_probability)
roc_eur <- select_ppv_cutoff(calibration_eur$chartlabel_strict, calibration_eur$komap_probability)
roc_trans <- select_ppv_cutoff(calibration_trans$chartlabel_trans, calibration_trans$komap_probability)

message(sprintf("AFR cutoff: %.6f (smoothed PPV %.4f)", roc_afr$cutoff, roc_afr$achieved_ppv))
message(sprintf("EUR cutoff: %.6f (smoothed PPV %.4f)", roc_eur$cutoff, roc_eur$achieved_ppv))
message(sprintf("Trans-ancestry cutoff: %.6f (smoothed PPV %.4f)", roc_trans$cutoff, roc_trans$achieved_ppv))

# Use ancestry-specific models for AFR/EUR and the trans-ancestry model for AMR.
# This avoids duplicate participant rows in the final combined phenotype file.
scores_with_ancestry <- komap_scores %>%
  left_join(select(covariates, id, ancestry), by = "id")

scored_afr <- scores_with_ancestry %>%
  filter(ancestry == "AFR") %>%
  mutate(
    komap_pred2 = predict(fit_afr, newdata = ., type = "response"),
    binary_ppv_90 = as.integer(komap_pred2 > roc_afr$cutoff)
  )

scored_eur <- scores_with_ancestry %>%
  filter(ancestry == "EUR") %>%
  mutate(
    komap_pred2 = predict(fit_eur, newdata = ., type = "response"),
    binary_ppv_90 = as.integer(komap_pred2 > roc_eur$cutoff)
  )

scored_amr <- scores_with_ancestry %>%
  filter(ancestry == "AMR") %>%
  mutate(
    komap_pred2 = predict(fit_trans, newdata = ., type = "response"),
    binary_ppv_90 = as.integer(komap_pred2 > roc_trans$cutoff)
  )

scored_participants <- bind_rows(scored_afr, scored_eur, scored_amr) %>%
  select(id, komap_pred2, binary_ppv_90) %>%
  mutate(id = as.character(id))

phecode_counts <- phecode_counts %>%
  select(id, phe_714_1) %>%
  mutate(id = as.character(id))

output <- covariates %>%
  mutate(id = as.character(id)) %>%
  left_join(scored_participants, by = "id") %>%
  left_join(phecode_counts, by = "id") %>%
  mutate(
    phe_714_1_binary2 = as.integer(phe_714_1 >= 2),
    # Exclude participants with exactly one RA PheCode or a missing count.
    binary_ppv_90_rm = if_else(
      phe_714_1 == 1 | is.na(phe_714_1),
      NA_integer_,
      coalesce(binary_ppv_90, 0L)
    )
  )

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(output, output_file, row.names = FALSE, sep = "\t", na = "NA")
