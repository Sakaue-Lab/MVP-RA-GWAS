[← Repository overview](../README.md)

# 01. RA phenotype calibration

`01_phenotype/komap_calib.R` calibrates KOMAP scores against manually reviewed
RA labels, selects the smoothed score threshold closest to 90% positive
predictive value (PPV), and applies the resulting models to the full cohort.
AFR and EUR use ancestry-specific models for which only a chart label of `Yes`
is positive. The trans-ancestry model is trained in AFR, EUR, and AMR and treats
`Yes`, `Probable`, and `Possible` as positive; it is used to score AMR
participants in the final combined phenotype file.

## Software

- R (version 4.0.3)
- R packages: `data.table`, `dplyr`, `MASS`, `glmpath`, and `glmnet`

Install the required packages once in R:

```r
install.packages(c("data.table", "dplyr", "MASS", "glmpath", "glmnet"))
```

## Input files

All four inputs must contain one row per participant. `id` must identify
the same participant consistently across files.

### 1. Chart-reviewed labels (`chart_labels.csv`)

Pipe-delimited text with a header, despite the `.csv` extension.

| Column | Description |
|---|---|
| `id` | Participant identifier |
| `chartlabel` | Chart-review classification: `Yes`, `No`, `Probable`, or `Possible` |

Example:

```text
id|chartlabel
sample_001|Yes
sample_002|No
```

### 2. KOMAP scores (`komap_scores.csv`)

Pipe-delimited text with a header. It should cover the cohort to which the
calibrated phenotype will be applied.

| Column | Description |
|---|---|
| `id` | Participant identifier |
| `komap_score` | Numeric, continuous KOMAP RA score |

### 3. Covariates (`covariates.txt`)

Tab-delimited text with a header. Additional covariates may be included and are
retained in the output.

| Column | Description |
|---|---|
| `id` | Participant identifier |
| `ancestry` | Genetic-ancestry group; this analysis expects `AFR`, `EUR`, or `AMR` |

### 4. RA PheCode counts (`phecode_counts.csv`)

Comma-delimited text with a header. Additional columns are permitted.

| Column | Description |
|---|---|
| `id` | Participant identifier |
| `phe_714_1` | Number of occurrences of RA PheCode 714.1 |

## Run the analysis

From the repository root:

```bash
Rscript 01_phenotype/komap_calib.R \
  data/chart_labels.csv \
  data/komap_scores.csv \
  data/covariates.txt \
  data/phecode_counts.csv \
  results/ra_phenotype.txt
```

The script prints the selected AFR, EUR, and trans-ancestry thresholds and their
smoothed PPVs to standard error.

## Output

`ra_phenotype.txt` is tab-delimited and contains all input covariate columns plus:

| Column | Description |
|---|---|
| `komap_pred2` | Calibrated probability from the applicable ancestry model |
| `binary_ppv_90` | `1` when `komap_pred2` exceeds the threshold closest to 90% PPV; otherwise `0` |
| `phe_714_1` | RA PheCode count copied from the PheCode input |
| `phe_714_1_binary2` | `1` for at least two RA PheCode occurrences; otherwise `0` |
| `binary_ppv_90_rm` | Final phenotype; missing for one or an unknown number of RA PheCode occurrences |

`01_phenotype/library_roc.R` contains the ROC smoothing helper functions and is
loaded automatically by `komap_calib.R`; it is not run separately.

