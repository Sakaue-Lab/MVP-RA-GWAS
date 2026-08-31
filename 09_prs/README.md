[← Repository overview](../README.md)

# 09. Polygenic risk scores

The MVP PGEN data are stored in multiple chunks per chromosome. To make LD
clumping computationally tractable, 1,000 randomly selected EUR participants
were retained in each chunk and the downsampled chunks were merged by
chromosome. These merged datasets are used only as LD references. PRSs are then
calculated in the full cohort from the original, un-downsampled PGEN chunks.

The complete order is:

1. `prs_prep.R`: prepare METAL summary statistics and remove the MHC;
2. `downsample.sh`: retain the same 1,000 EUR participants in every PGEN chunk;
3. `merge_pgen.sh`: merge the downsampled chunks within each chromosome;
4. `clump_merged.sh`: LD-clump at a selected p-value threshold;
5. `prepare_clumped_weights.R`: convert `.clumps` outputs into scoring weights;
6. `prs_calc.sh`: calculate partial scores in each full-cohort PGEN chunk; and
7. `prs_calc.R`: sum chunk scores by participant and evaluate held-out R2.

## Software and input conventions

- R packages: `data.table`, `dplyr`, and `boot`
- PLINK 2, available as `plink2` or supplied as the final shell-script argument
- PGEN chunks arranged as `PGEN_DIR/chr1/*.pgen`, ..., `PGEN_DIR/chr22/*.pgen`
- matching `.pvar` and `.psam` files for every `.pgen`
- variant IDs consistent between the METAL output and PGEN files

The downsampling keep file is a PLINK participant file with an `FID IID` header
and 1,000 randomly selected EUR participants. The same file must be used for
every chromosome and chunk. Creation of the random sample is treated as cohort
preprocessing and is not performed by these scripts.

## 1. Prepare summary statistics

The input is a METAL output containing `MarkerName`, `Effect`, `StdErr`,
`P-value`, `Allele1`, and `Allele2`. `MarkerName` must use
`chr:position:ref:alt` format. The script removes `chr6:25,000,000-35,000,000`
before clumping.

```bash
Rscript 09_prs/prs_prep.R \
  results/meta/trans_all_out1.txt \
  results/prs/snp_list_mhc_removed.csv \
  results/prs/gwas_weights_mhc_removed.csv
```

The GWAS output contains `ID`, `Beta`, `P`, `SE`, `A1`, and `A2`, where `A1`
is the METAL effect allele corresponding to `Beta`.

## 2. Downsample each chromosome's chunks

```bash
for chromosome in $(seq 1 22); do
  bash 09_prs/downsample.sh \
    "${chromosome}" \
    /path/to/full_chunked_pgen \
    data/clump_eur_1000_keep.txt \
    results/prs/downsampled \
    /path/to/plink2
done
```

This uses PLINK `--keep`, which retains participants rather than variants.

## 3. Merge the downsampled chunks by chromosome

```bash
for chromosome in $(seq 1 22); do
  bash 09_prs/merge_pgen.sh \
    "${chromosome}" \
    results/prs/downsampled \
    results/prs/merged_ld_reference \
    /path/to/plink2
done
```

The resulting prefixes are
`results/prs/merged_ld_reference/chrN/merged_chrN`.

## 4. Clump across a range of p-value thresholds

The study evaluated thresholds from 0.1 through `5e-8`, then selected the
best-performing threshold using Nagelkerke R2 in held-out participants. An
example threshold grid is:

```bash
thresholds=(
  1e-1 3e-2 1e-2 3e-3 1e-3 3e-4 1e-4
  3e-5 1e-5 3e-6 1e-6 3e-7 1e-7 5e-8
)

for threshold in "${thresholds[@]}"; do
  bash 09_prs/clump_merged.sh \
    results/prs/merged_ld_reference \
    results/prs/gwas_weights_mhc_removed.csv \
    "${threshold}" \
    results/prs/clumps \
    /path/to/plink2
done
```

Clumping uses an R2 threshold of 0.2 and the selected p-value as PLINK's
`--clump-p1` threshold.

## 5. Create per-chromosome scoring weights

This is the previously missing handoff between PLINK `.clumps` output and PRS
calculation. Run it once for each threshold:

```bash
for threshold in "${thresholds[@]}"; do
  Rscript 09_prs/prepare_clumped_weights.R \
    results/prs/gwas_weights_mhc_removed.csv \
    "results/prs/clumps/${threshold}" \
    "results/prs/scoring_weights/${threshold}"
done
```

Each threshold directory contains `1.gwasres.csv.gz` through
`22.gwasres.csv.gz`, with columns `ID`, `Beta`, `P`, and `A1`.

## 6. Calculate PRSs in the full chunked cohort

For every threshold and chromosome, score the original full-cohort PGEN chunks:

```bash
for threshold in "${thresholds[@]}"; do
  for chromosome in $(seq 1 22); do
    bash 09_prs/prs_calc.sh \
      "${chromosome}" \
      /path/to/full_chunked_pgen \
      "results/prs/scoring_weights/${threshold}/${chromosome}.gwasres.csv.gz" \
      "results/prs/scores/${threshold}" \
      /path/to/plink2
  done
done
```

Chunks without selected variants are skipped. Other chunks produce compressed
PLINK `.sscore.gz` files containing `#IID` and `SCORE1_SUM`.

## 7. Combine scores and evaluate held-out performance

The training-ID file must contain an `IID` column. The covariate file must
contain `id`, `ancestry`, `binary_ppv_90_rm_trans`, `age`, `sex`, and
`pc1`-`pc5`. Participants listed in the training file are excluded from
evaluation.

```bash
for threshold in "${thresholds[@]}"; do
  Rscript 09_prs/prs_calc.R \
    "results/prs/scores/${threshold}" \
    data/eur_training_ids.txt \
    data/covariates.txt \
    EUR \
    "results/prs/performance/${threshold}.csv" \
    1000
done
```

The script sums `SCORE1_SUM` across chunks by `#IID`, fits covariate-only and
PRS-plus-covariate logistic models, and reports the PRS coefficient, standard
error, p-value, Nagelkerke R2, and a percentile bootstrap 95% interval. Compare
the resulting performance files to select the best threshold.

