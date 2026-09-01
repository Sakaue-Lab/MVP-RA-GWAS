[← Repository overview](../README.md)

# BUHMBOX and cell-type PRS seropositivity analysis

This subfolder contains the following scripts:

1. `run_buhmbox.sh` tests whether a set of cell-type-associated disease-risk
   alleles shows excessive cross-locus correlation within cases, consistent
   with genetic heterogeneity.
2. `buhmbox_post.R` tests whether a cell-type-specific PRS is associated with
   seropositivity in held-out participants after adjustment for the full RA PRS
   and clinical/genetic covariates.

The second analysis does **not** use the BUHMBOX `.BBrst` statistic as a PRS
weight. Cell-type-specific PRSs must be constructed separately from the selected
cell-type variant sets and GWAS effect sizes before running `buhmbox_post.R`.

See the [official BUHMBOX page](https://software.broadinstitute.org/mpg/buhmbox/)
and [BUHMBOX manual](https://software.broadinstitute.org/mpg/buhmbox/buhmbox_manual.html)
for the method, original R implementation, argument definitions, and output
interpretation.

## Software

- Bash and `awk`
- PLINK 2
- R with `data.table` and `dplyr`
- The official `buhmbox.R` program downloaded from the BUHMBOX website

Make the shell script executable:

```bash
chmod +x 11_buhmbox_ctPRS/run_buhmbox.sh
```

## Part 1: run BUHMBOX

### Inputs

#### PGEN dataset

Provide the shared prefix for matching `.pgen`, `.pvar`, and `.psam` files.
Variant IDs must match the IDs in the cell-type SNP BED file.

#### Case-control phenotype file

Whitespace-delimited, without a header, with three columns:

| Column | Description |
|---|---|
| 1 | `FID` |
| 2 | `IID` |
| 3 | PLINK phenotype status: `1` control, `2` case |

```text
sample_001 sample_001 2
sample_002 sample_002 1
```

The script uses this file to create separate case and control keep files. Note that we kept only a random sample of controls: same as the number of cases, sampled from EUR controls. 

#### Cell-type SNP BED file

Headerless, tab-delimited, with at least seven columns:

| Column | Name | Description |
|---|---|---|
| 1 | `CHROM` | Chromosome |
| 2 | `START` | Zero-based interval start |
| 3 | `END` | Interval end |
| 4 | `ID` | Variant ID matching the PGEN data |
| 5 | `EFFECT_ALLELE` | GWAS effect/risk allele |
| 6 | `BETA` | GWAS log-odds effect estimate |
| 7 | `OR` | GWAS odds ratio |

```text
chr1 12344 12345 chr1:12345:A:G G 0.1823 1.20
```

The script converts this to BUHMBOX's four-column SNP format: variant ID,
effect/risk allele, allele frequency, and odds ratio. Because the BED input does
not include control allele frequency, the frequency is written as `NA`; BUHMBOX
then imputes it internally. Odds ratios below 1 are permitted by BUHMBOX v0.37
and later and are flipped internally.

#### Principal-component covariate file

Whitespace-delimited with a header and the following column order:

```text
FID IID sex age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
```

The script writes `FID`, `IID`, and `PC1`–`PC10` to the BUHMBOX PC file. Use `-`
instead of a covariate path to omit PC adjustment, although PC adjustment is
recommended. IDs are matched by `FID` and `IID`.

### Execute

The arguments are:

```text
run_buhmbox.sh CELL_TYPE PGEN_PREFIX PHENOTYPE SNP_BED COVARIATES OUTPUT_DIR [PLINK2] [BUHMBOX_R]
```

Example:

```bash
bash 11_buhmbox_ctPRS/run_buhmbox.sh \
  Tcell \
  data/genotypes/analysis_cohort \
  data/phenotypes/case_control.txt \
  data/cell_type_snps/Tcell_variants.bed \
  data/covariates.txt \
  results/buhmbox \
  /path/to/plink2 \
  /path/to/buhmbox.R
```

Run the command once per cell type. For multiple cell types:

```bash
for cell_type in Tcell Bcell Monocyte; do
  bash 11_buhmbox_ctPRS/run_buhmbox.sh \
    "${cell_type}" \
    data/genotypes/analysis_cohort \
    data/phenotypes/case_control.txt \
    "data/cell_type_snps/${cell_type}_variants.bed" \
    data/covariates.txt \
    results/buhmbox \
    /path/to/plink2 \
    /path/to/buhmbox.R
done
```

### Operations performed by the script

1. Convert the cell-type BED file to BUHMBOX SNP format.
2. Split the phenotype file into case and control keep files.
3. Estimate LD in controls and prune at `r2 < 0.1` using a 1 Mb window.
4. Export risk-allele dosages separately for cases and controls in PLINK `.raw`
   format.
5. Extract `PC1`–`PC10` when a covariate file is supplied.
6. Run BUHMBOX with arguments `YY Y Y`: use frequency and OR information,
   estimate the subgroup proportion, and run the risk-score/MR analysis.

### Outputs

For cell type `Tcell`, outputs are written under `OUTPUT_DIR/Tcell/` with prefix
`Tcell_buhmbox`:

- `Tcell_buhmbox.BBrst`: primary BUHMBOX result table;
- `Tcell_buhmbox.log`: BUHMBOX progress log; and
- `tmp/`: SNP lists, case/control lists, LD-pruning files, PC file, and exported
  dosage `.raw` files.

The `.BBrst` file can include the BUHMBOX p-value and Z-score, estimated mixture
proportions and likelihoods, and risk-score/MR coefficient statistics. A
significant BUHMBOX p-value indicates excessive positive correlation among the
tested loci in cases; it is evidence consistent with heterogeneity but is not,
by itself, proof of a specific clinical subgroup.

## Part 2: associate cell-type PRSs with seropositivity

`buhmbox_post.R` evaluates each cell-type PRS in held-out participants. It fits
the following logistic model separately for every cell type:

```text
seropositive ~ cell-type PRS + full RA PRS + sex + age + pc1 + ... + pc5
```

### Inputs

#### Cell-type list

Headerless text with one cell-type label per line. Labels determine both the
subdirectory and score filename:

```text
Tcell
Bcell
Monocyte
```

For `Tcell`, the script reads:

```text
CELL_SCORE_DIR/Tcell/scores_Tcell.csv
```

#### Cell-type PRS score files

Each score file must contain:

| Column | Description |
|---|---|
| `id` | Participant identifier |
| `raw_score` | Cell-type-specific PRS |

These scores must be created before this step. The general PRS workflow in the
parent folder can be adapted by restricting its weights to each cell-type SNP
set while retaining GWAS-derived effect sizes.

#### Full RA PRS file

| Column | Description |
|---|---|
| `id` | Participant identifier |
| `raw_score` | Full RA PRS, renamed internally to `raw_score_full` |

#### Training-ID file

Must contain an `IID` column. These participants are excluded so association
testing is performed in held-out samples.

#### Covariate file

Must contain `id`, `sex`, `age`, and `pc1`–`pc5`. The script does not filter an
ancestry column; provide a covariate file already restricted to the intended
analysis cohort if ancestry-specific evaluation is required.

#### Serology file

| Column | Description |
|---|---|
| `id` | Participant identifier |
| `ccflab` | Anti-CCP laboratory value; may be missing |
| `rflab` | Rheumatoid factor laboratory value; may be missing |
| `seropositive` | Binary seropositivity indicator |

Participants are included if either `ccflab` or `rflab` is observed. As in the
original analysis, a missing `seropositive` value among participants with an
observed serology laboratory value is assigned `0`.

### Execute

```bash
Rscript 11_buhmbox_ctPRS/buhmbox_post.R \
  data/cell_types.txt \
  results/cell_type_prs \
  results/prs/full_ra_scores.csv \
  data/eur_training_ids.txt \
  data/held_out_covariates.txt \
  data/serology.txt \
  results/buhmbox/prs_seropositivity_results.csv
```

### Output

The output CSV contains two rows per cell type:

| Column | Description |
|---|---|
| `cell_type` | Cell-type label |
| `variable` | `prs_ct` for the cell-type PRS or `prs_full` for the full RA PRS |
| `beta` | Logistic-regression coefficient |
| `se` | Standard error |
| `z` | Wald Z statistic |
| `p` | Two-sided Wald p-value |
| `n` | Participants included in the fitted model |

The `prs_ct` row is the primary test of whether the cell-type PRS is associated
with seropositivity independently of the full RA PRS and covariates.
