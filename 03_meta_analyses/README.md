[← Repository overview](../README.md)

# 03. Meta-analysis with METAL

`03_meta_analyses/run_metal.sh` performs a standard-error-weighted
meta-analysis using the 2011-03-25 release of METAL. The possible cohorts are
MVP_AFR, MVP_EUR, MVP_AMR, and RACI. Supply only the cohort files appropriate
for the analysis being reproduced. Seropositive analyses should use the
corresponding seropositive-only summary statistics rather than the overall RA
summary statistics.

## Software

- Bash
- `awk`, `gzip`, and standard Unix command-line utilities
- METAL, release 2011-03-25

METAL source and precompiled binaries are available from the [official METAL
download page](https://csg.sph.umich.edu/abecasis/metal/download/). After
unpacking and compiling the source distribution, the executable is commonly
located at `generic-metal/metal`.

## Input files

Each argument is a comma-delimited REGENIE result file, optionally compressed
with gzip. The script requires these columns by name:

All cohort files supplied to the METAL workflow must use the same genome build,
and each `ID` must agree with the corresponding chromosome, position, and
alleles. Genome-build harmonization is treated as upstream preprocessing and is
not performed by the scripts in this repository. The published analysis used
GRCh37 after harmonization with RACI.

| REGENIE column | METAL role | Description |
|---|---|---|
| `ID` | `SNP` | Variant identifier in `chr:position:ref:alt` format |
| `ALLELE1` | `EFFECT` | Tested/effect allele |
| `ALLELE0` | `OTHER` | Other/reference allele |
| `BETA` | `BETA` | Log-odds effect estimate for the effect allele |
| `SE` | `SE` | Standard error of `BETA` |
| `PVAL` | `PVALUE` | Association-test p-value |
| `N` | `N` | Per-variant sample size |
| `TEST` | filtering field | Only rows with `TEST=ADD` are included |

A complete input may also contain the other REGENIE columns shown below; they
are allowed but not used directly by this script:

```text
CHROM,GENPOS,ID,ALLELE0,ALLELE1,A1FREQ,INFO,N,TEST,BETA,SE,CHISQ,LOG10P,EXTRA,PVAL,ancestry,sex
```

Missing effect estimates, standard errors, p-values, or sample sizes are
excluded during conversion. Internally, the script converts each selected
cohort to the tab-delimited schema expected by METAL:

```text
SNP EFFECT OTHER BETA PVALUE SE N
```

## Run the analysis

Set the executable and output prefix, then list the desired cohorts as
positional arguments. For an overall trans-ancestry meta-analysis using all
four available cohorts:

```bash
export METAL_BIN=/path/to/generic-metal/metal
export OUTPUT_PREFIX="$PWD/results/meta/trans_all"

bash 03_meta_analyses/run_metal.sh \
  data/MVP_AFR.csv.gz \
  data/MVP_EUR.csv.gz \
  data/MVP_AMR.csv.gz \
  data/RACI.csv.gz
```

For an MVP-only analysis, omit RACI:

```bash
export OUTPUT_PREFIX="$PWD/results/meta/mvp_trans"

bash 03_meta_analyses/run_metal.sh \
  data/MVP_AFR.csv.gz \
  data/MVP_EUR.csv.gz \
  data/MVP_AMR.csv.gz
```

For a seropositive analysis, pass the seropositive summary-statistic files:

```bash
export OUTPUT_PREFIX="$PWD/results/meta/trans_seropositive"

bash 03_meta_analyses/run_metal.sh \
  data/MVP_AFR_seropositive.csv.gz \
  data/MVP_EUR_seropositive.csv.gz \
  data/MVP_AMR_seropositive.csv.gz \
  data/RACI_seropositive.csv.gz
```

The script accepts one or more cohorts, constructs the METAL instruction file,
runs METAL, and removes only its temporary converted inputs after completion.
METAL writes the meta-analysis results using `${OUTPUT_PREFIX}_out` as the
prefix. Its principal table includes `MarkerName`, `Allele1`, `Allele2`,
`Effect`, `StdErr`, `P-value`, `Direction`, and `N`.

