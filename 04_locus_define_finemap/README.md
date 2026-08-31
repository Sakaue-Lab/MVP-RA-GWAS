[← Repository overview](../README.md)

# 04. Lead-locus definition and fine-mapping

## Define lead loci

`04_locus_define_finemap/define_loci.R` identifies genome-wide-significant lead
variants from four METAL analyses:

1. overall trans-ancestry: MVP_AFR, MVP_EUR, MVP_AMR, and RACI;
2. overall European-ancestry: MVP_EUR and RACI;
3. seropositive trans-ancestry; and
4. seropositive European-ancestry.

All variant positions and distance calculations must use one consistent genome
build. The published analysis used GRCh37.
Loci are defined by physical distance rather than LD. Starting with the most
significant remaining variant, the script retains that variant and removes all
other variants within +/- 1 Mb on the same chromosome. It excludes the extended
MHC region (`chr6:25,000,000-35,000,000`) and chromosome X.

Overall trans-ancestry signals deliberately receive priority: lead variants are
first selected within the overall trans-ancestry analysis and placed ahead of
signals from the other three analyses. Thus, when signals from multiple
analyses fall within the same 1 Mb window, the overall trans-ancestry signal
defines the locus even if a secondary analysis has a larger test statistic.

### Software and inputs

- R (version 4.1 or newer recommended)
- R packages: `data.table` and `dplyr`
- four METAL output tables from the preceding workflow

Each METAL table must contain:

| Column | Description |
|---|---|
| `MarkerName` | Variant identifier in `chr:position:ref:alt` format, using the common analysis build |
| `Effect` | Meta-analysis log-odds effect estimate |
| `StdErr` | Standard error of `Effect` |
| `P-value` | Meta-analysis p-value |
| `Direction` | Cohort-specific effect directions from METAL |
| `Allele1` | Effect allele |
| `Allele2` | Other allele |

### Run the analysis

```bash
Rscript 04_locus_define_finemap/define_loci.R \
  results/meta/trans_all_out1.txt \
  results/meta/eur_all_out1.txt \
  results/meta/trans_seropositive_out1.txt \
  results/meta/eur_seropositive_out1.txt \
  results/loci/sig_alleles.csv
```

The exact numeric suffix in METAL filenames, such as `out1`, is assigned by
METAL; use the paths produced by the preceding commands.

The output is a comma-delimited table containing the lead variant's chromosome,
position, reference and alternate alleles, effect estimate, standard error,
p-value, direction, METAL alleles, chi-square statistic, and the analysis in
which the retained signal was identified (`ancestry`).

## Add ancestry-specific MVP statistics

`04_locus_define_finemap/annotate_ancestry_stats.R` annotates each lead variant
with allele frequency, minor allele frequency, imputation INFO score, and sample
size from the ancestry-specific MVP REGENIE results. These annotations come
from MVP_EUR, MVP_AFR, and MVP_AMR—not from METAL output, which does not retain
all of these cohort-level fields.

The lead-variant file and all three REGENIE annotation inputs must use the same
genome build. Coordinate conversion itself is outside the scope of this
repository. In each harmonized file, `GENPOS` is the position in the shared
analysis build. No filtering on `TEST` is performed, but each input must contain
at most one row for each combination of chromosome, position, reference allele,
and alternate allele.

### Required REGENIE columns

| Column | Description |
|---|---|
| `CHROM` | Chromosome |
| `GENPOS` | Position in the common analysis genome build |
| `ALLELE0` | Reference allele |
| `ALLELE1` | Alternate/tested allele |
| `A1FREQ` | Frequency of `ALLELE1` |
| `INFO` | Imputation INFO score |
| `N` | Per-variant sample size |

Both uncompressed `.csv` and gzip-compressed `.csv.gz` files are accepted by
`data.table::fread()`.

### Run the annotation

```bash
Rscript 04_locus_define_finemap/annotate_ancestry_stats.R \
  results/loci/sig_alleles.csv \
  data/MVP_EUR_harmonized.csv.gz \
  data/MVP_AFR_harmonized.csv.gz \
  data/MVP_AMR_harmonized.csv.gz \
  results/loci/sig_alleles_ancestry_annotated.csv
```

The script joins variants using chromosome, position in the common build,
reference allele, and alternate allele. It adds `EUR_AltFreq`, `EUR_MAF`,
`EUR_INFO`, and `EUR_N`, with corresponding fields for AFR and AMR. Minor allele
frequency is calculated as `min(A1FREQ, 1 - A1FREQ)`.

## Compare lead variants with a previous GWAS

`04_locus_define_finemap/annotate_previous_study.R` identifies the nearest
previously reported locus on the same chromosome and classifies a current lead
variant as novel when it is more than 1 Mb from that locus. The published
analysis used Ishigaki et al. because it was the most recent major RA GWAS at
the time, but the script is intentionally study-agnostic and can use another
appropriate reference GWAS.

The previous-study file should contain the previously reported significant or
lead variants—not every variant tested genome-wide—when the goal is to assess
locus novelty. The current and previous-study files must use the same genome
build. To compare against multiple earlier GWAS, combine their reported lead
variants into this one input table and, if useful, include a column identifying
the source study. A separate file or annotation pass is not required for each
study.

The script accepts a variant-ID column named any one of:

```text
VariantID
Variant.ID
ID
MarkerName
```

IDs must use `chr:position:ref:alt` format. All other columns are optional and
are carried into the output with a `previous_nearest_` prefix. Thus, users may
include a previous study's locus name, rsID, effect estimate, confidence
interval, p-value, population, serostatus, or allele frequencies without
changing the script.

```bash
Rscript 04_locus_define_finemap/annotate_previous_study.R \
  results/loci/sig_alleles_ancestry_annotated.csv \
  data/previous_ra_gwas_lead_variants.csv \
  results/loci/sig_alleles_previous_study_annotated.csv
```

The output adds parsed coordinates and metadata for the nearest previous-study
locus, `distance_to_previous_locus`, and the logical indicator
`novel_vs_previous_study`. A variant with no previous-study locus on the same
chromosome is also classified as novel.

## Construct 95% credible sets

`04_locus_define_finemap/fine_map.R` performs approximate-Bayes-factor
fine-mapping for one meta-analysis at a time. It is not tied to a particular
cohort combination: run it separately with the METAL output and corresponding
lead variants from any overall, ancestry-specific, or seropositive
meta-analysis.

For each lead variant, the script:

1. includes every variant within +/- 500 kb on the same chromosome, without
   additional MAF, INFO, or sample-size filtering;
2. computes Wakefield approximate Bayes factors using prior variance
   `w = 0.04`;
3. normalizes the Bayes factors to posterior inclusion probabilities (PIPs);
4. orders variants by decreasing PIP; and
5. retains the smallest set whose cumulative PIP reaches at least 95%.

The METAL file must contain `MarkerName`, `Effect`, `StdErr`, `P-value`,
`Allele1`, and `Allele2`. `MarkerName` must use `chr:position:ref:alt` format.
The lead-variant file must contain `Chr`, `Pos`, `RefAllele`, and `AltAllele`,
and it must use the same genome build as the METAL file.

```bash
Rscript 04_locus_define_finemap/fine_map.R \
  results/meta/trans_all_out1.txt \
  results/loci/trans_lead_variants.csv \
  results/finemapping/trans_95pct_credible_sets.csv
```

Repeat the command with another matched METAL/lead-variant pair to fine-map a
different meta-analysis. The output contains the variants in each 95% credible
set, their effect statistics, log approximate Bayes factors, PIPs, cumulative
PIPs, credible-set ranks and sizes, and identifying information for the lead
variant and locus.

