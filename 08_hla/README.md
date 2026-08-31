[← Repository overview](../README.md)

# 08. HLA analyses

The HLA workflow performs sequential amino-acid omnibus tests and estimates
effects for selected HLA haplotypes in MVP_AFR and MVP_EUR. The upstream HLA
imputation software is not prescribed by this repository; users must supply
the dosage, variant-metadata, and amino-acid dictionary files described below.

The scripts are ordered as follows:

1. `cond_haplo_kern.R`: initial, unconditioned omnibus tests across HLA genes;
2. `cond_haplo_kern_round2.R`: conditional tests within HLA-DRB1;
3. `cond_haplo_kern_round3.R`: tests of additional HLA genes conditional on
   selected HLA-DRB1 and HLA-B positions;
4. `cond_haplo_kern_round4.R`: a subsequent conditional omnibus-test round;
5. `est_haplo_effect.R`: effect estimation for selected haplotypes.

Conditioning positions were selected sequentially from the preceding round,
rather than specified completely before the first analysis.

The haplotype-effect figures and tables used haplotypes defined by the
EUR-identified HLA-DRB1 positions 11, 13, 71, and 74. The same EUR-defined
haplotypes were then evaluated in both EUR and AFR to enable a direct
cross-ancestry comparison. This effect-estimation step is distinct from the
AFR conditional analysis, in which position 85 was the ancestry-specific
secondary HLA-DRB1 signal.

The reported haplotype-effect models additionally adjust for HLA-B residue
`9:D` and HLA-DPB1 residue `9:F`, which were significant in the conditional
testing. An exploratory conditional step involving HLA-DQB1 position 37 is not
part of the documented final workflow.

## Shared input files

### HLA amino-acid dictionary

Comma-delimited text with one row for each HLA allele/amino-acid mapping:

| Column | Description |
|---|---|
| `gene` | HLA gene, such as `A`, `B`, or `DRB1` |
| `hla` | Imputed HLA allele represented by the row |
| `pos` | Amino-acid position; negative values denote positions before the mature protein sequence |
| `AA` | Amino-acid residue at that position for the indicated HLA allele |

### HLA dosage file

PLINK-style `.raw` text containing `FID`, `IID`, `PAT`, `MAT`, `SEX`, and
`PHENOTYPE`, followed by HLA allele-dosage columns. Dosage-column names must
correspond to the alleles represented in the HLA dictionary after PLINK's name
sanitization. The file contains patient-level data and is not distributed.

### Variant metadata (`.pvar`)

Tab-delimited PLINK 2 variant metadata containing:

```text
#CHROM POS ID REF ALT FILTER INFO
```

For imputed HLA variants, `INFO` must contain semicolon-delimited `AF`, `MAF`,
and `R2` entries. The conditional testing applies post-imputation thresholds of
MAF greater than 0.01 and R2 greater than 0.7.

### Covariate/phenotype file

Tab-delimited participant-level data containing:

```text
id ancestry binary_ppv_90_rm_trans sex age pc1 pc2 pc3 pc4 pc5
```

Only participants with `ancestry` equal to `AFR` or `EUR` are analyzed. The
logistic models use `binary_ppv_90_rm_trans` as the outcome and adjust for sex,
age, and the first five within-ancestry genetic principal components.

### Related-participant exclusion file

Comma-delimited text with one column named `IID`. Participants listed here are
removed before conditional testing.

The previously loaded MAF `.RDS` object was unused and is no longer required.

## Execution

All five scripts use the same positional arguments:

```text
1. HLA amino-acid dictionary CSV
2. phenotype/covariate TXT
3. HLA dosage RAW
4. related-participant exclusion CSV
5. variant-metadata PVAR
6. output CSV
```

Run the initial omnibus tests first:

```bash
Rscript 08_hla/cond_haplo_kern.R \
  data/hla_dictionary.csv \
  data/hla_covariates.txt \
  data/hla_dosages.raw \
  data/related_iids.csv \
  data/hla_variants.pvar \
  results/hla/round1_omnibus.csv
```

The subsequent scripts represent sequential conditional rounds. Conditioning
positions are encoded in each script and were selected after review of the
preceding round's results:

```bash
Rscript 08_hla/cond_haplo_kern_round2.R \
  data/hla_dictionary.csv data/hla_covariates.txt data/hla_dosages.raw \
  data/related_iids.csv data/hla_variants.pvar \
  results/hla/round2_conditional.csv

Rscript 08_hla/cond_haplo_kern_round3.R \
  data/hla_dictionary.csv data/hla_covariates.txt data/hla_dosages.raw \
  data/related_iids.csv data/hla_variants.pvar \
  results/hla/round3_conditional.csv

Rscript 08_hla/cond_haplo_kern_round4.R \
  data/hla_dictionary.csv data/hla_covariates.txt data/hla_dosages.raw \
  data/related_iids.csv data/hla_variants.pvar \
  results/hla/round4_conditional.csv
```

Each conditional-test output contains the ancestry, HLA gene, amino-acid
position, likelihood-ratio-test deviance, and omnibus p-value. Rounds 1 and 2
evaluate both EUR and AFR; the later ancestry-specific rounds implement the AFR
conditional analysis. The scripts print the current ancestry, gene, and
position to show progress.

Finally, estimate the EUR-defined haplotype effects in EUR and AFR:

```bash
Rscript 08_hla/est_haplo_effect.R \
  data/hla_dictionary.csv \
  data/hla_covariates.txt \
  data/hla_dosages.raw \
  data/related_iids.csv \
  data/hla_variants.pvar \
  results/hla/eur_defined_haplotype_effects.csv
```

The effect output contains ancestry, haplotype, log-odds estimate, odds ratio,
standard error, p-value, haplotype frequency, and McFadden pseudo-R-squared.
The same relatedness exclusion and post-imputation MAF/R2 thresholds are applied
in every round and in haplotype-effect estimation.

