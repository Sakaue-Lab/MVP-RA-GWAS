[← Repository overview](../README.md)

# 02. MVP genome-wide association analysis

The two scripts in `02_mvp_gwas/` run a binary-trait GWAS with REGENIE. Step 1
fits the null model using directly genotyped variants. Step 2 tests imputed
variants stored in chromosome- and chunk-specific BGEN files. The scripts retain
the historical software versions used for the analysis: REGENIE 3.1.3 for Step
1 and REGENIE 2.2.4 for Step 2.

The shell scripts are configured through exported environment variables so that
restricted filesystem paths do not need to be hard-coded in the repository.
They contain LSF `#BSUB` directives and were designed for submission with
`bsub`. Queue, project, memory, and log settings may require adaptation outside
the original MVP computing environment.

## Input files shared by both steps

### Phenotype file

Tab- or space-delimited text with a header and one row per participant:

| Column | Description |
|---|---|
| `FID` | Family identifier |
| `IID` | Individual identifier |
| `binary_ppv_90_rm` | Binary RA phenotype: `1` case, `0` control, missing for excluded participants |

```text
FID IID binary_ppv_90_rm
sample_001 sample_001 1
sample_002 sample_002 0
```

### Covariate file

Tab- or space-delimited text with a header and one row per participant:

| Column | Description |
|---|---|
| `FID` | Family identifier |
| `IID` | Individual identifier |
| `age` | Numeric age at the defined index date |
| `sex` | Numerically encoded sex covariate, using one consistent coding scheme |
| `pc1`–`pc5` | First five within-ancestry genetic principal components |

### Participant keep file

Tab- or space-delimited text containing the identifiers of participants in the
ancestry-specific analysis. The file has an `FID IID` header:

```text
FID IID
sample_001 sample_001
sample_002 sample_002
```

### Variant extract files

Plain text with one array variant identifier per line and no header. Variant IDs
in the original analysis followed the `AX-####` convention and must match the
IDs in the corresponding genotype data exactly.

```text
AX-000000001
AX-000000002
```

## Step 1: fit the null model

Additional inputs are a PLINK 1 binary genotype dataset (`.bed`, `.bim`, and
`.fam`) sharing one prefix, and an extract file containing the common variants
used to fit the null model.

Create the log directory, export the paths, and submit:

```bash
mkdir -p logs results/step1

export REGENIE_STEP1_BIN=/path/to/regenie-3.1.3/regenie
export BED_PREFIX=/path/to/genotypes/common_variants
export EXTRACT_FILE=/path/to/step1_variants.txt
export KEEP_FILE=/path/to/amr_keep.txt
export PHENO_FILE=/path/to/phenotype.txt
export COVAR_FILE=/path/to/covariates.txt
export RESULTS_DIR="$PWD/results/step1"
export ANCESTRY=AMR

bsub -env all < 02_mvp_gwas/regenie_step1.sh
```

`BED_PREFIX` is the shared path before `.bed`, `.bim`, and `.fam`, not the path
to any one of those files. The primary downstream output is
`RA.AMR.step1_pred.list`, which points to the LOCO prediction file(s) needed by
Step 2. REGENIE also writes a log file with the same output prefix.

## Step 2: association testing by BGEN chunk

The imputed data are arranged as:

```text
GENOTYPE_DIR/
  chr1/
    chr1.01.bgen
    chr1.01.sample
    chr1.02.bgen
    chr1.02.sample
  chr2/
    ...
```

The seed file has no header and contains exactly two whitespace-delimited fields
per row: the chromosome directory and the chunk prefix. The LSF array index
selects the corresponding row.

```text
chr10 chr10.01
chr10 chr10.02
chr10 chr10.03
```

Chromosome-specific extract files are expected by default to be named
`snplist.chr1.txt`, `snplist.chr2.txt`, and so forth. Set `EXTRACT_PREFIX` if a
different prefix is used.

```bash
mkdir -p logs results/step2

export REGENIE_STEP2_BIN=/path/to/regenie-2.2.4/regenie
export GENOTYPE_DIR=/path/to/imputed_bgen
export SEED_FILE=/path/to/chunk_seed.txt
export EXTRACT_DIR=/path/to/extract_files
export EXTRACT_PREFIX=snplist
export KEEP_FILE=/path/to/amr_keep.txt
export PHENO_FILE=/path/to/phenotype.txt
export COVAR_FILE=/path/to/covariates.txt
export PRED_FILE="$PWD/results/step1/RA.AMR.step1_pred.list"
export RESULTS_DIR="$PWD/results/step2"
export ANCESTRY=AMR

N_CHUNKS=$(wc -l < "${SEED_FILE}")
bsub -env all -J "RA_step2.amr[1-${N_CHUNKS}]" \
  < 02_mvp_gwas/regenie_step2.sh
```

For seed row `chr10 chr10.01`, the script reads
`GENOTYPE_DIR/chr10/chr10.01.bgen` and its matching `.sample` file, then writes
gzipped REGENIE results with prefix `AMR.chr10.01.step2`. Each array task checks
that all required inputs exist before starting REGENIE.

