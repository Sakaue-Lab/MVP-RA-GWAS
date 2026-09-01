[← Repository overview](../README.md)

# 02. MVP genome-wide association analysis

The scripts in `02_mvp_gwas/` run binary-trait analyses with REGENIE. Step 1
fits the null model using directly genotyped variants. The standard Step 2
script tests imputed variants stored in chromosome- and chunk-specific BGEN
files. A separate Step 2 script tests SNP-by-sex interactions for a targeted
set of variants. The scripts retain the historical software versions used for
the analysis: REGENIE 3.1.3 for Step 1, REGENIE 2.2.4 for the primary Step 2
GWAS, and REGENIE 3.4 for interaction testing.

The shell scripts are configured through exported environment variables so that
restricted filesystem paths do not need to be hard-coded in the repository.
They contain LSF `#BSUB` directives and were designed for submission with
`bsub`. Queue, project, memory, and log settings may require adaptation outside
the original MVP computing environment.

## Scripts and execution order

| Order | Script | Analysis | REGENIE version |
|---|---|---|---|
| 1 | `regenie_step1.sh` | Fit the ancestry-specific null logistic mixed model | 3.1.3 |
| 2a | `regenie_step2.sh` | Primary binary-trait single-variant GWAS | 2.2.4 |
| 2b | `regenie_step2_interaction.sh` | Targeted SNP-by-sex interaction testing | 3.4 |

### There are two separate Step 2 scripts

Run `regenie_step1.sh` first. Its prediction-list output is then used by **both**
of these independent Step 2 analyses:

1. `regenie_step2.sh` runs the standard genome-wide association analysis.
2. `regenie_step2_interaction.sh` runs the targeted SNP-by-sex interaction
   analysis.

The interaction script does not replace the standard Step 2 script, and it does
not use the standard Step 2 results as input. Both Step 2 scripts start from the
Step 1 prediction list and chunked BGEN files:

```text
regenie_step1.sh
  ├── regenie_step2.sh                 standard GWAS
  └── regenie_step2_interaction.sh     SNP-by-sex interaction tests
```

## General execution setup

Run the commands below from the repository root so the relative script and log
paths resolve correctly:

```bash
cd /path/to/MVP-RA-GWAS
chmod +x 02_mvp_gwas/regenie_step1.sh
chmod +x 02_mvp_gwas/regenie_step2.sh
chmod +x 02_mvp_gwas/regenie_step2_interaction.sh
mkdir -p logs
```

Each script stops if a required environment variable or input file is missing.
Export the variables in the same terminal used to submit the job, and retain
`-env all` in the `bsub` command so LSF passes those variables to the job.

## Input files shared by all three scripts

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

## 1. Step 1: fit the null model

Additional inputs are a PLINK 1 binary genotype dataset (`.bed`, `.bim`, and
`.fam`) sharing one prefix, and an extract file containing the common variants
used to fit the null model.

### Execute Step 1

Create the results directory, export the required settings, and submit one
Step 1 job. Step 1 is not an array job.

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
to any one of those files.

After submission, monitor the job and inspect its logs:

```bash
bjobs -J regenie_step1
less logs/regenie_step1.<JOB_ID>.stdout
less logs/regenie_step1.<JOB_ID>.stderr
```

For `ANCESTRY=AMR`, the output prefix is `RA.AMR.step1`. Confirm that
`RA.AMR.step1_pred.list` exists and points to the generated LOCO prediction
file(s). Supply this file to either Step 2 script through `PRED_FILE`.

## 2. Standard Step 2: association testing by BGEN chunk

Script: `regenie_step2.sh`

This is the first of the two separate Step 2 scripts. It performs the primary
GWAS and does not include an interaction term.

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

### Execute the standard Step 2 GWAS

1. Confirm that Step 1 completed and identify its `.pred.list` file.
2. Confirm that the seed file has one row for every BGEN chunk to test.
3. Export the genotype, extract, participant, phenotype, covariate,
   prediction-list, and output paths.
4. Count the seed-file rows and submit one LSF array task per row.

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

Monitor the standard Step 2 array and inspect task-specific logs:

```bash
bjobs -J 'RA_step2.amr[*]'
less logs/regenie_step2.<JOB_ID>.<ARRAY_INDEX>.stdout
less logs/regenie_step2.<JOB_ID>.<ARRAY_INDEX>.stderr
```

The results directory should contain one gzipped REGENIE result for each
successfully tested BGEN chunk. Combine only completed chunk outputs.

## 3. Interaction Step 2: SNP-by-sex interaction testing

`regenie_step2_interaction.sh` tests whether the effect of each targeted variant
differs by sex. It uses the same chunked BGEN data, seed file, phenotype,
covariates, ancestry-specific keep file, and Step 1 prediction list described
above. The published analysis used REGENIE 3.4 because interaction testing was
performed after the primary GWAS.

This is the second, independent Step 2 script. It uses the Step 1 `.pred.list`
file directly. The standard `regenie_step2.sh` analysis does not need to finish
first, and its association results are not inputs to the interaction analysis.

The targeted extract file is a single plain-text file with one variant ID per
line and no header. Unlike the standard Step 2 script, its path is supplied
directly through `EXTRACT_FILE`; it is not assumed to have separate files for
each chromosome. Variant IDs must match the BGEN data exactly.

Sex must be present in the covariate file as a categorical variable. By default,
the script passes `--catCovarList sex` and `--interaction 'sex[0]'` to REGENIE.
The exact meaning of level `0` depends on the coding in the covariate file, so
users must confirm which sex category is encoded as `0` before interpreting the
interaction coefficient. Age and `pc1`–`pc5` are included as adjustment
covariates.

### Execute Step 2 interaction testing

1. Confirm that Step 1 completed and identify its `.pred.list` file.
2. Prepare one targeted extract file containing all SNPs to test.
3. Confirm that `sex` is present in the covariate file and verify which category
   is encoded as `0`.
4. Export the required file paths and interaction settings.
5. Count the seed-file rows and submit one LSF array task per BGEN chunk.

```bash
mkdir -p logs results/interaction

export REGENIE_INTERACTION_BIN=/path/to/regenie-3.4/regenie
export GENOTYPE_DIR=/path/to/imputed_bgen
export SEED_FILE=/path/to/chunk_seed.txt
export EXTRACT_FILE=/path/to/target_variants.txt
export KEEP_FILE=/path/to/eur_keep.txt
export PHENO_FILE=/path/to/phenotype.txt
export COVAR_FILE=/path/to/covariates.txt
export PRED_FILE=/path/to/RA.EUR.step1_pred.list
export RESULTS_DIR="$PWD/results/interaction"
export ANCESTRY=EUR
export PHENO_COLUMN=binary_ppv_90_rm
export INTERACTION_COVARIATE=sex
export INTERACTION_LEVEL=0

N_CHUNKS=$(wc -l < "${SEED_FILE}")
bsub -env all -J "RA_step2_interaction.eur[1-${N_CHUNKS}]" \
  < 02_mvp_gwas/regenie_step2_interaction.sh
```

The optional variables `INTERACTION_COVARIATE` and `INTERACTION_LEVEL` default
to `sex` and `0`, respectively. For example, setting
`INTERACTION_LEVEL=1` changes the requested interaction term to `sex[1]`.

For seed row `chr10 chr10.01`, the default output prefix is
`EUR.chr10.01.sex_interaction.step2`. REGENIE writes gzipped association results
and a log using this prefix. Because the extract file contains only targeted
variants, array tasks for chunks containing no listed variants may finish
without association results; review the REGENIE logs before combining outputs.

Monitor the interaction array separately from the standard Step 2 array:

```bash
bjobs -J 'RA_step2_interaction.eur[*]'
less logs/regenie_step2_interaction.<JOB_ID>.<ARRAY_INDEX>.stdout
less logs/regenie_step2_interaction.<JOB_ID>.<ARRAY_INDEX>.stderr
```

Retain outputs containing the requested target variants. When combining the
results, select the interaction-test rows reported by REGENIE rather than the
variant main-effect rows.
