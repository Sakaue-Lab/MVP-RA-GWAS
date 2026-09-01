# MVP-RA-GWAS

This is a repository to host scripts used for analyzing the VA Million Veteran Program and perform genome-wide association study of rheumatoid arthritis (RA) followed by trans-ancestry meta-analyses and functional characterization of significant loci. Reference: Sakaue, Yang, Zhang et al. "Multi-ancestral GWAS with the VA Million Veteran Program enables functional interpretation of rheumatoid arthritis alleles" (https://www.medrxiv.org/content/10.64898/2026.04.22.26351423v1)

Each directory contains their README as well to guide the reproduction of our results.

### Table of Contents

- `01_phenotype/` : KOMAP Phenotyping of RA in the MVP
  - komap_calib.R: defining phenotype
- `02_mvp_gwas/` : GWAS script for the MVP using REGENIE
- `03_meta_analyses/` : Meta-analayses among MVP and between MVP and RACI
- `04_locus_define_finemap/` : Identify significant loci, loci annotation, fine mapping
- `05_ldsc_h2/`: Heritability and bias estimation of GWAS
- `06_gchromvar/` : Cell-type-specific enrichment of fine-mapped variants within ATAC peaks of synovial tissues at single-cell resolution
- `07_coloc/` : Colocalization of non-coding loci with TenK10K PBMC eQTL data
- `08_hla/`: HLA omnibus and conditional haplotype tests
- `09_prs/` : PRS training and testing
- `10_ldsc_gc/` : Genome-wide genetic correlation between RA and immune-mediated diseases (LDSC-rg)
- `11_buhmbox_ctPRS/` : Cell-type-specific heterogeneity of RA polygenic risk (BUHMBOX and cell-type-specific PRS vs. seropositivity)

