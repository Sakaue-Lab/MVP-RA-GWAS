# 11_buhmbox_ctPRS

Cell-type-specific heterogeneity of RA polygenic risk. We asked whether PRS
variants located in cell-type-specific active regulatory regions characterize
latent heterogeneity among RA cases.

## Workflow
1. `define_celltype_prs_variants.sh` : assign best-performing PRS variants to
   immune cell types by overlap with cell-type-specific active regulatory
   elements (single-cell ATAC-seq open-chromatin peaks from RA synovial tissue),
   extended by 250 bp on either end. Assignment is mutually non-exclusive.
2. `run_buhmbox.R` : run BUHMBOX (Han et al., Nat Genet 2016) on each
   cell-type-specific variant set, the union of cell types, and the full
   genome-wide PRS variant set, to detect a latent case subgroup from excess
   positive correlation among risk-increasing alleles. Writes per-set
   `*_buhmbox.weights` and `*_buhmbox.BBrst`, and `stats/buhmbox_summary.tsv`.
3. `ctPRS_seropos.R` : in held-out MVP RA cases (30% not used to train the
   overall PRS), fit `seropositive ~ ctPRS + overall_PRS + covariates` to estimate
   the additional contribution of each cell-type-specific PRS (ctPRS) on
   seropositivity, conditional on the overall PRS. Writes `prs_seropos_res.csv`.

## Notes
- BUHMBOX is an external method; set the `source()` path in `run_buhmbox.R` to your
  BUHMBOX installation (https://github.com/immunogenomics/BUHMBOX).
- Replace `/path_to/...` placeholders with your local paths.
- Cell-type ATAC peaks are the same synovial single-cell peaks used in
  `06_gchromvar/`.
