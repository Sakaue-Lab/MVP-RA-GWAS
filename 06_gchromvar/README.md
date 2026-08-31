[← Repository overview](../README.md)

# 06. Cell-type enrichment of fine-mapped variants (g-chromVAR)

Test whether fine-mapped RA variants are enriched within cell-type-specific
accessible chromatin of synovial tissue at single-cell resolution.

- `gchromvar.R` : pseudobulk single-cell ATAC-seq peak counts across cell
  states (states with >500 cells), add GC bias, import the fine-mapped RA
  variants (95% credible set, hg38) weighted by posterior inclusion probability,
  and run `gchromVAR` to compute a per-cell-state enrichment Z score. The
  weighted-deviation step is repeated (1,000 iterations) to derive the mean and
  standard error of the Z score per cell state.

Inputs are the single-cell ATAC peak-by-cell matrix and cell metadata (paths
marked `/path/...`) and the fine-mapped variant BED from `04_locus_define_finemap/`.
