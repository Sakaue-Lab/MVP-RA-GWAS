[← Repository overview](../README.md)

# 07. Colocalization with single-cell eQTL

Test whether RA GWAS signals share a causal variant with cell-type-specific
eQTL from peripheral blood immune cells (TenK10K Phase 1).

- `coloc.R` : for each of the 28 PBMC cell types, restrict to significant
  eGenes, and for each RA lead locus (across the four GWAS settings: trans /
  eur / trans-seropositive / eur-seropositive) extract the overlapping eQTL and
  GWAS summary statistics. Genes with more than 200 shared variants are tested
  with `coloc.abf` (approximate Bayes factor), modeling the GWAS as a
  case-control trait and the eQTL as a quantitative trait. Posterior
  probabilities of a shared causal variant (PP.H4) and of distinct causal
  variants (PP.H3) are recorded per locus–gene–cell-type; significant
  colocalization is defined as PP.H4 > 0.7.

GWAS locus summary statistics are lifted to hg38 before matching to the eQTL
(shared marker IDs). Replace the `/path_to/...` eQTL/eGene paths with your local
copies.
