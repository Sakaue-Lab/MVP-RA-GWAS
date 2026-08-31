# Define cell-type-specific subsets of the best-performing PRS variants by overlap
# with cell-type-specific active regulatory elements (single-cell ATAC-seq
# open-chromatin peaks from RA synovial tissue), extended by 250 bp on either end.
# Assignment is mutually non-exclusive (a variant may be assigned to >1 cell type).

# inputs:
#   stats/prs_best_variants.hg38.bed          : clumped PRS variants (hg38), col4 = varID
#   ../06_gchromvar/peaks/<celltype>.peaks.bed : cell-type ATAC peaks (hg38)

mkdir -p stats/celltype_variants

for celltype in Bcell DC endothelial fibroblast monocyte Tcell
do
  # extend peaks by 250 bp each side, then intersect with PRS variants
  bedtools slop -i ../06_gchromvar/peaks/${celltype}.peaks.bed \
                -g /path_to/hg38.chrom.sizes -b 250 \
  | bedtools intersect -a stats/prs_best_variants.hg38.bed -b - -wa \
  | cut -f4 | sort -u > stats/celltype_variants/${celltype}.variants.txt
  echo "${celltype}: $(wc -l < stats/celltype_variants/${celltype}.variants.txt) variants"
done

# union across all cell types (any regulatory element)
cat stats/celltype_variants/*.variants.txt | sort -u \
  > stats/celltype_variants/unioncells.variants.txt
