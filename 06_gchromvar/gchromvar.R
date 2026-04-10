options(stringsAsFactors=F)
#library(Signac)
library(Seurat)
library(ggplot2)
library(Matrix)
library(stringr)
library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
library(ggplot2)
library(dplyr)

meta<-readRDS("/path/metadata.rds")

atac<-readRDS("/path/atac.rds")

celltypes <- names(table(meta$mrna_fine_midRes_CT_symp)[table(meta$mrna_fine_midRes_CT_symp) > 500])

# pseudobulk atac matrix
allcells <- colnames(atac)
counts <- matrix(0, nrow=nrow(atac), ncol=length(celltypes) )
rownames(counts) <- rownames(atac)
colnames(counts) <- celltypes

for(cell in celltypes){
  cells <- subset(meta, mrna_fine_midRes_CT_symp==cell )$atacCN
  cells <- intersect(cells, allcells)
  d2 <- as.matrix(atac[, cells])
  # rowaverage <- rowSums(d2) / ncol(d2)
  # counts[, cell] <- rowaverage
  rowsum <- rowSums(d2)
  counts[, cell] <- rowsum
}

saveRDS(counts, "pseudobulk_atac.rds")

peaks<- GRanges(
    seqnames = str_split_fixed(rownames(atac),":",2)[,1],
    ranges = IRanges(start = as.numeric(str_split_fixed(str_split_fixed(rownames(atac),":",2)[,2],"-",2)[,1]), end = as.numeric(str_split_fixed(str_split_fixed(rownames(atac),":",2)[,2],"-",2)[,2]))
)

se <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = peaks,
  colData = DataFrame(cell = colnames(counts))
)

se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)

# load RA GWAS fine-mapping results at hg38
files<-c("/path/finemapped_our_gwas.hg38.bed")

# Computing weighted deviations
# many trials to derive SE

stats<-data.frame()

for(i in c(1:1000)){
  ra<-gchromVAR::importBedScore(rowRanges(se), files, colidx = 5)
  ra_wDEV <- computeWeightedDeviations(se, ra)
  zdf <- reshape2::melt(t(assays(ra_wDEV)[["z"]]))
  colnames(zdf) <- c("ct", "tr", "Zscore")
  zdf$gchromVAR_pvalue <- pnorm(zdf$Zscore, lower.tail = FALSE)
  zdf$iter<-i
  stats<-rbind(stats,zdf)
}

saveRDS(stats,"stats/gchromVAR.RA_summary.iter1000.rds")