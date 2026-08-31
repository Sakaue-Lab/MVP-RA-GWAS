library(data.table)
library(stringr)
library(dplyr)
library(coloc)

celltypes<-c("ASDC","B_intermediate","B_memory","B_naive","CD14_Mono","CD16_Mono","CD4_CTL","CD4_Naive","CD4_Proliferating","CD4_TCM","CD4_TEM","CD8_Naive","CD8_Proliferating","CD8_TCM","CD8_TEM","cDC1","cDC2","dnT","gdT","HSPC","ILC","MAIT","NK_CD56bright","NK","NK_Proliferating","pDC","Plasmablast","Treg")

egenes<-fread("/path_to/ST4.txt.gz")

res<-data.frame()

for(cell in celltypes){
  d<-fread(paste0("/path_to/TenK10K/",cell,"_common_all_cis_raw_pvalues_1000000bp.tsv.gz"))
  egene_list<-egenes[egenes$`Cell Type` == cell,]$`Gene (Ensembl gene ID)`
  d<-d[d$gene %in% egene_list,]
  for(categ in c("trans","eur","trans_seropos","eur_seropos")){
    markers<-read.table(paste0("stats/coloc/all_leadmarkers_",categ,".txt"))$V1
    for(marker in markers){
      gwas<-fread(paste0("stats/coloc/",marker,"_locus.sumstats_hg38.bed"))
      gwas$Marker<-paste(gsub("chr","",gwas$V1),gwas$V3,sep=":")
      gwas$alleles<-str_split_fixed(gwas$V4,":",3)[,3]
      gwas$Marker<-paste(gwas$Marker,gwas$alleles,sep=":")
      d2<-data.frame(d[d$MarkerID %in% gwas$Marker,])
      if(nrow(d2) > 200){
        n_snps<-d2 %>% group_by(gene) %>% summarise(n=n(), minP=min(p.value))
        n_snps<-n_snps[n_snps$n > 200,]
        if(nrow(n_snps) > 0){
          write(n_snps$gene,paste0("stats/coloc/",marker,"_locus.",cell,"_eQTL.genes.txt"))
          for(gene in n_snps$gene){ 
            eqtl<-d2[d2$gene == gene,]
            eqtl$MAF<-eqtl$AF_Allele2
            if(nrow(eqtl[eqtl$MAF>0.5,])>0){eqtl[eqtl$MAF>0.5,]$MAF<-1-eqtl[eqtl$MAF>0.5,]$MAF}
            eqtl_coloc <- list(varbeta = (eqtl$SE)^2, N=eqtl$N, MAF= eqtl$MAF, type = "quant", beta=eqtl$BETA, snp=as.character(eqtl$MarkerID),sdY=1)
            this_gwas<-gwas[gwas$Marker %in% eqtl$MarkerID ]
            gwas_coloc<-list(N=this_gwas$V9, type = "cc", beta=this_gwas$V7, varbeta=(this_gwas$V8)^2, snp=as.character(this_gwas$Marker))
            tryCatch({
coloc_res = coloc::coloc.abf(eqtl_coloc, gwas_coloc)
                                PP.H4.abf<-coloc_res$summary["PP.H4.abf"]
                                PP.H3.abf<-coloc_res$summary["PP.H3.abf"]
                                nsnps<-nrow(this_gwas)
                                tmp<-data.frame(cell=cell,marker=marker,gene=gene,PP.H3.abf=PP.H3.abf,PP.H4.abf=PP.H4.abf,nsnps=nsnps)
								res<-rbind(res,tmp)}, error = function(cond){
                        nsnps<-nrow(this_gwas)
                        tmp<-data.frame(cell=cell,marker=marker,gene=gene,PP.H3.abf=NA,PP.H4.abf=NA,nsnps=nsnps)
                        res<-rbind(res,tmp)
                        })}
                     }}}}
if(nrow(res)>0){saveRDS(res,paste0("stats/coloc/gwas_TenK10K_",cell,"_coloc_v2.rds"))}
}
