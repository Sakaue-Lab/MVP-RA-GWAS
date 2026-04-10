rm(list=ls())
library(dplyr)
library(stringr)
library(data.table)
library(pscl)

maf = readRDS(file="[FILE].RDS")
info = read.csv("[FILE].csv") # HLA dictionary
covar = fread("[FILE].txt")
covar$mvp001_id=as.character(covar$mvp001_id)
covar=data.frame(covar)

hla.two <- fread("[FILE].raw") # HLA AA dosage unphased
hla.two = data.frame(hla.two)
hla.two$IID=as.character(hla.two$IID)

  
#########################################################################
#####     Prepare DRB1 Data for Conditional Haplotype Testing     #######
#########################################################################
keep.col = colnames(hla.two)[grepl("HLA_DRB1",colnames(hla.two))]
hla.allele = hla.two[,c("IID",keep.col,"dosage_DRB1")]
hla.allele = hla.allele %>% 
  filter(dosage_DRB1 > 1.95,
         dosage_DRB1 < 2.05)
hla.allele = hla.allele %>%
  select_if(~!all(is.na(.)))
keep = intersect(covar$mvp001_id,hla.allele$IID)
tmp = covar %>% filter(mvp001_id %in% keep)
dat = left_join(tmp,
                hla.allele,
                by=c("mvp001_id"="IID"))
dat = as.data.frame(dat)

info.tmp = info[info$gene=="DRB1",]
info.tmp$tag=paste0(info.tmp$pos,":",info.tmp$AA)
allpos = sort( unique(info.tmp$pos) )
condpos=11
condpos2=13
condpos3=71
res=data.frame()
for(pos in allpos){
  x= info.tmp[info.tmp$pos %in% c(condpos), ]
  y= info.tmp[info.tmp$pos %in% c(condpos2), ]
  z= info.tmp[info.tmp$pos %in% c(condpos3), ]
  a= info.tmp[info.tmp$pos %in% c(pos), ]
  for(i in 1:length(unique(x$tag))){
  for(j in 1:length(unique(y$tag))){
  for(k in 1:length(unique(z$tag))){
  for(l in 1:length(unique(a$tag))){
    xtag = unique(x$tag)[i]
    ytag = unique(y$tag)[j]
    ztag = unique(z$tag)[k]
    atag = unique(a$tag)[l]
    hap = paste0(xtag,"_",ytag,"_",ztag,"_",atag)
    x_4d = subset(x,tag==xtag)$hla
    y_4d = subset(y,tag==ytag)$hla
    z_4d = subset(z,tag==ztag)$hla
    a_4d = subset(a,tag==atag)$hla
    hap_4d = intersect(x_4d,intersect(y_4d,intersect(z_4d,a_4d)))
    if(length(hap_4d)>0){
      out=data.frame(hap,hla=hap_4d,pos)
      res=rbind(res,out)
    } } }} }}

res.prev = res
add_pos="74"
thishaps = as.character(unique(subset(res.prev,
                                      pos==add_pos)$hap))
adopted.prev=NULL
for(thishap in thishaps){
  hlas=as.character(subset(res.prev,hap==thishap)$hla)
  hlas=hlas[hlas %in% colnames(dat)]
  if(length(hlas)>0){
    dat$thishap = rowSums(dat[hlas])
    colnames(dat)[ncol(dat)]=thishap
    adopted.prev=c(adopted.prev,thishap)
  }
}
adopted.prev
drb1.hap = dat[,c("mvp001_id",adopted.prev)]

#### Construct haplotypes for HLA-B position 9 
keep.col = colnames(hla.two)[grepl("HLA_B",colnames(hla.two))]
hla.allele = hla.two[,c("IID",keep.col,"dosage_B")]
hla.allele = hla.allele %>% 
  filter(dosage_B > 1.95,
         dosage_B < 2.05)
hla.allele = hla.allele %>%
  select_if(~!all(is.na(.)))
keep = intersect(covar$mvp001_id,hla.allele$IID)
tmp = covar %>% filter(mvp001_id %in% keep)
dat = left_join(tmp,
                hla.allele,
                by=c("mvp001_id"="IID"))
dat = as.data.frame(dat)

info.tmp = info[info$gene=="B",]
info.tmp$tag=paste0(info.tmp$pos,":",info.tmp$AA)
allpos = sort( unique(info.tmp$pos) )
res = data.frame()
for (pos in allpos){
  y = info.tmp[info.tmp$pos %in% c(pos), ]
  for (k in 1:length(unique(y$tag))){
    ytag = unique(y$tag)[k]
    hap = ytag
    y_4d = subset(y,tag==ytag)$hla
    hap_4d = y_4d
    if(length(hap_4d)>0){
      out=data.frame(hap,hla=hap_4d,pos)
      res=rbind(res,out)
    }
  }
}
res.prev = res
condpos="9"
thishaps = as.character(unique(subset(res.prev,
                                      pos==condpos)$hap))
adopted.prev=NULL
for(thishap in thishaps){
  hlas=as.character(subset(res.prev,hap==thishap)$hla)
  hlas=hlas[hlas %in% colnames(dat)]
  if(length(hlas)>0){
    dat$thishap = rowSums(dat[hlas])
    colnames(dat)[ncol(dat)]=thishap
    adopted.prev=c(adopted.prev,thishap)
  }
}
hlab.hap = dat[,c("mvp001_id",adopted.prev)]
cond.hap = left_join(hlab.hap,
                     drb1.hap)

#### Construct haplotypes for HLA-DQB1 position 37 
keep.col = colnames(hla.two)[grepl("HLA_DPB1",colnames(hla.two))]
hla.allele = hla.two[,c("IID",keep.col,"dosage_DPB1")]
hla.allele = hla.allele %>% 
  filter(dosage_DPB1 > 1.95,
         dosage_DPB1 < 2.05)
hla.allele = hla.allele %>%
  select_if(~!all(is.na(.)))
keep = intersect(covar$mvp001_id,hla.allele$IID)
tmp = covar %>% filter(mvp001_id %in% keep)
dat = left_join(tmp,
                hla.allele,
                by=c("mvp001_id"="IID"))
dat = as.data.frame(dat)

info.tmp = info[info$gene=="DPB1",]
info.tmp$tag=paste0(info.tmp$pos,":",info.tmp$AA)
allpos = sort( unique(info.tmp$pos) )
res = data.frame()
for (pos in allpos){
  y = info.tmp[info.tmp$pos %in% c(pos), ]
  for (k in 1:length(unique(y$tag))){
    ytag = unique(y$tag)[k]
    hap = ytag
    y_4d = subset(y,tag==ytag)$hla
    hap_4d = y_4d
    if(length(hap_4d)>0){
      out=data.frame(hap,hla=hap_4d,pos)
      res=rbind(res,out)
    }
  }
}
res.prev = res
condpos="9"
thishaps = as.character(unique(subset(res.prev,
                                      pos==condpos)$hap))
adopted.prev=NULL
for(thishap in thishaps){
  hlas=as.character(subset(res.prev,hap==thishap)$hla)
  hlas=hlas[hlas %in% colnames(dat)]
  if(length(hlas)>0){
    dat$thishap = rowSums(dat[hlas])
    colnames(dat)[ncol(dat)]=thishap
    adopted.prev=c(adopted.prev,thishap)
  }
}
hladqb1.hap = dat[,c("mvp001_id",adopted.prev)]
colnames(hladqb1.hap)[-1]=paste0("dpb1_",colnames(hladqb1.hap)[-1])
cond.hap = left_join(cond.hap,
                     hladqb1.hap)

keep = intersect(covar$mvp001_id,cond.hap$mvp001_id)
tmp = covar %>% filter(mvp001_id %in% keep)
dat = left_join(tmp,
                cond.hap)
dat = as.data.frame(dat)
write.csv(dat,
          file="/group/research/mvp001/hzhang/RA_GWAS/data/haplotype_eur_aa_seq_20240916.csv",
          row.names = F)

adopted.prev = colnames(cond.hap)[-1]

hap.effect=NULL
for(anc in c("AFR","EUR")){
  for(hap in adopted.prev){
    print(paste0(anc," ",hap))
    tryCatch({
  tmp = dat %>% filter(ancestry==anc)
  colnames(tmp)[which(colnames(tmp)==hap)]="haplotype"
  tmp$cond1 = tmp[,"9:D"]
  tmp$cond2 = tmp[,"dpb1_9:F"]
  obj2=glm(binary_95~haplotype+cond1+cond2+sex+age+pc1+pc2+pc3+pc4+pc5,
           data=tmp,family=binomial(link="logit"))
  psudor = pR2(obj2)["McFadden"]
  obj2=summary(obj2)
  hap.effect = rbind.data.frame(hap.effect,
                                data.frame("ancestry"=anc,
                                           "cond_gene"="out_drb1",
                                           "haplotype"=hap,
                                           "estimate"=obj2$coefficients["haplotype",1],
                                           "OR"=NA,
                                           "se"=obj2$coefficients["haplotype",2],
                                           "pval"=obj2$coefficients["haplotype",4],
                                           "maf"=sum(tmp[,"haplotype"],na.rm = T)/(2*(nrow(tmp)-sum(is.na(tmp[,"haplotype"])))),
                                           "psuedor"=psudor))
    },error=function(e){NA})
  }
}
hap.effect$OR = round(exp(hap.effect$estimate),digits = 2)
hap.effect = hap.effect %>% filter(maf >= 0.05)

write.csv(hap.effect,
          row.names = F,
          file="[FILE].csv")


tmp = dat %>% filter(ancestry=="AFR") %>%
  dplyr::select(-c(mvp001_id,ancestry,pc6,pc7,pc8,pc9,pc10,
                binary_99,binary_90,komap_score))
obj2=glm(binary_95 ~.,
         data=tmp,family=binomial(link="logit"))
psudor = pR2(obj2)["McFadden"]

tmp = dat %>% filter(ancestry=="EUR") %>%
  dplyr::select(-c(mvp001_id,ancestry,pc6,pc7,pc8,pc9,pc10,
                   binary_99,binary_90,komap_score))
obj2=glm(binary_95 ~.,
         data=tmp,family=binomial(link="logit"))
pR2(obj2)["McFadden"]

tmp = dat %>% filter(ancestry=="AFR") %>%
  dplyr::select(-c(mvp001_id,ancestry,pc6,pc7,pc8,pc9,pc10,
                   binary_99,binary_90,binary_95))
tmp$komap_score[is.na(tmp$komap_score)]=0
obj2=lm(komap_score ~.,
         data=tmp)
psudor = pR2(obj2)["McFadden"]

tmp = dat %>% filter(ancestry=="EUR") %>%
  dplyr::select(-c(mvp001_id,ancestry,pc6,pc7,pc8,pc9,pc10,
                   binary_99,binary_90,binary_95))
tmp$komap_score[is.na(tmp$komap_score)]=0
obj2=lm(komap_score ~.,
        data=tmp)
pR2(obj2)["McFadden"]

