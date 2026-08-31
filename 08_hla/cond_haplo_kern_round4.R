rm(list=ls())
library(dplyr)
library(stringr)
library(data.table)
args = commandArgs(trailingOnly = TRUE)
if(length(args) != 6) stop("Usage: Rscript cond_haplo_kern_round4.R <dictionary.csv> <covariates.txt> <dosages.raw> <exclude.csv> <variants.pvar> <output.csv>")
info = read.csv(args[1]) # HLA dictionary
covar = fread(args[2])
covar$id=as.character(covar$id)
covar=data.frame(covar)

hla.two <- fread(args[3]) # HLA AA dosage unphased
hla.two = data.frame(hla.two)
hla.two$IID=as.character(hla.two$IID)

exclude = fread(args[4]) # exclude related participants
covar = covar %>% filter(!(id %in% exclude$IID))

### Post-imputation QC, consistent with rounds 1-3
r2 = fread(args[5])
r2$r2 = unlist(lapply(r2$INFO,FUN=function(ll){
  junk = (strsplit(ll,";")[[1]][3])
  return(as.numeric(strsplit(junk,"=")[[1]][2]))
}))
r2$maf = unlist(lapply(r2$INFO,FUN=function(ll){
  junk = (strsplit(ll,";")[[1]][2])
  return(as.numeric(strsplit(junk,"=")[[1]][2]))
}))
r2 = r2 %>% filter(grepl("[*]",ID), !grepl("exon",ID))
varkeep = r2 %>% filter(maf > 0.01, r2 > 0.7, grepl(":",ID))
varkeep$format = paste0(gsub("[*]|:",".",varkeep$ID),"_T..A.")

######################################################
######   Implement conditional haplotype test ########
######################################################
rd4.export = NULL
for(HLA in unique(info$gene)){
  
  
#########################################################################
#####     Prepare DRB1 Data for Conditional Haplotype Testing     #######
#########################################################################
keep.col = colnames(hla.two)[grepl("HLA_DRB1",colnames(hla.two))]
keep.col = intersect(keep.col,varkeep$format)
hla.allele = hla.two[,c("IID",keep.col,"dosage_DRB1")]
hla.allele = hla.allele %>% 
  filter(dosage_DRB1 > 1.95,
         dosage_DRB1 < 2.05)
hla.allele = hla.allele %>%
  select_if(~!all(is.na(.)))
keep = intersect(covar$id,hla.allele$IID)
tmp = covar %>% filter(id %in% keep)
dat = left_join(tmp,
                hla.allele,
                by=c("id"="IID"))
dat = as.data.frame(dat)

info.tmp = info[info$gene=="DRB1",]
info.tmp$tag=paste0(info.tmp$pos,":",info.tmp$AA)
allpos = sort( unique(info.tmp$pos) )
condpos=11
res=data.frame()
for(pos in allpos){
  x= info.tmp[info.tmp$pos %in% c(condpos), ]
  y= info.tmp[info.tmp$pos %in% c(pos), ]
  for(i in 1:length(unique(x$tag))){
  for(k in 1:length(unique(y$tag))){
    xtag = unique(x$tag)[i]
    ytag = unique(y$tag)[k]
    hap = paste0(xtag,"_",ytag)
    x_4d = subset(x,tag==xtag)$hla
    y_4d = subset(y,tag==ytag)$hla
    hap_4d = intersect(x_4d,y_4d)
    if(length(hap_4d)>0){
      out=data.frame(hap,hla=hap_4d,pos)
      res=rbind(res,out)
    } } }}

res.prev = res
anc="AFR"
#anc="EUR"
add_pos="85"
condpos1 = "11"
condpos2 = "85"
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
drb1.hap = dat[,c("id",adopted.prev)]

#### Construct haplotypes for HLA-B position 9 
keep.col = colnames(hla.two)[grepl("HLA_B",colnames(hla.two))]
keep.col = intersect(keep.col,varkeep$format)
hla.allele = hla.two[,c("IID",keep.col,"dosage_B")]
hla.allele = hla.allele %>% 
  filter(dosage_B > 1.95,
         dosage_B < 2.05)
hla.allele = hla.allele %>%
  select_if(~!all(is.na(.)))
keep = intersect(covar$id,hla.allele$IID)
tmp = covar %>% filter(id %in% keep)
dat = left_join(tmp,
                hla.allele,
                by=c("id"="IID"))
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
hlab.hap = dat[,c("id",adopted.prev)]
cond.hap = left_join(hlab.hap,
                     drb1.hap)

#### Start haplotype testing 
info.tmp = info[info$gene==HLA,]
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

keep.col = colnames(hla.two)[grepl(paste0("HLA_",HLA),colnames(hla.two))]
keep.col = intersect(keep.col,varkeep$format)
hla.allele = hla.two[,c("IID",keep.col,colnames(hla.two)[grepl(paste0("dosage_",HLA),colnames(hla.two))])]
hla.allele = hla.allele[hla.allele[,colnames(hla.two)[grepl(paste0("dosage_",HLA),colnames(hla.two))]]>1.95,]
hla.allele = hla.allele[hla.allele[,colnames(hla.two)[grepl(paste0("dosage_",HLA),colnames(hla.two))]]<2.05,]
hla.allele = hla.allele %>%
  select_if(~!all(is.na(.)))
keep = intersect(covar$id,hla.allele$IID)
tmp = covar %>% filter(id %in% keep)
dat = left_join(tmp,
                hla.allele,
                by=c("id"="IID"))
dat = as.data.frame(dat)

keep = intersect(dat$id,cond.hap$id)
dat = dat %>% filter(id %in% keep)
dat = left_join(dat,
                cond.hap)
dat = as.data.frame(dat)

adopted.prev = colnames(cond.hap)[-1]

rd.cond.drb1=NULL
for(anc in c("AFR")){
pval_list = NULL
deviance_list = NULL
obj1=glm(binary_ppv_90_rm_trans~as.matrix(dat[dat$ancestry==anc,adopted.prev])+sex+age+pc1+pc2+pc3+pc4+pc5,data=subset(dat,ancestry==anc),family=binomial(link="logit"))
for(thispos in allpos){
  print(paste0(anc," ",HLA," ",which(allpos==thispos),"/",length(allpos)))

  thishaps = as.character(unique(subset(res,pos==thispos)$hap))
  adopted=NULL
  for(thishap in thishaps){
    hlas=as.character(subset(res,hap==thishap)$hla)
    hlas=hlas[hlas %in% colnames(dat)]
    if(length(hlas)>0){
      dat$thishap = rowSums(dat[hlas])
      colnames(dat)[ncol(dat)]=thishap
      adopted=c(adopted,thishap)
    }
  }
  obj2=glm(binary_ppv_90_rm_trans~as.matrix(dat[dat$ancestry==anc,c(adopted.prev,adopted)])+sex+age+pc1+pc2+pc3+pc4+pc5,data=subset(dat,ancestry==anc),family=binomial(link="logit"))
  chisqtest = anova(obj1,obj2,test="Chisq")
  pval = chisqtest$`Pr(>Chi)`[2]
  deviance = chisqtest$Deviance[2]
  pval_list=c(pval_list,pval)
  deviance_list=c(deviance_list,deviance)
  }
summary = data.frame(ancestry =anc,
                     position = allpos,
                     omnibus_deviance = deviance_list,
                     omnibus_pval = pval_list)
rd.cond.drb1=rbind.data.frame(rd.cond.drb1,
                     summary)
}

rd.cond.drb1$HLA=HLA
rd4.export = rbind.data.frame(rd4.export,
                              rd.cond.drb1)

}

dir.create(dirname(args[6]), recursive=TRUE, showWarnings=FALSE)
write.csv(rd4.export,
          row.names = F,
          file=args[6])
