rm(list=ls())
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

maf = readRDS(file="[FILE].RDS")
info = read.csv("[FILE].csv") # HLA dictionary
covar = fread("[FILE].txt")
covar$mvp001_id=as.character(covar$mvp001_id)
covar=data.frame(covar)

hla.two <- fread("[FILE].raw") # HLA AA dosage unphased
hla.two = data.frame(hla.two)
hla.two$IID=as.character(hla.two$IID)

exclude = fread("[FILE].csv") # exclude related
covar = covar %>% filter(!(mvp001_id %in% exclude$IID))

### Post-Impute QC
r2 = fread("[FILE].pvar")
r2$r2 = unlist(lapply(r2$INFO,FUN=function(ll){
  junk = (strsplit(ll,";")[[1]][3])
  return(as.numeric(strsplit(junk,"=")[[1]][2]))
}))
r2$maf = unlist(lapply(r2$INFO,FUN=function(ll){
  junk = (strsplit(ll,";")[[1]][2])
  return(as.numeric(strsplit(junk,"=")[[1]][2]))
}))
r2 = r2 %>% filter(grepl("[*]",ID),
                   !grepl("exon",ID))
r2$MAF_bin = unlist(lapply(r2$maf,FUN=function(ll){
  if(ll<.01){
    return(1)
  }else if(ll>=0.01 & ll<0.05){
    return(2)
  }else if (ll>=0.05 & ll<0.1){
    return(3)
  }else if (ll>=0.1){
    return(4)
  }
}))


varkeep = r2 %>% 
  filter(maf > 0.01,
         r2 > 0.7,
         grepl(":",ID))
varkeep$format = paste0(gsub("[*]|:",".",varkeep$ID),"_T..A.")

######################################################
######   Implement conditional haplotype test ########
######################################################
rd1.export = NULL
for(HLA in unique(info$gene)){

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


#########################################################################
#####    Prepare Data for Forward Conditional Haplotype Testing   #######
keep.col = colnames(hla.two)[grepl(paste0("HLA_",HLA),colnames(hla.two))]
keep.col = intersect(keep.col,varkeep$format)
hla.allele = hla.two[,c("IID",keep.col)]
junk = rowSums(hla.allele[,-1])
junk = junk < 2.10 & junk > 1.90
hla.allele = hla.allele[junk,]
hla.allele = hla.allele %>%
  select_if(~!all(is.na(.)))

revise = colnames(hla.allele)[-1]
revise = gsub("_T..A.","",revise)
revise = gsub("[.]","_",revise)
colnames(hla.allele)[-1]=revise

keep = intersect(covar$mvp001_id,hla.allele$IID)
tmp = covar %>% filter(mvp001_id %in% keep)
dat = left_join(tmp,
                hla.allele,
                  by=c("mvp001_id"="IID"))
dat = as.data.frame(dat)

rd1 = NULL
outcome="binary_ppv_90_rm_trans"
for(anc in c("AFR","EUR")){
  #for(anc in c("AFR")){
  pval_list=NULL
  deviance_list=NULL
  obj1 = glm(binary_ppv_90_rm_trans ~ sex + age + pc1+pc2+pc3+pc4+pc5,data = subset(dat,ancestry==anc),family=binomial(link="logit"))
  for(thispos in allpos){
    print(paste0(HLA," ",which(allpos==thispos),"/",length(allpos)))
    thishaps = as.character(unique(subset(res,pos==thispos)$hap))
    adopted=NULL
    for(thishap in thishaps){
      hlas = as.character(subset(res,hap==thishap)$hla)
      hlas=hlas[hlas %in% colnames(dat)]
    if(length(hlas>0)){
      dat$thishap = rowSums(dat[hlas])
      colnames(dat)[ncol(dat)] =thishap
      adopted = c(adopted,thishap)
    } }
    obj2 = glm(binary_ppv_90_rm_trans ~ as.matrix(dat[dat$ancestry==anc,adopted])+sex+ age + pc1+pc2+pc3+pc4+pc5,data = subset(dat,ancestry==anc),family=binomial(link="logit"))
    chisqtest = anova(obj1,obj2,test="Chisq")
    pval = chisqtest$`Pr(>Chi)`[2]
    deviance = chisqtest$Deviance[2]
    pval_list=c(pval_list,pval)
    deviance_list=c(deviance_list,deviance)
  }
  summary = data.frame("ancestry" =anc,
                       "outcome"=outcome,
                      "position" = allpos,
                       "omnibus_deviance" = deviance_list,
                       "omnibus_pval" = pval_list)
  rd1=rbind.data.frame(rd1,
                       summary)
}

rd1$HLA = HLA
rd1.export = rbind.data.frame(rd1.export,
                              rd1)

}

write.csv(rd1.export,
          row.names = F,
          file="[FILE].csv")

