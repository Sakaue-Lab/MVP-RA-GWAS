rm(list=ls())
library(dplyr)
library(stringr)
library(data.table)
maf = readRDS(file="[FILE].RDS")
info = read.csv("[FILE].csv") # HLA dictionary
covar = fread("[FILE].txt")
covar$mvp001_id=as.character(covar$mvp001_id)
covar=data.frame(covar)

hla.two <- fread("[FILE].raw") # HLA AA dosage unphased
hla.two = data.frame(hla.two)
hla.two$IID=as.character(hla.two$IID)

######################################################
######   Implement conditional haplotype test ########
######################################################
rd4.export = NULL
for(HLA in unique(info$gene)){
  
  
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
keep.col = colnames(hla.two)[grepl("HLA_DQB1",colnames(hla.two))]
hla.allele = hla.two[,c("IID",keep.col,"dosage_DQB1")]
hla.allele = hla.allele %>% 
  filter(dosage_DQB1 > 1.95,
         dosage_DQB1 < 2.05)
hla.allele = hla.allele %>%
  select_if(~!all(is.na(.)))
keep = intersect(covar$mvp001_id,hla.allele$IID)
tmp = covar %>% filter(mvp001_id %in% keep)
dat = left_join(tmp,
                hla.allele,
                by=c("mvp001_id"="IID"))
dat = as.data.frame(dat)

info.tmp = info[info$gene=="DQB1",]
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
condpos="37"
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
cond.hap = left_join(cond.hap,
                     hladqb1.hap)


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
hla.allele = hla.two[,c("IID",keep.col,colnames(hla.two)[grepl(paste0("dosage_",HLA),colnames(hla.two))])]
hla.allele = hla.allele[hla.allele[,colnames(hla.two)[grepl(paste0("dosage_",HLA),colnames(hla.two))]]>1.95,]
hla.allele = hla.allele[hla.allele[,colnames(hla.two)[grepl(paste0("dosage_",HLA),colnames(hla.two))]]<2.05,]
hla.allele = hla.allele %>%
  select_if(~!all(is.na(.)))
keep = intersect(covar$mvp001_id,hla.allele$IID)
tmp = covar %>% filter(mvp001_id %in% keep)
dat = left_join(tmp,
                hla.allele,
                by=c("mvp001_id"="IID"))
dat = as.data.frame(dat)

keep = intersect(dat$mvp001_id,cond.hap$mvp001_id)
dat = dat %>% filter(mvp001_id %in% keep)
dat = left_join(dat,
                cond.hap)
dat = as.data.frame(dat)

adopted.prev = colnames(cond.hap)[-1]

rd.cond.drb1=NULL
for(anc in c("AFR")){
pval_list = NULL
deviance_list = NULL
obj1=glm(binary_95_rm~as.matrix(dat[dat$ancestry==anc,adopted.prev])+sex+age+pc1+pc2+pc3+pc4+pc5,data=subset(dat,ancestry==anc),family=binomial(link="logit"))
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
  obj2=glm(binary_95_rm~as.matrix(dat[dat$ancestry==anc,c(adopted.prev,adopted)])+sex+age+pc1+pc2+pc3+pc4+pc5,data=subset(dat,ancestry==anc),family=binomial(link="logit"))
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

write.csv(rd4.export,
          row.names = F,
          file="[FILE].csv")

