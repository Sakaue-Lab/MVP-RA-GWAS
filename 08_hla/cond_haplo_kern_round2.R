rm(list=ls())
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
if(length(args) != 6) stop("Usage: Rscript cond_haplo_kern_round2.R <dictionary.csv> <covariates.txt> <dosages.raw> <exclude.csv> <variants.pvar> <output.csv>")
info = read.csv(args[1]) # HLA dictionary
covar = fread(args[2])
covar$id=as.character(covar$id)
covar=data.frame(covar)

hla.two <- fread(args[3]) # HLA AA dosage unphased
hla.two = data.frame(hla.two)
hla.two$IID=as.character(hla.two$IID)

exclude = fread(args[4]) # exclude related
covar = covar %>% filter(!(id %in% exclude$IID))

### Post-Impute QC
r2 = fread(args[5])
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
rd2.export = NULL
for(HLA in "DRB1"){
  
  
#########################################################################
#####     Prepare DRB1 Data for Conditional Haplotype Testing     #######
#########################################################################
keep.col = colnames(hla.two)[grepl(paste0("HLA_","DRB1"),colnames(hla.two))]
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

keep = intersect(covar$id,hla.allele$IID)
tmp = covar %>% filter(id %in% keep)
dat = left_join(tmp,
                hla.allele,
                by=c("id"="IID"))
dat = as.data.frame(dat)

HLA="DRB1"
info = info[info$gene==HLA,]
info$tag=paste0(info$pos,":",info$AA)
allpos = sort( unique(info$pos) )
res = data.frame()
for (pos in allpos){
  y = info[info$pos %in% c(pos), ]
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


### pulled from old
res.prev = res
add_pos="71"
condpos1 = "11"
condpos2 = "71"
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



allpos <- sort( unique(info$pos) )
allpos <- setdiff(allpos, condpos1) #other positions to analyse
allpos <- setdiff(allpos, condpos2)
res <- data.frame()
for( pos in allpos ){
  x <- info[ info$pos %in% c(condpos1), ] #haplotype information at the position 1 I want to condition on
  y <- info[ info$pos %in% c(condpos2), ] #haplotype information at the position 2 I want to condition on
  z <- info[ info$pos %in% c(pos), ] #haplotype information at the position I want to analyze
  for( i in 1:length( unique( x$tag )) ){
    for( j in 1:length( unique( y$tag )) ){
      for( k in 1:length( unique( z$tag )) ){
        xtag <- unique( x$tag )[ i ]
        ytag <- unique( y$tag )[ j ]
        ztag <- unique( z$tag )[ k ]
        hap <- paste0(xtag,"_",ytag)
        hap <- paste0(hap,"_",ztag)
        x_4d <- subset(x, tag == xtag )$hla
        y_4d <- subset(y, tag == ytag )$hla
        z_4d <- subset(z, tag == ztag )$hla
        hap_4d <- intersect( x_4d, y_4d )
        hap_4d <- intersect( hap_4d, z_4d )
        if( length(hap_4d) > 0 ){
          out <- data.frame(hap, hla = hap_4d, pos )
          res <- rbind(res, out)
        }
      }
    }
  }
}

drb1.hap = dat[,c("id",adopted.prev)]


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

keep = intersect(covar$id,hla.allele$IID)
tmp = covar %>% filter(id %in% keep)
dat = left_join(tmp,
                hla.allele,
                by=c("id"="IID"))
dat = as.data.frame(dat)
dat = left_join(dat,drb1.hap)


rd.cond.drb1=NULL
for(anc in c("EUR","AFR")){
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
rd2.export = rbind.data.frame(rd2.export,
                              rd.cond.drb1)
}

dir.create(dirname(args[6]), recursive=TRUE, showWarnings=FALSE)
write.csv(rd2.export,
          row.names = F,
          file=args[6])
