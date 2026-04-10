library(dplyr)
library(data.table)
library(tidyr)


# Prep for clumping
dat = fread("[PATH]/[FILE].csv")

dat = dat %>% 
  filter(!is.na(Chr),
         nchar(ID)>0) %>%
  select(ID)
colnames(dat)="ID"
write.csv(dat,
          "[PATH]/snp_list_trans.csv",
          row.names = F,
          quote = F)

dat = rbind(tmp1, tmp2)
dat = dat %>% 
  select(ID,Effect,PValue,StdErr,Allele1,Allele2) %>%
  filter(nchar(ID)>0)
colnames(dat)=c("ID","Beta","P","SE","Allele1","Allele2")
write.csv(dat,
          "[PATH]/gwas_trans.csv",
          row.names = F,
          quote = F)

fix = read.csv("[PATH]/snp_list_trans.csv")
colnames(fix) = "ID"
fix %>% write.csv("[PATH]/snp_list_trans.csv",
                  row.names = F,
                  quote = F)

fix = fread("[PATH]/gwas_trans.csv")
colnames(fix) = c("ID","Beta","P","SE","Allele1","Allele2")
fix %>% write.csv("[PATH]/gwas_trans.csv",
                  row.names = F,
                  quote = F)

#########
# Remove MHC
library(dplyr)
library(data.table)
library(tidyr)

snplist = fread("[PATH]/snp_list_trans.csv")
snplist = snplist %>% separate_wider_delim(ID,
                                           names=c("Chr", "Pos", "Ref", "Alt"),
                                           delim=":",
                                           too_many="merge",
                                           cols_remove=F)
snplist$Pos = as.numeric(snplist$Pos)
snplist = snplist %>% filter(!(Chr=="chr6" & Pos >= 25e6 & Pos <= 35e6))
snplist = snplist %>% select(ID)
snplist %>% write.csv("[PATH]/snp_list_trans_MHCrm.csv",
                      row.names = F,
                      quote = F)

gwas = fread("[PATH]/gwas_trans.csv")
gwas = gwas %>% filter(ID %in% snplist$ID)
gwas %>% write.csv("[PATH]/gwas_trans_MHCrm.csv",
                   row.names = F,
                   quote = F)

#########
# Downsample IDs for clumping
set.seed(123)
eur_split = fread("[PATH]/[FILE].txt") # IDs
eur_split = eur_split[sample(1:nrow(eur_split), 10000),]

write.table(eur_split, "clump_eur_keep.txt", 
            quote=F, sep="\t", row.names=F, col.names=F)

