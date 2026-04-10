library(dplyr)
library(data.table)
library(stringr)

# ALLELE0: reference
# ALLELE1: alternate 

# ALLELE1: effect
# ALLELE2: other

########################
# Running METAL
# 1. edit paths in run_metal.sh
# 2. ./run_metal.sh

#######################
# Find sig alleles
library(dplyr)
library(data.table)
library(tidyr)

# start by identifying independent signals in trans
trans = fread("[FILE].txt") # METAL output
trans = trans %>% separate_wider_delim(MarkerName, 
                                       names=c("chr", "pos", "ref_allele", "alt_allele"), delim=":", too_many = "merge")
dat = trans


dat$"P-value" = as.numeric(dat$"P-value")
dat = dat %>% select(chr, pos, ref_allele, alt_allele, Effect, StdErr, 'P-value', Direction, Allele1, Allele2)
colnames(dat) = c("Chr", "Pos", "RefAllele", "AltAllele", "Effect", "StdErr", "PValue", "Direction", "Allele1", "Allele2")
dat = dat %>% filter(PValue <= 5*10^(-8))
dat = dat %>% mutate(chisq_stat = (Effect / StdErr)^2)
dat = dat %>% arrange(desc(chisq_stat))
dat$Pos = as.numeric(dat$Pos)
dat = dat %>% filter(!(Chr == "chr6" & Pos >=25000000 & Pos <= 35000000)) # remove MHC

selected_variants <- data.frame()
remaining_variants <- dat
window_size <- 1e6

while(nrow(remaining_variants) > 0){
  
  top_variant <- remaining_variants[1, ]
  top_variant_chr = top_variant$Chr
  
  selected_variants <- rbind(selected_variants, top_variant)
  
  lower_bound <- top_variant$Pos - window_size
  upper_bound <- top_variant$Pos + window_size
  
  remaining_variants <- remaining_variants %>%
    filter((Chr==top_variant_chr & (Pos < lower_bound | Pos > upper_bound)) |
             Chr!=top_variant_chr )
}

trans = selected_variants # 137 autosomal 

trans$ancestry = rep("trans_combined", nrow(trans))

# Add other ancestries
eur = fread("[FILE].txt")
trans_seropos = fread("[FILE].txt")
eur_seropos = fread("[FILE].txt")

eur = eur %>% separate_wider_delim(MarkerName, 
                                   names=c("chr", "pos", "ref_allele", "alt_allele"), delim=":", too_many = "merge")
trans_seropos = trans_seropos %>% separate_wider_delim(MarkerName, 
                                                       names=c("chr", "pos", "ref_allele", "alt_allele"), delim=":", too_many = "merge")
eur_seropos = eur_seropos %>% separate_wider_delim(MarkerName, 
                                                   names=c("chr", "pos", "ref_allele", "alt_allele"), delim=":", too_many = "merge")

dat = eur
dat$"P-value" = as.numeric(dat$"P-value")
dat = dat %>% select(chr, pos, ref_allele, alt_allele, Effect, StdErr, 'P-value', Direction, Allele1, Allele2)
colnames(dat) = c("Chr", "Pos", "RefAllele", "AltAllele", "Effect", "StdErr", "PValue", "Direction", "Allele1", "Allele2")
dat = dat %>% filter(PValue <= 5*10^(-8))
dat = dat %>% mutate(chisq_stat = (Effect / StdErr)^2)
dat$Pos = as.numeric(dat$Pos)
dat = dat %>% filter(!(Chr == "chr6" & Pos >=25000000 & Pos <= 35000000))
eur = dat
eur$ancestry = rep("eur_combined", nrow(eur))

dat = trans_seropos
dat$"P-value" = as.numeric(dat$"P-value")
dat = dat %>% select(chr, pos, ref_allele, alt_allele, Effect, StdErr, 'P-value', Direction, Allele1, Allele2)
colnames(dat) = c("Chr", "Pos", "RefAllele", "AltAllele", "Effect", "StdErr", "PValue", "Direction", "Allele1", "Allele2")
dat = dat %>% filter(PValue <= 5*10^(-8))
dat = dat %>% mutate(chisq_stat = (Effect / StdErr)^2)
dat$Pos = as.numeric(dat$Pos)
dat = dat %>% filter(!(Chr == "chr6" & Pos >=25000000 & Pos <= 35000000))
trans_seropos = dat
trans_seropos$ancestry = rep("trans_seropositive", nrow(trans_seropos))

dat = eur_seropos
dat$"P-value" = as.numeric(dat$"P-value")
dat = dat %>% select(chr, pos, ref_allele, alt_allele, Effect, StdErr, 'P-value', Direction, Allele1, Allele2)
colnames(dat) = c("Chr", "Pos", "RefAllele", "AltAllele", "Effect", "StdErr", "PValue", "Direction", "Allele1", "Allele2")
dat = dat %>% filter(PValue <= 5*10^(-8))
dat = dat %>% mutate(chisq_stat = (Effect / StdErr)^2)
dat$Pos = as.numeric(dat$Pos)
dat = dat %>% filter(!(Chr == "chr6" & Pos >=25000000 & Pos <= 35000000))
eur_seropos = dat
eur_seropos$ancestry = rep("eur_seropositive", nrow(eur_seropos))


dat = rbind(eur, trans_seropos, eur_seropos)
dat = dat %>% arrange(desc(chisq_stat))

dat = rbind(trans, dat)
selected_variants <- data.frame()
remaining_variants <- dat
window_size <- 1e6

while(nrow(remaining_variants) > 0){
  
  top_variant <- remaining_variants[1, ]
  top_variant_chr = top_variant$Chr
  
  selected_variants <- rbind(selected_variants, top_variant)
  
  lower_bound <- top_variant$Pos - window_size
  upper_bound <- top_variant$Pos + window_size
  
  remaining_variants <- remaining_variants %>%
    filter((Chr==top_variant_chr & (Pos < lower_bound | Pos > upper_bound)) |
             Chr!=top_variant_chr )
}

selected_variants = selected_variants %>% filter(Chr != "chrX")
write.csv(selected_variants, "[FILE].csv", row.names = F)

# Check 1.5Mb window
remaining_variants = selected_variants
selected_variants = c()
window_size = 1.5e6
while(nrow(remaining_variants) > 0) {
  
  top_variant <- remaining_variants[1, ]
  top_variant_chr = top_variant$Chr
  
  lower_bound <- top_variant$Pos - window_size
  upper_bound <- top_variant$Pos + window_size
  
  remaining_variants = remaining_variants[-1,]
  
  selected_variants = remaining_variants %>%
    filter((Chr==top_variant_chr & (Pos >= lower_bound | Pos <= upper_bound)))
  
  remaining_variants <- remaining_variants %>%
    filter((Chr==top_variant_chr & (Pos < lower_bound | Pos > upper_bound)) |
             Chr!=top_variant_chr )
}

print(selected_variants)



##########################
# Annotate significant alleles: add ancestry-specific MAF, etc

dat = read.csv("sig_alleles.csv")
dat$Chr = gsub("chr", "", dat$Chr) %>% as.numeric()

sum_stat = fread("[FILE].csv")
sum_stat = sum_stat %>% select(CHROM, ALLELE0, ALLELE1, A1FREQ, INFO, N, GENPOS37)
colnames(sum_stat) = c("Chr", "RefAllele", "AltAllele", "EUR_AltFreq", "EUR_INFO", "EUR_N", "Pos") 
sum_stat = sum_stat %>% mutate(EUR_MAF = EUR_AltFreq * (EUR_AltFreq < 0.5) + 
                                 (1-EUR_AltFreq) * (EUR_AltFreq >= 0.5))
dat = dat %>% left_join(sum_stat, by=c("Chr", "Pos", "RefAllele", "AltAllele"))

sum_stat = fread("[FILE].csv")
sum_stat = sum_stat %>% select(CHROM, ALLELE0, ALLELE1, A1FREQ, INFO, N, GENPOS37)
colnames(sum_stat) = c("Chr", "RefAllele", "AltAllele", "AFR_AltFreq", "AFR_INFO", "AFR_N", "Pos") 
sum_stat = sum_stat %>% mutate(AFR_MAF = AFR_AltFreq * (AFR_AltFreq < 0.5) + 
                                 (1-AFR_AltFreq) * (AFR_AltFreq >= 0.5))
dat = dat %>% left_join(sum_stat, by=c("Chr", "Pos", "RefAllele", "AltAllele"))

sum_stat = fread("[FILE].csv")
sum_stat = sum_stat %>% select(CHROM, ALLELE0, ALLELE1, A1FREQ, INFO, N, GENPOS37)
colnames(sum_stat) = c("Chr", "RefAllele", "AltAllele", "AMR_AltFreq", "AMR_INFO", "AMR_N", "Pos") 
sum_stat = sum_stat %>% mutate(AMR_MAF = AMR_AltFreq * (AMR_AltFreq < 0.5) + 
                                 (1-AMR_AltFreq) * (AMR_AltFreq >= 0.5))
dat = dat %>% left_join(sum_stat, by=c("Chr", "Pos", "RefAllele", "AltAllele"))

dat %>% write.csv("sig_alleles.csv", row.names = F)


##########################
# Annotate significant alleles: add closest Ishigaki locus + identify *new* alleles
ish = read.csv("[FILE].csv")
ish_sum = read.csv("[FILE].csv")
ish = cbind(ish, select(ish_sum, Population, Serostatus, allele_freq_EAS, allele_freq_EUR, 
                        allele_freq_AFR, allele_freq_SAS))
ish = ish %>% separate_wider_delim(Variant.ID,
                                   names=c("Chr", "Pos", "RefAllele", "AltAllele"), 
                                   delim=":", too_many = "merge")
ish$Pos = as.numeric(ish$Pos)
ish = ish %>% select(Chr, Pos, RefAllele, AltAllele, Locus.ID, Rs.ID, OR, L95, U95, P.value, Direction, 
                     Population, Serostatus, allele_freq_EAS, allele_freq_EUR, 
                     allele_freq_AFR, allele_freq_SAS)
ish$Pos = as.numeric(ish$Pos)
colnames(ish) = paste0("Ish_nearest_", colnames(ish))
ish$Ish_nearest_Chr = as.integer(ish$Ish_nearest_Chr)

sig_alleles = read.csv("sig_alleles.csv")
sig_alleles$Chr = gsub("chr", "", sig_alleles$Chr) %>% as.numeric()

sig_dt <- as.data.table(select(sig_alleles, Chr, Pos))
ish_dt <- as.data.table(ish)
ish_dt$Pos = ish_dt$Ish_nearest_Pos
setnames(ish_dt, old = c("Ish_nearest_Chr"), new = c("Chr"))
setkey(ish_dt, Chr, Pos)
result <- ish_dt[sig_dt, on = .(Chr, Pos), roll = "nearest"]
sig_alleles = cbind(sig_alleles, select(result, -c(Chr,Pos)))


sig_alleles$dist_nearest_ish = abs(sig_alleles$Ish_nearest_Pos - sig_alleles$Pos)
sig_alleles$new = (sig_alleles$dist_nearest_ish > 1e6) %>% as.numeric()


write.csv(sig_alleles, "results/sig_alleles.csv", row.names = F)



##########################
# Check novel alleles against other previous studies
res = fread("results/sig_allele.csv")
res = filter(res, new == 1)

old = do.call(rbind, lapply(list.files("prev_studies", 
                                       pattern = "\\.tsv$", full.names = TRUE), 
                            read.delim, sep = "\t"))
old = mutate(old, ID = paste0(CHR_ID, "_", CHR_POS))
old = old %>% filter(P.VALUE < 5e-8)
new = select(res, Chr, Pos, Pos38)
new$Chr = as.character(new$Chr)

new = new %>% arrange(Chr, Pos38)
old = old %>% arrange(CHR_ID, CHR_POS)

within_1mb <- function(chr, pos, df2_chr, df2_pos) {
  # Subset df2 for same chromosome
  df2_subset <- df2_chr == chr
  # Get absolute distance
  any(abs(df2_pos[df2_subset] - pos) <= 1e6)
}

# Apply function row-wise
new$within_1mb <- mapply(
  within_1mb,
  new$Chr, new$Pos38,
  MoreArgs = list(df2_chr = old$CHR_ID, df2_pos = old$CHR_POS)
)

table(new$within_1mb)
new = filter(new, within_1mb)

res = fread("results/sig_alleles.csv")
res$new[paste0(res$Chr, ":", res$Pos) %in% paste0(new$Chr, ":", new$Pos)] = 2

write.csv(res, "results/sig_alleles.csv", row.names = F)



##########################
# Fine-mapping
dat = fread("[FILE].txt") # METAL output
dat = dat %>% separate_wider_delim(MarkerName, 
                                   names=c("chr", "pos", "ref_allele_old", "alt_allele_old"), 
                                   delim=":", too_many = "merge")
dat = filter(dat, chr != "chrX")
dat$chr = gsub("chr", "", dat$chr) %>% as.integer()
dat$pos = as.integer(dat$pos)
dat$'P-value' = as.numeric(dat$'P-value')

sig_alleles = read.csv("sig_alleles.csv")

w = 0.04
res = c()
for(i in (1:nrow(sig_alleles))) {
  print(i)
  sig_allele = sig_alleles[i,]
  d1 = filter(dat, chr == sig_allele$Chr, 
              pos >= sig_allele$Pos - 5e5, pos <= sig_allele$Pos + 5e5) %>%
    dplyr::select(chr, pos, ref_allele_old, alt_allele_old, Effect, StdErr, "P-value", Allele1, Allele2) %>%
    mutate(SNP = paste0(chr, "_", pos, "_", ref_allele_old, "_", alt_allele_old)) %>%
    dplyr::select(SNP, Effect, StdErr, "P-value", Allele1, Allele2)
  colnames(d1) <- c("SNP","Beta","SE","P", "Allele1", "Allele2")
  d1$Z <- qnorm(d1$P/2, lower.tail = FALSE)
  d1$abf <- sqrt( d1$SE^2 /( d1$SE^2 + w ) ) * exp( w * d1$Z^2 /( 2 *  (d1$SE^2 + w)) )
  #d1$abf <- sqrt( d1$SE^2 /( d1$SE^2 + w ) ) * exp( w * (d1$Beta)^2 /( 2 * d1$SE^2 * (d1$SE^2 + w)) )
  d1$pip <- d1$abf / sum(d1$abf)
  d1 <- d1[order(d1$pip,decreasing = T),]
  cumsum <- c()
  for(j in 1:nrow(d1)){
    cumsum <- c(cumsum, sum(d1$pip[1:j])) 
  }
  d1$cumsum <- cumsum
  lastrow <- (1:nrow(d1))[cumsum>0.95][1]
  d2 <- d1[1:lastrow, ] %>% cbind(data.frame(key = rep(i, lastrow),
                                             count = rep(lastrow, lastrow)))
  res = rbind(res, d2)
}

sig_alleles = mutate(sig_alleles, key = 1:nrow(sig_alleles))

sig_alleles = dplyr::select(sig_alleles, Chr, Pos, RefAllele, AltAllele, new, key)
colnames(sig_alleles) = c("lead_chr", "lead_pos", "lead_ref_allele", "lead_alt_allele", "new", "key")
res = left_join(res, sig_alleles, by="key")
length(unique(filter(res, new==1)$key))

fwrite(res, "results/[FILE].csv")


