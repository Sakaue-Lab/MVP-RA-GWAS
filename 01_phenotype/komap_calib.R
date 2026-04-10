library(dplyr)
library(data.table)

labels = read.csv("[PATH]/[FILE].csv", sep="|")

komap2 = read.csv("[PATH]/[FILE].csv", sep="|")
komap = read.table("[PATH]/[FILE].txt", 
                   sep="\t", header = T)

labels = labels %>% 
  left_join(komap, by="mvp001_id") %>% 
  dplyr::select(mvp001_id, ancestry, chartlabel) %>% 
  mutate(chartlabel_num = 1 * (chartlabel != "No"),
         chartlabel_num2 = 1 * (chartlabel == "Yes"))

labels = labels %>% left_join(dplyr::select(komap2, mvp001_id, komap_score), 
                              by = "mvp001_id")
labels = labels %>% filter(!is.na(komap_score))

# Model fitting for AFR patients
labels_afr = labels %>% filter(ancestry == "AFR") # 154 AFR
# komap_pred2 groups probable/possible under "No"
fit.calib.afr2 = glm(chartlabel_num2 ~ komap_score,
                     family="binomial", data=labels_afr)
labels_afr$komap_pred2 = predict(fit.calib.afr2, newdata = labels_afr,
                                 type="response")

# Model fitting for EUR patients
labels_eur = labels %>% filter(ancestry == "EUR") # 54 EUR
# komap_pred2 groups probable/possible under "No"
fit.calib.eur2 = glm(chartlabel_num2 ~ komap_score, 
                    family = "binomial", data = labels_eur)
labels_eur$komap_pred2 = predict(fit.calib.eur2, newdata = labels_eur, 
                                type="response")

# Trans-ancestry model fitting
labels_trans = labels
fit.calib.trans2 = glm(chartlabel_num ~ komap_score, 
                     family = "binomial", data = labels_trans)
labels_trans$komap_pred2 = predict(fit.calib.trans2, newdata = labels_trans, 
                                 type="response")

# CLEAN ROC calcs to get cutoffs -- requires functions in library_roc.R
# EUR
junk = ROC.Est.FUN(labels_eur$chartlabel_num2, labels_eur$komap_pred2, 
                   yy0 = 0.5, fpr0=seq(0,1,0.0001), yes.smooth=T)
junk[1]
junk = matrix(junk[-1], ncol = 6) %>% data.frame()
colnames(junk) = c("cut", "p.pos", "fpr", "tpr", "ppv", "npv")
View(junk)
cut_eur = NA # fill

# AFR
junk = ROC.Est.FUN(labels_afr$chartlabel_num2, labels_afr$komap_pred2, 
                   yy0 = 0.5, fpr0=seq(0,1,0.0001), yes.smooth=T)
junk[1]
junk = matrix(junk[-1], ncol = 6) %>% data.frame()
colnames(junk) = c("cut", "p.pos", "fpr", "tpr", "ppv", "npv")
View(junk)
cut_afr = NA # fill

# Trans-ancestry
junk = ROC.Est.FUN(labels_trans$chartlabel_num2, labels_trans$komap_pred2, 
                   yy0 = 0.5, fpr0=seq(0,1,0.0001), yes.smooth=T)
junk[1]
junk = matrix(junk[-1], ncol = 6) %>% data.frame()
colnames(junk) = c("cut", "p.pos", "fpr", "tpr", "ppv", "npv")
View(junk)
cut_trans = NA # fill

komap2 = komap2 %>% left_join(dplyr::select(komap, mvp001_id, ancestry), by = "mvp001_id")
komap2_afr$komap_pred2 = predict(fit.calib.afr2, newdata = komap2_afr, 
                                 type = "response")
komap2_eur$komap_pred2 = predict(fit.calib.eur2, newdata = komap2_eur, 
                                 type = "response")
komap2_trans$komap_pred2 = predict(fit.calib.trans2, newdata = komap2_trans,
                                   type = "response")
komap2_comb = rbind(komap2_afr, komap2_eur, komap2_trans)

# Change select here depending on komap_pred or komap_pred2
komap = komap %>% left_join(dplyr::select(komap2_comb, mvp001_id, komap_pred2), 
                            by = "mvp001_id")

# add phecode counts
pheno = fread("[PATH]/[FILE].csv") %>%
  dplyr::select(mvp001_id, phe_714_1)
pheno$mvp001_id = as.double(pheno$mvp001_id)

komap = komap %>% left_join(pheno, by = "mvp001_id")

covar = komap
covar = covar %>% mutate(phe_714_1_binary2 = 1*(phe_714_1 >= 2))

# Phenotypes based on PPV cutoffs
covar = covar %>% mutate(binary_ppv_90 = 1*(ancestry == "EUR")*(komap_pred2 > cut_eur) +
                           1*(ancestry=="AFR")*(komap_pred2>cut_afr) +
                           1*(!(ancestry %in% c("EUR", "AFR")))*(komap_pred2>cut_trans))


covar$binary_ppv_90[is.na(covar$binary_ppv_90)] = 0
covar$binary_ppv_95[is.na(covar$binary_ppv_95)] = 0

# Remove intermediate patients
rm_ind = (covar$phe_714_1 == 1) | is.na(covar$phe_714_1)
covar$binary_95_rm = covar$binary_95
covar$phe_714_1_binary2_rm = covar$phe_714_1_binary2
covar$binary_ppv_90_rm = covar$binary_ppv_90
covar$binary_ppv_95_rm = covar$binary_ppv_95

covar$binary_95_rm[rm_ind] = NA
covar$phe_714_1_binary2_rm[rm_ind] = NA
covar$binary_ppv_90_rm[rm_ind] = NA
covar$binary_ppv_95_rm[rm_ind] = NA


fwrite(covar, "[FILE].txt", row.names = F, sep = "\t")
