library(dplyr)
library(data.table)
library(rcompanion)

p = "5e-8"
prs_path = paste0("[PATH]", p, "/")

prs_files = list.files(path=prs_path,
                       pattern="^chr.*\\.sscore\\.gz$",
                       full.names = T,
                       recursive = T)
scores = lapply(prs_files, fread)

dat_scores = Reduce(`+`, lapply(scores, function(chunk) chunk[[2]]))
dat_scores = data.frame(mvp001_id = scores[[1]][,1], raw_score = dat_scores)
colnames(dat_scores) = c("mvp001_id", "raw_score")

eur_split = fread("[PATH]/[FILE].txt")
afr_split = fread("[PATH]/[FILE].txt")
amr_split = fread("[PATH]/[FILE].txt")
split = eur_split # change depending on ancestry being evaluated

covar = fread("[PATH]/[FILE].txt") # covariate file
covar = select(covar, mvp001_id, ancestry, 
               pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
               age, sex, binary_ppv_90_rm_trans)
covar = covar %>% filter(ancestry == "EUR") # change depending on ancestry being evaluated

dat_scores = dat_scores %>% filter(mvp001_id %in% covar$mvp001_id)
dat_scores = dat_scores %>% mutate(train = mvp001_id %in% split$IID)
dat_scores = dat_scores %>% left_join(covar, by = "mvp001_id")


# Evaluation
dat_scores = filter(dat_scores, !train)

model.null = glm(binary_ppv_90_rm_trans ~ sex + age +
              pc1 + pc2 + pc3 + pc4 + pc5,
            data = dat_scores,
            family = binomial())
model.full = glm(binary_ppv_90_rm_trans ~ raw_score + sex + age +
                   pc1 + pc2 + pc3 + pc4 + pc5,
                 data = dat_scores,
                 family = binomial())

nagelkerke(model.full, null = model.null)$Pseudo.R.squared.for.model.vs.null


######
# Bootstrap CI
library(dplyr)
library(data.table)
library(rcompanion)
library(boot)

nagelkerke_r2 <- function(model, null_model) {
  n = stats::nobs(model)
  ll_full = as.numeric(logLik(model))
  ll_red = as.numeric(logLik(null_model))
  
  r2_CS = 1 - exp((2/n) * (ll_red - ll_full))
  r2_CS / (1 - exp((2/n) * ll_red))
}

boot_nagelkerke <- function(model, null_model, data, R = 1000, conf = 0.95, seed = 123,
                            parallel = c("no", "multicore", "snow"), ncpus=1) {
  parallel = match.arg(parallel)
  frm = formula(model)
  frm_null = formula(null_model)
  fam = family(model)
  
  stat_fun <- function(d, i) {
    d2 = d[i, , drop = F]
    fit_b = glm(frm, data = d2, family = fam)
    fit_b_null = glm(frm_null, data=d2, family=fam)
    
    nagelkerke_r2(fit_b, fit_b_null)
  }
  
  set.seed(seed)
  b = boot::boot(data=data, statistic = stat_fun, R=R, 
                 parallel = if (parallel == "no") "no" else parallel, ncpus = ncpus)
  
  ci = boot::boot.ci(b, type="perc", conf = conf)
  
  list(point = nagelkerke(model, null=null_model),
       boot_object = b,
       ci_perc = if (!is.null(ci$percent)) c(lower = ci$percent[4], upper=ci$percent[5]) else NA)
}

res = boot_nagelkerke(model.full, model.null, dat_scores)
res$point
res$ci_perc

