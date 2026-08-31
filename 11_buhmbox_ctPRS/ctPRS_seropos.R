#!/usr/bin/env Rscript
# Cell-type-specific PRS (ctPRS) vs. seropositivity among RA cases.
# For each cell type we build a ctPRS from its cell-type-specific PRS variants in
# the held-out MVP participants (30% not used to train the overall PRS), then fit
#   seropositive ~ ctPRS + overall_PRS + covariates
# to estimate the additional contribution of the ctPRS on seropositivity (a known
# axis of RA heterogeneity), conditional on the overall genome-wide PRS.
#
# Output: prs_seropos_res.csv  (cell_type, variable[prs_ct/prs_full], beta, se, z, p)
suppressPackageStartupMessages({library(data.table)})

CELLTYPES <- c("Bcell", "DC", "endothelial", "fibroblast", "monocyte", "Tcell", "unioncells")

# held-out MVP RA cases: seropositivity + covariates + overall PRS
pheno <- fread("stats/heldout_cases.pheno.tsv")   # IID seropos age sex PC1..PC5 batch prs_full
# clumped PRS weights (varID, risk allele, weight) and per-sample dosages (held-out)
w    <- fread("stats/prs_best_variants.risk_alleles.txt")          # varID risk_allele beta
dos  <- fread("stats/heldout_cases.dosage.txt.gz")                 # rows varID ; cols IID (risk-allele dosage)

ct_prs <- function(vars){
  v  <- intersect(vars, w$varID)
  W  <- w[match(v, varID)]
  D  <- as.matrix(dos[match(v, varID), ..pheno$IID])              # variants x samples
  as.numeric(crossprod(D, W$beta))                                # weighted sum per sample
}
scale1 <- function(x) as.numeric(scale(x))

res <- rbindlist(lapply(CELLTYPES, function(ct){
  vars <- readLines(sprintf("stats/celltype_variants/%s.variants.txt", ct))
  d <- copy(pheno)
  d[, prs_ct   := scale1(ct_prs(vars))]
  d[, prs_full := scale1(prs_full)]
  fit <- glm(seropos ~ prs_ct + prs_full + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + batch,
             data = d, family = binomial())
  co <- summary(fit)$coefficients
  rbindlist(lapply(c("prs_ct", "prs_full"), function(v)
    data.table(cell_type = ct, variable = v,
               beta = co[v, "Estimate"], se = co[v, "Std. Error"],
               z = co[v, "z value"], p = co[v, "Pr(>|z|)"])))
}))

fwrite(res, "prs_seropos_res.csv")
cat("wrote prs_seropos_res.csv\n")
print(res[variable == "prs_ct", .(cell_type, beta = round(beta, 3), p = signif(p, 3))])
