## E4 - prognostic association of each of the six modules (HR per 1 SD, DerSimonian-Laird).
## Used to show how the two Barrier components behave separately, since in the ICI cohort
## they pull in opposite directions (CAF_ECM lower, TAM_SUP higher in responders).
suppressMessages({library(survival); library(metafor)})
R <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
res <- readRDS(file.path(R, "B2_scores.rds"))
MODS <- c("CAF_ECM", "TAM_SUP", "APM", "IFN", "TLS_B", "T_CELL")
uni <- do.call(rbind, lapply(names(res), function(t) {
  d <- res[[t]]; d <- d[!is.na(d$time) & d$time > 0, ]
  do.call(rbind, lapply(MODS, function(m) {
    f <- coxph(Surv(time, event) ~ scale(d[[m]])[, 1], data = d)
    data.frame(cohort = t, module = m, beta = coef(f), se = sqrt(vcov(f)[1, 1]))
  }))
}))
meta <- do.call(rbind, lapply(MODS, function(m) {
  u <- uni[uni$module == m, ]; r <- rma(yi = u$beta, sei = u$se, method = "DL")
  data.frame(module = m, HR = exp(r$b[1]), lo = exp(r$ci.lb), hi = exp(r$ci.ub), p = r$pval, I2 = r$I2)
}))
## correlation of the two Barrier components in bulk
rr <- sapply(res, function(d) cor(d$CAF_ECM, d$TAM_SUP, use = "pair"))
write.csv(meta, file.path(O, "E4_module_meta.csv"), row.names = FALSE)
print(format(meta, digits = 3)); cat("\nr(CAF_ECM, TAM_SUP) by cohort:\n"); print(round(rr, 2))
