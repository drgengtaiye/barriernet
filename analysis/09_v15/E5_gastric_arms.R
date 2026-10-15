## E5 (gastric) - prognosis of the two TAM_SUP arms and the CXCL9:SPP1 polarity score in the
## three gastric bulk cohorts, alongside CAF_ECM and TAM_SUP for reference.
## Peer-review addition; definitions and expectations fixed before outcomes were computed:
##   C1Q_arm  = C1QA, C1QB, C1QC                 - tracks immune engagement: protective or null (HR <= 1)
##   SPP1_arm = SPP1, TREM2, MARCO, APOE, LGALS3 - adverse or null (HR >= 1)
##   CXCL9_SPP1 = z(CXCL9) - z(SPP1)             - protective (HR < 1), reference
## Scoring as B2_scores.R: genes row-z within cohort, module = mean z (>= 3 genes present).
## Cox HR per SD, univariable; DerSimonian-Laird pooling.
suppressMessages({library(survival); library(metafor); library(matrixStats)})
R <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
MOD <- list(CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
            TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
            C1Q_arm = c("C1QA","C1QB","C1QC"), SPP1_arm = c("SPP1","TREM2","MARCO","APOE","LGALS3"))
zrow <- function(m) { s <- rowSds(m, na.rm = TRUE); s[s == 0 | is.na(s)] <- 1; (m - rowMeans(m, na.rm = TRUE)) / s }
S2 <- readRDS(file.path(R, "B2_scores.rds"))
rows <- list()
for (t in c("TCGA", "GSE84437", "GSE62254")) {
  d <- readRDS(file.path(R, paste0("B1_", t, ".rds"))); Z <- zrow(d$expr); s <- S2[[t]]
  idx <- match(s$sample, colnames(Z)); ok <- !is.na(idx) & !is.na(s$time) & s$time > 0 & !is.na(s$event)
  Z <- Z[, idx[ok]]; tt <- s$time[ok]; ee <- s$event[ok]
  sc <- lapply(MOD, function(G) { g <- intersect(G, rownames(Z)); if (length(g) >= 3) colMeans(Z[g, , drop = FALSE]) else NA })
  sc$CXCL9_SPP1 <- Z["CXCL9", ] - Z["SPP1", ]
  for (k in names(sc)) {
    f <- coxph(Surv(tt, ee) ~ scale(sc[[k]])[, 1])
    rows[[length(rows) + 1]] <- data.frame(cohort = t, score = k, n = length(tt),
      genes = if (k %in% names(MOD)) length(intersect(MOD[[k]], rownames(Z))) else 2,
      b = unname(coef(f)), se = sqrt(vcov(f)[1, 1]), p = summary(f)$coefficients[1, 5])
  }
}
per <- do.call(rbind, rows)
pool <- do.call(rbind, lapply(unique(per$score), function(k) { u <- per[per$score == k, ]
  m <- rma(yi = u$b, sei = u$se, method = "DL")
  data.frame(score = k, HR = exp(m$b[1]), lo = exp(m$ci.lb), hi = exp(m$ci.ub), p = m$pval, I2 = m$I2,
             per_cohort_HR = paste(sprintf("%.2f", exp(u$b)), collapse = " / ")) }))
write.csv(per, file.path(O, "E5_gastric_arms_per_cohort.csv"), row.names = FALSE)
write.csv(pool, file.path(O, "E5_gastric_arms_pooled.csv"), row.names = FALSE)
print(format(pool, digits = 3), row.names = FALSE)
