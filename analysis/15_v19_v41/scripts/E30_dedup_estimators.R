## E30 Q2 (gastric part; disclosure, no decision): every gastric pooled estimate reported in the main text from k <= 7
## cohorts - benchmark three, locked first three, extension four, locked seven - under DerSimonian-Laird (DL), REML,
## DL with Hartung-Knapp (DL-HK), REML with Hartung-Knapp (REML-HK) and the manuscript rule (DL estimate; wider of the
## DL-HK and DL intervals; larger P), with the GSE26253 duplicates (published per-cohort estimates) and without them
## (E30 primary: 8 arrays removed; and the 3-array sensitivity set).
## Per-cohort log HRs are those the published code pooled: E20_estimates_used_for_pooling.csv (benchmark from E8/E5,
## locked first three from E15, added four from E20, pre-registered duplicate exclusions), and the same file from the
## E30 re-runs (results/E30_dedup_runs/<mode>/). Output: results/E30_estimators_gastric.csv
suppressMessages(library(metafor))
A <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
SRC <- list(with_duplicates = file.path(A, "E20_estimates_used_for_pooling.csv"),
            without_duplicates_8 = file.path(A, "E30_dedup_runs/primary/E20_estimates_used_for_pooling.csv"),
            without_duplicates_3_sens = file.path(A, "E30_dedup_runs/sens3/E20_estimates_used_for_pooling.csv"))
BM <- c("TCGA", "GSE84437", "GSE62254"); L3 <- c("GSE15459", "GSE34942", "GSE26253")
NEW <- c("GSE13861", "GSE26901", "GSE26899", "GSE29272")
SETS <- list("benchmark 3" = BM, "locked first 3" = L3, "extension 4" = NEW, "locked 7" = c(L3, NEW))
SCORES <- c("NetImmune", "Immune", "Barrier", "CAF_ECM", "TAM_SUP", "IFN", "APM", "TLS_B", "T_CELL",
            "CXCL9_SPP1", "ESTIMATE", "CD8A", "C1Q_arm", "SPP1_arm")
## where the manuscript-rule value is reported in the main text (paragraph [n] / table [Tn]) of the v31 manuscript
MAIN <- c(
  "benchmark 3|NetImmune|unadjusted" = "[115] 0.74 (0.58-0.95; P = 0.036); [T1] Pooled; [T2]",
  "benchmark 3|CAF_ECM|unadjusted" = "[122] 1.29 (1.18-1.41; P = 0.001); [T2]",
  "benchmark 3|TAM_SUP|unadjusted" = "[122] 1.02 (0.81-1.27); [T2]",
  "benchmark 3|IFN|unadjusted" = "[122] 0.82 (0.60-1.12, P = 0.11); [T2]",
  "benchmark 3|Barrier|unadjusted" = "[123] 1.19 (1.05-1.35); [T2]",
  "benchmark 3|Immune|unadjusted" = "[T2] 0.84 (0.61-1.17)",
  "benchmark 3|CXCL9_SPP1|unadjusted" = "[T2] 0.84 (0.67-1.06)",
  "benchmark 3|C1Q_arm|unadjusted" = "[T2] 0.99 (0.80-1.22)",
  "benchmark 3|SPP1_arm|unadjusted" = "[T2] 1.04 (0.89-1.22)",
  "locked first 3|NetImmune|unadjusted" = "[14] P = 0.057; [120] 0.73 (0.53-1.02; P = 0.057); [154]; [272]; [281]; [285]; [T1] Pooled, first three",
  "locked first 3|NetImmune|adjusted" = "[T1] Pooled, first three: 0.85 (0.38-1.87); 0.46",
  "extension 4|NetImmune|unadjusted" = "[T1] Pooled, added four: 0.67 (0.53-0.85); 0.012",
  "extension 4|NetImmune|adjusted" = "[T1] Pooled, added four: 0.71 (0.59-0.85); 0.0098",
  "locked 7|NetImmune|unadjusted" = "[15] 0.70 (0.62-0.80); [120] 0.70 (0.62-0.80; P = 4.3e-4; I2 = 7%); [148]; [T1] Pooled, all seven locked",
  "locked 7|NetImmune|adjusted" = "[T1] Pooled, all seven locked: 0.75 (0.62-0.92); 0.012",
  "locked 7|CAF_ECM|unadjusted" = "[120] 1.34 (1.08-1.65)")
fit <- function(b, se) {
  m <- list(DL = rma(yi = b, sei = se, method = "DL"), REML = rma(yi = b, sei = se, method = "REML"),
            `DL-HK` = rma(yi = b, sei = se, method = "DL", test = "knha"),
            `REML-HK` = rma(yi = b, sei = se, method = "REML", test = "knha", control = list(maxiter = 1000)))
  out <- do.call(rbind, lapply(names(m), function(e) data.frame(estimator = e, logHR = m[[e]]$b[1], lo = m[[e]]$ci.lb,
           hi = m[[e]]$ci.ub, p = m[[e]]$pval, tau2 = m[[e]]$tau2, I2 = m[[e]]$I2, interval_used = e)))
  dl <- m$DL; hk <- m$`DL-HK`; w <- (hk$ci.ub - hk$ci.lb) >= (dl$ci.ub - dl$ci.lb)     # E15/E20 pool()
  rbind(out, data.frame(estimator = "manuscript rule", logHR = dl$b[1], lo = if (w) hk$ci.lb else dl$ci.lb,
                        hi = if (w) hk$ci.ub else dl$ci.ub, p = max(dl$pval, hk$pval), tau2 = dl$tau2, I2 = dl$I2,
                        interval_used = if (w) "DL-HK" else "DL"))
}
res <- list()
for (dv in names(SRC)) {
  E <- read.csv(SRC[[dv]], stringsAsFactors = FALSE)
  E <- E[grepl("^pre-registered", E$duplicate_rule) & E$group %in% c("benchmark", "locked (first 3)", "locked (added)"), ]
  for (st in names(SETS)) for (s in SCORES) for (m in c("unadjusted", "adjusted")) {
    a <- E[E$score == s & E$model == m & E$cohort %in% SETS[[st]], ]
    if (nrow(a) < length(SETS[[st]])) next                       # every cohort of the set must have the estimate
    a <- a[match(SETS[[st]], a$cohort), ]
    f <- fit(a$beta, a$se)
    key <- paste(st, s, m, sep = "|")
    res[[length(res) + 1]] <- data.frame(set = st, score = s, model = m, k = nrow(a), cohorts = paste(a$cohort, collapse = " "),
      duplicates = dv, f, HR = exp(f$logHR), HR_lo = exp(f$lo), HR_hi = exp(f$hi),
      reported_in_main_text = if (key %in% names(MAIN)) MAIN[[key]] else "", stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, res); rownames(R) <- NULL
R$t_crit_HK <- ifelse(grepl("HK", R$estimator), qt(0.975, R$k - 1), NA)
R <- R[, c("set", "score", "model", "k", "cohorts", "duplicates", "estimator", "interval_used", "HR", "HR_lo", "HR_hi", "p",
           "logHR", "lo", "hi", "tau2", "I2", "t_crit_HK", "reported_in_main_text")]
names(R)[names(R) %in% c("lo", "hi")] <- c("logHR_lo", "logHR_hi")
write.csv(R, file.path(A, "E30_estimators_gastric.csv"), row.names = FALSE)
## check: manuscript rule, with duplicates, reproduces the published pooled values (E15_pooled / E20_pooled)
P20 <- read.csv(file.path(A, "E20_pooled.csv"), stringsAsFactors = FALSE)
P20 <- P20[grepl("^pre-registered", P20$duplicate_rule), ]
map <- c("locked first 3" = "locked (first 3), E15 reference", "extension 4" = "locked (added 4)", "locked 7" = "locked 7 (primary)",
         "benchmark 3" = "benchmark 3 (reference)")
q <- R[R$estimator == "manuscript rule" & R$duplicates == "with_duplicates", ]
q$pub_HR <- mapply(function(st, s, m) { r <- P20[P20$set == map[[st]] & P20$score == s & P20$model == m, ]; if (nrow(r)) r$HR[1] else NA }, q$set, q$score, q$model)
q$pub_p <- mapply(function(st, s, m) { r <- P20[P20$set == map[[st]] & P20$score == s & P20$model == m, ]; if (nrow(r)) r$p[1] else NA }, q$set, q$score, q$model)
cat("manuscript rule vs E20_pooled.csv (with duplicates): max |HR diff|", signif(max(abs(q$HR - q$pub_HR), na.rm = TRUE), 3),
    "; max |P diff|", signif(max(abs(q$p - q$pub_p), na.rm = TRUE), 3), "over", sum(!is.na(q$pub_HR)), "estimates\n")
options(width = 220)
show <- R[R$reported_in_main_text != "" & R$duplicates != "without_duplicates_3_sens", ]
print(transform(show[, c("set", "score", "model", "duplicates", "estimator", "interval_used", "HR", "HR_lo", "HR_hi", "p", "tau2")],
                HR = round(HR, 3), HR_lo = round(HR_lo, 3), HR_hi = round(HR_hi, 3), p = signif(p, 2), tau2 = signif(tau2, 2)), row.names = FALSE)
