## E7b - Table 4 re-run on log2(TPM + 1) scores (E7a). Same comparison set and procedure as
## revision/scripts/A2_roc_full.R (pROC; stratified bootstrap CI, 2,000 resamples; DeLong vs
## NetImmune; BH across the 13 comparisons), except that 'TIDE' is replaced by what the
## locked script actually computed, the CAF/TAM exclusion proxy (inverted), and IPS is the
## Charoentong algorithm (raw aggregate score). Extra variants are reported outside the BH set.
suppressMessages(library(pROC)); set.seed(42)
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
d <- read.csv(file.path(O, "E7_scores_log2.csv"), row.names = 1)
y <- d$Responder
M <- list(`MSI/EBV` = d$MSI_EBV, IFN = d$IFN, Immune = d$Immune, GEP18 = d$GEP18,
          `PD-L1` = d$PDL1, APM = d$APM, NetImmune = d$NetImmune, T_CELL = d$T_CELL,
          TAM_SUP = d$TAM_SUP, `MSI-H` = d$MSI_H, TLS_B = d$TLS_B,
          `ExclusionProxy(neg)` = -d$ExclusionProxy, `Barrier(neg)` = -d$Barrier, IPS = d$IPS_AZ)
EXTRA <- list(`IPS (0-10 discrete)` = d$IPS, `GEP18 (z-mean)` = d$GEP18_zmean, `CAF_ECM(neg)` = -d$CAF_ECM)
auc_row <- function(nm, s) {
  r <- roc(y, s, quiet = TRUE, direction = "<")
  ci <- ci.auc(r, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
  list(r = r, row = data.frame(Marker = nm, AUC = as.numeric(auc(r)), CI_low = ci[1], CI_high = ci[3]))
}
A <- lapply(names(M), function(nm) auc_row(nm, M[[nm]])); names(A) <- names(M)
res <- do.call(rbind, lapply(A, `[[`, "row")); res <- res[order(-res$AUC), ]
dl <- do.call(rbind, lapply(setdiff(names(A), "NetImmune"), function(nm) {
  t <- roc.test(A$NetImmune$r, A[[nm]]$r, method = "delong")
  data.frame(Comparator = nm, AUC_comp = as.numeric(auc(A[[nm]]$r)),
             Delta = as.numeric(auc(A$NetImmune$r)) - as.numeric(auc(A[[nm]]$r)),
             Z = as.numeric(t$statistic), p = t$p.value)
}))
dl <- dl[order(dl$p), ]; dl$q <- p.adjust(dl$p, "BH")
ex <- do.call(rbind, lapply(names(EXTRA), function(nm) {
  a <- auc_row(nm, EXTRA[[nm]]); t <- roc.test(A$NetImmune$r, a$r, method = "delong")
  cbind(a$row, p_vs_NetImmune = t$p.value)
}))
write.csv(res, file.path(O, "E7_auc_ci_log2.csv"), row.names = FALSE)
write.csv(dl, file.path(O, "E7_delong_log2.csv"), row.names = FALSE)
write.csv(ex, file.path(O, "E7_auc_extra_log2.csv"), row.names = FALSE)
cat("=== AUC (log2), 2,000 stratified bootstrap ===\n"); print(format(res, digits = 3), row.names = FALSE)
cat("\n=== DeLong vs NetImmune, BH over 13 ===\n"); print(format(dl, digits = 3), row.names = FALSE)
cat("\n=== extra variants (not in BH set) ===\n"); print(format(ex, digits = 3), row.names = FALSE)
