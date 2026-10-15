## E1 - formal head-to-head comparison of prognostic scores.
##
## v13 compared scores by their meta-analysis p-values (NetImmune p=1.8e-7 vs CD8A p=0.094;
## NetImmune HR 0.74 vs Immune HR 0.84). A difference in significance is not a significant
## difference. Here each pair is compared directly: Harrell's C of each score, with its
## direction fixed a priori (no fitting), and the paired difference bootstrapped over
## patients. Per-cohort differences are pooled by inverse-variance weighting.
## ESTIMATE is the signature-level implementation used in B4 (ESTsig_diff); the official
## ssGSEA implementation fails on the GSE84437 array and is not used.
suppressMessages(library(survival))
set.seed(2026)
R <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
NB <- 2000
TAGS <- c("TCGA", "GSE84437", "GSE62254")
## risk orientation fixed before looking at outcome: high immune = low risk
RISK <- list(NetImmune = function(d) -d$NetImmune, Immune = function(d) -d$Immune,
             Barrier = function(d) d$Barrier, CD8A = function(d) -d$CD8A,
             ESTIMATE = function(d) -d$ESTsig_diff, NetImmune_ssGSEA = function(d) -d$NetImmune_ss,
             CAF_ECM = function(d) d$CAF_ECM,
             ## exploratory, defined after E4 showed TAM_SUP carries no prognostic signal
             Immune_minus_CAF = function(d) -(d$Immune - d$CAF_ECM))
PAIRS <- list(c("NetImmune", "CD8A"), c("NetImmune", "Immune"), c("NetImmune", "ESTIMATE"),
              c("NetImmune", "NetImmune_ssGSEA"), c("NetImmune", "Barrier"),
              c("NetImmune", "CAF_ECM"), c("Immune_minus_CAF", "NetImmune"))

hc <- function(tt, ee, r) concordance(Surv(tt, ee) ~ r, reverse = TRUE)$concordance

rows <- list(); cs <- list()
for (tag in TAGS) {
  d <- readRDS(file.path(R, "B2_scores.rds"))[[tag]]
  d <- d[!is.na(d$time) & d$time > 0 & !is.na(d$event), ]
  M <- sapply(RISK, function(f) f(d))
  C0 <- apply(M, 2, function(r) hc(d$time, d$event, r))
  B <- t(replicate(NB, {
    i <- sample(nrow(d), replace = TRUE)
    apply(M[i, ], 2, function(r) hc(d$time[i], d$event[i], r))
  }))
  cs[[tag]] <- data.frame(cohort = tag, score = names(C0), C = C0,
                          lo = apply(B, 2, quantile, 0.025), hi = apply(B, 2, quantile, 0.975))
  for (p in PAIRS) {
    dd <- B[, p[1]] - B[, p[2]]
    rows[[length(rows) + 1]] <- data.frame(cohort = tag, a = p[1], b = p[2], n = nrow(d),
      dC = C0[p[1]] - C0[p[2]], se = sd(dd), lo = quantile(dd, 0.025), hi = quantile(dd, 0.975),
      p = 2 * min(mean(dd <= 0), mean(dd >= 0)))
  }
  cat(tag, "done\n")
}
cidx <- do.call(rbind, cs); rownames(cidx) <- NULL
res <- do.call(rbind, rows); rownames(res) <- NULL
pool <- do.call(rbind, lapply(PAIRS, function(p) {
  s <- res[res$a == p[1] & res$b == p[2], ]
  w <- 1 / s$se^2; est <- sum(w * s$dC) / sum(w); se <- sqrt(1 / sum(w))
  data.frame(cohort = "pooled", a = p[1], b = p[2], n = sum(s$n), dC = est, se = se,
             lo = est - 1.96 * se, hi = est + 1.96 * se, p = 2 * pnorm(-abs(est / se)))
}))
res <- rbind(res, pool)
write.csv(cidx, file.path(O, "E1_cindex_by_score.csv"), row.names = FALSE)
write.csv(res, file.path(O, "E1_paired_dC.csv"), row.names = FALSE)
print(format(cidx, digits = 3)); cat("\n")
print(format(res, digits = 3))
