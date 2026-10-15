## E10 (step 4) - cross-trial synthesis of the score x checkpoint-arm interaction.
## Randomised comparisons (checkpoint-containing arm vs a comparator without it):
##   IMmotion150  atezolizumab +/- bevacizumab vs sunitinib  - PFS, interaction HR (E5_trials.csv)
##   I-SPY2       pembrolizumab + paclitaxel vs paclitaxel    - pCR, interaction OR (E10, HR-adjusted)
##   I-SPY2       durvalumab + olaparib + paclitaxel vs pac.  - pCR, interaction OR (E10, HR-adjusted)
## Sensitivity: + PRINCE (nivolumab-containing vs sotigalimab, all chemotherapy; both arms are
## immunotherapy), I-SPY2 without HR adjustment, and the two I-SPY2 comparisons alone (one scale).
## IMmotion150 atezolizumab monotherapy vs sunitinib is not added: it shares the sunitinib arm.
## Harmonisation: benefit-oriented log ratio per SD, y = -log(interaction HR) for survival and
## +log(interaction OR) for pCR, so y > 0 always means the score predicts more benefit from the
## checkpoint arm. HR and OR are different scales, so the pooled y is a direction-and-size
## summary, not an effect size; it is reported with a scale-free companion, the unweighted
## Stouffer combination of z = y / SE, and the count of comparisons with y > 0.
## Random effects: REML with the Hartung-Knapp adjustment (k - 1 df), as in the manuscript;
## fixed-effect inverse-variance estimates alongside.
suppressMessages(library(metafor))
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
E5 <- read.csv(file.path(O, "E5_trials.csv"))
E10 <- read.csv(file.path(O, "E10_ispy2_interactions.csv"))
ORDER <- c("CAF_ECM", "TAM_SUP", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1", "IFN", "APM", "TLS_B", "T_CELL",
           "Barrier", "Immune", "NetImmune")

surv <- function(pat, lab) {
  d <- E5[startsWith(E5$trial, pat), ]
  data.frame(trial = lab, score = d$score, endpoint = "PFS (HR)", y = -log(d$int_HR),
             se = (log(d$int_hi) - log(d$int_lo)) / (2 * qnorm(0.975)), p_trial = d$int_p)
}
pcr <- function(key, lab, tag = "HRadj") {
  d <- E10[E10$cohort == key, ]
  data.frame(trial = lab, score = d$score, endpoint = "pCR (OR)", y = d[[paste0(tag, "_logOR")]],
             se = d[[paste0(tag, "_se")]], p_trial = d[[paste0(tag, "_p")]])
}
IMM <- surv("IMmotion150: atezolizumab (+/- bev)", "IMmotion150")
PRI <- surv("PRINCE", "PRINCE")
PEM <- pcr("ISPY2_pembro", "I-SPY2 pembrolizumab"); DUR <- pcr("ISPY2_durva", "I-SPY2 durvalumab")
PEMu <- pcr("ISPY2_pembro", "I-SPY2 pembrolizumab", "overall"); DURu <- pcr("ISPY2_durva", "I-SPY2 durvalumab", "overall")
PEMs <- pcr("ISPY2_pembro", "I-SPY2 pembrolizumab", "strat"); DURs <- pcr("ISPY2_durva", "I-SPY2 durvalumab", "strat")
SETS <- list(`primary: IMmotion150 + 2 I-SPY2` = rbind(IMM, PEM, DUR),
             `+ PRINCE` = rbind(IMM, PEM, DUR, PRI),
             `I-SPY2 unadjusted` = rbind(IMM, PEMu, DURu),
             `I-SPY2 HR-stratified` = rbind(IMM, PEMs, DURs),
             `I-SPY2 only (pCR OR)` = rbind(PEM, DUR))
write.csv(do.call(rbind, lapply(names(SETS), function(s) cbind(set = s, SETS[[s]]))),
          file.path(O, "E10_meta_inputs.csv"), row.names = FALSE)

res <- do.call(rbind, lapply(names(SETS), function(s) do.call(rbind, lapply(ORDER, function(k) {
  d <- SETS[[s]][SETS[[s]]$score == k, ]
  re <- rma(yi = d$y, sei = d$se, method = "REML", test = "knha")
  fe <- rma(yi = d$y, sei = d$se, method = "FE")
  z <- d$y / d$se; zs <- sum(z) / sqrt(nrow(d))
  data.frame(set = s, score = k, k = nrow(d), k_benefit = sum(d$y > 0), k_benefit_p05 = sum(d$y > 0 & d$p_trial < 0.05),
             est = re$b[1], lo = re$ci.lb, hi = re$ci.ub, p = re$pval, tau2 = re$tau2, I2 = re$I2,
             fe_est = fe$b[1], fe_lo = fe$ci.lb, fe_hi = fe$ci.ub, fe_p = fe$pval,
             stouffer_z = zs, stouffer_p = 2 * pnorm(-abs(zs)))
}))))
write.csv(res, file.path(O, "E10_meta.csv"), row.names = FALSE)
options(width = 220)
for (s in names(SETS)) {
  cat("\n===", s, "=== (benefit-oriented log ratio per SD; REML + Hartung-Knapp)\n")
  r <- res[res$set == s, c("score", "k", "k_benefit", "est", "lo", "hi", "p", "I2", "fe_est", "fe_p", "stouffer_z", "stouffer_p")]
  r[, 4:12] <- round(r[, 4:12], 3); print(r, row.names = FALSE)
}
