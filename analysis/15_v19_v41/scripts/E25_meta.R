## E25 - pooling for Q1 (within-arm effects by arm) and Q2 (mutually adjusted Immune and CAF_ECM
## interactions), pre-registered in E25_prereg.md. Inputs from E25_trials.py. pool() is E21_meta.R's,
## unchanged (REML estimate; wider of Hartung-Knapp and DerSimonian-Laird interval; larger P).
suppressMessages(library(metafor))
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
W <- read.csv(file.path(O, "E25_within_arm_inputs.csv"))
J <- read.csv(file.path(O, "E25_joint_inputs.csv"))

pool <- function(d) {
  note <- ""
  re <- withCallingHandlers(rma(yi = d$y, sei = d$se, method = "REML", test = "knha"),
                            warning = function(w) { note <<- gsub("\n", " ", conditionMessage(w)); invokeRestart("muffleWarning") })
  dl <- rma(yi = d$y, sei = d$se, method = "DL")
  dlh <- rma(yi = d$y, sei = d$se, method = "DL", test = "knha")
  fe <- rma(yi = d$y, sei = d$se, method = "FE")
  hk_wider <- (re$ci.ub - re$ci.lb) >= (dl$ci.ub - dl$ci.lb)
  z <- d$y / d$se; zs <- sum(z) / sqrt(nrow(d))
  data.frame(k = nrow(d), k_benefit = sum(d$y > 0), k_benefit_p05 = sum(d$y > 0 & d$p_trial < 0.05),
             k_against_p05 = sum(d$y < 0 & d$p_trial < 0.05),
             est = re$b[1], lo = if (hk_wider) re$ci.lb else dl$ci.lb, hi = if (hk_wider) re$ci.ub else dl$ci.ub,
             p = max(re$pval, dl$pval), interval = if (hk_wider) "HK" else "DL",
             hk_lo = re$ci.lb, hk_hi = re$ci.ub, hk_p = re$pval, tau2_reml = re$tau2, I2 = re$I2,
             dl_est = dl$b[1], dl_lo = dl$ci.lb, dl_hi = dl$ci.ub, dl_p = dl$pval, tau2_dl = dl$tau2,
             dlhk_lo = dlh$ci.lb, dlhk_hi = dlh$ci.ub, dlhk_p = dlh$pval,
             fe_est = fe$b[1], fe_lo = fe$ci.lb, fe_hi = fe$ci.ub, fe_p = fe$pval,
             stouffer_z = zs, stouffer_p = 2 * pnorm(-abs(zs)), reml_warning = note)
}
SETS <- list(`all six` = c("breast", "renal"), breast = "breast", renal = "renal")

## Q1: within-arm effect, pooled separately over checkpoint and control arms
W$p_trial <- W$p
q1 <- do.call(rbind, lapply(names(SETS), function(s) do.call(rbind, lapply(unique(W$score), function(k)
  do.call(rbind, lapply(c("checkpoint", "control"), function(a) {
    d <- W[W$score == k & W$arm == a & W$cancer %in% SETS[[s]], ]
    cbind(question = "Q1 within-arm", set = s, score = k, arm = a, pool(d)) }))))))

## Q2: mutually adjusted interactions (and the unadjusted ones from the same run, for comparison)
mk <- function(y, se) data.frame(y = J[[y]], se = J[[se]], p_trial = 2 * pnorm(-abs(J[[y]] / J[[se]])), cancer = J$cancer)
Q2IN <- list(`Immune, adjusted for CAF_ECM` = mk("immune_adj_y", "immune_adj_se"),
             `CAF_ECM, adjusted for Immune` = mk("caf_adj_y", "caf_adj_se"),
             `Immune, unadjusted` = mk("immune_unadj_y", "immune_unadj_se"),
             `CAF_ECM, unadjusted` = mk("caf_unadj_y", "caf_unadj_se"))
q2 <- do.call(rbind, lapply(names(SETS), function(s) do.call(rbind, lapply(names(Q2IN), function(k) {
  d <- Q2IN[[k]]; d <- d[d$cancer %in% SETS[[s]], ]
  cbind(question = "Q2 interaction", set = s, score = k, arm = "interaction", pool(d)) }))))

res <- rbind(q1, q2)
write.csv(res, file.path(O, "E25_meta.csv"), row.names = FALSE)

show <- c("set", "score", "arm", "k", "k_benefit", "est", "lo", "hi", "p", "interval", "I2", "stouffer_p")
fmt <- function(r) { r[, c("est", "lo", "hi", "p", "I2", "stouffer_p")] <- lapply(r[, c("est", "lo", "hi", "p", "I2", "stouffer_p")], function(v) signif(v, 3)); r }
cat("\n=== Q1: within-arm effect per SD (benefit-oriented: + = better outcome with a higher score) ===\n")
print(fmt(q1[q1$score %in% c("CAF_ECM", "Immune", "NetImmune", "Barrier", "TAM_SUP", "ESTIMATE", "CD8A"), show]), row.names = FALSE)
cat("\n=== Q2: Immune and CAF_ECM interactions, mutually adjusted and unadjusted ===\n")
print(fmt(q2[, show]), row.names = FALSE)
cat("\nImmune-CAF_ECM correlation per comparison:", paste(sprintf("%s %.2f", J$trial, J$r_immune_caf), collapse = "; "), "\n")
