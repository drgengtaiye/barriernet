## E21 step 4 - cross-trial synthesis of the score x checkpoint-arm interaction (pre-registered:
## E17_E22_prereg.md, E21 and Deviations 2026-10-05). Extends E10_meta.R from 3 to 6 primary
## randomised comparisons; inputs and harmonisation are E10_meta.R's, unchanged:
##   IMmotion150         atezolizumab +/- bevacizumab vs sunitinib - PFS, interaction HR (E5_trials.csv)
##   I-SPY2              pembrolizumab + paclitaxel vs paclitaxel   - pCR, interaction OR, HR-adjusted (E10)
##   I-SPY2              durvalumab + olaparib + paclitaxel vs pac. - pCR, interaction OR, HR-adjusted (E10)
##   CheckMate 025       nivolumab vs everolimus                    - PFS, interaction HR (E21_trials.py)
##   JAVELIN Renal 101   avelumab + axitinib vs sunitinib           - PFS, interaction HR (E21_trials.py)
##   NeoTRIP             chemo + atezolizumab vs chemo              - pCR, interaction OR, batch-adjusted (E21)
## Sensitivity only: PRINCE (E5_trials.csv) and DUTRENEO (E21_trials.py).
## y = -log(interaction HR) or +log(interaction OR) per SD (benefit-oriented: y > 0 = the score predicts
## more benefit from the checkpoint-containing arm); HR and OR are different scales, so pooled y is a
## direction-and-size summary, reported with the unweighted Stouffer z of y / SE and the count with y > 0.
## Pooling (manuscript rule): random effects, point estimate REML (as E10_meta.R); reported interval =
## the wider of the REML Hartung-Knapp interval and the DerSimonian-Laird (Wald) interval, reported P =
## the larger of the two. This reproduces every published E10 value (checked below). Fixed-effect and
## DL-estimate (DL + DL-HK) versions are written alongside.
suppressMessages(library(metafor))
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
E5 <- read.csv(file.path(O, "E5_trials.csv"))
E10 <- read.csv(file.path(O, "E10_ispy2_interactions.csv"))
E21 <- read.csv(file.path(O, "E21_interactions.csv"))
ORDER <- c("CAF_ECM", "TAM_SUP", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1", "IFN", "APM", "TLS_B", "T_CELL",
           "Barrier", "Immune", "NetImmune")

## --- inputs (E10_meta.R functions unchanged, plus one for E21)
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
e21 <- function(key, ep, model, lab) {
  d <- E21[E21$key == key & E21$endpoint == ep & E21$model == model, ]
  stopifnot(nrow(d) == length(ORDER))
  data.frame(trial = lab, score = d$score, endpoint = paste0(ep, if (model == "Cox") " (HR)" else " (OR)"),
             y = d$y, se = d$int_se, p_trial = d$int_p)
}
IMM <- surv("IMmotion150: atezolizumab (+/- bev)", "IMmotion150")
PRI <- surv("PRINCE", "PRINCE")
PEM <- pcr("ISPY2_pembro", "I-SPY2 pembrolizumab"); DUR <- pcr("ISPY2_durva", "I-SPY2 durvalumab-olaparib")
CM <- e21("CM025", "PFS", "Cox", "CheckMate 025")
CMos <- e21("CM025", "OS", "Cox", "CheckMate 025 (OS)")
JR <- e21("JR101", "PFS", "Cox", "JAVELIN Renal 101")
NT <- e21("NeoTRIP", "pCR", "logistic + batch", "NeoTRIP")
NTu <- e21("NeoTRIP", "pCR", "logistic", "NeoTRIP (no batch)")
DT <- e21("DUTRENEO", "pCR", "logistic", "DUTRENEO")
DT1 <- e21("DUTRENEO_1M", "pCR", "logistic", "DUTRENEO (>= 1M)")
PRIM <- list(IMM, PEM, DUR, CM, JR, NT); names(PRIM) <- sapply(PRIM, function(d) d$trial[1])

SETS <- list(`reproduction: 3 published primary` = rbind(IMM, PEM, DUR),
             `primary: 6 comparisons` = do.call(rbind, PRIM))
for (nm in names(PRIM)) SETS[[paste("primary without", nm)]] <- do.call(rbind, PRIM[names(PRIM) != nm])
SETS <- c(SETS, list(
  `+ PRINCE + DUTRENEO (8)` = do.call(rbind, c(PRIM, list(PRI, DT))),
  `+ PRINCE (7)` = do.call(rbind, c(PRIM, list(PRI))),
  `+ DUTRENEO (7)` = do.call(rbind, c(PRIM, list(DT))),
  `renal only (IMmotion150, CheckMate 025, JAVELIN)` = rbind(IMM, CM, JR),
  `breast only (I-SPY2 x 2, NeoTRIP)` = rbind(PEM, DUR, NT),
  `3 new primary only` = rbind(CM, JR, NT),
  `primary, CheckMate 025 OS for PFS` = do.call(rbind, c(PRIM[names(PRIM) != "CheckMate 025"], list(CMos))),
  `primary, NeoTRIP without batch` = do.call(rbind, c(PRIM[names(PRIM) != "NeoTRIP"], list(NTu))),
  `+ PRINCE + DUTRENEO (>= 1M libraries)` = do.call(rbind, c(PRIM, list(PRI, DT1)))))
write.csv(do.call(rbind, lapply(names(SETS), function(s) cbind(set = s, SETS[[s]]))),
          file.path(O, "E21_meta_inputs.csv"), row.names = FALSE)

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
res <- do.call(rbind, lapply(names(SETS), function(s) do.call(rbind, lapply(ORDER, function(k)
  cbind(set = s, score = k, pool(SETS[[s]][SETS[[s]]$score == k, ]))))))
write.csv(res, file.path(O, "E21_meta.csv"), row.names = FALSE)
w <- res[res$reml_warning != "", c("set", "score", "tau2_reml", "tau2_dl", "est", "dl_est", "reml_warning")]
cat("REML warnings (metafor sets tau^2 = 0 where Fisher scoring stalls):\n"); print(w, row.names = FALSE)

## --- reproduction check against the published 3-comparison values (manuscript Table / E10_meta.csv)
pub <- data.frame(score = c("NetImmune", "CAF_ECM"), est = c(0.47, -0.37), lo = c(-0.04, -1.16), hi = c(0.97, 0.42), p = c(0.057, 0.18))
r3 <- res[res$set == "reproduction: 3 published primary", ]
e10 <- read.csv(file.path(O, "E10_meta.csv")); e10 <- e10[startsWith(e10$set, "primary"), ]
fe_w <- (e10$fe_hi - e10$fe_lo) > (e10$hi - e10$lo)                 # E10_fig_trials.py cons() rule
e10$c_lo <- ifelse(fe_w, e10$fe_lo, e10$lo); e10$c_hi <- ifelse(fe_w, e10$fe_hi, e10$hi); e10$c_p <- pmax(e10$p, e10$fe_p)
rep <- merge(r3[, c("score", "est", "lo", "hi", "p", "k_benefit", "stouffer_z")],
             setNames(e10[, c("score", "est", "c_lo", "c_hi", "c_p")], c("score", "E10_est", "E10_lo", "E10_hi", "E10_p")), by = "score")
rep$max_abs_diff_vs_E10 <- apply(abs(rep[, c("est", "lo", "hi", "p")] - rep[, c("E10_est", "E10_lo", "E10_hi", "E10_p")]), 1, max)
rep <- merge(rep, setNames(pub, c("score", "pub_est", "pub_lo", "pub_hi", "pub_p")), by = "score", all.x = TRUE)
rep$matches_published_2dp <- with(rep, ifelse(is.na(pub_est), NA, round(est, 2) == pub_est & round(lo, 2) == pub_lo &
                                                round(hi, 2) == pub_hi & abs(p - pub_p) < 0.0051))
rep <- rep[match(ORDER, rep$score), ]
write.csv(rep, file.path(O, "E21_reproduction.csv"), row.names = FALSE)
options(width = 230)
cat("=== reproduction of the published 3-comparison pooled values (conservative rule) ===\n")
print(rep[, c("score", "est", "lo", "hi", "p", "k_benefit", "stouffer_z", "max_abs_diff_vs_E10",
                                        "pub_est", "pub_lo", "pub_hi", "pub_p", "matches_published_2dp")], row.names = FALSE, digits = 3)

## --- expectations (pre-registered, pooled over the 6 primary comparisons)
EXP <- c(NetImmune = ">0", Immune = ">0", IFN = ">0", CAF_ECM = "<=0")
p6 <- res[res$set == "primary: 6 comparisons" & res$score %in% names(EXP), ]
p6$expected <- EXP[p6$score]
p6$held <- ifelse(p6$expected == ">0", p6$est > 0, p6$est <= 0)
p6$supported_conservative <- p6$held & p6$p < 0.05
write.csv(p6[, c("score", "expected", "k", "k_benefit", "est", "lo", "hi", "p", "interval", "stouffer_z", "stouffer_p",
                 "held", "supported_conservative")], file.path(O, "E21_expectations_pooled.csv"), row.names = FALSE)

show <- c("score", "k", "k_benefit", "k_benefit_p05", "est", "lo", "hi", "p", "interval", "I2", "fe_est", "fe_p",
          "stouffer_z", "stouffer_p")
for (s in names(SETS)) {
  cat("\n===", s, "=== (benefit-oriented log ratio per SD; REML, wider of HK / DL interval, larger P)\n")
  r <- res[res$set == s, show]; r[, c(5:7, 8, 10:14)] <- lapply(r[, c(5:7, 8, 10:14)], function(v) round(v, 3)); print(r, row.names = FALSE)
}
cat("\n=== leave-one-trial-out, NetImmune / CAF_ECM / Immune / IFN ===\n")
loo <- res[grepl("^primary", res$set) & res$score %in% c("NetImmune", "CAF_ECM", "Immune", "IFN"),
           c("set", "score", "k_benefit", "est", "lo", "hi", "p", "stouffer_z")]
loo[, 4:8] <- round(loo[, 4:8], 3); print(loo[order(loo$score, loo$set), ], row.names = FALSE)
