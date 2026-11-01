## E23 (post hoc, not pre-specified): are the gastro-oesophageal checkpoint cohorts consistent with the
## pan-cancer response estimates, and how much power did each have?
##   - 21 non-gastric checkpoint cohorts (E12 "all") + 3 gastro-oesophageal cohorts (PRJEB25780, PRJEB60680
##     primary processing, PERFECT); estimates taken unchanged from E12 / E22.
##   - subgroup test: DL random-effects meta-regression on a gastro-oesophageal indicator (z and Knapp-Hartung).
##   - consistency: each gastro-oesophageal estimate against the 21-cohort pooled estimate, with the variance of a
##     new cohort's estimate (its SE^2 + SE(pooled)^2 + tau^2, DL), i.e. a prediction-based test.
##   - power: two-sided alpha 0.05 to detect the 21-cohort pooled effect at each cohort's standard error.
## Pooling rule as elsewhere: DL estimate, the more conservative of the DL and Hartung-Knapp intervals.
suppressMessages(library(metafor))
R <- "results"
SC <- c("Immune", "IFN", "APM", "T_CELL", "TLS_B", "NetImmune", "CAF_ECM", "TAM_SUP", "CXCL9_SPP1")
e12 <- read.csv(file.path(R, "E12_icb_per_cohort.csv"))
pan <- subset(e12, set %in% c("E5", "new"), c(cohort, tumour, n, n_resp, score, logOR, se))
e22 <- read.csv(file.path(R, "E22_meta_inputs.csv"))
gej <- subset(e22, cohort %in% c("PRJEB25780", "PRJEB60680", "PERFECT_GSE165252"))
gej <- data.frame(cohort = gej$cohort, tumour = "gastro-oesophageal", n = gej$n, n_resp = gej$n_resp,
                  score = gej$score, logOR = gej$y, se = gej$se)
stopifnot(length(unique(pan$cohort)) == 21, length(unique(gej$cohort)) == 3)
d <- rbind(cbind(pan, gej_flag = 0), cbind(gej, gej_flag = 1))

cons <- function(m, h) {               # conservative of DL and HK intervals
  wide_h <- (h$ci.ub - h$ci.lb) > (m$ci.ub - m$ci.lb)
  c(est = m$b[1], lo = if (wide_h) h$ci.lb else m$ci.lb, hi = if (wide_h) h$ci.ub else m$ci.ub,
    p = max(m$pval, h$pval))
}
pool <- function(x) {
  m <- rma(yi = x$logOR, sei = x$se, method = "DL"); h <- rma(yi = x$logOR, sei = x$se, method = "DL", test = "knha")
  c(k = nrow(x), n = sum(x$n), n_resp = sum(x$n_resp), cons(m, h), I2 = m$I2, k_pos = sum(x$logOR > 0))
}
out <- do.call(rbind, lapply(SC, function(s) {
  x <- d[d$score == s, ]
  a <- pool(x); p21 <- pool(x[x$gej_flag == 0, ]); g3 <- pool(x[x$gej_flag == 1, ])
  mr <- rma(yi = logOR, sei = se, mods = ~ gej_flag, data = x, method = "DL")
  mrh <- rma(yi = logOR, sei = se, mods = ~ gej_flag, data = x, method = "DL", test = "knha")
  data.frame(score = s,
             all24_k = a["k"], all24_n = a["n"], all24_resp = a["n_resp"], all24_est = a["est"], all24_lo = a["lo"],
             all24_hi = a["hi"], all24_p = a["p"], all24_kpos = a["k_pos"], all24_I2 = a["I2"],
             pan21_est = p21["est"], pan21_lo = p21["lo"], pan21_hi = p21["hi"],
             gej3_est = g3["est"], gej3_lo = g3["lo"], gej3_hi = g3["hi"], gej3_p = g3["p"], gej3_kpos = g3["k_pos"],
             diff_gej_minus_pan = mr$b[2], diff_lo = mr$ci.lb[2], diff_hi = mr$ci.ub[2],
             p_subgroup = max(mr$QMp, mrh$QMp), row.names = NULL)
}))
write.csv(out, file.path(R, "E23_pooled_with_gastric.csv"), row.names = FALSE)

pw <- do.call(rbind, lapply(SC, function(s) {
  th <- out$pan21_est[out$score == s]
  p21 <- d[d$score == s & d$gej_flag == 0, ]
  m21 <- rma(yi = p21$logOR, sei = p21$se, method = "DL"); se_th <- m21$se; tau2 <- m21$tau2
  g <- d[d$score == s & d$gej_flag == 1, ]
  data.frame(score = s, cohort = g$cohort, n = g$n, n_resp = g$n_resp, logOR = g$logOR, se = g$se,
             pan21 = th,
             p_differs_from_pan21 = 2 * pnorm(-abs((g$logOR - th) / sqrt(g$se^2 + se_th^2 + tau2))),
             power_at_pan21 = pnorm(abs(th) / g$se - qnorm(0.975)) + pnorm(-abs(th) / g$se - qnorm(0.975)),
             n_for_80pct = ceiling(g$n * (g$se * (qnorm(0.975) + qnorm(0.8)) / th)^2))
}))
write.csv(pw, file.path(R, "E23_power_consistency.csv"), row.names = FALSE)

f <- function(x) sprintf("%+.2f", x)
sink(file.path(R, "E23_log.txt"))
cat("E23 post hoc: gastro-oesophageal consistency and power\n\n")
for (i in seq_len(nrow(out))) with(out[i, ], cat(sprintf(
  "%-10s 24 cohorts %s (%s to %s) P=%.2g, %d/%d positive | 21 pan %s | GEJ3 %s (%s to %s) %d/3 positive | GEJ-pan %s (%s to %s) P=%.2g\n",
  score, f(all24_est), f(all24_lo), f(all24_hi), all24_p, all24_kpos, all24_k, f(pan21_est), f(gej3_est), f(gej3_lo),
  f(gej3_hi), gej3_kpos, f(diff_gej_minus_pan), f(diff_lo), f(diff_hi), p_subgroup)))
cat("\n"); print(format(pw, digits = 3), row.names = FALSE)
sink()
cat(readLines(file.path(R, "E23_log.txt")), sep = "\n")
