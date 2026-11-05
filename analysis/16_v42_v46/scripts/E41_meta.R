# E41 analyses (protocol E41_prereg.md): random-effects pooling of the treatment x signature interaction from published
# randomized checkpoint trials. Input: results/E41_adjudicated.csv (one row per assessed comparison; the adjudicated
# table). Writes results/E41_inputs.csv, E41_meta.csv, E41_paired.csv, E41_expectations.csv.
suppressPackageStartupMessages(library(metafor))
R <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
A <- read.csv(file.path(R, "E41_adjudicated.csv"), stringsAsFactors = FALSE, na.strings = c("", "NA"))
Z <- qnorm(0.975)
yes <- function(x) !is.na(x) & tolower(x) == "yes"
logse <- function(lo, hi) (log(hi) - log(lo)) / (2 * Z)

# benefit-oriented log ratio: y > 0 = more checkpoint benefit in the high (or positive) group
pair_y <- function(h, hl, hh, l, ll, lh) {
  data.frame(y = -(log(h) - log(l)), se = sqrt(logse(hl, hh)^2 + logse(ll, lh)^2))
}
E <- A[yes(A$eligible), ]
S <- pair_y(E$hr_high, E$lo_high, E$hi_high, E$hr_low, E$lo_low, E$hi_low)
# a reported interaction (ratio of HRs, high over low) is used directly when the subgroup HRs are not both reported
use_int <- is.na(S$y) & !is.na(E$int_est) & E$int_scale == "ratio_high_over_low"
S$y[use_int] <- -log(E$int_est[use_int]); S$se[use_int] <- logse(E$int_lo[use_int], E$int_hi[use_int])
# per-SD scale: expected standard-normal score difference between the high and low groups at the reported split
# high = top fraction p_hi, low = bottom fraction p_lo (p_lo = 1 - p_hi unless a middle group was set aside, as for
# top vs bottom tertiles); for a complementary split this is phi(z_q) / (q (1 - q)) as in the protocol
p_hi <- ifelse(is.na(E$frac_high), 0.5, E$frac_high)
p_lo <- ifelse(is.na(E$frac_low), 1 - p_hi, E$frac_low)
dq <- dnorm(qnorm(1 - p_hi)) / p_hi + dnorm(qnorm(p_lo)) / p_lo
E$y <- S$y; E$se <- S$se; E$from_int <- use_int; E$split_assumed_median <- is.na(E$frac_high); E$d <- dq
E$y_sd <- E$y / dq; E$se_sd <- E$se / dq
E$survival <- E$endpoint %in% c("OS", "PFS")
E$ipd <- yes(E$in_ipd_six)
E$pdl1_y <- with(E, pair_y(pdl1_hr_pos, pdl1_lo_pos, pdl1_hi_pos, pdl1_hr_neg, pdl1_lo_neg, pdl1_hi_neg))$y
E$pdl1_se <- with(E, pair_y(pdl1_hr_pos, pdl1_lo_pos, pdl1_hi_pos, pdl1_hr_neg, pdl1_lo_neg, pdl1_hi_neg))$se
E$tmb_y <- with(E, pair_y(tmb_hr_high, tmb_lo_high, tmb_hi_high, tmb_hr_low, tmb_lo_low, tmb_hi_low))$y
E$tmb_se <- with(E, pair_y(tmb_hr_high, tmb_lo_high, tmb_hi_high, tmb_hr_low, tmb_lo_low, tmb_hi_low))$se
write.csv(E, file.path(R, "E41_inputs.csv"), row.names = FALSE)

out <- list()
fit <- function(name, d, y = "y", se = "se") {
  d <- d[!is.na(d[[y]]) & !is.na(d[[se]]), ]
  if (nrow(d) < 2) { out[[name]] <<- data.frame(analysis = name, k = nrow(d)); return(invisible(NULL)) }
  m <- rma(yi = d[[y]], sei = d[[se]], method = "REML", test = "knha")
  p <- predict(m)
  out[[name]] <<- data.frame(analysis = name, k = m$k, est = as.numeric(m$b), lo = m$ci.lb, hi = m$ci.ub, p = m$pval,
                             tau2 = m$tau2, I2 = m$I2, pi_lo = p$pi.lb, pi_hi = p$pi.ub,
                             trials = paste(d$trial, collapse = "; "))
  invisible(m)
}
P <- E[E$survival & !E$ipd & !is.na(E$y), ]
m1 <- fit("primary: survival, independent of the six", P)
fit("primary, per SD", P, "y_sd", "se_sd")
fit("all endpoints, independent of the six", E[!E$ipd, ])
fit("all endpoints, per SD", E[!E$ipd, ], "y_sd", "se_sd")
fit("survival, excluding numbers read from figures", P[!grepl("^(yes|estimated)", tolower(P$from_figure)), ])
fit("survival, excluding the figure-position estimate (post hoc)", P[!grepl("^estimated", tolower(P$from_figure)), ])
fit("survival, phase 3 only", P[grepl("3", P$phase), ])
fit("survival, adding the six", E[E$survival, ])
for (tt in names(which(table(P$tumor_group) >= 3))) fit(paste("survival, tumor:", tt), P[P$tumor_group == tt, ])
fit("PD-L1 IHC, survival, same comparisons", P, "pdl1_y", "pdl1_se")
fit("TMB, survival, same comparisons", P, "tmb_y", "tmb_se")
both <- P[!is.na(P$y) & !is.na(P$pdl1_y), ]
fit("signature, survival, comparisons that also report PD-L1", both)
fit("PD-L1, survival, comparisons that also report the signature", both, "pdl1_y", "pdl1_se")
if (!is.null(m1) && m1$k >= 2) {
  l1 <- leave1out(m1)
  out[["leave one out (range)"]] <- data.frame(analysis = "leave one out, survival (lo/hi = range of estimates)",
                                               k = m1$k - 1, lo = min(l1$estimate), hi = max(l1$estimate), p = max(l1$pval))
}
if (!is.null(m1) && m1$k >= 10) {
  eg <- regtest(m1)
  out[["Egger"]] <- data.frame(analysis = "Egger test, survival", k = m1$k, est = eg$zval, p = eg$pval)
}
M <- do.call(rbind, lapply(out, function(x) { for (c in c("est", "lo", "hi", "p", "tau2", "I2", "pi_lo", "pi_hi", "trials"))
  if (!c %in% names(x)) x[[c]] <- NA; x[, c("analysis", "k", "est", "lo", "hi", "p", "tau2", "I2", "pi_lo", "pi_hi", "trials")] }))
write.csv(M, file.path(R, "E41_meta.csv"), row.names = FALSE)
print(M[, 1:10], digits = 3)

both$diff <- both$y - both$pdl1_y
write.csv(both[, c("trial", "endpoint", "signature", "y", "pdl1_def", "pdl1_y", "diff")], file.path(R, "E41_paired.csv"),
          row.names = FALSE)
cat("paired with PD-L1:", nrow(both), "; signature larger in", sum(both$diff > 0), "\n")

g <- function(a) M[M$analysis == a, ]
k_surv <- nrow(P); pr <- g("primary: survival, independent of the six"); ps <- g("primary, per SD")
sg <- g("signature, survival, comparisons that also report PD-L1"); pl <- g("PD-L1, survival, comparisons that also report the signature")
n_int <- sum(yes(E$int_test_reported[!E$ipd])); n_el <- sum(!E$ipd)
X <- data.frame(
  id = c("M1", "M2", "M3", "M4", "M5"),
  expectation = c("at least 10 eligible comparisons with a survival endpoint",
                  "pooled y positive with 95% CI excluding 0",
                  "pooled per-SD estimate between +0.10 and +0.30",
                  "fewer than half of eligible comparisons report a formal interaction test",
                  "signature and PD-L1 pooled y differ by less than 0.20 (comparisons reporting both)"),
  observed = c(sprintf("%d (independent of the six)", k_surv),
               sprintf("%+.3f (%+.3f to %+.3f)", pr$est, pr$lo, pr$hi),
               sprintf("%+.3f (%+.3f to %+.3f)", ps$est, ps$lo, ps$hi),
               sprintf("%d of %d", n_int, n_el),
               if (nrow(sg) && !is.na(sg$est) && nrow(pl) && !is.na(pl$est)) sprintf("%+.3f vs %+.3f (k = %d)", sg$est, pl$est, sg$k) else "not estimable"),
  held = c(k_surv >= 10, isTRUE(pr$lo > 0), isTRUE(ps$est >= 0.10 & ps$est <= 0.30), n_int < n_el / 2,
           isTRUE(abs(sg$est - pl$est) < 0.20)))
write.csv(X, file.path(R, "E41_expectations.csv"), row.names = FALSE)
print(X)

# post hoc (not in the protocol), labelled as such: without the three comparisons whose signature class was in doubt
# (OAK AP/T score, ORIENT-11 and RATIONALE-309 clusters), and without the one covariate-adjusted result (IMbrave150)
ph <- list()
for (nm in c("without doubtful signature classes", "without the covariate-adjusted IMbrave150")) {
  d <- if (grepl("doubtful", nm)) P[!P$trial %in% c("OAK", "ORIENT-11", "RATIONALE-309"), ] else P[P$trial != "IMbrave150", ]
  m <- rma(yi = d$y, sei = d$se, method = "REML", test = "knha")
  ph[[nm]] <- data.frame(analysis = paste("post hoc:", nm), k = m$k, est = as.numeric(m$b), lo = m$ci.lb, hi = m$ci.ub,
                         p = m$pval, I2 = m$I2)
}
PH <- do.call(rbind, ph); write.csv(PH, file.path(R, "E41_posthoc.csv"), row.names = FALSE); print(PH, digits = 3)
