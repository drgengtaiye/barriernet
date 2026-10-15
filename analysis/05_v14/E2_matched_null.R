## E2 - structure-matched random-axis null for prognosis, with the SAME random gene set
## applied to all three cohorts.
##
## revision/scripts/v2_random_logrank.R drew a fresh gene set inside each cohort's loop, so
## random axis k was a different gene set in each cohort: its per-cohort ranks are valid, but
## its cross-cohort summaries ("significant in all three: 0.7%") combined unrelated axes.
## Here each random axis is one gene set: four random modules of 8/5/6/7 genes minus two of
## 11, from the 14,364 genes measured in all three cohorts, scored by the same row-z mean.
## A random axis is credited with its better direction, chosen once for all cohorts.
suppressMessages({library(survival); library(matrixStats)})
R <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
set.seed(2026)
NRAND <- 1000
SIZES <- c(APM = 8, IFN = 5, TLS_B = 6, T_CELL = 7, CAF_ECM = 11, TAM_SUP = 11)
TAGS <- c("TCGA", "GSE84437", "GSE62254")
zrow <- function(m) { s <- rowSds(m, na.rm = TRUE); s[s == 0 | is.na(s)] <- 1
                      (m - rowMeans(m, na.rm = TRUE)) / s }
D <- lapply(TAGS, function(t) readRDS(file.path(R, paste0("B1_", t, ".rds")))); names(D) <- TAGS
S <- readRDS(file.path(R, "B2_scores.rds"))
pool <- readLines(file.path(O, "bulk_gene_pool.txt"))
prep <- lapply(TAGS, function(t) {
  s <- S[[t]]; z <- zrow(D[[t]]$expr[pool, ])
  idx <- match(s$sample, colnames(z))
  ok <- !is.na(idx) & !is.na(s$time) & s$time > 0 & !is.na(s$event)
  list(z = z[, idx[ok]], time = s$time[ok], event = s$event[ok], net = s$NetImmune[ok])
}); names(prep) <- TAGS

fit <- function(p, x) {
  f <- coxph(Surv(p$time, p$event) ~ scale(x)[, 1])
  lr <- survdiff(Surv(p$time, p$event) ~ (x > median(x)))
  c(b = unname(coef(f)), se = sqrt(vcov(f)[1, 1]), p_lr = pchisq(lr$chisq, 1, lower.tail = FALSE))
}
br <- cumsum(SIZES); st <- c(1, head(br, -1) + 1)
axis_of <- function(z, g) {
  mods <- sapply(seq_along(br), function(j) colMeans(z[g[st[j]:br[j]], , drop = FALSE]))
  rowMeans(mods[, 1:4]) - rowMeans(mods[, 5:6])
}
obs <- sapply(TAGS, function(t) fit(prep[[t]], prep[[t]]$net))
rnd <- array(NA, c(NRAND, 3, 3), dimnames = list(NULL, TAGS, c("b", "se", "p_lr")))
for (k in seq_len(NRAND)) {
  g <- match(sample(pool, sum(SIZES)), pool)
  for (t in TAGS) rnd[k, t, ] <- fit(prep[[t]], axis_of(prep[[t]]$z, g))
}
fe <- function(b, se) { w <- 1 / se^2; sum(w * b) / sum(w) }
obs_pool <- fe(obs["b", ], obs["se", ])
rnd_pool <- apply(rnd, 1, function(m) fe(m[, "b"], m[, "se"]))
per <- do.call(rbind, lapply(TAGS, function(t) data.frame(
  cohort = t, HR_obs = exp(obs["b", t]),
  rank = sum(abs(rnd[, t, "b"]) >= abs(obs["b", t])) + 1, n_axes = NRAND + 1,
  p_emp = (sum(abs(rnd[, t, "b"]) >= abs(obs["b", t])) + 1) / (NRAND + 1),
  frac_random_logrank_sig = mean(rnd[, t, "p_lr"] < 0.05))))
same <- apply(sign(rnd[, , "b"]), 1, function(s) length(unique(s)) == 1)
sig3 <- apply(rnd[, , "p_lr"] < 0.05, 1, all)
joint <- data.frame(
  obs_pooled_HR = exp(obs_pool),
  rank_pooled = sum(abs(rnd_pool) >= abs(obs_pool)) + 1,
  p_emp_pooled = (sum(abs(rnd_pool) >= abs(obs_pool)) + 1) / (NRAND + 1),
  random_pooled_absHR_p95 = exp(quantile(abs(rnd_pool), 0.95)),
  random_same_direction_all3 = mean(same),
  random_sig_all3_same_direction = mean(sig3 & same))
write.csv(per, file.path(O, "E2_null_per_cohort.csv"), row.names = FALSE)
write.csv(joint, file.path(O, "E2_null_joint.csv"), row.names = FALSE)
saveRDS(list(obs = obs, rnd = rnd, rnd_pool = rnd_pool), file.path(O, "E2_null_draws.rds"))
print(format(per, digits = 3)); print(format(joint, digits = 3))
