## E6c - analyses requested by peer review of the ceiling claim.
##  (1) paired-bootstrap CI for dC = learner - NetImmune in each held-out cohort, pooled by inverse variance
##  (2) elastic net with NetImmune as an offset (beta from a training-cohort Cox fit), leave-one-cohort-out
##  (3) within-cohort ceiling: 10-fold cross-validated elastic net inside each cohort, 5 repeats
##  (4) rank-based inverse-normal transform of each gene within cohort, then leave-one-cohort-out EN
##  (5) reference: C of stage + age (apparent, from revision/results/B4_cindex.csv)
## Usage: Rscript E6_review.R <best learner key, e.g. rsf | xgb | screenridge>
source("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/E6_common.R")
library(parallel)
best_key <- commandArgs(TRUE)[1]; if (is.na(best_key)) best_key <- "rsf"
BEST <- list(rsf = function(...) fit_rsf(..., threads = 4), xgb = function(...) fit_xgbcox(..., threads = 4),
             screenridge = fit_screenridge, superpc = fit_superpc, pcacox = fit_pcacox, ridge = fit_glmnet(0))[[best_key]]
NB <- 2000
netz <- lapply(CT, function(t) as.numeric(scale(SC[[t]]$NetImmune))); names(netz) <- CT
INT <- lapply(Z, function(z) apply(z, 2, function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))))

## ---- held-out predictions (1, 2, 4) -----------------------------------------------------
fit_offset <- function(te) {
  tr <- setdiff(CT, te); d <- train_test(te); off_tr <- unlist(netz[tr])
  b <- coef(coxph(Surv(d$ttr, d$etr) ~ off_tr))           # NetImmune coefficient, training cohorts only
  set.seed(1)
  cv <- cv.glmnet(d$Xtr, Surv(d$ttr, d$etr), family = "cox", alpha = 0.5, nfolds = 5, standardize = FALSE, offset = b * off_tr)
  r <- as.numeric(predict(cv, d$Xte, s = "lambda.min", newoffset = b * netz[[te]]))
  attr(r, "nsel") <- sum(coef(cv, s = "lambda.min") != 0); r
}
fit_int <- function(te) {
  tr <- setdiff(CT, te)
  fit_glmnet(0.5)(do.call(rbind, INT[tr]), unlist(TT[tr]), unlist(EV[tr]), INT[[te]], seed = 1)
}
jobs <- expand.grid(what = c("en", "best", "offset", "int"), test = CT, stringsAsFactors = FALSE)
pred <- mclapply(seq_len(nrow(jobs)), function(i) {
  j <- jobs[i, ]; d <- train_test(j$test)
  switch(j$what, en = fit_glmnet(0.5)(d$Xtr, d$ttr, d$etr, d$Xte, seed = 1),
         best = BEST(d$Xtr, d$ttr, d$etr, d$Xte, seed = 1), offset = fit_offset(j$test), int = fit_int(j$test))
}, mc.cores = 12, mc.preschedule = FALSE)
names(pred) <- paste(jobs$what, jobs$test)
saveRDS(pred, file.path(O, "E6_review_predictions.rds"))

## paired bootstrap of dC (learner minus NetImmune) over held-out patients
boot_dc <- function(r, te, seed = 7) {
  t <- TT[[te]]; e <- EV[[te]]; n0 <- -netz[[te]]
  obs <- cidx(r, t, e) - cidx(n0, t, e); set.seed(seed)
  b <- replicate(NB, { i <- sample(length(t), replace = TRUE); cidx(r[i], t[i], e[i]) - cidx(n0[i], t[i], e[i]) })
  data.frame(test = te, C = cidx(r, t, e), C_net = cidx(n0, t, e), dC = obs, se = sd(b),
             lo = quantile(b, 0.025), hi = quantile(b, 0.975), p = 2 * min(mean(b <= 0), mean(b >= 0)))
}
dc <- do.call(rbind, lapply(c("en", "best", "offset", "int"), function(w) {
  x <- do.call(rbind, lapply(CT, function(te) boot_dc(pred[[paste(w, te)]], te)))
  wts <- 1 / x$se^2; est <- sum(wts * x$dC) / sum(wts); se <- sqrt(1 / sum(wts))
  x <- rbind(x, data.frame(test = "pooled", C = mean(x$C), C_net = mean(x$C_net), dC = est, se = se,
                           lo = est - 1.96 * se, hi = est + 1.96 * se, p = 2 * pnorm(-abs(est / se))))
  cbind(model = c(en = "Elastic net", best = paste("Best alternative:", best_key), offset = "NetImmune offset + elastic net",
                  int = "Elastic net, inverse-normal genes")[w], x)
}))
rownames(dc) <- NULL
write.csv(dc, file.path(O, "E6_review_dC.csv"), row.names = FALSE)
cat("\n=== held-out C and paired-bootstrap dC vs NetImmune (", NB, "resamples) ===\n"); print(format(dc, digits = 3), row.names = FALSE)
cat("\nnonzero coefficients in offset models:", sapply(CT, function(te) attr(pred[[paste("offset", te)]], "nsel")), "\n")

## ---- (3) within-cohort 10-fold CV elastic net, 5 repeats ----------------------------------
wjobs <- expand.grid(test = CT, rep = 1:5, fold = 1:10, stringsAsFactors = FALSE)
wres <- mclapply(seq_len(nrow(wjobs)), function(i) {
  j <- wjobs[i, ]; X <- Z[[j$test]]; t <- TT[[j$test]]; e <- EV[[j$test]]
  f <- folds_of(nrow(X), 10, seed = 100 + j$rep); a <- f != j$fold
  r <- fit_glmnet(0.5)(X[a, ], t[a], e[a], X[!a, , drop = FALSE], seed = j$rep)
  data.frame(test = j$test, rep = j$rep, idx = which(!a), risk = r)
}, mc.cores = 15, mc.preschedule = FALSE)
w <- do.call(rbind, wres)
within <- do.call(rbind, lapply(CT, function(te) {
  cs <- sapply(1:5, function(k) { s <- w[w$test == te & w$rep == k, ]; s <- s[order(s$idx), ]
                                   cidx(s$risk, TT[[te]][s$idx], EV[[te]][s$idx]) })
  data.frame(test = te, C_cv_mean = mean(cs), C_cv_min = min(cs), C_cv_max = max(cs), C_net = cidx(-netz[[te]], TT[[te]], EV[[te]]))
}))
write.csv(within, file.path(O, "E6_within_cohort_cv.csv"), row.names = FALSE)
cat("\n=== within-cohort 10-fold CV elastic net (5 repeats; C on pooled out-of-fold predictions) ===\n"); print(format(within, digits = 3), row.names = FALSE)

## ---- (5) clinical reference -------------------------------------------------------------
cl <- read.csv(file.path(P, "B4_cindex.csv")); cl <- unique(cl[cl$score == "NetImmune", c("cohort", "n", "C_clin", "C_full")])
cat("\n=== stage + age (apparent C, fitted in the same cohort; B4_cindex) ===\n"); print(cl, row.names = FALSE)
