## E6 - shared inputs and learners for the ceiling-robustness analyses.
## Inputs are exactly those of revision2/scripts/S11_ceiling.R: B1 expression matrices, genes
## common and finite in TCGA-STAD, GSE84437 and GSE62254 (14,386), z-scored within cohort.
## Every learner is tuned only inside the training cohorts and returns a risk score for the
## held-out cohort (higher = higher risk).
suppressMessages({library(survival); library(glmnet); library(matrixStats); library(ranger); library(xgboost)})
P <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
CT <- c("TCGA", "GSE84437", "GSE62254")
zrow <- function(m) { s <- rowSds(m, na.rm = TRUE); s[s == 0 | is.na(s)] <- 1; (m - rowMeans(m, na.rm = TRUE)) / s }
SC <- readRDS(file.path(P, "B2_scores.rds"))
E <- lapply(CT, function(t) readRDS(file.path(P, paste0("B1_", t, ".rds")))$expr); names(E) <- CT
g <- Reduce(intersect, lapply(E, rownames))
Z <- lapply(E, function(m) t(zrow(m[g, , drop = FALSE])))
keep <- Reduce(`&`, lapply(Z, function(z) apply(is.finite(z), 2, all)))
Z <- lapply(Z, function(z) z[, keep, drop = FALSE]); rm(E)
TT <- lapply(CT, function(t) SC[[t]]$time); names(TT) <- CT
EV <- lapply(CT, function(t) SC[[t]]$event); names(EV) <- CT
stopifnot(ncol(Z[[1]]) == 14386)

cidx <- function(risk, time, event) {
  np <- -as.numeric(risk)
  as.numeric(concordance(Surv(time, event) ~ np)$concordance)
}

## Cox score statistic at beta = 0 for every column of X at once (Breslow ties).
score_z <- function(X, time, event) {
  o <- order(-time); X <- X[o, , drop = FALSE]; t <- time[o]; d <- event[o]
  n_at <- seq_along(t)
  S1 <- colCumsums(X); S2 <- colCumsums(X^2)
  last <- !duplicated(t, fromLast = TRUE)   # risk set at a time includes all ties
  idx <- sapply(seq_along(t), function(i) max(which(t == t[i])))
  m1 <- S1[idx, , drop = FALSE] / n_at[idx]; m2 <- S2[idx, , drop = FALSE] / n_at[idx]
  ev <- d == 1
  U <- colSums(X[ev, , drop = FALSE] - m1[ev, , drop = FALSE])
  V <- colSums(m2[ev, , drop = FALSE] - m1[ev, , drop = FALSE]^2)
  U / sqrt(pmax(V, 1e-12))
}

folds_of <- function(n, k = 5, seed = 1) { set.seed(seed); sample(rep(seq_len(k), length.out = n)) }

## ---- learners: f(Xtr, ttr, etr, Xte, seed) -> risk for Xte -------------------------------
fit_glmnet <- function(alpha) function(Xtr, ttr, etr, Xte, seed = 1) {
  set.seed(seed)
  cv <- cv.glmnet(Xtr, Surv(ttr, etr), family = "cox", alpha = alpha, nfolds = 5, standardize = FALSE)
  as.numeric(predict(cv, Xte, s = "lambda.min"))
}

## (i) Cox on the top principal components of the training data; k by inner CV
fit_pcacox <- function(Xtr, ttr, etr, Xte, seed = 1, ks = c(1, 2, 3, 5, 10, 20, 40)) {
  mu <- colMeans(Xtr); sv <- svd(sweep(Xtr, 2, mu), nu = 0, nv = max(ks))
  Ftr <- sweep(Xtr, 2, mu) %*% sv$v; Fte <- sweep(Xte, 2, mu) %*% sv$v
  f <- folds_of(nrow(Xtr), 5, seed)
  cvC <- sapply(ks, function(k) mean(sapply(1:5, function(j) {
    a <- f != j; m <- coxph(Surv(ttr[a], etr[a]) ~ Ftr[a, 1:k, drop = FALSE])
    cidx(Ftr[!a, 1:k, drop = FALSE] %*% coef(m), ttr[!a], etr[!a]) })))
  k <- ks[which.max(cvC)]
  m <- coxph(Surv(ttr, etr) ~ Ftr[, 1:k, drop = FALSE])
  r <- as.numeric(Fte[, 1:k, drop = FALSE] %*% coef(m)); attr(r, "tuned") <- k; r
}

## (ii) univariate-screened ridge Cox; number of genes and lambda by inner CV (screening in-fold)
fit_screenridge <- function(Xtr, ttr, etr, Xte, seed = 1, ms = c(50, 100, 200, 500, 1000, 2000)) {
  f <- folds_of(nrow(Xtr), 5, seed); lam <- NULL
  res <- lapply(1:5, function(j) {
    a <- f != j; sz <- abs(score_z(Xtr[a, ], ttr[a], etr[a])); ord <- order(-sz)
    lapply(ms, function(m) {
      gi <- ord[1:m]
      fit <- glmnet(Xtr[a, gi], Surv(ttr[a], etr[a]), family = "cox", alpha = 0, standardize = FALSE,
                    lambda = exp(seq(log(50), log(0.05), length.out = 30)))
      pr <- predict(fit, Xtr[!a, gi])
      apply(pr, 2, cidx, time = ttr[!a], event = etr[!a])
    })
  })
  grid <- sapply(seq_along(ms), function(i) rowMeans(sapply(res, `[[`, i)))   # 30 lambdas x ms
  best <- which(grid == max(grid), arr.ind = TRUE)[1, ]
  lamv <- exp(seq(log(50), log(0.05), length.out = 30)); m <- ms[best[2]]; l <- lamv[best[1]]
  gi <- order(-abs(score_z(Xtr, ttr, etr)))[1:m]
  fit <- glmnet(Xtr[, gi], Surv(ttr, etr), family = "cox", alpha = 0, standardize = FALSE, lambda = lamv)
  r <- as.numeric(predict(fit, Xte[, gi], s = l)); attr(r, "tuned") <- sprintf("m=%d, lambda=%.3g", m, l); r
}

## (ii-b) supervised principal components (Bair & Tibshirani): first PC of screened genes
fit_superpc <- function(Xtr, ttr, etr, Xte, seed = 1, ms = c(20, 50, 100, 200, 500, 1000)) {
  f <- folds_of(nrow(Xtr), 5, seed)
  pc1 <- function(Xa, Xb, gi) { mu <- colMeans(Xa[, gi]); v <- svd(sweep(Xa[, gi], 2, mu), nu = 0, nv = 1)$v
    list(a = as.numeric(sweep(Xa[, gi], 2, mu) %*% v), b = as.numeric(sweep(Xb[, gi], 2, mu) %*% v)) }
  cvC <- sapply(ms, function(m) mean(sapply(1:5, function(j) {
    a <- f != j; gi <- order(-abs(score_z(Xtr[a, ], ttr[a], etr[a])))[1:m]
    p <- pc1(Xtr[a, ], Xtr[!a, ], gi); b <- coef(coxph(Surv(ttr[a], etr[a]) ~ p$a))
    cidx(p$b * b, ttr[!a], etr[!a]) })))
  m <- ms[which.max(cvC)]; gi <- order(-abs(score_z(Xtr, ttr, etr)))[1:m]
  p <- pc1(Xtr, Xte, gi); b <- coef(coxph(Surv(ttr, etr) ~ p$a))
  r <- p$b * b; attr(r, "tuned") <- m; r
}

## (iii-a) random survival forest (ranger) on in-fold univariate-screened genes;
##         min.node.size and screen size tuned on out-of-bag C within the training cohorts
fit_rsf <- function(Xtr, ttr, etr, Xte, seed = 1, ms = c(200, 1000), nodes = c(10, 30, 60), threads = 4) {
  sz <- abs(score_z(Xtr, ttr, etr)); best <- NULL
  for (m in ms) for (nd in nodes) {
    gi <- order(-sz)[1:m]
    fit <- ranger(x = Xtr[, gi], y = Surv(ttr, etr), num.trees = 500, min.node.size = nd,
                  mtry = floor(sqrt(m)), splitrule = "logrank", num.threads = threads, seed = seed)
    if (is.null(best) || fit$prediction.error < best$err) best <- list(err = fit$prediction.error, fit = fit, gi = gi, m = m, nd = nd)
  }
  ## note: out-of-bag C is optimistic for the screened genes (screening used all training
  ## patients) but that only affects tuning, not the held-out evaluation
  chf <- predict(best$fit, Xte[, best$gi], num.threads = threads)$chf
  r <- rowSums(chf); attr(r, "tuned") <- sprintf("m=%d, node=%d", best$m, best$nd); r
}

## (iii-b) gradient-boosted Cox (xgboost, survival:cox) on in-fold screened genes; rounds by inner CV
fit_xgbcox <- function(Xtr, ttr, etr, Xte, seed = 1, m = 1000, threads = 4) {
  gi <- order(-abs(score_z(Xtr, ttr, etr)))[1:m]
  lab <- ifelse(etr == 1, ttr, -ttr)
  d <- xgb.DMatrix(Xtr[, gi], label = lab)
  best <- NULL
  for (depth in c(1, 2)) {
    prm <- list(objective = "survival:cox", eval_metric = "cox-nloglik", eta = 0.02, max_depth = depth,
                subsample = 0.7, colsample_bytree = 0.3, min_child_weight = 5, nthread = threads)
    set.seed(seed)
    cv <- xgb.cv(params = prm, data = d, nrounds = 1500, nfold = 5, early_stopping_rounds = 100, verbose = 0)
    sc <- min(cv$evaluation_log$test_cox_nloglik_mean)
    nr <- which.min(cv$evaluation_log$test_cox_nloglik_mean)
    if (is.null(best) || sc < best$sc) best <- list(sc = sc, prm = prm, nr = nr, depth = depth)
  }
  set.seed(seed)
  fit <- xgb.train(params = best$prm, data = d, nrounds = best$nr, verbose = 0)
  r <- predict(fit, xgb.DMatrix(Xte[, gi])); attr(r, "tuned") <- sprintf("depth=%d, rounds=%d", best$depth, best$nr); r
}

train_test <- function(te) {
  tr <- setdiff(CT, te)
  list(Xtr = do.call(rbind, Z[tr]), ttr = unlist(TT[tr]), etr = unlist(EV[tr]),
       Xte = Z[[te]], tte = TT[[te]], ete = EV[[te]])
}
