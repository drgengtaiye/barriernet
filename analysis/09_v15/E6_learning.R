## E6b - learning curve: is out-of-cohort C still rising at the full training size?
## Usage: Rscript E6_learning.R <learner>   (learner = "en" or a key of LEARNERS below)
## For each held-out cohort, the two training cohorts are pooled and subsampled at n = 100,
## 200 and 400 (20 random subsamples each) and used whole (5 repeats with different tuning
## seeds); every fit is tuned inside its own subsample and evaluated on the held-out cohort.
source("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/E6_common.R")
library(parallel)
arg <- commandArgs(TRUE)[1]; if (is.na(arg)) arg <- "en"
LEARNERS <- list(en = fit_glmnet(0.5), ridge = fit_glmnet(0), pcacox = fit_pcacox, screenridge = fit_screenridge,
                 superpc = fit_superpc, rsf = function(...) fit_rsf(..., threads = 2), xgb = function(...) fit_xgbcox(..., threads = 2))
f <- LEARNERS[[arg]]
SIZES <- c(100, 200, 400); NREP <- 20; NFULL <- 5
jobs <- rbind(expand.grid(test = CT, n = SIZES, rep = seq_len(NREP), stringsAsFactors = FALSE),
              expand.grid(test = CT, n = NA, rep = seq_len(NFULL), stringsAsFactors = FALSE))
res <- mclapply(seq_len(nrow(jobs)), function(i) {
  j <- jobs[i, ]; d <- train_test(j$test)
  set.seed(1000 * j$rep + ifelse(is.na(j$n), 0, j$n))
  idx <- if (is.na(j$n)) seq_len(nrow(d$Xtr)) else sample(nrow(d$Xtr), j$n)
  r <- tryCatch(f(d$Xtr[idx, ], d$ttr[idx], d$etr[idx], d$Xte, seed = j$rep), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  data.frame(learner = arg, test = j$test, n = length(idx), full = is.na(j$n), rep = j$rep, C = cidx(r, d$tte, d$ete))
}, mc.cores = 14, mc.preschedule = FALSE)
out <- do.call(rbind, res)
write.csv(out, file.path(O, sprintf("E6_learning_%s.csv", arg)), row.names = FALSE)
out$size <- ifelse(out$full, "all", as.character(out$n))
print(aggregate(C ~ size + test, out, function(x) round(c(mean = mean(x), sd = sd(x)), 3)))
