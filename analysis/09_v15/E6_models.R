## E6a - more model classes for the prognostic ceiling (leave-one-cohort-out).
source("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/E6_common.R")
library(parallel)
MODELS <- list(
  "Elastic net (alpha 0.5)" = fit_glmnet(0.5), "Lasso" = fit_glmnet(1), "Ridge" = fit_glmnet(0),
  "PC-Cox" = fit_pcacox, "Screened ridge Cox" = fit_screenridge, "Supervised PC" = fit_superpc,
  "Random survival forest" = fit_rsf, "Boosted Cox (xgboost)" = fit_xgbcox)
jobs <- expand.grid(model = names(MODELS), test = CT, stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(jobs)), function(i) {
  j <- jobs[i, ]; d <- train_test(j$test); t0 <- Sys.time()
  r <- tryCatch(MODELS[[j$model]](d$Xtr, d$ttr, d$etr, d$Xte, seed = 1), error = function(e) { message(j$model, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  data.frame(test = j$test, model = j$model, C = cidx(r, d$tte, d$ete),
             tuned = ifelse(is.null(attr(r, "tuned")), "", as.character(attr(r, "tuned"))),
             minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}, mc.cores = 12, mc.preschedule = FALSE)
out <- do.call(rbind, res)
fixed <- do.call(rbind, lapply(CT, function(t) data.frame(test = t, model = c("NetImmune (fixed)", "ESTIMATE (fixed)"),
  C = c(cidx(-SC[[t]]$NetImmune, TT[[t]], EV[[t]]), cidx(-SC[[t]]$ESTsig_diff, TT[[t]], EV[[t]])), tuned = "", minutes = 0)))
out <- rbind(out, fixed)
write.csv(out, file.path(O, "E6_model_classes.csv"), row.names = FALSE)
W <- reshape(out[, c("test", "model", "C")], idvar = "model", timevar = "test", direction = "wide")
colnames(W) <- sub("^C\\.", "", colnames(W)); W$mean <- rowMeans(W[, CT]); W <- W[order(-W$mean), ]
print(format(W, digits = 3), row.names = FALSE)
print(out[, c("test", "model", "tuned", "minutes")])
