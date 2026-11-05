suppressMessages(library(metafor))
S <- "/groups/adv2105_gp/yichen/tmp/claude-3874/-groups-adv2105-gp-yichen-Yi-com-gastric/272d872d-f4a5-4666-a44a-998a271df7f6/scratchpad"
m <- read.csv(file.path(S, "e24_int_inputs.csv"), stringsAsFactors = FALSE)
pool <- function(d) {
  if (nrow(d) < 3) return(NULL)
  f  <- rma(yi = d$y, sei = d$se, method = "REML")
  hk <- tryCatch(rma(yi = d$y, sei = d$se, method = "REML", test = "knha"), error = function(e) NULL)
  dl <- rma(yi = d$y, sei = d$se, method = "DL")
  lo <- min(dl$ci.lb, if (!is.null(hk)) hk$ci.lb else Inf); hi <- max(dl$ci.ub, if (!is.null(hk)) hk$ci.ub else -Inf)
  p  <- max(dl$pval, if (!is.null(hk)) hk$pval else 0)
  data.frame(k = nrow(d), est = as.numeric(f$b), se = f$se, lo = lo, hi = hi, p = p,
             k_benefit = sum(d$y > 0))
}
out <- list()
for (s in unique(m$score)) for (st in c("all six", "renal", "breast")) {
  d <- m[m$score == s, ]; if (st != "all six") d <- d[d$cancer == st, ]
  r <- pool(d); if (!is.null(r)) out[[length(out) + 1]] <- cbind(score = s, set = st, r)
}
res <- do.call(rbind, out)
write.csv(res, file.path(S, "e24_int_pooled.csv"), row.names = FALSE)
n6 <- res[res$score == "NetImmune" & res$set == "all six", ]
cat(sprintf("check NetImmune all six: %+.3f (%.3f to %.3f) P = %.3f, %d/%d  [Table 4: +0.29 (0.03 to 0.55) 0.034 6/6]\n",
            n6$est, n6$lo, n6$hi, n6$p, n6$k_benefit, n6$k))
