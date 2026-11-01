## E30_dedup_E29b - copy of E29b_posthoc_adjusted.R for E30 Q1: sources E30_dedup_E29.R (same mode) instead of
## E29_dissociation.R; outputs go to results/E30_dedup_runs/<mode>/. No other change.
## ---------------------------------------------------------------- original header follows
## E29b (POST HOC, not in E29_prereg.md; added after the E29 run): pooled mutually adjusted coefficients of
## Immune and CAF_ECM, to describe what drives the E29 contrasts. Data loading reused verbatim from
## E29_dissociation.R (everything before its step 0). Same pooling rule.
src <- readLines("E30_dedup_E29.R"); end <- grep("step 0: reproduce", src)[1] - 1
eval(parse(text = src[setdiff(seq_len(end), grep("^sink\\(", src))]))
pool <- function(y, se) { m <- rma(yi = y, sei = se, method = "DL"); h <- rma(yi = y, sei = se, method = "DL", test = "knha")
  w <- (h$ci.ub - h$ci.lb) > (m$ci.ub - m$ci.lb)
  c(k = length(y), est = m$b[1], lo = if (w) h$ci.lb else m$ci.lb, hi = if (w) h$ci.ub else m$ci.ub, p = max(m$pval, h$pval), n_neg = sum(y < 0)) }
r <- list()
for (nm in c(E5C, NEWC)) { R_ <- RESP[[nm]]; if (all(is.na(R_$S$CAF_ECM))) next
  zi <- std(R_$S$Immune); zc <- std(R_$S$CAF_ECM); f <- suppressWarnings(glm(R_$y ~ zi + zc, family = binomial()))
  r[[nm]] <- data.frame(data = "response", cohort = nm, term = c("Immune", "CAF_ECM"), b = coef(f)[2:3], se = sqrt(diag(vcov(f)))[2:3]) }
for (nm in names(SURV)) { d <- SURV[[nm]]; d$zi <- as.numeric(scale(d$Immune)); d$zc <- as.numeric(scale(d$CAF_ECM))
  f <- coxph(Surv(time, event) ~ zi + zc, data = d)
  r[[nm]] <- data.frame(data = "survival", cohort = nm, term = c("Immune", "CAF_ECM"), b = coef(f), se = sqrt(diag(vcov(f)))) }
r <- do.call(rbind, r); rownames(r) <- NULL
out <- do.call(rbind, lapply(split(r, list(r$data, r$term)), function(q) data.frame(data = q$data[1], term = q$term[1], t(pool(q$b, q$se)))))
write.csv(r, file.path(O, "E29b_posthoc_adjusted_per_cohort.csv"), row.names = FALSE)
write.csv(out, file.path(O, "E29b_posthoc_adjusted_pooled.csv"), row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
