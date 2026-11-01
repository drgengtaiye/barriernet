## E31 Action (analysis/E31_prereg.md): every array meeting the duplicate rule outside GSE26253-GSE62254 is removed from
## the cohort analyzed later, and the main-text pooled survival estimates that use those cohorts are re-run and reported
## beside the current values. Within one group (the three benchmark cohorts) "later" is read as later in the protocol's
## cohort list (TCGA-STAD, GSE84437, GSE62254, GSE15459, GSE34942, GSE26253, GSE13861, GSE26901, GSE26899, GSE29272);
## the protocol does not order cohorts within a group, and this reading was fixed before anything below was run.
##
## Implementation: the patient-level scores of the current (E30 primary) analysis are re-used, as in E29: benchmark
## revision/results/B2_scores.rds; first three locked and extension: the E30_dedup_E15/E20 runs in mode primary, whose
## score files are kept in results/E30_dedup_runs/primary/scores/. Removed arrays are dropped and the Cox models refitted; scores are NOT recomputed without them (each
## cohort's z-scores move by the removal of at most 4 of 56-431 arrays). Step 0 must reproduce the current per-cohort
## and pooled values exactly before the removal is applied. Pooling: the manuscript rule for gastric pools
## (DerSimonian-Laird estimate; wider of Hartung-Knapp and DerSimonian-Laird interval; larger P), as in E20/E29.
## Random-axis null ranks and Harrell's C are not re-run. Deviation from the protocol's wording ("re-run with the E30
## machinery"): scores are not recomputed. Every pooled survival estimate of the main text is covered (v33 fact-check):
## benchmark 3, first three, added four, locked 7 and all 10; every module, axis and CD8A; unadjusted and, where every
## cohort of a set has it (E20 rule), adjusted for stage and age; Immune and CAF_ECM in one model.
suppressMessages({library(survival); library(metafor)})
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
O <- file.path(G, "2/analysis/results")
P <- file.path(O, "E30_dedup_runs/primary")
RDS <- file.path(P, "scores")
sink(file.path(O, "E31_rerun.log"), split = TRUE)
cat("E31 re-run", format(Sys.time()), "\n")

B2 <- readRDS(file.path(G, "revision/results/B2_scores.rds")); S15 <- readRDS(file.path(RDS, "E15_scores.rds"))
S20 <- readRDS(file.path(RDS, "E20_scores.rds"))$OS
SURV <- list(TCGA = B2$TCGA, GSE84437 = B2$GSE84437, GSE62254 = B2$GSE62254,
             GSE15459 = S15$GSE15459, GSE34942 = S15$GSE34942, GSE26253 = S15$GSE26253,
             GSE13861 = S20$GSE13861, GSE26901 = S20$GSE26901, GSE26899 = S20$GSE26899, GSE29272 = S20$GSE29272)
SETS <- list("benchmark 3" = names(SURV)[1:3], "locked, first three" = names(SURV)[4:6],
             "locked, added four" = names(SURV)[7:10], "locked 7" = names(SURV)[4:10], "all 10 gastric" = names(SURV))
SC <- c("NetImmune", "CAF_ECM", "TAM_SUP", "Barrier", "Immune", "APM", "IFN", "TLS_B", "T_CELL", "CD8A")

## ---- arrays to remove (rule i-iii met; GSE26253-GSE62254 stays under E30)
C <- read.csv(file.path(O, "E31_candidates.csv"), stringsAsFactors = FALSE)
C <- C[C$duplicate & !(C$cohort == "GSE62254" & C$comparator == "GSE26253"), ]
C$comparator[C$comparator == "TCGA-STAD"] <- "TCGA"; C$cohort[C$cohort == "TCGA-STAD"] <- "TCGA"
stopifnot(all(match(C$cohort, names(SURV)) < match(C$comparator, names(SURV))))   # comparator = later cohort
DROP <- split(C$best_match, C$comparator)
cat("pairs meeting the rule outside GSE26253-GSE62254:", nrow(C), "; arrays to remove:", length(unique(C$best_match)), "\n")
print(sapply(DROP, function(x) length(unique(x))))
already <- unlist(lapply(names(DROP), function(k) setdiff(DROP[[k]], SURV[[k]]$sample)))
cat("already outside the analysis set (excluded by the E20 rule):", if (length(already)) paste(already, collapse = " ") else "none", "\n")
DROP <- lapply(names(DROP), function(k) intersect(DROP[[k]], SURV[[k]]$sample)); names(DROP) <- sort(unique(C$comparator))

cox1 <- function(d, s, adj = FALSE) {                      # E20 cox1(): adjusted = + stage_ord + age where > 20 recorded
  d$x <- as.numeric(scale(d[[s]])); cov <- c()
  if (adj) { if (sum(!is.na(d$stage_ord)) > 20) cov <- c(cov, "stage_ord"); if (sum(!is.na(d$age)) > 20) cov <- c(cov, "age")
             if (!length(cov)) return(NULL) }
  dd <- d[complete.cases(d[, c("time", "event", "x", cov)]), ]
  f <- as.formula(paste("Surv(time, event) ~ x", if (length(cov)) paste("+", paste(cov, collapse = "+")) else ""))
  co <- summary(coxph(f, data = dd))$coefficients; c(b = co["x", 1], se = co["x", 3], n = nrow(dd), ev = sum(dd$event)) }
joint <- function(d, adj) {
  cov <- c()
  if (adj) { if (sum(!is.na(d$stage_ord)) > 20) cov <- c(cov, "stage_ord"); if (sum(!is.na(d$age)) > 20) cov <- c(cov, "age") }
  d$zi <- as.numeric(scale(d$Immune)); d$zc <- as.numeric(scale(d$CAF_ECM))
  dd <- d[complete.cases(d[, c("time", "event", "zi", "zc", cov)]), ]
  f <- coxph(as.formula(paste("Surv(time, event) ~ zi + zc", if (length(cov)) paste("+", paste(cov, collapse = "+")) else "")), data = dd)
  V <- vcov(f)
  c(bI = coef(f)[["zi"]], seI = sqrt(V["zi", "zi"]), bC = coef(f)[["zc"]], seC = sqrt(V["zc", "zc"]),
    D = -(coef(f)[["zi"]] + coef(f)[["zc"]]), seD = sqrt(V["zi", "zi"] + V["zc", "zc"] + 2 * V["zi", "zc"]))
}
pool <- function(b, se) {
  dl <- rma(yi = b, sei = se, method = "DL"); hk <- rma(yi = b, sei = se, method = "DL", test = "knha")
  w <- (hk$ci.ub - hk$ci.lb) >= (dl$ci.ub - dl$ci.lb)
  c(est = dl$b[1], lo = if (w) hk$ci.lb else dl$ci.lb, hi = if (w) hk$ci.ub else dl$ci.ub, p = max(dl$pval, hk$pval))
}
run <- function(S) {
  per <- do.call(rbind, lapply(names(S), function(k) {
    d <- S[[k]]
    rbind(do.call(rbind, lapply(c(paste(SC, "unadjusted"), paste(SC, "adjusted")), function(q) {
            s <- sub(" .*", "", q); e <- cox1(d, s, adj = grepl("adjusted$", q) && !grepl("unadjusted", q))
            if (is.null(e)) return(NULL)
            data.frame(cohort = k, quantity = q, b = e[["b"]], se = e[["se"]], n = e[["n"]], events = e[["ev"]]) })),
          { j <- joint(d, FALSE); data.frame(cohort = k, quantity = c("Immune | CAF_ECM", "CAF_ECM | Immune", "direct test D"),
                                             b = j[c("bI", "bC", "D")], se = j[c("seI", "seC", "seD")], n = nrow(d), events = sum(d$event)) },
          { j <- joint(d, TRUE); data.frame(cohort = k, quantity = "direct test D, adjusted", b = j[["D"]], se = j[["seD"]],
                                            n = nrow(d), events = sum(d$event)) })
  }))
  rownames(per) <- NULL
  pl <- list()
  for (s in names(SETS)) for (q in unique(per$quantity)) {
    a <- per[per$quantity == q & per$cohort %in% SETS[[s]], ]
    if (grepl(" adjusted$", q) && !all(SETS[[s]] %in% a$cohort)) next          # adjusted only where every cohort has it
    if (grepl("^direct test D, adjusted|^direct", q) && s != "all 10 gastric") next
    if (grepl("\\|", q) && s != "all 10 gastric") next
    pl[[length(pl) + 1]] <- data.frame(set = s, quantity = q, k = nrow(a), n = sum(a$n), events = sum(a$events), t(pool(a$b, a$se)))
  }
  list(per = per, pooled = do.call(rbind, pl))
}

## ---- step 0: reproduce the current values (every per-cohort and pooled row that has a published counterpart)
base <- run(SURV)
base$per$score <- sub(" .*", "", base$per$quantity); base$per$model <- sub("^[^ ]+ ", "", base$per$quantity)
base$pooled$score <- sub(" .*", "", base$pooled$quantity); base$pooled$model <- sub("^[^ ]+ ", "", base$pooled$quantity)
pubS <- read.csv(file.path(P, "E20_estimates_used_for_pooling.csv"), stringsAsFactors = FALSE)
pubS <- pubS[grepl("primary", pubS$duplicate_rule) & pubS$group %in% c("benchmark", "locked (first 3)", "locked (added)"), ]
g1 <- merge(base$per[base$per$model %in% c("unadjusted", "adjusted"), ], pubS, by = c("cohort", "score", "model"))
cat("step 0 per-cohort: rows", nrow(g1), "; max |diff| beta", signif(max(abs(g1$b - g1$beta)), 3), "se", signif(max(abs(g1$se.x - g1$se.y)), 3), "\n")
stopifnot(nrow(g1) >= 150, max(abs(g1$b - g1$beta)) < 1e-6, max(abs(g1$se.x - g1$se.y)) < 1e-6)
pl20 <- read.csv(file.path(P, "E20_pooled.csv"), stringsAsFactors = FALSE)
pl20 <- pl20[grepl("primary", pl20$duplicate_rule), ]
ref <- c("benchmark 3" = "benchmark 3 (reference)", "locked, first three" = "locked (first 3), E15 reference",
         "locked, added four" = "locked (added 4)", "locked 7" = "locked 7 (primary)", "all 10 gastric" = "all 10 gastric")
g2 <- do.call(rbind, lapply(seq_len(nrow(base$pooled)), function(i) { r <- base$pooled[i, ]
  q <- pl20[pl20$set == ref[[r$set]] & pl20$score == r$score & pl20$model == r$model, ]
  if (!nrow(q)) return(NULL); data.frame(set = r$set, quantity = r$quantity, d_est = abs(exp(r$est) - q$HR), d_p = abs(r$p - q$p)) }))
cat("step 0 pooled: rows", nrow(g2), "; max |diff| HR", signif(max(g2$d_est), 3), "P", signif(max(g2$d_p), 3), "\n")
stopifnot(nrow(g2) >= 50, max(g2$d_est) < 1e-6, max(g2$d_p) < 1e-6)
p29 <- read.csv(file.path(P, "E29_pooled.csv"), stringsAsFactors = FALSE)
q <- p29[p29$set == "Q2 survival, 10 gastric" & p29$immune == "Immune" & p29$method == "a", ]
r <- base$pooled[base$pooled$set == "all 10 gastric" & base$pooled$quantity == "direct test D", ]
cat("step 0 direct test: here", round(r$est, 6), "published", round(q$est, 6), "\n"); stopifnot(abs(r$est - q$est) < 1e-6)
p29b <- read.csv(file.path(P, "E29b_posthoc_adjusted_pooled.csv"), stringsAsFactors = FALSE)
for (t in c("Immune", "CAF_ECM")) {
  q <- p29b[p29b$data == "survival" & p29b$term == t, ]; r <- base$pooled[base$pooled$quantity == paste0(t, if (t == "Immune") " | CAF_ECM" else " | Immune"), ]
  cat("step 0 mutually adjusted", t, ": here", round(r$est, 6), "published", round(q$est, 6), "\n"); stopifnot(abs(r$est - q$est) < 1e-6)
}
cat("step 0 passed\n")

## ---- with the E31 removals
S2 <- SURV
for (k in names(DROP)) S2[[k]] <- S2[[k]][!S2[[k]]$sample %in% DROP[[k]], ]
new <- run(S2)
out <- merge(base$pooled[, names(new$pooled)], new$pooled, by = c("set", "quantity"), suffixes = c("_current", "_E31"), sort = FALSE)
hr <- !grepl("direct test", out$quantity)
f <- function(e, l, h, p, isHR) if (isHR) sprintf("%.2f (%.2f-%.2f); P = %.2g", exp(e), exp(l), exp(h), p) else
  sprintf("%+.2f (%.2f to %.2f); P = %.2g", e, l, h, p)
out$current <- mapply(f, out$est_current, out$lo_current, out$hi_current, out$p_current, hr)
out$E31 <- mapply(f, out$est_E31, out$lo_E31, out$hi_E31, out$p_E31, hr)
out$significance_changed <- (out$p_current < 0.05) != (out$p_E31 < 0.05)
out$abs_change <- ifelse(hr, abs(exp(out$est_E31) - exp(out$est_current)), abs(out$est_E31 - out$est_current))
options(width = 220)
print(out[, c("set", "quantity", "k_current", "n_current", "n_E31", "current", "E31", "abs_change", "significance_changed")], row.names = FALSE)
write.csv(out, file.path(O, "E31_rerun.csv"), row.names = FALSE)
write.csv(data.frame(cohort = C$comparator, array = C$best_match, matched_to = paste(C$cohort, C$sample)),
          file.path(O, "E31_removed_arrays.csv"), row.names = FALSE)
cat("pooled estimates compared:", nrow(out), "; largest change:", round(max(out$abs_change), 4), "(", out$set[which.max(out$abs_change)],
    out$quantity[which.max(out$abs_change)], ")\n")
cat("significance changed:", sum(out$significance_changed), ":", paste(out$set[out$significance_changed], out$quantity[out$significance_changed], collapse = "; "), "\n")
NP <- out$quantity %in% c("NetImmune unadjusted", "CAF_ECM unadjusted", "Immune unadjusted", "Barrier unadjusted", "NetImmune adjusted")
cat("largest change among NetImmune/CAF_ECM/Immune/Barrier unadjusted and NetImmune adjusted:", round(max(out$abs_change[NP]), 4), "\n")
sink()
