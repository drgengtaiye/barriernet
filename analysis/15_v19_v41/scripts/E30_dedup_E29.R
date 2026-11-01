## E30_dedup_E29 - copy of E29_dissociation.R for E30 Q1 (analysis/E30_prereg.md). No analytic change: the GSE26253
## sample filter enters through the E15 scores of the E30_dedup_E15.R run of the same mode (and E20 scores of the
## E30_dedup_E20.R run); the step-0 survival reference is that run's E20_estimates_used_for_pooling.csv. Outputs are
## redirected to results/E30_dedup_runs/<mode>/. Usage: E30_MODE=off|primary|sens3 E30_RDS=<scratch> Rscript E30_dedup_E29.R
## ---------------------------------------------------------------- original header follows
## E29 - direct test of Immune against CAF_ECM in the single-arm checkpoint cohorts (response) and the ten
## gastric cohorts (survival). Pre-registered in E29_prereg.md (sha256 in E29_prereg.sha256, written before
## this script was run). Scores are loaded and computed exactly as in E5_icb_split.R, E12_icb_expanded.R,
## E22_gastric_icb.R, E15/E20 (patient-level score files); step 0 must reproduce the published marginal
## estimates before anything else is computed.
suppressMessages({library(MultiAssayExperiment); library(SummarizedExperiment); library(matrixStats)
                  library(metafor); library(org.Hs.eg.db); library(survival)})
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
O_IN <- file.path(G, "2/analysis/results")                               # published results (read only)
MODE <- Sys.getenv("E30_MODE", "off"); stopifnot(MODE %in% c("off", "primary", "sens3"))
O <- file.path(O_IN, "E30_dedup_runs", MODE); stopifnot(file.exists(file.path(O, "E20_estimates_used_for_pooling.csv")))
ORDS <- file.path(Sys.getenv("E30_RDS"), MODE); stopifnot(nzchar(Sys.getenv("E30_RDS")))
RV <- file.path(G, "revision/results")
D5 <- file.path(G, "icb")
NB <- 2000
sink(file.path(O, "E29_log.txt"), split = TRUE)
cat("E29 run", format(Sys.time()), "\nprereg sha256:", readLines(file.path(G, "2/analysis/E29_prereg.sha256")), "\n")

MOD <- list(
  CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
  TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
  APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
  IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
  TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
  T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"))
IMM <- c("Immune", "APM", "IFN", "TLS_B", "T_CELL")       # Immune is primary; modules are S6
std <- function(v) (v - mean(v)) / sd(v)

## ------------------------------------------------------------------ response cohorts
## E5 loader and scorer (verbatim logic from E5_icb_split.R)
sym2ens <- suppressMessages(mapIds(org.Hs.eg.db, keys = unique(unlist(MOD)), keytype = "SYMBOL",
                                   column = "ENSEMBL", multiVals = "list"))
load_e5 <- function(nm) {
  x <- suppressWarnings(readRDS(file.path(D5, paste0(nm, ".rds"))))
  ex <- experiments(x); an <- names(ex)
  pick <- c("expr_gene_tpm", "expr")[c("expr_gene_tpm", "expr") %in% an][1]
  M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
  r <- as.character(cd$response); y <- ifelse(r == "R", 1, ifelse(r == "NR", 0, NA)); names(y) <- rownames(cd)
  keep <- intersect(colnames(M), names(y)[!is.na(y)]); M <- M[, keep, drop = FALSE]; y <- y[keep]
  rownames(M) <- sub("\\..*$", "", rownames(M)); M <- M[!duplicated(rownames(M)), , drop = FALSE]
  M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
  Z <- (M - rowMeans(M)) / rowSds(M)
  m <- sapply(names(MOD), function(k) {
    ids <- unlist(lapply(MOD[[k]], function(g) { e <- intersect(sym2ens[[g]], rownames(Z)); if (length(e)) e[1] else NULL }))
    if (length(ids) < 3) rep(NA_real_, ncol(Z)) else colMeans(Z[ids, , drop = FALSE]) })
  S <- as.data.frame(m); S$Immune <- rowMeans(S[, c("APM", "IFN", "TLS_B", "T_CELL")])
  list(S = S, y = unname(y))
}
## E12 scorer (verbatim logic from E12_icb_expanded.R)
score_e12 <- function(M) {
  M <- M[rowSds(M) > 0, , drop = FALSE]; Z <- (M - rowMeans(M)) / rowSds(M)
  S <- as.data.frame(sapply(names(MOD), function(k) {
    g <- intersect(MOD[[k]], rownames(Z)); if (length(g) < 3) rep(NA_real_, ncol(Z)) else colMeans(Z[g, , drop = FALSE]) }))
  S$Immune <- rowMeans(S[, c("APM", "IFN", "TLS_B", "T_CELL")]); S
}
pub <- read.csv(file.path(O_IN, "E12_icb_per_cohort.csv"), stringsAsFactors = FALSE)
E5C <- unique(pub$cohort[pub$set == "E5"]); NEWC <- unique(pub$cohort[pub$set == "new"])
stopifnot(length(E5C) == 11, length(NEWC) == 10)
CO <- readRDS(file.path(O_IN, "E12_cohorts.rds"))
RESP <- list()
for (nm in E5C) { d <- load_e5(nm); RESP[[nm]] <- list(S = d$S, y = d$y, set = "E5", tech = "WT") }
for (nm in NEWC) RESP[[nm]] <- list(S = score_e12(CO[[nm]]$M), y = CO[[nm]]$y, set = "new",
                                    tech = pub$tech[pub$cohort == nm][1])
RESP[["PRJEB25780"]] <- list(S = score_e12(CO$PRJEB25780$M), y = CO$PRJEB25780$y, set = "gastric", tech = "WT")
p6 <- read.csv(file.path(O_IN, "E22_scores_PRJEB60680.csv"), stringsAsFactors = FALSE)
RESP[["PRJEB60680"]] <- list(S = p6[, c(names(MOD), "Immune")], y = p6$y, set = "gastric", tech = "WT")

## ------------------------------------------------------------------ survival cohorts
B2 <- readRDS(file.path(RV, "B2_scores.rds")); S15 <- readRDS(file.path(ORDS, "E15_scores.rds"))
S20 <- readRDS(file.path(ORDS, "E20_scores.rds"))$OS
SURV <- list(TCGA = B2$TCGA, GSE84437 = B2$GSE84437, GSE62254 = B2$GSE62254,
             GSE15459 = S15$GSE15459, GSE34942 = S15$GSE34942, GSE26253 = S15$GSE26253,
             GSE13861 = S20$GSE13861, GSE26901 = S20$GSE26901, GSE26899 = S20$GSE26899, GSE29272 = S20$GSE29272)
SGRP <- setNames(rep(c("benchmark", "locked"), c(3, 7)), names(SURV))

## ------------------------------------------------------------------ step 0: reproduce published marginals
logit1 <- function(y, s) { f <- suppressWarnings(glm(y ~ std(s), family = binomial())); c(coef(f)[2], sqrt(vcov(f)[2, 2])) }
cox1 <- function(d, s) { d$x <- as.numeric(scale(d[[s]])); dd <- d[complete.cases(d[, c("time", "event", "x")]), ]
  co <- summary(coxph(Surv(time, event) ~ x, data = dd))$coefficients; c(co[1, 1], co[1, 3]) }
pubE22 <- read.csv(file.path(O_IN, "E22_PRJEB60680_response.csv"), stringsAsFactors = FALSE)
pubE22 <- pubE22[grepl("^primary", pubE22$analysis), ]
pubS <- read.csv(file.path(O, "E20_estimates_used_for_pooling.csv"), stringsAsFactors = FALSE)
pubS <- pubS[pubS$model == "unadjusted" & grepl("primary", pubS$duplicate_rule), ]
chk <- list()
for (nm in names(RESP)) for (k in c("Immune", "CAF_ECM")) {
  if (all(is.na(RESP[[nm]]$S[[k]]))) { cat("step 0:", nm, k, "not scorable (< 3 genes), as published\n"); next }
  e <- logit1(RESP[[nm]]$y, RESP[[nm]]$S[[k]])
  ref <- if (nm == "PRJEB60680") pubE22[pubE22$score == k, c("logOR", "se")] else
         pub[pub$cohort == nm & pub$score == k & pub$set %in% c("E5", "new", "gastric"), c("logOR", "se")]
  chk[[length(chk) + 1]] <- data.frame(data = "response", cohort = nm, score = k, est = e[1], se = e[2],
    published = ref$logOR[1], published_se = ref$se[1], diff = max(abs(e[1] - ref$logOR[1]), abs(e[2] - ref$se[1])))
}
for (nm in names(SURV)) for (k in c("Immune", "CAF_ECM")) {
  e <- cox1(SURV[[nm]], k); ref <- pubS[pubS$cohort == ifelse(nm == "TCGA", "TCGA", nm) & pubS$score == k, ]
  chk[[length(chk) + 1]] <- data.frame(data = "survival", cohort = nm, score = k, est = e[1], se = e[2],
    published = ref$beta[1], published_se = ref$se[1], diff = max(abs(e[1] - ref$beta[1]), abs(e[2] - ref$se[1])))
}
chk <- do.call(rbind, chk); rownames(chk) <- NULL
write.csv(chk, file.path(O, "E29_step0_reproduction.csv"), row.names = FALSE)
cat("\n== step 0: largest |difference| from published", signif(max(chk$diff), 3), "\n")
print(chk[chk$diff > 1e-6, ], row.names = FALSE)
if (any(chk$diff > 1e-6)) { cat("STEP 0 FAILED - stopping, no E29 estimate computed\n"); sink(); quit(status = 1) }

## ------------------------------------------------------------------ per-cohort contrasts
set.seed(2029)
resp_D <- function(y, sI, sC, boot = TRUE) {
  zi <- std(sI); zc <- std(sC)
  f <- suppressWarnings(glm(y ~ zi + zc, family = binomial())); V <- vcov(f)
  a <- c(bI = coef(f)[["zi"]], bC = coef(f)[["zc"]], D = coef(f)[["zi"]] + coef(f)[["zc"]],
         seD = sqrt(V["zi", "zi"] + V["zc", "zc"] + 2 * V["zi", "zc"]))
  mI <- unname(logit1(y, sI)); mC <- unname(logit1(y, sC)); b <- c(D = mI[1] + mC[1], seD = NA, cov = NA)
  if (boot) {
    idx1 <- which(y == 1); idx0 <- which(y == 0)
    dr <- t(replicate(NB, { i <- c(sample(idx1, replace = TRUE), sample(idx0, replace = TRUE))
      c(logit1(y[i], sI[i])[1], logit1(y[i], sC[i])[1]) }))
    cv <- cov(dr[, 1], dr[, 2]); b["cov"] <- cv; b["seD"] <- sqrt(mI[2]^2 + mC[2]^2 + 2 * cv)
  }
  list(a = a, b = b, r = cor(sI, sC, method = "spearman"))
}
surv_D <- function(d, kI, kC, adj = FALSE, boot = TRUE) {
  cov <- c()
  if (adj) { if (sum(!is.na(d$stage_ord)) > 20) cov <- c(cov, "stage_ord"); if (sum(!is.na(d$age)) > 20) cov <- c(cov, "age") }
  d$zi <- as.numeric(scale(d[[kI]])); d$zc <- as.numeric(scale(d[[kC]]))
  dd <- d[complete.cases(d[, c("time", "event", "zi", "zc", cov)]), ]
  f <- coxph(as.formula(paste("Surv(time, event) ~ zi + zc", if (length(cov)) paste("+", paste(cov, collapse = "+")) else "")), data = dd)
  V <- vcov(f)
  a <- c(bI = coef(f)[["zi"]], bC = coef(f)[["zc"]], D = -(coef(f)[["zi"]] + coef(f)[["zc"]]),
         seD = sqrt(V["zi", "zi"] + V["zc", "zc"] + 2 * V["zi", "zc"]))
  b <- c(D = NA, seD = NA, cov = NA)
  if (boot && !adj) {
    mI <- unname(cox1(d, kI)); mC <- unname(cox1(d, kC)); b["D"] <- -(mI[1] + mC[1])
    e1 <- which(dd$event == 1); e0 <- which(dd$event == 0)
    dr <- t(replicate(NB, { i <- c(sample(e1, replace = TRUE), sample(e0, replace = TRUE)); x <- dd[i, ]
      c(coef(coxph(Surv(time, event) ~ as.numeric(scale(x[[kI]])), data = x)),
        coef(coxph(Surv(time, event) ~ as.numeric(scale(x[[kC]])), data = x))) }))
    cv <- cov(dr[, 1], dr[, 2]); b["cov"] <- cv; b["seD"] <- sqrt(mI[2]^2 + mC[2]^2 + 2 * cv)
  }
  list(a = a, b = b, r = cor(d[[kI]], d[[kC]], method = "spearman", use = "complete.obs"), n = nrow(dd), ev = sum(dd$event))
}
rows <- list()
for (nm in names(RESP)) for (k in IMM) {
  R_ <- RESP[[nm]]; if (all(is.na(R_$S[[k]])) || all(is.na(R_$S$CAF_ECM))) next
  o <- resp_D(R_$y, R_$S[[k]], R_$S$CAF_ECM, boot = (k == "Immune"))
  rows[[length(rows) + 1]] <- data.frame(data = "response", cohort = nm, set = R_$set, tech = R_$tech, immune = k,
    n = length(R_$y), n_resp = sum(R_$y), rho = o$r, bI_a = o$a[["bI"]], bC_a = o$a[["bC"]], D_a = o$a[["D"]], seD_a = o$a[["seD"]],
    D_b = o$b[["D"]], seD_b = o$b[["seD"]], cov_b = o$b[["cov"]])
}
for (nm in names(SURV)) for (k in IMM) for (adj in c(FALSE, TRUE)) {
  if (adj && k != "Immune") next
  d <- SURV[[nm]]; d$cohort <- nm
  o <- surv_D(d, k, "CAF_ECM", adj = adj, boot = (k == "Immune"))
  rows[[length(rows) + 1]] <- data.frame(data = if (adj) "survival_adjusted" else "survival", cohort = nm, set = SGRP[[nm]],
    tech = NA, immune = k, n = o$n, n_resp = o$ev, rho = o$r, bI_a = o$a[["bI"]], bC_a = o$a[["bC"]], D_a = o$a[["D"]],
    seD_a = o$a[["seD"]], D_b = o$b[["D"]], seD_b = o$b[["seD"]], cov_b = o$b[["cov"]])
}
per <- do.call(rbind, rows); rownames(per) <- NULL
write.csv(per, file.path(O, "E29_per_cohort.csv"), row.names = FALSE)

## ------------------------------------------------------------------ pooling (manuscript rule) and sets
pool <- function(y, se, label) {
  ok <- is.finite(y) & is.finite(se); y <- y[ok]; se <- se[ok]
  m <- rma(yi = y, sei = se, method = "DL"); h <- rma(yi = y, sei = se, method = "DL", test = "knha")
  r <- rma(yi = y, sei = se, method = "REML")
  wide_h <- (h$ci.ub - h$ci.lb) > (m$ci.ub - m$ci.lb)
  data.frame(k = length(y), est = m$b[1], lo = if (wide_h) h$ci.lb else m$ci.lb, hi = if (wide_h) h$ci.ub else m$ci.ub,
             p = max(m$pval, h$pval), interval = if (wide_h) "HK" else "DL", I2 = m$I2, n_pos = sum(y > 0),
             reml_est = r$b[1], reml_lo = r$ci.lb, reml_hi = r$ci.ub, reml_p = r$pval)
}
SETS <- list(
  list("Q1 response, 21 cohorts", "response", c(E5C, NEWC)),
  list("S3 response, E5 11", "response", E5C),
  list("S3 response, new 10", "response", NEWC),
  list("S3 response, whole transcriptome", "response", c(E5C, NEWC[NEWC %in% per$cohort[per$tech == "WT"]])),
  list("Q2 survival, 10 gastric", "survival", names(SURV)),
  list("S5 survival, benchmark 3", "survival", names(SURV)[1:3]),
  list("S5 survival, locked 7", "survival", names(SURV)[4:10]),
  list("S5 survival, 10 gastric, adjusted", "survival_adjusted", names(SURV)))
out <- list()
for (s in SETS) for (k in IMM) for (meth in c("a", "b", "bI", "bC")) {
  q <- per[per$data == s[[2]] & per$cohort %in% s[[3]] & per$immune == k, ]
  if (!nrow(q)) next
  if (meth == "b" && (k != "Immune" || all(is.na(q$seD_b)))) next
  if (meth %in% c("bI", "bC") && k != "Immune") next
  y <- switch(meth, a = q$D_a, b = q$D_b, bI = q$bI_a, bC = q$bC_a)
  se <- switch(meth, a = q$seD_a, b = q$seD_b, bI = NA, bC = NA)
  if (meth %in% c("bI", "bC")) next                      # adjusted coefficients are reported per cohort only
  out[[length(out) + 1]] <- cbind(set = s[[1]], immune = k, method = meth, pool(y, se),
                                  rho_median = median(q$rho), rho_min = min(q$rho), rho_max = max(q$rho))
}
pooled <- do.call(rbind, out); rownames(pooled) <- NULL
## Q3 crossover (descriptive): D_r - D_s, independent data sets
q1 <- pooled[pooled$set == "Q1 response, 21 cohorts" & pooled$immune == "Immune" & pooled$method == "a", ]
q2 <- pooled[pooled$set == "Q2 survival, 10 gastric" & pooled$immune == "Immune" & pooled$method == "a", ]
se_from_ci <- function(r) (r$hi - r$lo) / (2 * qnorm(0.975))
x <- q1$est - q2$est; sx <- sqrt(se_from_ci(q1)^2 + se_from_ci(q2)^2)
cross <- data.frame(D_r = q1$est, D_s = q2$est, diff = x, lo = x - 1.96 * sx, hi = x + 1.96 * sx, p = 2 * pnorm(-abs(x / sx)))
write.csv(pooled, file.path(O, "E29_pooled.csv"), row.names = FALSE)
write.csv(cross, file.path(O, "E29_Q3_crossover.csv"), row.names = FALSE)

f <- function(v) sprintf("%+.2f", v)
cat("\n== pooled dominance contrasts (D > 0: Immune stronger for response; D < 0: CAF_ECM stronger for survival)\n")
for (i in seq_len(nrow(pooled))) with(pooled[i, ], cat(sprintf("%-38s %-7s (%s) k=%2d  D %s (%s to %s, %s) P=%.3g  pos %d/%d  REML %s P=%.3g  rho median %.2f [%.2f, %.2f]\n",
  set, immune, method, k, f(est), f(lo), f(hi), interval, p, n_pos, k, f(reml_est), reml_p, rho_median, rho_min, rho_max)))
cat("\n== Q3 crossover (descriptive; log OR vs log HR scales)\n"); print(cross, digits = 3)
cat("\n== S4 gastric checkpoint cohorts (not pooled)\n")
print(per[per$set == "gastric" & per$immune == "Immune", c("cohort", "n", "n_resp", "rho", "bI_a", "bC_a", "D_a", "seD_a", "D_b", "seD_b")], digits = 3, row.names = FALSE)
sink()
