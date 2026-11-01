## E28 Q1 and Q5 in the gastric benchmark cohorts (E28_prereg.md). Overall survival, univariable HR per SD and
## Harrell's C with the direction fixed in advance, pooled as in E8 (DerSimonian-Laird estimate; wider of the
## Hartung-Knapp and DerSimonian-Laird intervals, larger P).
##   Q1  official ESTIMATE with missing values set to the gene minimum (E28_estimate.R) against the z-mean
##       contrast of the manuscript (ESTsig_diff)
##   Q5  NetImmune and CAF_ECM scored without the other patients of the cohort: unit weights (module mean of
##       log expression), within-sample rank, and a frozen TCGA-STAD reference for the two array cohorts
## Usage: Rscript E28_gastric.R <E28_estimate.R workdir>
suppressMessages({library(survival); library(metafor)})
W <- commandArgs(TRUE)[1]
R <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
TAGS <- c("TCGA", "GSE84437", "GSE62254")
MOD <- list(CAF_ECM = c("COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"),
            TAM_SUP = c("SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"),
            APM = c("B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"),
            IFN = c("STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"),
            TLS_B = c("CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"),
            T_CELL = c("CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"))
G48 <- unique(unlist(MOD))
res <- readRDS(file.path(R, "B2_scores.rds"))
B1 <- lapply(TAGS, function(t) readRDS(file.path(R, paste0("B1_", t, ".rds")))); names(B1) <- TAGS

axes <- function(V) {
  m <- sapply(MOD, function(g) { g <- intersect(g, rownames(V)); if (length(g) >= 3) colMeans(V[g, , drop = FALSE], na.rm = TRUE) else NA })
  m <- as.data.frame(m); m$Barrier <- (m$CAF_ECM + m$TAM_SUP) / 2
  m$Immune <- (m$APM + m$IFN + m$TLS_B + m$T_CELL) / 4; m$NetImmune <- m$Immune - m$Barrier; m
}
e0 <- B1$TCGA$expr; e0 <- e0[intersect(G48, rownames(e0)), ]
mu0 <- rowMeans(e0, na.rm = TRUE); sd0 <- apply(e0, 1, sd, na.rm = TRUE)

hc <- function(tt, ee, r) concordance(Surv(tt, ee) ~ r, reverse = TRUE)$concordance
uni <- list(); cor_rows <- list()
for (t in TAGS) {
  d <- res[[t]]; e <- B1[[t]]$expr; idx <- match(d$sample, colnames(e))
  stopifnot(!anyNA(idx))
  e <- e[, idx]
  rk <- apply(e, 2, function(x) rank(x, na.last = "keep") / sum(!is.na(x)))
  rownames(rk) <- rownames(e)
  g0 <- intersect(names(mu0), rownames(e))
  fz <- (e[g0, , drop = FALSE] - mu0[g0]) / sd0[g0]
  off <- read.csv(file.path(W, paste0(t, "_official_scores.csv")))
  off <- off[match(d$sample, off$sample), ]
  S <- list(cohort_relative = d[, c("NetImmune", "CAF_ECM")],
            unit_weights = axes(e)[, c("NetImmune", "CAF_ECM")],
            within_sample_rank = axes(rk)[, c("NetImmune", "CAF_ECM")])
  if (t != "TCGA") S$frozen_TCGA_STAD <- axes(fz)[, c("NetImmune", "CAF_ECM")]
  for (v in names(S)) for (k in c("NetImmune", "CAF_ECM")) {
    s <- S[[v]][[k]]; ok <- !is.na(d$time) & d$time > 0 & !is.na(d$event) & is.finite(s)
    f <- coxph(Surv(time, event) ~ scale(s)[, 1], data = d[ok, ])
    sg <- if (k == "NetImmune") -1 else 1
    uni[[length(uni) + 1]] <- data.frame(cohort = t, scoring = v, score = k, beta = unname(coef(f)), se = sqrt(vcov(f)[1, 1]),
                                         C = hc(d$time[ok], d$event[ok], sg * s[ok]), r_cohort_relative = cor(s, d[[k]], use = "pair"))
  }
  est <- list(ESTIMATE_zmean = d$ESTsig_diff,
              ESTIMATE_official_na_min = as.numeric(scale(off$imm_na_min)) - as.numeric(scale(off$str_na_min)),
              ESTIMATE_official_as_is = as.numeric(scale(off$imm_as_is)) - as.numeric(scale(off$str_as_is)))
  for (v in names(est)) {
    s <- est[[v]]; ok <- !is.na(d$time) & d$time > 0 & !is.na(d$event) & is.finite(s)
    f <- coxph(Surv(time, event) ~ scale(s)[, 1], data = d[ok, ])
    uni[[length(uni) + 1]] <- data.frame(cohort = t, scoring = v, score = "ESTIMATE", beta = unname(coef(f)),
                                         se = sqrt(vcov(f)[1, 1]), C = hc(d$time[ok], d$event[ok], -s[ok]),
                                         r_cohort_relative = cor(s, d$ESTsig_diff, use = "pair"))
  }
}
U <- do.call(rbind, uni); U$HR <- exp(U$beta)
write.csv(U, file.path(O, "E28_gastric_per_cohort.csv"), row.names = FALSE)

pool <- function(u) {
  dl <- rma(yi = u$beta, sei = u$se, method = "DL"); hk <- rma(yi = u$beta, sei = u$se, method = "DL", test = "knha")
  w <- (hk$ci.ub - hk$ci.lb) >= (dl$ci.ub - dl$ci.lb)
  data.frame(k = nrow(u), HR = exp(dl$b[1]), lo = exp(if (w) hk$ci.lb else dl$ci.lb), hi = exp(if (w) hk$ci.ub else dl$ci.ub),
             p = max(dl$pval, hk$pval), I2 = dl$I2)
}
P <- do.call(rbind, lapply(split(U, list(U$scoring, U$score), drop = TRUE), function(u)
  cbind(scoring = u$scoring[1], score = u$score[1], cohorts = paste(u$cohort, collapse = "/"), pool(u),
        C = paste(sprintf("%.3f", u$C), collapse = " / "))))
rownames(P) <- NULL
write.csv(P, file.path(O, "E28_gastric_pooled.csv"), row.names = FALSE)
options(width = 250)
print(format(U, digits = 3), row.names = FALSE); print(format(P, digits = 3), row.names = FALSE)
