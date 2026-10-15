## E5 (cross-cancer) - do the two halves of Barrier behave outside gastric cancer as they did
## in PRJEB25780?
##
## Pre-stated expectation (written before running), from the gastric results:
##   TAM_SUP  - higher in responders (pooled log-OR per SD > 0), like the immune modules
##   CAF_ECM  - lower in responders or null (pooled log-OR <= 0)
##   IFN, Immune - higher in responders (log-OR > 0), the strongest of the six
##   Barrier  - attenuated towards 0, because its halves pull in opposite directions
##   NetImmune - > 0 but no larger than Immune
##
## Added at peer-review request, after the module-level results above had been seen but
## before any outcome was computed for these scores. Definitions and expected directions
## (Zhang et al., Cell 2020; Bill et al., Science 2023, CXCL9:SPP1 polarity):
##   C1Q_arm  = C1QA, C1QB, C1QC                        - tracks immune engagement: log-OR > 0
##   SPP1_arm = SPP1, TREM2, MARCO, APOE, LGALS3        - does not track response, or opposes it
##   CXCL9_SPP1 polarity = z(CXCL9) - z(SPP1) (both log-scale genes, row-z within cohort)
##                                                       - reference; expected log-OR > 0
## Size-matched nulls are drawn for 11-, 5- and 3-gene sets.
##
## Cohorts: the checkpoint cohorts of the companion benchmark that pass its entry rules
## (>= 25 labelled samples, >= 5 per class, >= 3,000 genes), minus ICB_Kim, which is the
## gastric discovery cohort (PRJEB25780), and minus ICB_VanDenEnde, whose matrix lacks the
## module genes (COL1A1, CD8A, STAT1 and the rest are absent), leaving 11 cohorts. Response = R vs NR as recorded in the objects.
## Scoring as in the manuscript: genes row-z scored within cohort, module = mean z of its
## genes present (>= 3), axes from the module means. Effect = log-odds of response per SD of
## the score (logistic regression); AUC oriented so that > 0.5 means higher in responders.
## Pooling: DerSimonian-Laird. Null: 1,000 random 11-gene sets drawn from genes measured in
## all cohorts, the same set in every cohort, pooled the same way and credited with its
## better direction (|pooled log-OR|).
suppressMessages({library(MultiAssayExperiment); library(SummarizedExperiment); library(matrixStats)
                  library(metafor); library(org.Hs.eg.db)})
set.seed(2026)
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
NRAND <- 1000
MOD <- list(
  CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
  TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
  APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
  IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
  TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
  T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"),
  C1Q_arm = c("C1QA","C1QB","C1QC"),
  SPP1_arm = c("SPP1","TREM2","MARCO","APOE","LGALS3"))
sym2ens <- suppressMessages(mapIds(org.Hs.eg.db, keys = unique(unlist(MOD)), keytype = "SYMBOL",
                                   column = "ENSEMBL", multiVals = "list"))
sym2ens <- c(sym2ens, suppressMessages(mapIds(org.Hs.eg.db, keys = "CXCL9", keytype = "SYMBOL",
                                              column = "ENSEMBL", multiVals = "list")))
COH <- c("ICB_Braun", "ICB_Gide", "ICB_Hugo", "ICB_Jung", "ICB_Limagne1", "ICB_Liu", "ICB_Mariathasan",
         "ICB_Miao1", "ICB_Puch", "ICB_Riaz", "ICB_Van_Allen")
TUMOUR <- c(ICB_Braun = "kidney", ICB_Miao1 = "kidney", ICB_Gide = "melanoma", ICB_Hugo = "melanoma",
            ICB_Liu = "melanoma", ICB_Riaz = "melanoma", ICB_Puch = "melanoma", ICB_Van_Allen = "melanoma",
            ICB_Jung = "lung", ICB_Limagne1 = "lung", ICB_Mariathasan = "bladder",
            ICB_VanDenEnde = "oesophageal adenocarcinoma")

load_cohort <- function(nm) {
  x <- suppressWarnings(readRDS(file.path(D, paste0(nm, ".rds"))))
  ex <- experiments(x); an <- names(ex)
  pick <- c("expr_gene_tpm", "expr")[c("expr_gene_tpm", "expr") %in% an][1]
  M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
  r <- as.character(cd$response); y <- ifelse(r == "R", 1, ifelse(r == "NR", 0, NA)); names(y) <- rownames(cd)
  keep <- intersect(colnames(M), names(y)[!is.na(y)]); M <- M[, keep, drop = FALSE]; y <- y[keep]
  rownames(M) <- sub("\\..*$", "", rownames(M)); M <- M[!duplicated(rownames(M)), , drop = FALSE]
  M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
  s <- rowSds(M); list(Z = (M - rowMeans(M)) / s, y = y, assay = pick)
}

scores_of <- function(Z) {
  used <- list()
  m <- sapply(names(MOD), function(k) {
    ids <- unlist(lapply(MOD[[k]], function(g) { e <- intersect(sym2ens[[g]], rownames(Z)); if (length(e)) e[1] else NULL }))
    used[[k]] <<- length(ids)
    if (length(ids) < 3) return(rep(NA_real_, ncol(Z)))
    colMeans(Z[ids, , drop = FALSE])
  })
  S <- as.data.frame(m)
  S$Barrier <- rowMeans(S[, c("CAF_ECM", "TAM_SUP")])
  S$Immune <- rowMeans(S[, c("APM", "IFN", "TLS_B", "T_CELL")])
  S$NetImmune <- S$Immune - S$Barrier
  gi <- function(g) { e <- intersect(sym2ens[[g]], rownames(Z)); if (length(e)) e[1] else NA }
  S$CXCL9_SPP1 <- if (!is.na(gi("CXCL9")) && !is.na(gi("SPP1"))) Z[gi("CXCL9"), ] - Z[gi("SPP1"), ] else NA_real_
  list(S = S, used = unlist(used))
}

logit_eff <- function(y, s) {
  s <- (s - mean(s)) / sd(s)
  f <- suppressWarnings(glm(y ~ s, family = binomial()))
  c(b = unname(coef(f)[2]), se = unname(sqrt(diag(vcov(f)))[2]))
}
auc <- function(y, v) { r <- rank(v); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
dl <- function(b, se) {             # DerSimonian-Laird pooled estimate (closed form, for the null)
  w <- 1 / se^2; fe <- sum(w * b) / sum(w); Q <- sum(w * (b - fe)^2); k <- length(b)
  t2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))); ws <- 1 / (se^2 + t2)
  sum(ws * b) / sum(ws)
}

SC <- c("CAF_ECM", "TAM_SUP", "IFN", "APM", "TLS_B", "T_CELL", "Barrier", "Immune", "NetImmune",
        "C1Q_arm", "SPP1_arm", "CXCL9_SPP1")
L <- list(); rows <- list()
for (nm in COH) {
  d <- load_cohort(nm); sc <- scores_of(d$Z); L[[nm]] <- d
  for (k in SC) {
    v <- sc$S[[k]]; if (all(is.na(v))) next
    e <- logit_eff(d$y, v)
    rows[[length(rows) + 1]] <- data.frame(cohort = nm, tumour = TUMOUR[[nm]], n = length(d$y), n_resp = sum(d$y),
      score = k, genes = if (k %in% names(sc$used)) sc$used[[k]] else NA, logOR = e[["b"]], se = e[["se"]],
      AUC = auc(d$y, v))
  }
  cat(sprintf("%-16s %-26s n=%3d resp=%3d  genes used: %s\n", nm, TUMOUR[[nm]], length(d$y), sum(d$y),
              paste(names(sc$used), sc$used, sep = "=", collapse = " ")))
}
per <- do.call(rbind, rows); rownames(per) <- NULL
write.csv(per, file.path(O, "E5_icb_per_cohort.csv"), row.names = FALSE)

pool_one <- function(p) {
  m <- rma(yi = p$logOR, sei = p$se, method = "DL")
  data.frame(score = p$score[1], k = nrow(p), logOR = m$b[1], lo = m$ci.lb, hi = m$ci.ub, p = m$pval,
             I2 = m$I2, n_pos = sum(p$logOR > 0), median_AUC = median(p$AUC))
}
pooled <- do.call(rbind, lapply(SC, function(k) pool_one(per[per$score == k, ])))
noVDE <- do.call(rbind, lapply(SC, function(k) { r <- pool_one(per[per$score == k & per$tumour != "kidney", ])
  r$set <- "without kidney (Braun, Miao1)"; r }))
pooled$set <- "all 11"
loo <- do.call(rbind, lapply(c("CAF_ECM", "TAM_SUP", "C1Q_arm", "SPP1_arm"), function(k) {
  p <- per[per$score == k, ]
  est <- sapply(seq_len(nrow(p)), function(i) rma(yi = p$logOR[-i], sei = p$se[-i], method = "DL")$b[1])
  data.frame(score = k, loo_min = min(est), loo_max = max(est))
}))

## size-matched null: the same random 11-gene set in every cohort
pool_ids <- Reduce(intersect, lapply(L, function(d) rownames(d$Z)))
cat("null gene pool (measured in all 11 cohorts):", length(pool_ids), "\n")
draw_null <- function(size) replicate(NRAND, {
  ids <- sample(pool_ids, size)
  e <- t(sapply(L, function(d) logit_eff(d$y, colMeans(d$Z[ids, , drop = FALSE]))))
  dl(e[, "b"], e[, "se"])
})
NULLS <- list(`11` = draw_null(11), `5` = draw_null(5), `3` = draw_null(3))
null <- NULLS[["11"]]
obs <- setNames(pooled$logOR, pooled$score)
SIZE <- c(CAF_ECM = "11", TAM_SUP = "11", SPP1_arm = "5", C1Q_arm = "3")
nullsum <- do.call(rbind, lapply(names(SIZE), function(k) { nd <- NULLS[[SIZE[[k]]]]
  data.frame(score = k, set_size = as.integer(SIZE[[k]]), obs = obs[[k]], null_median_abs = median(abs(nd)),
             null_p95_abs = unname(quantile(abs(nd), 0.95)), p_emp = (sum(abs(nd) >= abs(obs[[k]])) + 1) / (NRAND + 1)) }))
write.csv(rbind(pooled, noVDE), file.path(O, "E5_icb_pooled.csv"), row.names = FALSE)
write.csv(loo, file.path(O, "E5_icb_loo.csv"), row.names = FALSE)
write.csv(nullsum, file.path(O, "E5_icb_null_summary.csv"), row.names = FALSE)
write.csv(data.frame(draw = seq_len(NRAND), size11 = NULLS[["11"]], size5 = NULLS[["5"]], size3 = NULLS[["3"]]),
          file.path(O, "E5_icb_null_draws.csv"), row.names = FALSE)
cat("\n=== pooled (DerSimonian-Laird), log-odds of response per SD ===\n")
print(format(rbind(pooled, noVDE), digits = 3), row.names = FALSE)
cat("\n=== leave-one-cohort-out range ===\n"); print(format(loo, digits = 3), row.names = FALSE)
cat("\n=== size-matched null (|pooled log-OR|, 11-, 5- and 3-gene sets) ===\n"); print(format(nullsum, digits = 3), row.names = FALSE)
