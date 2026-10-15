## E12 - more checkpoint cohorts: does the E5 split (immune modules track response, CAF_ECM does
## not) replicate in cohorts that E5 did not use?
##
## Expected directions, written before any E12 outcome was computed (from the 11-cohort E5 result):
##   CAF_ECM              ~ 0 (pooled log-OR per SD within about +/-0.15; no consistent sign)
##   IFN, Immune, APM, TLS_B, T_CELL   > 0 (higher in responders)
##   NetImmune            > 0, no larger than Immune
##   TAM_SUP              leaning > 0, weaker than the immune modules
##   Barrier              between CAF_ECM and TAM_SUP, near 0
##   C1Q_arm, SPP1_arm    both leaning > 0 (as in E5; no opposition expected any more)
##   CXCL9_SPP1           > 0
## The replication test is the new cohorts pooled on their own: the immune modules should have
## a positive pooled log-OR and CAF_ECM should not.
##
## Entry rules, fixed before scoring (cohort facts in results/E12_cohort_audit.csv, from
## E12a_cohort_audit.R, which never relates a score to response):
##   - checkpoint antibody given without chemotherapy; one pre-treatment sample per patient;
##     response = CR/PR against SD/PD (or R/NR as recorded in the compendium objects);
##   - >= 5 responders and >= 5 non-responders (E5 also asked for >= 25 labelled samples; that
##     stricter rule is kept as a sensitivity set);
##   - not the same patients as an included cohort (identifier match, or expression twins:
##     a best cross-cohort Spearman partner standing > 0.03 clear of the second best);
##   - a module is scored where >= 3 of its genes are measured; axes need all their modules.
## Chemotherapy-plus-checkpoint cohorts (GSE207422, GSE235919) enter only a sensitivity set.
## PRJEB25780 (gastric discovery) is not in any pooled set; it is reported beside the others.
## Scoring exactly as E5: genes row-z within cohort on the log scale, module = mean z of its genes
## present, Barrier = mean(CAF_ECM, TAM_SUP), Immune = mean(APM, IFN, TLS_B, T_CELL),
## NetImmune = Immune - Barrier, CXCL9_SPP1 = z(CXCL9) - z(SPP1); effect = logistic log-odds of
## response per SD. The 11 E5 cohorts are taken from results/E5_icb_per_cohort.csv unchanged.
## Pooling: DerSimonian-Laird random effects (metafor, Wald CI) and the same model with the
## Hartung-Knapp adjustment (test = "knha").
suppressMessages({library(matrixStats); library(metafor)})
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"

MOD <- list(
  CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
  TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
  APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
  IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
  TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
  T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"),
  C1Q_arm = c("C1QA","C1QB","C1QC"),
  SPP1_arm = c("SPP1","TREM2","MARCO","APOE","LGALS3"))
SC <- c("CAF_ECM", "TAM_SUP", "IFN", "APM", "TLS_B", "T_CELL", "Barrier", "Immune", "NetImmune",
        "C1Q_arm", "SPP1_arm", "CXCL9_SPP1")

CO <- readRDS(file.path(O, "E12_cohorts.rds"))
aud <- read.csv(file.path(O, "E12_cohort_audit.csv"), stringsAsFactors = FALSE)

## ---- inclusion, applied mechanically to the audit
DUP <- c(ICB_Jerby_Arnon = "same patients as ICB_Liu (91/96 identifiers match; each is its Liu twin's best expression partner, median gap 0.067; labels agree 90/90)",
         mel_dfci_2019 = "same samples as ICB_Liu (121/121 identifiers); no response field")
GENES <- c(ICB_VanDenEnde = "compendium matrix retains 1,365 protein-coding genes (mostly pseudogene/small-RNA rows); 3 of 48 module genes; no module scorable; not rescuable offline")
why <- sapply(aud$id, function(id) {
  a <- aud[aud$id == id, ]
  if (id %in% names(GENES)) return(GENES[[id]])
  if (id %in% names(DUP)) return(DUP[[id]])
  if (id == "PRJEB25780") return("gastric discovery cohort (reported separately)")
  if (id == "ICB_Padron") return("randomised PRINCE trial, analysed as a trial in v15; 3 responders")
  if (a$n_resp < 5 || a$n_nonresp < 5) return(sprintf("%d responders / %d non-responders (< 5 in a class)", a$n_resp, a$n_nonresp))
  if (a$role == "chemoICB") return("included (chemo-ICB sensitivity only)")
  "included"
})
aud$status <- unname(why)
NEW <- aud$id[aud$status == "included"]; CHEMO <- aud$id[aud$status == "included (chemo-ICB sensitivity only)"]
cat("new cohorts:", NEW, "\nchemo-ICB:", CHEMO, "\n")

## ---- per-cohort effects for the new cohorts (and PRJEB25780, and chemo-ICB)
scores_sym <- function(M) {
  M <- M[rowSds(M) > 0, , drop = FALSE]; Z <- (M - rowMeans(M)) / rowSds(M)
  used <- sapply(MOD, function(g) sum(g %in% rownames(Z)))
  S <- as.data.frame(sapply(names(MOD), function(k) {
    g <- intersect(MOD[[k]], rownames(Z)); if (length(g) < 3) rep(NA_real_, ncol(Z)) else colMeans(Z[g, , drop = FALSE]) }))
  S$Barrier <- rowMeans(S[, c("CAF_ECM", "TAM_SUP")])
  S$Immune <- rowMeans(S[, c("APM", "IFN", "TLS_B", "T_CELL")])
  S$NetImmune <- S$Immune - S$Barrier
  S$CXCL9_SPP1 <- if (all(c("CXCL9", "SPP1") %in% rownames(Z))) Z["CXCL9", ] - Z["SPP1", ] else NA_real_
  list(S = S, used = used)
}
logit_eff <- function(y, s) {
  s <- (s - mean(s)) / sd(s)
  f <- suppressWarnings(glm(y ~ s, family = binomial()))
  c(b = unname(coef(f)[2]), se = unname(sqrt(diag(vcov(f)))[2]))
}
auc <- function(y, v) { r <- rank(v); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
TUM <- c(ICB_Nathanson = "melanoma", ICB_Snyder = "urothelial", ICB_Limagne2 = "lung", ICB_Shiuan = "kidney",
         GSE176307 = "urothelial", GSE161537 = "lung", GSE159067 = "head and neck", GSE215868 = "melanoma",
         GSE93157 = "mixed (lung, HNSCC, melanoma)", brca_iatlas_anders_2022 = "breast",
         GSE207422 = "lung", GSE235919 = "colorectal", PRJEB25780 = "gastric")
rows <- list()
for (nm in c(NEW, CHEMO, "PRJEB25780")) {
  d <- CO[[nm]]; sc <- scores_sym(d$M); a <- aud[aud$id == nm, ]
  for (k in SC) {
    v <- sc$S[[k]]; if (all(is.na(v))) next
    e <- logit_eff(d$y, v)
    rows[[length(rows) + 1]] <- data.frame(cohort = nm, tumour = TUM[[nm]], n = length(d$y), n_resp = sum(d$y),
      score = k, genes = if (k %in% names(sc$used)) sc$used[[k]] else NA, logOR = e[["b"]], se = e[["se"]],
      AUC = auc(d$y, v), tech = a$tech, set = if (nm %in% NEW) "new" else if (nm %in% CHEMO) "chemoICB" else "gastric")
  }
}
## gastric sensitivity: the double log2 used by E7a
d <- CO$PRJEB25780; sc <- scores_sym(log2(d$M + 1))
for (k in SC) { v <- sc$S[[k]]; e <- logit_eff(d$y, v)
  rows[[length(rows) + 1]] <- data.frame(cohort = "PRJEB25780 (E7a double log)", tumour = "gastric", n = length(d$y),
    n_resp = sum(d$y), score = k, genes = if (k %in% names(sc$used)) sc$used[[k]] else NA, logOR = e[["b"]],
    se = e[["se"]], AUC = auc(d$y, v), tech = "WT", set = "gastric_sens") }
new <- do.call(rbind, rows)
old <- read.csv(file.path(O, "E5_icb_per_cohort.csv"), stringsAsFactors = FALSE)
old$tech <- "WT"; old$set <- "E5"
per <- rbind(old[, names(new)], new); rownames(per) <- NULL
write.csv(per, file.path(O, "E12_icb_per_cohort.csv"), row.names = FALSE)

## module coverage on each included panel / matrix
cov <- aud[, c("id", "cancer", "treatment", "platform", "tech", "role", "genes", "n_labelled", "n_resp", "n_nonresp",
               names(MOD), "modules_ok", "status")]
write.csv(cov, file.path(O, "E12_cohort_table.csv"), row.names = FALSE)

## ---- pooling
pool_one <- function(p, set) {
  if (nrow(p) < 2) return(NULL)
  m <- rma(yi = p$logOR, sei = p$se, method = "DL")
  h <- rma(yi = p$logOR, sei = p$se, method = "DL", test = "knha")
  data.frame(set = set, score = p$score[1], k = nrow(p), n = sum(p$n), n_resp = sum(p$n_resp),
             logOR = m$b[1], lo_DL = m$ci.lb, hi_DL = m$ci.ub, p_DL = m$pval,
             lo_HK = h$ci.lb, hi_HK = h$ci.ub, p_HK = h$pval, I2 = m$I2, tau2 = m$tau2,
             n_pos = sum(p$logOR > 0), median_AUC = median(p$AUC))
}
E5C <- unique(old$cohort)
STRICT <- aud$id[aud$id %in% NEW & aud$n_labelled >= 25]
SETS <- list(
  "E5 11 (reference)" = E5C,
  "new only (replication)" = NEW,
  "new only, E5 entry rule (>= 25 labelled)" = STRICT,
  "all" = c(E5C, NEW),
  "all, whole transcriptome" = c(E5C, aud$id[aud$id %in% NEW & aud$tech == "WT"]),
  "all, panel only" = aud$id[aud$id %in% NEW & aud$tech == "panel"],
  "all, without kidney" = setdiff(c(E5C, NEW), c("ICB_Braun", "ICB_Miao1", "ICB_Shiuan")),
  "all + chemo-ICB" = c(E5C, NEW, CHEMO),
  "new, whole transcriptome" = aud$id[aud$id %in% NEW & aud$tech == "WT"])
pooled <- do.call(rbind, lapply(names(SETS), function(s) do.call(rbind, lapply(SC, function(k)
  pool_one(per[per$score == k & per$cohort %in% SETS[[s]], ], s)))))
write.csv(pooled, file.path(O, "E12_icb_pooled.csv"), row.names = FALSE)

## ---- leave-one-cohort-out
loo <- do.call(rbind, lapply(c("new only (replication)", "all"), function(s) do.call(rbind, lapply(SC, function(k) {
  p <- per[per$score == k & per$cohort %in% SETS[[s]], ]
  est <- t(sapply(seq_len(nrow(p)), function(i) {
    m <- rma(yi = p$logOR[-i], sei = p$se[-i], method = "DL"); h <- rma(yi = p$logOR[-i], sei = p$se[-i], method = "DL", test = "knha")
    c(b = m$b[1], pDL = m$pval, pHK = h$pval) }))
  data.frame(set = s, score = k, k = nrow(p), loo_min = min(est[, "b"]), loo_max = max(est[, "b"]),
             max_p_DL = max(est[, "pDL"]), max_p_HK = max(est[, "pHK"]),
             most_influential = p$cohort[which.max(abs(est[, "b"] - rma(yi = p$logOR, sei = p$se, method = "DL")$b[1]))])
}))))
write.csv(loo, file.path(O, "E12_icb_loo.csv"), row.names = FALSE)

## ---- console summary
f <- function(x) sprintf("%+.2f", x)
for (s in c("new only (replication)", "all", "all, whole transcriptome", "all, panel only", "all, without kidney",
            "new only, E5 entry rule (>= 25 labelled)", "all + chemo-ICB", "E5 11 (reference)")) {
  cat("\n==", s, "\n")
  q <- pooled[pooled$set == s, ]
  for (i in seq_len(nrow(q))) with(q[i, ], cat(sprintf("  %-11s k=%2d  %s (DL %s..%s, P=%.3g; HK %s..%s, P=%.3g)  I2=%2.0f%%  pos %d/%d\n",
    score, k, f(logOR), f(lo_DL), f(hi_DL), p_DL, f(lo_HK), f(hi_HK), p_HK, I2, n_pos, k)))
}
cat("\n== leave-one-out\n"); print(loo, digits = 3)
cat("\n== gastric PRJEB25780\n"); print(per[grepl("PRJEB", per$cohort), c("cohort", "score", "logOR", "se", "AUC")], digits = 3)
